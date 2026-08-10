import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import torch
import numpy as np
import pandas as pd
import h5py
from torch.utils.data import DataLoader

from models import ConductivityNet, PotentialNet
from physics_informer import PhysicsInformer
from pytorch_dataset import ERTDataset
from train import train_pinn
from simpeg_data_inversion import SimPEGDataForward, train_data_only

try:
    import wandb
except ImportError: 
    wandb = None


def main():
    parser = argparse.ArgumentParser(description="Entrenamiento PINN para ERT 3D")
    parser.add_argument("--w_data", type=float, default=5.0, help="Peso para el data loss")
    parser.add_argument("--w_pde", type=float, default=0.1, help="Peso para el residual PDE")
    parser.add_argument("--w_bc", type=float, default=0.1, help="Peso para condiciones de frontera")
    parser.add_argument("--w_reg", type=float, default=1e-4, help="Peso de regularizacion espacial")
    parser.add_argument("--w_flux", type=float, default=0.1, help="Peso para conservacion de flujo")
    parser.add_argument("--use_wandb", action="store_true", help="Activar logging en Weights & Biases")
    parser.add_argument("--wandb_project", type=str, default="ERT_PINN_3D")
    parser.add_argument("--wandb_name", type=str, default="baseline_training_run")
    parser.add_argument("--csv", type=str, default="dataset_output_test/measurements.csv",
                        help="CSV de mediciones generado junto a campaign.h5")
    parser.add_argument("--epochs_adam", type=int, default=1000)
    parser.add_argument("--epochs_lbfgs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--background_weight", type=float, default=0.1,
                        help="Peso del prior global estimado desde la mediana de Rho_a")
    parser.add_argument("--residual_weight", type=float, default=1.0,
                        help="Peso de la diferencia respecto al fondo homogéneo")
    parser.add_argument("--lr_sigma", type=float, default=None)
    parser.add_argument("--warmup_epochs", type=int, default=300,
                        help="Epocas para incorporar gradualmente las restricciones fisicas")
    parser.add_argument("--legacy_pinn", action="store_true",
                        help="Usar la PINN secundaria antigua en lugar del forward diferenciable")
    args = parser.parse_args()

    if args.use_wandb:
        if wandb is None:
            raise ImportError("wandb no esta instalado, pero --use_wandb fue solicitado.")
        api_key = os.getenv("WANDB_API_KEY")
        if api_key:
            wandb.login(key=api_key)
        wandb.init(
            project=args.wandb_project,
            name=args.wandb_name,
            config={
                "w_data": args.w_data,
                "w_pde": args.w_pde,
                "w_bc": args.w_bc,
                "w_reg": args.w_reg,
                "w_flux": args.w_flux,
            },
        )

    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print(f"Iniciando entrenamiento en: {device}")

    repo_root = Path(__file__).resolve().parents[2]
    csv_filepath = Path(args.csv)
    if not csv_filepath.is_absolute():
        csv_filepath = repo_root / csv_filepath
    gamma = 4.0

    weights = {
        "w_data": args.w_data,
        "w_pde": args.w_pde,
        "w_bc": args.w_bc,
        "w_reg": args.w_reg,
        "w_flux": args.w_flux,
    }

    print("Cargando dataset y generando puntos de colocacion fisicos...")
    dataset = ERTDataset(
        csv_filepath=csv_filepath,
        n_pde=500,
        n_bc_surf=100,
        n_bc_inf=100,
        n_flux=50,
        epsilon=gamma,
    )
    current_I = dataset.current_I
    
    # IMPORTANTE: Para una inversión PINN, solo debemos entrenar sobre UN conjunto 
    # de mediciones (un "survey" específico). Extraemos la primera muestra (idx=0).
    subset = torch.utils.data.Subset(dataset, [0])
    dataloader = DataLoader(subset, batch_size=1, shuffle=True)

    campaign_df = pd.read_csv(csv_filepath)
    rho_scale = float(np.median(np.abs(campaign_df["Rho_a"].to_numpy())))
    rho_values = np.abs(campaign_df["Rho_a"].to_numpy())
    # Keep the inverse field within a robust data-derived range.  The old
    # factor of five allowed 700+ ohm-m artifacts although the campaign spans
    # roughly 70-140 ohm-m in apparent resistivity.
    sigma_min = 1.0 / max(2.0 * rho_scale, 1e-6)
    sigma_max = min(1.0 / max(0.4 * float(rho_values.min()), 1e-6), 0.1)
    sigma_init = 1.0 / max(rho_scale, 1e-6)
    sigma_net = ConductivityNet(
        hidden_layers=4,
        hidden_dim=128,
        sigma_init=sigma_init,
        sigma_min=sigma_min,
        sigma_max=max(sigma_max, sigma_init * 1.5),
    ).to(device)

    if not args.legacy_pinn:
        h5_path = csv_filepath.with_name("campaign.h5")
        with h5py.File(h5_path, "r") as handle:
            mesh_centers = np.asarray(handle["mesh/cell_centers"])
        forward_model = SimPEGDataForward(
            campaign_df, cell_size=5.0, padding=15.0, mesh_centers=mesh_centers
        )
        trained_sigma_net = train_data_only(
            sigma_net=sigma_net,
            forward=forward_model,
            df=campaign_df,
            epochs=args.epochs_adam,
            lr=args.lr,
            rho_background=rho_scale,
            background_weight=args.background_weight,
            refresh_jacobian=True,
            device=device,
        )
        torch.save({
            "sigma_net_state_dict": trained_sigma_net.state_dict(),
            "model_config": {
                "hidden_layers": 4,
                "hidden_dim": 128,
                "sigma_init": sigma_init,
                "sigma_min": sigma_min,
                "sigma_max": max(sigma_max, sigma_init * 1.5),
            },
            "campaign_csv": str(csv_filepath),
        }, "sigma_net.pth")
        print("Inversión forward diferenciable completada. Pesos guardados en sigma_net.pth")
        if args.use_wandb and wandb is not None:
            wandb.finish()
        return

    pot_net = PotentialNet(conductivity_net=sigma_net).to(device)
    informer = PhysicsInformer(sigma_net, pot_net, source_radius=gamma)

    print("Iniciando entrenamiento PINN con todo el dataset...")
    trained_pot_net, trained_sigma_net = train_pinn(
        u_net=pot_net,
        sigma_net=sigma_net,
        informer=informer,
        dataloader=dataloader,
        weights=weights,
        current_I=current_I,
        gamma=gamma,
        num_epochs_adam=args.epochs_adam,
        num_epochs_lbfgs=args.epochs_lbfgs,
        lr=args.lr,
        lr_sigma=args.lr_sigma,
        warmup_epochs=args.warmup_epochs,
        device=device,
        use_wandb=args.use_wandb,
    )

    print("Entrenamiento completado. Guardando pesos...")
    torch.save({
        "sigma_net_state_dict": trained_sigma_net.state_dict(),
        "model_config": {
            "hidden_layers": 4,
            "hidden_dim": 128,
            "sigma_init": sigma_init,
            "sigma_min": sigma_min,
            "sigma_max": max(sigma_max, sigma_init * 1.5),
        },
        "campaign_csv": str(csv_filepath),
    }, "sigma_net.pth")
    torch.save(trained_pot_net.state_dict(), "pot_net.pth")
    if args.use_wandb and wandb is not None:
        wandb.finish()


if __name__ == "__main__":
    main()
