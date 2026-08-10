"""Data-only conductivity inversion using SimPEG sensitivities.

The model is a free conductivity field evaluated by ConductivityNet.  No
anomaly metadata is read here; SimPEG is used only for the survey physics and
the data Jacobian.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import discretize
from tqdm import tqdm
from simpeg import maps
import simpeg.electromagnetics.static.resistivity as dc


def build_simpeg_survey(df):
    sources = []
    ordered = []
    currents = df["I"].to_numpy(float)
    if not np.allclose(currents, currents[0], rtol=1e-5, atol=1e-7):
        raise ValueError("La campaña debe tener una corriente I constante")
    injected_current = float(currents[0])
    for a_idx in np.unique(df["A_idx"].to_numpy()):
        rows = df[df["A_idx"] == a_idx]
        m = rows[["M_x", "M_y", "M_z"]].to_numpy(float)
        n = rows[["N_x", "N_y", "N_z"]].to_numpy(float)
        a = rows[["A_x", "A_y", "A_z"]].iloc[0].to_numpy(float)
        sources.append(
            dc.sources.Pole(
                [dc.receivers.Dipole(m, n)],
                location=a,
                current=injected_current,
            )
        )
        ordered.extend(rows.index.to_list())
    return dc.Survey(sources), ordered


class SimPEGDataForward:
    def __init__(self, df, cell_size=5.0, padding=15.0, mesh_centers=None):
        electrodes = np.vstack([
            df[["A_x", "A_y", "A_z"]].to_numpy(float),
            df[["M_x", "M_y", "M_z"]].to_numpy(float),
            df[["N_x", "N_y", "N_z"]].to_numpy(float),
        ])
        if mesh_centers is not None:
            centers = np.asarray(mesh_centers, dtype=float)
            xs, ys, zs = (np.unique(centers[:, i]) for i in range(3))
            dx = float(np.median(np.diff(xs)))
            dy = float(np.median(np.diff(ys)))
            dz = float(np.median(np.diff(zs)))
            self.mesh = discretize.TensorMesh(
                [np.full(xs.size, dx), np.full(ys.size, dy), np.full(zs.size, dz)],
                x0=np.array([xs.min() - dx / 2, ys.min() - dy / 2, zs.min() - dz / 2]),
            )
        else:
            xmin, ymin = electrodes[:, :2].min(axis=0) - padding
            xmax, ymax = electrodes[:, :2].max(axis=0) + padding
            zmax = max(50.0, float(electrodes[:, 2].max()) + 50.0)
            nx = int(np.ceil((xmax - xmin) / cell_size))
            ny = int(np.ceil((ymax - ymin) / cell_size))
            nz = int(np.ceil(zmax / cell_size))
            self.mesh = discretize.TensorMesh(
                [np.full(nx, cell_size), np.full(ny, cell_size), np.full(nz, cell_size)],
                x0=np.array([xmin, ymin, 0.0]),
            )
        self.centers = np.asarray(self.mesh.cell_centers, dtype=np.float32)
        self.n_cells = self.mesh.n_cells
        survey, ordered = build_simpeg_survey(df)
        self.ordered = np.asarray(ordered, dtype=int)
        try:
            from pymatsolver import Pardiso as Solver
        except ImportError:
            from simpeg.utils import SolverLU as Solver
        self.simulation = dc.Simulation3DCellCentered(
            self.mesh, survey=survey, sigmaMap=maps.IdentityMap(self.mesh), solver=Solver
        )
        # SimPEG's dpred is the measured potential difference for the survey
        # current. Train against the actual voltage column, not resistance.
        self.target = df.loc[self.ordered, "V"].to_numpy(float)

    def predict_and_jacobian(self, sigma):
        sigma = np.asarray(sigma, dtype=float).reshape(-1)
        predicted = np.asarray(self.simulation.dpred(sigma), dtype=float)
        jacobian = self.simulation.getJ(sigma)
        if hasattr(jacobian, "toarray"):
            jacobian = jacobian.toarray()
        return predicted, np.asarray(jacobian, dtype=float)


def train_data_only(
    sigma_net,
    forward,
    df,
    epochs=100,
    lr=2e-4,
    rho_background=100.0,
    background_weight=1e-3,
    # The sensitivity must follow the current conductivity model.  Reusing the
    # Jacobian of the initial homogeneous model drives the optimization toward
    # shallow, overly smooth artifacts.
    refresh_jacobian=True,
    device="cpu",
):
    sigma_net.to(device)
    coords = torch.tensor(forward.centers, dtype=torch.float32, device=device)
    optimizer = torch.optim.Adam(sigma_net.parameters(), lr=lr)
    target = forward.target
    target_scale = max(float(np.median(np.abs(target))), 1e-6)
    fixed_jacobian = None
    for epoch in tqdm(range(epochs), desc="SimPEG-data", unit="epoch"):
        optimizer.zero_grad(set_to_none=True)
        sigma = sigma_net(coords)
        sigma_np = sigma.detach().cpu().numpy().ravel()
        if fixed_jacobian is None or refresh_jacobian:
            predicted, jacobian = forward.predict_and_jacobian(sigma_np)
            if not refresh_jacobian:
                fixed_jacobian = jacobian
        else:
            predicted = np.asarray(forward.simulation.dpred(sigma_np), dtype=float)
            jacobian = fixed_jacobian

        pred_abs = np.maximum(np.abs(predicted), 1e-8)
        target_abs = np.maximum(np.abs(target), 1e-8)
        log_residual = np.log(pred_abs) - np.log(target_abs)
        data_loss = float(np.mean(log_residual ** 2))

        # d mean((log|d_pred|-log|d_obs|)^2) / d sigma.
        dloss_dd = 2.0 * log_residual / len(target)
        dlog_dpred = dloss_dd * np.sign(predicted) / pred_abs
        grad_sigma = jacobian.T @ dlog_dpred

        log_sigma = torch.log(torch.clamp(sigma, min=1e-8))
        bg_loss = torch.mean((log_sigma - np.log(1.0 / rho_background)) ** 2)
        grad_tensor = torch.tensor(grad_sigma, dtype=sigma.dtype, device=device).reshape_as(sigma)
        torch.autograd.backward(
            (sigma, bg_loss),
            (grad_tensor, torch.tensor(background_weight, dtype=sigma.dtype, device=device)),
        )
        torch.nn.utils.clip_grad_norm_(sigma_net.parameters(), 1.0)
        optimizer.step()

        if epoch % 10 == 0 or epoch == epochs - 1:
            print(f"SimPEG-data epoch {epoch}: log-data={data_loss:.4e}")
    return sigma_net
