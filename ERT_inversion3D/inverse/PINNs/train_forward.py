"""Train conductivity directly through the differentiable ERT forward solve."""

import torch
from tqdm import tqdm


def train_forward_inversion(
    sigma_net,
    forward_model,
    dataloader,
    epochs_adam=1000,
    epochs_lbfgs=100,
    lr=1e-3,
    background_weight=0.1,
    residual_weight=1.0,
    device="cpu",
):
    sigma_net.to(device)
    optimizer = torch.optim.Adam(sigma_net.parameters(), lr=lr)
    batch = next(iter(dataloader))["data"]
    batch = {k: (v[0] if isinstance(v, torch.Tensor) and v.ndim > 2 else v) for k, v in batch.items()}
    r_m = batch["r_m"].to(device)
    r_n = batch["r_n"].to(device)
    source = batch["source"].to(device)
    target = batch["delta_v"].to(device)
    scale = target.abs().mean().clamp_min(1e-3)
    grid = forward_model.cell_coords
    with torch.no_grad():
        background_sigma = torch.full_like(grid[:, :1], sigma_net.sigma_init)
        background_prediction = forward_model.predict(background_sigma, source, r_m, r_n)
        residual_scale = (target - background_prediction).abs().mean().clamp_min(1e-4)

    def objective():
        sigma = sigma_net(grid)
        prediction = forward_model.predict(sigma, source, r_m, r_n)
        # ERT voltages span a wide dynamic range.  Log residuals prevent the
        # largest voltages from dominating the inversion; the sign term keeps
        # the measured polarity.
        log_prediction = torch.log(torch.abs(prediction).clamp_min(1e-6))
        log_target = torch.log(torch.abs(target).clamp_min(1e-6))
        data_loss = torch.mean((log_prediction - log_target) ** 2)
        sign_loss = torch.mean(torch.relu(-prediction * target) / (scale ** 2))
        data_loss = data_loss + 0.1 * sign_loss
        residual_prediction = prediction - background_prediction
        residual_target = target - background_prediction
        residual_loss = torch.mean(((residual_prediction - residual_target) / residual_scale) ** 2)
        sigma_grid = sigma.reshape(forward_model.nx, forward_model.ny, forward_model.nz)
        log_sigma_grid = torch.log(torch.clamp(sigma_grid, min=1e-8))
        background_loss = torch.mean((log_sigma_grid - torch.log(torch.tensor(
            sigma_net.sigma_init, device=device))) ** 2)
        smooth_loss = (
            (log_sigma_grid[1:] - log_sigma_grid[:-1]).abs().mean()
            + (log_sigma_grid[:, 1:] - log_sigma_grid[:, :-1]).abs().mean()
            + (log_sigma_grid[:, :, 1:] - log_sigma_grid[:, :, :-1]).abs().mean()
        ) / 3.0
        total_loss = (
            data_loss
            + residual_weight * residual_loss
            + background_weight * background_loss
            + 1e-2 * smooth_loss
        )
        return total_loss, data_loss, prediction

    for epoch in tqdm(range(epochs_adam), desc="Adam-forward"):
        optimizer.zero_grad()
        loss, _, _ = objective()
        if not torch.isfinite(loss):
            raise ValueError("La pérdida forward contiene NaN o Inf")
        loss.backward()
        torch.nn.utils.clip_grad_norm_(sigma_net.parameters(), 1.0)
        optimizer.step()

    if epochs_lbfgs > 0:
        lbfgs = torch.optim.LBFGS(
            sigma_net.parameters(), lr=0.1, max_iter=1, history_size=20,
            line_search_fn="strong_wolfe",
        )
        for _ in tqdm(range(epochs_lbfgs), desc="LBFGS-forward"):
            def closure():
                lbfgs.zero_grad()
                loss, _, _ = objective()
                loss.backward()
                return loss
            lbfgs.step(closure)

    os.makedirs("checkpoints", exist_ok=True)
    torch.save({
        "sigma_net_state_dict": sigma_net.state_dict(),
        "model_config": {
            "hidden_layers": len(sigma_net.mlp.blocks),
            "hidden_dim": sigma_net.mlp.output_layer.in_features,
            "sigma_init": sigma_net.sigma_init,
            "sigma_min": sigma_net.sigma_min,
            "sigma_max": sigma_net.sigma_max,
        },
        "forward_bounds": forward_model.domain_bounds,
        "forward_grid": (forward_model.nx, forward_model.ny, forward_model.nz),
    }, "checkpoints/final_forward_checkpoint.pth")
    return sigma_net
