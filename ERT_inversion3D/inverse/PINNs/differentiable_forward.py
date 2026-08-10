"""Differentiable finite-volume forward model for pole ERT data."""

from __future__ import annotations

import torch


class DifferentiableERTForward:
    """Solve the variable-conductivity DC equation on a regular cell grid.

    The only trainable input is sigma on the cells.  Source and receiver
    locations come from the survey; no anomaly geometry is used.
    """

    def __init__(self, bounds, grid=(24, 24, 12), padding=(50.0, 50.0, 50.0), device="cpu"):
        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = bounds
        self.nx, self.ny, self.nz = grid
        self.device = device
        self.domain_bounds = tuple(bounds)
        px, py, pz = padding
        self.x0, self.x1 = self.x_min - px, self.x_max + px
        self.y0, self.y1 = self.y_min - py, self.y_max + py
        self.z0, self.z1 = self.z_min, self.z_max + pz
        self.dx = (self.x1 - self.x0) / self.nx
        self.dy = (self.y1 - self.y0) / self.ny
        self.dz = (self.z1 - self.z0) / self.nz
        x = torch.linspace(self.x0 + self.dx / 2, self.x1 - self.dx / 2, self.nx)
        y = torch.linspace(self.y0 + self.dy / 2, self.y1 - self.dy / 2, self.ny)
        z = torch.linspace(self.z0 + self.dz / 2, self.z1 - self.dz / 2, self.nz)
        X, Y, Z = torch.meshgrid(x, y, z, indexing="ij")
        self.cell_coords = torch.stack([X, Y, Z], dim=-1).reshape(-1, 3).to(device)
        self.n_cells = self.cell_coords.shape[0]
        self._build_stencil()

    def _build_stencil(self):
        edges_i, edges_j, edge_coeff = [], [], []
        boundary_i, boundary_coeff = [], []
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    idx = ix * self.ny * self.nz + iy * self.nz + iz
                    if ix + 1 < self.nx:
                        edges_i.append(idx)
                        edges_j.append((ix + 1) * self.ny * self.nz + iy * self.nz + iz)
                        edge_coeff.append(self.dy * self.dz / self.dx)
                    if iy + 1 < self.ny:
                        edges_i.append(idx)
                        edges_j.append(ix * self.ny * self.nz + (iy + 1) * self.nz + iz)
                        edge_coeff.append(self.dx * self.dz / self.dy)
                    if iz + 1 < self.nz:
                        edges_i.append(idx)
                        edges_j.append(ix * self.ny * self.nz + iy * self.nz + iz + 1)
                        edge_coeff.append(self.dx * self.dy / self.dz)
                    if ix in (0, self.nx - 1):
                        boundary_i.append(idx)
                        boundary_coeff.append(self.dy * self.dz / (self.dx / 2))
                    if iy in (0, self.ny - 1):
                        boundary_i.append(idx)
                        boundary_coeff.append(self.dx * self.dz / (self.dy / 2))
                    if iz == self.nz - 1:
                        boundary_i.append(idx)
                        boundary_coeff.append(self.dx * self.dy / (self.dz / 2))
        self.edge_i = torch.tensor(edges_i, dtype=torch.long, device=self.device)
        self.edge_j = torch.tensor(edges_j, dtype=torch.long, device=self.device)
        self.edge_coeff = torch.tensor(edge_coeff, dtype=torch.float32, device=self.device)
        self.boundary_i = torch.tensor(boundary_i, dtype=torch.long, device=self.device)
        self.boundary_coeff = torch.tensor(boundary_coeff, dtype=torch.float32, device=self.device)

    def _nearest_cells(self, coords):
        scaled = torch.stack(
            [
                (coords[:, 0] - self.x0) / self.dx,
                (coords[:, 1] - self.y0) / self.dy,
                (coords[:, 2] - self.z0) / self.dz,
            ],
            dim=1,
        )
        ijk = torch.round(scaled - 0.5).long()
        ijk[:, 0].clamp_(0, self.nx - 1)
        ijk[:, 1].clamp_(0, self.ny - 1)
        ijk[:, 2].clamp_(0, self.nz - 1)
        return ijk[:, 0] * self.ny * self.nz + ijk[:, 1] * self.nz + ijk[:, 2]

    def _matrix(self, sigma):
        """Build -div(sigma grad) with top Neumann and distant grounded sides."""
        n = self.n_cells
        A = torch.zeros((n, n), dtype=sigma.dtype, device=sigma.device)
        sigma = sigma.reshape(-1)
        si = sigma[self.edge_i]
        sj = sigma[self.edge_j]
        conductance = self.edge_coeff.to(sigma.dtype) * (2 * si * sj / (si + sj + 1e-12))
        flat_indices = torch.cat([
            self.edge_i * n + self.edge_i,
            self.edge_j * n + self.edge_j,
            self.edge_i * n + self.edge_j,
            self.edge_j * n + self.edge_i,
        ])
        flat_values = torch.cat([conductance, conductance, -conductance, -conductance])
        A.view(-1).index_add_(0, flat_indices, flat_values)

        # Fix the arbitrary potential reference of the all-Neumann system.
        A[0, :] = 0.0
        A[:, 0] = 0.0
        A[0, 0] = 1.0
        return A

    def predict(self, sigma, source_coords, r_m, r_n, current=1.0):
        """Return predicted voltage differences for every survey row."""
        A = self._matrix(sigma)
        source_A = source_coords[:, :3]
        unique_A, inverse = torch.unique(source_A, dim=0, return_inverse=True)
        source_cells = self._nearest_cells(unique_A)
        B = torch.zeros((self.n_cells, unique_A.shape[0]), dtype=sigma.dtype, device=sigma.device)
        source_columns = torch.arange(unique_A.shape[0], device=sigma.device)
        B[source_cells, source_columns] = current

        # A pole source returns at infinity.  Approximate that return with a
        # charge-balanced flux over the remote outer boundary.
        ix = torch.arange(self.nx, device=sigma.device).view(-1, 1, 1)
        iy = torch.arange(self.ny, device=sigma.device).view(1, -1, 1)
        iz = torch.arange(self.nz, device=sigma.device).view(1, 1, -1)
        outer = ((ix == 0) | (ix == self.nx - 1) |
                 (iy == 0) | (iy == self.ny - 1) |
                 (iz == self.nz - 1)).reshape(-1)
        outer[source_cells] = False
        sink_cells = torch.nonzero(outer, as_tuple=False).flatten()
        if sink_cells.numel() > 0:
            B[sink_cells, :] -= current / sink_cells.numel()
        potentials = torch.linalg.solve(A, B)
        m_cells = self._nearest_cells(r_m)
        n_cells = self._nearest_cells(r_n)
        return (potentials[m_cells, inverse] - potentials[n_cells, inverse]).unsqueeze(-1)
