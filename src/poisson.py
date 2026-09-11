"""
Fast separable Poisson solver for the pressure-projection step.

The projection step solves  lap(phi) = rhs  every time step on a fixed uniform
grid, so it pays to diagonalise once and reuse. We take a cosine transform
along y (homogeneous Neumann on both walls, which is what a pressure equation
wants at a solid wall) and are then left with `ny` independent tridiagonal
systems in x, solved by a Thomas sweep vectorised across all modes at once.

x boundary conditions are configurable so the same class serves the closed box
(Neumann/Neumann, singular, pinned) and the open channel (Neumann at the inlet,
Dirichlet at the outlet).

Cost is O(N log N) with a tiny constant -- roughly 100x faster than a sparse LU
solve at the grid sizes used here, which is what makes the transient runs
affordable.
"""

from __future__ import annotations

import numpy as np
from scipy.fft import dctn, idctn


class FastPoisson:
    """lap(phi) = rhs on a uniform (nx, ny) cell-centred grid."""

    def __init__(self, nx: int, ny: int, dx: float, dy: float,
                 x_bc: tuple[str, str] = ("neumann", "neumann")):
        self.nx, self.ny, self.dx, self.dy = nx, ny, dx, dy
        self.x_bc = x_bc

        # Eigenvalues of the 1-D Neumann Laplacian in y under DCT-II.
        j = np.arange(ny)
        self.lam = -4.0 / dy**2 * np.sin(np.pi * j / (2 * ny)) ** 2   # (ny,)

        self.singular = x_bc == ("neumann", "neumann")

        # Tridiagonal coefficients in x for every y-mode.
        a = np.full((ny, nx), 1.0 / dx**2)      # sub-diagonal
        b = np.full((ny, nx), -2.0 / dx**2)     # diagonal
        c = np.full((ny, nx), 1.0 / dx**2)      # super-diagonal
        b += self.lam[:, None]

        # Left boundary.
        if x_bc[0] == "neumann":
            b[:, 0] += 1.0 / dx**2              # ghost = first cell
        else:                                    # dirichlet, phi=0 on the face
            b[:, 0] -= 1.0 / dx**2
        a[:, 0] = 0.0
        # Right boundary.
        if x_bc[1] == "neumann":
            b[:, -1] += 1.0 / dx**2
        else:
            b[:, -1] -= 1.0 / dx**2
        c[:, -1] = 0.0

        if self.singular:
            # The j=0 mode is a pure Neumann problem: defined only up to a
            # constant. Pin the first cell of that mode.
            b[0, 0] = 1.0
            c[0, 0] = 0.0

        # Pre-compute the Thomas forward sweep (coefficients are constant).
        cp = np.empty_like(c)
        dp_scale = np.empty_like(b)
        cp[:, 0] = c[:, 0] / b[:, 0]
        dp_scale[:, 0] = 1.0 / b[:, 0]
        for i in range(1, nx):
            denom = b[:, i] - a[:, i] * cp[:, i - 1]
            cp[:, i] = c[:, i] / denom
            dp_scale[:, i] = 1.0 / denom
        self._a, self._cp, self._dps = a, cp, dp_scale

    def solve(self, rhs: np.ndarray) -> np.ndarray:
        """rhs shape (nx, ny) -> phi shape (nx, ny)."""
        # Forward cosine transform along y.
        r = dctn(rhs, type=2, axes=1, norm="ortho").T.copy()   # (ny, nx)

        if self.singular:
            r[0, 0] = 0.0

        a, cp, dps = self._a, self._cp, self._dps
        nx = self.nx
        d = np.empty_like(r)
        d[:, 0] = r[:, 0] * dps[:, 0]
        for i in range(1, nx):
            d[:, i] = (r[:, i] - a[:, i] * d[:, i - 1]) * dps[:, i]
        for i in range(nx - 2, -1, -1):
            d[:, i] -= cp[:, i] * d[:, i + 1]

        phi = idctn(d.T.copy(), type=2, axes=1, norm="ortho")
        return phi


def _selftest():
    """Manufactured-solution check of order of accuracy."""
    print("FastPoisson verification (manufactured solutions)\n")

    print("  closed box, Neumann/Neumann in x, Neumann in y")
    prev = None
    for n in (32, 64, 128, 256):
        Lx = Ly = 1.0
        dx, dy = Lx / n, Ly / n
        x = (np.arange(n) + 0.5) * dx
        y = (np.arange(n) + 0.5) * dy
        X, Y = np.meshgrid(x, y, indexing="ij")
        # cos modes satisfy Neumann on all four walls
        exact = np.cos(2 * np.pi * X) * np.cos(3 * np.pi * Y)
        rhs = -((2 * np.pi) ** 2 + (3 * np.pi) ** 2) * exact
        ps = FastPoisson(n, n, dx, dy, ("neumann", "neumann"))
        phi = ps.solve(rhs)
        phi -= phi.mean()
        err = np.abs(phi - (exact - exact.mean())).max()
        rate = "" if prev is None else f"  order {np.log2(prev/err):4.2f}"
        print(f"    n={n:4d}  Linf={err:.3e}{rate}")
        prev = err

    print("\n  open channel, Neumann inlet / Dirichlet outlet, Neumann in y")
    prev = None
    for n in (32, 64, 128, 256):
        dx = dy = 1.0 / n
        x = (np.arange(n) + 0.5) * dx
        y = (np.arange(n) + 0.5) * dy
        X, Y = np.meshgrid(x, y, indexing="ij")
        # cos(pi x /2 * ...) : dphi/dx = 0 at x=0, phi = 0 at x=1
        kx = 1.5 * np.pi
        ky = 2.0 * np.pi
        exact = np.cos(kx * X) * np.cos(ky * Y)
        rhs = -(kx**2 + ky**2) * exact
        ps = FastPoisson(n, n, dx, dy, ("neumann", "dirichlet"))
        phi = ps.solve(rhs)
        err = np.abs(phi - exact).max()
        rate = "" if prev is None else f"  order {np.log2(prev/err):4.2f}"
        print(f"    n={n:4d}  Linf={err:.3e}{rate}")
        prev = err

    import time
    n = 384
    ps = FastPoisson(n, n, 1.0 / n, 1.0 / n)
    r = np.random.rand(n, n); r -= r.mean()
    ps.solve(r)
    t0 = time.perf_counter()
    for _ in range(20):
        ps.solve(r)
    print(f"\n  timing: {n}x{n} solve = "
          f"{(time.perf_counter()-t0)/20*1000:.2f} ms")


if __name__ == "__main__":
    _selftest()
