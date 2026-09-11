"""
Grid-convergence index following Roache's procedure (ASME V&V 20 / Celik et al.,
J. Fluids Eng. 130, 2008).

The three grids used here do not share a constant refinement ratio, so the
observed order p is not available in closed form and is found by fixed-point
iteration on

    p = |ln|eps32/eps21| + q(p)| / ln(r21),
    q(p) = ln((r21^p - s) / (r32^p - s)),   s = sign(eps32/eps21)

which reduces to the textbook formula when r21 == r32.
"""

from __future__ import annotations

import math


def gci(h1: float, h2: float, h3: float, dx1: float, dx2: float, dx3: float,
        Fs: float = 1.25):
    """
    h1 = finest grid value, h3 = coarsest. dx1 < dx2 < dx3.
    Returns observed order, extrapolated value and the fine-grid GCI in percent.
    """
    r21 = dx2 / dx1
    r32 = dx3 / dx2
    e21 = h2 - h1
    e32 = h3 - h2
    if e21 == 0:
        return math.nan, h1, 0.0
    s = 1.0 if (e32 / e21) > 0 else -1.0

    p = 2.0
    for _ in range(200):
        q = math.log((r21**p - s) / (r32**p - s)) if s != 1 or r21 != r32 else 0.0
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-12:
            p = p_new
            break
        p = 0.5 * p + 0.5 * p_new

    h_ext = (r21**p * h1 - h2) / (r21**p - 1.0)
    ea = abs(e21 / h1)
    gci_fine = Fs * ea / (r21**p - 1.0) * 100.0
    return p, h_ext, gci_fine


if __name__ == "__main__":
    import json
    d = json.load(open("results/data/cfd_grid.json"))
    rec = sorted(d["records"], key=lambda r: r["dx_mm"])
    dx = [r["dx_mm"] for r in rec]
    h = [r["h_quasi_steady"] for r in rec]
    print("grid  dx (mm)   h (W/m2K)")
    for a, b in zip(dx, h):
        print(f"       {a:5.2f}   {b:8.1f}")
    p, he, g = gci(h[0], h[1], h[2], dx[0], dx[1], dx[2])
    print(f"\nthree finest grids: observed order p = {p:.2f}")
    print(f"Richardson value at dx -> 0 : {he:.1f} W/m2-K")
    print(f"fine-grid GCI               : {g:.1f} %")
    d["extra"] = {"p": p, "h_extrapolated": he, "GCI_pct": g}
    json.dump(d, open("results/data/cfd_grid.json", "w"), indent=1, default=float)
    print("\n(written back into results/data/cfd_grid.json)")
