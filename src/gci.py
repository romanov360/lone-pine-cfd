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
    Returns the observed order, the extrapolated value, the fine-grid GCI in
    percent, and whether convergence was OSCILLATORY.

    That last flag matters. Roache's iteration still returns numbers when the
    two differences have opposite signs, but they describe a sequence that is
    bouncing rather than converging, and quoting the extrapolated value as
    though an asymptotic range had been established would be wrong. Callers
    must check it.
    """
    r21 = dx2 / dx1
    r32 = dx3 / dx2
    e21 = h2 - h1
    e32 = h3 - h2
    if e21 == 0:
        return math.nan, h1, 0.0, False
    s = 1.0 if (e32 / e21) > 0 else -1.0
    oscillatory = s < 0

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
    return p, h_ext, gci_fine, oscillatory


if __name__ == "__main__":
    import json
    d = json.load(open("results/data/cfd_grid.json"))
    rec = sorted(d["records"], key=lambda r: r["dx_mm"])
    dx = [r["dx_mm"] for r in rec]
    h = [r["h_quasi_steady"] for r in rec]
    print("grid  dx (mm)   h (W/m2K)")
    for a, b in zip(dx, h):
        print(f"       {a:5.2f}   {b:8.1f}")
    p, he, g, osc = gci(h[0], h[1], h[2], dx[0], dx[1], dx[2])
    print(f"\nthree finest grids ({dx[0]:.2f}, {dx[1]:.2f}, {dx[2]:.2f} mm):")
    print(f"  observed order p          : {p:.2f}")
    print(f"  Richardson value, dx -> 0 : {he:.1f} W/m2-K")
    print(f"  fine-grid GCI             : {g:.1f} %")
    if osc:
        print("\n  *** OSCILLATORY CONVERGENCE ***")
        print("  The two grid differences have opposite signs, so this sequence")
        print("  is bouncing rather than converging. No asymptotic range has")
        print("  been established and the extrapolated value above must NOT be")
        print("  quoted as a converged answer. Here the cause is physical: the")
        print("  plume is genuinely unsteady at Ra ~ 1e9, and successive runs")
        print("  differ by more than the discretisation does.")
    d["extra"] = {"p": p, "h_extrapolated": he, "GCI_pct": g,
                  "oscillatory": osc,
                  "grids_mm": [dx[0], dx[1], dx[2]]}
    json.dump(d, open("results/data/cfd_grid.json", "w"), indent=1, default=float)
    print("\n(written back into results/data/cfd_grid.json)")
