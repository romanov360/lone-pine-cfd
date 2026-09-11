"""
Figures for the study. Writes PNGs into results/figures/.

Palette is the validated categorical set (blue / orange / aqua / yellow /
magenta), checked with the data-viz validator for CVD separation on a light
surface. Colour always encodes the same thing across every figure: BLUE is the
still box, ORANGE is the moving creek.
"""

from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

FIG = "results/figures"
os.makedirs(FIG, exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW, MAGENTA = ("#2a78d6", "#eb6834", "#1baf7a",
                                       "#eda100", "#e87ba4")
VIOLET, RED, GREEN = "#4a3aa7", "#e34948", "#008300"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
SURF = "#fcfcfb"
GRID = "#e6e5e1"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": GRID, "axes.linewidth": 1.0,
    "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 11.5, "axes.titleweight": "600",
    "axes.titlelocation": "left", "axes.titlepad": 10,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "xtick.major.size": 0, "ytick.major.size": 0,
    "grid.color": GRID, "grid.linewidth": 0.8,
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
    "figure.dpi": 150,
})


def _style(ax, xlabel=None, ylabel=None, title=None, grid="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    if grid:
        ax.grid(True, axis=grid, zorder=0)
        ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return ax


def _minutes(ax):
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v/60:g}"))


def save(fig, name):
    path = f"{FIG}/{name}.png"
    fig.savefig(path, bbox_inches="tight", dpi=170)
    plt.close(fig)
    print(f"  {path}")
    return path


# ---------------------------------------------------------------------------


def fig_resistance_chain():
    """Where the heat transfer resistance actually sits."""
    import sys
    sys.path.insert(0, "src")
    from correlations import h_external, h_internal, resistance_split
    from geometry import BOTTLE_500 as B

    Ts, Tinf = 22.0, 10.0
    hi = h_internal(B.D_inner, B.H_inner, 30.0, Ts)
    cases = [("Still box\n(no current)", 0.0),
             ("Slow creek\n0.05 m/s", 0.05),
             ("Typical creek\n0.35 m/s", 0.35),
             ("Fast creek\n1.5 m/s", 1.5)]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    labels, ins, wall, outs, Us = [], [], [], [], []
    for name, U in cases:
        ho = h_external(B.D_outer, B.H_outer, U, Ts, Tinf)
        sp = resistance_split(B, ho, hi)
        labels.append(name)
        ins.append(sp["inside"] * 100)
        wall.append(sp["wall"] * 100)
        outs.append(sp["outside"] * 100)
        Us.append(1.0 / (sp["R_total"] * B.A_outer))

    y = np.arange(len(cases))[::-1]
    h = 0.55
    ax.barh(y, ins, h, color=AQUA, zorder=3, label="Water inside the bottle")
    ax.barh(y, wall, h, left=ins, color=VIOLET, zorder=3,
            label="Glass wall")
    ax.barh(y, outs, h, left=np.array(ins) + np.array(wall), color=ORANGE,
            zorder=3, label="Outside film (what the creek changes)")
    for i, yy in enumerate(y):
        c = [ins[i], wall[i], outs[i]]
        x = 0
        for j, v in enumerate(c):
            if v > 7:
                ax.text(x + v / 2, yy, f"{v:.0f}%", ha="center", va="center",
                        color="white", fontsize=9, fontweight="600", zorder=4)
            x += v
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_xlim(0, 100)
    _style(ax, xlabel="share of total thermal resistance  (%)", grid="x")
    ax.set_title("Only the orange slice is the creek's to change")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)

    ax2.barh(y, Us, h, color=[BLUE, ORANGE, ORANGE, ORANGE], zorder=3)
    for i, yy in enumerate(y):
        ax2.text(Us[i] + 2, yy, f"{Us[i]:.0f}", va="center", fontsize=9.5,
                 color=INK, fontweight="600")
    ax2.set_yticks(y); ax2.set_yticklabels([])
    ax2.set_xlim(0, max(Us) * 1.22)
    _style(ax2, xlabel="overall U  (W/m$^2$K)", grid="x")
    ax2.set_title("Still water to a fast creek: +35% on U")
    fig.subplots_adjust(wspace=0.08)
    return save(fig, "01_resistance_chain")


def fig_cooling_curves():
    import sys
    sys.path.insert(0, "src")
    from lumped import Scenario, simulate

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3), sharey=True)

    runs = [
        ("Box, 5 L", dict(U=0.0, bath_volume=0.005), BLUE, "-"),
        ("Box, 20 L", dict(U=0.0, bath_volume=0.020), BLUE, "--"),
        ("Box, 200 L", dict(U=0.0, bath_volume=0.200), BLUE, ":"),
        ("Creek, 0.35 m/s", dict(U=0.35, bath_volume=None), ORANGE, "-"),
    ]
    for name, kw, c, ls in runs:
        sc = Scenario(name, **kw)
        r = simulate(sc, t_end=5400)
        ax.plot(r["t"], r["T_bottle"], color=c, ls=ls, label=name, zorder=3)
        if kw["bath_volume"]:
            ax.plot(r["t"], r["T_bath"], color=c, ls=ls, lw=1.0, alpha=0.45,
                    zorder=2)
    ax.axhline(10, color=INK3, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(5300, 10.5, "creek / bath at 10 $^\\circ$C", ha="right",
            fontsize=8.5, color=INK3)
    ax.set_xlim(0, 5400); ax.set_ylim(8, 41)
    _minutes(ax)
    _style(ax, xlabel="time (minutes)", ylabel="bottle temperature ($^\\circ$C)",
           title="A small box stalls; the creek does not")
    ax.legend(loc="upper right")

    for U, c, ls, lab in ((0.0, BLUE, "-", "No current (still)"),
                          (0.02, ORANGE, ":", "Creek 0.02 m/s"),
                          (0.1, ORANGE, "--", "Creek 0.1 m/s"),
                          (0.35, ORANGE, "-", "Creek 0.35 m/s"),
                          (1.5, ORANGE, "-.", "Creek 1.5 m/s")):
        sc = Scenario("", U=U, bath_volume=None)
        r = simulate(sc, t_end=5400)
        ax2.plot(r["t"], r["T_bottle"], color=c, ls=ls, label=lab, zorder=3)
    ax2.axhline(10, color=INK3, lw=1, ls=(0, (4, 3)), zorder=1)
    ax2.set_xlim(0, 3600)
    _minutes(ax2)
    _style(ax2, xlabel="time (minutes)",
           title="With an unlimited bath, velocity barely matters")
    ax2.legend(loc="upper right")
    return save(fig, "02_cooling_curves")


def fig_h_vs_velocity():
    import sys
    sys.path.insert(0, "src")
    from correlations import (_FORCED, h_external, h_forced_crossflow,
                              h_free_vertical)
    from geometry import BOTTLE_500 as B

    Ts, Tinf = 22.0, 10.0
    U = np.logspace(-3, 0.55, 200)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))

    # Band = spread of the six correlations AFTER each is blended with
    # buoyancy, i.e. the uncertainty in the quantity actually plotted. The raw
    # forced-convection spread would dive away from the curve at low velocity,
    # where every pure-forced correlation tends to zero but the real transfer
    # is floored by natural convection.
    from correlations import h_mixed
    hfree0 = h_free_vertical(B.H_outer, Ts, Tinf, D=B.D_outer)
    band = np.array([[h_mixed(h_forced_crossflow(B.D_outer, u, Ts, Tinf, m),
                              hfree0) for u in U] for m in _FORCED])
    ax.fill_between(U, band.min(0), band.max(0), color=ORANGE, alpha=0.18,
                    lw=0, zorder=2,
                    label="spread of six correlations")
    hmix = np.array([h_external(B.D_outer, B.H_outer, u, Ts, Tinf) for u in U])
    ax.plot(U, hmix, color=ORANGE, zorder=4,
            label="mixed convection (what a creek gives)")
    hfree = h_free_vertical(B.H_outer, Ts, Tinf, D=B.D_outer)
    ax.axhline(hfree, color=BLUE, zorder=3,
               label=f"still box, buoyancy only ({hfree:.0f})")
    ax.axhline(B.U_wall, color=VIOLET, ls=(0, (5, 3)), lw=1.6, zorder=3)
    ax.text(4.3e-3, B.U_wall * 1.06,
            f"the glass wall alone is worth only {B.U_wall:.0f} W/m$^2$K",
            fontsize=8.5, color=VIOLET)

    import os
    cpath = "results/data/cfd_creek.json"
    if not os.path.exists(cpath):
        cpath = "results/data/cfd_creek_partial.json"
    if os.path.exists(cpath):
        rec = json.load(open(cpath))["records"]
        us = np.array([max(r["U"], 1.1e-3) for r in rec])
        hs = np.array([r["h_quasi_steady"] for r in rec])
        sd = np.array([r.get("h_std", 0.0) for r in rec])
        ax.errorbar(us, hs, yerr=sd, fmt="o", color=INK, ms=6.5, zorder=6,
                    mec=SURF, mew=1.5, capsize=3, elinewidth=1.2,
                    label="this study's CFD")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(4e-3, 3.5); ax.set_ylim(250, 7000)
    _style(ax, xlabel="creek velocity (m/s)",
           ylabel="outside film coefficient  $h$  (W/m$^2$K)",
           title="The creek wins the film coefficient easily -- and it still barely helps",
           grid="both")
    ax.legend(loc="upper left", fontsize=8.6)
    return save(fig, "03_h_vs_velocity")


def fig_aeration():
    import sys
    sys.path.insert(0, "src")
    from aeration import h_aerated
    from geometry import BOTTLE_500 as B

    voids = np.linspace(0, 0.30, 40)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3))

    for U, c, ls, lab in ((0.0, BLUE, "-", "still box"),
                          (0.05, AQUA, "--", "slow creek 0.05 m/s"),
                          (0.35, ORANGE, "-", "creek 0.35 m/s"),
                          (1.0, RED, "-.", "fast creek 1.0 m/s")):
        net = [h_aerated(B.D_outer, B.H_outer, U, a, 22, 10)["net_effect_pct"]
               for a in voids]
        ax.plot(voids * 100, net, color=c, ls=ls, label=lab, zorder=3)
    ax.axhline(0, color=INK3, lw=1, zorder=2)
    _style(ax, xlabel="entrained air, void fraction (%)",
           ylabel="change in $h$ versus bubble-free water (%)",
           title="Bubbles help still water and hurt fast water")
    ax.legend(loc="upper right")

    r = [h_aerated(B.D_outer, B.H_outer, 0.35, a, 22, 10) for a in voids]
    ax2.plot(voids * 100, [x["property_effect_pct"] for x in r], color=VIOLET,
             label="property penalty (lower $\\rho c_p$ and $k$)", zorder=3)
    ax2.plot(voids * 100, [x["agitation_effect_pct"] for x in r], color=AQUA,
             label="agitation benefit (bubble-driven stirring)", zorder=3)
    ax2.plot(voids * 100, [x["net_effect_pct"] for x in r], color=ORANGE,
             lw=2.6, label="net", zorder=4)
    ax2.axhline(0, color=INK3, lw=1, zorder=2)
    _style(ax2, xlabel="entrained air, void fraction (%)",
           ylabel="change in $h$ (%)",
           title="The two effects, separated (creek at 0.35 m/s)")
    ax2.legend(loc="lower left")
    return save(fig, "04_aeration")


def fig_dissolved_oxygen():
    import sys
    sys.path.insert(0, "src")
    from props import bubbly_mixture, do_saturation, water, water_with_dissolved_o2

    base = water(10.0)
    o2 = water_with_dissolved_o2(10.0, saturation=1.0)
    names = ["density\n$\\rho$", "heat capacity\n$c_p$",
             "conductivity\n$k$", "viscosity\n$\\mu$"]
    do_eff = [abs(getattr(o2, f) / getattr(base, f) - 1) * 100
              for f in ("rho", "cp", "k", "mu")]
    b1 = bubbly_mixture(10.0, 0.01)
    bub = [abs(getattr(b1, f) / getattr(base, f) - 1) * 100
           for f in ("rho", "cp", "k", "mu")]

    fig, ax = plt.subplots(figsize=(8.0, 4.3))
    x = np.arange(4)
    ax.bar(x - 0.19, np.maximum(do_eff, 1e-7), 0.36, color=AQUA, zorder=3,
           label="dissolved oxygen, 0 -> 100% saturated (11.3 mg/L)")
    ax.bar(x + 0.19, np.maximum(bub, 1e-7), 0.36, color=ORANGE, zorder=3,
           label="entrained air bubbles, 1% by volume")
    for xi, v in zip(x - 0.19, do_eff):
        ax.text(xi, max(v, 1e-7) * 1.5, f"{v:.1e}%" if v > 0 else "0",
                ha="center", fontsize=8, color=INK2, rotation=0)
    for xi, v in zip(x + 0.19, bub):
        ax.text(xi, max(v, 1e-7) * 1.5, f"{v:.2f}%", ha="center", fontsize=8,
                color=INK2)
    ax.set_yscale("log"); ax.set_ylim(1e-7, 30)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9)
    _style(ax, ylabel="change in the property (%, log scale)",
           title="Dissolved oxygen is ~10,000x too small to matter; bubbles are not")
    ax.legend(loc="upper left")
    return save(fig, "05_dissolved_oxygen")


def fig_validation():
    """The three benchmarks the solver was checked against."""
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.9))

    ax = axes[0]
    try:
        lid = json.load(open("results/data/validation_cavity.json"))["lid_cavity"]
        for rec, c in zip(lid, (BLUE, ORANGE, AQUA)):
            ax.plot(rec["u"], rec["y"], color=c, zorder=3,
                    label=f"Re = {rec['Re']:.0f}")
            ax.plot(rec["ghia"], rec["y"], "o", color=c, ms=4.5, mfc="none",
                    mew=1.3, zorder=4)
        ax.plot([], [], "o", color=INK3, ms=4.5, mfc="none", mew=1.3,
                label="Ghia et al. (1982)")
    except FileNotFoundError:
        pass
    ax.set_xlim(-0.5, 1.05); ax.set_ylim(0, 1)
    _style(ax, xlabel="$u$ on the vertical centreline", ylabel="$y/L$",
           title="Lid-driven cavity", grid="both")
    ax.legend(loc="upper left", fontsize=8.4)

    ax = axes[1]
    try:
        hc = json.load(open("results/data/validation_heated.json"))["heated_cavity"]
        Ra = [r["Ra"] for r in hc]
        ax.plot(Ra, [r["ref"] for r in hc], "o", color=INK3, ms=7, mfc="none",
                mew=1.5, zorder=3, label="de Vahl Davis (1983)")
        ax.plot(Ra, [r["Nu"] for r in hc], "x", color=ORANGE, ms=8, mew=2.2,
                zorder=4, label="this solver")
        for r in hc:
            ax.annotate(f"{r['error_pct']:+.2f}%", (r["Ra"], r["Nu"]),
                        textcoords="offset points", xytext=(7, -11),
                        fontsize=8, color=INK2)
    except FileNotFoundError:
        pass
    ax.set_xscale("log"); ax.set_yscale("log")
    _style(ax, xlabel="Rayleigh number", ylabel="average Nusselt number",
           title="Buoyancy: heated cavity", grid="both")
    ax.legend(loc="upper left", fontsize=8.4)

    ax = axes[2]
    import os
    cp = "results/data/validation_cyl.json"
    if not os.path.exists(cp):
        cp = "results/data/validation_cyl_partial.json"
    try:
        cyl = json.load(open(cp))["cylinder"]
        Re = [r["Re"] for r in cyl]
        ax.plot(Re, [r["Nu_ChurchillBernstein"] for r in cyl], "-o", color=INK3,
                ms=6, mfc="none", mew=1.5, zorder=3,
                label="Churchill-Bernstein")
        ax.plot(Re, [r["Nu"] for r in cyl], "x", color=ORANGE, ms=8, mew=2.2,
                zorder=4, label="this solver")
        for r in cyl:
            ax.annotate(f"{r['Nu_error_pct']:+.0f}%", (r["Re"], r["Nu"]),
                        textcoords="offset points", xytext=(6, -12),
                        fontsize=8, color=INK2)
    except (FileNotFoundError, KeyError):
        pass
    ax.set_xscale("log")
    _style(ax, xlabel="Reynolds number", ylabel="Nusselt number",
           title="Forced convection: cylinder (Pr = 0.7)", grid="both")
    ax.legend(loc="upper left", fontsize=8.4)
    fig.subplots_adjust(wspace=0.32)
    return save(fig, "06_validation")


def fig_grid_convergence():
    """
    Honest grid study: convergence here is not monotone, because the flow
    itself is unsteady at Ra ~ 1e9. Plotting a Richardson extrapolation would
    imply an asymptotic range this data does not have.
    """
    try:
        d = json.load(open("results/data/cfd_grid.json"))
    except FileNotFoundError:
        return None
    rec = sorted(d["records"], key=lambda r: -r["dx_mm"])
    dx = np.array([r["dx_mm"] for r in rec])
    h = np.array([r["h_quasi_steady"] for r in rec])
    sd = np.array([float(np.nanstd(np.asarray(r["series"]["h_eff"])[
        np.asarray(r["series"]["t"]) > 0.6 * r["t_end"]])) for r in rec])

    import sys
    sys.path.insert(0, "src")
    from correlations import h_free_vertical
    from geometry import BOTTLE_500 as B
    hc = h_free_vertical(B.H_outer, 25.0, 10.0, D=B.D_outer)

    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    fine = dx <= 1.5
    lo = (h[fine] - sd[fine]).min()
    hi = (h[fine] + sd[fine]).max()
    ax.axhspan(lo, hi, color=ORANGE, alpha=0.13, lw=0, zorder=2,
               label=f"1.5 mm and finer: {h[fine].mean():.0f} "
                     f"$\\pm$ {(hi-lo)/2:.0f} W/m$^2$K")
    ax.errorbar(dx, h, yerr=sd, fmt="o-", color=BLUE, ms=7, capsize=4,
                elinewidth=1.4, lw=2.0, zorder=4,
                label="CFD, bars = unsteadiness of the plume")
    ax.axhline(hc, color=INK3, ls=":", zorder=3,
               label=f"Churchill-Chu, infinite medium ({hc:.0f})")

    for x, y, s_ in zip(dx, h, sd):
        ax.annotate(f"{x:.2f} mm", (x, y), textcoords="offset points",
                    xytext=(0, -(s_ + 16)), fontsize=8.2, color=INK2,
                    ha="center")
    ax.annotate("under-resolved", (dx[0], h[0]), textcoords="offset points",
                xytext=(14, 2), fontsize=8.5, color=INK2, ha="left")

    ax.set_xlim(2.25, 0.55)
    ax.set_ylim(340, 545)
    _style(ax, xlabel="cell size $\\Delta x$ (mm)   --   finer to the right",
           ylabel="quasi-steady $h$  (W/m$^2$K)",
           title="The mesh is not the limit -- the flow's own unsteadiness is")
    ax.legend(loc="upper left", fontsize=8.6)
    ax.text(0.60, hc - 34,
            "the remaining gap to the correlation\nis confinement, not "
            "discretisation error",
            fontsize=8.3, color=INK3, ha="left", va="top")
    return save(fig, "07_grid_convergence")


def fig_crossover():
    import sys
    sys.path.insert(0, "src")
    from lumped import Scenario, simulate, time_to

    vols = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5, np.inf]
    Us = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
    Z = np.full((len(vols), len(Us)), np.nan)
    for i, V in enumerate(vols):
        for j, U in enumerate(Us):
            sc = Scenario("", U=U, bath_volume=None if np.isinf(V) else V)
            t = time_to(simulate(sc, t_end=3 * 3600), 15.0)
            Z[i, j] = t / 60.0 if np.isfinite(t) else np.nan

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "seq", ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#256abf",
                "#184f95", "#0d366b"])
    im = ax.imshow(Z, cmap=cmap, aspect="auto", origin="lower",
                   vmin=np.nanmin(Z), vmax=np.nanmax(Z))
    for i in range(len(vols)):
        for j in range(len(Us)):
            v = Z[i, j]
            txt = "never" if not np.isfinite(v) else f"{v:.0f}"
            frac = 0 if not np.isfinite(v) else (v - np.nanmin(Z)) / (np.nanmax(Z) - np.nanmin(Z))
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.5,
                    color="white" if frac > 0.55 else INK,
                    fontweight="600" if not np.isfinite(v) else "normal")
    ax.set_xticks(range(len(Us)))
    ax.set_xticklabels([f"{u:g}" for u in Us])
    ax.set_yticks(range(len(vols)))
    ax.set_yticklabels(["inf" if np.isinf(v) else f"{v*1000:g} L" for v in vols])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("current past the bottle (m/s)")
    ax.set_ylabel("bath volume")
    ax.set_title("Minutes to cool 500 mL from 40 $^\circ$C to 15 $^\circ$C")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("minutes", color=INK2)
    cb.outline.set_visible(False)
    return save(fig, "08_crossover")


def fig_fields():
    """
    Temperature fields from the head-to-head runs.

    The two rows have different domains -- a closed 24 x 30 cm box and an
    open 45 x 26 cm channel -- so the row heights are set from the true
    aspect ratios rather than forcing them into one grid, which would
    misrepresent the geometry.
    """
    import glob

    import matplotlib.gridspec as gridspec

    files = sorted(glob.glob("results/data/field_*.npz"),
                   key=lambda f: 0 if "still_box" in f else 1)
    if not files:
        return None

    rows = []
    for f in files:
        d = np.load(f)
        name = "Still box" if "still_box" in f else "Creek, 0.02 m/s"
        dx = float(d["dx"])
        # Only the matched snapshot times, so the two rows are directly
        # comparable. The final-state field is dropped: the two runs ended at
        # different times (360 s for the box, 300 s for the creek), so putting
        # them side by side would invite exactly the wrong comparison.
        panels = [(t, d[f"snapT_{t}"]) for t in (20, 120, 300)
                  if f"snapT_{t}" in d]
        rows.append({"name": name, "dx": dx, "outer": d["outer"],
                     "panels": panels})

    ncol = max(len(r["panels"]) for r in rows)
    aspects = [r["panels"][0][1].shape[1] / r["panels"][0][1].shape[0]
               for r in rows]          # height / width of each domain
    W = 3.35 * ncol
    heights = [W / ncol * a for a in aspects]
    fig = plt.figure(figsize=(W + 1.1, sum(heights) + 1.5))
    gs = gridspec.GridSpec(len(rows), ncol, figure=fig,
                           height_ratios=aspects, hspace=0.42, wspace=0.28,
                           left=0.055, right=0.9, top=0.9, bottom=0.06)

    vmin, vmax = 10.0, 15.5
    im = None
    for i, r in enumerate(rows):
        for j, (t, T) in enumerate(r["panels"]):
            ax = fig.add_subplot(gs[i, j])
            nx, ny = T.shape
            dx = r["dx"]
            Tm = np.ma.masked_where(r["outer"], T)
            im = ax.imshow(Tm.T, origin="lower", cmap="RdYlBu_r",
                           extent=[0, nx * dx * 100, 0, ny * dx * 100],
                           vmin=vmin, vmax=vmax, interpolation="bilinear")
            ax.contour(np.linspace(0, nx * dx * 100, nx),
                       np.linspace(0, ny * dx * 100, ny), Tm.T,
                       levels=np.arange(10.15, vmax, 0.35), colors="k",
                       linewidths=0.25, alpha=0.3)
            xs, ys = np.where(r["outer"])
            ax.add_patch(plt.Rectangle(
                (xs.min() * dx * 100, ys.min() * dx * 100),
                (xs.max() - xs.min() + 1) * dx * 100,
                (ys.max() - ys.min() + 1) * dx * 100,
                fc="#cfcec9", ec=INK, lw=1.0, zorder=5))
            ax.set_title(f"t = {t} s", fontsize=9.5, color=INK2)
            if j == 0:
                ax.set_ylabel(f"{r['name']}\ncm", fontsize=9.5, color=INK)
            ax.set_xlabel("cm", fontsize=8)
            ax.tick_params(labelsize=7.5)
            ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)

    cax = fig.add_axes([0.915, 0.12, 0.016, 0.72])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("water temperature ($^\\circ$C)", color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle("The still box warms the water the bottle sits in; "
                 "the creek carries it away",
                 x=0.055, y=0.965, ha="left", fontsize=12.5,
                 fontweight="600", color=INK)
    return save(fig, "09_fields")


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["all"]
    reg = {"resistance": fig_resistance_chain, "cooling": fig_cooling_curves,
           "hvel": fig_h_vs_velocity, "aeration": fig_aeration,
           "do": fig_dissolved_oxygen, "validation": fig_validation,
           "grid": fig_grid_convergence, "crossover": fig_crossover,
           "fields": fig_fields}
    for k, fn in reg.items():
        if "all" in which or k in which:
            fn()
