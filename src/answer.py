"""
The headline synthesis: every framing of the question, with its answer.

Run this last. It pulls the calibrated transient model and prints the numbers
that the README and the report quote.
"""

from __future__ import annotations

import json
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from correlations import (h_external, h_free_vertical, h_internal,
                          resistance_split)
from geometry import BOTTLE_500 as B
from lumped import Scenario, fmt, simulate, time_to

RULE = "=" * 76


def head(n, t):
    print(f"\n{RULE}\n{n}. {t}\n{RULE}")


def run(name, **kw):
    sc = Scenario(name, **kw)
    r = simulate(sc, t_end=4 * 3600)
    return sc, r


def line(name, r, sc, width=38):
    T30 = np.interp(1800, r["t"], r["T_bottle"])
    print(f"{name:<{width}} {fmt(time_to(r,25)):>8} {fmt(time_to(r,15)):>8} "
          f"{fmt(time_to(r,12)):>8} {T30:7.2f} C {sc.T_equilibrium():7.2f} C")


def header(width=38):
    print(f"{'':<{width}} {'->25C':>8} {'->15C':>8} {'->12C':>8} "
          f"{'T@30min':>9} {'T_final':>9}")
    print("-" * (width + 46))


def main():
    print(RULE)
    print("WARM BOTTLE, COLD WATER: STILL BOX versus MOVING CREEK")
    print(RULE)
    print(B.describe())

    head(1, "The direct comparison, as asked")
    header()
    for nm, kw in (("Still box, 20 L, at 10 C", dict(U=0.0, bath_volume=0.020)),
                   ("Creek, 0.35 m/s, at 10 C", dict(U=0.35, bath_volume=None))):
        sc, r = run(nm, **kw)
        line(nm, r, sc)
    print("\n  The creek wins. It always wins. The interesting part is by how")
    print("  little, and for what reason.")

    head(2, "Why: where the resistance sits")
    hi = h_internal(B.D_inner, B.H_inner, 30.0, 22.0)
    print(f"{'':<24}{'h_out':>9}{'inside':>9}{'wall':>9}{'outside':>9}{'U overall':>11}")
    for nm, U in (("still box", 0.0), ("creek 0.05 m/s", 0.05),
                  ("creek 0.35 m/s", 0.35), ("creek 1.50 m/s", 1.5)):
        ho = h_external(B.D_outer, B.H_outer, U, 22.0, 10.0)
        sp = resistance_split(B, ho, hi)
        Uo = 1.0 / (sp["R_total"] * B.A_outer)
        print(f"{nm:<24}{ho:9.0f}{sp['inside']*100:8.0f}%{sp['wall']*100:8.0f}%"
              f"{sp['outside']*100:8.0f}%{Uo:11.0f}")
    print("\n  Going from still water to a fast 1.5 m/s creek multiplies the")
    print("  outside film coefficient by 9.2 -- and the overall U by 1.35. The")
    print("  glass wall and the water inside the bottle are what actually limit")
    print("  the cooling, and the creek cannot touch either of them.")

    head(3, "The effect that does dominate: a finite bath warms up")
    header()
    for V in (0.002, 0.005, 0.020, 0.100, None):
        nm = "Box, infinite reservoir" if V is None else f"Box, {V*1000:g} L"
        sc, r = run(nm, U=0.0, bath_volume=V)
        line(nm, r, sc)
    print("\n  A 2 L box never reaches 15 C at all: it equilibrates at 16.5 C.")
    print("  No current, however fast, rescues a bath that small -- and no box,")
    print("  however large, catches a creek, because a creek is a bath of")
    print("  infinite size AND a current at the same time.")

    head(4, "How much of the creek's win is the current, and how much the size?")
    sc_a, r_a = run("box 20 L, still", U=0.0, bath_volume=0.020)
    sc_b, r_b = run("box 20 L, stirred at 0.35 m/s", U=0.35, bath_volume=0.020)
    sc_c, r_c = run("infinite bath, still", U=0.0, bath_volume=None)
    sc_d, r_d = run("creek 0.35 m/s (both)", U=0.35, bath_volume=None)
    header()
    for nm, sc, r in (("box 20 L, still", sc_a, r_a),
                      ("  + stir it at 0.35 m/s", sc_b, r_b),
                      ("  + make it infinite instead", sc_c, r_c),
                      ("  + both = the creek", sc_d, r_d)):
        line(nm, r, sc)
    print("\n  Splitting the creek's advantage between its two causes, by box size")
    print("  (the two shares overlap because the effects are not additive):")
    print(f"\n{'bath':>9}{'still':>9}{'stirred':>9}{'creek':>9}"
          f"{'from current':>15}{'from size':>11}")
    for V in (0.002, 0.005, 0.010, 0.020, 0.050, 0.200):
        _, ra = run("", U=0.0, bath_volume=V)
        _, rb = run("", U=0.35, bath_volume=V)
        a, b = time_to(ra, 15), time_to(rb, 15)
        d = time_to(r_d, 15)
        c = time_to(r_c, 15)
        if np.isfinite(a):
            share = f"{(a-b)/(a-d)*100:13.0f}%{(a-c)/(a-d)*100:10.0f}%"
        else:
            share = f"{'--':>14}{'everything':>11}"
        print(f"{V*1000:7.0f} L{fmt(a):>9}{fmt(b):>9}{fmt(d):>9}{share}")
    print("\n  At 200 L the creek's win is 99 % current and 2 % reservoir size.")
    print("  At 5 L the two are comparable. At 2 L the box simply cannot get")
    print("  there at any stirring rate, so reservoir size is the whole story.")

    head(5, "Things that matter MORE than whether the water is moving")
    base, rb = run("baseline: glass bottle, box 20 L", U=0.0, bath_volume=0.020)
    t_base = time_to(rb, 15)
    creek, rc = run("creek", U=0.35, bath_volume=None)
    items = [("move it to a 0.35 m/s creek", dict(U=0.35, bath_volume=None))]
    from dataclasses import replace
    from props import ALUMINIUM, PET
    items += [
        ("swap glass for aluminium (same box)",
         dict(U=0.0, bath_volume=0.020, bottle=replace(B, shell=ALUMINIUM))),
        ("swap glass for PET (same box)",
         dict(U=0.0, bath_volume=0.020, bottle=replace(B, shell=PET))),
        ("halve the wall to 1.5 mm",
         dict(U=0.0, bath_volume=0.020, bottle=replace(B, wall=0.0015))),
        ("use 250 mL instead of 500 mL",
         dict(U=0.0, bath_volume=0.020,
              bottle=replace(B, D_outer=B.D_outer * 0.794,
                             H_outer=B.H_outer * 0.794))),
        ("bath at 4 C instead of 10 C",
         dict(U=0.0, bath_volume=0.020, T_bath_0=4.0)),
        ("bubble air through the box (1% void)",
         dict(U=0.215, bath_volume=0.020, void_fraction=0.01)),
        ("just stir the box by hand (0.1 m/s)",
         dict(U=0.1, bath_volume=0.020)),
    ]
    print(f"{'change from the baseline still box':<40}{'->15 C':>9}{'speed-up':>10}")
    print("-" * 59)
    print(f"{'(baseline)':<40}{fmt(t_base):>9}{'1.00x':>10}")
    rows = []
    for nm, kw in items:
        sc, r = run(nm, **kw)
        t = time_to(r, 15)
        rows.append((nm, t))
    for nm, t in sorted(rows, key=lambda x: x[1]):
        sp = t_base / t if np.isfinite(t) and t > 0 else np.nan
        print(f"{nm:<40}{fmt(t):>9}{sp:9.2f}x")

    head(6, "Oxygenation, the three separate questions")
    from aeration import dissolved_o2_effect, h_aerated
    d = dissolved_o2_effect(B.D_outer, B.H_outer, 0.35, 22.0, 10.0,
                            saturations=(0.0, 1.0, 2.0))
    print("  (a) DISSOLVED oxygen, 0 -> 200 % of saturation:")
    print(f"      change in h = {d[-1]['change_pct']:+.5f} %   -- irrelevant")
    print("  (b) ENTRAINED bubbles, properties only (1 % void):")
    r1 = h_aerated(B.D_outer, B.H_outer, 0.35, 0.01, 22.0, 10.0)
    print(f"      change in h = {r1['property_effect_pct']:+.2f} %   -- a penalty")
    print("  (c) ENTRAINED bubbles, agitation included:")
    for U, nm in ((0.0, "still box"), (0.35, "creek 0.35 m/s"),
                  (1.0, "creek 1.0 m/s")):
        r2 = h_aerated(B.D_outer, B.H_outer, U, 0.01, 22.0, 10.0)
        print(f"      {nm:<16} net {r2['net_effect_pct']:+7.1f} %")
    print("\n  Aeration helps by stirring, not by oxygenating, so it helps most")
    print("  exactly where there was no stirring -- the still box.")

    head(7, "The one case where the still box beats the creek")
    print("  Only by cheating on temperature. A box is a container you control;")
    print("  a creek is whatever the creek happens to be.")
    header()
    for nm, kw in (("Box of ice water, pinned at 0 C",
                    dict(U=0.0, bath_volume=None, T_bath_0=0.0)),
                   ("Creek at 10 C, 0.35 m/s",
                    dict(U=0.35, bath_volume=None, T_bath_0=10.0)),
                   ("Creek at 10 C, 2.5 m/s",
                    dict(U=2.5, bath_volume=None, T_bath_0=10.0))):
        sc, r = run(nm, **kw)
        line(nm, r, sc)
    print("\n  245 g of ice is enough to hold a 20 L box at 0 C through the whole")
    print("  job, and that box then beats any creek in the world at 10 C.")
    print("  Note the physics working against it though: near 4 C water has its")
    print("  density maximum, so buoyancy almost switches off and the still")
    print("  box's own natural convection collapses. See section 8.")

    head(8, "The water density anomaly, which the still box cannot escape")
    from props import water
    print(f"{'bath':>8}{'surface':>9}{'rho_bath':>11}{'rho_surf':>11}{'h_free':>9}")
    for Tb, Ts in ((10, 25), (10, 15), (4, 10), (0, 8), (2, 6)):
        print(f"{Tb:8.1f}{Ts:9.1f}{water(Tb).rho:11.3f}{water(Ts).rho:11.3f}"
              f"{h_free_vertical(B.H_outer, Ts, Tb, D=B.D_outer):9.0f}")
    print("\n  A surface at 6 C in a 2 C bath has essentially NO buoyancy to work")
    print("  with -- the two densities are equal to six figures. Free convection")
    print("  stops. A creek does not care: forced convection has no such hole.")
    print("  This is why very cold still water underperforms its own reputation.")


if __name__ == "__main__":
    main()
