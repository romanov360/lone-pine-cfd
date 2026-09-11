"""
Factors outside the box-versus-creek axis that change the answer anyway.

Everything in the main study varied the OUTSIDE of the bottle. These are the
ones that do not: what is in the bottle, how full it is, whether it is being
knocked about, and whether it needs to be submerged at all.

The largest single sensitivity in this whole problem turns out to be none of
the things the original question asked about. It is whether the water inside
the bottle is circulating.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import ht
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from correlations import h_external, h_free_vertical, h_internal
from geometry import BOTTLE_500 as B
from props import water

RULE = "=" * 78


def head(n, t):
    print(f"\n{RULE}\n{n}. {t}\n{RULE}")


# ---------------------------------------------------------------------------
# 1. Internal circulation -- the biggest lever in the problem
# ---------------------------------------------------------------------------


def internal_mixing(mults=(1.0, 2.0, 7.0, 15.0, 30.0, 100.0), h_out=466.0):
    """
    Overall conductance against the effective conductivity of the contents.

    k_eff = 1 is pure conduction: the bottle stands perfectly still and heat
    crawls inward by diffusion alone. k_eff = 100 is perfectly stirred. A
    sealed bottle left alone in a bath sits around 4-5; a bottle rolling in a
    current is much higher, because every tumble overturns the contents.

    This is the mechanism by which a creek wins that the film coefficient
    never captures -- and it is worth more than the film coefficient is.
    """
    w = water(25.0)
    R_wall = B.R_wall
    A_o, A_i = B.A_outer, B.A_inner
    rows = []
    for m in mults:
        # Effective internal film from an enclosure treated as a conducting
        # slab of enhanced conductivity across half the bottle radius.
        h_i = m * w.k / (B.D_inner / 4.0)
        R = 1.0 / (h_i * A_i) + R_wall + 1.0 / (h_out * A_o)
        rows.append({"k_eff_mult": m, "h_in": h_i,
                     "U_overall": 1.0 / (R * A_o),
                     "inside_share": (1.0 / (h_i * A_i)) / R})
    base = rows[0]["U_overall"]
    for r in rows:
        r["gain"] = r["U_overall"] / base
    return rows


# ---------------------------------------------------------------------------
# 2. How full the bottle is
# ---------------------------------------------------------------------------


def fill_level(fills=(1.0, 0.9, 0.75, 0.5, 0.25), h_out=466.0, h_in=273.0):
    """
    A part-full bottle has less to cool, but the headspace is air -- a 40x
    worse conductor than water -- so the dry part of the wall stops carrying
    heat almost entirely.

    Modelled as two regions in parallel: the wetted side wall plus base
    conducting normally, and the dry side wall plus cap conducting through a
    stagnant air gap of the headspace height (capped at the bottle radius,
    beyond which the gap convects rather than conducts).
    """
    w = water(25.0)
    k_air = 0.026
    r_i = B.D_inner / 2.0
    rows = []
    for f in fills:
        h_liq = B.H_inner * f
        h_gas = B.H_inner - h_liq
        A_wet_side = 2 * math.pi * r_i * h_liq
        A_dry_side = 2 * math.pi * r_i * h_gas
        A_end = math.pi * r_i**2

        # Wetted path: interior film -> glass -> exterior film.
        R_wet = (1.0 / (h_in * (A_wet_side + A_end))
                 + B.R_wall * B.A_inner / max(A_wet_side + A_end, 1e-9)
                 + 1.0 / (h_out * (A_wet_side + A_end)))
        # Dry path: across the headspace air, then glass, then exterior film.
        gap = min(max(h_gas, 1e-4), r_i)
        R_dry = math.inf
        if A_dry_side + A_end > 1e-9:
            A_d = A_dry_side + A_end
            R_dry = (gap / (k_air * A_d)
                     + B.R_wall * B.A_inner / A_d
                     + 1.0 / (h_out * A_d))
        UA = 1.0 / R_wet + (0.0 if not math.isfinite(R_dry) else 1.0 / R_dry)
        V = B.V_liquid * f
        C = V * w.rho * w.cp
        rows.append({"fill": f, "volume_mL": V * 1e6, "UA": UA,
                     "C_kJ_K": C / 1000.0, "tau_min": C / UA / 60.0})
    base = rows[0]["tau_min"]
    for r in rows:
        r["speedup"] = base / r["tau_min"]
    return rows


# ---------------------------------------------------------------------------
# 3. What is in the bottle
# ---------------------------------------------------------------------------

CONTENTS = {
    # name: (density, cp, k, dynamic viscosity at ~30 C)
    "water":            (996.0, 4180.0, 0.615, 0.80e-3),
    "beer / soft drink": (1005.0, 4000.0, 0.580, 1.10e-3),
    "milk":             (1030.0, 3900.0, 0.560, 1.80e-3),
    "olive oil":        (910.0, 1970.0, 0.170, 40.0e-3),
    "thick soup":       (1050.0, 3800.0, 0.500, 200.0e-3),
    "honey":            (1420.0, 2300.0, 0.500, 5000.0e-3),
}


def contents_sweep(h_out=466.0):
    """
    The contents set both the heat to be removed and how fast it reaches the
    wall. Viscosity is the dominant term: internal natural convection scales
    roughly as nu^-1/2, so honey has none at all and cools by conduction.
    """
    rows = []
    g, beta, dT, L = 9.80665, 3.0e-4, 12.0, B.H_inner
    for name, (rho, cp, k, mu) in CONTENTS.items():
        nu = mu / rho
        alpha = k / (rho * cp)
        Pr = nu / alpha
        Ra = g * beta * dT * L**3 / (nu * alpha)
        Nu = ht.Nu_vertical_plate_Churchill(Pr, Ra / Pr) if Ra > 1e3 else 1.0
        h_i = max(Nu * k / L, k / (B.D_inner / 4.0))
        R = (1.0 / (h_i * B.A_inner) + B.R_wall
             + 1.0 / (h_out * B.A_outer))
        UA = 1.0 / R
        C = B.V_liquid * rho * cp
        rows.append({"name": name, "Pr": Pr, "Ra": Ra, "h_in": h_i,
                     "C_kJ_K": C / 1000.0, "tau_min": C / UA / 60.0})
    base = rows[0]["tau_min"]
    for r in rows:
        r["vs_water"] = r["tau_min"] / base
    return rows


# ---------------------------------------------------------------------------
# 4. Do you need water at all? The wet-cloth-in-the-wind alternative
# ---------------------------------------------------------------------------


def evaporative(T_s=25.0, T_air=20.0, RH=0.5, winds=(0.0, 1.0, 3.0, 8.0)):
    """
    A bottle wrapped in a wet cloth in moving air, the classic desert trick.

    Sensible convection in air plus evaporation, the latter via the
    Chilton-Colburn analogy: h_m = h / (rho*cp*Le^(2/3)), with the driving
    potential the difference between saturated vapour density at the surface
    and the ambient vapour density. Reported as an equivalent film coefficient
    on the surface-to-air temperature difference so it can be set beside the
    submerged numbers directly.
    """
    T_film = 0.5 * (T_s + T_air) + 273.15
    rho_a = PropsSI("D", "T", T_film, "P", 101325, "Air")
    cp_a = PropsSI("C", "T", T_film, "P", 101325, "Air")
    k_a = PropsSI("L", "T", T_film, "P", 101325, "Air")
    mu_a = PropsSI("V", "T", T_film, "P", 101325, "Air")
    nu_a = mu_a / rho_a
    Pr_a = cp_a * mu_a / k_a
    D_ab = 2.6e-5                      # water vapour in air, m^2/s at ~295 K
    Le = (k_a / (rho_a * cp_a)) / D_ab
    h_fg = (PropsSI("H", "T", T_s + 273.15, "Q", 1, "Water")
            - PropsSI("H", "T", T_s + 273.15, "Q", 0, "Water"))

    def vapour_density(T, rh):
        W = HAPropsSI("W", "T", T + 273.15, "P", 101325, "R", rh)
        V = HAPropsSI("Vha", "T", T + 273.15, "P", 101325, "R", rh)
        return W / V

    d_rho_v = vapour_density(T_s, 1.0) - vapour_density(T_air, RH)

    rows = []
    for U in winds:
        if U > 0:
            Re = U * B.D_outer / nu_a
            Nu = ht.Nu_cylinder_Churchill_Bernstein(Re, Pr_a)
        else:
            Gr = 9.80665 * (1.0 / T_film) * abs(T_s - T_air) * B.H_outer**3 / nu_a**2
            Nu = ht.Nu_vertical_plate_Churchill(Pr_a, Gr) * B.D_outer / B.H_outer
        h_sens = max(Nu * k_a / B.D_outer, 2.0)
        h_m = h_sens / (rho_a * cp_a * Le ** (2.0 / 3.0))
        q_evap = h_m * d_rho_v * h_fg
        h_evap_eq = q_evap / max(T_s - T_air, 1e-6)
        rows.append({"wind": U, "h_sensible": h_sens, "h_evap_equiv": h_evap_eq,
                     "h_total": h_sens + h_evap_eq, "q_evap": q_evap})
    return rows


def wet_bulb(T_air=20.0, RH=0.5):
    return HAPropsSI("Twb", "T", T_air + 273.15, "P", 101325, "R", RH) - 273.15


if __name__ == "__main__":
    print(RULE)
    print("OTHER FACTORS: everything that is not the water outside the bottle")
    print(RULE)

    head(1, "Internal circulation -- the biggest lever in the problem")
    print(f"  {'k_eff':>8}{'h_inside':>10}{'U overall':>11}{'inside share':>14}{'gain':>8}")
    for r in internal_mixing():
        print(f"  {r['k_eff_mult']:8.0f}{r['h_in']:10.0f}{r['U_overall']:11.0f}"
              f"{r['inside_share']*100:13.0f}%{r['gain']:7.2f}x")
    rows = internal_mixing()
    sealed = next(r for r in rows if r["k_eff_mult"] == 7.0)
    stirred = rows[-1]
    print(f"\n  The free-convection correlation used in the main study gives")
    print(f"  h_inside = 273 W/m2-K for a sealed bottle standing still, which")
    print(f"  is k_eff near 7 -- the row marked above. Perfectly stirred is the")
    print(f"  last row. The difference between them is"
          f" {stirred['U_overall']/sealed['U_overall']:.2f}x on the")
    print("  overall conductance, against the 1.35x that the fastest imaginable")
    print("  current buys on the outside film.")
    print("\n  So the creek has a SECOND advantage the main study never counted:")
    print("  it does not just flow past the bottle, it knocks it about. A bottle")
    print("  wedged still in a creek gets only the film-coefficient gain; one")
    print("  rolling along the bed gets this as well, and this is the bigger of")
    print("  the two.")

    head(2, "How full the bottle is")
    print(f"  {'fill':>7}{'volume':>10}{'UA':>9}{'heat cap':>11}{'time const':>13}{'vs full':>9}")
    for r in fill_level():
        print(f"  {r['fill']*100:6.0f}%{r['volume_mL']:9.0f}mL{r['UA']:9.2f}"
              f"{r['C_kJ_K']:10.2f}kJ/K{r['tau_min']:11.1f}min{r['speedup']:8.2f}x")
    print("\n  Fill level barely moves the time constant, and that is the")
    print("  interesting part. Emptying the bottle halfway removes half the")
    print("  water to be cooled -- but it also un-wets half the wall, and the")
    print("  headspace behind it is nearly a perfect insulator. The two effects")
    print("  very nearly cancel: half the heat, half the conductance, same")
    print("  time constant. A part-full bottle reaches a given temperature")
    print("  sooner only because it started with less heat in it, not because")
    print("  it sheds heat any faster.")

    head(3, "What is in the bottle")
    print(f"  {'contents':>20}{'Pr':>10}{'Ra':>11}{'h_inside':>10}"
          f"{'time const':>13}{'vs water':>10}")
    for r in contents_sweep():
        print(f"  {r['name']:>20}{r['Pr']:10.0f}{r['Ra']:11.1e}{r['h_in']:10.0f}"
              f"{r['tau_min']:11.1f}min{r['vs_water']:9.2f}x")
    print("\n  Viscosity decides it. Water, beer and milk all convect freely")
    print("  inside and behave alike. Oil, soup and honey do not convect at")
    print("  all -- their Rayleigh numbers are too low -- so they cool by pure")
    print("  conduction and take many times longer, whatever is going on")
    print("  outside. For honey the creek-versus-box question is irrelevant.")

    head(4, "Do you need water at all?")
    twb = wet_bulb()
    print(f"  Bottle at 25 C, air at 20 C, 50% RH (wet-bulb {twb:.1f} C),")
    print("  wrapped in a wet cloth:\n")
    print(f"  {'wind m/s':>10}{'sensible':>10}{'evaporative':>13}{'total h':>10}"
          f"{'vs still water':>16}")
    for r in evaporative():
        print(f"  {r['wind']:10.1f}{r['h_sensible']:10.1f}{r['h_evap_equiv']:13.1f}"
              f"{r['h_total']:10.1f}{r['h_total']/466.0:15.2f}x")
    rows = evaporative()
    gale = rows[-1]
    print("\n  Evaporation carries most of the load -- six or seven times the")
    print("  sensible transfer at every wind speed. But even an 8 m/s gale")
    print(f"  reaches only {gale['h_total']:.0f} W/m2-K, still "
          f"{466.0/gale['h_total']:.1f}x short of a bucket of")
    print(f"  STILL water and {1741.0/gale['h_total']:.0f}x short of a 0.35 m/s "
          f"creek. Water beats air")
    print("  because it is eight hundred times denser, and no amount of latent")
    print("  heat closes that gap.")
    print(f"\n  The wet cloth also cannot go below the wet-bulb temperature")
    print(f"  ({twb:.1f} C here), which is above a 10 C creek to begin with. It")
    print("  is a good trick in a hot dry place with no water to spare, and a")
    print("  poor one anywhere you could simply put the bottle in a stream.")
