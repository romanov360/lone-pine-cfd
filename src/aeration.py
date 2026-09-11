"""
What aeration actually does to the heat transfer.

"Aerated creek water" bundles together three physically distinct things, and
they do not push in the same direction:

  1. DISSOLVED oxygen. A true solute at ~11 mg/L, i.e. a mass fraction of
     about 1e-5. Its effect on rho, cp, k and mu is in the sixth decimal
     place. It is thermally irrelevant, and props.py measures exactly how
     irrelevant.

  2. ENTRAINED gas -- real bubbles occupying real volume. This changes the
     mixture properties in proportion to the void fraction, and it does so in
     the WRONG direction: bubbles lower rho*cp and lower k, both of which
     reduce the film coefficient. A 5% void fraction costs about 5% of the
     volumetric heat capacity and 7% of the conductivity.

  3. BUBBLE-INDUCED AGITATION. Rising bubbles drag liquid with them. This is
     the effect that actually helps, and in a still tank it is the dominant
     one: bubbling air through a tank is a well-known way to raise heat
     transfer precisely because the induced liquid circulation replaces a
     0.05 m/s buoyant plume with a 0.2-0.4 m/s recirculation.

So the honest statement is that aeration helps by stirring and hurts through
its properties, and which wins depends entirely on whether the water was
already moving. This module quantifies both sides.
"""

from __future__ import annotations

import numpy as np

from correlations import h_external, h_free_vertical
from props import bubbly_mixture, film, water

G = 9.80665


def induced_liquid_velocity(j_g: float, L: float, C: float = 1.4) -> float:
    """
    Liquid circulation velocity driven by a rising bubble swarm.

    Standard plume scaling for gas-agitated vessels (Sahai & Guthrie, 1982):
    the buoyant gas delivers power per unit mass proportional to g*j_g, and
    dimensional analysis over a plume of height L gives

        U_l ~ C * (g * j_g * L)^(1/3)

    with C around 1.2-1.6 fitted across ladle-metallurgy and bubble-column
    data. `j_g` is the superficial gas velocity, m^3 of gas per m^2 per s.
    """
    if j_g <= 0.0:
        return 0.0
    return C * (G * j_g * L) ** (1.0 / 3.0)


def void_to_superficial(void_fraction: float, u_slip: float = 0.23) -> float:
    """
    Convert a void fraction to a superficial gas velocity.

    In a bubbly swarm the gas rises relative to the liquid at roughly its
    terminal slip velocity; 0.23 m/s is the well-known plateau for 2-6 mm
    air bubbles in water, where the rise velocity is almost independent of
    size.
    """
    return void_fraction * u_slip


def h_aerated(D: float, L: float, U_stream: float, void_fraction: float,
              T_s: float, T_inf: float, include_agitation: bool = True,
              include_properties: bool = True):
    """
    Film coefficient in aerated water, separating the two competing effects.

    Returns a dict so the property penalty and the agitation benefit can be
    read off independently rather than only as a net number.
    """
    f_clean = film(T_s, T_inf)
    f_bubbly = (bubbly_mixture(0.5 * (T_s + T_inf), void_fraction)
                if include_properties and void_fraction > 0 else f_clean)

    h_base = h_external(D, L, U_stream, T_s, T_inf, fl=f_clean)
    h_prop = h_external(D, L, U_stream, T_s, T_inf, fl=f_bubbly)

    U_eff = U_stream
    if include_agitation and void_fraction > 0:
        j_g = void_to_superficial(void_fraction)
        U_ind = induced_liquid_velocity(j_g, L)
        # The induced circulation and the mean stream add in quadrature: they
        # are uncorrelated contributions to the near-wall velocity scale.
        U_eff = float(np.hypot(U_stream, U_ind))
    h_full = h_external(D, L, U_eff, T_s, T_inf, fl=f_bubbly)

    return {"void": void_fraction,
            "h_clean": h_base,
            "h_property_only": h_prop,
            "h_full": h_full,
            "U_effective": U_eff,
            "property_effect_pct": (h_prop / h_base - 1) * 100 if h_base else 0.0,
            "agitation_effect_pct": (h_full / h_prop - 1) * 100 if h_prop else 0.0,
            "net_effect_pct": (h_full / h_base - 1) * 100 if h_base else 0.0}


def dissolved_o2_effect(D: float, L: float, U: float, T_s: float, T_inf: float,
                        saturations=(0.0, 0.5, 1.0, 1.5, 2.0)):
    """Film coefficient as a function of dissolved-oxygen saturation."""
    from props import water_with_dissolved_o2
    out = []
    Tf = 0.5 * (T_s + T_inf)
    for sat in saturations:
        f = water_with_dissolved_o2(Tf, saturation=sat)
        out.append({"saturation": sat,
                    "h": h_external(D, L, U, T_s, T_inf, fl=f)})
    h0 = out[0]["h"]
    for o in out:
        o["change_pct"] = (o["h"] / h0 - 1) * 100
    return out


if __name__ == "__main__":
    from geometry import BOTTLE_500 as B

    D, L = B.D_outer, B.H_outer
    Ts, Tinf = 22.0, 10.0

    print("A. DISSOLVED oxygen: effect on the film coefficient")
    print(f"   {'saturation':>11} {'DO mg/L':>9} {'h W/m2K':>10} {'change':>10}")
    from props import do_saturation
    for o in dissolved_o2_effect(D, L, 0.35, Ts, Tinf):
        print(f"   {o['saturation']*100:10.0f}% "
              f"{o['saturation']*do_saturation(0.5*(Ts+Tinf)):9.2f} "
              f"{o['h']:10.1f} {o['change_pct']:+9.5f}%")

    print("\nB. ENTRAINED gas: properties vs agitation")
    for U in (0.0, 0.05, 0.35, 1.0):
        label = "still box" if U == 0 else f"stream {U} m/s"
        print(f"\n   {label}")
        print(f"   {'void':>6} {'U_eff':>7} {'h_clean':>9} {'h_prop':>9} "
              f"{'h_full':>9} {'prop':>8} {'agit':>8} {'net':>8}")
        for a in (0.0, 0.002, 0.01, 0.05, 0.15, 0.30):
            r = h_aerated(D, L, U, a, Ts, Tinf)
            print(f"   {a*100:5.1f}% {r['U_effective']:7.3f} "
                  f"{r['h_clean']:9.1f} {r['h_property_only']:9.1f} "
                  f"{r['h_full']:9.1f} {r['property_effect_pct']:+7.2f}% "
                  f"{r['agitation_effect_pct']:+7.1f}% {r['net_effect_pct']:+7.1f}%")
