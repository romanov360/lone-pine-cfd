"""
Convection coefficients, built on the `ht` correlation library.

Everything here returns a film coefficient h in W/m^2-K. Three families:

  * free      -- buoyancy only (the box)
  * forced    -- imposed flow only (a fast creek)
  * mixed     -- both, blended by Churchill's asymptotic rule

The mixed form matters more than it looks. A "creek" at U -> 0 must reduce
continuously to the still-water case, otherwise any comparison between the two
setups is an artifact of using two unrelated correlations.
"""

from __future__ import annotations

import math

import ht

from props import Fluid, Gr_density, Gr_from_densities, film


# ---------------------------------------------------------------------------
# Free convection
# ---------------------------------------------------------------------------


def _conduction_floor(k: float, L: float, D: float | None) -> float:
    """
    Lower bound on the film coefficient when buoyancy vanishes.

    A body suspended in a stagnant infinite medium still loses heat by pure
    conduction. The classical limit is Nu_D -> 2 for a sphere, so the floor is
    taken as 2*k/D_eq with D_eq the diameter of a sphere of the same surface
    area as the body. Without this floor the correlations return zero right at
    the water density maximum, which would say a bottle in 4 C water never
    cools at all -- wrong, and numerically fatal.
    """
    if D is None:
        return 2.0 * k / L
    A = math.pi * D * L + 2.0 * math.pi * (D / 2.0) ** 2
    D_eq = math.sqrt(A / math.pi)
    return 2.0 * k / D_eq


def h_free_vertical(L: float, T_s: float, T_inf: float,
                    D: float | None = None, fl: Fluid | None = None,
                    rho_pair: tuple[float, float] | None = None) -> float:
    """
    Free convection from a vertical cylinder of height L (and diameter D).

    Churchill-Chu for the vertical plate, with the Popiel-Churchill curvature
    correction when D is given -- a 70 mm bottle is slender enough for the
    boundary layer thickness to be a noticeable fraction of the radius.
    """
    f = fl or film(T_s, T_inf)
    floor = _conduction_floor(f.k, L, D)
    if abs(T_s - T_inf) < 1e-9:
        return floor
    # `rho_pair` is (rho at the bulk temperature, rho at the surface
    # temperature) for whatever the bath actually is. Without it the bath is
    # assumed to be water.
    Gr = (Gr_from_densities(rho_pair[0], rho_pair[1], f, L)
          if rho_pair is not None else Gr_density(T_s, T_inf, L, fl=f))
    if Gr <= 0.0 or not math.isfinite(Gr):
        return floor
    if D is None:
        Nu = ht.Nu_vertical_plate_Churchill(f.Pr, Gr)
    else:
        Nu = ht.Nu_vertical_cylinder(Pr=f.Pr, Gr=Gr, L=L, D=D,
                                     Method="Popiel & Churchill")
    return max(Nu * f.k / L, floor)


def h_free_horizontal(D: float, T_s: float, T_inf: float,
                      fl: Fluid | None = None,
                      rho_pair: tuple[float, float] | None = None) -> float:
    """Free convection from a horizontal cylinder (bottle lying down)."""
    f = fl or film(T_s, T_inf)
    floor = _conduction_floor(f.k, D, D)
    if abs(T_s - T_inf) < 1e-9:
        return floor
    Gr = (Gr_from_densities(rho_pair[0], rho_pair[1], f, D)
          if rho_pair is not None else Gr_density(T_s, T_inf, D, fl=f))
    if Gr <= 0.0 or not math.isfinite(Gr):
        return floor
    Nu = ht.Nu_horizontal_cylinder_Churchill_Chu(f.Pr, Gr)
    return max(Nu * f.k / D, floor)


# ---------------------------------------------------------------------------
# Forced convection
# ---------------------------------------------------------------------------

_FORCED = {
    "Churchill-Bernstein": lambda Re, Pr, f: ht.Nu_cylinder_Churchill_Bernstein(Re, Pr),
    "Zukauskas":           lambda Re, Pr, f: ht.Nu_cylinder_Zukauskas(Re, Pr, Prw=Pr),
    "Whitaker":            lambda Re, Pr, f: ht.Nu_cylinder_Whitaker(Re, Pr, mu=f.mu, muw=f.mu),
    "Sanitjai-Goldstein":  lambda Re, Pr, f: ht.Nu_cylinder_Sanitjai_Goldstein(Re, Pr),
    "McAdams":             lambda Re, Pr, f: ht.Nu_cylinder_McAdams(Re, Pr),
    "Fand":                lambda Re, Pr, f: ht.Nu_cylinder_Fand(Re, Pr),
}


def h_forced_crossflow(D: float, U: float, T_s: float, T_inf: float,
                       method: str = "Churchill-Bernstein",
                       fl: Fluid | None = None) -> float:
    """Cross-flow over a cylinder of diameter D at approach velocity U."""
    if U <= 0.0:
        return 0.0
    f = fl or film(T_s, T_inf)
    Re = f.Re(D, U)
    Nu = _FORCED[method](Re, f.Pr, f)
    return Nu * f.k / D


def forced_spread(D: float, U: float, T_s: float, T_inf: float,
                  fl: Fluid | None = None) -> dict[str, float]:
    """All six correlations at once, to size the correlation uncertainty."""
    return {m: h_forced_crossflow(D, U, T_s, T_inf, m, fl) for m in _FORCED}


# ---------------------------------------------------------------------------
# Mixed convection
# ---------------------------------------------------------------------------


def h_mixed(h_forced: float, h_free: float, n: float = 4.0,
            assisting: bool = True) -> float:
    """
    Churchill's asymptotic blend of a forced and a free limit.

    n = 4 is Churchill's recommendation for cross-flow over a cylinder with
    the buoyant plume transverse to the oncoming stream, which is exactly the
    geometry of an upright bottle standing in a creek. `assisting=False`
    (opposing flow) subtracts, and is the one case where adding a slow current
    can *reduce* the transfer coefficient below the still-water value.
    """
    if assisting:
        return (h_forced**n + h_free**n) ** (1.0 / n)
    val = abs(h_forced**n - h_free**n)
    return val ** (1.0 / n)


def h_external(D: float, L: float, U: float, T_s: float, T_inf: float,
               orientation: str = "vertical",
               method: str = "Churchill-Bernstein",
               fl: Fluid | None = None,
               n: float = 4.0,
               rho_pair: tuple[float, float] | None = None) -> float:
    """
    The external film coefficient for either setup.

    U = 0 gives the box; U > 0 gives the creek. One function, one continuous
    curve, so the two setups are never compared across a modelling seam.
    """
    f = fl or film(T_s, T_inf)
    if orientation == "vertical":
        hf = h_free_vertical(L, T_s, T_inf, D=D, fl=f, rho_pair=rho_pair)
    else:
        hf = h_free_horizontal(D, T_s, T_inf, fl=f, rho_pair=rho_pair)
    hc = h_forced_crossflow(D, U, T_s, T_inf, method=method, fl=f)
    return h_mixed(hc, hf, n=n)


# ---------------------------------------------------------------------------
# Internal convection (inside the bottle)
# ---------------------------------------------------------------------------


def h_internal(D_i: float, H_i: float, T_bulk: float, T_wall: float,
               fl: Fluid | None = None) -> float:
    """
    Free convection inside the bottle, from the bulk liquid to the wall.

    Standard treatment for a cooling liquid-filled vessel: the downward
    boundary layer on the cold inner wall is a vertical-plate free-convection
    layer driven by the bulk-to-wall difference. Uses the same Churchill-Chu
    form with curvature correction as the outside.
    """
    return h_free_vertical(H_i, T_wall, T_bulk, D=D_i, fl=fl)


def U_overall(bottle, h_out: float, h_in: float) -> float:
    """
    Overall bottle-to-environment coefficient referred to the OUTSIDE area.

    Series chain: interior film -> shell conduction -> exterior film.
    """
    A_o, A_i = bottle.A_outer, bottle.A_inner
    R = 1.0 / (h_in * A_i) + bottle.R_wall + 1.0 / (h_out * A_o)
    return 1.0 / (R * A_o)


def resistance_split(bottle, h_out: float, h_in: float) -> dict[str, float]:
    """Where the resistance actually sits, as a fraction of the total."""
    R_in = 1.0 / (h_in * bottle.A_inner)
    R_w = bottle.R_wall
    R_out = 1.0 / (h_out * bottle.A_outer)
    tot = R_in + R_w + R_out
    return {"inside": R_in / tot, "wall": R_w / tot, "outside": R_out / tot,
            "R_total": tot}


if __name__ == "__main__":
    from geometry import BOTTLE_500 as B

    Ts, Tinf = 25.0, 10.0
    print(f"Bottle D={B.D_outer*1000:.0f} mm  H={B.H_outer*1000:.0f} mm, "
          f"surface {Ts} C in {Tinf} C water\n")

    hfv = h_free_vertical(B.H_outer, Ts, Tinf, D=B.D_outer)
    hfh = h_free_horizontal(B.D_outer, Ts, Tinf)
    print(f"free convection, upright   h = {hfv:7.1f} W/m^2-K")
    print(f"free convection, on its side h = {hfh:7.1f} W/m^2-K\n")

    print("forced cross-flow, correlation spread:")
    print(f"{'U m/s':>7} {'Re':>9} " + " ".join(f"{m[:9]:>10}" for m in _FORCED))
    f = film(Ts, Tinf)
    for U in (0.02, 0.05, 0.1, 0.2, 0.35, 0.6, 1.0, 2.0):
        sp = forced_spread(B.D_outer, U, Ts, Tinf, f)
        print(f"{U:7.2f} {f.Re(B.D_outer,U):9.0f} "
              + " ".join(f"{sp[m]:10.0f}" for m in _FORCED))

    print("\nmixed (Churchill n=4) vs pure-forced, upright bottle:")
    print(f"{'U m/s':>7} {'h_forced':>9} {'h_free':>8} {'h_mixed':>9} {'Ri':>8}")
    for U in (0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.6, 1.0, 2.0):
        hc = h_forced_crossflow(B.D_outer, U, Ts, Tinf, fl=f)
        hm = h_mixed(hc, hfv)
        Ri = f.Gr(B.D_outer, Ts - Tinf) / max(f.Re(B.D_outer, U), 1e-9) ** 2
        print(f"{U:7.2f} {hc:9.1f} {hfv:8.1f} {hm:9.1f} {Ri:8.2f}")

    print("\nresistance chain (upright, creek U=0.35 m/s):")
    hi = h_internal(B.D_inner, B.H_inner, 25.0, 24.0)
    for label, U in (("box  U=0", 0.0), ("creek U=0.35", 0.35)):
        ho = h_external(B.D_outer, B.H_outer, U, Ts, Tinf, fl=f)
        rs = resistance_split(B, ho, hi)
        Uo = U_overall(B, ho, hi)
        print(f"  {label:14s} h_out={ho:6.0f}  h_in={hi:5.0f}  U_overall={Uo:6.1f}"
              f"   inside {rs['inside']*100:4.1f}% | wall {rs['wall']*100:4.1f}%"
              f" | outside {rs['outside']*100:4.1f}%")
