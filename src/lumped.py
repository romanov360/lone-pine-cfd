"""
Transient cooling models at the 0-D / 1-node level.

Two bodies: the bottle (water core + glass shell) and the bath. The bath is
either a finite box that warms up as it absorbs heat, or a creek whose
temperature never moves. That single difference -- finite vs infinite heat
capacity -- turns out to dominate the whole comparison, so it is worth having
a clean model of it before any CFD is run.

Properties are re-evaluated at every step, and the film coefficients are
recomputed from the instantaneous temperature difference, so the exponential
decay here is not a constant-coefficient approximation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from correlations import (
    U_overall,
    h_external,
    h_internal,
    resistance_split,
)
from geometry import BOTTLE_500, Bottle
from props import brine_nacl, brine_rho, bubbly_mixture, film, water


@dataclass
class Scenario:
    """One cooling experiment."""

    name: str
    bottle: Bottle = field(default_factory=lambda: BOTTLE_500)
    T_bottle_0: float = 40.0        # C
    T_bath_0: float = 10.0          # C
    U: float = 0.0                  # m/s, 0 = still box
    bath_volume: float | None = None  # m^3; None = infinite (creek)
    orientation: str = "vertical"
    void_fraction: float = 0.0      # entrained air in the bath
    bath_kind: str = "water"        # 'water' | 'brine'
    brine_w: float = 0.23           # NaCl mass fraction when bath_kind='brine'
    shelter: float = 1.0            # fraction of the nominal current the
                                    # bottle actually sees; < 1 for a bottle
                                    # tucked behind a rock or lying in the
                                    # bed boundary layer
    engaged_fraction: float = 1.0   # fraction of the box that actually takes
                                    # part on the timescale of interest; < 1
                                    # represents thermal stratification
    T_room: float | None = None     # if set, the box also exchanges with the
                                    # room through its walls
    UA_room: float = 0.0            # W/K for that exchange
    forced_method: str = "Churchill-Bernstein"
    h_out_scale: float = 1.0        # multiplier for sensitivity studies
    h_in_scale: float = 1.0

    @property
    def infinite_bath(self) -> bool:
        return self.bath_volume is None

    def C_bath(self) -> float:
        if self.infinite_bath:
            return np.inf
        f = (brine_nacl(self.T_bath_0, self.brine_w)
             if self.bath_kind == "brine" else water(self.T_bath_0))
        return self.engaged_fraction * self.bath_volume * f.rho * f.cp

    def T_equilibrium(self) -> float:
        """Where the two bodies end up once nothing else can happen."""
        if self.infinite_bath:
            return self.T_bath_0
        Cb = self.bottle.C_total(self.T_bottle_0)
        Cw = self.C_bath()
        return (Cb * self.T_bottle_0 + Cw * self.T_bath_0) / (Cb + Cw)


def _bath_fluid(sc: Scenario, T_s: float, T_bath: float):
    """Film properties of the bath, including any entrained air."""
    if sc.bath_kind == "brine":
        return brine_nacl(0.5 * (T_s + T_bath), sc.brine_w)
    if sc.void_fraction > 0.0:
        return bubbly_mixture(0.5 * (T_s + T_bath), sc.void_fraction)
    return film(T_s, T_bath)


def _rho_pair(sc: Scenario, T_s: float, T_bath: float):
    """
    Densities at the bulk and surface temperatures, for the buoyancy term.

    Only needed when the bath is not water; the correlations assume water
    otherwise.
    """
    if sc.bath_kind == "brine":
        return (brine_rho(T_bath, sc.brine_w), brine_rho(T_s, sc.brine_w))
    return None


def coefficients(sc: Scenario, T_bottle: float, T_bath: float):
    """
    Solve the series chain  interior film -> shell -> exterior film  for the
    outer wall temperature, and return the resulting coefficients.

    Both films depend on their own temperature drop, so the wall temperature
    is the root of a heat-balance residual rather than something that can be
    written down. It is found with a bracketed Brent solve on

        residual(T_s) = q_in(T_s) - q_out(T_s)

    which is guaranteed to bracket a root: at T_s = T_bath nothing leaves the
    outside but heat still flows in from the core, and at T_s = T_bottle the
    reverse. A damped fixed-point iteration was tried first and is not good
    enough -- near the water density maximum h_out is so steep a function of
    T_s that the iterate never settles to a repeatable value, which leaves the
    ODE integrator unable to estimate its own error and it stalls.
    """
    B = sc.bottle
    A_o, A_i, R_w = B.A_outer, B.A_inner, B.R_wall

    def films(T_s):
        f_bath = _bath_fluid(sc, T_s, T_bath)
        h_o = sc.h_out_scale * h_external(
            B.D_outer, B.H_outer, sc.U * sc.shelter, T_s, T_bath,
            orientation=sc.orientation, method=sc.forced_method, fl=f_bath,
            rho_pair=_rho_pair(sc, T_s, T_bath))
        q_out = h_o * A_o * (T_s - T_bath)
        T_wi = T_s + q_out * R_w
        h_i = sc.h_in_scale * h_internal(B.D_inner, B.H_inner, T_bottle, T_wi)
        q_in = h_i * A_i * (T_bottle - T_wi)
        return h_o, h_i, q_out, q_in, T_wi

    if abs(T_bottle - T_bath) < 1e-9:
        h_o, h_i, *_ = films(T_bath)
        return {"h_out": h_o, "h_in": h_i, "U": U_overall(B, h_o, h_i),
                "T_surface": T_bath,
                "split": resistance_split(B, h_o, h_i)}

    lo, hi = min(T_bath, T_bottle), max(T_bath, T_bottle)
    eps = 1e-7 * max(1.0, hi - lo)

    def residual(T_s):
        _, _, q_out, q_in, _ = films(T_s)
        return q_in - q_out

    try:
        T_s = brentq(residual, lo + eps, hi - eps, xtol=1e-9, rtol=1e-12)
    except ValueError:
        T_s = 0.5 * (lo + hi)

    h_o, h_i, _, _, _ = films(T_s)
    return {"h_out": h_o, "h_in": h_i, "U": U_overall(B, h_o, h_i),
            "T_surface": T_s,
            "split": resistance_split(B, h_o, h_i)}


def simulate(sc: Scenario, t_end: float = 14400.0, n_out: int = 2000,
             stop_at: float | None = None):
    """
    Integrate the two-body system. Returns a dict of arrays.

    `stop_at` terminates the run when the bottle reaches that temperature.
    Its purpose is sub-zero baths: this model has no phase change, so left to
    itself it would happily carry 500 mL of water down to -16 C as a liquid.
    Stopping at 0 C and reporting the freeze separately (see freeze_time) is
    both faster and honest about what the model actually knows.
    """
    B = sc.bottle
    A_o = B.A_outer
    C_bath = sc.C_bath()

    def rhs(t, y):
        Tb, Tw = y
        C_b = B.C_total(Tb)
        c = coefficients(sc, Tb, Tw)
        q = c["U"] * A_o * (Tb - Tw)           # W, bottle -> bath
        dTb = -q / C_b
        if np.isinf(C_bath):
            dTw = 0.0
        else:
            q_room = 0.0
            if sc.T_room is not None and sc.UA_room > 0.0:
                q_room = sc.UA_room * (Tw - sc.T_room)
            dTw = (q - q_room) / C_bath
        return [dTb, dTw]

    events = None
    if stop_at is not None:
        def hit_floor(t, y):
            return y[0] - stop_at
        hit_floor.terminal = True
        hit_floor.direction = -1
        events = [hit_floor]

    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(rhs, (0.0, t_end), [sc.T_bottle_0, sc.T_bath_0],
                    t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10,
                    max_step=t_end / 50.0, events=events)

    T_bottle, T_bath = sol.y[0], sol.y[1]
    if stop_at is not None and sol.t_events and len(sol.t_events[0]):
        keep = sol.t <= sol.t_events[0][0]
        sol.t, T_bottle, T_bath = sol.t[keep], T_bottle[keep], T_bath[keep]
    coef = [coefficients(sc, a, b) for a, b in
            zip(T_bottle[::max(1, n_out // 200)], T_bath[::max(1, n_out // 200)])]
    return {
        "t": sol.t,
        "T_bottle": T_bottle,
        "T_bath": T_bath,
        "scenario": sc,
        "h_out": np.array([c["h_out"] for c in coef]),
        "h_in": np.array([c["h_in"] for c in coef]),
        "U": np.array([c["U"] for c in coef]),
        "t_coef": sol.t[::max(1, n_out // 200)],
    }


def time_to(res, T_target: float) -> float:
    """Time in seconds for the bottle to reach T_target, or inf if never."""
    t, T = res["t"], res["T_bottle"]
    if T[-1] > T_target:
        return np.inf
    i = int(np.argmax(T <= T_target))
    if i == 0:
        return 0.0
    t0, t1, T0, T1 = t[i - 1], t[i], T[i - 1], T[i]
    return t0 + (T0 - T_target) * (t1 - t0) / (T0 - T1)


def half_life(res) -> float:
    """Time to close half the initial gap to the *bath's initial* temperature."""
    sc = res["scenario"]
    return time_to(res, 0.5 * (sc.T_bottle_0 + sc.T_bath_0))


def freeze_time(sc: Scenario, T_start: float = 0.0) -> tuple[float, float]:
    """
    Time to cool the contents to 0 C, and the extra time to then freeze them
    solid, seconds.

    The second number is the latent-heat plateau: 500 mL of water gives up
    167 kJ turning to ice, against roughly 2.1 kJ per kelvin of sensible heat,
    so freezing through takes far longer than the whole cooling run that
    preceded it. Anyone using a salt-ice bath to chill a drink has a wide
    window before the bottle is in danger -- and a glass bottle full of ice is
    a burst bottle, so the window matters.
    """
    from props import ICE_LATENT, water

    r = simulate(sc, t_end=4 * 3600, stop_at=T_start + 0.25)
    t_to_zero = time_to(r, T_start + 0.5)
    if not np.isfinite(t_to_zero):
        return np.inf, np.inf

    # Heat rate at the moment freezing starts, held for the plateau.
    c = coefficients(sc, T_start, r["T_bath"][-1])
    q = c["U"] * sc.bottle.A_outer * (T_start - r["T_bath"][-1])
    if q <= 0:
        return t_to_zero, np.inf
    m = sc.bottle.V_liquid * water(T_start).rho
    return t_to_zero, m * ICE_LATENT / q


def fmt(t: float) -> str:
    if not np.isfinite(t):
        return "  never"
    if t < 90:
        return f"{t:5.1f} s"
    if t < 5400:
        return f"{t/60:5.1f} m"
    return f"{t/3600:5.2f} h"


if __name__ == "__main__":
    box20 = Scenario("box, 20 L still", U=0.0, bath_volume=0.020)
    creek = Scenario("creek, 0.35 m/s", U=0.35, bath_volume=None)

    print(f"{'scenario':22s} {'T_eq':>6} {'t->25C':>8} {'t->15C':>8} "
          f"{'t->12C':>8} {'T@30min':>8} {'T@2h':>7}")
    for sc in (box20, creek):
        r = simulate(sc, t_end=4 * 3600)
        T30 = np.interp(1800, r["t"], r["T_bottle"])
        T2h = np.interp(7200, r["t"], r["T_bottle"])
        print(f"{sc.name:22s} {sc.T_equilibrium():6.2f} "
              f"{fmt(time_to(r,25)):>8} {fmt(time_to(r,15)):>8} "
              f"{fmt(time_to(r,12)):>8} {T30:7.2f}C {T2h:6.2f}C")
