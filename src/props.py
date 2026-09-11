"""
Thermophysical property layer for the bottle-cooling study.

Three things live here:

1. Pure liquid water properties from CoolProp (IAPWS-95 backed).
2. The effect of *dissolved* oxygen on those properties (the thing people
   usually mean by "aerated creek water").
3. The effect of *entrained* gas -- actual bubbles -- which is a different
   and much larger effect than dissolved gas.

The distinction in (2) vs (3) turns out to matter a great deal, see
docs/oxygenation.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from CoolProp.CoolProp import PropsSI

P_ATM = 101325.0  # Pa

# ---------------------------------------------------------------------------
# 1. Pure water
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Fluid:
    """Frozen bundle of the properties a convection correlation needs."""

    T: float        # K
    rho: float      # kg/m^3
    cp: float       # J/kg-K
    k: float        # W/m-K
    mu: float       # Pa-s
    beta: float     # 1/K  (volumetric thermal expansion)

    @property
    def nu(self) -> float:
        """Kinematic viscosity, m^2/s."""
        return self.mu / self.rho

    @property
    def alpha(self) -> float:
        """Thermal diffusivity, m^2/s."""
        return self.k / (self.rho * self.cp)

    @property
    def Pr(self) -> float:
        return self.cp * self.mu / self.k

    def Gr(self, L: float, dT: float, g: float = 9.80665) -> float:
        """Grashof number on length L for temperature difference dT."""
        return g * abs(self.beta) * abs(dT) * L**3 / self.nu**2

    def Ra(self, L: float, dT: float, g: float = 9.80665) -> float:
        return self.Gr(L, dT, g) * self.Pr

    def Re(self, L: float, U: float) -> float:
        return U * L / self.nu


@lru_cache(maxsize=200_000)
def _water_cached(T_mK: int, P: float) -> Fluid:
    T = T_mK / 1000.0
    return Fluid(
        T=T,
        rho=PropsSI("D", "T", T, "P", P, "Water"),
        cp=PropsSI("C", "T", T, "P", P, "Water"),
        k=PropsSI("L", "T", T, "P", P, "Water"),
        mu=PropsSI("V", "T", T, "P", P, "Water"),
        beta=PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", P, "Water"),
    )


def water(T_C: float, P: float = P_ATM) -> Fluid:
    """
    Liquid water properties at temperature `T_C` in degrees Celsius.

    An ice-water bath sits exactly at 0 C, which is the melting point and so
    just outside the liquid range CoolProp will evaluate. The temperature is
    nudged 10 mK into the liquid region; no property changes measurably over
    that interval.
    """
    # Quantised to 1 mK so the lookup can be cached. Every property used here
    # varies by far less than its own correlation uncertainty over 1 mK, and
    # the transient solvers call this tens of millions of times.
    T_mK = int(round(max(T_C + 273.15, 273.16) * 1000.0))
    return _water_cached(T_mK, P)


def film(T_surface_C: float, T_bulk_C: float) -> Fluid:
    """Properties evaluated at the film temperature, the usual convention."""
    return water(0.5 * (T_surface_C + T_bulk_C))


def Gr_from_densities(rho_bulk: float, rho_surface: float, fl: "Fluid",
                      L: float, g: float = 9.80665) -> float:
    """
    Grashof number from an explicit pair of densities.

    The general form behind Gr_density, so a bath that is not pure water --
    brine, say -- can supply its own density pair.
    """
    return g * abs(rho_bulk - rho_surface) / fl.rho * L**3 / fl.nu**2


def Gr_density(T_surface_C: float, T_bulk_C: float, L: float,
               g: float = 9.80665, fl: "Fluid | None" = None) -> float:
    """
    Grashof number from the actual density difference rather than from a
    linearised beta*dT.

    This matters near an ice bath. Water has a density maximum at 3.98 C, so
    beta passes through zero there and is NEGATIVE below it. The usual
    Gr = g*beta*dT*L^3/nu^2 then returns a negative number and the free
    convection correlations produce NaN -- which is the algebra complaining
    about a real physical effect, not a bug to be clamped away.

    Using Gr = g*|rho_bulk - rho_surface|/rho_film * L^3/nu^2 is exact for any
    temperature pair, handles the anomaly correctly, and collapses to the
    linearised form away from 4 C. Near 4 C it correctly predicts that
    buoyancy almost switches off: a surface at 6 C in a 0 C bath drives far
    weaker convection than a naive beta*dT would suggest, because the two
    temperatures sit on opposite sides of the density maximum and their
    densities are nearly equal.
    """
    f = fl if fl is not None else film(T_surface_C, T_bulk_C)
    rho_s = water(T_surface_C).rho
    rho_b = water(T_bulk_C).rho
    # The thermal density difference comes from the liquid; `f` supplies the
    # kinematic viscosity and reference density, which is what changes when
    # the bath is carrying bubbles.
    return g * abs(rho_b - rho_s) / f.rho * L**3 / f.nu**2


# ---------------------------------------------------------------------------
# 2. Dissolved oxygen
# ---------------------------------------------------------------------------


def do_saturation(T_C: float, salinity_ppt: float = 0.0,
                  elevation_m: float = 0.0) -> float:
    """
    Dissolved-oxygen solubility in mg/L at 100% saturation.

    Benson & Krause (1984) as adopted by APHA Standard Methods / USGS. This
    is the curve every stream gauge in the country is calibrated against.
    """
    T = T_C + 273.15
    ln_c = (
        -139.34411
        + 1.575701e5 / T
        - 6.642308e7 / T**2
        + 1.243800e10 / T**3
        - 8.621949e11 / T**4
    )
    # Salinity correction (negligible for a freshwater creek, kept for rigor).
    if salinity_ppt:
        ln_c -= salinity_ppt * (
            0.017674 - 10.754 / T + 2140.7 / T**2
        )
    c = math.exp(ln_c)  # mg/L at 1 atm
    # Barometric correction for elevation.
    if elevation_m:
        c *= math.exp(-elevation_m / 8400.0)
    return c


def water_with_dissolved_o2(T_C: float, do_mg_per_L: float | None = None,
                            saturation: float | None = None,
                            P: float = P_ATM) -> Fluid:
    """
    Water properties corrected for dissolved oxygen content.

    Oxygen dissolved in water is a true solute, so it perturbs the mixture
    properties in proportion to its *mass fraction*. At 10 C, 100% saturated
    freshwater holds ~11.3 mg/L, i.e. a mass fraction of ~1.1e-5. The
    resulting property shifts are parts-per-million.

    We model the shift with ideal mixing on a mass basis for rho, cp and k,
    using the partial-molar volume of dissolved O2 (~31 cm^3/mol, Millero)
    rather than gas-phase density -- dissolved O2 is not a bubble.
    """
    base = water(T_C, P)
    if do_mg_per_L is None:
        sat = 1.0 if saturation is None else saturation
        do_mg_per_L = sat * do_saturation(T_C)

    # Mass fraction of O2 in the solution: mg/L -> kg/m^3 -> kg per kg.
    w = (do_mg_per_L * 1e-3) / base.rho

    M_O2 = 0.0319988          # kg/mol
    V_partial = 31.0e-6       # m^3/mol, partial molar volume of aqueous O2
    rho_o2_aq = M_O2 / V_partial   # ~1032 kg/m^3 -- close to water, not to gas

    # Volume-additive density.
    v_mix = (1.0 - w) / base.rho + w / rho_o2_aq
    rho = 1.0 / v_mix

    # Mass-weighted cp. Aqueous O2 partial molar heat capacity ~ 200 J/mol-K
    # (Hnedkovsky & Wood), i.e. ~6250 J/kg-K.
    cp = (1.0 - w) * base.cp + w * 6250.0

    # Conductivity: dilute solute, Maxwell/Nernst-type shift is far below the
    # uncertainty of the base correlation. Treated as unchanged but the
    # sensitivity is retained explicitly so the magnitude is auditable.
    k = base.k * (1.0 + 0.0 * w)

    # Viscosity: aqueous O2 raises viscosity by a Jones-Dole B-coefficient
    # of roughly 0.1 L/mol for a nonelectrolyte gas.
    c_molar = (do_mg_per_L * 1e-3) / M_O2 / 1000.0   # mol/L
    mu = base.mu * (1.0 + 0.10 * c_molar)

    return Fluid(T=base.T, rho=rho, cp=cp, k=k, mu=mu, beta=base.beta)


# ---------------------------------------------------------------------------
# 3. Entrained gas (bubbles) -- the effect people actually mean
# ---------------------------------------------------------------------------


def bubbly_mixture(T_C: float, void_fraction: float, P: float = P_ATM) -> Fluid:
    """
    Homogeneous-mixture properties for water carrying a gas void fraction.

    `void_fraction` is volumetric (0.01 = 1% of the volume is air). A riffle
    or plunge pool runs 0.1%-5%; whitewater can exceed 30% locally.

    Density and heat capacity mix by volume/mass. Conductivity uses the
    Maxwell-Eucken bound for dispersed non-conducting spheres in a continuous
    liquid, which is the right model for dilute bubbles. Viscosity uses
    Einstein's dilute-suspension result with the Taylor correction for
    circulating (non-rigid) inclusions -- for gas bubbles the Taylor factor
    tends to 1.0, so mu_eff = mu_l * (1 + a).
    """
    if not 0.0 <= void_fraction < 1.0:
        raise ValueError("void_fraction must be in [0, 1)")

    liq = water(T_C, P)
    a = void_fraction
    T = T_C + 273.15

    rho_g = PropsSI("D", "T", T, "P", P, "Air")
    cp_g = PropsSI("C", "T", T, "P", P, "Air")
    k_g = PropsSI("L", "T", T, "P", P, "Air")

    rho = (1.0 - a) * liq.rho + a * rho_g
    # Mass-weighted cp so that rho*cp is the correct volumetric capacity.
    cp = ((1.0 - a) * liq.rho * liq.cp + a * rho_g * cp_g) / rho

    # Maxwell-Eucken, continuous phase = liquid.
    kr = k_g / liq.k
    k = liq.k * (2.0 + kr - 2.0 * a * (1.0 - kr)) / (2.0 + kr + a * (1.0 - kr))

    mu = liq.mu * (1.0 + 1.0 * a)   # Taylor limit for gas bubbles

    # Buoyancy: the *thermal* expansion coefficient of the liquid is what
    # drives natural convection; bubbles add a separate, much stronger
    # buoyancy source handled in the CFD, not here.
    return Fluid(T=liq.T, rho=rho, cp=cp, k=k, mu=mu, beta=liq.beta)


# ---------------------------------------------------------------------------
# Brines -- baths that stay liquid below 0 C
# ---------------------------------------------------------------------------


# Ice liquidus of the water-NaCl system: mass percent NaCl against the
# temperature at which ice first forms. Tabulated rather than fitted -- a
# cubic through the dilute end missed the eutectic by 4 K, and the eutectic is
# the whole point of a salt-ice bath.
_NACL_LIQUIDUS = (
    (0.0, 0.0), (2.0, -1.13), (4.0, -2.35), (6.0, -3.63), (8.0, -4.97),
    (10.0, -6.56), (12.0, -8.18), (14.0, -10.0), (16.0, -11.9),
    (18.0, -14.0), (20.0, -16.5), (22.0, -19.2), (23.3, -21.1),
)


def nacl_freezing_point(w: float) -> float:
    """
    Temperature at which ice first forms in aqueous NaCl, degrees C, for mass
    fraction `w`. Linear interpolation of the tabulated ice liquidus.

    Clamped at the eutectic (23.3 wt%, -21.1 C). Past that the solid that
    forms is the salt dihydrate rather than ice and the liquidus turns back
    upward, so more salt buys nothing -- which is exactly why a salt-ice bath
    bottoms out around -21 C however much salt is thrown at it.
    """
    x = min(max(w, 0.0), 0.233) * 100.0
    xs = [p[0] for p in _NACL_LIQUIDUS]
    ts = [p[1] for p in _NACL_LIQUIDUS]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            f = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ts[i] + f * (ts[i + 1] - ts[i])
    return ts[-1]


NACL_EUTECTIC_W = 0.233
NACL_EUTECTIC_T = -21.1


@lru_cache(maxsize=4096)
def _brine_freeze_cached(w_ppm: int, P: float) -> float:
    return PropsSI("T_FREEZE", "T", 273.15, "P", P,
                   f"INCOMP::MNA[{w_ppm / 1e6:.4f}]") - 273.15


@lru_cache(maxsize=400_000)
def _brine_cached(T_mK: int, w_ppm: int, P: float) -> Fluid:
    T = T_mK / 1000.0          # absolute, K
    w = w_ppm / 1e6
    fluid = f"INCOMP::MNA[{w:.4f}]"
    rho = PropsSI("D", "T", T, "P", P, fluid)

    # Thermal expansion coefficient by finite difference. The incompressible
    # backend does not expose beta, and it refuses any evaluation outside the
    # solution's liquid range -- so right at the freezing point a centred
    # stencil steps outside and throws. Pick the stencil that fits.
    T_lo = _brine_freeze_cached(w_ppm, P) + 273.15
    T_hi = 313.15
    dT = 0.5
    if T - dT <= T_lo:
        a, b_ = T, min(T + 2 * dT, T_hi)
    elif T + dT >= T_hi:
        a, b_ = max(T - 2 * dT, T_lo), T
    else:
        a, b_ = T - dT, T + dT
    rho_a = PropsSI("D", "T", a, "P", P, fluid)
    rho_b = PropsSI("D", "T", b_, "P", P, fluid)
    beta = -(rho_b - rho_a) / (b_ - a) / rho

    return Fluid(
        T=T, rho=rho,
        cp=PropsSI("C", "T", T, "P", P, fluid),
        k=PropsSI("L", "T", T, "P", P, fluid),
        mu=PropsSI("V", "T", T, "P", P, fluid),
        beta=beta,
    )


def brine_freeze_C(w: float, P: float = P_ATM) -> float:
    """
    CoolProp's own freezing temperature for the NaCl solution, degrees C.

    This is the backend's validity limit, so it is what the property calls
    have to respect. It agrees with the tabulated liquidus above to a few
    tenths of a kelvin, which is a useful cross-check on both.
    """
    w = min(max(w, 1e-4), 0.23)
    return _brine_freeze_cached(int(round(w * 1e6)), P)


# CoolProp's NaCl solution is correlated over 173.15-313.15 K.
_BRINE_T_MAX_C = 313.15 - 273.15 - 0.02


def _brine_T(T_C: float, w: float, P: float) -> float:
    """Clamp to just inside CoolProp's valid range for this solution."""
    return min(max(T_C, brine_freeze_C(w, P) + 0.02), _BRINE_T_MAX_C)


def brine_nacl(T_C: float, w: float = 0.23, P: float = P_ATM) -> Fluid:
    """
    Aqueous NaCl brine properties from CoolProp's incompressible solutions.

    A salt-and-ice bath is the classic way to get a still container well below
    0 C, so it is the natural candidate for "when does the box win". Note what
    it costs: at 23 wt% the viscosity is five times water's and the
    conductivity is 11 % lower, both of which work against the film
    coefficient. The bath wins on temperature, not on transport.
    """
    # CoolProp's NaCl solution is correlated to 23 wt%, which is essentially
    # the eutectic (23.3 wt%); clamp rather than extrapolate. Quantised to
    # 1 mK and 1 ppm so the lookup can be cached -- the transient solvers call
    # this millions of times and an uncached CoolProp incompressible-solution
    # call is slow enough to stall the ODE integrator outright.
    w = min(max(w, 1e-4), 0.23)
    T_mK = int(round((_brine_T(T_C, w, P) + 273.15) * 1000.0))
    return _brine_cached(T_mK, int(round(w * 1e6)), P)


def brine_rho(T_C: float, w: float, P: float = P_ATM) -> float:
    return brine_nacl(T_C, w, P).rho


# ---------------------------------------------------------------------------
# Ice
# ---------------------------------------------------------------------------

ICE_K = 2.22          # W/m-K at 0 C
ICE_RHO = 917.0       # kg/m^3
ICE_LATENT = 333.5e3  # J/kg


# ---------------------------------------------------------------------------
# Solids
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Solid:
    name: str
    rho: float
    cp: float
    k: float

    @property
    def alpha(self) -> float:
        return self.k / (self.rho * self.cp)


GLASS = Solid("soda-lime glass", rho=2500.0, cp=840.0, k=1.05)
PET = Solid("PET", rho=1380.0, cp=1200.0, k=0.24)
ALUMINIUM = Solid("aluminium", rho=2700.0, cp=900.0, k=205.0)
STAINLESS = Solid("304 stainless", rho=8000.0, cp=500.0, k=16.2)


if __name__ == "__main__":
    for T in (5, 10, 20, 25, 40):
        w = water(T)
        print(f"T={T:5.1f} C  rho={w.rho:7.2f}  cp={w.cp:7.1f}  k={w.k:.4f} "
              f" mu={w.mu*1e6:7.2f}e-6  beta={w.beta*1e4:6.3f}e-4  Pr={w.Pr:5.2f}")
    print()
    for T in (5, 10, 20, 25):
        print(f"DO saturation at {T:4.1f} C = {do_saturation(T):5.2f} mg/L")
    print()
    base = water(10.0)
    o2 = water_with_dissolved_o2(10.0, saturation=1.0)
    print("dissolved O2 at 100% sat, relative change vs pure water:")
    for f in ("rho", "cp", "k", "mu"):
        b, v = getattr(base, f), getattr(o2, f)
        print(f"   {f:4s}: {(v/b - 1)*100:+.6f} %")
    print()
    for a in (0.001, 0.01, 0.05, 0.20):
        m = bubbly_mixture(10.0, a)
        print(f"void={a*100:5.1f}%  rho={m.rho:7.1f} ({(m.rho/base.rho-1)*100:+6.2f}%) "
              f" k={m.k:.4f} ({(m.k/base.k-1)*100:+6.2f}%) "
              f" rho*cp={(m.rho*m.cp)/1e6:.3f} MJ/m3K ({((m.rho*m.cp)/(base.rho*base.cp)-1)*100:+6.2f}%)")
