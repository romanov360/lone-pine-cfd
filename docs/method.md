# Method

The question is simple to ask and awkward to answer with one tool, because it
spans three regimes that no single model covers well:

* a **boundary-layer** problem (what is the film coefficient?),
* a **conjugate conduction** problem (the glass wall and the water inside),
* a **system energy balance** problem (does the bath warm up?).

So the study uses four levels of model, arranged so that each one calibrates or
checks the next. Every level is run over the same cases, and where they
disagree the disagreement is reported rather than averaged away.

## Level 0 — Properties (`src/props.py`)

Liquid water properties come from CoolProp's IAPWS-95 backend, evaluated at the
film temperature and re-evaluated at every step of every transient. Three
additions:

**Dissolved oxygen** (`water_with_dissolved_o2`). Saturation from Benson &
Krause (1984), the curve used by USGS stream gauges. The solute is mixed on a
mass basis using the *partial molar volume* of aqueous O₂ (≈31 cm³/mol) rather
than gas-phase density — dissolved oxygen is not a bubble, and treating it as
one would overstate its effect by three orders of magnitude.

**Entrained gas** (`bubbly_mixture`). Homogeneous mixture: volume-weighted
density, mass-weighted c_p so that ρc_p is the correct volumetric capacity,
Maxwell–Eucken conductivity for dispersed non-conducting spheres, and Einstein
viscosity in the Taylor limit for circulating inclusions.

**Density-difference Grashof number** (`Gr_density`). The usual
Gr = g·β·ΔT·L³/ν² fails near an ice bath, because water's density maximum at
3.98 °C makes β negative below it and the correlations return NaN. That is not
a numerical nuisance to be clamped — it is real physics. The code uses

    Gr = g · |ρ(T_bulk) − ρ(T_surface)| / ρ_film · L³ / ν²

which is exact for any temperature pair, reduces to the linearised form away
from 4 °C, and correctly predicts that buoyancy nearly vanishes when the two
temperatures straddle the density maximum.

A **conduction floor** is imposed on every free-convection coefficient: a body
in a stagnant infinite medium still loses heat, with Nu → 2 on the
equal-area-sphere diameter. Without it the correlations would claim a bottle in
4 °C water never cools at all.

## Level 1 — Correlations (`src/correlations.py`, on the `ht` library)

* Free convection, vertical: Churchill–Chu with the Popiel–Churchill curvature
  correction (a 70 mm bottle is slender enough for it to matter).
* Free convection, horizontal: Churchill–Chu for a horizontal cylinder.
* Forced convection: six independent correlations — Churchill–Bernstein,
  Žukauskas, Whitaker, Sanitjai–Goldstein, McAdams, Fand — carried in parallel
  so the *correlation uncertainty* (±20 %) is visible next to every result.
* **Mixed convection**: Churchill's asymptotic blend, n = 4 for a buoyant plume
  transverse to the oncoming stream, which is exactly an upright bottle in a
  creek.

The mixed form matters more than it looks. A creek at U → 0 must reduce
continuously to the still-water case; otherwise the whole comparison is an
artifact of using two unrelated correlations on the two sides. `h_external()` is
one function producing one continuous curve from U = 0 upward, so the box and
the creek are never compared across a modelling seam.

## Level 2 — Transient system model (`src/lumped.py`)

Two bodies: the bottle (water core plus glass shell) and the bath. The bath is
either a finite box that warms as it absorbs heat, or a creek whose temperature
never moves. Integrated with SciPy's LSODA, properties and film coefficients
recomputed at every step.

The outer wall temperature is not known in advance — both films depend on their
own temperature drop — so each evaluation solves

    residual(T_s) = q_in(T_s) − q_out(T_s) = 0

by bracketed Brent iteration. A damped fixed-point iteration was tried first and
abandoned: near the density maximum h_out is so steep a function of T_s that the
iterate never settles to a *repeatable* value, which leaves the ODE integrator
unable to estimate its own error, and it stalls indefinitely. The bracketing
root solve is deterministic and the stall disappears.

This level carries all the parameter sweeps, because it is the only one cheap
enough to run hundreds of times.

## Level 3 — Axisymmetric conjugate model (`src/axisym_fipy.py`, FiPy)

A true (r, z) finite-volume model on FiPy's `CylindricalGrid2D`, where the cell
volumes carry the 2πr weighting and the geometry is the real cylinder rather
than a Cartesian slab. Water core, glass shell, and a thin outer annulus whose
conductivity is set to h·t so it represents the external film exactly.

Its job is to price the *internal* resistance, which the lumped model can only
assume. Internal circulation is represented as an effective conductivity
k_eff = Nu_i · k_water, and the model is run across that multiplier:

| k_eff multiplier | → 25 °C | → 15 °C | centre-to-edge ΔT at 5 min |
|---|---|---|---|
| 1 (pure conduction, no stirring) | 9.2 min | 31.2 min | 16.9 K |
| 5 | 4.5 min | 12.4 min | 10.5 K |
| 20 | 3.4 min | 8.8 min | 2.6 K |
| 100 (perfectly mixed) | 3.0 min | 7.7 min | 0.4 K |
| *lumped 0-D model, for comparison* | *5.1 min* | *15.6 min* | *(assumes 0)* |

Two things follow. First, whether the water inside the bottle circulates is
worth a **factor of three** — considerably more than the entire box-versus-creek
question. Second, the 0-D lumped model corresponds to k_eff ≈ 4–5, which is a
reasonable value for buoyant circulation inside a sealed vessel, so the two
levels are consistent.

## Level 4 — CFD (`src/nsolver.py`, `src/case_bottle.py`, `src/campaign.py`)

A purpose-built 2-D incompressible solver:

* Staggered MAC grid — u on x-faces, v on y-faces, p and T at centres.
* van Leer TVD advection. Cell Reynolds numbers reach O(50) here, well past
  where central differencing stays stable.
* Explicit viscous and conductive terms (at these spacings the diffusive step
  limit is far above the CFL limit, so implicit treatment would buy nothing).
* Chorin projection with a **separable Poisson solver**: cosine transform in y,
  a Thomas sweep vectorised across all y-modes in x. O(N log N), about 8 ms at
  384², roughly 100× faster than a sparse LU — which is what makes the
  transient runs affordable at all.
* The immersed solid by **Brinkman volume penalisation**.
* Conjugate energy equation on a single domain with harmonic-mean face
  conductivity, which enforces interface flux continuity exactly.

### Two modelling decisions worth stating plainly

**The bottle interior is not resolved as a free fluid.** It was tried. A 3 mm
glass wall is two or three cells thick, and no immersed-boundary treatment
holds a membrane that thin: hard-masking the solid faces destroyed the
divergence-free condition (max|div| went from 10⁻¹⁴ to O(1)), while penalisation
left ~2 × 10⁻³ m/s inside the wall — enough to carry the bottle's contents
bodily into the bath within twenty seconds. Both failure modes were caught by
the energy-conservation diagnostic. The interior is therefore carried as a
conducting region at k_eff, with the multiplier taken from Level 3 rather than
guessed, and the sensitivity reported.

**The fast-creek film coefficients come from correlations, not CFD.** A creek at
0.35 m/s past a 70 mm bottle is Re ≈ 23,000: three-dimensional, turbulent, and
with a thermal boundary layer around 0.2 mm. Resolving that in 2-D would not be
more truthful than the correlations, it would just be more expensive and wrong
in a different way — two-dimensional turbulence cascades the wrong direction.
The CFD is run where a 2-D simulation is defensible (Re ≲ 3,000, and the
buoyancy-driven box at Ra ≈ 10⁹), validated against the correlations there, and
the correlations carry the extrapolation.

### What the CFD is actually for

Three things the correlations cannot tell you, all of which turned out to
matter:

1. **Confinement.** Textbook free-convection correlations assume an infinite
   quiescent medium. A real box is not one.
2. **Stratification.** A still box does not warm uniformly; the bottle's own
   plume pools warm water at the top, around the bottle.
3. **The mixed-convection regime** at slow creek velocities, where neither the
   free nor the forced correlation is clearly in charge.

## Verification and validation

Every numerical component is checked before use; see
[validation.md](validation.md). In order: the Poisson solver against
manufactured solutions, the momentum solver against Ghia et al.'s lid-driven
cavity, the coupled buoyancy-energy solver against de Vahl Davis's heated
cavity, the immersed boundary and forced convection against cylinder cross-flow
data, and the grid-convergence index following Roache/Celik — with the observed
order found by iteration, because the grid sequence has a non-constant
refinement ratio and the closed-form version returns a negative (meaningless)
order on this data.
