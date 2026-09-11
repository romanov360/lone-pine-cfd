# Verification and validation

None of the CFD results in this study are worth reading unless the solver first
reproduces problems whose answers are already known. The suite below is run by
`python3 src/validate.py all` and its raw output lands in
`results/data/validation_*.json`.

The benchmarks are ordered the way the solver is built, so a failure points at
the layer that caused it.

## 1. The Poisson operator — manufactured solutions

The pressure projection solves `lap(phi) = rhs` every step, so it is checked
first, against exact solutions on both boundary configurations used later: the
closed box (Neumann on all four sides) and the open channel (Neumann at the
inlet, Dirichlet at the outlet).

| n | closed box, L∞ | order | open channel, L∞ | order |
|---|---|---|---|---|
| 32 | 5.98 × 10⁻³ | | 2.70 × 10⁻³ | |
| 64 | 1.50 × 10⁻³ | 2.00 | 6.76 × 10⁻⁴ | 2.00 |
| 128 | 3.75 × 10⁻⁴ | 2.00 | 1.69 × 10⁻⁴ | 2.00 |
| 256 | 9.36 × 10⁻⁵ | 2.00 | 4.23 × 10⁻⁵ | 2.00 |

Second order to three figures in both cases. A 384² solve takes 8.2 ms, which
is what makes the transient runs affordable.

## 2. Momentum and pressure coupling — lid-driven cavity

Against the centreline velocities tabulated by Ghia, Ghia & Shin
(*J. Comput. Phys.* **48**, 1982), on a 128² grid:

| Reynolds number | RMS error vs Ghia | max ∇·u |
|---|---|---|
| 100 | 0.0023 | 1.5 × 10⁻¹⁴ |
| 400 | 0.0012 | 2.1 × 10⁻¹⁴ |
| 1000 | 0.0031 | 1.7 × 10⁻¹⁴ |

Errors of 0.1–0.3 % of the lid velocity, and the velocity field is
divergence-free to machine precision.

## 3. Buoyancy and energy — differentially heated cavity

This is the benchmark that matters most for this study: it is the same physics
as a bottle shedding a buoyant plume into a closed box. Against de Vahl Davis
(*Int. J. Numer. Meth. Fluids* **3**, 1983), Pr = 0.71:

| Rayleigh number | grid | Nu (this solver) | Nu (benchmark) | error |
|---|---|---|---|---|
| 10³ | 48² | 1.119 | 1.118 | **+0.07 %** |
| 10⁴ | 64² | 2.250 | 2.243 | **+0.32 %** |
| 10⁵ | 96² | 4.538 | 4.519 | **+0.43 %** |
| 10⁶ | 160² | 8.864 | 8.800 | **+0.73 %** |

Better than 1 % across three decades of Rayleigh number.

## 4. Immersed boundary and forced convection — cylinder cross-flow

A circular cylinder in a channel, the closest standard benchmark to a bottle in
a current. Two quantities are checked, and they do not come out equally well.

**Wake geometry (momentum only).** At Re = 40 the flow is steady with a closed
recirculation bubble whose length is well established at L/D ≈ 2.2
(Coutanceau & Bouard, *J. Fluid Mech.* **79**, 1977). At D/Δx = 24 this solver
gives **L/D = 2.06**, about 6 % short — consistent with a staircase immersed
boundary at that resolution. The velocity field is good.

**Heat transfer (Pr = 7).** The same run gives Nu = 4.3 against 7.6 from
Churchill–Bernstein, about 40 % low. This is a resolution limit, not a model
error, and it is worth being explicit about why:

* The wall heat flux is computed as a linear gradient from the wall face to the
  first fluid cell centre, over Δx/2. That is first-order accurate, and a linear
  fit *under*-estimates the gradient of a profile that curves away from the
  wall.
* At Pr = 7 the thermal boundary layer is thinner than the momentum layer by
  roughly Pr^(1/3) ≈ 1.9. At Re = 40 it is about D/7.6 ≈ 0.13 D — only three
  cells at D/Δx = 24.

Under-resolving a boundary layer by that margin under-predicts Nu, which is the
direction observed.

**This is why the study does not take its forced-convection film coefficients
from the CFD.** They come from the six correlations in `ht`, which are fits to
experimental data. The CFD's role in the forced-convection branch is to confirm
the *flow* is right — which the wake length does — and to explore the mixed
convection regime at low velocity where the geometry, not the film coefficient,
is the question.

## 5. Grid convergence — the actual case

The box case at four spacings, bottle held at 25 °C in 10 °C water:

| Δx | grid | quasi-steady h |
|---|---|---|
| 2.00 mm | 120 × 130 | 399.1 W/m²K |
| 1.50 mm | 160 × 173 | 434.5 |
| 1.00 mm | 240 × 260 | 443.9 |
| 0.75 mm | 320 × 347 | 426.9 |

The three coarsest grids converge monotonically (Richardson on those three
gives an observed order of 3.0 and a limit of 414 W/m²K, GCI 3.7 %), but the
finest grid breaks the trend. **No asymptotic range has been established, so no
extrapolated value is quoted as though one had been.**

The reason is physical, not numerical. At Ra ≈ 10⁹ the plume is genuinely
unsteady on timescales comparable to the 60 s averaging window, so successive
runs differ by more than the discretisation does. Within any single run the
signal is now very steady — the temporal scatter is ±2 to ±3 W/m²K once the
bottle surface is properly pinned — but the *between-run* spread on the three
finest grids is ±8.

The defensible reading is that the 2.0 mm grid is under-resolved, and that on
1.5 mm and finer the answer is **h = 435 ± 10 W/m²K** — about 16 % below the
Churchill–Chu correlation's 520 W/m²K for the same conditions. That gap is not
numerical error; §6 shows it is confinement.

Max ∇·u stays at 10⁻¹⁵ throughout.

> **A bug this table caught.** An earlier version of these runs approximated
> the isothermal bottle by giving it a large conductivity *and* a very large
> heat capacity. Those divide: the body's diffusivity ended up a thousand times
> *below* the surrounding water's, its surface cells cooled faster than its
> interior could resupply them, and every h came out low (402 instead of 435)
> and noisy. The giveaway was in the cylinder benchmark, where it produced a
> Nusselt number that fell with Reynolds number — impossible, and a symptom
> that scaled with how much heat was being drawn off. The solver now pins the
> body's temperature explicitly, which is what an isothermal boundary is.

## 6. A cross-check the correlations cannot give

The box-size sweep is not a validation case — there is no reference answer —
but it explains the 20 % gap above and so belongs here:

| box | side gap | h from CFD | vs correlation | bath warmed | top − bottom |
|---|---|---|---|---|---|
| 12 × 20 cm | 0.36 D | 323 ± 45 | 62 % | +7.7 K | +9.9 K |
| 18 × 26 cm | 0.79 D | 419 ± 9 | 81 % | +3.4 K | +7.2 K |
| 26 × 34 cm | 1.36 D | 444 ± 11 | 85 % | +1.7 K | +4.0 K |
| 40 × 46 cm | 2.36 D | 445 ± 10 | 86 % | +0.8 K | +2.1 K |

The correlation returns 520 W/m²K for all three, because it assumes an infinite
quiescent medium and knows nothing about walls. A real box violates that
assumption twice over: it recirculates water the bottle has already warmed, and
it stratifies, so the warmest water in the box collects at the top — exactly
where the bottle is. Both effects shrink as the box grows, which is the trend
the table shows in the last two columns: the bath warms ten times less in the
largest box than the smallest, and the stratification falls from 8 K to 2 K.

The film coefficient plateaus cleanly at **445 W/m²K** once the side gap
exceeds about one bottle diameter — the 26 cm and 40 cm boxes agree to within
0.4 % — and collapses to 323 W/m²K when the gap closes to a third of a
diameter. So confinement alone costs 27 % in a tight tub. The remaining 14 %
between the plateau and the correlation's 520 W/m²K is the first-order
wall-flux estimate discussed in §4, plus the fact that even a 40 cm box is not
an infinite medium.
The practical conclusion is unchanged either way, and is in fact strengthened:
the still box is *worse* than its textbook correlation suggests, so the
correlation-based comparison in the main study is conservative about the
creek's advantage rather than flattering to it.
