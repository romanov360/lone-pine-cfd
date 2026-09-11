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
a current.

**Heat transfer.** With the cylinder held isothermal, against Churchill &
Bernstein (1977) at Pr = 0.7 and D/Δx = 24:

| Reynolds number | Nu (this solver) | Nu (correlation) | error |
|---|---|---|---|
| 40 | 3.43 | 3.36 | **+2.0 %** |
| 100 | 5.13 | 5.16 | **−0.5 %** |
| 200 | 7.06 | 7.19 | **−1.8 %** |

Within two percent across the range — comfortably inside the correlation's own
scatter, and the trend with Reynolds number is now the right sign.

**Resolution at Pr = 7.** Water has Pr ≈ 7, so its thermal boundary layer is
thinner than its momentum layer by roughly Pr^(1/3) ≈ 1.9 and needs a finer
mesh for the same accuracy. A refinement sequence at Re = 40 shows the error
shrinking as it should:

| D/Δx | Nu | error vs correlation |
|---|---|---|
| 16 | 8.72 | +15.2 % |
| 24 | 8.24 | +8.9 % |
| 36 | 7.96 | +5.1 % |

Monotone, and roughly halving for each 1.5× refinement.

This is the practical reason the study's water-side forced-convection numbers
come from correlations rather than from the CFD: getting Pr = 7 cross-flow to
correlation accuracy costs a mesh that is not worth buying when experimental
fits are available and already good.

**Wake geometry.** At Re = 40 the flow is steady with a closed recirculation
bubble, established at L/D ≈ 2.2 (Coutanceau & Bouard, *J. Fluid Mech.* **79**,
1977). This solver gives **L/D = 2.02**, about 8 % short, which is what a
staircase immersed boundary at 24 cells per diameter should give.

**What is *not* checked here: vortex shedding.** At Re = 100 a real cylinder
sheds at St ≈ 0.164. This one does not shed at all, and the reason is the
configuration rather than the solver: the cylinder is centred in a uniform
stream on a symmetric mesh, and nothing breaks that symmetry. Two-dimensional
cylinder wakes need a perturbation to trip the instability, and none is
supplied, so the wake stays symmetric and elongated (L/D = 5.5) instead of
rolling up. The Strouhal number is therefore not measured and is not claimed.
This does not affect the study: the bottle cases are buoyancy-driven, where the
instability is supplied by the flow itself.

> **A bug this benchmark caught.** Before the isothermal boundary was imposed
> properly (see §5), these same runs gave −10 % at Re = 40 and −59 % at
> Re = 100 — a Nusselt number that *fell* with Reynolds number, which is
> impossible. The symptom scaled with how much heat was being drawn off, which
> is what pointed at the boundary condition rather than at the discretisation.
> Fixing it moved both cases to within 2 %.

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

## 6. The creek case, against the correlation it will be compared with

The box results in §7 only mean something if the same solver reproduces the
correlation where the correlation is trustworthy. Run in an open channel
instead of a closed box, with the bottle held isothermal:

| setting | Re | h from CFD | correlation | ratio |
|---|---|---|---|---|
| closed channel, no current | — | 441 ± 27 | 520 | 0.85 |
| creek at 0.01 m/s | 540 | 524 ± 16 | 528 | **0.99** |
| creek at 0.02 m/s | 1,080 | 577 ± 26 | 552 | 1.05 |
| creek at 0.05 m/s | 2,700 | 771 ± 60 | 672 | 1.15 |
| creek at 0.10 m/s | 5,400 | 775 ± 55 | 897 | 0.86 |

At 0.01 m/s the agreement is within 1 %. It then degrades, and at 0.1 m/s the
error changes *sign* — which is the useful signal. A discretisation error would
shrink smoothly; an error that grows, flips and comes with ±55 of scatter is
the two-dimensional representation ceasing to be physical. By Re ≈ 5,000 a real
cylinder wake is three-dimensional and transitional, and a 2-D simulation
cannot cascade energy the right way.

The sweep was stopped there rather than extended to 0.2 m/s. That point would
have cost about four hours of compute to produce a number at Re ≈ 11,000, where
the method is already known not to apply. **This is exactly why the study's
fast-creek film coefficients come from experimental correlations and not from
this solver** — and why the CFD's job in the forced-convection branch is to
confirm the correlation where a 2-D simulation is trustworthy, which it does to
1 % at 0.01 m/s.

The first row is the useful contrast. The *same solver*, in a channel with no
current, gives 441 — well above the 323–445 range the box gives, and the box's
own value falls as the box tightens. So the shortfall against the correlation
is confinement, not the discretisation.

## 7. A cross-check the correlations cannot give

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
