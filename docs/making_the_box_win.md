# When does the box win?

The main study found that a creek beats a still box of the same temperature,
always, and that the margin is smaller than intuition suggests. This is the
inverse question: what would have to be true for the box to win?

## The structural obstacle

Everything in this study blends free and forced convection with Churchill's
asymptotic rule,

    h_mixed = (h_forced ⁿ + h_free ⁿ)^(1/n)

which is **bounded below by the free-convection limit for any velocity**. So at
equal bath temperature, with an unlimited bath and a current that is transverse
to or assisting the plume, the creek cannot lose. That is not a result — it is a
property of the model, and it is worth stating so the real answers are not
confused with it.

Every genuine box win therefore has to break one of those three assumptions.
There turn out to be five ways to do it, and they are very unequal.

| | lever | worth |
|---|---|---|
| 1 | **Bath temperature** | 2–11 K of bath temperature ≈ the whole creek |
| 2 | **Phase change** (ice) | unlimited reservoir *and* a lower temperature |
| 3 | Shelter | erodes the creek's margin 1.55× → 1.07×, never reverses it |
| 4 | Opposed flow | can halve h, but only in a narrow velocity band |
| 5 | Scale and drive | creek's edge falls 5.0× → 2.4× across bottle sizes |

## Lever 1 — the exchange rate between temperature and current

How cold must a **still** bath be to match a creek at 10 °C?

| creek speed | to reach 25 °C | to reach 20 °C | to reach 15 °C |
|---|---|---|---|
| 0.05 m/s | 3.6 K colder | 3.2 K | 2.4 K |
| 0.20 m/s | 9.2 K | 6.8 K | 4.4 K |
| 0.35 m/s | 10.4 K | 8.5 K | 5.1 K |
| 1.00 m/s | 11.1 K | 10.4 K | 6.0 K |
| 2.50 m/s | 11.5 K | 10.7 K | 6.8 K |

**Two to eleven kelvin.** Swapping tap water for fridge water — 10 °C to 4 °C —
is worth about as much as finding a brisk creek.

Note the trend *reverses* with target temperature: matching a 0.35 m/s creek
down to 25 °C needs a bath 10.4 K colder, but matching it down to 15 °C needs
only 5.1 K. The early phase is governed by how fast heat crosses the film, where
the creek is strong. The last few degrees are governed by the bath temperature
itself, where it is not.

## Lever 2 — the cold-bath ladder

| bath | → 25 °C | → 15 °C | → 10 °C | → 5 °C | → 1 °C |
|---|---|---|---|---|---|
| Creek, 0.35 m/s at 10 °C | 3.7 min | 10.7 min | **never** | never | never |
| Creek, 1.5 m/s at 10 °C | 3.5 min | 10.1 min | **never** | never | never |
| Still box, tap water 10 °C | 5.2 min | 16.6 min | never | never | never |
| Still box, fridge water 4 °C | 4.1 min | 10.5 min | 18.6 min | 3.0 h | never |
| Ice-water bath, 0 °C | 3.8 min | 9.5 min | 14.5 min | 25.3 min | 54.0 min |
| Salt-ice, 10 wt% (−6.6 °C) | 2.6 min | 5.6 min | 8.1 min | 12.8 min | 18.1 min |
| Salt-ice, 20 wt% (−16.5 °C) | 2.1 min | 4.2 min | 5.8 min | 8.6 min | 11.1 min |
| Salt-ice, 23 wt% (−20.5 °C) | **1.9 min** | **3.8 min** | **5.3 min** | **7.7 min** | **9.7 min** |

A 10 °C creek cannot reach 10 °C at any velocity — that *is* its temperature.
Everything from the ice bath down goes somewhere no 10 °C creek can follow, and
a salt-ice bath reaches 1 °C in under ten minutes while sitting perfectly still
in a bucket.

**It is also safe.** Freezing 500 mL of water takes 167 kJ against roughly
2.1 kJ per kelvin of sensible heat, so the latent plateau is far longer than the
cooling that preceded it:

| bath | margin before the contents freeze solid |
|---|---|
| Salt-ice, 10 wt% | 20.3 h |
| Salt-ice, 20 wt% | 8.0 h |
| Salt-ice, 23 wt% | 6.4 h |

The pinned-temperature assumption is physical here for a specific reason: a
salt-**ice** bath is a slurry. Near the warm bottle the ice melts, absorbing
latent heat locally, which holds the bath at its liquidus however hard the
bottle pushes. A bath of cold brine with no ice in it would warm up like any
other finite reservoir.

## Lever 3 — the creek the bottle actually feels

A bottle wedged behind a rock, or lying in the bed boundary layer, does not see
the depth-averaged velocity:

| shelter | velocity felt | → 15 °C | vs the 20 L box |
|---|---|---|---|
| 100 % | 0.350 m/s | 10.7 min | 1.55× |
| 50 % | 0.175 | 11.2 min | 1.48× |
| 30 % | 0.105 | 11.7 min | 1.42× |
| 15 % | 0.052 | 12.6 min | 1.32× |
| 5 % | 0.017 | 14.3 min | 1.16× |
| 0 % | 0.000 | 15.6 min | 1.07× |

Shelter erodes the margin but never reverses it, because a zero-velocity creek
is still an *infinite* bath, and that alone beats a 20 L one. Shelter is a
multiplier on the other levers, not a lever on its own.

From the log law over a rough gravel bed, a 70 mm bottle lying on the bed sees
56–79 % of the depth-averaged velocity — a real effect, but not a large one.
The big shelter factors need an actual obstruction.

## Lever 4 — the one direction where a current hurts

A bottle in a downwelling — the throat of a plunge pool — meets water moving
*down* past a plume trying to rise. Churchill's blend subtracts in that case:

| velocity | h forced | h free | assisting | **opposing** |
|---|---|---|---|---|
| 0.010 m/s | 260 | 466 | 477 | 454 |
| 0.020 | 370 | 466 | 507 | 411 |
| **0.030** | 456 | 466 | 548 | **251** |
| 0.050 | 595 | 466 | 645 | 529 |
| 0.100 | 862 | 466 | 879 | 842 |

Near 0.03 m/s the two mechanisms nearly cancel and moving water transfers *less*
heat than still water. This is the only configuration in the whole study where a
current is a liability.

Treat the magnitude with caution: the asymptotic blend is least reliable exactly
at the crossover it is predicting, and a real opposing flow separates rather
than cancelling cleanly. The *sign* is a well-documented effect — opposed mixed
convection does dip below both pure limits — but the depth of the dip here is
the model's, not a measurement.

## Lever 5 — scale and drive

Free convection is nearly size-independent once the Rayleigh number is large
(Nu ~ Ra^⅓ makes the length cancel), while forced convection falls as D^−½. And
free convection grows as ΔT^⅓ while forced convection does not scale with ΔT at
all. Both push the same way:

| bottle | h free | h forced | creek's advantage |
|---|---|---|---|
| 20 mL vial | 517 | 2588 | **5.00×** |
| 500 mL bottle | 466 | 1738 | 3.73× |
| 8.9 L jerrycan | 436 | 1212 | 2.79× |
| 37.7 L drum | 426 | 1027 | **2.43×** |

| bottle temperature | h free | h mixed | creek's advantage |
|---|---|---|---|
| 15 °C (barely warm) | 238 | 1670 | **7.03×** |
| 40 °C | 494 | 1751 | 3.55× |
| 95 °C (near boiling) | 877 | 1932 | **2.20×** |

A big hot drum drives its own plume hard enough that a creek adds little. A
small tepid vial has almost no plume of its own and the creek is worth
everything.

## What CFD says about designing the box

Two things a person can actually do, both measured with the bottle held
isothermal at 25 °C in 10 °C water.

### Keep the bottle off the top

A still box stratifies, and the warm layer is one the bottle built itself. Where
it sits in that layer matters (box 24 × 34 cm):

| bottle height | h (W/m²K) | stratification it sits in |
|---|---|---|
| 30 % below centre | 434 ± 16 | +2.0 K |
| 15 % below centre | 455 ± 14 | +3.1 K |
| centre | 443 ± 8 | +4.1 K |
| 15 % above centre | 410 ± 12 | +7.2 K |
| 30 % above centre | **332 ± 11** | **+9.3 K** |

Putting the bottle high in the box costs **25 %**. Putting it low buys nothing
over the centre — the three lower positions are within each other's scatter. The
rule is one-sided: *don't let it float near the surface.*

### Give it a diameter of clearance, and some headroom

Box shape at constant volume (720 cm² in section):

| box | side gap | head room | h (W/m²K) |
|---|---|---|---|
| 14 × 51 cm | 0.50 D | 1.09 H | 406 ± 23 |
| 19 × 38 cm | 0.86 D | 0.67 H | 443 ± 15 |
| 24 × 30 cm | 1.21 D | 0.43 H | 440 ± 13 |
| 31 × 23 cm | 1.71 D | 0.22 H | 439 ± 6 |
| 40 × 18 cm | 2.36 D | 0.06 H | **392 ± 3** |

A plateau with a penalty at each end. Too narrow and the plume has no return
path down the sides; too shallow and it has nowhere to rise and turn over. The
three middle shapes are indistinguishable. So the requirement is a threshold,
not an optimum: **roughly one bottle diameter of side clearance and a fifth of a
bottle height above it**, after which shape stops mattering and only volume
does.

## The verdict

The box wins whenever you use the one thing it has that a creek does not: **you
choose what is in it.**

* Ice water at 0 °C beats any 10 °C creek outright, at any velocity, and 245 g
  of ice holds a 20 L box there for the whole job.
* Salt and ice takes it to −21 °C and reaches 1 °C in under ten minutes, with
  hours of margin before the bottle is in danger.
* A warm summer creek loses to a bucket of tap water.

The box loses whenever the comparison is made at equal temperature with an
unlimited bath, because then the only thing left to compete on is the flow
field, and a still box has none.
