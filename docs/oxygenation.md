# Does oxygenation matter?

You asked specifically about the difference in oxygenation between aerated
creek water and still box water. It is worth separating carefully, because the
phrase "aerated water" covers three physically different things that push in
different directions, and only one of them matters.

## 1. Dissolved oxygen: no

Dissolved oxygen is a genuine solute. Its concentration is set by Henry's law
and is strongly temperature dependent — `props.do_saturation()` implements the
Benson & Krause (1984) curve that every USGS stream gauge is calibrated
against:

| water temperature | 100 %-saturated dissolved O₂ |
|---|---|
| 0 °C | 14.62 mg/L |
| 5 °C | 12.77 mg/L |
| 10 °C | 11.29 mg/L |
| 20 °C | 9.09 mg/L |
| 25 °C | 8.26 mg/L |

A fully aerated creek at 10 °C therefore carries about 11.3 mg of O₂ per litre.
That is a **mass fraction of 1.1 × 10⁻⁵**. Treating it as an ideal solution with
the measured partial molar volume of aqueous O₂ (≈ 31 cm³/mol, so the dissolved
oxygen has an effective density near 1030 kg/m³ — it is *not* a gas once
dissolved) and its partial molar heat capacity, the property shifts are:

| property | change from 0 % to 100 % saturation |
|---|---|
| density ρ | +3.6 × 10⁻⁵ % |
| heat capacity c_p | +5.5 × 10⁻⁴ % |
| conductivity k | below the noise floor of the base correlation |
| viscosity μ | +3.5 × 10⁻³ % |

Propagated through the convection correlations, going from completely
deoxygenated water to supersaturated water (200 %) changes the film
coefficient by **−0.001 %**. For comparison, the six standard cross-flow
correlations in the `ht` library disagree with each other by ±20 % on the same
case. Dissolved oxygen is four to five orders of magnitude below the
uncertainty of anything else in the problem.

There is no thermal mechanism by which dissolved oxygen could matter. It does
not change the phase behaviour (the bottle is submerged, so there is no
evaporation), it does not change the buoyancy meaningfully, and it carries no
latent heat. Its biological importance in a creek is real and large; its
thermal importance here is nil.

## 2. Entrained bubbles: yes, and they *hurt* the fluid properties

What people usually see and call "aerated" is not dissolved gas — it is
entrained gas. Actual bubbles, occupying actual volume. A riffle or the tail of
a plunge pool runs a void fraction of roughly 0.1 % to 5 %; whitewater exceeds
30 % locally.

Bubbles displace water, and water is what carries the heat. Using
Maxwell–Eucken for the conductivity of dispersed non-conducting spheres and
volume mixing for the capacity (`props.bubbly_mixture`):

| void fraction | ρ | k | ρc_p (volumetric capacity) |
|---|---|---|---|
| 0.1 % | −0.10 % | −0.14 % | −0.10 % |
| 1 % | −1.00 % | −1.40 % | −1.00 % |
| 5 % | −4.99 % | −6.86 % | −5.00 % |
| 20 % | −19.98 % | −25.68 % | −19.99 % |

Every one of these is the wrong direction for heat transfer. On properties
alone, aeration is a **penalty**, and at 1 % void it is already about 3,000
times larger than the entire dissolved-oxygen effect.

## 3. Bubble-driven agitation: yes, and this is the one that helps

Rising bubbles drag liquid with them. This is why bubbling air through a tank
is a standard way to improve heat transfer in it. The liquid circulation
velocity follows the plume scaling used for gas-agitated vessels (Sahai &
Guthrie 1982),

    U_liquid ≈ 1.4 (g · j_g · L)^(1/3)

with j_g the superficial gas velocity. A 1 % void fraction of 3 mm bubbles
rising at their 0.23 m/s slip velocity gives j_g ≈ 2.3 mm/s and an induced
liquid velocity of about **0.22 m/s** past a 160 mm bottle.

That is the whole story. Aeration does not help by adding oxygen; it helps by
adding *motion*. And so its benefit depends entirely on how much motion was
already there:

| | still box | creek 0.05 m/s | creek 0.35 m/s | creek 1.0 m/s |
|---|---|---|---|---|
| void 0.2 % | **+111 %** | +59 % | +3 % | +0.1 % |
| void 1 % | **+178 %** | +104 % | +8 % | −0.4 % |
| void 5 % | **+251 %** | +156 % | +14 % | −5 % |
| void 15 % | **+263 %** | +163 % | +8 % | −19 % |
| void 30 % | +207 % | +123 % | −11 % | **−38 %** |

Read the corners. Bubbling air into the **still box** nearly quadruples its
film coefficient, because 0.2 m/s of induced circulation is enormous compared
with the 0.05 m/s buoyant plume it replaces. Bubbling air into a **fast creek**
makes things worse, because the creek was already supplying more velocity than
the bubbles can add, and all the bubbles do is dilute the water.

## The conclusion that actually matters

Aeration is a red herring for the question as posed, for a reason that has
nothing to do with any of the above: **the outside film is not what limits this
bottle.** At a creek velocity of 0.35 m/s the external film is only about 10 %
of the total resistance — the glass wall is 53 % and the water inside the
bottle is 37 %. Tripling the external film coefficient by aerating changes the
overall U by a few percent.

So the honest ranking of the oxygenation question is:

1. Dissolved oxygen — utterly irrelevant (10⁻³ % effects).
2. Bubble properties — a real but modest penalty (−1 % per 1 % void).
3. Bubble agitation — a large benefit in still water, a penalty in fast water.
4. All of it together — still second-order for a *glass* bottle, because the
   bottle's own wall and contents dominate. It would matter much more for a
   thin aluminium can, where the outside film is the controlling resistance.
