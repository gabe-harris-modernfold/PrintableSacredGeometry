# The human thermal plume

A thermal lattice-Boltzmann study of the natural-convection plume that rises off a standing
person, driven by **measured** per-segment skin temperatures, and analysed with the same
vortex-identification pipeline the experimental PIV groups use.

The question this project exists to answer precisely: *does a human body produce a Rankine
vortex from its infrared emission?* The short answer is no, and the interesting part is why the
question is malformed in two separate places at once.

## The two confusions, separated

**Infrared does not stir the air.** Air is very nearly transparent at the 9–10 µm wavelengths
skin emits, so that radiation crosses the room and warms walls, furniture and other people. It
is not a negligible channel — over the torso it is the *larger* one, 60% of the combined
coefficient — but almost none of it ends up as air motion. The air moves because of conduction
across the few millimetres touching skin, followed by buoyancy. Those are different physics with
different governing numbers.

**A Rankine vortex is not a cause.** It is a kinematic description of a velocity field that
already rotates: solid-body rotation inside a core radius, irrotational flow outside it. The
same model is applied to tornadoes, bathtub drains and optical vortices in partially coherent
light. Saying a plume "is a Rankine vortex" says nothing about what set it turning; it is a
ruler, not a mechanism.

What *is* real, measured and reproducible is the buoyant plume, and the fact that it sheds
discrete rotating structures above the head every few seconds. That is what this project
simulates and measures.

## Files

| File | What it does |
|---|---|
| `body_model.py` | Silhouette geometry plus the measured temperature and heat-transfer tables. Run it directly to print the dimensionless groups. |
| `plume_lbm.py` | The solver. D2Q9 momentum (TRT collision, Guo forcing) + D2Q5 temperature, Boussinesq coupling. |
| `vortex_analysis.py` | Q-criterion, swirling strength λ_ci, core detection, Rankine and Oseen profile fits, shedding spectrum. |
| `figures.py` | Renders every figure and writes `summary.json`. |
| `build_body_asset.py` | Traces the body outline from `Thermal-Infrared-Full-Body-Scan.jpg` into an RGBA sprite plus a packed lattice collision mask, the scan's own recovered colour ramp, and per-segment label anchors — all emitted to `body_mask.js` and pasted into the exhibit. |
| `refresh_shedding.py` | Re-runs only the spectral analysis and its figure. |
| `refresh_rankine.py` | Re-renders the Rankine/Oseen figure from `core_fits.json` without refitting. |

```bash
python body_model.py                                              # dimensionless groups
python plume_lbm.py --calibrate                                   # stability sweep
python plume_lbm.py --ra 1e7 --dx 0.015 --seconds 90 --smag 0 --out run_ra1e7
python figures.py --run run_ra1e7
```

Run these from inside this directory — the scripts read and write bare relative filenames.

## Boundary conditions, and where they come from

Nothing in the wall condition is a modelling convenience. Each segment carries a measured skin
temperature, and the resulting buoyancy is deliberately non-monotonic up the body:

| Segment | Skin °C | ΔT vs 24.4 °C ambient | h_r | h_c | Radiative share |
|---|---:|---:|---:|---:|---:|
| Anterior neck | 35.0 | 10.6 | 4.1 | 3.6 | 53% |
| Head / face | 34.8 | 10.4 | 4.1 | 3.6 | 53% |
| Hand (palm) | 33.8 | 9.4 | 4.1 | 4.1 | 50% |
| Chest * | 33.6 | 9.2 | 4.5 | 3.0 | 60% |
| Abdomen * | 32.8 | 8.4 | 4.4 | 2.9 | 60% |
| Lower leg (shin) | 32.7 | 8.3 | 5.3 | 4.1 | 56% |
| Upper arm / forearm | 32.6 | 8.2 | 5.2 / 4.9 | 2.9 / 3.7 | 64% / 57% |
| Pelvis * | 32.2 | 7.8 | 4.2 | 3.4 | 55% |
| Thigh | 31.9 | 7.5 | 4.3 | 4.1 | 51% |
| Sole of foot | 30.0 | 5.6 | 3.9 | 5.1 | 43% |

`*` Chest, abdomen and pelvis are not among the 15 exposed sites of the supine SST study; their
values follow the separate mapping result that the trunk runs warmest overall, and are flagged
`derived=True` in the code.

Two features of this table are worth stating out loud, because both are counterintuitive and
both matter:

- **The shin is warmer than the thigh above it** (32.7 vs 31.9 °C). Buoyancy along the leg does
  not increase monotonically with height.
- **Feet have the lowest radiative and the highest convective coefficient; the chest has the
  opposite.** Thin distal limbs shed heat convectively far better than the broad torso, which
  grows a thick, insulating boundary layer. The foot/chest convective ratio is 5.1/3.0 ≈ 1.7.

## Dimensionless groups

With the anterior neck's 10.6 K excess as the reference and body height as the length scale:

```
Ra over body height     4.88e9
free-fall velocity      0.765 m/s      sqrt(g beta dT L)
buoyant time scale      2.22 s         L / U_ff
measured shed interval  ~5 s           = 2.25 buoyant times, St = f L / U_ff ~ 0.44
```

The reported laminar→turbulent transition at 0.80 m of body height corresponds to a **local**
Rayleigh number of 3.8×10⁸ for the 306 K / 298 K conditions it was quoted at, which sits in the
accepted 10⁸–10⁹ band for a vertical surface. That number, not the height, is the transferable
quantity.

## Method

Thermal lattice Boltzmann, all quantities nondimensional inside the solver:

- **Momentum** — D2Q9, two-relaxation-time collision with the magic parameter Λ = 3/16, Guo
  forcing for buoyancy.
- **Temperature** — D2Q5 advection–diffusion, relaxation set by Pr = 0.71.
- **Body** — halfway bounce-back for momentum, halfway anti-bounce-back for temperature, which
  pins each wall node at its own segment's measured skin temperature.
- **Domain** — no-slip cool floor, open sides and top, with a sponge layer below the top
  boundary well above every probe.
- **No turbulence closure.** A Smagorinsky term was implemented and tested; at these reduced
  Rayleigh numbers it changes peak velocity by about 3% and costs 40% of the throughput, so the
  production run is a reduced-Ra laminar DNS and is described as one.

### Two implementation details that are load-bearing

Both were found by the simulation failing, and both are easy to get wrong silently:

1. **Guo forcing must use τ⁻, not τ⁺.** Momentum is an odd moment, so under TRT it relaxes at
   τ⁻. Only that choice makes the per-step momentum increment equal F. With τ⁺ ≈ 0.51 the
   forcing delivers roughly 5% of the intended buoyancy — the solver stays perfectly stable and
   the plume simply never develops, which reads as a physics result rather than a bug.
2. **Reflect a population only if it lands on fluid.** Applying bounce-back or anti-bounce-back
   to populations aimed at another solid cell lets adjacent solid cells reflect off each other.
   Anti-bounce-back has eigenvalue −1, so this grows an unbounded alternating mode in
   temperature — first visible in the gap between the feet and along the foot/floor junction.

The analysis pipeline is unit-tested against a synthetic Oseen vortex of known circulation and
core radius on a sheared, drifting background: it recovers Γ to 0.3% and the core radius to
0.1%, locates the centre to 0.02 cells, and λ_ci correctly ignores the background shear that
vorticity lights up on (2×10⁻⁷ against 8.5×10⁻⁵).

## Rankine versus Oseen

Each detected core's azimuthal profile is fitted twice, with the same number of free parameters:

```
Rankine   v(r) = Gam/(2 pi a^2) r     r < a ;   Gam/(2 pi r)   r > a
Oseen     v(r) = Gam/(2 pi r) [1 - exp(-1.2564 r^2/a^2)]
```

Rankine has a corner at r = a. A viscous core does not. On the synthetic test above — a pure
Oseen vortex — the Rankine fit tops out at R² = 0.87 while Oseen reaches 1.0000, and that gap is
structural rather than a matter of tuning. The literature's use of "Rankine-like" for buoyant
plume cores is defensible as shorthand for *near-solid-body rotation inside, decaying
circulation outside*, but the corner itself is fiction.

## Results

Production run: `Ra = 1e7`, `dx = 15 mm`, 133×213 lattice, 91,840 steps covering 90.0 s of room
time in 86 minutes, no turbulence closure, no divergence. The first 36 s are spin-up and are
excluded from every number below.

| Quantity | Simulated | Measured | Note |
|---|---:|---:|---|
| Shedding period, 0.5 m above the head | 3.6 s | ≈5 s | peak-count cross-check gives 3.71 s |
| Strouhal number `f L / U_ff` | 0.62 | 0.44 | the comparable form; ~40% fast |
| Slow meander of the whole column | 9.0 s | — | separate mode, not the shedding |
| Peak plume speed at the probe | 0.41 m/s | 0.2–0.3 m/s | time-mean of the row maximum |
| Per-segment conductance, 11 segments | r = 0.79 | — | Spearman 0.63 |
| … with the foot excluded | r = 0.48 | — | one point carries much of it |
| Foot / chest convective ratio | 30× | 1.7× | ordering right, magnitude wrong |
| Mean fitted core radius | 0.11 m | — | no direct measured counterpart |

### Rankine versus Oseen, over 695 fitted cores

| | Rankine | Oseen |
|---|---:|---:|
| mean R² | 0.737 | 0.774 |
| median R² | 0.803 | 0.798 |
| wins head-to-head | **61.2%** | 37.6% |

Only **9.9%** of cores reach R² = 0.9 under *either* model, and the best single fit is 0.98.

This inverts the clean-theory expectation and is the study's most interesting result. On a
synthetic Oseen vortex the fitter is decisive — Oseen 1.000, Rankine 0.866 — so the corner really
is fiction for an isolated viscous vortex. But on cores the plume actually sheds, Rankine is the
better fit on nearly two thirds of them. Oseen keeps the higher *mean* only because it fails less
badly on the worst cores; the medians are effectively tied. Neither is the shape of a shed core,
because a shed core is being stretched, tilted and advected while it is measured. The literature's
"Rankine-like" comes out of this better than the textbook argument predicts.

### Per-segment ordering

Simulated: foot > hand > head > lower leg > upper arm > forearm > neck > thigh > chest > abdomen > pelvis
Measured: foot > hand > thigh > lower leg > forearm > head > neck > pelvis > chest > abdomen > upper arm

The extremes agree — foot highest, chest and abdomen lowest, which is the counterintuitive part —
while the middle of the ranking scrambles. The pelvis comes out slightly *negative* in the
simulation, meaning the air arriving there is warmer than the skin: air heated by the 32.7 °C shin
reaches a 32.2 °C pelvis. That is physically correct for a 2D slice with no fresh entrainment
between segments, and it is one reason the thigh is badly misplaced.

## What the study earns, and what it does not

| Claim | Status | Why |
|---|---|---|
| Structure along the body: thin layer at the ankles, thickening, going unsteady, detaching above the head | **Reproduced** | Sequence and ordering match schlieren and PIV imaging |
| Discrete vortex shedding above the head | **Reproduced** | St = 0.62 against a measured 0.44; two independent estimators agree to 3% |
| Per-segment convective ordering, feet above chest | **Partly** | Extremes right, middle scrambles, correlation substantially foot-driven, contrast exaggerated 18× |
| On-body laminar→turbulent transition at 0.80 m | **Not reproduced** | Needs local Ra ≈ 3.8×10⁸; the run is 488× below full scale and never reaches it |
| Absolute boundary-layer thickness in cm | **Not claimed** | δ scales as Ra^(−1/4); at reduced Ra the layer is roughly 4.7× too thick, and the figure shows both the raw and the rescaled curve |
| Which vortex idealisation is right | **Neither** | Only 10% of 695 cores reach R² = 0.9 under either |

Two further caveats that apply to every figure here:

- **Two dimensions flatter vortices.** 2D turbulence has no vortex stretching, so cores stay
  coherent longer and read more cleanly than in a real room. Their persistence is not a
  measurement.
- **The silhouette is anterior.** The width the flow sees at the torso is shoulder width rather
  than chest depth, which overstates the torso's characteristic length. It was chosen because
  the narrow ankles and broad trunk are exactly what drives the h_c contrast being tested.
- **The browser exhibit traces a real scan; the study's solver does not.** `plume_lbm.py` uses
  the parametric silhouette in `body_model.py` (arms down, legs apart). The exhibit uses the
  photographed outline (arms out, legs together, 0.93 m across). They are different geometries
  and their numbers are not interchangeable — the exhibit is illustrative, the study is the
  measurement.

## Note on the scan image

`Thermal-Infrared-Full-Body-Scan.jpg` was supplied for this project. `build_body_asset.py`
isolates the subject as the largest connected non-black component, which also removes the
"MEDICAL INFRARED" watermark sitting on the background. **Confirm the image is cleared for
publication before deploying the site** — the study and the solver do not depend on it, and the
exhibit falls back to a shaded mask when the sprite cannot be loaded.

Its palette is uncalibrated false colour for one individual and is never read as data: it
supplies shape and appearance only, while wall temperature comes from the measured population
values in `body_model.SEGMENTS`. In this pose the upper arms are contiguous with the torso above
the armpit, so those cells take the chest value rather than the upper-arm value.

## Sources

- Regional skin surface temperature, 30 subjects, 24.4 °C ambient, 46% RH, 15 exposed sites,
  supine; ten distinct thermal zones, 29.8–35.0 °C — <https://pmc.ncbi.nlm.nih.gov/articles/PMC7992731/>
- Whole-body thermal mapping (Chudecka et al.): trunk warmest, distal lower limb coolest,
  right-side dominance, Y-shaped anterior cool zone.
- 16-segment articulated thermal manikin, boundary-layer wind tunnel: per-segment h_r and h_c,
  standing, still air.
- Human thermal plume structure and velocity — <https://pmc.ncbi.nlm.nih.gov/articles/PMC7817963/>
- Time-resolved PIV above a manikin head: intermittent shedding at roughly 5 s intervals,
  isolated with vorticity, Q and λ_ci.
- Exhaled breath interacting with the plume: counter-rotating vortex pairs —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC7118971/>
- Rankine-like tangential profiles in swirling buoyant plumes — <https://www.mdpi.com/2076-3417/10/10/3428>
- Rankine vortex definition — <https://glossary.ametsoc.org/wiki/rankine-vortex/>

The browser exhibit built on this study is `web/human-thermal-plume.html`.
