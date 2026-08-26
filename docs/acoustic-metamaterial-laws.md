# Acoustic Metamaterial & Passive Hologram — Foundational Laws

**Status:** living document. Rev 4 — 2026-08-23.
**Scope:** governing relations for passive (unpowered) acoustic structures, reduced to closed-form
laws and bound to this shop's actual hardware (Bambu H2S, 0.4 mm nozzle, PETG, 320 x 320 x 320 mm
envelope).

**Two branches, one physics.** Both are unpowered printed geometry acting on an incident field;
they differ only in the **source**, and every divergence below traces back to that.

| | **Ultrasonic branch** (Parts I–X) | **Audible branch** (Parts XI–XVI) |
|---|---|---|
| Source | 40 kHz transducer (air) / MHz emitter (water) | a raised human voice, ~100 Hz – 8 kHz |
| Coherence | single frequency, stable phase front | broadband, ~6 octaves, no phase reference |
| `lambda_air` | 8.575 mm @ 40 kHz | 343 mm @ 1 kHz — **40x larger** |
| Governing limit | index contrast and wrap geometry | conservation of power and Q vs. intelligibility |
| Output | shaped fields, traps, radiation force | directivity, tonal shaping, projection |

Shared across both: horn/impedance matching (Part II), space-coiling (Part V), and every
printability floor in Part IV. Where the audible branch reuses a relation, the ultrasonic Part is
authoritative and the audible Part cites it.

**How to extend:** add new laws under the numbered Part that owns the physics, register any new
material or hardware value in Part 0, and append to the Revision Log. Do not silently change a
`[PARAM]` — supersede it with a dated row so downstream numbers can be re-derived.

## Legend

| Tag | Meaning |
|---|---|
| `[LAW]` | Closed-form relation derived here. Holds for any parameter set. |
| `[PARAM]` | A number. Always carries provenance: `given` / `assumed` / `derived` / `measured`. |
| `[RESULT]` | A numeric conclusion for a specific parameter set. Recompute if a `[PARAM]` moves. |
| `[DEAD]` | Route investigated and closed, with the mechanism that closes it. Recorded so it is not re-tried. |
| `[OPEN]` | Unresolved. Ranked in Part XVII. |

## Contents

**Front matter**

- [Part 0 — Parameter Registry](#part-0--parameter-registry)

**Ultrasonic branch — coherent single-frequency source**

- [Part I — The Master Relation](#part-i--the-master-relation)
- [Part II — Impedance, Transmission, and Etalon Ripple](#part-ii--impedance-transmission-and-etalon-ripple)
- [Part III — Dispersion and Bandwidth](#part-iii--dispersion-and-bandwidth)
- [Part IV — Sampling, Quantization, and Printability](#part-iv--sampling-quantization-and-printability)
- [Part V — Space-Coiling Metasurfaces](#part-v--space-coiling-metasurfaces)
- [Part VI — Field Topologies and Trapping](#part-vi--field-topologies-and-trapping)
- [Part VII — Phase Retrieval](#part-vii--phase-retrieval)
- [Part VIII — Radiation Pressure and Asymmetric Patches](#part-viii--radiation-pressure-and-asymmetric-patches)
- [Part IX — Material Routes Closed](#part-ix--material-routes-closed)
- [Part X — Design Envelope for This Shop](#part-x--design-envelope-for-this-shop)

**Audible branch — incoherent broadband source (the human voice)**

- [Part XI — Scope Decision: What the Pivot Discards](#part-xi--scope-decision-what-the-pivot-discards)
- [Part XII — Output Taxonomy Across Both Branches](#part-xii--output-taxonomy-across-both-branches)
- [Part XIII — The Governing Constraint: No Passive Gain](#part-xiii--the-governing-constraint-no-passive-gain)
- [Part XIV — Architecture of a Voice Projector](#part-xiv--architecture-of-a-voice-projector)
- [Part XV — Speech-Band Tuning and the Q Budget](#part-xv--speech-band-tuning-and-the-q-budget)
- [Part XVI — Build and Measurement Sequence](#part-xvi--build-and-measurement-sequence)

**Back matter**

- [Part XVII — Open Questions, Ranked](#part-xvii--open-questions-ranked)
- [Part XVIII — Source Lineage and Unverified Claims](#part-xviii--source-lineage-and-unverified-claims)
- [Revision Log](#revision-log)

---

## Part 0 — Parameter Registry

### 0.1 Media

| Symbol | Value | Provenance |
|---|---|---|
| `c_air` | 343 m/s | standard, 20 C |
| `rho_air` | 1.2 kg/m^3 | standard |
| `Z_air` | 413 Rayl | derived, `Z = rho c` |
| `c_water` | 1500 m/s | standard |
| `rho_water` | 1000 kg/m^3 | standard |
| `Z_water` | **1.5 MRayl** | given (source material) |
| `mu_air` | 1.81e-5 Pa s | standard |
| `gamma` | 1.4 | standard |
| `Pr` | 0.71 | standard |

Reference wavelengths, ultrasonic branch: `lambda_air(40 kHz) = 8.575 mm`,
`lambda_water(1 MHz) = 1.50 mm`, `lambda_water(2 MHz) = 0.75 mm`.

Reference wavelengths, audible branch (`lambda = c_air / f`):

| f | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz | 4 kHz | 8 kHz |
|---|---|---|---|---|---|---|---|
| `lambda_air` | 2744 mm | 1372 mm | 686 mm | 343 mm | 172 mm | 86 mm | 43 mm |

> **The 320 mm bed is a wavelength constraint, not just a size constraint.** In the ultrasonic
> branch a 100 mm plate is tens of wavelengths across; in the audible branch the largest printable
> aperture is `0.9 lambda` at 1 kHz and `0.23 lambda` at 250 Hz. This single row governs Part XIV.

### 0.2 Plate materials

| Material | rho (kg/m^3) | c_L (m/s) | Z (MRayl) | Provenance |
|---|---|---|---|---|
| Photopolymer resin (SLA) | ~1150 *assumed* | 2174 *derived* | **2.50** *given* | Source text gives Z only; c back-solved as `Z/rho`. |
| **PETG (filament bulk)** | **1230–1330, ~1270** *given* | **UNKNOWN** | derived once c is known | This is the shop material. |
| PETG (printed solid) | ~0.95–0.99 x bulk | lower than bulk | — | Inter-bead voids + moisture. See §4.6. |

> **Critical asymmetry.** For resin the source supplies `Z`, so `c = Z/rho`. For PETG we have
> `rho` but not `c`, so `Z = rho c` and **`c_L` is the free variable that drives everything.**
> Working assumption throughout: **`c_L = 2200 m/s`** → `Z_PETG = 2.794 MRayl`. Every `[RESULT]`
> tagged *(c2200)* inherits this and must be recomputed on measurement. See §17.1. It is the free
> variable for **both** branches — every printed-structure number depends on it.

`[PARAM]` *given (measured externally — Zvoníček et al. 2023, §18.1), audible band:* FDM PET-G,
0.4 mm nozzle, impedance tube 200–3200 Hz. Mean sound reflection `beta_m` = **0.963** for
two-perimeter shell over 40% Grid infill at 0.3 mm layers (best structure tested); **0.852** for
100% solid at 0.2 mm; 0.632–0.878 for Cubic infill (erratic — avoid in reflecting bodies).
Surface roughness across 0.1–0.5 mm layer heights: **Sa 9.7–13.9 μm, Sz 76–111 μm** (flat
surfaces; vertical layer-ridge walls unmeasured). Normal incidence on flat discs — anchors
material/structure choices, not device transfer functions.

### 0.3 Hardware

| Parameter | Value | Provenance |
|---|---|---|
| Printer | Bambu H2S, FDM/FFF | given |
| Nozzle | **0.4 mm** | given |
| Lateral feature floor `d` | 0.4 mm (one extrusion) | derived from nozzle |
| Layer height range | 0.08–0.28 mm typical for 0.4 mm nozzle | standard practice |
| Build envelope | 320 x 320 x 320 mm | repo constraint (CLAUDE.md) |
| Material | PETG | repo constraint |

### 0.4 The voice as a source

Registry for the audible branch. All values are textbook speech science (lineage in §18.1); they
are population means, not measurements of any one talker — see the caveat in §18.3.

| Parameter | Value | Provenance |
|---|---|---|
| `f0` spoken, adult male | 85–180 Hz | given, standard |
| `f0` spoken, adult female | 165–255 Hz | given, standard |
| `f0` sung, operatic | 82 Hz (bass low E2) to 1047 Hz (soprano high C6); 1397 Hz (F6, Queen of the Night) | given, standard |
| Vocal tract length | ~170 mm male, ~145 mm female | given, standard |
| Epilaryngeal tube | 25–35 mm effective; decouples as its own resonator when its area is <= ~1/6 of the pharynx | given (Sundberg) |
| Singer's formant (*squillo*) | 2.4–3.2 kHz cluster; lower for basses, higher for tenors | given (Sundberg) |
| Vibrato | rate 5–7 Hz; extent ±0.5 to ±1 semitone = ±3–6% | given, standard |
| Stop-consonant burst | 5–20 ms | given, standard |
| Syllable rate, running speech | 4–7 per second | given, standard |
| Ear canal | ~25 mm closed pipe → quarter-wave near 3.4 kHz | given, standard |
| Glottal source resistance | 30–100 cgs acoustic ohm = 3–10 MPa·s/m³ during phonation (the voice is a stiff but finite flow source) | given, coarse |

Vowel formants, adult male means (Peterson & Barney 1952) — consumed by §15.10:

| vowel | F1 (Hz) | F2 (Hz) |
|---|---|---|
| ee (heed) | 270 | 2290 |
| ih (hid) | 390 | 1990 |
| eh (head) | 530 | 1840 |
| ae (had) | 660 | 1720 |
| ah (father) | 730 | 1090 |
| aw (bought) | 570 | 840 |
| oo (boot) | 300 | 870 |

> Note three lengths in this table — epilarynx 25–35 mm, ear canal ~25 mm — and hold them next to
> Law 14.4's 1–4 kHz mouth passband. §15.11 shows the coincidence is the design.

---

## Part I — The Master Relation

### 1.1 Phase from thickness

`[LAW 1.1]`

```
Delta_phi(x,y) = 2 pi f ( 1/c_p - 1/c_m ) h(x,y)
```

Assumes each column propagates as an independent plane wave. Validity is bounded by §4.5.

### 1.2 The 2-pi thickness is a pure material property

`[LAW 1.2]` — set `|Delta_phi| = 2 pi`, normalize by `lambda_m = c_m/f`. Frequency cancels:

```
  h_2pi          1                              c_m
  ------  =  ---------- ,      with   n  =  ---------
  lambda_m    | n - 1 |                        c_p
```

`n` is the acoustic refractive index of the plate relative to the medium. This is the optical
phase-plate result. **Required thickness in wavelengths is frequency-independent.**

`[RESULT]` *(c2200 for PETG)*

| medium / plate | n | h_2pi / lambda | h_2pi @ 40 kHz air | h_2pi @ 1 MHz water | h_2pi @ 2 MHz water |
|---|---|---|---|---|---|
| air / resin | 0.158 | 1.19 | 10.2 mm | — | — |
| air / PETG | 0.156 | 1.18 | 10.2 mm | — | — |
| water / resin | 0.690 | 3.23 | — | 4.84 mm | 2.42 mm |
| water / PETG | 0.682 | 3.14 | — | **4.71 mm** | 2.36 mm |

> **Water holograms are ~2.7x thicker in wavelengths than air ones would be**, purely because
> water's index contrast against thermoplastics is weak (n ~ 0.68 vs 0.16). Low contrast is the
> single root cause of plate thickness, which in turn is the root cause of the column-approximation
> breakdown (§4.5) and the etalon ripple count (§2.3). These are not three problems.

### 1.3 Sensitivity of thickness to plate speed

`[LAW 1.3]` — elasticity of `h_2pi` with respect to `c_p`:

```
  d(h_2pi)/h_2pi              n
  --------------  =  -  -------------
    d(c_p)/c_p             1 - n
```

`[RESULT]`

| n | elasticity | interpretation |
|---|---|---|
| 0.158 (air/PETG) | −0.19 | 1% speed error → 0.2% thickness error. Immune. |
| 0.652 | −1.87 | |
| 0.682 (water/PETG, c2200) | **−2.15** | 1% speed error → 2.1% thickness error. |
| 0.714 | −2.50 | |
| 0.750 | −3.00 | |

The `1/(1-n)` divergence means water-immersed designs are hypersensitive to plate speed and air
designs are not. **This law is why §17.1 is the top-ranked open item.**

### 1.4 PETG sensitivity sweep (water)

`[RESULT]` — `rho = 1270`, sweeping the unknown `c_L`:

```
 c_L    Z=rho*c   n=1500/c   h_2pi/lam   h_2pi@1MHz   r=Z/1.5   T_min   p_min   ripple cyc/wrap
2000     2.540      0.750       4.00       6.00 mm     1.693    0.767   0.876       6.00
2100     2.667      0.714       3.50       5.25 mm     1.778    0.730   0.854       5.00
2200     2.794      0.682       3.14       4.71 mm     1.863    0.695   0.833       4.29
2300     2.921      0.652       2.88       4.31 mm     1.947    0.661   0.813       3.75
```

Sweeping `rho` across its full 1230–1330 band at fixed `c` moves `h_2pi` **not at all** (Law 1.2
depends only on `c`) and moves `T_min` by ~3 points. **Speed dominates; density is nearly inert.**

---

## Part II — Impedance, Transmission, and Etalon Ripple

### 2.1 Single interface

`[LAW 2.1]`

```
              4 Z1 Z2
  T  =  -----------------      (intensity)
          ( Z1 + Z2 )^2
```

`[RESULT]`

| interface | T | dB / interface | dB through plate (x2) |
|---|---|---|---|
| resin ↔ water | 0.9375 | −0.28 | −0.56 |
| PETG ↔ water *(c2200)* | 0.9092 | −0.41 | **−0.83** |
| resin ↔ air | 6.61e-4 | −31.8 | −63.6 |
| **PETG ↔ air** *(c2200)* | 5.91e-4 | −32.3 | **−64.6** |

### 2.2 Air is closed to monolithic plates

`[DEAD]` A solid monolithic plate in air is ~64 dB opaque. The half-wave resonance escape
(`T = 1` at `h = m lambda_p / 2`) does not rescue it. Half-power occurs at
`sin(k_p h) = 2/(r - 1/r)` with `r = Z_p/Z_air = 6765` *(PETG, c2200)*:

```
  sin(k_p h) = 2.96e-4
  k_p        = 2 pi f / c_p = 114.2 rad/m       (40 kHz)
  Delta_h    < 2.96e-4 / 114.2  =  2.6 um       [thickness tolerance]
  fractional bandwidth ~ 1e-4
```

Both are orders of magnitude beyond any printer.

> **Consequence — this reframes space-coiling entirely.** Coiling in air is *not* a thinner
> alternative to a monolithic hologram; it is **the only available mechanism**, because it keeps
> the propagating medium as **air inside the channels** and uses geometry for delay rather than
> bulk transmission through solid. This is why the airborne literature describes "metamaterial
> delay-line structures" and why every printed monolithic hologram in the source material is
> **water-immersed**. The impedance table predicts the split in the literature.

### 2.3 Etalon ripple — amplitude modulation locked to the phase map

Each column is a different-thickness Fabry-Perot etalon.

`[LAW 2.3a]`

```
                        1                                     Z_p
  T(h)  =  -----------------------------  ,      r  =  ---------
            1 + (1/4)(r - 1/r)^2 sin^2(k_p h)              Z_m
```

`[LAW 2.3b]` — ripple cycles per full 2-pi phase wrap:

```
   h_2pi           2 c_m
  --------  =  -------------
  lambda_p/2     c_p - c_m
```

`[RESULT]` *(PETG, c2200, water)* — `r = 1.863`:

```
  T   in [0.695, 1.000]   intensity
  p   in [0.833, 1.000]   pressure amplitude      = +/- 9% modulation
  4.29 ripple cycles per phase wrap
```

For resin (`r = 1.667`) it is `T in [0.779, 1.000]`, 4.45 cycles. **PETG's higher impedance makes
the etalon problem measurably worse than the resin baseline in the source material.**

This is not noise. It is a deterministic modulation spatially locked to the thickness map, on a
period finer than the zone spacing, and it therefore feeds *coherently* into sidelobes.

### 2.4 The two knobs are separable

`[LAW 2.4]` Since `Z = rho c`:

- **Ripple depth** is set by impedance contrast `Z_p/Z_m`.
- **Ripple count and plate thickness** are set by speed contrast `c_p/c_m`.

They can only be decoupled through **density**. This motivated the low-density routes in Part IX
(all closed) and, failing those, the computational route in §7.3.

---

## Part III — Dispersion and Bandwidth

### 3.1 An unwrapped monolithic plate is achromatic

`[LAW 3.1]` A monolithic plate imposes a **true time delay** `tau(x,y) = K h(x,y)` with
`K = |1/c_p - 1/c_m|`, independent of `f`.

Proof of achromaticity: the phase required to focus at `F` is
`phi_req(r) = -(2 pi f / c_m)(sqrt(r^2 + F^2) - F)`, which scales as `f`. The delivered phase
`2 pi f K h` also scales as `f`. Setting delivered = required at `f0(1+eps)` gives
`sqrt(r^2 + F^2) - F = sqrt(r^2 + F0^2) - F0`, satisfied by `F = F0` for all `eps`.

> **There is no chromatic focal shift in an unwrapped monolithic hologram.** Delay-based
> beamforming is inherently broadband. This is the opposite of the usual intuition about
> "narrowband metamaterials," and it identifies the real culprit below.

### 3.2 Wrapping destroys it — the Fresnel penalty

`[LAW 3.2]` Cliffs are exactly `2 pi` only at `f0`. At `f = f0(1+eps)` each cliff errs by
`2 pi eps`, accumulating across Fresnel zones. Outermost-zone error is `2 pi N_z eps`. Requiring
`<= lambda/4` wavefront error:

```
  Delta_f          1
  -------   <=   -------
    f             4 N_z
```

with zone count

```
                sqrt( (D/2)^2 + F^2 )  -  F
  N_z   =   -----------------------------------
                        lambda_m
```

`[RESULT]` D = 50 mm, F = 50 mm:

| frequency | lambda | N_z | Delta_f/f | unwrapped slab thickness |
|---|---|---|---|---|
| water 1 MHz | 1.50 mm | **3.93** | **6.4%** | 18.5 mm *(c2200)* |
| water 2 MHz | 0.75 mm | 7.87 | 3.2% | 18.5 mm |

> **Thickness buys bandwidth one-for-one with Fresnel zone count.** A wrapped plate is a Fresnel
> lens and inherits Fresnel-lens chromatics; the thick monolithic hologram is the achromatic one.
> Note the unwrapped slab thickness is identical at both frequencies — it is `N_z * h_2pi`, and
> both factors scale inversely with lambda.

### 3.3 Coiled cells are strictly worse

`[LAW 3.3]` Coiled cells are wrapped by construction (each cell is mod 2 pi), so they inherit the
`1/(4 N_z)` ceiling **and add resonance dispersion**: each cell has a different path length `L`,
hence a Fabry-Perot resonance at a different frequency (spacing `c/2L`), hence a different
`d(phi)/df` and a different transmitted amplitude off design. The aperture develops
**cell-dependent amplitude speckle and non-uniform phase drift** — the map is corrupted, not
merely shifted.

Ranking:

```
  unwrapped monolithic  >  wrapped monolithic  >  coiled
      (broadband)           (1 / 4 N_z)          (1 / 4 N_z, minus resonance)
```

Horn-like coiling structures matter here specifically because they make amplitude a *designed*
variable rather than a resonance artifact.

---

## Part IV — Sampling, Quantization, and Printability

### 4.1 Nozzle diameter is a Nyquist limit

`[LAW 4.1]` Grating-lobe-free encoding requires pixel pitch `d <= lambda_m / 2`, hence

```
                  c_m
  f_Nyquist  =  -------
                 2 d
```

`[RESULT]` `d = 0.4 mm`:

```
  water:   f_Nyquist  =  1500 / (2 * 0.4e-3)  =  1.875 MHz
  air:     f_Nyquist  =   343 / (2 * 0.4e-3)  =    429 kHz
```

| medium / freq | lambda/2 | d/lambda | margin | verdict |
|---|---|---|---|---|
| air 40 kHz | 4.29 mm | 1/21 | 10.7x | fine |
| water 1 MHz | 0.75 mm | 1/3.75 | **1.88x** | fine |
| water 2 MHz | 0.375 mm | 1/1.88 | 0.94x | aliased |
| water 5 MHz | 0.15 mm | 1/0.75 | 0.38x | dead |

> **A 0.4 mm nozzle caps water-immersed monolithic holograms just under 1.9 MHz.** Air at 40 kHz
> has 10x headroom.

### 4.2 The historical "100x" resolution claim, recomputed

Source claim: replacing a transducer array with a printed plate improved resolution ~100x because
the aperture element dropped from transducer diameter to printer feature size. Recovered exactly
by `10 mm pitch / 0.1 mm feature` (SLA-class printer). **For an H2S: `10 / 0.4 = 25x`.**

`[LAW 4.2]` This does not weaken the conclusion. Lateral resolution is a **threshold**
(`d <= lambda/2`), not a gradient — focal spot is aperture-limited at `w ~ lambda F / D`
regardless of pixel count. Once Nyquist is cleared, extra pixels buy nothing.

Concretely, in wavelength units at 40 kHz in air:

```
  array   d = 10 mm  = 1.17 lambda  ->  grating lobe at sin(theta) = lambda/d = 0.857, theta = 59 deg
  plate   d = 0.4 mm = lambda/21    ->  grating-lobe free at all angles
```

The array is spatially *undersampled*; the plate is oversampled. **You do not need a resin printer
for this work. You need to stay under 1.9 MHz.**

### 4.3 Axial quantization — N-step blazed grating, not Strehl

`[LAW 4.3]` Layer stepping on a blazed profile is a staircase, so use grating-order efficiency:

```
          [ sin(pi/N) ]^2                       h_2pi
  eta  =  [ --------- ]      ,     N  =  ------------------
          [   pi/N    ]                    layer height
```

`[RESULT]` *(c2200)*

| | h_2pi | 0.08 mm layers | 0.20 mm layers |
|---|---|---|---|
| water 1 MHz | 4.71 mm | N=59 → 99.9% | N=24 → **99.4%** |
| water 2 MHz | 2.36 mm | N=30 → 99.6% | N=12 → 97.7% |

**Layer height is not the limiting error term anywhere in the usable band.** The lost few percent
goes into spurious diffraction orders rather than heat, but it is second-order. Print at
0.16–0.20 mm and spend effort elsewhere.

### 4.4 Zone width, numerical aperture, and the printability bound

`[LAW 4.4a]` Rim zone width for a focusing plate:

```
                lambda_m                                    r
  Delta_r  =  ----------- ,       NA  =  sin(theta)  =  -------------------
                  NA                                    sqrt(r^2 + F^2)
```

`[LAW 4.4b]` Requiring at least 4 extrusion widths per zone (`Delta_r >= 4 d`):

```
              lambda_m
  NA   <=   ------------
               4 d
```

`[RESULT]` `d = 0.4 mm`:

```
  water 1 MHz  ->  NA <= 0.94    (effectively unconstrained)
  water 2 MHz  ->  NA <= 0.47    (binding; D=50/F=50 sits at NA = 0.447, right at the edge)
```

> After the frequency ceiling, the second thing the nozzle costs is **numerical aperture at the
> top of the band.**

### 4.5 Wrap-cliff aspect ratio — the dominant monolithic error

`[LAW 4.5]` Dividing Law 1.2 by Law 4.4a, lambda cancels:

```
  cliff height           NA
  ------------  =  -------------          (frequency-independent)
   zone width        | n - 1 |
```

`[RESULT]`

| case | NA | n | aspect |
|---|---|---|---|
| water / PETG *(c2200)* | 0.447 | 0.682 | **1.41 : 1** — walls taller than wide |
| air / PETG | 0.447 | 0.156 | 0.53 : 1 — walls wider than tall |

At 1 MHz that is a 4.71 mm wall every 3.36 mm, i.e. 3.1 wavelengths tall with 2.2-wavelength
spacing. **Shadowing and inter-column diffraction are first-order, not a correction term.** This
invalidates Law 1.1's independent-column assumption exactly where monolithic plates are most
useful — and note it is created by the mod-2-pi step, not by the smooth ramp between cliffs.

Once again the root cause is low index contrast: the air row shows the problem simply does not
exist at `n = 0.16`.

### 4.6 Print orientation and printed-vs-bulk density

**Orientation.** Print **flat face down, structured face up**. Cliffs become vertical walls (FDM's
best case, zero overhang); ramps become upward-facing slopes. At 1 MHz the ramp is 35 deg from
vertical, inside the 45 deg limit even if flipped. No support anywhere.

**Density.** The registered 1270 kg/m^3 is *filament bulk*. Printed solid differs:

- **Inter-bead voids.** "100% infill" FDM typically lands at 95–99% of bulk. At 97%: `rho -> 1232`,
  and stiffness drops faster than density (~6–8% off E for 3% porosity), netting **c down ~2%**.
  Through Law 1.3 (elasticity −2.15) that is a **~4.5% thickness error** — 0.21 mm on a 4.71 mm
  wrap, about one layer, and systematic.
- **Moisture.** PETG is hygroscopic; wet filament produces voids and poor interlayer bonding —
  the same failure mode, amplified and made non-uniform. Dry the spool.

> Both effects argue the same thing: **measure `c_L` on a printed block from the production
> profile, never from a datasheet.** The measurement captures voids, layer adhesion, and
> orientation anisotropy in one number, and it calibrates out the systematic offset.

---

## Part V — Space-Coiling Metasurfaces

### 5.1 Coiling factor

`[LAW 5.1]` An air-filled folded channel of physical thickness `t` and path length `L`, against a
free reference path `t`:

```
  Delta_phi  =  2 pi (L - t) / lambda        ->     full 2 pi requires   L = t + lambda

                       L               lambda
  n_eff       =       ---     =   1 + --------
                       t                  t
```

`[RESULT]` at 40 kHz in air (`lambda = 8.575 mm`):

| t | n_eff | L | passes at a = lambda/2 |
|---|---|---|---|
| lambda/4 = 2.14 mm | 5.0 | 10.72 mm | ~3 |
| lambda/2 = 4.29 mm | 3.0 | 12.86 mm | ~3 |

The "under a quarter-wavelength thick" target in the literature is exactly a **5x coiling factor**,
in a unit cell that must also stay `<= lambda/2 = 4.29 mm` laterally (Law 4.1).

### 5.2 Wall thickness budget — the 0.4 mm nozzle penalty

With 3 serpentine passes and 0.4 mm dividers, `3h + 2(0.4) = t`:

```
  t = lambda/4 = 2.14 mm   ->   h = 0.448 mm
  t = lambda/2 = 4.29 mm   ->   h = 1.163 mm
```

**0.4 mm walls consume 37% of a quarter-wavelength cell's thickness at 40 kHz.**

### 5.3 Viscothermal loss

`[LAW 5.3]`

```
                                        omega     delta_v  [      gamma - 1  ]
  delta_v = sqrt(2 mu / rho omega),   alpha = ----- * ------- * [ 1 + --------- ]
                                          c        2 h     [      sqrt(Pr)  ]
```

`[RESULT]` air, 40 kHz — `delta_v = 11 um`:

| t | channel h | alpha | alpha L | loss |
|---|---|---|---|---|
| lambda/4 | 0.448 mm | 13.2 Np/m | 0.142 Np | **1.23 dB** |
| lambda/2 | 1.163 mm | 5.09 Np/m | 0.065 Np | **0.57 dB** |

> **Half the loss for twice the thickness.** Since air holograms have no impedance or absorption
> problem, `t = lambda/2` is the better operating point on this machine unless compactness is the
> actual goal. Note `alpha ~ 1/h`: aggressive coiling is paid for in amplitude, linearly.

`[OPEN]` FDM channel walls run ~50–100 um Ra against an 11 um viscous boundary layer, so roughness
dominates the wall interaction. **Treat 1.23 dB as a floor; real loss is plausibly 1.5–3x higher.**

### 5.4 Closed form, valid in both branches

`[LAW 5.4]` Substituting `delta_v = sqrt(2 mu / rho omega)` into Law 5.3 and collapsing constants:

```
                1.474         2 mu            sqrt(omega)
  alpha  =  ----------- sqrt( ------ )  =  K --------------  ,     K = 1.180e-5 (air, SI)
               2 h c           rho                  h
```

```
  alpha  ~  sqrt(f) / h
```

Two consequences that matter more in the audible branch than the ultrasonic one:

1. **Loss per metre falls as frequency falls**, but audible coiled paths are ~50x longer
   (0.5 m vs 10 mm), so total `alpha L` can still be larger.
2. **The loss is a `sqrt(f)` tilt, not a flat insertion loss.** A coiled path is a low-pass filter
   whose slope is set entirely by channel height.

`[RESULT] 5.5` Audible branch, 0.5 m coiled path (verifies against Part 5.3: at 40 kHz,
`h = 0.448 mm`, Law 5.4 returns 13.2 Np/m — exact match):

| channel h | 1 kHz | 4 kHz | tilt across 1–4 kHz |
|---|---|---|---|
| 1 mm | 4.06 dB | 8.13 dB | **−4.1 dB** |
| 3 mm | 1.35 dB | 2.71 dB | −1.4 dB |
| 5 mm | 0.81 dB | 1.63 dB | −0.8 dB |
| 10 mm | 0.41 dB | 0.81 dB | −0.4 dB |

> **This is the quantitative reason coiling is dangerous for speech.** The tilt lands directly on
> the 1–4 kHz intelligibility band (§15.6). Narrow channels do not merely lose level — they
> selectively remove the consonants. **Keep audible-branch channels at or above ~5 mm**, and treat
> the roughness caveat above as making these figures optimistic.

---

## Part VI — Field Topologies and Trapping

### 6.1 The vortex parity argument

Vortex beams carry OAM via `exp(i l theta)`; distinct `l` are orthogonal since
`integral exp(i(l1-l2)theta) dtheta = 0`, which is what enables OAM multiplexing.

`[LAW 6.1]` Under `l -> -l`:

```
  trapping force    F   ~  grad |p|^2        ->   EVEN in sign of l
  radiation torque  tau ~  l P_abs / omega   ->   ODD  in sign of l
```

Torque spins the trapped particle up until centrifugal demand exceeds the radial trap:

```
  m Omega^2 R  >  F_trap,max        =>   ejection
```

**Virtual vortices** — pulse sequences of equal helicity, alternating chirality — time-average the
odd term to zero while leaving the even term at full strength. The decoupling of trap force from
rotation rate is **parity cancellation**, not a tuned effect.

### 6.2 Static plates cannot do it, and the obvious workaround fails

`[DEAD]` A frozen plate has one `l` and a nonzero time-averaged torque. Spatial superposition is
not a substitute for temporal alternation:

```
  exp(+i l theta) + exp(-i l theta)  =  2 cos(l theta)
```

That is a `2l`-lobed **standing** pattern with zero OAM — the vortex ring is destroyed and replaced
by discrete lobes. There is no static encoding of alternating chirality.

### 6.3 Resulting particle-size bound

`[LAW 6.3]` Mie / wavelength-order means `ka ~ 1`; Gorkov (Rayleigh) validity is `ka < 0.3`.

```
  a_Mie       =  lambda / (2 pi)
  a_Rayleigh  <  0.3 lambda / (2 pi)
```

`[RESULT]`

| | a_Mie | a_Rayleigh |
|---|---|---|
| air 40 kHz | 1.37 mm | < 0.41 mm |
| water 1 MHz | 0.24 mm | < 0.072 mm |
| water 2 MHz | 0.12 mm | < 0.036 mm |

The static-plate vortex trap is confined to roughly **the bottom decade of particle size** a phased
array can handle.

### 6.4 Bottle beams — the escape route

Bottle beams (null-pressure pocket inside a closed high-pressure shell) were first demonstrated
with **no metamaterial at all** — phase engineering across a conventional array in homogeneous
space — showing genuine negative (pulling) force on a rigid ball inside, plus obstacle-bending and
self-healing. Later work ported the topology to phase-modulating holographic lenses, including
**multi-bottle configurations from a single static water-immersed hologram**.

`[LAW 6.4]` Bottle confinement is a **pure spatial-topology** trap: no OAM, therefore no torque,
therefore nothing for the parity argument of §6.1 to break.

> **Design rule: static plates should use bottles, not vortices.** The correct reading of the
> source material is not "static plates cannot hold large particles" — it is that they must reach
> them through the cage topology rather than the helical one.

---

## Part VII — Phase Retrieval

### 7.1 The GS family

```
  GS:  p_n = Prop[ A_n e^{i phi_n} ];  enforce |p| = p_target;  back-propagate;  enforce |A| = const
```

| Method | Relaxation / update | Buys |
|---|---|---|
| **GS** | amplitude constrained on both planes | baseline, monotonic error decrease |
| **GS-PAT** | drops target-point amplitude constraint entirely | large speed gain — the constrained system collapses to the `N_trap x N_trap` subspace rather than the full `M_transducer` problem |
| **WGS** | `w_k <- w_k * (p_bar / abs(p_k))` each iteration | minimizes pressure variance across simultaneous traps |

### 7.2 Unified objective

`[LAW 7.2]` All three minimize

```
  J  =  sum_k  w_k ( abs(p_k) - p_k_target )^2   +   lambda * C(hardware)
```

with different `(w, C)`. This is the "iterative backpropagation" framing.

### 7.3 Why static plates need a different optimizer

`[LAW 7.3]` The two-step route (retrieve `phi` → wrap → `h = phi / 2 pi f K`) **manufactures the
geometry that invalidates its own approximation**: wrapping creates the cliffs of §4.5, while the
smooth ramp between them is gentle. Direct 3D thickness optimization under a differentiable
forward model (the "End-to-End Homogeneous Physics" approach) can trade smoothness against ideal
phase and place discontinuities where they hurt least — a degree of freedom the linear conversion
structurally cannot express.

**The same optimizer absorbs the etalon ripple.** `T(h)` from Law 2.3a is a closed-form
deterministic function of the very thickness map being optimized. It is not noise; it slots
directly into `C(h)`:

```
  J  =  sum_k  w_k ( abs(p_k) - p_k_target )^2   +   lambda * C(h)      <-  C(h) carries T(h)
```

> Both dominant static-plate error terms — wrap-cliff diffraction and etalon ripple — are
> consequences of geometry that the two-step route treats as invisible, and both are addressable
> by the same method. This is the strongest single argument in the document for direct 3D
> thickness optimization.

---

## Part VIII — Radiation Pressure and Asymmetric Patches

### 8.1 Mechanism

Metamaterial patches with sawtooth surface features let a single external speaker exert directional
force on an attached object, by reflecting incident sound differently across asymmetric facets.
Demonstrated pushing, pulling, and rotating floating objects (wood, wax, foam) and fully submerged
objects in 3D.

**Distinguish sharply** from rotating asymmetric propeller blades (a separate mechanical steering
approach for UUVs): that is fluid-dynamic thrust vectoring with moving parts; this is
radiation-pressure asymmetry with none. The real advantage is **zero embedded electronics on the
manipulated object** — a passive patch retrofits existing hardware.

### 8.2 Force scaling

`[LAW 8.2]`

```
  I      =  p0^2 / (2 rho c)
  P_rad  =  2 I / c                        (perfect reflector, normal incidence)
  F      =  P_rad * A
  F_lat  ~  (0.2 - 0.5) * F                (sawtooth asymmetry fraction)
```

`[RESULT]` water, `p0 = 100 kPa`, 30 x 30 mm patch:

```
  I      = 3.33 kW/m^2
  P_rad  = 4.44 Pa
  F      = 4.0 mN
  F_lat  ~ 1 - 2 mN                =  0.1 - 0.2 gram-force
  v_term = sqrt(1.5e-3 / 0.45)     =  ~6 cm/s        (C_d = 1, neutrally buoyant)
```

> **This actuates grams.** Drug-delivery capsules and small assembly components are in range;
> UUV propulsion is not, absent a very large insonified aperture. The source material lists force
> magnitude as unquantified — this is an order-of-magnitude bound, derived here, not sourced.

`[OPEN]` Angular steering resolution remains unbounded; it depends on facet-scale scattering detail
the source does not supply.

---

## Part IX — Material Routes Closed

Recorded so they are not re-attempted. All target the same goal: null the §2.3 ripple by bringing
`Z_p` to `Z_water = 1.5 MRayl`.

### 9.1 Naive target: low-density solid

`[DEAD]` Initial target was `rho ~ 500 kg/m^3` at `c ~ 3000 m/s` (giving `Z = 1.5`, `n = 0.5`,
`h_2pi = 2 lambda`). Correct in isolation, unreachable by foaming — see 9.2.

### 9.2 Bending-dominated (Gibson-Ashby) foam

`[DEAD]` For open-cell foam, `E*/E_s = rho_rel^2`, hence

```
  c*  =  sqrt(E*/rho*)  =  c_s * sqrt(rho_rel)
```

`[RESULT]` PETG *(c2200)* at `rho_rel = 0.40`:

```
  c*     = 1391 m/s
  n      = 1500/1391 = 1.078       <- index INVERTS past 1
  h_2pi  = 12.8 lambda             <- 4x thicker
  Z*     = 0.707 MRayl -> T_min = 0.59    <- ripple WORSE
```

Worse, there is a hard singularity inside the useful porosity range:

`[LAW 9.2]`

```
  c* = c_water   at   rho_rel = ( c_water / c_s )^2  =  0.465

  =>  n = 1 exactly, h_2pi -> infinity, the plate stops working
```

At that point `Z* = 0.886 MRayl` — **the index singularity is reached before the impedance match.**
Foaming toward `Z = 1.5 MRayl` walks straight into it.

### 9.3 Axially-aligned sealed-channel lattice

Structurally the correct idea. Voigt bound for aligned prisms gives `E_axial = E_s * rho_rel`
exactly, so:

```
  c_eff  =  c_s                     (unchanged)
  Z_eff  =  rho_rel * Z_s           (scales linearly)

  Z_eff = 1.5 MRayl   =>   rho_rel = 1.5 / 2.794 = 0.537
```

`[RESULT]` A remarkable coincidence at 1 MHz in water, where the required cell period is
`lambda/2 = 0.75 mm`:

```
  minimum lamellar period  =  0.4 mm wall + 0.35 mm gap  =  0.75 mm  =  lambda/2 exactly
  solid fraction           =  0.4 / 0.75  =  0.533
  Z_eff                    =  0.533 * 2.794  =  1.49 MRayl     vs. water 1.50 MRayl
```

**0.7% from a perfect impedance match, using the minimum feature the nozzle can produce.** The
printer's resolution floor and the acoustic matching condition intersect at 1 MHz in water.

At 2 MHz the required period is 0.375 mm — **below a single extrusion**. The structure cannot exist
on this machine above ~1.9 MHz, which is the same ceiling as Law 4.1.

`[DEAD]` **Bar-mode hazard closes it anyway.** Thin lattice walls have free lateral surfaces, so
the propagating mode is the bar speed, not bulk longitudinal:

`[LAW 9.3]`

```
    c_L                 (1 - nu)
  -------  =  sqrt( ------------------ )
   c_bar             (1+nu)(1-2nu)
```

| nu | ratio | c_bar *(c_L = 2200)* | resulting n |
|---|---|---|---|
| 0.30 | 1.160 | 1896 m/s | 0.791 |
| 0.35 | 1.267 | 1737 m/s | 0.864 |
| 0.40 | 1.464 | **1503 m/s** | **1.00 — at the singularity** |

And the hazard is inescapable: avoiding the lattice's own diffraction orders requires cell period
`< lambda/2`, so walls are necessarily thin compared to lambda. **The impedance-matching lattice
drags `c` toward exactly the value where the hologram collapses**, and where it lands depends on
PETG's Poisson ratio — spanning "workable at 4.8 lambda" to "nonfunctional."

### 9.4 Resolution

With the material routes closed, ripple correction belongs in the End-to-End thickness optimizer
(§7.3), where `T(h)` is a known closed-form term rather than an error to be suppressed.

---

## Part X — Design Envelope for This Shop

Everything above converges on one corner of the space.

```
  medium        water immersion            (air is -64.6 dB closed to monolithic; Law 2.2)
  frequency     1 MHz                      (nozzle Nyquist = 1.875 MHz; 1.88x margin; Law 4.1)
  material      PETG, dried                (see 4.6)
  orientation   flat face down             (zero overhang; see 4.6)
  aperture      D = 50 mm, F = 50 mm       (NA 0.447)
  layers        0.16 - 0.20 mm             (eta >= 99.4%; Law 4.3)
  nozzle        0.4 mm

  ---- derived, all tagged (c2200) - RECOMPUTE ON MEASUREMENT ----
  h_2pi         = 4.71 mm
  N_z           = 3.93 zones
  total relief  = 18.5 mm   (+ substrate => ~23 mm, ~120 layers at 0.20 mm)
  bandwidth     = Delta_f/f <= 6.4%   =>  +/- 64 kHz
  rim zone      = 3.36 mm  =  8.4 extrusion widths
  cliff aspect  = 1.41 : 1            <- dominant error term (Law 4.5)
  ripple        = p in [0.833, 1.000], 4.29 cycles per wrap (Law 2.3)
  transmission  = -0.83 dB from impedance alone (Law 2.1)
```

Trivially inside the 320^3 envelope.

**Airborne variant** (40 kHz, space-coiled, per Part V): `t = lambda/2 = 4.29 mm` cells,
`n_eff = 3`, `a <= 4.29 mm` lateral, 3 serpentine passes with 0.4 mm dividers giving 1.16 mm
channels, ~0.57 dB per cell before roughness correction.

---

## Part XI — Scope Decision: What the Pivot Discards

Parts I–X assume a **powered, coherent, single-frequency source**. Replace it with a person
speaking loudly into an aperture and three things change immediately. They drive every choice in
Parts XII–XVI.

1. **Broadband, not monochromatic.** An ultrasonic plate is designed at one frequency; speech spans
   ~6 octaves with fast transients. **Every narrowband mechanism in Parts I–X becomes a
   *coloration* here rather than a *function*.** Law 3.2's `Delta_f/f <= 1/(4 N_z)` is not a
   tolerance to meet in this branch — at `N_z >= 1` it is already violated across the speech band
   by two orders of magnitude.
2. **No phase reference.** Holography, vortex generation, and trapping all assume a coherent
   single-frequency source. A voice has no stable phase front to sculpt.
   `[DEAD]` for this branch: holograms and kinoforms (Part I), vortex torque (Part VI), bottle and
   twin traps (§6.4), particle sorting, and the entire phase-retrieval apparatus of Part VII.
3. **Wavelengths are 40–800x larger.** `lambda_air(1 kHz) = 343 mm` against `8.575 mm` at 40 kHz.
   Aperture behaviour that was trivially satisfied by a 100 mm ultrasonic plate now collides with
   the 320 mm print envelope. See Law 14.4 — this is the binding constraint of the whole branch.

**What survives the pivot:** horn / impedance matching (the Part II mechanism, applied to radiation
rather than transmission), space-coiling as a *packaging* trick only (Part V, and note Law 5.4's
`sqrt(f)/h` penalty), Helmholtz and quarter-wave resonators as *tuning* elements, reflector and
lens geometry for directivity, and every printability floor in Part IV.

> **Symmetry worth noting.** The ultrasonic branch fails from *too much* geometric detail (wrap
> cliffs, §4.5) and the audible branch fails from *too little* aperture (Law 14.4). Both are the
> same ratio — feature size against wavelength — read from opposite ends.

---

## Part XII — Output Taxonomy Across Both Branches

The map of the field, with each row marked for branch applicability. The unifying chain:

```
printed geometry
  -> local phase / amplitude response
    -> pressure and velocity field
      -> momentum flux
        -> force, torque, or imaging output
```

### 12.1 Wavefield outputs

| Design | Output | Typical result | Audible branch |
|---|---|---|---|
| Thickness-coded hologram / kinoform | Prescribed 3D pressure field | Focuses one source into a point, line, letter, logo, or distributed intensity pattern | `[DEAD]` — needs coherence |
| Spiral phase plate | Acoustic vortex | Hollow-core beam, helical phase, OAM, torque on suitable particles | `[DEAD]` |
| Bottle-beam hologram | 3D low-pressure cage | Null region inside a high-pressure shell; contactless trapping | `[DEAD]` |
| Twin-trap metasurface | Two opposed pressure lobes | Pinch trap; levitate or translate a small object; basis of levitated-particle displays | `[DEAD]` |
| Phase-gradient metasurface | Beam steering / anomalous reflection | Redirects a plane wave off the specular angle | Partial — works on a *narrow band* of the voice spectrum only |
| Space-coiling slab | Large phase delay in a thin panel | Refraction, focusing, bending without bulk thickness | **Retained as packaging** — §14.3 |

### 12.2 Mechanical outputs

Field gradients and momentum flux, not merely "shaped sound." All require ultrasonic intensities;
all are `[DEAD]` at voice SPL. Part VIII quantifies why — `F = 2IA/c` yields ~4 mN at 100 kPa, and
a raised voice at 1 m is roughly six orders of magnitude below that intensity.

- **Axial push or pull** — tailored scattering gives positive radiation force, or negative
  (tractor / bottle configurations).
- **Lateral translation** — phase-gradient or asymmetric scattering converts axial momentum to
  sideways momentum (Part VIII).
- **Torque and rotation** — vortex beams transfer orbital angular momentum (§6.1).
- **Particle sorting** — density, compressibility, size, and acoustic contrast factor give
  different force responses; separates cells, polymers, metals.
- **Multi-particle assembly** — holographic fields impose many discrete traps or arbitrary paths;
  the basis of acoustic bioprinting and micro-assembly work.

### 12.3 Imaging and information outputs

- **Subwavelength imaging** — negative- or near-zero-index designs manipulate evanescent or
  strongly refracted components; resolution beyond ordinary diffraction-limited behaviour in
  particular configurations.
- **Acoustic cloaking** — transformation-acoustics structures route sound around an object and
  reconstruct the field beyond it. Demonstrated in 3D with perforated plastic networks.
- **Encoded acoustic images** — metaholograms that reconstruct a symbol only at the correct
  frequency, incidence, or decoding configuration (butterfly, sunglasses, text demos).
- **OAM channels** — distinct topological charges are orthogonal (§6.1); candidate multiplexed
  channels for underwater or airborne acoustic communication.

All `[DEAD]` for the audible branch on the coherence argument (Part XI, point 2).

### 12.4 Noise-control outputs

The one family **fully live at voice frequencies**, because it is defined in the audible band to
begin with.

- **Low-frequency absorption in thin panels** — labyrinthine channels, Helmholtz cavities, or
  membrane resonators dissipate far below the frequency the panel thickness would suggest.
- **Reflective low-frequency barriers** — high transmission loss without mass-law thickness.
- **Underwater anechoic coatings** — structured inclusions plus viscoelastic material to cut
  sonar-relevant reflection. Note this is the Part 18.2 stealth family: absorption, not projection.

`[PARAM]` *given, unverified:* published designs report near-perfect absorption at **~125 Hz** from
a thin perforated plate plus coiled cavity — a depth where porous foam would need far more. Useful
as the existence proof that sub-wavelength audible-band control is real. At 125 Hz,
`lambda = 2744 mm`, so a panel of a few tens of mm is `~lambda/100`.

---

## Part XIII — The Governing Constraint: No Passive Gain

`[LAW 13.1] — Conservation.` A passive structure adds **zero acoustic power**. Total radiated
power is bounded above by what the voice supplies, minus viscothermal (Law 5.4) and structural
loss.

This is the audible branch's analogue of Law 2.2: a hard closure that determines what the rest of
the branch is even allowed to attempt.

| Mechanism | What it buys | What it costs |
|---|---|---|
| Impedance matching (horn) | More of the throat's high-impedance energy actually radiates | Length and mouth area |
| Directivity | Higher SPL on-axis | SPL everywhere off-axis |
| Resonant storage | Higher local pressure at `f0` | Bandwidth, and ring-down time |
| Path routing | Sound delivered around a corner or barrier | Wall loss, dispersion |

`[LAW 13.2] — The gain ledger.` Any reported "gain" is a redistribution. Read every dB figure by
asking what was traded: bandwidth, directivity, stored energy, or radiation elsewhere.

`[PARAM]` *given, unverified:* metamaterial localization cavities have measured up to **13 dB local
SPL gain at resonance**. Treat this as proof of compact pressure localization, **not** as a
broadband speech-amplifier specification — §15.5 shows why a 13 dB peak in a 20–50 Hz-wide band is
nearly useless for speech.

---

## Part XIV — Architecture of a Voice Projector

### 14.1 Signal path

```
Raised voice
  -> oval mouthpiece / throat
    -> short smooth transition
      -> coiled or folded channels   (optional; packaging only)
        -> large flared mouth
          -> directional spoken-word projection
```

### 14.2 Geometry `[PARAM]` *assumed, first build*

| Element | Value | Rationale |
|---|---|---|
| Mouthpiece (oval) | 35–50 mm x 60–90 mm | Couples to the mouth without sealing it |
| Transition length | 100–150 mm | Smooth, no abrupt steps |
| Profile | Exponential or tractrix | Standard low-reflection impedance taper |
| Output mouth | 180–300 mm | Bounded above by the 320 mm bed |

### 14.3 Functional decomposition

Do **not** tile a uniform cavity array across the whole path. Split by function:

- **Main horn** — broadband transmission and forward radiation. This is the acoustic engine.
- **Three or four tunable side cavities** — selectively alter problem bands. Threaded caps or
  sliding volume plugs so resonance can be swept experimentally instead of reprinted.
- **Coiled segments** — compress physical path length only. Smooth walls, generous cross-section.
  **Governed by Law 5.4:** keep channel height at or above ~5 mm, or the `sqrt(f)/h` tilt strips
  the consonant band (§5.5, §15.6).
- **Removable output plate** — swap straight channels, diffraction slots, shallow phase contours
  without reprinting the body.

`[LAW 14.4] — Mouth-to-wavelength bound.` A horn controls directivity only where the mouth is a
meaningful fraction of `lambda = c/f`.

| f | `lambda_air` | 250 mm mouth is... |
|---|---|---|
| 250 Hz | 1372 mm | 0.18 lambda — no directivity |
| 500 Hz | 686 mm | 0.36 lambda — weak |
| 1 kHz | 343 mm | 0.73 lambda — usable |
| 2 kHz | 172 mm | 1.5 lambda — good |
| 4 kHz | 86 mm | 2.9 lambda — strong |

`[RESULT] 14.5` A 200–300 mm mouth meaningfully controls the 1–4 kHz intelligibility band and does
essentially nothing for the 100–300 Hz fundamental range. **The listener will perceive added
clarity and presence before they perceive added chest weight.** Design the device around clarity,
not bass.

> This is the branch's hard ceiling and it is set by the printer bed, not by acoustics. The
> ultrasonic branch's equivalent ceiling — Law 4.1's 1.875 MHz Nyquist — is set by the nozzle.
> **Both branches are bounded at opposite ends of the same hardware.**

### 14.6 Helmholtz sizing

```
  f_H  =  ( c / 2 pi ) * sqrt( A / ( V * L_eff ) )
```

`A` = neck area, `V` = cavity volume, `L_eff` = neck length including end correction. Build for
adjustability; the end correction is the term you will get wrong on paper.

### 14.7 Horn profiles — the Webster laws

Part 14.2 named "exponential or tractrix" without mathematics. Here is the mathematics.

`[LAW 14.7] — Exponential cutoff.` One-dimensional propagation in a duct of slowly varying area
`S(x)` obeys the Webster horn equation; for an exponential flare it has a hard cutoff:

```
  p'' + (S'/S) p' + k^2 p  =  0    ,      S(x) = S_t e^{m x}

  k_axial  =  sqrt( k^2 - m^2/4 )      ->    propagates only above    f_c  =  m c / (4 pi)
```

Below `f_c` the axial wavenumber goes imaginary: the horn does not work *less*, it becomes a
reactive plug that stores and returns energy instead of radiating it. For the ideal infinite horn
cutoff is a cliff; a finite printed horn softens it into a steep slope with mouth-reflection ripple.

`[RESULT]` Printable straight horn — throat 40 x 70 mm oval (mid-range of §14.2,
`A_t = pi a b = 2199 mm^2`), mouth Ø250 (`A_m = 49 087 mm^2`), axial length 280 mm inside the bed:

```
  m    =  ln( A_m / A_t ) / L  =  ln(22.3) / 0.28  =  11.1 m^-1
  f_c  =  m c / 4 pi           =  303 Hz
```

`[LAW 14.8] — Mouth matching.` A mouth of radius `r_m` radiates cleanly only where `k r_m >= 1`:

```
  f_match  =  c / ( 2 pi r_m )          Ø250  ->  437 Hz
```

Below `f_match` the mouth is an impedance step and reflects; the horn interior develops standing
waves and response ripple. **This is Law 2.3's etalon reappearing in the audible branch** — the
mismatch is radiation impedance at the mouth rather than material impedance at a face, but the
physics (partial reflection at a boundary, interference locked to geometry) is identical.

> **The printable low edge is overdetermined.** Cutoff (303 Hz), mouth matching (437 Hz), and
> directivity (Law 14.4: nothing below ~500 Hz) are three independent calculations, and all three
> pile onto the same edge. Nothing is gained by fighting any one of them — the other two remain.

`[RESULT 14.9] — Tractrix on this bed.` The tractrix (the curve whose tangent segment to the axis
has constant length `a`) makes a horn whose mouth radius equals `a`, with `a = c / (2 pi f_c)`:

```
  bed-maximum mouth  a = 160 mm   ->   f_c = 341 Hz
  a 300 Hz tractrix needs  a = 182 mm  ->  Ø364 mm  — exceeds the 320 mm bed
```

So on this machine the tractrix caps at ~340 Hz — essentially the same floor as the exponential.
*(Geometry note: the tractrix revolved about its axis is Beltrami's pseudosphere, the surface of
constant negative curvature `K = -1/a^2`. The horn mouth is a printed slice of hyperbolic
geometry; its cutoff frequency and its curvature are the same number read in different units.)*

`[DEAD]` **Bass horn.** A 100 Hz cutoff needs `m = 3.66 m^-1`, i.e. `L = 0.85 m` at the same area
ratio — and, independently, `k r_m >= 1` demands a Ø1.09 m mouth. Length is fixable (§14.10);
the mouth is not. Closed on aperture, permanently, by the bed.

### 14.10 Folding — you can coil the path; you cannot coil the mouth

`[LAW 14.10]` The natural fold for a flaring duct is the logarithmic spiral `r = r0 e^{k theta}`
(constant pitch angle, so the channel never pinches). Its arc length between radii is closed-form:

```
  s  =  ( sqrt(1 + k^2) / k ) ( r2 - r1 )
```

`[RESULT]` The dead bass horn's 0.85 m path folds comfortably: from `r = 30 mm` to `r = 140 mm`,
`k = 0.131` (pitch angle 7.4 deg), sweep `= ln(140/30)/k = 11.8 rad` — **1.9 turns**, flat inside
the bed. The length constraint dissolves. The mouth constraint does not move at all: 0.93 m^2 of
radiating aperture must exist as physical boundary on the outside of the print, and no folding of
the interior creates exterior surface.

> Coiling is a transformation of the path integral, not of the boundary condition. A spiral
> defeats `m`; nothing printable defeats `k r_m >= 1`. **This is Law 14.4 re-derived from the
> opposite side, and it is why Part V calls coiling "packaging only."**

### 14.11 The precedent that beats the bound — the cochlea

There is one known passive structure that is simultaneously broadband (~10 octaves), audible-band,
and radically sub-wavelength: a coiled tapered channel, 2.75 turns, 35 mm uncoiled, handling 20 Hz
where the wavelength in its own fluid is ~75 m — the organ is `lambda/2000` at its low edge. It
does **not** use cavity resonance. The basilar membrane grades in stiffness by ~4 orders of
magnitude from base to apex; the fast pressure wave couples into a **slow structural traveling
wave** whose local wavelength collapses as it approaches its characteristic place, and the energy
of each frequency is deposited at a different position (von Bekesy, Nobel 1961). The
frequency-place map is log-periodic:

`[PARAM]` *given (Greenwood 1990), human:*

```
  f(x)  =  165.4 ( 10^{0.06 x} - 0.88 )  Hz        x in mm from apex

  f(0) = 20 Hz ,   f(35) = 20.7 kHz     ->    ~3.5 mm per octave, a logarithmic ruler
```

The lesson for this document: the winning passive architecture for broadband sub-wavelength
audible work is a **graded dispersive line** — a slow wave on a stiffness gradient — not a
resonant cavity. It evades Law 15.13's bandwidth/ring-time trap because each position is broadband
by itself; selectivity comes from *where* energy lands, not *how long* it rings. This is exactly
the trick "acoustic rainbow" metamaterial ridge arrays copy, and it is the same family as §12.4's
125 Hz coiled absorber. Registered as open question 11.

---

## Part XV — Speech-Band Tuning and the Q Budget

### 15.1 Tuning targets

| Perceptual aim | Band | Passive feature |
|---|---|---|
| Reduce rumble / plosives | 80–180 Hz | High-pass mouth geometry; avoid oversized sealed cavities |
| Add vocal body | 200–500 Hz | Large cavity or gentle horn loading |
| Reduce boxiness | 400–700 Hz | Tuned side-branch notch / cavity damping |
| Improve intelligibility | 1–4 kHz | Horn directivity plus mild resonant support |
| Avoid harshness | 2.5–5 kHz | Shallow output diffuser or damped cavity |
| Preserve air / sibilance | 5–8 kHz | Short smooth final outlet; no narrow labyrinths |

Target **broad 2–5 dB shaping across several bands**, never a single 13 dB peak.

### 15.2 The voice is already a resonant system

Vocal folds produce a harmonic-rich source; the tract (larynx, pharynx, mouth, tongue, lips, nasal
path) filters it into **formants**, and the moving formant pattern is what distinguishes vowels and
identifies the speaker.

```
vocal-fold source -> many harmonics -> vocal-tract resonances -> speech
                          F1         F2         F3
                           ^          ^          ^
spectrum:              ___/ \______/ \______/ \_____
```

`[LAW 15.3] — Fixed resonance is phoneme-blind.` An external cavity peaking at, say, 990 Hz boosts
whatever passes near 990 Hz regardless of whether the talker is saying "ah," "ee," "s," or "k." It
competes with, and can dominate, the naturally moving formant pattern.

The device is a multiplier on the voice spectrum:

```
  P_out(f)  =  H_structure(f) * P_voice(f)
```

A projector wants `H_structure(f)` comparatively smooth. A high-`Q` cavity instead gives:

```
Gain
 ^
 |                  /\
 |                 /  \        narrow, high resonance
 |________________/    \_______________________  Frequency
                 990 Hz
```

"Hollow" is the perceptual name for this: a few frequencies strongly reinforced, the gaps between
them relatively weak, so the ear hears pipe, helmet, tunnel, or box rather than a person.

### 15.4 The five failure modes

| Mechanism | Acoustic result | What the listener hears |
|---|---|---|
| Narrowband boost | One formant region dominates | Honky, nasal, hollow |
| Narrowband rejection nearby | Detail between resonances weakened | Missing consonants, poor word recognition |
| Long ring-down | Previous phonemes persist after the mouth has moved on | Smearing, echo, blurred syllables |
| Multiple reflections | Frequency-dependent addition and cancellation | Comb-filtered, phasey |
| Uneven directivity | Tone varies with angle | Clear at one spot, muffled elsewhere |

Ring-down is the central one. A resonator does not stop when the voice changes:

```
  Q  ~=  2 pi f0 * E_stored / E_lost_per_cycle          delta_f  ~=  f0 / Q
```

### 15.5 `[RESULT]` Q budget at f0 = 990 Hz

| Q | Bandwidth | Speech outcome |
|---|---|---|
| 5 | 198 Hz | Broad coloration — usable |
| 10 | 99 Hz | Obvious vowel/formant emphasis |
| 25 | 40 Hz | Strong one-note, tube-like |
| 50 | 20 Hz | Very selective ringing; poor natural speech |
| 100 | 10 Hz | A tone resonator, not a speech projector |

This resolves the 13 dB claim in §13.2: 13 dB centred in a 20–50 Hz band touches a thin spectral
slice of speech. Dramatic at a probe microphone, near-irrelevant to a listener.

> **Cross-branch note.** The ultrasonic branch *wants* the high-Q end — Law 3.3's coiled-cell
> resonance is a design tool there. Here the same physics is the primary failure mode. `Q` is the
> clearest single parameter separating the two branches.

### 15.6 Consonants are the vulnerable class

Vowels are sustained and tolerate tonal coloring. Consonants are short and broadband — stops
(p, t, k, b, d, g), fricatives (s, f, sh, th), affricates (ch, j). A ~1 kHz cavity can make vowels
sound substantial while failing to radiate the higher-frequency transients that separate
*seat/sheet*, *tin/kin*, *pin/bin*. Intelligibility depends disproportionately on **1–4 kHz**:
removing content below 500 Hz leaves speech broadly understandable, while low-passing near 1 kHz
degrades it sharply. Added reverberation consistently reduces intelligibility by obscuring temporal
detail; listeners adapt somewhat to a consistent signature, but the degradation remains.

**This band is attacked from two independent directions** — by resonance (§15.5) and by coiled-
channel `sqrt(f)` tilt (§5.5). Both must be budgeted, and neither shows up in an overall dBA
figure.

### 15.7 `[LAW]` Comb filtering — the subtle failure

Two paths to the listener (coiled path plus direct leakage, or two outlets) sum as

```
  P_total(f)  =  P_0(f) * [ 1 + alpha * exp( -j 2 pi f tau ) ]
```

giving regularly spaced peaks and nulls at spacing `delta_f ~= 1/tau`.

```
Level
  ^
  |    /\    /\    /\    /\
  |___/  \__/  \__/  \__/  \___ Frequency
      peak null peak null
```

`[RESULT] 15.8`

| Path mismatch | `tau` | `delta_f_comb` | Verdict |
|---|---|---|---|
| 10 cm | 0.29 ms | 3.43 kHz | One broad tilt — tolerable |
| 50 cm | 1.46 ms | 686 Hz | Comb sits inside the formant region — destructive |

Worse, the mismatch changes with listener head position, so the timbre swims. **Rule:** if you use
multiple channels, either match path lengths to a small fraction of the shortest wavelength of
interest, or separate the outlets far enough that they never recombine at one listener.

### 15.9 When high Q *is* the goal

A high-`Q` resonator is not a bad device, it is a different device. Choose it deliberately for a
ritual or sculptural voice-coloring chamber, a sustained-note resonator, an acoustic instrument, a
vowel/formant demonstrator, a frequency-selective voice gate, or a local pressure-field experiment
read by a microphone. Do not choose it for a megaphone.

### 15.10 The tract is a quarter-wave instrument

`[LAW 15.10]` A tube closed at the glottis and open at the lips resonates at odd quarter-waves:

```
  F_n  =  ( 2n - 1 ) c / ( 4 L_tract )
```

`[RESULT]` `L = 170 mm` gives **504, 1513, 2522 Hz** — the neutral vowel (schwa), the sound the
tract makes when nobody is steering it. Every vowel is a perturbation of this stack: constrictions
and lip rounding pull individual resonances off their quarter-wave posts. A 145 mm female tract
scales the whole stack up 17% (591, 1774, 2957 Hz).

The vowel plane, from the §0.4 registry (F1 down the page = jaw opening; F2 across = tongue
front/back):

```
              front  <----------- F2 -----------> back

  close    ee (270, 2290)                  oo (300,  870)
             |                               |
             ih (390, 1990)                aw (570,  840)
             |                               |
  open       ae (660, 1720)                ah (730, 1090)

                                     ^^^^^^^^^^^^^^^^^^^
                                     the §15.3 cavity at 990 Hz lands here:
                                     above every F1, below every front-vowel F2,
                                     inside back-vowel F2 territory (840-1090 Hz)
```

**Law 15.3 now has its mechanism.** A fixed 990 Hz peak sits above every male F1 and below every
front-vowel F2 — squarely inside the F2 band of *ah / aw / oo*. The ear parses it as a spurious
back-vowel formant welded onto every phoneme: "ee" acquires a ghost of "aw." That is precisely the
timbre listeners call hollow, honky, speaking-into-a-jug — §15.3's "hollow" is a phantom back
vowel, and the vowel plane says which one. And because female formants run 15–20% higher, the same
cavity collides with a different part of each talker's vowel space: **the coloration is not even
the same insult to every voice.**

> The tract is an instrument retuned every few tens of milliseconds, and that movement *is* the
> message. A printed cavity is one key held down through the whole performance.

### 15.11 The singer's formant — the opera house as existence proof

An unamplified voice fills a 2000-seat hall over a hundred instruments. Law 13.1 binds the singer
exactly as it binds the horn: no gain exists, only placement. Where the trained voice places it:

- The orchestra's long-term average spectrum peaks near 450 Hz and falls steeply above — the
  masker is loudest in the bass and quiet in the treble.
- The ear is most sensitive near 3–4 kHz (the ear canal's own quarter-wave, §0.4).
- Between the falling masker and the peaking receiver there is a spectral window, and the trained
  voice clusters resonances F3–F5 into it — the **singer's formant**, the ring Italians call
  *squillo* — by narrowing the epilaryngeal tube until (area ratio ≲ 1:6, §0.4) it decouples from
  the pharynx and resonates as its own quarter-wave pipe:

```
  f_epi  =  c / ( 4 L_epi )  =  343 / ( 4 x 0.028 )  ~  3.1 kHz      (25-35 mm -> 2.4-3.4 kHz)
```

`[RESULT 15.11]` Three quarter-wave resonators align on one band, each for an unrelated reason:

| resonator | scale | band | why it is there |
|---|---|---|---|
| epilarynx (source) | 25–35 mm tube | 2.4–3.4 kHz | trained: put energy where the orchestra is quiet |
| horn mouth (device) | Ø200–300 mm aperture | 1–4 kHz (Law 14.4) | printed: the only band the bed can steer |
| ear canal (receiver) | ~25 mm tube | ~3.4 kHz | evolved: the ear's own horn throat |

The larynx, the printer bed, and the ear converge on the same octave by three independent
constraints. §14.5's conclusion — *design around clarity, not bass* — is the same discovery every
conservatory makes, and it upgrades from advice to strategy: **the device is a prosthetic
epilarynx.** It does by aperture what the singer does by tube-narrowing — adds nothing, places
everything. Opera is Law 13.2's ledger with a two-century experimental record.

### 15.12 Vibrato through a resonator — FM becomes AM

Operatic vibrato modulates pitch at 5–7 Hz with ±3–6% excursion (§0.4). Feed that through a
Lorentzian resonance and the frequency motion is transcribed into amplitude motion:

`[LAW 15.12]`

```
  |H(delta)|^2  =  1 / ( 1 + (2 Q delta / f0)^2 )         delta = detuning from f0

  AM depth  =  10 log10( 1 + (2 Q e)^2 )   dB    ,    e = fractional pitch excursion
```

A harmonic centered on the peak crosses it twice per vibrato cycle — flutter at 10–14 Hz. A
harmonic parked on the skirt rides the slope — the resonator is an **FM discriminator**, and the
tremolo comes out at the vibrato rate itself, 5–7 Hz.

`[RESULT]` at ±1 semitone (`e = 0.06`):

| Q | AM depth |
|---|---|
| 5 | 1.3 dB |
| 10 | 3.8 dB |
| 25 | 9.9 dB |
| 50 | 15.6 dB |
| 100 | 21.5 dB |

(±0.5 semitone halves `e`: Q=25 → 5.0 dB, Q=50 → 9.8 dB. Tracking lag is no escape: the amplitude
time constant `Q / (pi f0)` is 16 ms at Q=50, f0=990 Hz, against an 83 ms vibrato half-period —
the resonator follows faithfully, so the full depth develops.)

A 5–7 Hz amplitude wobble is, note for note, the fault singing teachers reject as a broken vibrato
— the *wobble*. Prosodic pitch movement in speech does the same at 1–4 Hz. **A high-Q cavity does
not merely color a moving voice; it transcribes pitch motion into loudness motion at exactly the
rates the ear reads as a damaged instrument.** This is a sixth failure mode for §15.4's table, and
it is invisible to any swept-sine measurement, because a swept sine never sits still and shakes.

### 15.13 Ring-down against the phoneme clock

A resonator's envelope decays as `exp( - pi f0 t / Q )`.

`[LAW 15.13a]`

```
  t_60  =  ln(10^3) Q / ( pi f0 )  =  2.2 Q / f0
```

`[LAW 15.13b]` — and since the half-power bandwidth is `Delta_f = f0 / Q`:

```
  Delta_f * t_60  =  2.2
```

Bandwidth and ring time are one number wearing two units — Gabor's uncertainty relation in lab
clothes. Demanding the tail die inside a stop-consonant closure (~50 ms, §0.4) therefore demands
**at least 44 Hz of bandwidth at every center frequency**, and a Q ceiling that rises linearly
with `f0`:

`[RESULT]` `Q_max = f0 * 0.050 / 2.2`:

| f0 | Q_max | t_60 at Q_max |
|---|---|---|
| 250 Hz | 6 | 50 ms |
| 500 Hz | 11 | 50 ms |
| 1 kHz | 23 | 50 ms |
| 2 kHz | 45 | 50 ms |
| 4 kHz | 91 | 50 ms |

Two cross-checks land exactly where they should: §15.5 judged Q=25 at 990 Hz "strong one-note" —
it sits on this ceiling to within rounding. And §15.1's "mild resonant support" being licensed at
1–4 kHz but nowhere lower is now a law rather than a taste: **sharpness is a privilege of the
treble.** A 250 Hz cavity is permitted almost no selectivity at all before it starts speaking over
the talker.

One more clock. Running speech delivers syllables at 4–7 per second — the same band as operatic
vibrato, and the band where the ear's sensitivity to amplitude fluctuation peaks (~4 Hz). The ear
is tuned to the syllable rate; vibrato deliberately lives in that same window. That single
coincidence is why a smeared syllable (§15.4) and a wobbling tone (§15.12) are both so mercilessly
audible: the failure modes land on the ear's most-watched channel.

### 15.14 Spoken vs sung — a line spectrum samples the structure

A voice is not noise. It is a harmonic comb at spacing `f0` sampling `H_structure(f)` (§15.3) at
discrete points.

`[LAW 15.14]` A structural feature of bandwidth `Delta_f` is *guaranteed* audible only if
`Delta_f >= f0`. Narrower features are sampled by luck, re-drawn at every pitch change.

- **Spoken male voice**, `f0 ~ 120 Hz`: a comb every 120 Hz. Every feature wider than ~120 Hz is
  always struck; the transfer function is fully interrogated and every Q budget above applies at
  full force.
- **Soprano at A5** (880 Hz): nine harmonics below 8 kHz. A Q=25 notch at 2 kHz is 80 Hz wide —
  it falls silently between harmonics, or lands on one and deletes it outright, and which of the
  two happens is re-decided at every note. The instrument the singer experiences is a lottery
  drawn once per pitch.

This resolves §15.9 from the other side. The resonator device that ruins speech can be a
legitimate *instrument* for the sung voice — but only with a player in the loop, because sopranos
already play this game against their own anatomy: at pitches where `f0` climbs above a vowel's F1
they retune, opening the jaw to raise F1 up to meet the fundamental. The tract follows the note.
A printed cavity cannot follow — but §14.3's threaded caps and sliding plugs are exactly a jaw for
the device. Score the deliberate high-Q build (§15.9) like an instrument part: tunings within
reach of a hand mid-phrase, one resonance per harmonic of the intended note, and no pretense of
serving speech.

---

## Part XVI — Build and Measurement Sequence

1. **Plain horn first.** Uncoiled, smooth. Measure against an unmodified voice source at 0, 30,
   60, and 90 degrees off-axis.
2. **Swept sine before voice.** Drive the throat with a small driver (or the existing 12 in
   speaker / 275 mm membrane rig), 100 Hz – 8 kHz, and capture the transfer response. Look for
   deep nulls, narrow high peaks, and extended ring-down in a spectrogram or waterfall.
3. **One branch at a time.** Add a single removable Helmholtz or quarter-wave branch tuned near
   500 Hz, 1 kHz, or 2–3 kHz. Compare spectra.
4. **Coil last.** Only after the plain horn is shown not to be excessively lossy. Coiling buys
   compactness; Law 5.4 says narrow channels spend the energy you are projecting, preferentially
   at the top of the intelligibility band.
5. **Map SPL and spectrum** with a calibrated USB measurement mic at fixed distances. Report
   overall dBA *and* narrowband FFT / third-octave — louder and more intelligible do not track
   together.
6. **Test words, not tones.** Consonant-rich phrases ("Peter Piper picked a peck...") expose the
   failures in §15.4 within seconds.

### 16.1 Candidate deliverables

- **Directional speech beam** — on-axis clarity, off-axis rolloff.
- **Localized listening point** — horn plus lens/reflector maximizing SPL at a chosen position.
- **Voice-coloring resonator** — cathedral, metal throat, robotic vowel, produced acoustically
  (the deliberate §15.9 case).
- **Frequency-selective voice gate** — output strong only when specific formant bands excite it.
- **Multi-exit voice sculpture** — one mouthpiece, several coiled paths of different length,
  spatially separated outputs with different delays and timbres (deliberate comb filtering,
  §15.7 inverted into a feature).
- **Passive communication tube/panel** — folded path across a barrier or around a corner.
- **Spoken-word focal reflector** — curved panel or array aimed at a zone; architectural acoustics
  rather than loudspeaker.

### 16.2 Measurables that connect geometry to a digital twin

Microphone pressure maps, FFT-based insertion loss, beam angle, on/off-axis SPL difference,
ring-down time (hence `Q`), and word-recognition score.

---

## Part XVII — Open Questions, Ranked

### 17.0 Shared across both branches

1. **`c_L` of printed PETG.** Drives everything through Law 1.3 (elasticity −2.15). A 1%
   measurement buys 2.1% on thickness. **Method:** pulse-echo time-of-flight through a printed
   ~20 mm calibration block from the production profile — not filament, not a datasheet (§4.6).
   Also yields `Z = rho c`, hence ripple depth, and the margin to the `n = 1` singularity.
   *(Numbered §17.1 throughout the document.)*
2. **Real viscothermal loss in FDM channels** (§5.3, §5.5), where surface roughness (~50–100 um Ra)
   exceeds the boundary layer by 5–10x at 40 kHz and by ~1x at 1 kHz. **The roughness-to-`delta_v`
   ratio inverts between branches**, so one correction factor will not serve both — measure per
   band. *Partially bounded (Rev 4):* Zvoníček et al. 2023 measure FDM PET-G at Sa 9.7–13.9 μm /
   Sz 76–111 μm (§0.2) — an order of magnitude below the Ra assumed above, on flat surfaces. If
   vertical channel walls are comparable, the §5.3 "1.5–3x worse" caveat is conservative in the
   audible band and roughly right at 40 kHz only if ridge peaks (Sz-scale) dominate. Vertical-wall
   measurement still needed.

### 17.A Ultrasonic branch

3. **MHz absorption in printed PETG.** Unquantified, plausibly larger than the 0.83 dB impedance
   loss. Thermoplastic attenuation runs on the order of dB/mm at MHz and scales ~`f^1.5-2`. If
   large, it independently forces 1 MHz over 2 MHz regardless of the Nyquist ceiling. **Second
   thing to measure.**
4. **Can End-to-End 3D optimization recover unwrapped achromaticity at wrapped thickness?** Laws
   3.2 and 4.5 share a root cause — the mod-2-pi discontinuity — so one optimizer targeting cliff
   placement should improve bandwidth and diffraction error together. Not tested in the source
   material.
5. **Poisson ratio of printed PETG.** Decides whether §9.3's lattice is merely marginal or
   categorically dead. Low priority while §9.4 stands.
6. **Multi-bottle static holograms vs. phased-array vortex traps at matched aperture.** §6.4 says
   static plates are viable for large particles via bottles; no head-to-head on trap stiffness or
   capacity exists in the source material.
7. **Sawtooth angular steering resolution** (§8.2). Force magnitude is bounded; resolution is not.

### 17.B Audible branch

8. **On-axis SPL delta of a 250 mm-mouth printed horn** over an unaided raised voice, at 2 m, in
   the 1–4 kHz band. Law 13.1 bounds it from above by conservation and Law 14.4 says the mechanism
   is directivity rather than gain, but **nothing in this document predicts the number.** It must
   be measured (Part XVI step 1). This is the branch's single most important unknown — every
   deliverable in §16.1 is contingent on it being non-trivial.
9. **Structure-borne re-radiation.** Whether a printed PETG wall is stiff enough to avoid
   re-radiating at speech levels, or whether it needs mass loading / damping. Note the tension
   with §9.x: mass loading raises `rho` and therefore `Z`, the opposite of what the ultrasonic
   branch wants from the same material.
10. **Do the two 1–4 kHz attack vectors compound or partially cancel?** Resonance loss (§15.5) and
    coiled-channel `sqrt(f)` tilt (§5.5) both target the intelligibility band. §5.5 is a smooth
    tilt and §15.5 is a peaked response, so a cavity tuned into the tilt might partially
    compensate. Untested, and worth a swept-sine experiment before committing to a coiled build.
11. **A printable graded dispersive line** (§14.11). The cochlea proves broadband sub-wavelength
    passive audible-band control is physical — via a slow wave on a stiffness gradient, not a
    cavity. Can an FDM PETG rib array (graded rib height/width along a channel) produce a
    measurable frequency-place separation across the 320 mm bed in the speech band? If yes, it is
    a passive spectral analyzer and the natural mechanism for §16.1's frequency-selective voice
    gate — and it sidesteps Law 15.13's bandwidth/ring-time trap entirely.
12. **Vibrato-FM audibility of side branches** (§15.12). Law 15.12 predicts flutter depth from Q
    in closed form, and no stationary swept sine can observe it. Verify with a sung glide or a
    synthesized 6 Hz vibrato tone through one printed Helmholtz branch, mic at the mouth: measured
    AM depth vs `10 log10(1 + (2Qe)^2)` is a two-point experiment.

---

## Part XVIII — Source Lineage and Unverified Claims

### 18.1 Lineage

- **2015** — Bristol/Sussex, *Nature Communications*: first single-sided acoustic tractor beam,
  64-element loudspeaker array, generating tweezer, vortex, and cage force fields.
- **2016** — Max Planck Institute for Intelligent Systems: array replaced by a single large
  transducer plus a 3D-printed contoured block; ~100x finer resolution as the aperture element
  dropped from transducer diameter to printer feature size.
- **Subsequently** — Asier Marzo (Bristol, later UPNA/Navarre): open-source, hobbyist-buildable
  "sonic tractor beam"; published fabrication files and an Instructables guide for the
  metamaterial delay-line structure.
- **2018** — *Nature Communications*: horn-like space-coiling structures modulating phase and
  amplitude simultaneously.
- **LBNL** — bottle beams demonstrated with no metamaterial, via phase engineering across a
  conventional array (§6.4).
- **ASA 188th meeting** — underwater asymmetric sawtooth patches (§8.1).
- **2025 arXiv** — unified iterative-backpropagation framing plus an open-source Python full-stack
  pipeline for propagation modeling, phase retrieval, and hardware control (§7.2).
- **2023** — Zvoníček, Vašina, Pata & Smolka, *Polymers* 15:2025 (MDPI, open access; local copy in
  repo root): ISO 10534-2 impedance-tube sound reflection of FDM PET-G/ASA/PLA, Grid/Gyroid/Cubic
  infill, 0.1–0.5 mm layers, 200–3200 Hz, plus Zygo areal roughness. The only entry in this
  lineage measured on **this shop's material and process in the audible band**; parameters
  registered in §0.2, consumed by §17.0-2 and by the voice-horn POC study
  ([voice-horn-pocs.md](voice-horn-pocs.md), Rev 5, findings F1–F5).

The audible branch additionally stands on a speech-science spine, cited once here and consumed by
§0.4, §14.11, §15.10–15.14:

- **1952** — Peterson & Barney: vowel formant means across talkers (§0.4, §15.10).
- **1961** — von Bekesy, Nobel in Physiology or Medicine: passive traveling-wave cochlear
  mechanics (§14.11).
- **1974** — Sundberg: the singer's formant explained as epilaryngeal quarter-wave clustering
  (§15.11).
- **1990** — Greenwood: the human frequency-place map (§14.11).

**Organizing thesis:** every result above follows from *where the phase map lives* — in electronics
(reconfigurable in `t`, coarse in `x`) or in geometry (frozen in `t`, fine in `x`). The audible
branch obeys the same thesis one level up: the singer's phase map lives in muscle (reconfigurable
in `t`), the horn's in PETG (frozen in `t`), and Part XV is the accounting of what freezing costs.

### 18.2 Adjacent, not core — naval impedance work

Hard inclusions in viscoelastic matrices as anechoic hull coatings, impedance-matched to water while
damping structural radiation. **Design goal is absorption (stealth) — the opposite sign of the
objective here.** Submarine propeller-shaft radiation studies confirm impedance mismatch as the
general control principle for structure-borne sound coupling into the water column. Relevance to
Part II is methodological — same governing quantity, inverted target. **Do not cite as supporting
evidence for manipulation plates.**

### 18.3 Claims not independently verified

Carried from the source material without confirmation. **Quarantined here deliberately — do not
promote into a law section without a citation.**

*Ultrasonic branch:*

- The 188th ASA meeting attribution for the sawtooth patch work.
- The 2025 arXiv paper and its open-source pipeline.
- The "roughly 100x" resolution figure — at least self-consistent, since §4.2 recovers it exactly
  from `10 mm pitch / 0.1 mm feature`.
- Photopolymer resin at 2.5 MRayl (used only as the resin baseline; PETG numbers do not depend
  on it).

*Audible branch:*

- 13 dB local SPL gain in metamaterial localization cavities (§13.2). Bounded in interpretation by
  §15.5, not verified as a measurement.
- Near-perfect absorption at ~125 Hz from a thin perforated plate plus coiled cavity (§12.4).
- 3D acoustic cloaking demonstrated with perforated plastic networks (§12.3).
- Encoded metahologram demos — butterfly, sunglasses, text (§12.3).
- Speech-science values in §0.4 (formants, vibrato, singer's formant band, epilarynx geometry,
  Greenwood map) are population means from the textbook literature (§18.1 spine). Sound as
  physics; but no individual talker matches them. Treat as design centers, never as measurements
  of the actual user's voice — Part XVI step 2 measures the real one.

---

## Revision Log

| Rev | Date | Change |
|---|---|---|
| 1 | 2026-08-21 | Initial. Laws 1.1–1.3, 2.1–2.4, 3.1–3.3, 4.1–4.5, 5.1–5.3, 6.1–6.4, 7.2–7.3, 8.2, 9.2–9.3. Parameters bound to Bambu H2S / 0.4 mm nozzle / PETG. Material routes 9.1–9.3 closed. Design envelope Part X set at 1 MHz water. |
| 2 | 2026-08-21 | **Audible branch merged in** from the former `passive-voice-projection.md`, which is now deleted — this document is the single source. New Parts XI–XVI (scope pivot, output taxonomy, no-passive-gain, horn architecture, Q/intelligibility, build sequence); Laws 13.1–13.2, 14.4, 15.3, 15.7 renumbered from that document. Back matter renumbered XI→XVII, XII→XVIII, and both branches' open questions and unverified claims merged. **New derived content:** Law 5.4 (closed-form `alpha ~ sqrt(f)/h`, verified against the 40 kHz case) and Result 5.5 (audible coiled-path loss table), which quantify the former doc's open question on audible-band channel loss and give Part XIV its ~5 mm channel floor. Audible reference wavelengths added to Part 0. |
| 3 | 2026-08-21 | **Horn geometry and voice-source physics.** New §0.4 (voice registry: f0 ranges, Peterson–Barney formants, vibrato, epilarynx, ear canal, phoneme clocks). Part XIV: Law 14.7 (Webster exponential cutoff — printable straight horn lands at 303 Hz), Law 14.8 (mouth matching `k r_m >= 1` at 437 Hz; identified as the audible-branch etalon, same physics as Law 2.3), Result 14.9 (tractrix caps at 341 Hz on this bed; pseudosphere note), Law 14.10 (log-spiral fold arc length — *coil the path, not the mouth*; bass horn closed on aperture), §14.11 (cochlea precedent, Greenwood map, graded-dispersive-line lesson). Part XV: Law 15.10 + vowel plane (Law 15.3's "hollow" mechanized: a 990 Hz cavity is a phantom back-vowel F2), Result 15.11 (singer's formant; epilarynx / horn mouth / ear canal align on 2–4 kHz — the device is a prosthetic epilarynx), Law 15.12 (vibrato FM→AM transcription with closed-form depth; new failure mode invisible to swept sine), Laws 15.13a/b (`Delta_f * t_60 = 2.2`; 44 Hz minimum bandwidth; Q ceiling linear in f0 — cross-checked against §15.5's Q=25 verdict), Law 15.14 (line-spectrum sampling; sung-voice lottery; §15.9 rescored as a playable instrument). Open questions 11–12 added; speech-science lineage added to §18.1; textbook-values caveat added to §18.3. |
| 4 | 2026-08-23 | **First measured audible-band wall data.** Zvoníček et al. (*Polymers* 2023) registered: §0.2 gains a measured `[PARAM]` block (PET-G `beta_m` by infill structure and layer height; Sa/Sq/Sz roughness), §17.0-2's roughness assumption partially bounded (measured Sa an order of magnitude below the assumed Ra; vertical walls still open), §18.1 lineage entry added. Design consequences applied in voice-horn-pocs.md Rev 5 (F1–F5: shell+Grid sandwich walls, 0.3 mm body layers, Cubic ban, RULE W3 measured ceiling). |
| 5 | 2026-08-23 | Glottal source resistance registered in §0.4 (30–100 cgs ohm; the voice as a finite-impedance flow source). Consumed by the voice-horn Forge cycle 2 (voice-horn-pocs.md Rev 6), where it closed the ideal-source tiny-throat artifact and showed maximum power transfer wants `S_t ≈ rho c / Z_s` — smaller than anatomy permits, so the throat pins to the anatomical floor. |
