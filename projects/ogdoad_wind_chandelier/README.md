# The Hymn of the Ogdoad — a wind-blown Helmholtz chandelier

Nine pigeon-whistle (*geling*) voices, tuned to the nine Solfeggio tones and hung as a geocentric
cosmos. **Earth** (174 Hz) is the big gourd at the centre. A wheel carries the seven planetary
spheres in the Chaldean order the Hermetic *Poimandres* climbs (Moon, Mercury, Venus, Sun, Mars,
Jupiter, Saturn), and at the top the eighth sphere, the **Ogdoad** (963 Hz). A dove-tail vane keeps
every intake facing the wind. A crown takes the swivel, and a sand-filled plumb-bob hangs below as
ballast.

Every voice shares one whistle cut-up (window W = 2.4 mm), so each speaks in its own band of wind
speed, in order of pitch. A breath of air wakes Earth, and only a fresh breeze reaches the Ogdoad.

Exposé page: [web/hymn-of-the-ogdoad.html](../../web/hymn-of-the-ogdoad.html), which is staged in `site/`.

| | |
|---|---|
| ![Three-quarter](preview_three_quarter.png) | ![Front](preview_front.png) |
| ![Side](preview_side.png) | ![Head section](preview_head_section.png) |

**Status:** designed and checked in CAD on 2026-10-20. Not yet printed or heard.

## How a voice works

Each gourd is a Helmholtz cavity, and the whistle head glued to its neck holds the mouth. Wind enters a
bellmouth and leaves a 0.7 mm windway as a flat jet. The jet crosses the 2.4 mm window (the
cavity's mouth) and splits on a 30° labium. The jet's waving pumps air in and out of the window
at the cavity's resonance. The head's tall block shelters the window, so the jet crosses slower air
and keeps the shear that makes it wave.

- **Pitch:** the voice is solved as two coupled resonators, `outside —M_mouth— V_upper —M_waist— V_lower`.
  A bottle gourd's waist is itself a short neck. The lower root is the voice, and the upper sits 7–8×
  higher. The mouth's end correction is the radiation mass of the actual W × b slit (closed-form
  rectangle self-potential), not a round-hole formula.
- **Wind bands (model):** jet speed U = 0.85 × wind. A voice speaks while the jet Strouhal number
  St = f·W/U is in 0.10–0.40, best at 0.22 (jet-delay criterion; recorder and ocarina data put the
  best response at 0.15–0.25). With W fixed, a voice's band scales with its pitch.

## The nine voices

| # | Voice | Tone | Gourd Ø × H (mm) | Cavity (cm³) | Window b (mm) | Wind on / best / overblow (m/s) | Angle on wheel | PETG |
|---|---|---|---|---|---|---|---|---|
| 0 | Earth | 174 Hz | 115 × 200 | 955.7 | 29.9 | 1.2 / 2.2 / 4.9 | centre, in front of the hub | 225 g |
| 1 | Moon | 285 Hz | 75 × 133 | 252.6 | 19.4 | 2.0 / 3.7 / 8.0 | 247.5° (left) | 81 g |
| 2 | Mercury · UT | 396 Hz | 57 × 102 | 105.0 | 14.7 | 2.8 / 5.1 / 11.2 | 292.5° (right) | 45 g |
| 3 | Venus · RE | 417 Hz | 54 × 98 | 91.6 | 14.1 | 2.9 / 5.4 / 11.8 | 337.5° (right) | 41 g |
| 4 | Sun · MI | 528 Hz | 45 × 82 | 49.2 | 11.6 | 3.7 / 6.8 / 14.9 | 202.5° (left) | 29 g |
| 5 | Mars · FA | 639 Hz | 39 × 72 | 30.3 | 10.2 | 4.5 / 8.2 / 18.0 | 22.5° (right) | 24 g |
| 6 | Jupiter · SOL | 741 Hz | 36 × 65 | 20.9 | 9.2 | 5.2 / 9.5 / 20.9 | 157.5° (left) | 21 g |
| 7 | Saturn · LA | 852 Hz | 33 × 59 | 14.8 | 8.5 | 6.0 / 10.9 / 24.1 | 67.5° (right) | 19 g |
| 8 | Ogdoad · SI | 963 Hz | 30 × 55 | 11.0 | 7.9 | 6.8 / 12.4 / 27.2 | 112.5° (left) | 17 g |

Angles are measured counter-clockwise from +x, seen from upwind. The rays climb in four tiers. Each
tier's pair is split left/right to balance the wheel, so the climb zig-zags across the axis
(the caduceus).

**One known clash:** 396 and 417 Hz beat at 21 Hz, the roughest possible spacing for pure tones.
They are also neighbours on the wind ladder. Every other pair is mild, because pure sines only clash
when they are close in pitch.

## Checks (`python chandelier.py`)

- Every built cavity solid matches its design volume, so every pitch is within 0.1 cents of its target.
- Mating parts do not interpenetrate. Non-mating parts stay at least 3 mm apart, and the closest pair
  (Earth and the crown) is 22 mm apart.
- All 24 STLs are watertight and single-bodied, fit the 320 mm bed, and pass `lib/validate_mesh.py`
  (PRINT-READY, 0 self-intersections).
- Printed mass is 831 g. The centre of mass is at (−0.7, +4.5, −6.0) mm, so hang from the crown hole
  at y = +6 mm (toward the vane).
- Assembled size: about 388 mm wide × 557 mm deep × 543 mm tall, including the bob on a 190 mm cord.
- The web page's `window.__selfTest()` re-derives every pitch, end correction and band in JS
  (56/56 pass).

## Print

| Part | Orientation | Notes |
|---|---|---|
| 9 gourds | upright | Supports not needed: every surface, inside and out, is within 45° of vertical. 3 perimeters. **Must be airtight.** |
| 9 heads | on their side, D-flat up | 0.4 mm nozzle, 0.12 mm layers, slow outer walls: keep the 0.7 mm windway open and the labium sharp. |
| wheel | rear face down | Hub, 8 spokes, ring, offset foot pegs, crown and pendant lugs. |
| crown, fin, stopper | flat or upright | As exported. |
| vane boom | lying | 274 mm long. |
| plumb-bob | flat top down, point up | Fill with dry sand (~0.4 kg) or steel shot (~1.1 kg). |

Hardware: one M3×20 bolt and nut, a ball-bearing fishing swivel and cord, a steel split ring for the
bob, and epoxy or gap-filling CA.

## Assemble, tune, hang

1. Test and tune each voice **before** gluing: blow across the intake (hair-dryer on cool, or a fan) and
   read a tuner app.
   - Flat: drop plasticine or glass beads into the gourd. 1 % less volume raises the pitch about 9 cents.
   - Sharp: tape across one end of the window. 1 % less width lowers the pitch about 7 cents.
   - Pitch follows air temperature: about −45 cents at 5 °C and +43 cents at 35 °C, relative to 20 °C.
2. Glue the heads on. The D-flat only lets the intake face upwind, and the glue seals the joint.
3. Glue the planets onto the wheel's foot pegs. The pegs are offset, so each gourd goes on only one way.
4. Push Earth's rear boss onto the hub's front peg, the boom into the hub, and the fin onto the boom.
5. Bolt the crown over the top lug and hang the swivel from the +6 mm hole. Move one hole either way
   until the wheel hangs plumb.
6. Hang the filled bob about 190 mm below the bottom lug.
7. Hang the chandelier in open air, ideally under an eave: the upward-facing mouths collect rain.

**Field test:** note each voice's start and stop wind speed with a pocket anemometer. Those numbers
replace the model's ST_ON / ST_OPT / ST_OFF.

## Files

- `voices.py` — acoustic design: gourd profiles, exact volumes, end corrections, 2-DOF pitch solve,
  wind bands, roughness. Writes `voices.json`.
- `chandelier.py` — build123d CAD for every part, cavity re-measurement, layout and balance,
  assembly clearance checks, STL export (`stl/`), previews. `--quick` skips the slow checks and
  `--no-preview` skips the renders.
- `sync_page.py` — injects `voices.json` into the web page and stages the page, previews and STL zip
  in `site/`.

Parametric: change `W_CUT`, `BETA`, `RF_MIN` or the head constants, then run `chandelier.py` followed
by `sync_page.py`.
