# Passive Incense Tornado

A fully passive helical smoke vortex — **no fans, no electronics**. The burning
cone's own heat drives the updraft (stack effect in a clear 100 mm tube), and
twelve tangential slots in the base drum force the incoming makeup air into
rotation (the classic slotted-chimney fire-tornado principle, miniaturized). As
the swirling air converges from the chamber wall (r ≈ 45 mm) to the smoke core
(r ≈ 6 mm), conservation of angular momentum amplifies the spin ~7×, wrapping
the plume into a standing helix. One slot carries an external **mouthpiece
nozzle** — a gentle puff injects extra angular momentum and visibly tightens
the column.

Swirl direction: **clockwise viewed from above.**

## Parts

| File | Part | Print notes |
|---|---|---|
| `vortex_base.stl` | Slotted swirl drum, 176 × 176 × 57 mm | Upright as-is. PETG/ABS, 0.2 mm layers, 3 perimeters, ~15 % infill. No supports (slot roofs are 10 mm bridges; optionally add a dab of support under the mouthpiece tip). |
| `vortex_top_collar.stl` | Exit nozzle cap, Ø106 × 28 mm | Flange down, as oriented. Overhangs ≤ 45°, no supports. |

**Tube:** the build is sized for **100 mm OD × 2 mm wall clear acrylic tube,
15.5" (393.7 mm) long**. There is no printable-tube STL at this size — Ø100 ×
394 mm exceeds the 320 mm bed (the `part = 3` module remains in the .scad for
reference only). A truly transparent wall is what makes the vortex worth
watching; acrylic wins anyway.

## Bill of materials

- 1× printed base, 1× printed collar
- 100 mm OD × 2 mm wall clear acrylic tube, 393.7 mm (15.5")
- 1× metal bottle cap or 2–3 layers of aluminum foil, ≤30 mm — sits in the
  pedestal recess as a heat shield under the cone (**required** — PETG softens
  at ~80 °C)
- Standard upward-burning incense cones

## Assembly & use

1. Drop the metal cap / foil disc into the pedestal recess.
2. Light the cone **with the tube removed**, place it on the disc, wait for a
   steady smoke ribbon.
3. Lower the tube into the base groove; set the collar on top.
4. The column spins up on its own within ~20 s as the tube warms. For a tighter
   helix, give short, *gentle* puffs into the side mouthpiece — it feeds a
   tangential slot, so breath converts directly to swirl.
5. Keep the device out of cross-drafts (HVAC, open windows) — the vortex core
   is only a few millibars of structure.

## Tuning (edit `incense_vortex.scad`, re-export)

- `vane_angle` (68°): higher = more swirl per unit draft, more intake
  resistance. 60–72° is the useful band; too high starves the draft, too low
  gives a lazy column.
- `tube_len`: taller tube = stronger stack draft = faster, more stable vortex.
- `exit_r` (35 mm): smaller exit tightens the core at the outlet but adds
  resistance. 31–42 mm reasonable.
- `slot_w`, `slot_h`, `n_slots`: total slot area is ~0.76× the bore area;
  keep the ratio ≥ 0.6 or the chimney will choke.

## Maintenance & safety

- Incense tar coats the tube interior and slot edges; it changes the
  aerodynamics surprisingly fast. Everything separates without tools — wash
  PETG parts and acrylic with isopropyl alcohol.
- Never leave a burning cone unattended; the printed parts are fuel-adjacent
  plastic, not fireproof. The heat-shield disc is not optional.
- Re-light or swap cones only with the tube lifted off.

---

# Wind Edition (`incense_vortex_wind.scad`)

Outdoor variant optimised to **catch the wind and use its power** — still zero
moving parts. Wind is harvested at both ends of the tube, from any direction:

**Intake pinwheel (bottom).** Twelve tall guide fins radiate from the drum,
each extending one tangential slot's guide wall outward — a fixed turbine
stator. From any wind azimuth, the windward fins form converging funnels
(~4:1 area contraction) that accelerate the breeze and inject it as
same-handed tangential jets. Ram pressure of a 2–3 m/s breeze (~5 Pa) is
roughly **8× the stack-effect draft** of the 394 mm tube, so the vortex
carries far more angular momentum. Leeward slots sit in separated low-pressure
flow, so net circulation stays clockwise.

**Venturi cowl (top).** A stacked-disc chimney-cowl: crosswind is squeezed
through the converging gap between a 45° skirt and a flat hat, dropping static
pressure right over the exit throat (Bernoulli) and sucking the column upward
from any wind direction. The hat also blocks rain and downdrafts. The
gap-perimeter area equals the throat bore area, so the cowl never chokes the
chimney in calm air.

| File | Part | Print notes |
|---|---|---|
| `wind_base.stl` | Pinwheel drum + staked skirt, 242 × 242 × 63 mm | Upright, no supports. Three countersunk holes in the skirt take M4 screws or tent stakes — anchor it; a 450 mm tower in wind needs it. |
| `wind_venturi_head.stl` | Venturi skirt + posts, Ø124 × 44 mm | As oriented (flange down). All overhangs ≤ 45°. |
| `wind_venturi_hat.stl` | Hat disc, Ø140 × 10 mm | Flat side down. Press-fits onto the head's three posts (dab of CA glue optional). |

Tube: same 100 mm interface — the same acrylic tube as the indoor version.

**Behavior:** happiest in a steady 1–4 m/s breeze (porch, open window,
sheltered garden). In gusty wind the core will intermittently burst into
turbulence and re-form — that's vortex breakdown, not a defect. In dead calm
it behaves like the indoor version, running on cone heat alone. Blowing at the
fins from any side replaces the indoor version's mouthpiece.

---

# Vortex Dust Plate (`vortex_dust_plate.scad`)

A flat, open-topped **dust-devil arena** (Ø200 × 79 mm assembled) that lifts a
tracer powder into a visible free-standing vortex using only ambient wind. No
tube, no moving parts.

Designed around **Borozin powder (zinc stearate)** — a micronized air-flow
tracer made to float and follow air currents like smoke, which makes it ideal
here: it entrains at very low velocities, so the vortex core picks it up and
renders itself visible. Any similar micronized powder works.

**How it works.** Twelve guide fins form a converging crown that turns wind
from any azimuth into same-handed tangential jets (58° from radial, ~9 mm
throats, ~4:1 contraction). The flow spirals inward and spins up; the core's
pressure minimum sits over a Ø40 × 3 mm saucer dimple at dead center where a
teaspoon of powder pools. The arena floor is dished (rim high, center low) so
centrifuged powder slides back for re-entrainment. A removable roof ring keys
onto pins on the fin tops, keeping wind in the channels instead of spilling
over, with a Ø136 mm sky opening for the devil to rise through.

| File | Part | Print notes |
|---|---|---|
| `dust_plate_arena.stl` | Arena + fin crown, 200 × 200 × 76 mm | Upright, no supports. **Print in a dark color**: sun-warmed floor adds the thermal updraft real dust devils feed on. |
| `dust_plate_roof.stl` | Roof ring, Ø200 × 5 mm | Flat, no supports. Drops onto the 12 fin-top pins. |

**XL280 variant** (`dust_plate_arena_XL280.stl` + `dust_plate_roof_XL280.stl`):
Ø280 × 99 mm assembled — the practical ceiling of the ≤100 mm flat form.
Roughly 2× the wind-capture area, ~4× the gust resistance (stored angular
momentum grows as D⁴), a taller crown sampling faster wind, and 13 mm throats
that lose proportionally less to wall friction. Same powder saucer principle
(Ø54 dish). Needs a 300-class bed; if yours is 220-class, ask for the
bolted-halves split. Rendered from the same parametric source via `-D`
overrides (see `vortex_dust_plate.scad` parameters).

**Honest expectations.** With no tube to confine it, this produces *dancing,
intermittent* dust whirls in a steady 2–5 m/s breeze — bursts of a few seconds
that wander around the dish, die, and re-form — not a permanent standing
column. Fluff the powder loosely (don't pack it), use a level surface, stake
it down, and favor sunny days with a steady breeze. Powder is consumable:
some escapes through the channels each session; the channel floors (6 mm
below the arena floor) act as gutters that catch much of it for reuse.

**Safety:** zinc stearate is low-toxicity but don't breathe the dust cloud —
and like all metal-soap powders, a dispersed stearate dust cloud is a
**combustible dust**. Never use the dust plate near the burning incense
devices or any open flame.

## Validation

All STLs (all editions) verified watertight and 2-manifold with pymeshlab
(`validate_meshes.py`); each base's second shell is the intentional sealed
air-gap cavity inside the pedestal (a conduction heat break under the cone).
