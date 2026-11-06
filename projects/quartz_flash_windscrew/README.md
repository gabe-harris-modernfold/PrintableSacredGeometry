# Quartz Flash Windscrew

A desk-mounted Archimedes windscrew for a windy ridge. It spins an icosahedron of loose quartz
stones until they flash, and it turns a communication disc that comes to rest under a pointer
when the wind drops: YES · ABOVE · MAYBE · NO · BELOW · NOT HUMAN, with a ring of numbers
0–9 around the rim.

It is styled as if it came out of 1926: a scientific instrument from the year after the Paris
*Exposition des Arts Décoratifs*, with the spirit-board craze of the 1920s on its dial.

![assembly](preview_assembly.png)

| file | what it does |
|---|---|
| `windscrew.py` | parametric generator → `stl/*.stl` (print orientation), `assembly.npz`, part table, clearance check (fails the build if anything touches) |
| `deco_motifs.py` | Art Deco motif library: sunbursts, rays, chevrons, lozenges, clipped-corner plans, Egyptian eye, sun and moon, curved and straight type |
| `physics.py` | the numbers below → `physics.png` |
| `preview.py` | `preview_*.png` renders |

Run from this directory: `python windscrew.py && python physics.py && python preview.py`.

## How quartz makes light when it is struck

Mechanoluminescence covers any light made by mechanical action. *Triboluminescence* is light
from rubbing or crushing; *fractoluminescence* is light from fracture itself. Quartz does both,
and the light is not the stone getting hot as a whole. The bulk stays cool. In the windscrew
two different processes are at work:

1. **Fracture (non-thermal).** A crack opening across quartz leaves the two new faces with
   opposite surface charges. The field across the gap is strong enough to break down the air,
   so the crack lights up as a microscopic gas discharge. Nitrogen's ultraviolet-to-violet
   lines are the standard signature of fracture light in air. Defect centres in the quartz
   lattice add their own luminescence, which is excited by electrons and by the discharge.
   Piezoelectricity is not needed: fused silica glass does this too.
2. **Sliding hot spots (thermal, at micron scale).** When quartz slides on quartz, the real
   contact is a few asperities, each a few µm across. Friction heat lands on them faster
   than quartz can conduct it away. On these spots the Tian–Kennedy flash-temperature model
   reaches the silica softening ceiling (~1950 K) once sliding passes about **1–2 m/s**.
   At **hand-rubbing speed (~0.2 m/s)** it only reaches 500–970 K. On the way up the spot
   crosses the α–β inversion at 573 °C, and the volume jump there cracks the surface, which
   feeds channel 1.

The literature has not settled how much of the orange glow from rubbed quartz is thermal.
Chapman and Walton (*J. Appl. Phys.* 54, 5961, 1983) and Chapman's 1981 thesis compared
continuum spectra with blackbody curves. Other work attributes a 300–600 nm band to
triboelectric gas discharge. `physics.py` makes the thermal question testable instead of
assuming an answer. **Clear quartz is transparent**, so a micron-thick hot skin of it is a very
poor emitter. To clear the Hecht–Shlaer–Pirenne threshold (54–148 photons at the cornea) from
30 cm, the hot skin needs a visible emissivity of about **6×10⁻⁴**. That works out to an
absorption of about **2 /cm**. Clear rock crystal does not absorb that much. Milky quartz,
iron-stained quartz or skins carrying wear debris plausibly do.

**An experiment you can do with your own stones:** load clear crystal, then milky quartz, then
iron-stained quartz, and compare them on a dark night. If the clear stones show only blue-violet
crackle and the stained ones glow orange, the thermal channel is real and depends on absorption.

## The tumbler problem, and the paddle

The rotor runs at 180–1300 rpm. At those speeds loose stones in a container on the shaft are
**pinned to the walls**. The Froude number ω²r/g passes 1 at 128 rpm, which is a 2.5 m/s
wind, and the rotor's cut-in is 3.7 m/s. So a plain tumbler would be silent whenever it turns
at all: the stones would just ride round in the vertex pockets at 2–100 g.

The fix is a **stationary paddle** inside the icosahedron. It hangs from a fixed sleeve
that comes out of the downwind bearing block, around the shaft. Every pinned stone hits it once
per revolution. A stone knocked out of its pocket flies across and hits the stones pinned in the
other pockets, which gives the quartz-on-quartz collisions that flash. One of your stones is
epoxied *through* the paddle near its tip, so it presents quartz on both faces and works in
either spin direction. Quartz on PETG gives no flash, because the plastic softens long before
the contact gets hot. The icosahedron's faces and vertices give the stones somewhere to collect
and corners to rattle off.

The sleeve tip carries its own 608 bearing. That supports the shaft right at the centre of the
container, so the icosahedron and disc are not hanging off the riser bearing alone.

## Numbers (`physics.py`)

| wind m/s | rpm | strikes/s | stone speed m/s | energy/strike mJ | Hertz peak force N | rim tension / strength |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 179 | 18 | 1.0 | 10 | 600 | 9× |
| 8 | 399 | 40 | 2.3 | 48 | 1560 | 13× |
| 12 | 610 | 61 | 3.5 | 111 | 2590 | 15× |
| 20 | 1025 | 103 | 5.9 | 314 | 4830 | 19× |

- **The tip-speed ratio does not depend on wind speed.** The paddle's demand goes as ω³ and
  the wind's power as v³, so v cancels. The rotor settles at λ ≈ 0.8 in every wind, against a
  runaway λ of 2.0. Brightness and strike rate rise smoothly with wind, with no stall and no
  threshold.
- **Every strike exceeds the Hertz cone-crack threshold by about 10–20×.** So each strike chips
  quartz: that is where the fracture light comes from, and the stones *wear*. Expect grit, and
  top the stones up over time.
- **Stone sliding speeds are 1–6 m/s**, 5–30× faster than hand rubbing. Every scrape between two
  quartz faces sits on the flash-temperature ceiling.
- **Cut-in is 3.7 m/s.** That comes from bearing drag plus the stones resting on the paddle.

![physics](physics.png)

**Assumptions and uncertainties:**
- **Rotor.** Ct(λ) = 0.2 (1 − λ/2) gives Cp_max = 0.10. That is conservative for a 300 mm
  Liam-type rotor. Independent tests of the full-size rotor report 0.1–0.2, not the
  manufacturer's figure.
- **Strike energy.** Taken as η·½mv² with η = 1.5, which is the strike itself plus the
  friction loss of re-accelerating the stone.
- **Quartz properties** are room-temperature directional averages. Quartz's conductivity falls
  as it heats, so the flash temperatures are lower bounds, and they are capped at the ceiling
  anyway.
- **Visibility.** The photon count is for one scrape seen by a dark-adapted rod system. Colour
  vision needs more light than that.

## The 1926 styling

![dial](preview_disc.png)

- **Dial.** Dark brown dial with brass relief, the look of Bakelite and brass. Sun beside YES
  and crescent moon beside NO, as on 1920s talking boards. Twin chevrons for ABOVE and BELOW, a
  sunrise for MAYBE, and an Egyptian-revival eye for NOT HUMAN; the Tutankhamun find in 1922
  made that eye everywhere by 1926. Around them sit a railroad-track number ring, sector
  dividers set with lozenges, a 36-ray sunburst and a stepped brass hub under a spire finial.
- **Type.** Century Gothic Bold, the geometric Futura/Twentieth Century family, with wide Deco
  letter-spacing. I tested Broadway, the archetypal Deco face, and rejected it: even at 14 mm
  caps its hairlines lose about 12% of the glyph area to a 0.8 mm minimum stroke, so they
  would not print. Every ornament stroke is at least 0.9 mm.
- **Frame.**
  - The base is a three-tier ziggurat with clipped corners, a sunburst and a MCMXXVI (1926)
    plaque.
  - The post is a fluted column with stepped flutes and a stepped capital.
  - The yaw hub is a stepped octagon, with ziggurat corbels and a skirt hiding the nut.
  - The risers step back like skyscraper setbacks under the bearing blocks, with rondels at
    the knees, triple banding and chevrons.
  - The tail fin is a scalloped rising-sun fan.
  - The pointer has a pierced quarter-sun bracket and a stepped crown.
- **Icosahedron.** Beaded edges and vertex studs make it read as cut crystal in a leaded frame.

The ornament is raised 1 mm on the show face (+y) of the yoke and fin, and engraved 0.8 mm into
the side that rests on the bed. That costs at most 4% of any cross-section, so there are no
weakening cuts.

## Parts

| part | filament | notes |
|---|---|---|
| `rotor` | clear PETG | One piece, 288 × 275 × 303 mm, printed **nose-down with no supports**: blade rulings are at 50°, and the generator measures 0.000 of the area as unsupported overhang. 4 perimeters, 2.0 mm blades. M8 nut traps at both hub ends. |
| `ico_inboard`, `ico_outboard` | clear PETG | Hollow icosahedron, circumradius 70 mm, with beaded edges clipped flush to the face it prints on. Joined by six tabs with M3 × 16 bolts. Print outer face down. The first ring of faces sits at 41.8°, just past the 45° rule, so print the outer walls slowly or support above 50°. |
| `comm_disc` | dark brown, **swap to brass/gold PETG at Z = 4.0 mm** | Ø220 mm. Everything above the dial face is relief, so the swap turns the words, symbols, ticks, sunburst, bezel and stepped hub brass. |
| `finial` ×2 | brass/gold PETG | Stepped spires that press over the two axle-end M8 nuts. The collar clears the stationary face by 1 mm. |
| `pointer` | dark brown | Print lying flat. The arrow tip is shaped across the thickness, so it reads as a point from the front, and stops just outside the number ring. Two M3 × 20 screws go through the plinth into the seat on `block_mouth`. |
| `block_nose`, `block_mouth` | dark brown | Octagonal 608 blocks with stepped caps. Print standing on the rotor-side face, where the bearing pocket opens. `block_mouth` carries the fixed sleeve, which has a hex tip for the paddle and a 608 in its tip. |
| `race_spacer` | dark brown | 12 × 36 mm tube inside `block_nose`. It carries the outboard nut's clamp to the bearing's inner race. |
| `paddle` | dark brown | Slides onto the sleeve hex and is held by an M3 set screw. 16.5 mm through-hole for the paddle stone (epoxy). |
| `yoke_nose`, `yoke_mouth` | dark brown | Print flat on the −y face. `yoke_mouth` carries the fan fin, which lies flush with that face so it prints on the bed. |
| `yaw_hub` | dark brown | Two 608s on the vertical rod. Print body-down. The ziggurat corbels under the pad use 4 × 4 mm steps, and the skirt has a 45° inner cone. M5 nuts drop into the hex channels from below. |
| `post`, `base` | dark brown | Ø40 × 280 mm fluted post with a Ø52 capital. Three-tier 200 mm base; the four Ø6.4 holes in the bottom tier's exposed corners take desk screws. |

**Hardware:**
- M8 threaded rod: 550 mm for the axle and 360 mm for the post, grade 8.8 or A2. Bending stress
  in the axle at a 25 m/s gust is about 170 MPa.
- 5 × 608-2RS bearings, plus M8 nuts and washers. One extra M8 nut acts as a jam nut under the
  disc finial.
- 4 × M5 × 30 bolts with nuts, for the yoke feet to the yaw hub.
- 4 × M5 × 35 screws, self-tapped into Ø4.2 pilots in the risers, for the bearing blocks.
- 6 × M3 × 16 (icosahedron), 2 × M3 × 20 (pointer), 1 × M3 × 10 (paddle set screw).
- Epoxy for the paddle stone.

**Stones:** 6–8 tumbled quartz pieces, about 15–25 mm, for the icosahedron, plus one 16 mm piece
for the paddle. Match their masses roughly so the container runs balanced.

**Assembly order (downwind end):**
1. Fit `block_mouth`.
2. Slide `ico_inboard` over the sleeve. Its open side then reaches the sleeve tip.
3. Fit the paddle and its set screw.
4. Load the stones.
5. Bolt on `ico_outboard`. Fit a nut inside the hub and a nut in the outer trap.
6. Fit `comm_disc` and its nut.
7. Fit the jam nut and press on the finial.

At the nose: nut, bearing, `race_spacer`, washer, nut, then the finial.

**Clearances** (swept surfaces of revolution, checked by `windscrew.py`, which asserts > 0.5 mm):
- Rotor to yoke: 15 mm.
- Rotor to nose block: 9.8 mm.
- Icosahedron to riser: 12 mm.
- Icosahedron to pointer: 8 mm.
- Icosahedron to sleeve: 2 mm (the designed annulus).
- Paddle to walls: about 4 mm.
- Dial relief to pointer: 3 mm.
- Nose finial to block: 1 mm.

Axle height is 554 mm and the rotor top is at 704 mm.

**Not yet printed or tested.** Before trusting the model, measure the cut-in wind and the rpm
with a phone tachometer.

**Storm limit:** at 25 m/s the model gives about 1300 rpm and 12 N of stone pinning force.
That is safe for the parts, but the stones will chew themselves to grit. Tie the rotor off in
gales.

**Sources:**
- [Chapman, *Image intensifier studies of low intensity triboluminescence…* (Open Univ. thesis)](https://oro.open.ac.uk/64651)
- [AVS 2012, gas discharge in quartz–diamond friction](https://www2.avs.org/symposium2012/Papers/Paper_TR-TuP2.html)
- Tian & Kennedy, *J. Tribology* 116:167 (1994)
- Hecht, Shlaer & Pirenne, *J. Gen. Physiol.* 25:819 (1942)
