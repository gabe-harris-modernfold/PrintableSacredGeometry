# Dodecahedron Mirror Ball

A 12-piece printable frame that turns twelve **10 cm hexagonal mirror tiles**
(the common 20 × 17.32 × 10 cm spec — those three numbers are across-corners,
across-flats, and side of one regular hexagon) into a mirrored regular
dodecahedron, Ø300.6 mm point-to-point.

The hexagons must be cut down to pentagons, and that is not a limitation of
this design — it is forced by the geometry. **No polyhedron can be built from
hexagons alone.** A regular hexagon's interior angle is 120°, so three around a
vertex sum to 360°: flat, the honeycomb tiling, no corner. Descartes' theorem
requires the angular defects at the vertices of a convex polyhedron to total
720°, and hexagons contribute 0° each. Euler says the same thing: for F
hexagonal faces, E = 3F and V ≤ 2F, so χ ≤ 0, never the 2 a sphere needs.
Twelve hexagons land at χ = 0 exactly — torus topology, which is the shape they
actually want to be.

Twelve is still the magic number here, just for pentagons: each pentagon vertex
in a soccer ball carries 360 − (108 + 120 + 120) = 12° of defect, and
12 pentagons × 5 vertices × 12° = 720°. That is why *every* hexagon-based
ball — soccer ball, geodesic dome, fullerene — has exactly twelve pentagons no
matter how large it grows. Here the pentagons are the whole solid.

Each mirror is inset 5 mm from its face edge, so the frame reads as a 10 mm rib
between facets. That inset is deliberate: 4 mm glass cut to the exact face
pentagon cannot butt at a 116.565° dihedral without bevelled edges, because two
square-cut slabs meeting there overlap by 180 − 116.565 = 63.435°.

## Parts

| File | Part | Print notes |
|---|---|---|
| `dodeca-mirror-tile.stl` | Bezel tile, 173.6 × 165.1 × 14 mm — **print 12** | **As oriented, pocket up.** Support-free: both walls lean out only 31.7° from vertical and the ledge grows off the bed rather than bridging. First layer is 3839 mm² over an 8.65 mm-wide ring, so **no brim needed** — bed adhesion is generous. |
| `dodeca-mirror-pin.stl` | Alignment pin, Ø4 × 30 mm — **print 11, not 30** | Or substitute 4 mm dowel or brass rod. Only 11 of the 30 joints can take a pin — see Assembly. Groove is Ø4.3, so 0.15 mm radial; FDM channels print undersize, so ease the pins with a file if they bind. |
| `dodeca-mirror-assembly.stl` | Finished-ball preview, frames + glass | **Reference only, do not print.** |
| `pentagon-cut-template.pdf` | 1:1 cutting template, A4/Letter landscape | Print at **100% / Actual Size**, not "fit to page". Verify the 100 mm scale bar with a ruler before cutting glass. |

Tiles mate on mitre faces at **58.283°**, half the 116.565° dihedral. All twelve
are the same part: alignment is a half-round groove down each mitre face, so two
mating tiles form one Ø4 mm bore. That works because a dodecahedron has a 2-fold
rotation axis through every edge midpoint, and a groove parallel to the edge and
centred on the mitre plane maps onto itself under that rotation.

Bed fit: 173.6 mm bare clears a 180 mm bed by 6.4 mm and anything larger
comfortably. A brim would need 200 mm, but the footprint makes one unnecessary.

## Bill of materials

- 12 printed tiles, 11 pins — 553 cm³ solid, so ~686 g at 100% infill but
  realistically **300–350 g** at 3 perimeters and 15% infill
- 12 hexagonal mirror tiles, 10 cm side. Order 15: twelve finished pentagons
  means 60 scores, and mirror backing likes to flake at the corners
- Neutral-cure silicone to bed the glass (**required** — see below)
- Glass cutter, straightedge, safety glasses

Finished ball: 238.9 mm face-to-face, 300.6 mm point-to-point. Glass weighs
**2.06 kg at 4 mm**, 1.55 kg at 3 mm.

There is **no hanging or mounting provision** — this is a bare shell. The
central openings give interior access if you want to bond an eye bolt to a
tile's inner face; otherwise plan on a stand.

## Cutting the mirrors

The largest regular pentagon that fits in a 100 mm regular hexagon has a
**102.37 mm** side, achieved with a pentagon corner pointed at a hexagon corner
(24° of rotation). The frame is designed for **100 mm**, and the pocket accepts
anything from 100.00 up to **101.74 mm** — that +1.74 mm is the margin for
hand-cutting, and it is the tightest tolerance in the whole design, because all
five edges have to fit at once. Aim slightly under rather than over; undersize
just leaves a marginally wider gap, oversize means the mirror will not seat. If
one edge fouls, grind it back with a diamond file rather than recutting.

Trace the template, score all five chords edge-to-edge, snap each one. You keep
about 69% of each hexagon's area.

If the tiles have a factory bevel or polished edge, it will read as a groove line
that your fresh cuts won't match. This design is unbothered — all 60 edges are
new cuts, so the finish stays uniform.

## Assembly

1. Bed each mirror on silicone in its pocket, reflective side out, and let it
   cure. There is **no retaining lip** over the front, and there cannot be: the
   mirror is larger than the tile's central opening, so it can only be loaded
   from the front. Silicone is the right choice against glass regardless — it
   absorbs the differential expansion that rigid adhesive would transfer as
   stress into the pane. With 1.2 mm of pocket clearance the glass can shift, so
   centre each pane by eye before the silicone skins over.
2. **The frame self-jigs, so pins are optional.** Every tile sits in a tapered
   socket: measured against the assembled ball, a tile lifts straight out
   radially with zero interference at any distance, and pushing it inward
   interferes immediately (4.0 mm³ at 5 µm). It can only enter from outside,
   and it wedges on the mitre faces when it seats. The last tile is a keystone.
3. **A pin only works on a joint you close by pure normal translation** —
   two tiles brought face-to-face on the bench. Verified clean at every
   separation along the mitre normal.
4. **A pin cannot be used to seat a tile into the standing assembly.** The
   wedge motion is radial, the groove axes are tangential, and the two are
   exactly perpendicular, so the grooves shear across the pins: binding starts
   at 0.25 mm of travel and is solid by 0.5 mm. Since a tile must travel many
   millimetres to seat, pins are unusable there. This is fundamental rather
   than a flaw in the groove — *any* key crossing a mitre face resists sliding
   in that face, and wedge closure is exactly that sliding. A key aligned to
   the insertion direction would not survive the edge's 2-fold rotation, so it
   would need mirrored left/right tiles instead of 12 identical ones.
5. So: pick a build order where each new tile touches exactly one already-placed
   tile — a spanning tree, **11 of the 30 joints** — and pin those, unrolled
   like a papercraft net. The remaining 19 joints close on the final fold, glue
   only. Or skip pins entirely and let the taper do the work; the grooves then
   serve as glue reservoirs, which is worth having either way.
6. Glue the mitre faces progressively, closing the ball last. The final tile is
   a captive fit — leave it dry-fitted if you want the interior accessible.

## Sizing constraints (edit the parameter block in `dodeca_frame.py`, re-run)

Two of these are coupled and will silently produce an unprintable part if you
change one alone. The generator prints and grades both on every run.

- `MIRROR_THK` (4.0) — **verify this against your actual glass.** Many
  decorative hex tiles are 3 mm, which would leave every mirror standing 1 mm
  proud of the rib.
- `CLEAR` (1.2) — radial pocket clearance, and therefore the cutting tolerance.
  Was 0.6, which allowed only +0.87 mm of oversize; too tight for glass cut by
  hand to a traced template, where ±1 mm per edge is normal.
- `RIB` (5.0) — visible frame width per side; the rib between facets is 2×.
  **Must satisfy `RIB > CLEAR + MIRROR_THK × face_inr/r_i + ~1.2`.** The mitre
  plane leans inward at 0.618 mm per mm of depth, so the pocket wall is a
  *wedge*: 3.80 mm at the mirror face but only 1.33 mm where it meets the ledge,
  having lost 2.47 mm over the 4 mm pocket. That narrow end is what the slicer
  prints first. At `RIB` 4.0 with `CLEAR` 1.2 it fell to 0.33 mm — under one
  extrusion width, so the first traces of the rim would be dropped or broken.
- `THK` (14.0) — frame thickness. Below ~12 the groove walls get thin.
- `MIRROR_SIDE` (100.0) — rescales the entire ball.

The central ledge cut is a pyramid from the body centre, not a prism, so it
tapers at the same rate as the shell and holds a constant ledge width at every
depth. A prism gets overtaken by the mitre taper and the ledge vanishes.

## Validation

`dodeca_frame.py` asserts on every run; `verify.py` does the pin-fit test;
`audit.py`, `audit2.py` and `audit3.py` are the geometric, assembly and
dimensional audits.

| Check | Result |
|---|---|
| Tile watertight, 1 body, genus 1 | True / 1 / 1 |
| Frame ∩ mirror | 0.00 mm³ |
| Tile ∩ neighbouring tile | 0.010 mm³ |
| Gap between mating mitre faces | < 5 µm (contact is exact) |
| Radial withdrawal, 0.25–10 mm | 0.00 mm³ at every step |
| Adjacent mirror clearance | 4.30 mm |
| Groove ∩ mirror pocket | 0.00 mm³; 2.85 mm of wall each side |
| Pocket wall at the ledge | 1.33 mm = 3.2 traces |
| Mirror support ledge | 8.14 mm |
| First-layer footprint | 3839 mm² |
| Pin ∩ any solid, seated | 0.000 mm³ |
| Tiles gripping each pin | exactly 2, balanced 152/152 mm³ |
| Pinned pair separated along −n | clean at every distance |
| Pin vs radial insertion | **binds at 0.25 mm — see Assembly** |
| Symmetry residual | 2.1e-13 mm |

Three traps worth recording, since each produced a confident-looking pass over a
real defect:

**Summing signed volumes proves nothing about overlap.** For separate closed
bodies it returns the sum either way. An early version reported a perfect
`12 × (tile + mirror)` match while pin bosses were driven 2086 mm³ into the
mirror pockets. Only real boolean intersections catch interference.

**Boolean union is not a gap test.** Unioning the 12 tiles reports 34 components
and euler 46 — not a valid closed-surface value, so it is a boolean artifact on
coincident faces, not a gap. Manifold declines to fuse face-touching solids at
all: a pairwise union of two mating tiles also stays 2 components, at exactly
2× the tile volume. Contact was instead proved by pushing one tile inward and
watching interference rise linearly from zero (4.0 mm³ at 5 µm), which bounds
the gap below 5 µm.

**A wall's nominal width is not the width that gets printed.** Every wall bounded
by a mitre face is a wedge, because the mitre leans in at 0.618 mm/mm. Quoting
`RIB − CLEAR` for the pocket wall overstated it by 2.47 mm and hid a
0.33 mm sliver. Cross-sectioning through the print height is what exposed it:
the ring collapsed from 4068 mm² to 230 mm² across the pocket floor.
