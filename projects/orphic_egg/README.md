# orphic_egg

The Orphic Egg: a serpent wound three times round the world-egg, its head raised over the crown and
a skeleton key in its jaws.
Sculpted after the reference pendant (`reference_front.webp`, `reference_head_on.webp`,
`reference_above.webp`) **without its hanger loop**, so the crown of the egg is smooth. Because it
no longer hangs, the serpent's last turn slides off the egg onto the bed and its tail lies round the
egg's tip as a foot.

Three sculptures, each printable bare or fused onto its own plinth:

| File | What | Size |
| --- | --- | --- |
| `orphic_egg.stl` | the serpent on its egg | 103 × 103 × 143 mm, 514 cm³ |
| `orphic_serpent.stl` | the same serpent with no egg: the empty coil, as if the egg has hatched | 103 × 103 × 143 mm, 165 cm³ |
| `orphic_shedding.stl` | the empty coil shedding its skin | 108 × 110 × 143 mm, 197 cm³ |
| `orphic_egg_plinth.stl` | the egg sculpture on its moulded drum, one print | 149 × 149 × 165 mm, 875 cm³ |
| `orphic_serpent_plinth.stl` | the empty coil on its moulded drum, one print | 149 × 149 × 165 mm, 527 cm³ |
| `orphic_shedding_plinth.stl` | the shedding serpent on its drum, the Seed of Life cut 1 mm into the top | 149 × 149 × 165 mm, 556 cm³ |

```bash
python orphic_egg.py                # orphic_egg.stl, about 3 min
python orphic_egg.py --serpent      # orphic_serpent.stl, about 3 min
python orphic_shedding.py           # orphic_shedding.stl, about 14 min (--draft: 0.42 mm voxels, 2 min)
python plinth.py [name ...]         # *_plinth.stl from the STLs above (default all three), about a minute each
python orphic_egg.py --draft        # orphic_egg_draft.stl at 0.42 mm voxels, about 20 s
python preview_orphic_egg.py                                                 # orphic_egg_views.png, orphic_egg_eye.png
python preview_orphic_egg.py orphic_serpent.stl orphic_serpent_views.png -   # serpent views only
python preview_shedding.py          # orphic_shedding_views.png, orphic_shedding_details.png
```

## Proportions

Sized from the shared grid (`lib/sacred_grid.py`, R = 15 mm) like the rest of the matched set:

| Part | Size |
| --- | --- |
| Egg | 8R × 5R = 120 × 75 mm, a Fibonacci 8:5 (near φ), Hügelschäffer profile, blunt end up |
| Serpent | body Ø R = 15 mm, about 950 mm long; coils 2.25R apart along the meridian |
| Scales | 26 dorsal rows of imbricate diamonds (about 3.6 × 4.2 mm, 0.7 mm relief), belly plates on the egg side |
| Head | about 35 mm long and 21 mm wide, carried 20 mm above the crown |
| Eyes | 6 × 5.5 mm domed eyeball set 1 mm proud in a recessed socket, under a swept brow shelf the socket cuts into a lid edge |
| Pupils | an exact vesica piscis, 2.5 × 4.3 mm (√3 : 1), its horns 0.6 mm deep; its inscribed circle (Ø 2.5 mm) is a vortex, a disc falling to its centre with a 1.5-turn spiral groove (0.45 mm wide, 0.83 mm pitch) winding inward to a 1 mm-deep pit. The two eyes spin in mirror image |
| Key | 45 mm skeleton key bitten crosswise: bow on the serpent's right (toward the front view), a 16 mm medallion of the Tripod of Life (three circles, each through the others' centres) in a ring tangent to its lobes; Ø 3.6 mm shaft with a throat collar; warded bit hanging on the left; plates 2.6 mm thick |
| Plinth | a moulded round drum, top radius 4.5R (15 mm clear of the serpent's widest reach, 8 mm clear of the shedding skin's), 1.5R tall, 149 mm across at the foot. From the bed up: a chamfered foot, a plain band, a cove, a bead and a top band with a rounded edge |
| Seed of Life | on the shedding serpent's plinth only: the set's figure (seven circles, each through its neighbours' centres) at circle radius 2R, so it spans 120 mm and its outer arcs stop 7.5 mm inside the top edge; grooves 1.8 mm wide (the set's stroke), 1 mm deep, left uncut under the sculpture's footprint |

## How it is built

One signed-distance field: the egg (exact 2-D distance to its meridian), the serpent as a tube
round a centreline, with scale relief added as a displacement, and a head sculpted from ellipsoids,
with eyes, brow ridges, nostrils and a mouth slit, plus the key (extruded 2-D fields with rounded
edges). They are joined by smooth unions, so every coil is filleted onto the egg, the neck grows
out of the skull and the lips close round the key's shaft. The centreline is a helix at constant
meridian pitch, offset from the egg so the coil sinks 0.9 mm into it. Where the helix would pass
below the bed, a smooth max with the bed takes over. A Hermite curve carries the neck up to the
head; the tightest visible bend is 1.6 body radii. `--serpent` leaves the egg out of the field but
still uses it to shape the coils, so the serpent is identical in both prints.

Each eye is built in its own skin frame: a flattened ellipsoid eyeball, an ellipsoidal orbit carved
round it, and a brow shelf swept over it. The orbit is carved after the brow, so its rim becomes a
crisp upper lid. The vortex pupils are too fine for the body's 0.2 mm voxel, so each is meshed on
its own at 0.04 mm as a cutter solid and subtracted from the finished print with manifold3d.

It is meshed by marching cubes in z-slabs. The exact field is evaluated only within 2.4 mm of the
surface (a 1 mm coarse field, trilinearly interpolated, fills in elsewhere). The slab seams are
welded, then 8 passes of Taubin smoothing take the marching-cubes stair-steps off the cut edges
(pupil, lids, mouth, scale steps) without shrinking the volume (0.01 %). Finally the mesh is
quadric-decimated with pymeshlab to keep each STL under GitHub's 50 MB warning. Marching cubes
leaves a few hundred sub-voxel crumbs off the sharp scale steps; the build drops them and reports
the largest.

## The shedding serpent

`orphic_shedding.py` builds the empty coil again and lays its old skin over it in the same field.
The front half, from the head to about 1.6 turns down, shows the new scales. At the peel front,
450 mm along the body, the old skin has been pushed back: a rolled collar, then concertina folds,
and ahead of it a cuff torn into four flared, ruffled tongues. Behind the front the skin lies loose
over the rest of the body, crumpled, with sharper creases, and the crumpling grows where the tail lies on the
ground. Eighteen tears show the new scales underneath, each with a curled-up rim, and four of them
have their torn piece still hinged on and curling up like a lid. Three drapes of skin hang off
the coils' outer flanks: one onto the coil below, one onto the tail, and one down to the ground,
where its foot splays out. Four crumpled sheets lie on the ground beside the tail.

The old skin is built in the serpent's own tube coordinates (s along the centreline, φ round it,
ρ out from it) as an outer radius and, where it lifts, an inner one 1.3 mm in. So it is solid and
fused to the body where it lies down and a shell where it lifts. Each radius is divided by its
slope so the field stays a true distance on the steep folds. The drapes are shells in cylindrical
coordinates about the coil's axis, starting and ending inside a tube. The ground sheets are height
fields. All of the old skin carries the net of scale outlines a shed skin keeps (grooves 0.35 mm
deep on the same diamond lattice as the scales), so it reads differently from the overlapping
scales of the new skin.

The plinth is revolved exactly from its profile. Each sculpture's STL is set 0.5 mm into its top,
so no faces are coplanar, and unioned with manifold3d at full detail (48 MB each). The union
leaves one zero-volume scrap where the tail tip grazes the top plane; the build drops it.

## Printing

- **Upright, as exported.** `orphic_egg.stl` stands on the egg's tip (a disc 21 mm across) plus a
  1.14r-wide strip along the tail: 791 mm² on the bed within a 2745 mm² footprint. It leans about
  21° before tipping. `orphic_serpent.stl` stands on the tail alone: 460 mm² on the bed, also about
  21°. The plinth prints stand on the drum's full base.
- **Supports: tree.** About 18 % of the egg print's surface overhangs more than 45°: the underside
  of each coil, the head and chin, and the key's bow and bit, which hang at least 12 mm clear of the
  coils and egg. A small pocket between the egg's tip and the tail ring also needs support. The
  fillet where coil meets egg helps the coils' inner edge self-support.
- **The empty coil needs more.** Without the egg every coil overhangs along its whole underside
  (22 % of the surface), and the upper coils sit over the lower ones. Let tree supports land on the
  model, not just the build plate. On the plinth prints they can also land on the plinth top.
- **The shedding serpent** stands on the tail, the ground sheets and the foot of one drape: 1896 mm²
  on the bed within a 6000 mm² footprint, leaning about 28° before it tips. Only 16 % of its
  surface overhangs more than 45° (22 % for the bare coil), because the drapes are near-vertical
  walls that tie the lower coils to each other and to the ground. The cuff's tongues and the
  peeled flaps lift off the body with thin, 1.3 mm edges and need tree supports beneath them. Remove
  those supports gently. On the plinth print the Seed of Life grooves need nothing, and the
  drape's foot lies across them.
- **The key joins only through the jaws** (a Ø 3.6 mm shaft), so leave its supports on until the
  print has cooled.
- **The plinth needs no support.** The cove only steps inward, the bead's underside is a 45° bevel
  and the foot has a 0.6 mm chamfer against elephant's foot. It is 362 cm³ solid; 10 % infill is
  plenty there, and it makes the print far steadier on the bed than the sculpture's own foot.
- The vortex pupils are at the edge of what FDM can draw: the spiral's groove and land are both
  about one nozzle wide. On a 0.4 mm nozzle at 0.12 mm layers, with a slow outer wall, the vesica
  and its circle read clearly and the spiral reads as a soft whorl. A 0.2 mm nozzle or resin
  shows it crisp.
- Scales are 0.7 mm relief on 3.6 mm diamonds: fine on a 0.4 mm nozzle at 0.12–0.16 mm layers. On
  the top faces of the coils they will show layer stepping.
- Use 10–15 % infill and 3 walls; the scale relief lives entirely in the walls. Without the egg
  the serpent is a free coil spring of Ø 15 mm PETG; a coil-spring estimate puts its sag under its
  own weight below 1 mm.
