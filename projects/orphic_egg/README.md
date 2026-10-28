# orphic_egg

The Orphic Egg: a serpent wound three times round the world-egg, its head raised over the crown and
a skeleton key in its jaws.
Sculpted after the reference pendant (`reference_front.webp`, `reference_head_on.webp`,
`reference_above.webp`) **without its hanger loop**, so the crown of the egg is smooth. Because it
no longer hangs, the serpent's last turn slides off the egg onto the bed and its tail lies round the
egg's tip as a foot.

```bash
python orphic_egg.py            # orphic_egg.stl, about 3 min
python orphic_egg.py --draft    # orphic_egg_draft.stl at 0.42 mm voxels, about 20 s
python preview_orphic_egg.py    # orphic_egg_views.png
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
| Print | 103 × 103 × 143 mm, one watertight body, about 935 k triangles (about 47 MB) |

## How it is built

One signed-distance field: the egg (exact 2-D distance to its meridian), the serpent as a tube
round a centreline, with scale relief added as a displacement, and a head sculpted from ellipsoids,
with eyes, brow ridges, nostrils and a mouth slit, plus the key (extruded 2-D fields with rounded
edges). They are joined by smooth unions, so every coil is filleted onto the egg, the neck grows
out of the skull and the lips close round the key's shaft. The centreline is a helix at constant
meridian pitch, offset from the egg so the coil sinks 0.9 mm into it. Where the helix would pass
below the bed, a smooth max with the bed takes over. A Hermite curve carries the neck up to the
head; the tightest visible bend is 1.6 body radii.

Each eye is built in its own skin frame: a flattened ellipsoid eyeball, an ellipsoidal orbit carved
round it, and a brow shelf swept over it. The orbit is carved after the brow, so its rim becomes a
crisp upper lid. The vortex pupils are too fine for the body's 0.2 mm voxel, so each is meshed on
its own at 0.04 mm as a cutter solid and subtracted from the finished print with manifold3d.

It is meshed by marching cubes in z-slabs. The exact field is evaluated only within 2.4 mm of the
surface (a 1 mm coarse field, trilinearly interpolated, fills in elsewhere). The slab seams are
welded, then 8 passes of Taubin smoothing take the marching-cubes stair-steps off the cut edges
(pupil, lids, mouth, scale steps) without shrinking the volume (0.01 %). Finally the mesh is
quadric-decimated with pymeshlab to keep the STL under GitHub's 50 MB warning. Marching cubes leaves a few hundred sub-voxel crumbs off the sharp scale steps; the build
drops them and reports the largest.

## Printing

- **Upright, as exported.** The flat at z = 0 is the egg's tip (a disc 21 mm across) plus a 1.14r-wide
  strip along the tail: 787 mm² on the bed within a 2747 mm² footprint. It leans about 21° before
  tipping.
- **Supports: tree, from the build plate only.** About 18 % of the surface overhangs more than 45°:
  the underside of each coil, the head and chin, and the key's bow and bit, which hang at least
  12 mm clear of the coils and egg. The key joins the piece only through the jaws (a Ø 3.6 mm
  shaft), so let its supports stay on until the print has cooled. A small pocket between the egg's tip and the
  tail ring also needs support. The fillet where coil meets egg helps the coils' inner edge
  self-support.
- The vortex pupils are at the edge of what FDM can draw: the spiral's groove and land are both
  about one nozzle wide. On a 0.4 mm nozzle at 0.12 mm layers, with a slow outer wall, the vesica
  and its circle read clearly and the spiral reads as a soft whorl. A 0.2 mm nozzle or resin
  shows it crisp.
- Scales are 0.7 mm relief on 3.6 mm diamonds: fine on a 0.4 mm nozzle at 0.12–0.16 mm layers. On
  the top faces of the coils they will show layer stepping.
- Solid it is about 513 cm³. Use 10–15 % infill and 3 walls; the scale relief lives entirely in the
  walls.
