# From One Circle — the matched set

Every piece of the sacred-geometry set is sized from one base circle, **R = 15 mm**, so the prints
match each other. The geometry lives in [lib/sacred_grid.py](../../lib/sacred_grid.py); the pieces are
built by their own projects; this directory gathers them for the site page
[web/from-one-circle.html](../../web/from-one-circle.html).

## The three sizing rules

| Family | Rule | Result |
| --- | --- | --- |
| Plates | Circles of radius R on one triangular lattice of spacing R, one direction vertical; common 170 mm ring, 3 pins at 0°/120°/180° | Every plate stacks in register; Metatron's Cube is 10R = 150 mm |
| Solids | Seen down a three-fold axis, the outermost corners land on the circle of radius 4R = 60 mm through Metatron's Cube's outer centres | Tetrahedron, cube, octahedron: every corner on a Metatron centre. Dodecahedron 8/20, icosahedron 6/12 (φ) |
| Spheres | Radius R on a cubic grid of spacing R·√(3/2) = 18.37 mm | Corner-on shadow is the plate lattice: 2³ casts the Seed, 3³ the Flower |

Stellations grow from the grid solids (apex where neighbour face planes cross), so their cores are
exactly the solids; the tetrahedron and cube get their dual compounds. The five tetrahedra sit on the
grid dodecahedron's corners.

## Where each piece is built

| Project | Script | Files |
| --- | --- | --- |
| [sacred_geometry_plates](../sacred_geometry_plates) | `plates.py` | `01_circle.stl` … `08_tree_of_life.stl`, `plate_pin.stl` |
| [platonic_solid_frames](../platonic_solid_frames) | `frames.py` | `<solid>_frame.stl` ×5 |
| [stellated_platonic_solids](../stellated_platonic_solids) | `stellations.py` | `stellated_<solid>_frame.stl` ×5 |
| [five_intersecting_tetrahedra](../five_intersecting_tetrahedra) | `five_intersecting_tetrahedra.py` | woven + rigid |
| [flower_of_life_spheres](../flower_of_life_spheres) | `spheres.py` | `egg_of_life.stl`, `flower_of_life_spheres.stl` |

Every script checks its output: one body, watertight, consistent winding, inside the 320 mm bed.
The superseded files of each project are in its `_legacy/` folder.

## Rebuilding the page assets

Run each project script from its own directory, then:

```bash
cd projects/from_one_circle && python web_assets.py
```

That writes `site/assets/one-circle/` (viewer meshes `.nsm`, studio renders `.webp`, `geometry.js` with
the figures and corner sets for the page's SVGs), copies the STLs to `site/stl/one-circle/`, and zips
the whole set with printing notes into `site/stl/one-circle-set.zip`.

The menu preview `site/previews/from-one-circle.webp` was cropped from a browser screenshot of the
solids section: `make-previews.ps1` screenshots `web/` as `file://`, where `site/assets/` does not
resolve, and headless Edge ignores the `#anchor`.
