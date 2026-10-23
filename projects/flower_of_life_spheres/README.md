# flower_of_life_spheres

The Egg of Life and the Flower of Life as clusters of overlapping spheres, for the
[matched set](../from_one_circle).

```bash
python spheres.py
```

Spheres of radius R = 15 mm sit on a cubic grid of spacing R·√(3/2) = 18.37 mm. Seen down the body
diagonal, that grid's shadow is the plates' triangular lattice of spacing R, so:

| Output | Spheres | Casts | Print |
| --- | --- | --- | --- |
| `egg_of_life.stl` | 2 × 2 × 2 | the Seed of Life (7 circles) | 48 × 48 × 47 mm |
| `flower_of_life_spheres.stl` | 3 × 3 × 3 | the Flower of Life (19 circles) | 67 × 67 × 65 mm |

At this spacing each cubic cell's centre sits just outside all eight spheres, a sealed 0.9 mm bubble;
a small filler sphere closes it, so each cluster is one solid body. Seated on a cube face with a 1.5 mm
footprint; supports under the lower spheres.

`_legacy/` holds the earlier exploratory generator (hex / fcc / cubic modes) and its 4 × 4 × 4 print.
