# sacred_geometry_plates

The eight stackable plates of the [matched set](../from_one_circle): circle, vesica piscis, tripod of
life, seed, flower and fruit of life, Metatron's Cube, Tree of Life. Every circle has radius R = 15 mm
and every centre sits on one triangular lattice, so stacked plates put their lines on one another.

```bash
python plates.py
```

| Output | |
| --- | --- |
| `01_circle.stl` … `08_tree_of_life.stl` | 170 × 170 × 3 mm, one body each, flat, no supports |
| `plate_pin.stl` | D-section 3 mm pin, 24 mm long; print three, lying on the flat |
| `preview_plates.png` | all eight from above |

The common ring frame runs from 5R = 75 mm (where the Fruit of Life's outer circles end) to 85 mm. Where
a figure stops short of the ring, spokes along the lattice rays tie it in. Pin holes are 3.4 mm at 0°,
120° and 180°: no rotation or flip maps that set onto itself, so a plate seats one way only.

The Tree of Life is the Golden Dawn / Kircher layout (10 sephirot, 22 paths) on the Flower of Life
lattice with Tiferet at the centre; Daath is drawn as a ring.

`_legacy/` holds the earlier OpenSCAD generator and its two generations of mixed-size plates and 80 mm
solids, with their own README.
