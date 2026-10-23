# five_intersecting_tetrahedra

The chiral compound of five tetrahedra as woven, flat-faced beam frames (Francis Ow's 60° unit), on
the grid dodecahedron of the [matched set](../from_one_circle): tetrahedron edge 103.9 mm, the grid
cube's face diagonal.

```bash
python five_intersecting_tetrahedra.py
```

| Output | |
| --- | --- |
| `five_intersecting_tetrahedra.stl` | five free interlocked frames, 0.8 mm clearance, beam width L/14.0 |
| `five_intersecting_tetrahedra_rigid.stl` | the same compound welded into one solid, beam width L/9.5 |
| `five_intersecting_tetrahedra.png` | preview, one colour per frame |

Both are 121 × 125 × 101 mm, standing on a dodecahedron face (five tips on the bed). The jamming
width m* = L/12.09 is solved by LP and matches the origami unit's L/12. `_legacy/` holds the 150 mm
originals.
