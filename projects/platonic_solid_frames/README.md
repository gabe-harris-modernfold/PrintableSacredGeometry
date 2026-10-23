# platonic_solid_frames

The five Platonic solids as edge frames (5 mm struts, 7 mm corner nodes), sized on the shared grid of
the [matched set](../from_one_circle): seen down a three-fold axis, each solid's outermost corners land
on the 120 mm circle through Metatron's Cube's outer centres.

```bash
python frames.py
```

| Solid | Edge | Corner sphere | Corners on a Metatron centre | Print |
| --- | --- | --- | --- | --- |
| tetrahedron | 103.9 mm | 127.3 mm | 4 / 4 | 111 × 97 × 91 mm |
| cube | 73.5 mm | 127.3 mm | 8 / 8 | 111 × 111 × 80 mm |
| octahedron | 103.9 mm | 147.0 mm | 6 / 6 | 111 × 127 × 91 mm |
| dodecahedron | 45.4 mm | 127.3 mm | 8 / 20 | 127 × 132 × 108 mm |
| icosahedron | 64.2 mm | 122.2 mm | 6 / 12 | 111 × 127 × 103 mm |

The cube and octahedron at these sizes are exact duals (edges cross at midpoints). The dodecahedron
shares the cube's 8 corners. The dodecahedron's and icosahedron's other corners cannot reach the
lattice: they carry φ, and Metatron's Cube has rings only at 2R and 4R.

Each STL is seated on its widest face with a 0.6 mm flat trimmed off for bed contact. Print with tree
supports and a brim.
