# stellated_platonic_solids

First stellation of each Platonic solid, as an edge frame grown from the grid-sized solid in
[platonic_solid_frames](../platonic_solid_frames), so every star's core is exactly that solid.

```bash
python stellations.py
```

| From | Result | Kind | Tips on | Print |
| --- | --- | --- | --- | --- |
| tetrahedron | stella octangula | dual compound | 127.3 mm | 97 × 97 × 80 mm |
| cube | cube + octahedron | dual compound | 147.0 mm | 134 × 127 × 110 mm |
| octahedron | stella octangula | true stellation | 254.6 mm | 215 × 215 × 153 mm |
| dodecahedron | small stellated dodecahedron | true stellation | 226.2 mm | 211 × 228 × 186 mm |
| icosahedron | small triambic icosahedron | true stellation | 130.2 mm | 133 × 135 × 124 mm |

The method is the one derived in `_legacy/stellated_platonics.py`: the pyramid on each face has its
apex where the neighbouring face planes cross its normal, and a solid whose dihedral angle is 90° or less
(tetrahedron, cube) has no stellation, so it gets its midsphere dual compound instead. Collinear struts
are merged, so the stars read as long straight lines through their core corners.

`_legacy/` holds the earlier generator, its hollow-shell solids at 150 mm, and the combined plate.
