# sacred_geometry_plates

The numbered 2D plate series, from a single circle up to the icosahedron. All shapes share one base
radius `R` so the plates overlap exactly when stacked.

`Sacred_Geometry_Generator.scad` is the source. It is **not** a batch exporter: pick a shape by
uncommenting its module call at the bottom of the file, render with F6, then
`File → Export → Export as STL`. Key parameters at the top: `R` (base circle radius, 10 mm),
`T` (extrusion thickness, 2 mm), `$fn` (curve resolution, 64).

`00_alignment_preview.png` shows the whole series on a common grid.

## Note on the numbering

The prefixes come from two different generator runs, so a few names repeat at different numbers
(`05_fruit_of_life.stl` / `06_fruit_of_life.stl`, `06_metatrons_cube.stl` / `07_metatrons_cube.stl`,
`07_tree_of_life.stl` / `08_tree_of_life.stl`). The paired files are **not** identical meshes — both
generations were kept rather than guessing which is current. The `_2d` suffix marks the flat
extruded plate version of a solid that also exists as a 3D form.
