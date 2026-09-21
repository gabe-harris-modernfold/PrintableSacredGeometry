# Nickel-box carry handle

Slip-on end cap for a bank box of nickels (240 x 88 x 100 mm, ~10 kg for a $100 box), with a
flared D-handle on top. Print **two**; one on each end turns the box into a two-handled carry /
row / farmer's-walk weight. The caps just push on; nothing clamps.

| file | what |
|---|---|
| `nickel_box_handle.scad` | parametric source (OpenSCAD) -- box size, clearance, wall, handle sizes at the top |
| `nickel_box_handle.stl` | one cap, already in print orientation (end plate on the bed) |
| `render_previews.py` | regenerates the three preview PNGs from the STL in seconds |
| `preview_end.png`, `preview_part.png`, `preview_assembly.png` | renders |

## Design

The part is an extruded end-view profile: rectangular sleeve + D-loop, closed by a 4 mm end plate.
The handle is flush with the end face and 30 mm thick, so in print orientation every feature rises
straight off the bed -- **no supports, no bridges**. In-use loads act in the layer plane.

- cavity 91 x 103 mm (1.5 mm clearance per side), sleeve reaches 70 mm along the box, 3.6 mm walls
- every outside edge rounded 3 mm (stepped-offset fillet, 6 x 0.5 mm slices: OpenSCAD 2021's CGAL
  minkowski() took 10 min on this profile, this takes ~2.5 min), 1.5 mm 45-degree lead-in at the mouth
- four scalloped finger grooves (22 mm pitch, 3 mm deep) under the 24 mm grip bar
- two paracord ears at the lower corners of the hoop, 7 mm holes running along the box: lash the two
  caps together along both sides so they cannot walk off the box
- finger opening 44 mm tall, 112 mm wide under a 22 mm grip bar; uprights flare so the grip is
  wider than the box
- the 30 mm handle section is a full hoop, 6 mm thicker than the sleeve, running all the way round the box, so the grip load is carried around and under the box, not just through the top wall
- gusset tapers the loop back 25 mm along the top wall

Part envelope 190 x 154 x 70 mm, watertight, ~483 cm^3 (~610 g PETG at 100% walls; use 3-4
perimeters and 30% infill and it drops well under 300 g).

## Rebuild

```
openscad -o nickel_box_handle.stl nickel_box_handle.scad
python render_previews.py
```

Measure your actual box first and set `BOX_W` / `BOX_H` (end face width / height); swap them if you
want the handle on the other face.
