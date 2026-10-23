#!/usr/bin/env python3
"""The Egg of Life and the Flower of Life as sphere clusters on the shared grid.

Spheres of radius R = 15 mm sit on a cubic grid of spacing R*sqrt(3/2). Seen straight down
the cube's body diagonal, that grid's shadow is the plates' triangular lattice of spacing R,
and each sphere's shadow is a circle of radius R, so

    2 x 2 x 2  =  8 spheres  ->  casts the Seed of Life (7 circles; two spheres share the centre)
    3 x 3 x 3  = 27 spheres  ->  casts the Flower of Life (19 circles)

Set either on its plate, corner down, and its outline covers the plate's circles. The
spheres overlap their neighbours by 39% of a diameter, so each cluster is one solid.
Geometry lives in lib/sacred_grid.py; the exploratory hex/fcc/cubic generator this
replaces is in _legacy/.

Run from this directory: python spheres.py
  egg_of_life.stl, flower_of_life_spheres.stl   seated on a cube face, 1.5 mm footprint
  preview_spheres.png                           from above the diagonal, and at rest
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import sacred_grid as sg

PIECES = {"egg_of_life": 2, "flower_of_life_spheres": 3}


def main():
    ok, items = True, []
    for name, n in PIECES.items():
        man = sg.spheres(n)
        mesh = sg.to_trimesh(sg.print_pose(man, sg.sphere_centres(n), flat=1.5))
        path = f"{name}.stl"
        ok &= sg.check(mesh, path)
        mesh.export(path)
        items.append((f"{name}: down the diagonal", sg.to_trimesh(sg.grid_view(name))))
        items.append((f"{name}: at rest", mesh))
    print(f"    spheres R = {sg.R} mm, spacing {sg.SPHERE_STEP:.2f} mm")
    sg.contact_sheet([items[0], items[1], items[2], items[3]], "preview_spheres.png", cols=4,
                     az=-90, el=89.9, tile=360)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
