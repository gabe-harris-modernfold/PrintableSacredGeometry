#!/usr/bin/env python3
"""The five Platonic solids as edge frames, sized on the shared grid.

Each solid is scaled so that, seen down a three-fold axis, its outermost corners land on
the circle through Metatron's Cube's outer six centres (radius 4R = 60 mm). The
tetrahedron, cube and octahedron then sit exactly on lines of the Metatron's Cube plate;
the dodecahedron and icosahedron carry the golden ratio and only touch that circle.
Struts 5 mm, corner nodes 7 mm. Geometry lives in lib/sacred_grid.py.

Run from this directory: python frames.py
  tetrahedron_frame.stl ... icosahedron_frame.stl   seated on a face, 0.6 mm footprint
  preview_frames.png
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import numpy as np

import sacred_grid as sg


def main():
    ok, items = True, []
    for name in sg.SOLIDS:
        V, E = sg.solid(name)
        mesh = sg.to_trimesh(sg.print_pose(sg.frame(V, E), V))
        path = f"{name}_frame.stl"
        ok &= sg.check(mesh, path)
        mesh.export(path)
        r = np.linalg.norm(V, axis=1).max()
        edge = min(np.linalg.norm(V[i] - V[j]) for i, j in E)
        print(f"    edge {edge:.1f} mm, corners on a {2 * r:.1f} mm sphere, "
              f"shadow circle {2 * np.hypot(V[:, 0], V[:, 1]).max():.1f} mm")
        items.append((name, sg.to_trimesh(sg.frame(V, E))))
    sg.contact_sheet(items, "preview_frames.png", cols=5, az=20, el=18, tile=360)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
