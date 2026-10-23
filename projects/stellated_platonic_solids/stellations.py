#!/usr/bin/env python3
"""First stellation of each Platonic solid, as an edge frame on the shared grid.

Each one is grown from the grid-sized solid in projects/platonic_solid_frames, so the
core of every star is exactly its solid. The method is the one derived in
_legacy/stellated_platonics.py: the pyramid on each face has its apex where the planes of
the face's edge-neighbours cross its normal.

    octahedron    -> stella octangula              true stellation
    dodecahedron  -> small stellated dodecahedron  true stellation
    icosahedron   -> small triambic icosahedron    true stellation
    tetrahedron   -> stella octangula              dual compound (it does not stellate)
    cube          -> cube + octahedron             dual compound (it does not stellate)

Struts that continue one another in a straight line are merged, so the star forms read as
long lines through their core corners. Geometry lives in lib/sacred_grid.py.

Run from this directory: python stellations.py
  stellated_<solid>_frame.stl   seated on its widest hull face, 0.6 mm footprint
  preview_stellations.png
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import numpy as np

import sacred_grid as sg


def main():
    ok, items = True, []
    for name in sg.SOLIDS:
        V, E, true = sg.stellation(name)
        man = sg.frame(V, E)
        mesh = sg.to_trimesh(sg.print_pose(man, V))
        path = f"stellated_{name}_frame.stl"
        ok &= sg.check(mesh, path)
        mesh.export(path)
        print(f"    {sg.STELLATED[name]}; tips on a {2 * np.linalg.norm(V, axis=1).max():.1f} mm sphere, "
              f"{len(E)} struts")
        items.append((f"{name} -> {sg.STELLATED[name].split(' (')[0]}", sg.to_trimesh(man)))
    sg.contact_sheet(items, "preview_stellations.png", cols=5, az=20, el=18, tile=360)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
