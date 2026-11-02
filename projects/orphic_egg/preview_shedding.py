#!/usr/bin/env python3
"""Headless renders of the shedding serpent: the standard five views, and a sheet of close-ups
of the peel front, the drapes and the skin slumped round the tail.

Run from this directory:  python preview_shedding.py [orphic_shedding.stl]
    orphic_shedding_views.png, orphic_shedding_details.png
"""

import sys

import numpy as np
import pyvista as pv

import orphic_egg as oe
import orphic_shedding as osh
from preview_orphic_egg import VIEWS, shoot


def main():
    stl = sys.argv[1] if len(sys.argv) > 1 else "orphic_shedding.stl"
    mesh = pv.read(stl)
    b = mesh.bounds
    shots = []
    for title, d, dist, focus in VIEWS:
        c = np.array([(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, b[4] + focus * (b[5] - b[4])])
        if title == "head":
            c = mesh.points[mesh.points[:, 2] > b[5] - 14].mean(axis=0) - [0, 0, 3]
        shots.append((title, c, d, dist))
    shoot(mesh, shots, "orphic_shedding_views.png")

    body = oe.Piece(with_egg=False).body
    up = np.array([0.0, 0.0, 1.0])

    def on_body(s):
        i = int(np.argmin(np.abs(body.s - s)))
        c = body.C[i]
        return c, body.T[i], np.r_[c[:2] / np.linalg.norm(c[:2]), 0.0]

    c, t, out = on_body(osh.FRONT)
    details = [("peel front, from ahead", c, -t + 0.6 * out + 0.5 * up, 75),
               ("peel front, from behind", c, t + 0.6 * out + 0.5 * up, 75)]
    for (s0, s1, to), title in zip(osh.CURTAINS, ("drape onto the coil below", "drape onto the tail",
                                                   "drape to the ground")):
        c, _, out = on_body((s0 + s1) / 2)
        c = c.copy()
        c[2] = c[2] / 2 if to == "ground" else c[2] - 17.0
        details.append((title, c, out + 0.15 * up, 120))
    details.append(("tail and sheets", np.array([-20.0, 0.0, 8.0]), np.array([-1.0, 0.0, 0.45]), 150))
    shoot(mesh, details, "orphic_shedding_details.png", size=480, height=480)


if __name__ == "__main__":
    main()
