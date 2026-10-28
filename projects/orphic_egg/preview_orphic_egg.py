#!/usr/bin/env python3
"""Headless renders of the Orphic Egg: the reference pendant's three viewpoints (front,
head-on, from above) plus the back and the head, and a second sheet of close-ups of the
eye and its vortex pupil.

Run from this directory:  python preview_orphic_egg.py [orphic_egg.stl] [views.png] [eye.png]
"""

import sys

import numpy as np
import pyvista as pv

import orphic_egg as oe

STL = sys.argv[1] if len(sys.argv) > 1 else "orphic_egg.stl"
OUT = sys.argv[2] if len(sys.argv) > 2 else "orphic_egg_views.png"
EYE_OUT = sys.argv[3] if len(sys.argv) > 3 else "orphic_egg_eye.png"

#          title                 camera direction        distance  focus (fraction of height)
VIEWS = [("front",               (0.0, -1.0, 0.12),      420, 0.50),
         ("head-on",             (0.92, -0.38, 0.25),    420, 0.50),
         ("from above",          (-0.55, 0.65, 0.75),    420, 0.50),
         ("back",                (0.2, 1.0, 0.15),       420, 0.50),
         ("head",                (0.35, -1.0, 0.35),     150, 0.90)]


def shoot(mesh, shots, out, size=560, height=820):
    """shots: (title, focal point, camera direction, distance)."""
    pl = pv.Plotter(off_screen=True, shape=(1, len(shots)), window_size=(size * len(shots), height))
    for i, (title, c, d, dist) in enumerate(shots):
        pl.subplot(0, i)
        pl.set_background("white")
        pl.add_mesh(mesh, color="#8fdc86", smooth_shading=True, split_sharp_edges=True,
                    feature_angle=60, specular=0.25, specular_power=18, ambient=0.18)
        pl.add_text(title, font_size=10, color="#333333")
        d = np.array(d, float) / np.linalg.norm(d)
        pl.camera.position = (np.asarray(c) + d * dist).tolist()
        pl.camera.focal_point = list(c)
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30
    pl.screenshot(out)
    print("wrote", out)


mesh = pv.read(STL)
b = mesh.bounds
shots = []
for title, d, dist, focus in VIEWS:
    c = np.array([(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, b[4] + focus * (b[5] - b[4])])
    if title == "head":
        c = mesh.points[mesh.points[:, 2] > b[5] - 14].mean(axis=0) - [0, 0, 3]
    shots.append((title, c, d, dist))
shoot(mesh, shots, OUT)

hd = oe.Piece().head
eye = hd.O + oe.HEAD_SCALE * hd.eye_at @ hd.M
look = hd.eye_axes[2] @ hd.M                               # straight into the left eye
mid = hd.O + 0.45 * 35 * hd.F
shoot(mesh, [("vortex pupil", eye, look, 12),
             ("eye", eye, look + 0.25 * hd.F, 26),
             ("head", mid, np.cross(hd.U, hd.F) + 0.6 * hd.F + 0.35 * hd.U, 90),
             ("head-on, with the key", mid + 4 * hd.M[1], hd.F + 0.15 * hd.U, 120)], EYE_OUT, height=560)
