#!/usr/bin/env python3
"""Each sculpture on its own plinth, fused into a single print: the Orphic Egg, and the empty
coil (the same serpent after the hatching).

The plinth is a moulded round drum. From the bed up: a chamfered foot, a plain band, a cove,
a bead and a top band with a rounded edge. Every face is printable without support: the cove
only steps inward and the bead's underside is a 45-degree bevel. Sized from the shared grid:
the top face has radius 4.5R, 15 mm clear of the serpent's widest reach, and the drum is 1.5R
tall.

Each sculpture is unioned onto its drum with manifold3d at full detail.

Run from this directory, after orphic_egg.py and orphic_egg.py --serpent:
    python plinth.py      orphic_egg_plinth.stl, orphic_serpent_plinth.stl, orphic_plinths_views.png
"""

import math
import os
import sys

import manifold3d as m3
import numpy as np
import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import sacred_grid as sg

R = sg.R
TOP_R = 4.5 * R             # radius of the plinth's top face
HEIGHT = 1.5 * R
SINK = 0.5                  # sculpture set this far into the top, so no faces are coplanar
SEGMENTS = 360
PIECES = ("orphic_egg", "orphic_serpent")


def arc(cd, cz, r, a0, a1, n=12):
    """Profile points on a circle centred at inset cd, height cz, angles in degrees measured
    from the outward direction."""
    a = np.radians(np.linspace(a0, a1, n))
    return np.c_[cd - r * np.cos(a), cz + r * np.sin(a)]


def profile():
    """(inset from the base radius, height) from the bed to the top face."""
    p = [(0.6, 0.0), (0.0, 0.6), (0.0, 7.0), (0.8, 7.8)]          # chamfered foot, plain band
    p += arc(0.8, 11.8, 4.0, -90, -180, 14)[1:].tolist()          # cove: steps in, never out
    p += [(4.8, 12.2)]
    p += arc(5.0, 15.4, 2.4, -45, 90, 16).tolist()                # bead, 45-degree underside
    p += [(5.6, 17.8), (5.6, 21.3)]                               # top band
    p += arc(6.8, 21.3, 1.2, 0, 90, 8)[1:].tolist()               # rounded top edge
    p = np.array(p)
    assert abs(p[-1, 1] - HEIGHT) < 1e-9
    return p


def plinth():
    """The drum, revolved from its profile: rings of SEGMENTS points, capped top and bottom."""
    prof = profile()
    base = TOP_R + prof[-1, 0]
    a = np.linspace(0, 2 * math.pi, SEGMENTS, endpoint=False)
    rings = [np.c_[(base - d) * np.cos(a), (base - d) * np.sin(a), np.full(SEGMENTS, z)] for d, z in prof]
    n, k = SEGMENTS, len(rings)
    V = np.r_[np.concatenate(rings), [[0, 0, 0], [0, 0, HEIGHT]]]
    i = np.arange(n)
    j = (i + 1) % n
    F = []
    for r in range(k - 1):
        lo, hi = r * n, (r + 1) * n
        F += [np.c_[lo + i, lo + j, hi + j], np.c_[lo + i, hi + j, hi + i]]
    F += [np.c_[np.full(n, k * n), j, i], np.c_[np.full(n, k * n + 1), (k - 1) * n + i, (k - 1) * n + j]]
    return trimesh.Trimesh(V, np.concatenate(F), process=False)


def manifold(m):
    return m3.Manifold(m3.Mesh(vert_properties=np.asarray(m.vertices, np.float32),
                               tri_verts=np.asarray(m.faces, np.uint32)))


def main():
    base = plinth()
    assert base.is_watertight and base.is_volume, "plinth mesh is not a closed solid"
    print(f"    plinth {base.extents[0]:.1f} mm across, {HEIGHT:.1f} mm tall, {base.volume / 1000:.0f} cm3")
    ok, meshes = True, []
    for name in PIECES:
        piece = manifold(trimesh.load(f"{name}.stl")).translate((0.0, 0.0, HEIGHT - SINK))
        m = sg.to_trimesh(manifold(base) + piece)
        m.merge_vertices(digits_vertex=5)
        m.update_faces(m.nondegenerate_faces())
        m.update_faces(m.unique_faces())
        m.remove_unreferenced_vertices()
        shells = sorted(m.split(only_watertight=False), key=lambda p: -len(p.faces))
        m = shells[0]                  # zero-volume scraps where the tail tip grazes the top
        out = f"{name}_plinth.stl"
        m.export(out)
        back = trimesh.load(out)
        ok &= sg.check(back, out)
        print(f"        {len(shells) - 1} scraps dropped, {os.path.getsize(out) / 1e6:.1f} MB")
        meshes.append(back)

    import pyvista as pv
    c = np.array([0.0, 0.0, 0.45 * meshes[0].extents[2]])
    shots = []
    for m, name in zip(meshes, ("egg", "empty coil")):
        pm = pv.wrap(m)
        shots += [(f"{name}: front", pm, c, (0.0, -1.0, 0.15)),
                  (f"{name}: three-quarter", pm, c, (0.8, -0.8, 0.4))]
    sheet(shots, "orphic_plinths_views.png")
    sys.exit(0 if ok else 1)


def sheet(shots, out, size=560, height=760):
    """shots: (title, mesh, focal point, camera direction) -- one mesh per tile."""
    import pyvista as pv
    pl = pv.Plotter(off_screen=True, shape=(1, len(shots)), window_size=(size * len(shots), height))
    for i, (title, mesh, c, d) in enumerate(shots):
        pl.subplot(0, i)
        pl.set_background("white")
        pl.add_mesh(mesh, color="#8fdc86", smooth_shading=True, split_sharp_edges=True,
                    feature_angle=60, specular=0.25, specular_power=18, ambient=0.18)
        pl.add_text(title, font_size=10, color="#333333")
        d = np.array(d, float) / np.linalg.norm(d)
        pl.camera.position = (c + d * 470).tolist()
        pl.camera.focal_point = c.tolist()
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30
    pl.screenshot(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
