#!/usr/bin/env python3
"""The dish the uterus candle holder stands in: 155 mm across, with a shallow wax well, a
rolled lip, and a raised round boss in the centre. The
boss carries a socket shaped like the sculpture's base cross-section, so the sculpture drops
in and can only face one way.

Every face prints without support: the outer wall leans out at about 30 degrees from vertical,
and the boss and inner wall slope at 45 degrees or gentler.

Run from this directory, after cervix_candle.py:   python dish.py      dish.stl
"""

import manifold3d as m3
import numpy as np
import trimesh
from shapely.affinity import translate
from shapely.geometry.polygon import orient

DIAMETER = 155.0
FLOOR = 3.0             # wax-well floor thickness
RIM_H = 18.0
BOSS_H = 10.0           # top of the centre boss
SOCKET_DEPTH = 5.0
SOCKET_FLOOR_Z = BOSS_H - SOCKET_DEPTH
CLEAR = 0.4             # socket clearance around the sculpture's base
SEGMENTS = 360


def base_polygon(path="cervix_candle_upright.stl", at=1.0):
    m = trimesh.load(path)
    sec = m.section(plane_origin=(0, 0, at), plane_normal=(0, 0, 1))
    planar, _ = sec.to_2D(to_2D=np.eye(4))          # identity: planar x, y are world x, y
    return max(planar.polygons_full, key=lambda p: p.area)


def centre_offset():
    """Shift that puts the sculpture's footprint centroid on the dish's axis."""
    c = base_polygon().centroid
    return -c.x, -c.y


def profile(boss_r):
    R = DIAMETER / 2
    return [(0, 0), (R - 8, 0), (R, 14), (R, 16), (R - 1, 17.4), (R - 2.5, RIM_H),
            (R - 4, RIM_H), (R - 5.5, 16.5), (R - 18.5, FLOOR), (boss_r + BOSS_H - FLOOR, FLOOR),
            (boss_r, BOSS_H), (0, BOSS_H)]


def build():
    dx, dy = centre_offset()
    base = translate(base_polygon(), dx, dy)
    socket = base.buffer(CLEAR, join_style="round")
    boss_r = max(np.hypot(*np.asarray(socket.exterior.coords).T)) + 8.0
    dish = m3.Manifold.revolve(m3.CrossSection([profile(boss_r)]), SEGMENTS)
    ring = np.asarray(orient(socket, 1.0).exterior.coords)[:-1]      # counter-clockwise
    cutter = m3.Manifold.extrude(m3.CrossSection([ring]), SOCKET_DEPTH + 1.0) \
        .translate((0, 0, SOCKET_FLOOR_Z))
    out = (dish - cutter).to_mesh()
    mesh = trimesh.Trimesh(out.vert_properties[:, :3], out.tri_verts)
    mesh.export("dish.stl")
    return mesh, boss_r, base


if __name__ == "__main__":
    mesh, boss_r, base = build()
    print("base footprint", np.round(base.bounds, 1), "boss r", round(boss_r, 1))
    print("dish extents", np.round(mesh.extents, 1), "watertight", mesh.is_watertight,
          "volume cm3", round(mesh.volume / 1000, 1))
