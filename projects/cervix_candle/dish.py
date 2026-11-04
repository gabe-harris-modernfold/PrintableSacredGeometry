#!/usr/bin/env python3
"""The dish the uterus candle holder stands in: a vesica piscis, 207.8 mm point to point and
120 mm across, with a shallow wax well, a rolled lip, and a raised boss in the centre. The
boss is a smaller vesica on the same axis and carries a socket shaped like the sculpture's base
cross-section, so the sculpture drops in and can only face one way. The points lie left and
right, along the sculpture's width.

The vesica is two circles of radius R_VESICA, each through the other's centre. Every contour of
the dish (the rim, the lip, the well floor, the boss) is the lens between two circles on those
same two centres, so each one is an exact parallel offset of the outline and the points stay
sharp all the way down.

Every face prints without support: the outer wall leans out at about 30 degrees from vertical,
and the boss and inner wall slope at 45 degrees or gentler.

Run from this directory, after cervix_candle.py:   python dish.py      dish.stl
"""

import manifold3d as m3
import numpy as np
import trimesh
from shapely.affinity import translate
from shapely.geometry.polygon import orient

R_VESICA = 120.0        # radius of the two circles; the lens is R wide and R*sqrt(3) long
FLOOR = 3.0             # wax-well floor thickness
RIM_H = 18.0
BOSS_H = 10.0           # top of the centre boss
BOSS_MARGIN = 8.0       # boss top beyond the socket, measured square to the boss edge
SOCKET_DEPTH = 5.0
SOCKET_FLOOR_Z = BOSS_H - SOCKET_DEPTH
CLEAR = 0.4             # socket clearance around the sculpture's base
SEGMENTS = 1440         # rays round the axis; a multiple of 4 so both points are sampled

# The wall, as (inset from the outline, z) from the foot to the edge of the well floor
WALL = [(8, 0), (0, 14), (0, 16), (1, 17.4), (2.5, RIM_H), (4, RIM_H), (5.5, 16.5),
        (18.5, FLOOR)]


def base_polygon(path="cervix_candle_upright.stl", at=1.0):
    m = trimesh.load(path)
    sec = m.section(plane_origin=(0, 0, at), plane_normal=(0, 0, 1))
    planar, _ = sec.to_2D(to_2D=np.eye(4))          # identity: planar x, y are world x, y
    return max(planar.polygons_full, key=lambda p: p.area)


def centre_offset():
    """Shift that puts the sculpture's footprint centroid on the dish's axis."""
    c = base_polygon().centroid
    return -c.x, -c.y


def lens_ring(c, rho, z, theta):
    """Points of the lens between circles of radius rho centred at (0, +-c), on rays from the
    origin at angles theta, at height z."""
    u = np.stack([np.cos(theta), np.sin(theta)], axis=1)
    t = np.inf
    for s in (1, -1):
        uc = u[:, 1] * s * c
        t = np.minimum(t, uc + np.sqrt(uc ** 2 - c ** 2 + rho ** 2))
    return np.column_stack([u * t[:, None], np.full(len(theta), z)])


def boss_radius(socket):
    """Smallest boss vesica whose outline clears the socket by BOSS_MARGIN all round."""
    p = np.asarray(socket.exterior.coords)
    lo, hi = 1.0, R_VESICA
    for _ in range(60):
        rb = 0.5 * (lo + hi)
        fits = all(np.hypot(p[:, 0], p[:, 1] - s * rb / 2).max() <= rb - BOSS_MARGIN
                   for s in (1, -1))
        lo, hi = (lo, rb) if fits else (rb, hi)
    return hi


def lens_dish(boss_r):
    theta = np.linspace(0, 2 * np.pi, SEGMENTS, endpoint=False)
    c = R_VESICA / 2
    rings = [lens_ring(c, R_VESICA - d, z, theta) for d, z in WALL]
    # the boss: a smaller vesica, its flank a 45-degree offset from the well floor to its top
    rings += [lens_ring(boss_r / 2, boss_r + BOSS_H - z, z, theta) for z in (FLOOR, BOSS_H)]
    n = SEGMENTS
    verts = np.vstack([[0, 0, 0]] + rings + [[[0, 0, BOSS_H]]])
    top = len(verts) - 1
    k = np.arange(n)
    k1 = (k + 1) % n
    faces = [np.column_stack([np.zeros(n, int), 1 + k1, 1 + k])]               # underside
    for i in range(len(rings) - 1):
        a, b = 1 + i * n, 1 + (i + 1) * n
        faces += [np.column_stack([a + k, a + k1, b + k1]),
                  np.column_stack([a + k, b + k1, b + k])]
    last = 1 + (len(rings) - 1) * n
    faces.append(np.column_stack([np.full(n, top), last + k, last + k1]))       # boss top
    mesh = trimesh.Trimesh(verts, np.vstack(faces), process=True)
    if mesh.volume < 0:
        mesh.invert()
    return mesh


def build():
    dx, dy = centre_offset()
    base = translate(base_polygon(), dx, dy)
    socket = base.buffer(CLEAR, join_style="round")
    boss_r = boss_radius(socket)
    blank = lens_dish(boss_r)
    dish = m3.Manifold(m3.Mesh(blank.vertices.astype(np.float32),
                               blank.faces.astype(np.uint32)))
    ring = np.asarray(orient(socket, 1.0).exterior.coords)[:-1]      # counter-clockwise
    cutter = m3.Manifold.extrude(m3.CrossSection([ring]), SOCKET_DEPTH + 1.0) \
        .translate((0, 0, SOCKET_FLOOR_Z))
    out = (dish - cutter).to_mesh()
    mesh = trimesh.Trimesh(out.vert_properties[:, :3], out.tri_verts)
    mesh.export("dish.stl")
    return mesh, boss_r, base


if __name__ == "__main__":
    mesh, boss_r, base = build()
    c = R_VESICA / 2
    floor_half = (R_VESICA - WALL[-1][0]) - c, np.sqrt((R_VESICA - WALL[-1][0]) ** 2 - c ** 2)
    boss_half = boss_r / 2 + BOSS_H - FLOOR, np.sqrt((boss_r + BOSS_H - FLOOR) ** 2 - (boss_r / 2) ** 2)
    print("base footprint", np.round(base.bounds, 1), "boss vesica r", round(boss_r, 1),
          "boss top", round(np.sqrt(boss_r ** 2 - (boss_r / 2) ** 2) * 2, 1), "x",
          round(boss_r, 1))
    print("well floor width: across", round(floor_half[0] - boss_half[0], 1),
          "point-ward", round(floor_half[1] - boss_half[1], 1))
    print("dish extents", np.round(mesh.extents, 1), "watertight", mesh.is_watertight,
          "volume cm3", round(mesh.volume / 1000, 1))
