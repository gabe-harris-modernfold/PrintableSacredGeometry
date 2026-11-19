#!/usr/bin/env python3
"""Ratcheted dodecahedra: Dan Winter's "DNA helix is ratcheted dodecahedra" as a printable lattice.

Stand a dodecahedron on a face. Its twenty corners lie in four rings of five: a small face ring
top and bottom, and a wide zig-zag girdle round the waist. Choose the edge so the girdle radius is
phi^3. Then everything in Winter's graphic falls out exactly:

    girdle radius            phi^3          ("radius of turn")
    girdle side, seen top    phi^2          ("horizontal increment of turn"; 2 phi^3 sin 18 = phi^2)
    girdle zig-zag rise      phi            ("vertical increment of turn")
    face-ring radius         phi^2
    face-ring height         phi^2 + phi/2  (so a face edge climbs phi^2 while stepping out phi)
    edge                     2 phi^2 sin 36

Now "ratchet": copy the dodecahedron, turn it 36 degrees about its axis and lift it phi. The copy's
lower girdle lands exactly on the original's upper girdle -- the two share five corners -- so the
stack locks together like a ratchet, one tooth (36 deg) per step, ten steps to the full turn. Each
girdle edge that climbs +36 deg joins up with the next copy's, and the chain of them is a helix:
radius phi^3, 36 deg and phi per step, pitch 10 phi. Five such right-handed strands thread the
column (72 deg apart) and five left-handed ones mirror them, because a 36-deg turn of a dodecahedron
is the same as flipping it end-for-end. A double helix is two of the right-handed strands; strands
144 deg apart leave a narrow and a wide groove (144 / 216 deg), as DNA's minor and major grooves do.

The print is one full turn -- ten dodecahedra -- as a strut lattice, with the two helix strands in
heavier struts. Seen down the axis the girdle decagon touches the 4R = 60 mm circle of the matched
sacred-geometry set (lib/sacred_grid.py).

Printability, by design rather than by slicer supports:
  * struts have a diamond section with a vertical diagonal, so even the shallow girdle edges
    (31.7 deg) present 53-deg faces underneath, and level face edges are 45-deg ridges that bridge
    corner to corner (43.6 mm);
  * every face-ring corner of every dodecahedron sits on one of ten vertical lines (radius phi^2,
    every 36 deg), two dodecahedra apart. Posts along those lines carry each otherwise floating
    face pentagon, so no corner starts in mid-air. They read as the column's inner decagonal cage;
    set POSTS = False for the bare dodecahedra (then print with tree supports).
  * no base: it stands on the bottom face pentagon, its ridges trimmed flat by BED_FLAT, plus
    five post feet. Use a slicer brim for adhesion; it peels off.

Concept and original graphic: Dan Winter, goldenmean.info/dnamanifesto ("DNA Helix is Wratcheted
Dodecahedra"). This script builds his picture as an exact, printable solid.

Run from this directory:  python ratcheted_dodecahedra.py
Writes ratcheted_dodecahedra.stl and ratcheted_dodecahedra_views.png.
"""

import math
import os
import sys
import time
from itertools import combinations

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

PHI = (1 + 5 ** 0.5) / 2

# ------------------------------------------------------------------ parameters
TURN = 10                       # dodecahedra in the print: one full turn of the helix
STEP_DEG = 36.0                 # the ratchet tooth
GRID_R = 15.0                   # the matched set's one measure (lib/sacred_grid.R)
S = 4 * GRID_R / PHI ** 3       # mm per unit: girdle radius phi^3 -> 4R = 60 mm
STRANDS = (0.0, 144.0)          # phase of the two highlighted strands (multiples of 72)
STRUT = 6.0                     # lattice strut, diamond diagonal (mm)
STRAND = 9.0                    # helix strand strut, diamond diagonal (mm)
POST = 4.5                      # post diagonal (mm)
POSTS = True
BED_FLAT = 1.2                  # trimmed off the bottom struts' ridges to give them a footprint
BED = 320.0
OUT = "ratcheted_dodecahedra"


# ------------------------------------------------------------------ exact geometry (units)

def polar(r, deg, z):
    a = math.radians(deg)
    return np.array([r * math.cos(a), r * math.sin(a), z])


R_GIRDLE = PHI ** 3
R_FACE = PHI ** 2
Z_GIRDLE = PHI / 2
Z_FACE = math.sqrt(PHI ** 6 + Z_GIRDLE ** 2 - PHI ** 4)       # same circumsphere as the girdle
EDGE = 2 * PHI ** 2 * math.sin(math.radians(36))


def dodecahedron(k):
    """Corners and edges of the k-th ratcheted copy: turned 36k deg, lifted phi*k."""
    t, z = STEP_DEG * k, PHI * k
    top = [polar(R_FACE, t + 72 * m, z + Z_FACE) for m in range(5)]
    up = [polar(R_GIRDLE, t + 72 * m, z + Z_GIRDLE) for m in range(5)]
    lo = [polar(R_GIRDLE, t + 36 + 72 * m, z - Z_GIRDLE) for m in range(5)]
    bot = [polar(R_FACE, t + 36 + 72 * m, z - Z_FACE) for m in range(5)]
    V = top + up + lo + bot                       # indices 0-4, 5-9, 10-14, 15-19
    E = []
    for m in range(5):
        n = (m + 1) % 5
        E += [(m, n), (m, 5 + m),                 # top face, top -> upper girdle
              (5 + m, 10 + m), (10 + m, 5 + n),   # girdle zig-zag
              (10 + m, 15 + m), (15 + m, 15 + n)] # lower girdle -> bottom, bottom face
    return V, E


def key(p):
    return tuple(np.round(p, 6))


def build_lattice(n=TURN):
    """Merged corners and edges of n ratcheted copies, plus the two helix strands' edges."""
    idx, V, E = {}, [], set()

    def vid(p):
        k = key(p)
        if k not in idx:
            idx[k] = len(V)
            V.append(np.asarray(p, float))
        return idx[k]

    for k in range(n):
        Vk, Ek = dodecahedron(k)
        ids = [vid(p) for p in Vk]
        E |= {tuple(sorted((ids[i], ids[j]))) for i, j in Ek}

    strands = []
    for c in STRANDS:
        pts = [polar(R_GIRDLE, c + STEP_DEG * j, PHI * j + Z_GIRDLE) for j in range(-1, n)]
        strands.append([tuple(sorted((idx[key(a)], idx[key(b)]))) for a, b in zip(pts, pts[1:])])
    return np.array(V), sorted(E), strands


def posts(V):
    """Vertical struts joining face-ring corners on the ten radius-phi^2 lines, plus feet."""
    out, feet = [], []
    r = np.hypot(V[:, 0], V[:, 1])
    on = np.where(np.abs(r - R_FACE) < 1e-6)[0]
    ang = np.round(np.degrees(np.arctan2(V[on, 1], V[on, 0])), 3) % 360   # 360 -> 0
    zmin = V[:, 2].min()
    for a in np.unique(ang):
        col = on[ang == a]
        col = col[np.argsort(V[col, 2])]
        out += [(int(i), int(j)) for i, j in zip(col, col[1:])]
        if V[col[0], 2] > zmin + 1e-6:
            feet.append(int(col[0]))
    return out, feet


def verify(V, E, strands):
    """Every number in the graphic, asserted."""
    for k in range(3):
        Vk, Ek = dodecahedron(k)
        L = [np.linalg.norm(Vk[i] - Vk[j]) for i, j in Ek]
        assert len(Ek) == 30 and np.allclose(L, EDGE), "not a regular dodecahedron"
        # 12 pentagonal faces <=> every corner has degree 3 and the solid is convex-equal-edge
        deg = np.bincount(np.array(Ek).ravel(), minlength=20)
        assert (deg == 3).all()
        R = [np.linalg.norm(p - np.array([0, 0, PHI * k])) for p in Vk]
        assert np.allclose(R, R[0]), "corners not on one sphere"
    assert math.isclose(Z_FACE, PHI ** 2 + PHI / 2), "face ring height"
    assert math.isclose(2 * R_GIRDLE * math.sin(math.radians(18)), PHI ** 2), "horizontal increment"
    # ratchet: copy k's upper girdle == copy k+1's lower girdle
    A, _ = dodecahedron(0)
    B, _ = dodecahedron(1)
    assert {key(p) for p in A[5:10]} == {key(p) for p in B[10:15]}, "girdles do not lock"
    assert len(V) == 20 * TURN - 5 * (TURN - 1)
    Es = set(E)
    for s in strands:
        assert set(s) <= Es, "strand leaves the lattice"
        a, b = V[s[0][0]], V[s[-1][1]]
        p0 = V[[i for i in s[0] if i not in s[1]][0]]
        rise = (V[list(s[-1])][:, 2].max() - p0[2]) / len(s)
        assert math.isclose(rise, PHI), "vertical increment"
    print(f"verified: edge {EDGE:.4f} = 2 phi^2 sin36, girdle r phi^3 = {R_GIRDLE:.4f}, "
          f"top-view side phi^2 = {PHI**2:.4f}, rise phi = {PHI:.4f}, face ring z = phi^2+phi/2")
    print(f"lattice: {len(V)} corners, {len(E)} edges in {TURN} dodecahedra; "
          f"strands of {len(strands[0])} edges = {len(strands[0]) * STEP_DEG:.0f} deg")


def closest(p, q, r, s):
    d1, d2, w = q - p, s - r, p - r
    a, e, f = d1 @ d1, d2 @ d2, d2 @ w
    b, c = d1 @ d2, d1 @ w
    den = a * e - b * b
    t = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-12 else 0.0
    u = np.clip((b * t + f) / e, 0, 1)
    t = np.clip((b * u - c) / a, 0, 1)
    return np.linalg.norm((p + d1 * t) - (r + d2 * u))


# ------------------------------------------------------------------ solids (mm)

def section(p, q, h):
    """Diamond cross-section corners (vertical diagonal) for the strut p->q."""
    d = (q - p) / np.linalg.norm(q - p)
    u = np.cross(d, [0, 0, 1.0])
    if np.linalg.norm(u) < 1e-9:                   # vertical: face the diamond tangentially
        u = np.cross([0, 0, 1.0], [p[0], p[1], 0])
        u = u / (np.linalg.norm(u) or 1) if np.linalg.norm(u) > 1e-9 else np.array([1.0, 0, 0])
    u /= np.linalg.norm(u)
    w = np.cross(u, d)
    return [h * w, h * u, -h * w, -h * u]


def hull_solid(pts):
    from build123d import Face, Shell, Solid, Wire
    pts = np.asarray(pts, float)
    hull = ConvexHull(pts)
    faces, seen = [], {}
    for simplex, eq in zip(hull.simplices, hull.equations):
        k = tuple(np.round(eq, 5))
        seen.setdefault(k, set()).update(simplex.tolist())
    for eq, ids in seen.items():
        P = pts[sorted(ids)]
        c = P.mean(axis=0)
        n = np.array(eq[:3])
        a = P[0] - c
        a /= np.linalg.norm(a)
        b = np.cross(n, a)
        order = np.argsort(np.arctan2((P - c) @ b, (P - c) @ a))
        faces.append(Face(Wire.make_polygon([tuple(x) for x in P[order]], close=True)))
    return Solid(Shell(faces))


def build_solid(V, E, strands, post_edges, feet):
    Vm = V * S
    Vm[:, 2] -= Vm[:, 2].min()                     # bottom face centre-lines on z = 0
    strand_set = {e for s in strands for e in s}
    h_of = {e: (STRAND if e in strand_set else STRUT) / 2 for e in E}
    for e in post_edges:
        h_of[e] = POST / 2
    ends = {}                                      # node -> list of section corners
    solids = []

    def strut(i, j, h, p=None, q=None):
        p = Vm[i] if p is None else p
        q = Vm[j] if q is None else q
        sec = section(p, q, h)
        pts = [p + s for s in sec] + [q + s for s in sec]
        solids.append(hull_solid(pts))
        if i is not None:
            ends.setdefault(i, []).extend(p + s for s in sec)
        if j is not None:
            ends.setdefault(j, []).extend(q + s for s in sec)

    for e, h in h_of.items():
        strut(e[0], e[1], h)
    z_cut = -STRUT / 2 + BED_FLAT                  # bed plane, through the bottom face's struts
    for i in feet:                                 # posts down to the bed (trimmed with the rest)
        foot = Vm[i].copy()
        foot[2] = -STRUT
        strut(None, i, POST / 2, p=foot)
    for i, corners in ends.items():                # joint = hull of all incident end sections
        solids.append(hull_solid(np.vstack(corners + [Vm[i][None, :]])))

    from build123d import Box, Pos
    t0 = time.time()
    body = solids[0].fuse(*solids[1:]).clean()
    body = (body - Pos(0, 0, z_cut - 500) * Box(1000, 1000, 1000)).clean()
    body = Pos(0, 0, -z_cut) * body                # bed plane -> z = 0
    Vm[:, 2] -= z_cut
    print(f"fused {len(solids)} solids in {time.time() - t0:.1f} s")
    return body, Vm


def strand_mesh(strands, V):
    """Just the two strands, for colouring the preview."""
    import trimesh
    Vm = V * S
    Vm[:, 2] += STRUT / 2 - BED_FLAT - Vm[:, 2].min()      # same placement as the print
    meshes = []
    for s in strands:
        tri = []
        for i, j in s:
            sec = section(Vm[i], Vm[j], STRAND / 2 + 0.15)
            pts = [Vm[i] + c for c in sec] + [Vm[j] + c for c in sec]
            tri.append(trimesh.convex.convex_hull(np.array(pts)))
        meshes.append(trimesh.util.concatenate(tri))
    return meshes


def main():
    V, E, strands = build_lattice()
    verify(V, E, strands)
    post_edges, feet = posts(V) if POSTS else ([], [])
    if POSTS:
        print(f"posts: {len(post_edges)} on 10 lines, {len(feet)} feet to the bed")

    # strut clearance: centre-line distance between struts that share no corner
    allE = [tuple(e) for e in E] + post_edges
    gaps = sorted(closest(V[a] * S, V[b] * S, V[c] * S, V[d] * S)
                  for (a, b), (c, d) in combinations(allE, 2) if not {a, b} & {c, d})
    print(f"closest non-adjacent struts: {gaps[0]:.2f} mm centre to centre "
          f"({sum(g < STRUT for g in gaps)} pairs under one strut width -> fused crossings)")

    from build123d import export_stl
    body, Vm = build_solid(V, E, strands, post_edges, feet)
    export_stl(body, OUT + ".stl", tolerance=0.02, angular_tolerance=0.1)

    import trimesh
    m = trimesh.load(OUT + ".stl")
    ext = m.extents
    print(f"{OUT}.stl: {len(m.faces)} tris, extents {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm, "
          f"watertight {m.is_watertight}, volume {m.volume / 1000:.1f} cm^3, "
          f"bodies {len(m.split(only_watertight=False))}")
    assert m.is_watertight and (ext <= BED).all()

    from render_mesh import MAT, render
    MAT["strandA"] = {"col": (214, 92, 70), "amb": 0.42, "spec": 0.25, "shin": 22.0}
    MAT["strandB"] = {"col": (64, 128, 214), "amb": 0.42, "spec": 0.25, "shin": 22.0}
    sa, sb = strand_mesh(strands, V)
    layers = [(m.vertices, m.faces, "frame"), (sa.vertices, sa.faces, "strandA"),
              (sb.vertices, sb.faces, "strandB")]
    from PIL import Image
    tiles = []
    for name, az, el in (("side", 0, 0), ("oblique", 30, 22), ("top", 0, 90)):
        f = f"_view_{name}.png"
        render(layers, az, el, w=700, h=1100 if el < 80 else 700, out=f)
        tiles.append(Image.open(f))
    W = sum(t.width for t in tiles)
    sheet = Image.new("RGB", (W, 1100), (242, 242, 240))
    x = 0
    for t in tiles:
        sheet.paste(t, (x, (1100 - t.height) // 2))
        x += t.width
    sheet.save(OUT + "_views.png")
    for t in tiles:
        os.remove(t.filename)
    print("wrote", OUT + "_views.png")


if __name__ == "__main__":
    main()
