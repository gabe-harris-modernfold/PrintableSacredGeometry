#!/usr/bin/env python3
"""The shared grid behind the matched sacred-geometry set ("From One Circle").

Every piece is sized from one base circle of radius R, so the prints share one scale:

  plates   circles of radius R on the triangular lattice of spacing R, one lattice
           direction vertical. Every figure sits on the same lattice, so the plates stack
           in register on three pins through a common ring frame.
  solids   edge frames, scaled so that seen down a three-fold axis the outermost corners
           land on the circle through Metatron's Cube's outer six centres (radius 4R). The
           tetrahedron, cube, octahedron and star tetrahedron then lie exactly on lines of
           the plate; the dodecahedron and icosahedron carry phi and only touch the circle.
           Cube and octahedron at this size are exact duals (edges cross at midpoints).
  stellations  grown from those same solids by the method of
           projects/stellated_platonic_solids/_legacy/stellated_platonics.py: each face's
           pyramid apex is where its neighbours' face planes cross. The tetrahedron and cube
           do not stellate and get their dual compound instead (stella octangula; cube +
           octahedron). Collinear struts are merged, so stars read as long straight lines.
  spheres  spheres of radius R on a cubic grid of spacing R*sqrt(3/2), whose shadow down
           the body diagonal is the plate lattice: 2x2x2 casts the Seed of Life, 3x3x3 the
           Flower of Life.

Geometry is built with manifold3d (exact unions, always watertight) and shapely (plate
outlines). All lengths in mm."""

import math
from itertools import combinations

import manifold3d as m3
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from shapely.geometry import LineString, Point
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

PHI = (1 + 5 ** 0.5) / 2
R = 15.0                    # the one measure
BED = 320.0                 # printer bed, cube

# plates
PLATE_T = 3.0
STROKE = 1.8                # circle line width
LINE = 1.6                  # straight lines of Metatron's Cube
PATH = 2.2                  # Tree of Life paths
NODE_R = 0.27 * R           # Tree of Life sephirah disc
SPOKE = 3.0
FRAME_IN = 5 * R            # the Fruit of Life's outer circles reach exactly 5R
FRAME_OUT = 5 * R + 10.0
PIN_AT = (FRAME_IN + FRAME_OUT) / 2
PIN_HOLE = 3.4
PIN_D = 3.0
PIN_ANGLES = (0.0, 120.0, 180.0)   # no rotation or flip maps this set onto itself
QUAD = 32                   # shapely segments per quarter circle
RAYS = (30.0, 90.0, 150.0, 210.0, 270.0, 330.0)   # the six lattice directions

# frames and spheres
STRUT_D = 5.0
NODE_D = 7.0
SEGS = 24
BED_FLAT = 0.6              # trimmed off the bottom so round parts get a footprint

PLATES = ("circle", "vesica_piscis", "tripod_of_life", "seed_of_life", "flower_of_life",
          "fruit_of_life", "metatrons_cube", "tree_of_life")


# ---------------------------------------------------------------- the plate lattice

def ray(deg, r):
    a = math.radians(deg)
    return (r * math.cos(a), r * math.sin(a))


def lattice(i, j):
    """Lattice point i*A + j*B, A at 30 deg and B vertical, both of length R."""
    ax, ay = ray(30.0, R)
    return (i * ax, i * ay + j * R)


def hex_points(n):
    """Lattice points within n steps of the centre (1, 7, 19, 37 ...)."""
    return [lattice(i, j) for i in range(-n, n + 1) for j in range(-n, n + 1)
            if max(abs(i), abs(j), abs(i + j)) <= n]


SEPHIROT = {   # the Tree of Life on the Flower of Life lattice, Tiferet at the centre
    "Kether": (0.0, 2 * R), "Chokhmah": ray(60.0, math.sqrt(3) * R), "Binah": ray(120.0, math.sqrt(3) * R),
    "Daath": (0.0, R), "Chesed": ray(30.0, R), "Gevurah": ray(150.0, R), "Tiferet": (0.0, 0.0),
    "Netzach": ray(330.0, R), "Hod": ray(210.0, R), "Yesod": (0.0, -R), "Malkuth": (0.0, -2 * R),
}
PATHS = [("Kether", "Chokhmah"), ("Kether", "Binah"), ("Kether", "Tiferet"), ("Chokhmah", "Binah"),
         ("Chokhmah", "Chesed"), ("Chokhmah", "Tiferet"), ("Binah", "Gevurah"), ("Binah", "Tiferet"),
         ("Chesed", "Gevurah"), ("Chesed", "Tiferet"), ("Chesed", "Netzach"), ("Gevurah", "Tiferet"),
         ("Gevurah", "Hod"), ("Tiferet", "Netzach"), ("Tiferet", "Hod"), ("Tiferet", "Yesod"),
         ("Netzach", "Hod"), ("Netzach", "Yesod"), ("Netzach", "Malkuth"), ("Hod", "Yesod"),
         ("Hod", "Malkuth"), ("Yesod", "Malkuth")]


def figure(name):
    """A plate figure: circles (x, y, r), lines (p, q, width), discs (x, y, r)."""
    c, lines, discs = [], [], []
    seed = [(0.0, 0.0)] + [ray(a, R) for a in RAYS]
    fruit = [(0.0, 0.0)] + [ray(a, 2 * R) for a in RAYS] + [ray(a, 4 * R) for a in RAYS]
    if name == "circle":
        c = [(0.0, 0.0, R)]
    elif name == "vesica_piscis":          # the second circle drawn on the first one's rim
        c = [(0.0, 0.0, R), (0.0, R, R)]
    elif name == "tripod_of_life":
        c = [(0.0, 0.0, R), (0.0, R, R), (*ray(30.0, R), R)]
    elif name == "seed_of_life":
        c = [(x, y, R) for x, y in seed]
    elif name == "flower_of_life":
        c = [(x, y, R) for x, y in hex_points(2)] + [(0.0, 0.0, 3 * R)]
    elif name in ("fruit_of_life", "metatrons_cube"):
        c = [(x, y, R) for x, y in fruit]
        if name == "metatrons_cube":
            lines = [(p, q, LINE) for p, q in combinations(fruit, 2)]
    elif name == "tree_of_life":
        s = SEPHIROT
        lines = [(s[a], s[b], PATH) for a, b in PATHS]
        discs = [(x, y, NODE_R) for k, (x, y) in s.items() if k != "Daath"]
        c = [(*s["Daath"], NODE_R - PATH / 2)]      # the hidden sephirah, drawn as a ring
    else:
        raise ValueError(name)
    return {"circles": c, "lines": lines, "discs": discs}


def _ring(x, y, r, w):
    return Point(x, y).buffer(r + w / 2, quad_segs=QUAD).difference(Point(x, y).buffer(r - w / 2, quad_segs=QUAD))


def _art(name):
    f = figure(name)
    w = PATH if name == "tree_of_life" else STROKE
    parts = [_ring(x, y, r, w) for x, y, r in f["circles"]]
    parts += [LineString([p, q]).buffer(lw / 2, quad_segs=8) for p, q, lw in f["lines"]]
    parts += [Point(x, y).buffer(r, quad_segs=QUAD) for x, y, r in f["discs"]]
    return unary_union(parts)


def spokes(name):
    """Segments tying the figure to the frame along the lattice rays it does not reach."""
    art, out = _art(name), []
    for a in RAYS:
        hit = LineString([(0, 0), ray(a, FRAME_IN)]).intersection(art)
        far = max((math.hypot(*xy) for g in getattr(hit, "geoms", [hit]) for xy in g.coords), default=0.0)
        if far < FRAME_IN - 0.5:
            out.append((ray(a, far - SPOKE), ray(a, FRAME_IN + 1.0)))
    return out


def plate_outline(name):
    """Figure + spokes + common frame, minus the pin holes, as one shapely polygon."""
    frame = Point(0, 0).buffer(FRAME_OUT, quad_segs=QUAD * 2).difference(
        Point(0, 0).buffer(FRAME_IN, quad_segs=QUAD * 2))
    bars = [LineString(s).buffer(SPOKE / 2, cap_style="flat") for s in spokes(name)]
    shape = unary_union([_art(name), frame] + bars)
    pins = unary_union([Point(*ray(a, PIN_AT)).buffer(PIN_HOLE / 2, quad_segs=QUAD) for a in PIN_ANGLES])
    return shape.difference(pins)


def _cross_section(geom):
    contours = []
    for poly in getattr(geom, "geoms", [geom]):
        poly = orient(poly, 1.0)
        contours.append(np.asarray(poly.exterior.coords)[:-1])
        contours += [np.asarray(i.coords)[:-1] for i in poly.interiors]
    return m3.CrossSection(contours, m3.FillRule.Positive)


def plate(name):
    return m3.Manifold.extrude(_cross_section(plate_outline(name)), PLATE_T)


def pin():
    """A D-shaped pin that prints lying on its flat: through all eight plates."""
    L = len(PLATES) * PLATE_T
    rod = m3.Manifold.cylinder(L, PIN_D / 2, PIN_D / 2, 48).rotate((0.0, 90.0, 0.0))
    return rod.trim_by_plane((0.0, 0.0, 1.0), -PIN_D / 2 + 0.4).translate((0.0, 0.0, PIN_D / 2 - 0.4))


# ---------------------------------------------------------------- solids on the grid

def _rot_between(a, b):
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    v, c = np.cross(a, b), float(a @ b)
    if np.linalg.norm(v) < 1e-12:
        if c > 0:
            return np.eye(3)
        p = np.cross(a, [1.0, 0.0, 0.0] if abs(a[0]) < 0.9 else [0.0, 1.0, 0.0])
        p /= np.linalg.norm(p)
        return 2 * np.outer(p, p) - np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k / (1 + c)


def _rz(deg):
    t = math.radians(deg)
    return np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1.0]])


def _grid_rotation():
    """Cube coordinates -> grid: body diagonal (1,1,1) up, corner (1,1,-1) seen at 90 deg."""
    q = _rot_between((1, 1, 1), (0, 0, 1))
    p = q @ np.array([1.0, 1.0, -1.0])
    return _rz(90.0 - math.degrees(math.atan2(p[1], p[0]))) @ q


Q = _grid_rotation()
SHADOW = math.sqrt(8 / 3)            # cube corner's distance from the axis, in cube units
UNIT = 4 * R / SHADOW                # one cube unit in mm: half the grid cube's edge


def _edges_at(V, d, tol=1e-6):
    return [(i, j) for i, j in combinations(range(len(V)), 2)
            if abs(np.linalg.norm(V[i] - V[j]) - d) < tol * max(1.0, d)]


def _distance_classes(V):
    ds = sorted({round(float(np.linalg.norm(V[i] - V[j])), 6) for i, j in combinations(range(len(V)), 2)})
    return ds


def _cyclic(t):
    x, y, z = t
    return [(x, y, z), (z, x, y), (y, z, x)]


def _signs(t):
    out = set()
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                out.add((sx * t[0], sy * t[1], sz * t[2]))
    return sorted(out)


def _cube_v():
    return np.array(_signs((1.0, 1.0, 1.0)))


def _tetra_v():
    return np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], float)


def _octa_v():
    return np.array([(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)], float)


def _dodeca_v():
    v = list(_signs((1.0, 1.0, 1.0)))
    for p in _cyclic((0.0, 1 / PHI, PHI)):
        v += _signs(p)
    return np.array(sorted(set(v)))


def _icosa_raw():
    v = []
    for p in _cyclic((0.0, 1.0, PHI)):
        v += _signs(p)
    return np.array(sorted(set(v)))


def _to_grid(V):
    return (Q @ np.asarray(V, float).T).T * UNIT


def _spin_to_rays(G):
    """Spin about the axis so the outermost shadow corners sit on lattice rays."""
    r = np.hypot(G[:, 0], G[:, 1])
    outer = G[np.isclose(r, r.max(), rtol=1e-6)]
    ang = math.degrees(math.atan2(outer[0, 1], outer[0, 0]))
    best = min(RAYS, key=lambda a: abs((ang - a + 180) % 360 - 180))
    return (_rz(best - ang) @ G.T).T


def _icosa_grid():
    G = _to_grid(_icosa_raw())
    G *= 4 * R / np.hypot(G[:, 0], G[:, 1]).max()
    return _spin_to_rays(G)


def _five_tetrahedra(D):
    """Exact cover of the dodecahedron's 20 corners by 5 regular tetrahedra (one hand)."""
    d = 2 * math.sqrt(2) * UNIT
    near = {i: {j for j in range(len(D)) if j != i and abs(np.linalg.norm(D[i] - D[j]) - d) < 1e-6 * d}
            for i in range(len(D))}
    tets = sorted({tuple(sorted(t)) for i in range(len(D)) for t in
                   ([i, a, b, c] for a, b, c in combinations(sorted(near[i]), 3)
                    if b in near[a] and c in near[a] and c in near[b])})

    def cover(chosen, used):
        if len(chosen) == 5:
            return chosen
        first = min(set(range(len(D))) - used)
        for t in tets:
            if first in t and not used & set(t):
                got = cover(chosen + [t], used | set(t))
                if got:
                    return got
        return None
    return cover([], set())


def solid(name):
    """(V, E) in mm, grid orientation: three-fold axis along z, centred on the origin."""
    if name == "icosahedron":
        I = _icosa_grid()
        return I, _edges_at(I, _distance_classes(I)[0])
    C, T, O, D = (_to_grid(f()) for f in (_cube_v, _tetra_v, _octa_v, _dodeca_v))
    if name == "tetrahedron":
        return T, _edges_at(T, 2 * math.sqrt(2) * UNIT)
    if name == "cube":
        return C, _edges_at(C, 2 * UNIT)
    if name == "octahedron":
        return O, _edges_at(O, 2 * math.sqrt(2) * UNIT)
    if name == "dodecahedron":
        return D, _edges_at(D, 2 / PHI * UNIT)
    if name == "five_tetrahedra":                   # drawing data; the print is woven beams
        return D, [e for t in _five_tetrahedra(D) for e in combinations(t, 2)]
    raise ValueError(name)


def _polygon_faces(V):
    """Convex hull of V as polygons: (loop, outward unit n, d) with loop CCW from outside.
    Same construction as stellated_platonics.polygonal_faces: group the hull's triangles
    by plane, then refit each normal from all of the face's vertices."""
    hull = ConvexHull(V)
    groups = {}
    for eq, simplex in zip(hull.equations, hull.simplices):
        groups.setdefault(tuple(np.round(eq / np.linalg.norm(eq[:3]), 6)), set()).update(int(i) for i in simplex)
    out = []
    for key, idx in groups.items():
        idx = sorted(idx)
        pts = V[idx]
        c = pts.mean(axis=0)
        n = np.linalg.svd(pts - c)[2][-1]
        if n @ np.array(key[:3]) < 0:
            n = -n
        u = pts[0] - c
        u -= (u @ n) * n
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        loop = [idx[k] for k in np.argsort(np.arctan2((pts - c) @ w, (pts - c) @ u))]
        out.append((loop, n, float(n @ c)))
    return out


def _merge_collinear(V, E):
    """Join struts that continue one another along a straight line."""
    groups = {}
    for i, j in E:
        d = V[j] - V[i]
        d = d / np.linalg.norm(d)
        if d[np.argmax(np.abs(d))] < 0:
            d = -d
        foot = V[i] - (V[i] @ d) * d
        key = tuple(np.round(np.hstack([d, foot]), 4))
        groups.setdefault(key, (d, set()))[1].update((i, j))
    out = []
    for d, idx in groups.values():
        idx = sorted(idx, key=lambda k: V[k] @ d)
        out.append((idx[0], idx[-1]))
    return out


def stellation(name):
    """First stellation of a grid solid as (V, E, true), or its dual compound (true=False)."""
    V, E = solid(name)
    faces = _polygon_faces(V)
    owner = {}
    for fi, (loop, _, _) in enumerate(faces):
        for a, b in zip(loop, loop[1:] + loop[:1]):
            owner.setdefault(frozenset((a, b)), []).append(fi)
    nb = [set() for _ in faces]
    for pair in owner.values():
        nb[pair[0]].add(pair[1])
        nb[pair[1]].add(pair[0])
    apexes = []
    for fi, (loop, n, _) in enumerate(faces):     # apex = c + h n, h from each neighbour plane
        c = V[loop].mean(axis=0)
        hs = [(faces[fj][2] - c @ faces[fj][1]) / (n @ faces[fj][1]) if abs(n @ faces[fj][1]) > 1e-9 else -1.0
              for fj in nb[fi]]
        assert max(hs) - min(hs) < 1e-6 * max(1.0, abs(hs[0])), "neighbour planes disagree"
        apexes.append(c + hs[0] * n if hs[0] > 1e-9 else None)
    n0 = len(V)
    if all(a is not None for a in apexes):
        W = np.vstack([V, apexes])
        lateral = [(i, n0 + fi) for fi, (loop, _, _) in enumerate(faces) for i in loop]
        return W, _merge_collinear(W, E + lateral), True
    rho = float(np.mean([np.linalg.norm((V[i] + V[j]) / 2) for i, j in E]))   # midradius
    dual = np.array([(rho * rho / d) * n for _, n, d in faces])          # the midsphere dual
    dual_e = [(n0 + a, n0 + b) for a in range(len(faces)) for b in nb[a] if a < b]
    return np.vstack([V, dual]), E + dual_e, False


SOLIDS = ("tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron")
STELLATED = {"tetrahedron": "stella octangula (dual compound; the tetrahedron does not stellate)",
             "cube": "compound of cube and octahedron (dual compound; the cube does not stellate)",
             "octahedron": "stella octangula (true first stellation)",
             "dodecahedron": "small stellated dodecahedron (true first stellation)",
             "icosahedron": "small triambic icosahedron (true first stellation)"}


def _strut(p, q, r):
    p, q = np.asarray(p, float), np.asarray(q, float)
    L = float(np.linalg.norm(q - p))
    m = np.hstack([_rot_between((0, 0, 1), q - p), p[:, None]])
    return m3.Manifold.cylinder(L, r, r, SEGS).transform(m)


def frame(V, E):
    parts = [m3.Manifold.sphere(NODE_D / 2, SEGS).translate(tuple(v)) for v in V]
    parts += [_strut(V[i], V[j], STRUT_D / 2) for i, j in E]
    return m3.Manifold.batch_boolean(parts, m3.OpType.Add)


def _closest(p, q, r, s):
    """Closest points between segments pq and rs."""
    d1, d2, w = q - p, s - r, p - r
    a, e, f = d1 @ d1, d2 @ d2, d2 @ w
    b, c = d1 @ d2, d1 @ w
    den = a * e - b * b
    t = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-12 else 0.0
    u = np.clip((b * t + f) / e, 0, 1)
    t = np.clip((b * u - c) / a, 0, 1)
    return p + d1 * t, r + d2 * u


def strut_gaps(V, E):
    """(centre-line distance, point, point) for every pair of struts sharing no corner."""
    out = []
    for (i, j), (k, l) in combinations(E, 2):
        if not {i, j} & {k, l}:
            p, q = _closest(V[i], V[j], V[k], V[l])
            out.append((float(np.linalg.norm(p - q)), p, q))
    return out


def min_strut_gap(V, E):
    return min((g for g, _, _ in strut_gaps(V, E)), default=math.inf)



# ---------------------------------------------------------------- sphere clusters

SPHERE_STEP = R * math.sqrt(1.5)     # cubic spacing whose diagonal shadow has spacing R


def sphere_centres(n):
    g = (np.arange(n) - (n - 1) / 2) * SPHERE_STEP
    return np.array([(x, y, z) for x in g for y in g for z in g])


def spheres(n):
    """n^3 spheres of radius R, cube-aligned (print orientation). Each cubic cell's centre
    lies sqrt(3)/2 * SPHERE_STEP = 1.06 R from its eight spheres, leaving a sealed 0.9 mm
    bubble; a small filler sphere there closes it without touching the outside."""
    parts = [m3.Manifold.sphere(R, 96).translate(tuple(c)) for c in sphere_centres(n)]
    cells = sphere_centres(n - 1) if n > 1 else []
    parts += [m3.Manifold.sphere(0.12 * R, 24).translate(tuple(c)) for c in cells]
    return m3.Manifold.batch_boolean(parts, m3.OpType.Add)


# ---------------------------------------------------------------- printing and checks

def to_trimesh(man):
    mm = man.to_mesh()
    return trimesh.Trimesh(np.asarray(mm.vert_properties)[:, :3], np.asarray(mm.tri_verts), process=True)


def grid_view(name):
    """The piece as it sits on its plate: three-fold axis up, centred."""
    if name in ("egg_of_life", "flower_of_life_spheres"):
        man = spheres(2 if name == "egg_of_life" else 3)
        return man.transform(np.hstack([Q, np.zeros((3, 1))]))
    if name.startswith("stellated_"):
        V, E, _ = stellation(name[len("stellated_"):])
    else:
        V, E = solid(name)
    return frame(V, E)


def print_pose(man, pts, flat=BED_FLAT):
    """Rest the piece on the largest facet of its corners' hull (pts: corner or sphere
    centres, in the manifold's frame), then trim a footprint off the bottom."""
    pts = np.asarray(pts, float)

    def area(loop):
        p = pts[loop] - pts[loop].mean(axis=0)
        return 0.5 * float(np.linalg.norm(np.cross(p, np.roll(p, -1, axis=0)).sum(axis=0)))
    loop, n, _ = max(_polygon_faces(pts), key=lambda f: area(f[0]))
    rot = _rot_between(n, (0, 0, -1))
    man = man.transform(np.hstack([rot, np.zeros((3, 1))]))
    lo = man.bounding_box()[2]
    man = man.translate((0.0, 0.0, -lo))
    if flat > 0:
        man = man.trim_by_plane((0.0, 0.0, 1.0), flat).translate((0.0, 0.0, -flat))
    b = man.bounding_box()
    return man.translate((-(b[0] + b[3]) / 2, -(b[1] + b[4]) / 2, 0.0))


def decimated(mesh, faces):
    if len(mesh.faces) <= faces:
        return mesh
    import pyvista as pv
    f = np.hstack([np.full((len(mesh.faces), 1), 3), mesh.faces]).ravel()
    pd = pv.PolyData(mesh.vertices, f).decimate(1.0 - faces / len(mesh.faces))
    return trimesh.Trimesh(pd.points, pd.faces.reshape(-1, 4)[:, 1:])


def contact_sheet(items, out, cols=3, az=35, el=24, tile=420, faces=12000):
    """Studio renders of (label, mesh) pairs on one sheet, via lib/render_mesh."""
    import os
    import tempfile
    from PIL import Image, ImageDraw
    import render_mesh as rm
    rows = (len(items) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * tile, rows * (tile + 28)), rm.BG)
    draw = ImageDraw.Draw(sheet)
    tmp = os.path.join(tempfile.gettempdir(), "sacred_grid_tile.png")
    for k, (label, mesh) in enumerate(items):
        m = decimated(mesh, faces)
        rm.render([(m.vertices, m.faces, "frame")], az, el, w=tile, h=tile, out=tmp)
        x, y = (k % cols) * tile, (k // cols) * (tile + 28)
        sheet.paste(Image.open(tmp), (x, y))
        draw.text((x + 12, y + tile + 6), label, fill=(60, 64, 72))
    sheet.save(out)
    print("wrote", out)


def check(mesh, label):
    bodies = len(mesh.split(only_watertight=False))
    ok = mesh.is_watertight and mesh.is_winding_consistent and bodies == 1 and bool((mesh.extents <= BED).all())
    print(f"{'ok ' if ok else 'BAD'} {label:44s} {' x '.join(f'{e:6.1f}' for e in mesh.extents)} mm  "
          f"bodies={bodies} watertight={mesh.is_watertight} vol={mesh.volume / 1000:7.1f} cm3  tris={len(mesh.faces)}")
    return ok
