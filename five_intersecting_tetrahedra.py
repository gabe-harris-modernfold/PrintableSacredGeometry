#!/usr/bin/env python3
"""Five intersecting tetrahedra -- printable 150 mm compound of woven tetrahedral frames.

Geometry
--------
The 20 vertices of a dodecahedron split into five sets of four, each set a regular
tetrahedron; the five tetrahedra together are the classic chiral compound. Turning each
solid tetrahedron into a *frame* (six beams along its edges) lets the five pass through
one another without intersecting -- the woven object built in modular origami from Francis
Ow's 60-degree unit.

Each beam is the intersection of six half-spaces: the four face planes of its tetrahedron
(the two faces meeting at the beam's edge form its outer surfaces, the two opposite faces
miter its ends) plus one inner plane. That gives a triangular-section wedge with its ridge
on the tetrahedron edge and a flat inner face -- the same cross-section the folded paper
unit produces. The single free parameter is `m`, the width of the strip the beam occupies
on each face.

`m` is solved for, not assumed: an LP over every pair of beams belonging to different
tetrahedra finds the width at which the frames jam. That limit is m* = L/12.09 (L = the
tetrahedron edge), which is why Ow's unit folded from 1x3 paper -- frame width L/12 --
is the width that works. Backing off from m* by a printer clearance gives the woven STL;
pushing past it fuses the five frames into one rigid solid.

Outputs
-------
five_intersecting_tetrahedra_150mm.stl        five free interlocked frames (print-in-place)
five_intersecting_tetrahedra_150mm_rigid.stl  the same compound welded into one solid
five_intersecting_tetrahedra.png              preview
"""

import itertools
import os
import tempfile

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.spatial.transform import Rotation

# ---------------------------------------------------------------- parameters

SIZE = 150.0        # largest bounding-box dimension of the finished object, mm
GAP = 0.8           # clearance between neighbouring frames in the woven version, mm
WELD = 1.6          # interference between frames in the rigid version, mm
FIVEFOLD_UP = True  # stand the compound on a dodecahedron face (five vertices on the bed)

PHI = (1 + 5 ** 0.5) / 2
EDGES = list(itertools.combinations(range(4), 2))
FACES = list(itertools.combinations(range(4), 3))

# The two chiral partitions of the 20 dodecahedron vertices into five regular tetrahedra
# (enumerated in five_tetrahedra combinatorics; there are exactly ten such tetrahedra and
# exactly two partitions, mirror images of one another). Index into dodeca_vertices().
TETRAHEDRA = [(0, 3, 5, 6), (1, 8, 12, 19), (2, 9, 16, 17),
              (4, 10, 11, 18), (7, 13, 14, 15)]

FRAME_COLORS = ['#d64545', '#3f8f4f', '#3b6fb6', '#d99b2c', '#8b5cb8']


# ---------------------------------------------------------------- geometry

def dodeca_vertices():
    """The 20 vertices of a dodecahedron: circumradius sqrt(3), bounding box 2*phi."""
    v = [(a, b, c) for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)]
    for s in (-1, 1):
        for t in (-1, 1):
            v.append((0, s / PHI, t * PHI))
            v.append((s / PHI, t * PHI, 0))
            v.append((s * PHI, 0, t / PHI))
    return np.array(v, float)


def tet_face_halfspaces(P):
    """Outward unit normals and offsets of the four faces of tetrahedron P (A x <= b)."""
    centre = P.mean(0)
    A, b = [], []
    for tri in FACES:
        p0, p1, p2 = P[list(tri)]
        n = np.cross(p1 - p0, p2 - p0)
        n /= np.linalg.norm(n)
        if n @ (p0 - centre) < 0:
            n = -n
        A.append(n)
        b.append(n @ p0)
    return A, b


def beam_halfspaces(P, i, j, m):
    """Half-spaces of the frame beam along edge (i, j) of tetrahedron P.

    Four tetrahedron faces plus one inner plane. The inner plane is parallel to the edge
    and passes through the two lines lying at in-plane distance `m` from the edge on the
    two adjacent faces, so the section is a triangle of ridge-to-inner-face depth
    m * cos(dihedral/2) = 0.8165 m.
    """
    A, b = tet_face_halfspaces(P)
    e = P[j] - P[i]
    e = e / np.linalg.norm(e)
    perps = []
    for k in range(4):                       # the two faces containing both i and j
        if k in (i, j):
            continue
        u = P[k] - P[i]
        u = u - (u @ e) * e                  # in-plane, perpendicular to the edge, inward
        perps.append(u / np.linalg.norm(u))
    n = perps[0] + perps[1]
    n /= np.linalg.norm(n)
    A.append(n)
    b.append(n @ P[i] + m * (n @ perps[0]))
    return np.array(A), np.array(b)


def all_beams(V, m):
    """(frame index, A, b) for all 30 beams of the compound."""
    return [(t, *beam_halfspaces(V[list(tet)], i, j, m))
            for t, tet in enumerate(TETRAHEDRA) for (i, j) in EDGES]


# ---------------------------------------------------------------- width solver

def pair_slack(A1, b1, A2, b2):
    """max s over x with a.x + s <= b for both polytopes, all |a| = 1.

    s > 0 means the two beams overlap; s < 0 means they are separated and the true
    distance between them is at least -2s (offsetting half-spaces outward by t grows a
    polytope to at least its radius-t Minkowski sum, so the bound is conservative).
    """
    A = np.vstack([A1, A2])
    b = np.concatenate([b1, b2])
    res = linprog(c=[0, 0, 0, -1], A_ub=np.hstack([A, np.ones((len(A), 1))]), b_ub=b,
                  bounds=[(None, None)] * 4)
    if res.status != 0:
        raise RuntimeError(f"clearance LP failed: {res.message}")
    return res.x[3]


def worst_gap(V, m):
    """Signed clearance between the closest pair of beams from different frames.

    Positive: guaranteed gap. Negative: depth of interference.
    """
    beams = all_beams(V, m)
    s = max(pair_slack(A1, b1, A2, b2)
            for (ta, A1, b1), (tb, A2, b2) in itertools.combinations(beams, 2)
            if ta != tb)
    return -2 * s


def solve_width(V, target_gap, lo=None, hi=None, iters=34):
    """Beam width `m` whose closest cross-frame clearance equals `target_gap`."""
    L = np.linalg.norm(V[TETRAHEDRA[0][0]] - V[TETRAHEDRA[0][1]])
    lo = L / 60 if lo is None else lo
    hi = L / 6 if hi is None else hi
    for _ in range(iters):
        mid = (lo + hi) / 2
        if worst_gap(V, mid) > target_gap:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------- solid construction

def polytope_faces(A, b, tol=1e-7):
    """Vertices and outward-wound face loops of the bounded polytope {A x <= b}."""
    res = linprog(c=[0, 0, 0, -1], A_ub=np.hstack([A, np.ones((len(A), 1))]), b_ub=b,
                  bounds=[(None, None)] * 4)
    centre, inradius = res.x[:3], res.x[3]
    if inradius <= tol:
        raise RuntimeError(f"degenerate beam (inradius {inradius:.3g})")
    verts = HalfspaceIntersection(np.hstack([A, -b[:, None]]), centre).intersections

    loops = []
    for n, off in zip(A, b):
        on = np.flatnonzero(np.abs(verts @ n - off) < 1e-6)
        if len(on) < 3:
            continue
        u = np.cross(n, [0, 0, 1.0])
        if np.linalg.norm(u) < 0.1:
            u = np.cross(n, [1.0, 0, 0])
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        d = verts[on] - verts[on].mean(0)
        loops.append(on[np.argsort(np.arctan2(d @ w, d @ u))])   # CCW seen from outside
    return verts, loops


def make_solid(A, b):
    """One beam as an exact OCCT solid built from its planar faces."""
    from build123d import Face, Shell, Solid, Vector, Wire
    verts, loops = polytope_faces(A, b)
    faces = [Face(Wire.make_polygon([Vector(*verts[k]) for k in loop], close=True))
             for loop in loops]
    solid = Solid(Shell(faces))
    if solid.volume <= 0:
        raise RuntimeError("beam solid has non-positive volume")
    return solid


def build_frames(V, m):
    """One fused OCCT solid per tetrahedral frame."""
    frames = []
    for tet in TETRAHEDRA:
        P = V[list(tet)]
        beams = [make_solid(*beam_halfspaces(P, i, j, m)) for (i, j) in EDGES]
        frame = beams[0]
        for beam in beams[1:]:
            frame = frame + beam
        frames.append(frame.clean())
    return frames


# ---------------------------------------------------------------- placement

def orientation():
    """Rotation standing the compound on a dodecahedron face.

    A face normal is a five-fold axis, so pointing one straight down puts that face's five
    vertices -- five frame tips, one from each tetrahedron -- on the bed in a pentagon.
    All twelve faces are equivalent under the compound's rotation group, so any will do.
    """
    if not FIVEFOLD_UP:
        return np.eye(3)
    V = dodeca_vertices()
    normals = np.unique(np.round(ConvexHull(V).equations[:, :3], 6), axis=0)
    assert len(normals) == 12, f"expected 12 dodecahedron faces, got {len(normals)}"
    R, _ = Rotation.align_vectors([[0, 0, -1.0]], [normals[0]])
    return R.as_matrix()


def placed_vertices():
    """Dodecahedron vertices rotated to the print orientation and scaled to SIZE."""
    V = dodeca_vertices() @ orientation().T
    V *= SIZE / (V.max(0) - V.min(0)).max()
    return V


# ---------------------------------------------------------------- output

def export(frames, path, fuse):
    """Export frames as one STL; `fuse` welds them into a single connected solid."""
    import trimesh
    from build123d import export_stl

    meshes = []
    with tempfile.TemporaryDirectory() as tmp:
        parts = frames
        if fuse:
            whole = parts[0]
            for p in parts[1:]:
                whole = whole + p
            parts = [whole.clean()]
        for k, part in enumerate(parts):
            f = os.path.join(tmp, f"part{k}.stl")
            export_stl(part, f, tolerance=1e-3, angular_tolerance=0.1)
            meshes.append(trimesh.load(f))
    combined = trimesh.util.concatenate(meshes)
    combined.export(path)
    return meshes, combined


def render(frame_meshes, path):
    """Two-view preview, one colour per tetrahedron (z-buffered, so the weave reads)."""
    import pyvista as pv
    pv.OFF_SCREEN = True

    centre = np.concatenate([m.vertices for m in frame_meshes]).mean(0)
    d = 3.1 * SIZE
    views = [("the woven compound", np.array([0.86, -0.98, 0.62])),
             ("down a five-fold axis", np.array([0.0, 0.0, 1.0]))]

    plotter = pv.Plotter(shape=(1, 2), window_size=(1560, 820),
                         off_screen=True, border=False)
    for col, (title, direction) in enumerate(views):
        plotter.subplot(0, col)
        for mesh, color in zip(frame_meshes, FRAME_COLORS):
            plotter.add_mesh(pv.wrap(mesh), color=color, smooth_shading=False,
                             specular=0.25, specular_power=12, ambient=0.22)
        eye = centre + d * direction / np.linalg.norm(direction)
        up = (0, 1, 0) if direction[2] > 0.99 else (0, 0, 1)
        plotter.camera_position = [tuple(eye), tuple(centre), up]
        plotter.camera.zoom(1.35)
        plotter.add_text(title, position='lower_edge', font_size=13, color='#333333')
    plotter.set_background('white')
    plotter.screenshot(path)
    plotter.close()


# ---------------------------------------------------------------- interlock check

def linking_number(loop_a, loop_b, view):
    """Linking number of two closed polylines, from signed crossings in projection."""
    u = np.cross(view, [1.0, 0, 0])
    if np.linalg.norm(u) < 0.3:
        u = np.cross(view, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    basis = np.array([u, np.cross(view, u)])

    def segments(loop):
        return [(loop[i], loop[(i + 1) % len(loop)]) for i in range(len(loop))]

    def cross2(a, b):
        return a[0] * b[1] - a[1] * b[0]

    total = 0.0
    for p0, p1 in segments(loop_a):
        for q0, q1 in segments(loop_b):
            a0, a1, b0, b1 = (basis @ p for p in (p0, p1, q0, q1))
            r, s = a1 - a0, b1 - b0
            den = cross2(r, s)
            if abs(den) < 1e-9:
                raise RuntimeError("degenerate projection direction")
            t, v = cross2(b0 - a0, s) / den, cross2(b0 - a0, r) / den
            if not (1e-9 < t < 1 - 1e-9 and 1e-9 < v < 1 - 1e-9):
                continue
            za = (p0 + t * (p1 - p0)) @ view
            zb = (q0 + v * (q1 - q0)) @ view
            total += np.sign(den) if za > zb else -np.sign(den)
    return total / 2


def check_interlock(V):
    """Every pair of frames must be linked, or the compound would fall apart.

    Each frame carries four triangular cycles (its faces). Two frames are inseparable if
    any cycle of one has non-zero linking number with a cycle of the other.
    """
    view = np.array([0.31, 0.57, 0.76])
    view /= np.linalg.norm(view)
    cycles = [(t, V[list(tri)]) for t, tet in enumerate(TETRAHEDRA)
              for tri in itertools.combinations(tet, 3)]
    linked = {}
    for (ta, la), (tb, lb) in itertools.combinations(cycles, 2):
        if ta != tb and linking_number(la, lb, view) != 0:
            linked[ta, tb] = linked.get((ta, tb), 0) + 1
    pairs = list(itertools.combinations(range(len(TETRAHEDRA)), 2))
    missing = [p for p in pairs if p not in linked]
    if missing:
        raise RuntimeError(f"frames not interlocked: {missing}")
    return len(pairs), min(linked.values())


# ---------------------------------------------------------------- main

def main():
    V = placed_vertices()
    L = np.linalg.norm(V[TETRAHEDRA[0][0]] - V[TETRAHEDRA[0][1]])
    m_max = solve_width(V, 0.0)
    n_pairs, n_cycles = check_interlock(V)
    print(f"tetrahedron edge      {L:.3f} mm")
    print(f"jamming beam width    m* = {m_max:.3f} mm = L/{L / m_max:.3f}"
          f"   (origami 1x3 unit: L/12)")
    print(f"interlock             all {n_pairs} frame pairs linked"
          f" ({n_cycles}+ linked cycle pairs each)")

    for tag, target, fuse in [("woven", GAP, False), ("rigid", -WELD, True)]:
        m = solve_width(V, target)
        gap = worst_gap(V, m)
        print(f"\n--- {tag} ---")
        print(f"beam width            {m:.3f} mm on each face  = L/{L / m:.2f}"
              f"   (depth {0.81650 * m:.3f} mm)")
        print(f"closest frame pair    {gap:+.3f} mm")

        frames = build_frames(V, m)
        suffix = "" if tag == "woven" else "_rigid"
        path = f"five_intersecting_tetrahedra_{SIZE:.0f}mm{suffix}.stl"
        meshes, combined = export(frames, path, fuse)

        lo = combined.bounds[0]
        combined.apply_translation([0, 0, -lo[2]])           # sit on the bed
        combined.export(path)
        ext = combined.extents
        print(f"{path}")
        print(f"  bodies {len(combined.split(only_watertight=False))}"
              f"  triangles {len(combined.faces)}"
              f"  watertight {combined.is_watertight}")
        print(f"  extents {ext[0]:.2f} x {ext[1]:.2f} x {ext[2]:.2f} mm"
              f"   volume {combined.volume / 1000:.1f} cm^3")
        feet = combined.vertices[combined.vertices[:, 2] < 0.05]
        span = np.linalg.norm(feet[:, :2] - feet[:, :2].mean(0), axis=1)
        print(f"  bed contact: {len(feet)} frame tips on a"
              f" {2 * span.max():.1f} mm circle")
        if tag == "woven":
            for k, fm in enumerate(meshes):
                assert fm.is_watertight, f"frame {k} not watertight"
            print(f"  every frame watertight, volume/frame"
                  f" {meshes[0].volume / 1000:.2f} cm^3")
            render([fm.copy().apply_translation([0, 0, -lo[2]]) for fm in meshes],
                   "five_intersecting_tetrahedra.png")
            print("  preview -> five_intersecting_tetrahedra.png")

        bed = 320.0
        assert (ext <= bed).all(), f"does not fit the {bed:.0f} mm bed"


if __name__ == "__main__":
    main()
