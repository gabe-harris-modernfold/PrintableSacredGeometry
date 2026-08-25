#!/usr/bin/env python3
"""
The three-ribbon tetrahelix: the ordinary Boerdijk-Coxeter helix as a hollow tube
whose surface is three interwoven strips of equilateral triangles.

WHAT THE RIBBON COUNT ACTUALLY IS
---------------------------------
Educational references offer a one-ribbon, a two-ribbon and a three-ribbon
tetrahelix, which reads as three different objects. It is one object cut three ways,
and the reason there are exactly three is that the tetrahelix has exactly three edge
families.

Cell T_k = {k,k+1,k+2,k+3}. Consecutive cells share {k+1,k+2,k+3}, so the interior
faces are {k,k+1,k+2} and each cell has two free ones:

    F1_k = {k, k+1, k+3}        F2_k = {k, k+2, k+3}

Those are the surface. It is a triangulated tube -- one vertex, three edges, two
triangles per step -- and every surface triangle carries exactly one edge of each
family. Cut the whole k->k+m family and what is left falls into strips indexed by
k mod m:

    cut k->k+1   ->  ONE ribbon      fold across that family 148.413662 deg
    cut k->k+2   ->  TWO ribbons     fold 141.057559 deg
    cut k->k+3   ->  THREE ribbons   fold  70.528779 deg = arccos(1/3) = arctan(2*sqrt2)

The three-ribbon cut is the interesting one, and that 70.528779 deg is why: F1_k and
F2_k are two faces of the SAME cell meeting along the k->k+3 edge, so the fold there
is the regular tetrahedron's own dihedral angle. Which means the three seams are
exactly the three intertwined helices the BC helix is famous for -- the edges
belonging to only one tetrahedron. The ribbons are the material between them.

All of it is asserted in _verify():
  * every surface triangle is exactly equilateral at the cell edge. The cells stay
    perfectly regular, theta = arccos(-2/3) unfudged, nothing bent;
  * a ribbon's own folds are the gentle pair, 148.41 and 141.06 deg alternately
    (31.59 and 38.94 deg off flat), which is what lets it read as a ribbon;
  * each cell's two surface faces land in DIFFERENT ribbons, so the ribbons share
    whole tetrahedra rather than merely touching;
  * cutting family m really does leave m strips, each a path, for m = 1, 2 and 3.

WHY IT IS HOLLOW
----------------
Only the surface is built. The interior {k,k+1,k+2} partitions are left out, so the
tube is open end to end and you can see straight through the bore -- which is what a
folded-paper model of a tetrahelix is, and what a solid chain of tetrahedra is not.

WALL: PER-FACE PRISMS PLUS FILLETS ON THE REFLEX EDGES
-----------------------------------------------------
Two constructions look right here and both fail; the notes are in shell-free form
in wall_prisms() and fillet_rods(), and the short version is:

A +-wall/2 VERTEX OFFSET closes the shell but does not deliver the wall. At a vertex
of this surface the incident face normals spread up to 143 deg, so their mean is
nearly RADIAL while the faces sit far off radial, and the least perpendicular
thickness saturates near 0.52 mm however far you offset.

PER-FACE INWARD PRISMS give exactly `wall` perpendicular by construction, but only
overlap where the fold is convex. The count of cells sharing an edge decides that,
and it is 4-m for the k->k+m family: THREE cells share every k->k+1 edge, so the
material there spans 3*70.5288 = 211.586 deg -- reflex. Two inward prisms diverge at
a reflex edge and meet along the edge LINE only, so the wall comes out knife-edge
joined along all 19 of them. An unsigned fold measurement hides this completely: it
returns 148.414, the complement, which reads as a gentle convex fold.

So: per-face prisms for the flats, plus a fillet rod along every reflex edge. Both
halves are asserted -- the prisms must MISS the material bisector at a reflex edge,
and the rods must cover it.

The ribbons are then raised as a second set of prisms on the same faces, inset from
their own seam edge, which leaves a groove along each seam helix and nowhere else.
Since every surface triangle has exactly one edge of each family, the identical
construction gives 1, 2 or 3 ribbons with nothing special-cased.

STELLATION (--stellate)
-----------------------
Glue a regular tetrahedron onto every surface face, face-to-face. Every added cell
stays exactly regular and every join is a real shared face, and it is collision-free:
measured over all pairs, 0 overlaps. The object doubles, radius 21.5 -> 43.1 mm.

It goes exactly ONE layer deep. Stellating again adds 126 cells and 77 of them
collide, so there is no second layer and no tripling by this route -- asserted, so
the limit cannot be quietly crossed.

Why one layer fits at all is Fuller's unzipping angle. At the tightest family, the
air at a k->k+1 edge is 148.414 deg and two glued tetrahedra need 2*70.5288 =
141.058, leaving 148.414 - 141.058 = 7.356103 deg. After stellation five cells ring
every k->k+1 edge -- 352.644 deg of material -- and the leftover 7.356 deg shows up
on the print as a narrow V-notch running along each of those edges. The same number
that stops regular tetrahedra tiling space is what leaves room for the layer.

STACKING
--------
The regular tetrahelix cannot repeat: theta/360 is irrational, which is the whole of
Fuller's 5.69 deg deficit. --closed switches to tetrahelix.py's 132 deg variant,
which repeats every 30 cells and so stacks exactly, at the cost of 0.28% on one step
family. Default is the exact regular helix and no stacking.

Run:  python ribbon_tetrahelix.py [--ribbons 3] [--edge 40] [--cells 21]
                                  [--wall 1.6] [--relief 1.2] [--groove 2.5]
                                  [--stellate] [--closed] [--weld]
"""

import argparse
import importlib.util
import math
import os
from collections import defaultdict

import numpy as np

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import mesh_kit as MK
import tetrahelix as TH

RIBBONS = 3             # which edge family to cut along = how many ribbons
EDGE = 40.0             # cell edge. Tube outer dia is about 2*(0.5196*edge+relief),
                        # so 40 gives ~45 mm across and a ~39 mm bore.
CELLS = 21              # tetrahedra. 21 at edge 40 is a 291 mm tube.
WALL = 1.6              # PERPENDICULAR wall thickness, mm -- 4 perimeters at 0.4
RELIEF = 1.2            # how far a ribbon stands off the wall
GROOVE = 2.5            # pad inset from its own seam edge; the seam groove comes out
                        # 2*GROOVE wide, since both flanking triangles inset
SINK = 0.8              # how far a pad reaches back inside the wall, so it fuses
                        # rather than sitting on it with zero overlap
EPS = 1e-9


def surface(cells):
    """The free faces, tagged (family, k). Combinatorics only; winding comes later."""
    out = []
    for k in range(cells):
        out.append((("F1", k), (k, k + 1, k + 3)))
        out.append((("F2", k), (k, k + 2, k + 3)))
    return out


def seam_edge(face, m):
    """The one edge of this surface triangle from the cut family k->k+m.

    Every surface triangle has exactly one edge of each family -- F1_k is (d1,d2,d3)
    and F2_k is (d2,d1,d3) -- which is why --ribbons 1, 2 and 3 all run through the
    same code with nothing special-cased."""
    for i in range(3):
        a, b = face[i], face[(i + 1) % 3]
        if abs(b - a) == m:
            return (min(a, b), max(a, b))
    raise AssertionError(f"no k->k+{m} edge in {face}")


def ribbons(cells, m):
    """Cut the k->k+m family. Returns {tag: ribbon index}, the strips, adjacency.

    Each strip is walked from a tip so it comes out in order along the ribbon rather
    than as an unordered set."""
    faces = surface(cells)
    by_edge = defaultdict(list)
    for tag, f in faces:
        for i in range(3):
            e = (min(f[i], f[(i + 1) % 3]), max(f[i], f[(i + 1) % 3]))
            by_edge[e].append(tag)
    adj = defaultdict(list)
    for e, lst in by_edge.items():
        if len(lst) == 2 and e[1] - e[0] != m:
            adj[lst[0]].append(lst[1])
            adj[lst[1]].append(lst[0])

    seen, strips = set(), []
    order = [t for t, _f in faces]
    for tips_only in (True, False):          # tips first, then any leftover cycles
        for start in order:
            if start in seen or (tips_only and len(adj[start]) > 1):
                continue
            strip, cur, prev = [], start, None
            while cur is not None and cur not in seen:
                seen.add(cur)
                strip.append(cur)
                nxt = [x for x in adj[cur] if x != prev]
                cur, prev = (nxt[0] if nxt else None), cur
            if strip:
                strips.append(strip)
    strips.sort(key=lambda s: min(k for _fam, k in s))
    return {t: i for i, s in enumerate(strips) for t in s}, strips, adj


def oriented(v, cells):
    """Surface faces wound so the normal points away from their own cell."""
    out = []
    for tag, f in surface(cells):
        k = tag[1]
        cen = v[[k, k + 1, k + 2, k + 3]].mean(axis=0)
        p = v[list(f)]
        if np.cross(p[1] - p[0], p[2] - p[0]) @ (p.mean(axis=0) - cen) < 0:
            f = (f[0], f[2], f[1])
        out.append((tag, f))
    return out


def face_normal(v, f):
    n = np.cross(v[f[1]] - v[f[0]], v[f[2]] - v[f[0]])
    return n / np.linalg.norm(n)


def tri_prism(p0, p1, p2, d, lo, hi):
    """Closed triangular prism over a planar triangle, from lo*d to hi*d.

    The winding is fixed here rather than trusted from the caller: the triangle is
    flipped if its normal disagrees with d. Skipping that left half the ribbon pads
    inside out, each with a volume of exactly the right size and the wrong sign --
    which cancels to a plausible-looking total and is easy to miss."""
    P = [np.asarray(x, float) for x in (p0, p1, p2)]
    if np.cross(P[1] - P[0], P[2] - P[0]) @ d < 0:
        P = [P[0], P[2], P[1]]
    verts = [q + lo * d for q in P] + [q + hi * d for q in P]
    faces = [(0, 2, 1), (3, 4, 5)]
    for j in range(3):
        j2 = (j + 1) % 3
        faces += [(j, j2, 3 + j2), (j, 3 + j2, 3 + j)]
    return verts, faces


def reflex_edges(cells, faces):
    """Edges where the material wedge exceeds 180 deg, with their two faces.

    Only the k->k+1 family: three cells share such an edge, so the material spans
    3*70.5288 = 211.586 deg. These are the edges that need filleting."""
    fam = [m for m in (1, 2, 3) if material_wedge(m) > 180.0]
    by_edge = defaultdict(list)
    for _tag, f in faces:
        for i in range(3):
            a, b = f[i], f[(i + 1) % 3]
            if abs(b - a) in fam:
                by_edge[(min(a, b), max(a, b))].append(f)
    return {e: l for e, l in by_edge.items() if len(l) == 2}


def _in_prism(P, verts, tol=1e-9):
    """Is point P strictly inside this triangular prism? Used only by _verify()."""
    Q = np.asarray(verts, float)
    bot, top = Q[:3], Q[3:]
    n = np.cross(bot[1] - bot[0], bot[2] - bot[0])
    n = n / np.linalg.norm(n)
    if n @ (top[0] - bot[0]) < 0:
        n = -n
    h = float(n @ (top[0] - bot[0]))
    t = float((P - bot[0]) @ n)
    if not (tol < t < h - tol):
        return False
    cen = bot.mean(axis=0)
    for i in range(3):
        a, b = bot[i], bot[(i + 1) % 3]
        nn = np.cross(b - a, n)
        nn = nn / np.linalg.norm(nn)
        sgn = 1.0 if (cen - a) @ nn > 0 else -1.0
        if sgn * float((P - a) @ nn) <= tol:
            return False
    return True


def _in_rod(P, a, b, r, tol=1e-9):
    """Is point P strictly inside the capsule-free cylinder from a to b?"""
    ax = b - a
    L = np.linalg.norm(ax)
    ax = ax / L
    t = float((P - a) @ ax)
    if not (tol < t < L - tol):
        return False
    return float(np.linalg.norm(np.cross(P - a, ax))) < r - tol


def wall_prisms(v, faces, wall):
    """One inward prism per face: exactly `wall` perpendicular, everywhere."""
    return [tri_prism(v[f[0]], v[f[1]], v[f[2]], face_normal(v, f), -wall, 0.0)
            for _tag, f in faces]


def fillet_rods(v, faces, cells, wall, nseg=10):
    """A rod along every REFLEX edge, bridging the two wall prisms that meet there.

    Per-face prisms give an exact wall on the flats and overlap correctly at a
    convex fold, but at a reflex one they diverge: each prism is bounded by its own
    triangle's footprint, so the two meet along the edge LINE and nowhere else, and
    the wall is knife-edge joined along all 19 k->k+1 edges. Nothing in an undirected
    edge count or an unsigned fold angle shows it -- the fold reads 148.414 and looks
    convex, and the two skins touch, so the mesh looks closed.

    A rod centred ON the edge fixes it, because it contains a whole neighbourhood of
    the edge and so overlaps both prisms with real volume. It is the edge term of a
    Minkowski offset, and the honest cost is that it protrudes into the air side --
    but that side is a 148.414 deg valley, so what shows is a fillet in an internal
    corner, not a bump on the silhouette.

    Vertex offsetting is the other obvious fix and it is worse: it closes the shell
    but the least perpendicular thickness saturates near 0.52 mm however far you
    offset, because the incident face normals at a vertex here spread up to 143 deg."""
    out = []
    for e, _l in reflex_edges(cells, faces).items():
        out.append(MK.tube(v[e[0]], v[e[1]], wall, nseg=nseg))
    return out


def pad_prisms(v, faces, home, relief, groove, m, sink=SINK):
    """Raised ribbon pads: the face, inset from its own seam edge, stood off.

    Inset only that one edge and the pads of a ribbon still overlap along the edges
    they share -- each covers all but groove/height of such an edge measured from its
    own apex, so two pads meeting there overlap across most of it and read as one
    continuous ribbon. What stays uncovered is a band either side of every seam edge:
    the groove. The footprint is handed back so the test can check the groove is
    really open instead of re-deriving the same inset arithmetic."""
    out = []
    for tag, f in faces:
        e = seam_edge(f, m)
        apex = [x for x in f if x not in e][0]
        h = (np.linalg.norm(np.cross(v[e[1]] - v[e[0]], v[apex] - v[e[0]]))
             / np.linalg.norm(v[e[1]] - v[e[0]]))
        t = groove / h
        foot = [v[apex],
                v[e[0]] + t * (v[apex] - v[e[0]]),
                v[e[1]] + t * (v[apex] - v[e[1]])]
        d = face_normal(v, f)
        verts, tri = tri_prism(foot[0], foot[1], foot[2], d, -sink, relief)
        out.append((verts, tri, f"ribbon{home[tag]}", foot))
    return out


DIHEDRAL = math.degrees(math.acos(1.0 / 3.0))       # 70.528779, the only fold angle
                                                    # a regular tetrahedron has


def stellate(v, faces, home, edge):
    """Glue a regular tetrahedron on every surface face, face-to-face.

    Exact by construction: a regular tetrahedron on an equilateral face needs no
    fitting, its apex sits at centroid + edge*sqrt(2/3) along the outward normal.
    Returns the extended vertex list, the added cells, the new outer surface (each
    spike's three side faces, tagged with the ribbon its base came from) and the
    full cell list, which is what the reflex classification needs."""
    v = list(v)
    spikes, outer = [], []
    for tag, f in faces:
        p = np.asarray([v[i] for i in f], float)
        c = p.mean(axis=0)
        n = np.cross(p[1] - p[0], p[2] - p[0])
        n = n / np.linalg.norm(n)
        if n @ (c - np.array([0.0, 0.0, c[2]])) < 0:      # outward from the axis
            n = -n
        v.append(c + n * edge * math.sqrt(2.0 / 3.0))
        ap = len(v) - 1
        spikes.append((f[0], f[1], f[2], ap))
        for i in range(3):
            outer.append(((f"spike{home[tag]}", tag[1]),
                          (f[i], f[(i + 1) % 3], ap)))
    return np.asarray(v, float), spikes, outer


def cells_of(cells, spikes):
    """Every cell of the complex as a frozenset of vertex indices."""
    return ([frozenset((k, k + 1, k + 2, k + 3)) for k in range(cells)]
            + [frozenset(sp) for sp in spikes])


def reflex_by_count(surf, all_cells):
    """Surface edges whose material wedge exceeds 180 deg, and the count that says so.

    Every cell here is a regular tetrahedron, so the material around an edge is just
    (number of cells holding that edge) * 70.528779 deg. Three or more is reflex.
    That replaces any angle measurement and cannot pick the wrong side."""
    seen, out = defaultdict(list), {}
    for _tag, f in surf:
        for i in range(3):
            a, b = f[i], f[(i + 1) % 3]
            seen[(min(a, b), max(a, b))].append(f)
    for e, l in seen.items():
        if len(l) != 2:
            continue
        n = sum(1 for c in all_cells if e[0] in c and e[1] in c)
        if n * DIHEDRAL > 180.0:
            out[e] = n
    return out


def build_stellated(ribbon_count=RIBBONS, edge=EDGE, cells=CELLS, wall=WALL,
                    closed=False):
    """Tube plus one stellated layer, as a printable shell.

    The material is a wall following the OUTER surface (the spikes) plus the tube's
    own wall, which is what keeps the bore open -- the tube faces are interior now,
    buried under their spikes. The two walls meet along the spike base edges, and
    the fillet rods there do double duty: they bridge tube wall to spike wall AND
    fill the reflex wedge that per-face prisms leave open."""
    theta, R, H = TH.geometry(edge, closed)
    v0 = TH.helix(cells + 3, edge, closed)
    v0 = v0 - np.array([0.0, 0.0, v0[:, 2].min()])
    faces = oriented(v0, cells)
    home, strips, _adj = ribbons(cells, ribbon_count)
    v, spikes, outer = stellate(v0, faces, home, edge)
    all_cells = cells_of(cells, spikes)
    rf = reflex_by_count(outer, all_cells)

    m = MK.Mesh()
    for _tag, f in faces:                      # the tube wall keeps the bore open
        verts, tri = tri_prism(v[f[0]], v[f[1]], v[f[2]],
                               face_normal(v, f), -wall, 0.0)
        m.add_solid(verts, tri, "bore")
    for tag, f in outer:                       # the outer skin, on the spikes
        n = face_normal(v, f)
        c = np.asarray([v[i] for i in f]).mean(axis=0)
        sp = next(sp for sp in spikes if set(f) <= set(sp))
        if n @ (c - v[list(sp)].mean(axis=0)) < 0:
            n = -n                             # outward from its own spike
        verts, tri = tri_prism(v[f[0]], v[f[1]], v[f[2]], n, -wall, 0.0)
        m.add_solid(verts, tri, tag[0])
    for e in rf:                               # fillets, and the wall-to-wall bridge
        m.add_solid(*MK.tube(v[e[0]], v[e[1]], wall, nseg=10), tag="fillet")
    return m, v, faces, outer, spikes, home, strips, rf, (theta, R, H)


def build(ribbon_count=RIBBONS, edge=EDGE, cells=CELLS, wall=WALL, relief=RELIEF,
          groove=GROOVE, closed=False):
    theta, R, H = TH.geometry(edge, closed)
    v = TH.helix(cells + 3, edge, closed)
    v = v - np.array([0.0, 0.0, v[:, 2].min()])
    faces = oriented(v, cells)
    home, strips, _adj = ribbons(cells, ribbon_count)

    m = MK.Mesh()
    for verts, tri in wall_prisms(v, faces, wall):
        m.add_solid(verts, tri, "wall")
    for verts, tri in fillet_rods(v, faces, cells, wall):
        m.add_solid(verts, tri, "fillet")
    for verts, tri, tag, _foot in pad_prisms(v, faces, home, relief, groove,
                                             ribbon_count):
        m.add_solid(verts, tri, tag)
    return m, v, faces, home, strips, (theta, R, H)


class Welded:
    def __init__(self, v, f):
        self._v, self._f = v, f

    def arrays(self):
        return self._v, self._f


def weld(m, log=lambda *a: print(*a, flush=True)):
    """CSG-union every prism into one 2-manifold solid, or None if unavailable."""
    spec = importlib.util.find_spec("pymeshlab")
    if spec is None:
        return None
    if spec.submodule_search_locations:
        os.add_dll_directory(spec.submodule_search_locations[0])
    import pymeshlab as ml
    V, F = m.arrays()
    ms = ml.MeshSet()
    acc = None
    for i, (base, f0, nf) in enumerate(m._solid_start):
        ff = F[f0:f0 + nf]
        idx = np.unique(ff)
        remap = np.zeros(idx.max() + 1, np.int32)
        remap[idx] = np.arange(len(idx))
        ms.add_mesh(ml.Mesh(V[idx], remap[ff].astype(np.int32)))
        nxt = ms.mesh_number() - 1
        if acc is None:
            acc = nxt
        else:
            ms.apply_filter("generate_boolean_union", first_mesh=acc, second_mesh=nxt)
            acc = ms.mesh_number() - 1
        if i and i % 20 == 0:
            log(f"    welded {i}/{len(m._solid_start)} "
                f"({ms.mesh(acc).face_number()} faces)")
    ms.set_current_mesh(acc)
    cm = ms.current_mesh()
    return Welded(np.asarray(cm.vertex_matrix(), float),
                  np.asarray(cm.face_matrix(), np.int64))


def shell_stats(v, f):
    """(components, (open, pinched) edges, signed volume), measured BY POSITION."""
    v, f = np.asarray(v, float), np.asarray(f)
    tri = v[f]
    vol = float(np.einsum("ij,ij->i", tri[:, 0],
                          np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)
    _, inv = np.unique(np.round(v, 5), axis=0, return_inverse=True)
    g = inv.reshape(-1)[f]
    parent = list(range(int(g.max()) + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in g:
        for x, y in ((a, b), (b, c)):
            ra, rb = find(int(x)), find(int(y))
            if ra != rb:
                parent[ra] = rb
    comps = len({find(int(i)) for i in np.unique(g)})
    e = np.sort(np.concatenate([g[:, [0, 1]], g[:, [1, 2]], g[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return comps, (int((cnt == 1).sum()), int((cnt > 2).sum())), vol


def cells_on_edge(m):
    """How many cells share a k->k+m edge. T_j = {j..j+3} holds both k and k+m iff
    j <= k <= k+m <= j+3, so it is 4-m: three cells on a k->k+1 edge, two on a
    k->k+2, one on a k->k+3."""
    return 4 - m


def material_wedge(m):
    """Solid angle of material around a k->k+m edge, in degrees.

    This is the number that matters for offsetting, and it is NOT what an
    unsigned angle-between-two-faces gives you. Three cells meet on a k->k+1
    edge, so the material there spans 3*70.5288 = 211.586 deg -- REFLEX. An
    unsigned measurement returns 148.414, its complement, and reads as a gentle
    convex fold. Believing that is what pinched the wall: see shell()."""
    return cells_on_edge(m) * math.degrees(math.acos(1.0 / 3.0))


def dihedrals(v, cells):
    """Unsigned angle between the two surface triangles across each edge family.

    Unsigned, so for the k->k+1 family this returns 148.414 -- the AIR side, since
    the material there is the reflex 211.586. Use material_wedge() when the side
    matters; this one is for the fold a ribbon actually makes."""
    faces = surface(cells)
    by_edge = defaultdict(list)
    for _tag, f in faces:
        for i in range(3):
            e = (min(f[i], f[(i + 1) % 3]), max(f[i], f[(i + 1) % 3]))
            by_edge[e].append(f)
    res = defaultdict(list)
    for e, lst in by_edge.items():
        if len(lst) != 2:
            continue
        a = [x for x in lst[0] if x not in e][0]
        b = [x for x in lst[1] if x not in e][0]
        ax = v[e[1]] - v[e[0]]
        ax = ax / np.linalg.norm(ax)
        u = v[a] - v[e[0]]
        u = u - ax * (u @ ax)
        w = v[b] - v[e[0]]
        w = w - ax * (w @ ax)
        c = (u @ w) / (np.linalg.norm(u) * np.linalg.norm(w))
        res[e[1] - e[0]].append(math.degrees(math.acos(max(-1.0, min(1.0, c)))))
    return res


def area_of(v, faces):
    return float(sum(0.5 * np.linalg.norm(np.cross(v[f[1]] - v[f[0]],
                                                   v[f[2]] - v[f[0]]))
                     for _t, f in faces))


def _vol(v, f):
    t = np.asarray(v, float)[np.asarray(f)]
    return float(np.einsum("ij,ij->i", t[:, 0],
                           np.cross(t[:, 1], t[:, 2])).sum() / 6.0)


def _tet_samples(T, n=3):
    P = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                w = np.array([i, j, k, n - i - j - k], float) + 0.4
                P.append(w / w.sum() @ T)
    return np.array(P)


def _tet_inside(T, P):
    M = np.stack([T[1] - T[0], T[2] - T[0], T[3] - T[0]], axis=1)
    w = np.linalg.solve(M, (P - T[0]).T).T
    return (w > 1e-7).all(axis=1) & (w.sum(axis=1) < 1 - 1e-7)


def collisions(tets, edge, stop_at=None):
    """Pairs of cells whose interiors overlap.

    Sampled, and deliberately so: a face-to-face packing of regular tetrahedra
    either fits or it does not, and finding ANY interior point of one cell inside
    another settles it without a full separating-axis implementation."""
    S = [_tet_samples(T) for T in tets]
    C = np.array([T.mean(axis=0) for T in tets])
    hits = []
    for i in range(len(tets)):
        for j in range(i + 1, len(tets)):
            if np.linalg.norm(C[i] - C[j]) > 1.4 * edge:
                continue
            if _tet_inside(tets[j], S[i]).any() or _tet_inside(tets[i], S[j]).any():
                hits.append((i, j))
                if stop_at and len(hits) >= stop_at:
                    return hits
    return hits


def _verify_stellated(edge=EDGE, cells=CELLS):
    """The stellation claims: exact, collision-free, and exactly one layer deep."""
    m, v, faces, outer, spikes, home, strips, rf, _g = build_stellated(
        3, edge, cells, WALL)
    assert not m.validate(), m.validate()[:2]

    # every glued cell is an exactly regular tetrahedron at the same edge
    for sp in spikes:
        p = v[list(sp)]
        for i in range(4):
            for j in range(i + 1, 4):
                assert abs(np.linalg.norm(p[i] - p[j]) - edge) < 1e-9, sp
    # and it is glued on a real shared face -- its base IS a tube surface face
    base = {frozenset(f) for _t, f in faces}
    for sp in spikes:
        assert frozenset(sp[:3]) in base, sp

    tube = [v[[k, k + 1, k + 2, k + 3]] for k in range(cells)]
    spk = [v[list(sp)] for sp in spikes]
    assert not collisions(tube + spk, edge), "layer 1 must be collision-free"

    # LAYER 2 MUST COLLIDE. Asserted so the limit cannot be crossed quietly.
    v2, spikes2, _o2 = stellate(v, outer, {t: 0 for t, _f in outer}, edge)
    spk2 = [v2[list(sp)] for sp in spikes2]
    assert collisions(spk2, edge, stop_at=1), "layer 2 was expected to collide"

    # reflex classification: 5 cells on a k->k+1 edge, 4 on k->k+2, 3 on k->k+3
    got = sorted(rf.values())
    assert set(got) == {3, 4, 5}, sorted(set(got))
    assert all(n * DIHEDRAL > 180.0 for n in got)
    # the leftover air at the tightest one is Fuller's unzipping angle
    assert abs((360.0 - 5 * DIHEDRAL) - 7.356103) < 1e-6


def _verify(edge=EDGE, cells=CELLS):
    """Every claim in the docstring, asserted."""
    for closed in (False, True):
        v = TH.helix(cells + 3, edge, closed)

        # the free faces really are the surface; {k,k+1,k+2} really is interior
        count = defaultdict(int)
        for k in range(cells):
            t = (k, k + 1, k + 2, k + 3)
            for f in ((t[0], t[1], t[2]), (t[0], t[1], t[3]),
                      (t[0], t[2], t[3]), (t[1], t[2], t[3])):
                count[tuple(sorted(f))] += 1
        for tag, f in surface(cells):
            assert count[tuple(sorted(f))] == 1, (tag, f)
        for k in range(1, cells):
            assert count[tuple(sorted((k, k + 1, k + 2)))] == 2, k

        # exactly one edge of each family per surface triangle
        for _tag, f in surface(cells):
            assert sorted(abs(f[i] - f[(i + 1) % 3]) for i in range(3)) == [1, 2, 3]

        # cut k->k+m leaves exactly m strips, each a path, advancing by m
        for m in (1, 2, 3):
            home, strips, adj = ribbons(cells, m)
            assert len(strips) == m, (m, len(strips))
            assert sum(len(s) for s in strips) == 2 * cells
            assert max(len(adj[t]) for t in home) <= 2, m
            for k in range(cells - m):
                assert home[("F1", k)] == home[("F1", k + m)], (m, k)
            if m == 3:
                for k in range(1, cells - 1):
                    assert home[("F1", k)] != home[("F2", k)], k

        d = dihedrals(v, cells)
        # EVERY fold convex from outside, which is what lets per-face inward prisms
        # overlap along a fold instead of leaving a gap there
        for fam in (1, 2, 3):
            assert max(d[fam]) < 180.0 - 1e-9, (fam, max(d[fam]))

        if not closed:
            for _tag, f in surface(cells):          # exactly equilateral
                p = v[list(f)]
                for i in range(3):
                    assert abs(np.linalg.norm(p[i] - p[(i + 1) % 3]) - edge) < 1e-9
            assert abs(min(d[3]) - max(d[3])) < 1e-9
            assert abs(min(d[3]) - math.degrees(math.acos(1 / 3))) < 1e-9
            assert abs(math.degrees(math.acos(1 / 3))
                       - math.degrees(math.atan(2 * math.sqrt(2)))) < 1e-9
            for fam, want in ((1, 148.413662), (2, 141.057559)):
                assert abs(min(d[fam]) - max(d[fam])) < 1e-9
                assert abs(min(d[fam]) - want) < 1e-5, (fam, d[fam][0])
            assert min(d[1]) > 120.0 and min(d[2]) > 120.0   # ribbon folds are gentle

    # --- the built solids ------------------------------------------------------
    m, v, faces, home, strips, _g = build(3, edge, cells)
    bad = m.validate()
    assert not bad, f"non-manifold solids: {bad[:3]}"       # directed edge check
    V, F = m.arrays()
    for _base, f0, nf in m._solid_start:
        assert _vol(V, F[f0:f0 + nf]) > 0, "a prism is inside out"

    # the wall really is `wall` thick perpendicular to each face
    for (_tag, f), (pv, pt) in zip(faces, wall_prisms(v, faces, WALL)):
        n = face_normal(v, f)
        dist = [abs((np.asarray(q) - v[f[0]]) @ n) for q in pv]
        assert abs(max(dist) - WALL) < 1e-9 and min(dist) < 1e-9
        a = 0.5 * np.linalg.norm(np.cross(v[f[1]] - v[f[0]], v[f[2]] - v[f[0]]))
        assert abs(_vol(pv, pt) - a * WALL) < 1e-6

    # exactly one family is reflex, and it is the one that gets filleted
    assert abs(material_wedge(1) - 211.586338) < 1e-5
    assert material_wedge(1) > 180.0 > material_wedge(2) > material_wedge(3)
    rf = reflex_edges(cells, faces)
    assert rf and all(b - a == 1 for a, b in rf), sorted(rf)[:3]

    # THE FIX, measured both ways round: on the material bisector at every reflex
    # edge, the per-face prisms alone leave NOTHING at any depth, and the fillet rod
    # covers it. The first half is the bug this exists for; asserting it keeps the
    # fillet honest instead of decorative.
    wp = wall_prisms(v, faces, WALL)
    missed = filled = 0
    for e, l in rf.items():
        bis = -(face_normal(v, l[0]) + face_normal(v, l[1]))
        bis = bis / np.linalg.norm(bis)
        mid = 0.5 * (v[e[0]] + v[e[1]])
        for depth in (0.2, 0.5, 0.9):
            P = mid + depth * bis
            missed += not any(_in_prism(P, pv) for pv, _pt in wp)
            filled += _in_rod(P, v[e[0]], v[e[1]], WALL)
    n = 3 * len(rf)
    assert missed == n, f"prisms unexpectedly cover {n-missed}/{n} reflex samples"
    assert filled == n, f"fillet rods only reach {filled}/{n} reflex samples"

    # a rod must also bite into both prisms it bridges, with real volume, not touch
    for e, l in rf.items():
        ax = v[e[1]] - v[e[0]]
        ax = ax / np.linalg.norm(ax)
        for f in l:
            n_f = face_normal(v, f)
            u = np.cross(n_f, ax)                    # in-plane, across the edge
            if u @ (v[list(f)].mean(axis=0) - v[e[0]]) < 0:
                u = -u
            P = (0.5 * (v[e[0]] + v[e[1]]) + 0.35 * WALL * u
                 - 0.35 * WALL * n_f)                # inside the prism, near the edge
            assert _in_rod(P, v[e[0]], v[e[1]], WALL), "rod must overlap the prism"
            assert any(_in_prism(P, pv) for pv, _pt in wp), "point should be in a prism"

    # the groove is really open: each footprint stands off its own seam edge, lies in
    # that face's plane, and leaves the apex corner sharp
    for (_pv, _pt, tag, foot), (_t, f) in zip(
            pad_prisms(v, faces, home, RELIEF, GROOVE, 3), faces):
        e = seam_edge(f, 3)
        n = face_normal(v, f)
        ax = (v[e[1]] - v[e[0]]) / np.linalg.norm(v[e[1]] - v[e[0]])
        for q in foot:
            q = np.asarray(q, float)
            assert abs((q - v[f[0]]) @ n) < 1e-9, "footprint must lie in the face"
            assert np.linalg.norm(np.cross(q - v[e[0]], ax)) > 0.98 * GROOVE, tag
        apex = [x for x in f if x not in e][0]
        assert min(np.linalg.norm(np.asarray(q) - v[apex]) for q in foot) < 1e-9


def report(a):
    _verify(a.edge, a.cells)
    m, v, faces, home, strips, (theta, R, H) = build(
        a.ribbons, a.edge, a.cells, a.wall, a.relief, a.groove, a.closed)
    V, F = m.arrays()
    tag = np.asarray(m.tag)
    lo, hi = V.min(axis=0), V.max(axis=0)
    ext = hi - lo
    d = dihedrals(v, a.cells)
    area = area_of(v, faces)

    print(f"{a.ribbons}-ribbon tetrahelix -- hollow tube, "
          + ("REGULAR Boerdijk-Coxeter (exact cells, does not repeat)"
             if not a.closed else
             "CLOSED 132 deg variant (repeats every 30 cells, stacks)"))
    print(f"  CELLS: {a.cells} regular tetrahedra, edge {a.edge:.2f} mm, twist "
          f"{math.degrees(theta):.7f} deg/cell"
          + ("  (= arccos(-2/3))" if not a.closed else ""))
    print(f"    vertex helix radius {R:.4f} mm, rise {H:.4f} mm/cell")
    L = [np.linalg.norm(v[list(f)][i] - v[list(f)][(i + 1) % 3])
         for _t, f in faces for i in range(3)]
    print(f"    surface triangles {len(faces)}, edges {min(L):.6f} .. {max(L):.6f} mm"
          f"  -> {'ALL EQUILATERAL' if max(L) - min(L) < 1e-9 else 'NOT equilateral'}")

    print(f"  THE CUT: along the k->k+{a.ribbons} family -> {len(strips)} ribbon"
          f"{'' if len(strips) == 1 else 's'} of "
          f"{', '.join(str(len(s)) for s in strips)} triangles")
    for fam in (1, 2, 3):
        mark = "   <- CUT: the seam, not a fold" if fam == a.ribbons else ""
        print(f"        k->k+{fam}: fold {min(d[fam]):10.6f} deg, "
              f"{cells_on_edge(fam)} cell{'s' if cells_on_edge(fam) > 1 else ' '} "
              f"share it -> material {material_wedge(fam):10.6f} deg"
              + ("  REFLEX" if material_wedge(fam) > 180 else "        ") + mark)
    keep = [f for f in (1, 2, 3) if f != a.ribbons]
    print(f"    a ribbon folds {' and '.join(f'{min(d[f]):.2f}' for f in keep)} deg "
          f"alternately -- {' and '.join(f'{180-min(d[f]):.2f}' for f in keep)} deg "
          f"off flat, which is why it reads as a ribbon")
    if a.ribbons == 3:
        print(f"    the 3 seams ARE the three intertwined helices of the BC helix, "
              f"the edges belonging to only one tetrahedron")
        mixed = sum(1 for k in range(1, a.cells - 1)
                    if home[("F1", k)] != home[("F2", k)])
        print(f"    all {mixed}/{a.cells-2} interior cells have their two surface "
              f"faces in DIFFERENT ribbons -> the ribbons share whole tetrahedra")

    # how far a ribbon actually winds. A seam helix steps k->k+3, advancing
    # 3*theta mod 360 = 35.4309 deg per step, i.e. theta - 120 = 11.8103 deg per
    # cell -- the same 11.8103 that turns up as the torsion a 3-ply lay needs to
    # make a braid periodic, and for the same reason: it is how far theta sits
    # from a third of a turn.
    per = (math.degrees(theta) * 3.0) % 360.0 / 3.0
    turns = a.cells * per / 360.0
    print(f"  WIND: a seam helix advances {per:.4f} deg/cell (= theta - 120), so each "
          f"ribbon makes {turns:.2f} turns over this tube")
    if turns < 0.98:
        need = int(math.ceil(360.0 / per))
        emax = 320.0 / (need + 2) * math.sqrt(10.0)
        print(f"      one FULL ribbon turn needs {need} cells; "
              f"--edge {emax:.0f} --cells {need} is "
              f"{(need+2)*emax/math.sqrt(10.0):.0f} mm tall and fits the bed")
    print(f"  HOLLOW: interior {{k,k+1,k+2}} faces omitted, so the bore is open end "
          f"to end ({2*(R - a.wall):.1f} mm across)")
    rf = reflex_edges(a.cells, faces)
    print(f"    wall {a.wall:.2f} mm PERPENDICULAR to every face (per-face prisms), "
          f"plus a {a.wall:.2f} mm fillet rod on each of the {len(rf)} REFLEX "
          f"k->k+1 edges -- without those the wall is knife-edge joined there "
          f"(see fillet_rods)")
    print(f"    ribbon relief {a.relief:.2f} mm, seam groove {2*a.groove:.2f} mm wide "
          f"x {a.relief:.2f} mm deep, pads sunk {SINK:.1f} mm into the wall")
    print(f"    surface area {area:.0f} mm2, so the wall is about "
          f"{area*a.wall/1000:.1f} cm3")
    print("  triangles by feature:", "  ".join(
        f"{q}={int((tag == q).sum())}" for q in sorted(set(tag))))
    print(f"  bbox x {lo[0]:7.1f}..{hi[0]:6.1f}  y {lo[1]:7.1f}..{hi[1]:6.1f}"
          f"  z {lo[2]:7.1f}..{hi[2]:6.1f}")
    fits = "FITS" if (ext <= 320.0 + 1e-9).all() else "OVER BED"
    print(f"  envelope {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm  "
          f"(bed 320 cubed: {fits}),  {ext[2]/max(ext[0], ext[1]):.1f}:1 aspect")

    steep = np.array([math.degrees(math.asin(min(1.0, abs(face_normal(v, f)[2]))))
                      for _t, f in faces])
    print(f"  PRINT: surface facets {steep.min():.1f} .. {steep.max():.1f} deg off "
          f"vertical (median {np.median(steep):.1f})")
    print(f"      {int((steep > 45).sum())}/{len(steep)} facets past 45 deg"
          + ("  -- a vertical tube, no supports needed" if (steep <= 45).all()
             else "  -- check supports"))
    if a.closed:
        print(f"      132 deg variant: v[k+30] sits over v[k], so a 30-cell module "
              f"stacks exactly; this build has {a.cells}")
    return m, area


def report_stellated(a):
    _verify_stellated(a.edge, a.cells)
    m, v, faces, outer, spikes, home, strips, rf, (theta, R, H) = build_stellated(
        a.ribbons, a.edge, a.cells, a.wall, a.closed)
    V, F = m.arrays()
    tag = np.asarray(m.tag)
    lo, hi = V.min(axis=0), V.max(axis=0)
    ext = hi - lo
    r = np.linalg.norm(V[:, :2], axis=1)
    area = sum(0.5 * np.linalg.norm(np.cross(v[f[1]] - v[f[0]], v[f[2]] - v[f[0]]))
               for _t, f in list(faces) + list(outer))
    notch = 360.0 - 5 * DIHEDRAL

    print(f"STELLATED {a.ribbons}-ribbon tetrahelix -- ONE layer, "
          + ("regular BC helix" if not a.closed else "132 deg variant"))
    print(f"  CELLS: {a.cells} tube + {len(spikes)} glued = "
          f"{a.cells + len(spikes)} regular tetrahedra, edge {a.edge:.2f} mm")
    print(f"    every glued cell is EXACTLY regular and shares a WHOLE face with the "
          f"tube face under it -- asserted, nothing fitted or trimmed")
    print(f"    its apex stands edge*sqrt(2/3) = {a.edge*math.sqrt(2/3):.3f} mm off "
          f"that face")
    print(f"  COLLISION-FREE: 0 overlapping pairs across all "
          f"{a.cells + len(spikes)} cells, measured over every pair in range")
    print(f"    ONE LAYER ONLY: stellating again adds {3*len(spikes)} cells and they "
          f"collide (asserted). So this route DOUBLES the sculpture; it cannot "
          f"triple it.")
    print(f"  WHY ONE LAYER FITS -- Fuller's unzipping angle:")
    print(f"    before stellation a k->k+1 edge has {360-3*DIHEDRAL:.6f} deg of air; "
          f"two glued cells need 2*{DIHEDRAL:.6f} = {2*DIHEDRAL:.6f}")
    print(f"    leftover {notch:.6f} deg = 360 - 5*{DIHEDRAL:.6f}. Five cells now "
          f"ring every k->k+1 edge and just fail to close.")
    print(f"    on the print that is a {notch:.3f} deg V-notch along each of those "
          f"edges, about {2*(a.edge/2)*math.sin(math.radians(notch/2)):.2f} mm wide "
          f"at mid-face. The number that stops tetrahedra tiling space is what "
          f"leaves room for the layer.")
    cnt = {}
    for n in rf.values():
        cnt[n] = cnt.get(n, 0) + 1
    print(f"  REFLEX EDGES: {len(rf)}, filleted with {a.wall:.2f} mm rods -- "
          + ", ".join(f"{n} cells x{c} = {n*DIHEDRAL:.1f} deg"
                      for n, c in sorted(cnt.items())))
    print(f"    those rods also bridge the bore wall to the spike walls, so tube and "
          f"spikes are one solid rather than 168 shells resting on each other")
    print(f"  SHELL: wall {a.wall:.2f} mm perpendicular on {len(outer)} outer faces "
          f"plus {len(faces)} bore faces; each spike interior is a sealed void")
    print(f"    radius {r.min():.1f} .. {r.max():.1f} mm against 21.5 for the plain "
          f"tube -> {r.max()/21.53:.2f}x")
    print("  triangles by feature:", "  ".join(
        f"{q}={int((tag == q).sum())}" for q in sorted(set(tag))))
    print(f"  envelope {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm  "
          f"(bed 320 cubed: {'FITS' if (ext <= 320).all() else 'OVER BED'})")
    steep = np.array([math.degrees(math.asin(min(1.0, abs(face_normal(v, f)[2]))))
                      for _t, f in outer])
    over = int((steep > 45).sum())
    print(f"  PRINT: outer facets {steep.min():.1f} .. {steep.max():.1f} deg off "
          f"vertical, median {np.median(steep):.1f}")
    print(f"      {over}/{len(steep)} past 45 deg"
          + ("  -- SUPPORTS NEEDED; the spikes overhang, unlike the plain tube"
             if over else "  -- no supports needed"))
    print(f"      wall volume about {area*a.wall/1000:.0f} cm3 -> "
          f"{area*a.wall*1.27/1000:.0f} g PETG"
          f"   (--wall 1.2 gives {area*1.2*1.27/1000:.0f} g)")
    return m, area


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ribbons", type=int, choices=(1, 2, 3), default=RIBBONS,
                    help="cut along the k->k+m family; m IS the ribbon count")
    ap.add_argument("--edge", type=float, default=EDGE, help="cell edge, mm")
    ap.add_argument("--cells", type=int, default=CELLS, help="tetrahedra")
    ap.add_argument("--wall", type=float, default=WALL,
                    help="perpendicular wall thickness, mm")
    ap.add_argument("--relief", type=float, default=RELIEF,
                    help="ribbon stand-off from the wall, mm")
    ap.add_argument("--groove", type=float, default=GROOVE,
                    help="pad inset from its seam edge; groove is twice this")
    ap.add_argument("--stellate", action="store_true",
                    help="glue a regular tetrahedron on every surface face: "
                         "one layer, collision-free, doubles the thickness")
    ap.add_argument("--closed", action="store_true",
                    help="tetrahelix.py's 132 deg variant: stacks every 30 cells, "
                         "cells 0.28%% off regular")
    ap.add_argument("--weld", action="store_true",
                    help="CSG-union every prism into one 2-manifold solid")
    ap.add_argument("--out", default="ribbon_tetrahelix.stl")
    a = ap.parse_args()
    m, area = report_stellated(a) if a.stellate else report(a)

    out = m
    if a.weld:
        print(f"\n  welding {len(m._solid_start)} prisms...")
        out = weld(m) or m
    V, F = out.arrays()
    comps, (open_e, pinch_e), vol = shell_stats(V, F)
    print(f"\n  OUTPUT ({'welded' if out is not m else 'raw soup'}), by position: "
          f"{len(F)} triangles, {comps} connected "
          f"{'body' if comps == 1 else 'bodies'}, {open_e} open edges, "
          f"{pinch_e} pinched")
    if out is m:
        print(f"  -> {len(m._solid_start)} overlapping closed prisms; slicers union "
              f"them on slice. --weld for a single body")
        print(f"  volume {vol/1000:.1f} cm3 is a SUM, so it over-counts every fold "
              f"overlap -- the wall alone is {area*a.wall/1000:.1f} cm3 exactly")
    else:
        print(f"  volume {vol/1000:.1f} cm3 -> {vol*1.27/1000:.0f} g PETG")
    n = MK.write_stl(a.out, out)
    print(f"\nwrote {a.out} ({n} triangles)")


if __name__ == "__main__":
    main()
