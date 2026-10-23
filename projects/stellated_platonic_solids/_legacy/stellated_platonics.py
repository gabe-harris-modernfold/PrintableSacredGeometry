#!/usr/bin/env python3
"""First stellation of the Platonic solids -- printable STLs.

What "one level of stellation" means here
-----------------------------------------
Stellation extends a polyhedron's *face planes* until they close off new cells,
then keeps the core plus the first shell of cells. For a Platonic solid every
face is equivalent, so the first shell is one cell per face, and that cell is the
pyramid standing on face F bounded by the planes of F's edge-neighbours. Its base
is exactly F (plane(F) meets plane(G) along their shared edge), so the apex sits
on F's normal axis where all the neighbour planes cross it:

    apex = c_F + h * n_F ,   h = (d_G - c_F . n_G) / (n_F . n_G)

with plane(G) = {x : x . n_G = d_G}. The same h comes out for every neighbour G --
the script asserts that rather than assuming it.

That single formula also decides *whether* a stellation exists. Writing D for the
dihedral angle, n_F . n_G = -cos(D), so

    D > 90 deg  ->  h > 0    a real cell            octahedron, dodecahedron, icosahedron
    D = 90 deg  ->  h = inf  planes never meet      cube            (D = 90.000)
    D < 90 deg  ->  h < 0    the cell is behind F   tetrahedron     (D = 70.529)

So the tetrahedron and the cube have **no stellations at all** -- a standard result
(4 planes in general position bound one cell; the cube's 6 planes are 3 parallel
pairs and bound one cell) that this script rederives numerically instead of
hard-coding. The three that do stellate give the three classical solids:

    octahedron   -> stella octangula            (compound of two tetrahedra)
    dodecahedron -> small stellated dodecahedron  (Kepler)
    icosahedron  -> small triambic icosahedron    (first stellation of the icosahedron)

For the two degenerate cases the script still emits the conventional spiky
stand-in so the printed set is complete: the **dual compound**, the union of the
solid with its midsphere dual. That is a compound, not a stellation, and is
labelled as such. It gives the stella octangula for the tetrahedron (the "Merkaba"
that gets printed under the name "stellated tetrahedron" -- note this is congruent
to the octahedron's stellation, because the stella octangula is both) and the
compound of cube and octahedron for the cube.

Mesh
----
Every original face is completely covered by its own pyramid base, so the surface
of the result is nothing but the pyramid flanks: one triangle fan per face, apex to
base edge. That is closed and 2-manifold by construction -- no boolean union is
needed, which matters because this box has no CSG backend for trimesh.

Outputs (150 mm max extent, seated on the widest face of the convex hull)
-------------------------------------------------------------------------
stellated_octahedron_150mm.stl        stella octangula                TRUE stellation
stellated_dodecahedron_150mm.stl      small stellated dodecahedron    TRUE stellation
stellated_icosahedron_150mm.stl       small triambic icosahedron      TRUE stellation
stellated_tetrahedron_150mm.stl       stella octangula                dual compound
stellated_cube_150mm.stl              cube + octahedron compound      dual compound
stellated_platonics_plate.stl         all five at 90 mm on one 320 mm bed

Run:  python stellated_platonics.py
"""

import math

import numpy as np
import trimesh
from scipy.spatial import ConvexHull

# ---------------------------------------------------------------- parameters

SIZE = 150.0          # max bounding-box dimension of each single-solid STL, mm
WALL = 4.0            # shell wall thickness, mm -- measured off this repo's own
                      # Platonic solids: 09..13_*_2d.stl are line-work 4.00 mm wide
                      # by 4.00 mm tall, the same on all five plates
PLATE_SIZE = 90.0     # per-solid size on the combined print plate, mm
PLATE_GAP = 12.0      # clearance between solids on the plate, mm
BED = 320.0           # printer bed, mm (320 x 320 x 320, PETG)

PHI = (1.0 + 5.0 ** 0.5) / 2.0
PLANE_KEY = 7         # decimals used to group coplanar hull facets

# ------------------------------------------------------------------- solids


def platonic(name):
    """Vertices of the five Platonic solids, centred on the origin."""
    if name == "tetrahedron":
        return np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], float)
    if name == "cube":
        return np.array([(x, y, z) for x in (-1, 1) for y in (-1, 1)
                         for z in (-1, 1)], float)
    if name == "octahedron":
        return np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                         (0, -1, 0), (0, 0, 1), (0, 0, -1)], float)
    if name == "dodecahedron":
        a, b = 1.0 / PHI, PHI
        v = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
        v += [(0, s * a, t * b) for s in (-1, 1) for t in (-1, 1)]
        v += [(s * a, t * b, 0) for s in (-1, 1) for t in (-1, 1)]
        v += [(s * b, 0, t * a) for s in (-1, 1) for t in (-1, 1)]
        return np.array(v, float)
    if name == "icosahedron":
        v = []
        for s in (-1, 1):
            for t in (-1, 1):
                v += [(0, s, t * PHI), (s, t * PHI, 0), (s * PHI, 0, t)]
        return np.array(v, float)
    raise ValueError(name)


def polygonal_faces(V):
    """Convex hull of V as polygons, not triangles.

    Returns [(loop, n, d)] with `loop` the vertex indices ordered counter-clockwise
    seen from outside and the face plane {x : x . n = d}, n outward unit.
    """
    hull = ConvexHull(V)
    groups = {}
    for eq, simplex in zip(hull.equations, hull.simplices):
        key = tuple(np.round(eq, PLANE_KEY))
        groups.setdefault(key, set()).update(int(i) for i in simplex)

    faces = []
    for key, idxs in groups.items():
        idxs = sorted(idxs)
        pts = V[idxs]
        c = pts.mean(axis=0)
        # Re-fit the plane from all of the face's vertices. Qhull's per-simplex
        # equation is only good to ~1e-8 on a 5-gon, and the apex solve below
        # divides by (n_F . n_G) -- it wants the normal to full precision.
        n = np.linalg.svd(pts - c)[2][-1]
        if n @ np.array(key[:3], float) < 0:
            n = -n
        d = float(n @ c)
        u = pts[0] - c
        u -= (u @ n) * n
        u /= np.linalg.norm(u)
        w = np.cross(n, u)                      # (u, w, n) right-handed
        ang = np.arctan2((pts - c) @ w, (pts - c) @ u)
        loop = [idxs[i] for i in np.argsort(ang)]
        faces.append((loop, n, d))
    return faces


def edge_neighbours(faces):
    """face index -> list of face indices sharing an edge."""
    owner = {}
    for fi, (loop, _n, _d) in enumerate(faces):
        for a, b in zip(loop, loop[1:] + loop[:1]):
            owner.setdefault(frozenset((a, b)), []).append(fi)
    nb = [set() for _ in faces]
    for pair in owner.values():
        assert len(pair) == 2, "non-manifold edge in the base solid"
        nb[pair[0]].add(pair[1])
        nb[pair[1]].add(pair[0])
    return [sorted(s) for s in nb]


# --------------------------------------------------------------- stellation


def dihedral_deg(faces, nb):
    a, b = faces[0][1], faces[nb[0][0]][1]
    return math.degrees(math.pi - math.acos(float(np.clip(a @ b, -1, 1))))


def stellation_apexes(V, faces, nb):
    """Apex per face from the neighbouring face planes, or (None, reason)."""
    apexes = np.zeros((len(faces), 3))
    for fi, (loop, n, _d) in enumerate(faces):
        c = V[loop].mean(axis=0)
        hs = []
        for fj in nb[fi]:
            _lj, nj, dj = faces[fj]
            den = float(n @ nj)
            if abs(den) < 1e-9:
                return None, "adjacent face planes are parallel to the face normal"
            hs.append((dj - float(c @ nj)) / den)
        spread = max(hs) - min(hs)
        assert spread < 1e-9 * max(1.0, abs(hs[0])), "neighbour planes disagree"
        h = float(np.mean(hs))
        if h <= 1e-9:
            return None, "the cell closes behind the face, not in front of it"
        apexes[fi] = c + h * n
    return apexes, None


def dual_compound(V, faces):
    """Fallback for the two solids that do not stellate: union with the dual.

    P and its dual P* share a midsphere, so P and P* cut each other along the
    rectified solid R = hull(edge midpoints of P) -- octahedron for the tetrahedron,
    cuboctahedron for the cube. Everything of P or P* sticking out past R is one
    corner, and that corner is the pyramid on the face of R that faces it, with its
    apex at the vertex of P or P* extreme in that direction. So

        P union P*  =  Kleetope(R)

    with no boolean needed. Tetrahedron -> stella octangula (the Merkaba);
    cube -> compound of cube and octahedron.

    Returns (R, faces_of_R, apexes) ready for spiked_mesh().
    """
    seen, mids = set(), []
    for loop, _n, _d in faces:
        for a, b in zip(loop, loop[1:] + loop[:1]):
            key = frozenset((a, b))
            if key not in seen:
                seen.add(key)
                mids.append(0.5 * (V[a] + V[b]))
    R = np.array(mids)

    rho = float(np.linalg.norm(R[0]))                 # midradius
    inr = faces[0][2]                                 # inradius of P
    dual_verts = np.array([(rho * rho / inr) * n for _loop, n, _d in faces])
    tips = np.vstack([V, dual_verts])

    rfaces = polygonal_faces(R)
    apexes = np.array([tips[int(np.argmax(tips @ n))] for _loop, n, _d in rfaces])
    return R, rfaces, apexes


def spiked_mesh(V, faces, apexes):
    """Surface = pyramid flanks only; the original faces are buried under them."""
    verts = np.vstack([V, apexes])
    base = len(V)
    tris = []
    for fi, (loop, _n, _d) in enumerate(faces):
        a = base + fi
        for i, j in zip(loop, loop[1:] + loop[:1]):
            tris.append((i, j, a))
    return verts, np.array(tris, np.int64)


def check_planes_agree(V, faces, apexes, nb):
    """Every pyramid flank must lie in a face plane of the ORIGINAL solid.

    That is the defining property of a stellation and the reason the flanks fuse
    into the big star faces (pentagrams, triambi) of the classical solids.
    """
    worst = 0.0
    for fi, (loop, _n, _d) in enumerate(faces):
        for i, j in zip(loop, loop[1:] + loop[:1]):
            tri = np.array([V[i], V[j], apexes[fi]])
            best = min(float(np.abs(tri @ faces[fj][1] - faces[fj][2]).max())
                       for fj in nb[fi])
            worst = max(worst, best)
    return worst


# ------------------------------------------------------------ orient / scale


def seat(verts):
    """Rotate the widest convex-hull face down, scale later, drop onto z = 0.

    For a star solid the hull is the hull of the spike tips, so the widest hull
    face is the most stable tripod/pentapod of tips to print on.
    """
    hull_faces = polygonal_faces(verts)
    best, best_area, best_n = None, -1.0, 0
    for loop, n, _d in hull_faces:
        pts = verts[loop]
        c = pts.mean(axis=0)
        area = 0.5 * float(np.linalg.norm(
            np.cross(pts - c, np.roll(pts - c, -1, axis=0)).sum(axis=0)))
        if area > best_area:
            best, best_area, best_n = n, area, len(loop)

    down = -np.asarray(best, float)
    target = np.array([0.0, 0.0, -1.0])
    v = np.cross(down, target)
    s, c = float(np.linalg.norm(v)), float(down @ target)
    if s < 1e-12:
        R = np.eye(3) if c > 0 else -np.eye(3)
    else:
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + K + K @ K * ((1 - c) / (s * s))
    out = verts @ R.T
    return out, R, best_area, best_n


def to_size(verts, size):
    """Scale to `size` across, still centred on the origin (hollow() needs that)."""
    ext = verts.max(axis=0) - verts.min(axis=0)
    return verts * (size / float(ext.max()))


def drop(verts):
    v = verts.copy()
    v[:, :2] -= 0.5 * (v[:, :2].min(axis=0) + v[:, :2].max(axis=0))
    v[:, 2] -= v[:, 2].min()
    return v


# ----------------------------------------------------------------- shelling


def hollow(verts, tris, wall):
    """Shell to a perpendicular wall thickness of `wall` by scaling, not offsetting.

    Every flank of a Platonic stellation lies in a face plane of the core solid, and
    those planes are all the same distance d from the centre. Pushing all of them in
    by t therefore lands on the plane set of the SAME stellation at inradius d - t --
    that is, on a uniform scale by k = 1 - t/d. So the inner surface is the outer one
    scaled about the centre, and the wall is exactly t on every flank, with no offset
    machinery and no boolean.

    The dual compounds are the one case with two plane distances (the cube's faces and
    the octahedron's). Scaling keys off the nearer set, so the thinnest wall is still
    exactly t and the other family comes out proportionally thicker.

    Nesting the reversed inner surface inside the outer one gives a closed shell around
    a sealed void -- two components, which is what a hollow part is supposed to be.
    """
    tri = verts[tris]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    d = np.einsum("ij,ij->i", n, tri[:, 0])          # every plane, origin to plane
    dmin = float(d.min())
    if wall >= dmin:
        raise ValueError(f"wall {wall} mm exceeds the {dmin:.2f} mm inradius")
    k = 1.0 - wall / dmin
    v = np.vstack([verts, verts * k])
    f = np.vstack([tris, tris[:, ::-1] + len(verts)])
    walls = d * (1.0 - k)
    return v, f, {"k": k, "inradius": dmin,
                  "wall_min": float(walls.min()), "wall_max": float(walls.max())}


# ------------------------------------------------------------------- driver

SOLIDS = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"]

TRUE_NAMES = {
    "octahedron": "stella octangula (compound of two tetrahedra)",
    "dodecahedron": "small stellated dodecahedron",
    "icosahedron": "small triambic icosahedron",
}
FALLBACK_NAMES = {
    "tetrahedron": "stella octangula = tetrahedron + dual (Merkaba)",
    "cube": "compound of cube and octahedron",
}
FALLBACK_BASE = {"tetrahedron": "octahedron", "cube": "cuboctahedron"}


def build(name):
    V = platonic(name)
    faces = polygonal_faces(V)
    nb = edge_neighbours(faces)
    D = dihedral_deg(faces, nb)
    inradius = faces[0][2]

    apexes, why = stellation_apexes(V, faces, nb)
    true_stellation = apexes is not None
    if true_stellation:
        planar_err = check_planes_agree(V, faces, apexes, nb)
        label = TRUE_NAMES[name]
        base_name, base_V, base_faces = name, V, faces
    else:
        base_V, base_faces, apexes = dual_compound(V, faces)
        planar_err = float("nan")
        label = FALLBACK_NAMES[name]
        base_name = FALLBACK_BASE[name]

    verts, tris = spiked_mesh(base_V, base_faces, apexes)
    ratio = float(np.linalg.norm(apexes, axis=1).max()) / inradius

    verts, R, contact_area, contact_pts = seat(verts)
    sides = sorted({len(loop) for loop, _n, _d in base_faces})
    return {
        "name": name, "label": label, "true": true_stellation, "why": why,
        "dihedral": D, "n_faces": len(base_faces), "base_name": base_name,
        "sides": "+".join(str(k) for k in sides),
        "ratio": ratio, "planar_err": planar_err,
        "verts": verts, "tris": tris, "rot": R,
        "core_V": V @ R.T, "core_faces": faces,
        "contact_area": contact_area, "contact_pts": contact_pts,
    }


def scale_factor(rec, size):
    ext = rec["verts"].max(axis=0) - rec["verts"].min(axis=0)
    return size / float(ext.max())


def as_mesh(rec, size, wall=WALL):
    v = to_size(rec["verts"], size)
    if wall and wall > 0.0:
        v, f, info = hollow(v, rec["tris"], wall)
    else:
        f, info = rec["tris"], None
    return trimesh.Trimesh(vertices=drop(v), faces=f, process=False), info


def core_triangles(rec):
    """Fan-triangulated core solid, in the same frame as rec["verts"]."""
    tris = []
    for loop, _n, _d in rec["core_faces"]:
        for k in range(1, len(loop) - 1):
            tris.append((loop[0], loop[k], loop[k + 1]))
    return rec["core_V"], np.array(tris, np.int64)


def main():
    recs = [build(n) for n in SOLIDS]

    print("=" * 96)
    print("FIRST STELLATION OF THE PLATONIC SOLIDS")
    print("=" * 96)
    print(f"{'solid':<13}{'dihedral':>10}{'stellates':>11}  {'result':<44}{'tip/in':>7}")
    print("-" * 96)
    for r in recs:
        print(f"{r['name']:<13}{r['dihedral']:>9.3f}d{('yes' if r['true'] else 'NO'):>11}"
              f"  {r['label']:<44}{r['ratio']:>7.4f}")
    print("-" * 96)
    for r in recs:
        if not r["true"]:
            print(f"  {r['name']}: no stellation exists -- {r['why']} "
                  f"(dihedral {r['dihedral']:.3f} deg <= 90). Emitting the dual-compound "
                  f"cumulation instead.")
    print()

    print(f"Shelled to WALL = {WALL:.1f} mm, measured off this repo's own Platonic")
    print(f"solids (09..13_*_2d.stl are 4.00 mm line-work on all five plates).")
    print()

    for r in recs:
        m, sh = as_mesh(r, SIZE)
        solid, _ = as_mesh(r, SIZE, wall=0.0)
        out = f"stellated_{r['name']}_{int(SIZE)}mm.stl"
        m.export(out)
        ext = m.extents
        euler = len(m.vertices) - len(m.edges_unique) + len(m.faces)
        fits = all(e <= BED + 1e-6 for e in ext)
        print(f"{out}")
        print(f"    {r['n_faces']} pyramids on the {r['sides']}-gon faces of the "
              f"{r['base_name']} -> {len(m.faces)} triangles, {len(m.vertices)} "
              f"vertices, Euler {euler}")
        print(f"    extents {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm  "
              f"{'FITS' if fits else 'TOO BIG'} the {BED:.0f} mm bed")
        print(f"    wall {sh['wall_min']:.2f} - {sh['wall_max']:.2f} mm  "
              f"(inner shell = outer x {sh['k']:.4f}, core inradius {sh['inradius']:.1f} mm)")
        print(f"    watertight {m.is_watertight}   winding consistent {m.is_winding_consistent}"
              f"   volume {m.volume / 1000.0:.1f} cm^3 ({'outward' if m.volume > 0 else 'INVERTED'})")
        print(f"    hollowing saves {100.0 * (1.0 - m.volume / solid.volume):.0f}% "
              f"of the {solid.volume / 1000.0:.0f} cm^3 solid "
              f"({1.27 * m.volume / 1000.0:.0f} g of PETG at 100% infill)")
        print(f"    seats on {r['contact_pts']} spike tips")
        if r["true"]:
            print(f"    flanks lie in the original face planes to {r['planar_err']:.2e} "
                  f"(exact stellation)")
        print()

    # ------------------------------------------------ one plate, all five
    tiles = []
    for r in recs:
        tiles.append(as_mesh(r, PLATE_SIZE)[0])
    cols = [3, 2]
    pitch = PLATE_SIZE + PLATE_GAP
    placed, k = [], 0
    rows = len(cols)
    for ri, ncol in enumerate(cols):
        y = (ri - (rows - 1) / 2.0) * pitch
        for ci in range(ncol):
            x = (ci - (ncol - 1) / 2.0) * pitch
            t = tiles[k].copy()
            t.apply_translation([x, y, 0.0])
            placed.append(t)
            k += 1
    plate = trimesh.util.concatenate(placed)
    plate.export("stellated_platonics_plate.stl")
    pe = plate.extents
    print(f"stellated_platonics_plate.stl")
    print(f"    all five at {PLATE_SIZE:.0f} mm, {len(plate.faces)} triangles")
    print(f"    footprint {pe[0]:.1f} x {pe[1]:.1f} x {pe[2]:.1f} mm  "
          f"{'FITS' if pe[0] <= BED and pe[1] <= BED and pe[2] <= BED else 'TOO BIG'} "
          f"the {BED:.0f} mm bed")
    print(f"    watertight {plate.is_watertight}  "
          f"{plate.body_count} surfaces (5 objects x outer+inner)  "
          f"volume {plate.volume / 1000.0:.1f} cm^3")
    print()
    print("Print notes")
    print("  - Each file is a closed shell around a sealed void, so it reads as two")
    print("    surfaces per object. That is what a hollow part looks like; slicers")
    print("    fill the wall and leave the cavity empty. Set infill to 0%: the wall is")
    print("    already the 4 mm you asked for, and infill would only fill the tips.")
    print("  - Every spike overhangs, so slice with supports on. The seated orientation")
    print("    puts the widest hull facet on the bed, the flattest rest each solid has.")


if __name__ == "__main__":
    main()
