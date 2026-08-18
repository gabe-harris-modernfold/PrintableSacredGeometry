#!/usr/bin/env python3
"""E8 horn torus lattice -- the 240 roots of E8 as a printable 150 mm strut lattice.

The construction
----------------
The Coxeter element c of E8 (the product of its eight simple reflections) has order 30, and
it acts on R^8 by rotating four mutually orthogonal 2-planes at once, through 2*pi*m/30 for
the eight exponents m = 1, 7, 11, 13, 17, 19, 23, 29 (conjugate in pairs, so four planes).
It permutes the 240 roots in exactly 8 orbits of 30.

Project a root alpha onto the m = 1 plane and onto the m = 7 plane. Writing those two
projections as complex numbers z1, z7, applying c multiplies them by exp(2*pi*i/30) and
exp(14*pi*i/30). So |z1| is constant along an orbit while arg(z1) and arg(z7) each advance
by a fixed step -- which is to say the pair (arg z1, arg z7) maps the root system onto a
2-torus, and each Coxeter orbit becomes a closed coil on it.

That is the object built here:

  toroidal angle theta = arg(z1)     the Coxeter (Petrie) plane angle, step 2*pi/30
  poloidal angle phi   = arg(z7)     step 7 * 2*pi/30
  tube radius          = |z1|        the Petrie ring radius, 8 distinct values

Each of the 8 orbits therefore becomes a curve that runs once around the torus while winding
7 times around the tube -- a (1, 7) coil -- and the 8 coils sit at 8 nested radii inside the
tube. The cross-section of the tube is the eight-ring Petrie diagram, which is what makes
this a torus and not just a ring of struts.

Which torus: a horn torus, tube radius equal to major radius, so the hole closes to a single
point on the axis rather than staying open. TUBE_FRAC picks this out of the family and horn is
not an arbitrary choice -- see the note on that parameter. One consequence is structural and
worth expecting in the print: with a = R the outer shell's inner equator collapses onto the
axis, so one root sits exactly at the origin and the outer coil passes through that same point
once per poloidal turn. Seven strands and a node meet there in a 7-fold hub, which is the
densest and strongest point in the object.

The eight radii come out as four exact golden-ratio pairs:

  0.2091 : 0.3383    0.4158 : 0.6728    0.5028 : 0.8135    0.6180 : 1.0000

each right-hand value phi = 1.618... times its partner. That is the E8 -> H4 folding showing
up in the print: E8's roots project to two 600-cells, one phi times the size of the other,
so the shells of this torus pair off the same way. The script asserts it rather than trusting
it. (0.6180 : 1.0000 is 1/phi : 1, and five of the radii are 2*sin(k*6 deg), k = 1..5.)

Struts
------
coils   the 8 orbits, swept as continuous closed tubes through their own 30 nodes. The frame
        used for the sweep is the poloidal radial direction, which is exactly perpendicular
        to the coil tangent (asserted) and 2*pi-periodic, so the tube closes on itself with
        no accumulated twist and no seam.
braces  root pairs at 60 degrees -- inner product 1, meaning alpha - beta is also a root --
        which is the real adjacency of the E8 root system. See "How much of E8 is here".
nodes   240 spheres, one per root, sitting on the coil sample points they fall on.

How much of E8 is here
----------------------
All 240 roots, and the placement is faithful: the map alpha -> (arg z1, arg z7, |z1|) is
injective on them, so no two roots collide and nothing about their image in the
(m=1) + (m=7) four-space is discarded. |z7| is not thrown away either -- it is constant on
each ring and the |z7| radii are the same eight numbers as the |z1| radii, permuted, so it
carries no information the shell index does not already have. What a 3D object must lose is
the m = 11 and m = 13 planes.

The edges are where a choice has to be made. Every root has 56 neighbours at 60 degrees;
that is 6720 edges, 335 m of strut, about 2.7 kg of PETG at any sane diameter, and it would
print as a solid lump rather than a lattice. So a subset is unavoidable -- but a length
cutoff is the wrong way to pick it, because the Coxeter element does not act on this
embedding as a rigid motion, so a metric rule cuts across the symmetry and leaves some roots
with no braces at all.

Instead: c permutes the 6720 edges into exactly 224 orbits of 30, and whole orbits are kept,
shortest first, until every root carries MIN_BRACE_DEGREE of them. That subset is
Coxeter-equivariant, and it forces the brace degree to be constant on each of the eight root
orbits -- an edge orbit gives every root of the orbits it touches the same number of ends. At
the default degree 12 that comes to 63 orbits, 1890 braces, 28.1% of E8's edge set, degree 12
to 18 per shell. MIN_BRACE_DEGREE is the density knob.

Targeting the degree rather than a fixed orbit count matters because the ordering of the
orbits by length depends on the embedding: the 58 orbits that give every root 12 braces on a
ring torus leave a shell with only 9 on this one.

Why the object has no rotational symmetry: c acts here as 12 degrees about the z axis plus
84 degrees poloidally, and a poloidal rotation is not a rigid motion of R^3. It would become
one if some power of c had zero poloidal shift, i.e. if 7k = 0 mod 30 for some k < 30 -- but
E8's exponents are exactly the units mod 30, so no power of c works and the honest embedding
is forced to be asymmetric. The 30-fold order is visible in the angles, not in the solid.

Printing
--------
Torus axis vertical, as generated: 147 x 150 x 79 mm, ~199 g of PETG, standing on the seven
poloidal low points of the outermost coil -- seven contact patches spread round a 76 mm
circle, so it is stable on the plate but wants a brim. 37% of the 52 m of centreline lies
within 20 degrees of horizontal, with shallow runs up to 48 mm, so this needs supports;
organic/tree supports suit a lattice this open far better than grid. Braces cross and fuse
where they cross, which is what a space frame does and is where its stiffness comes from, but
no two coils come closer than 3.5 mm, so the eight shells stay legible as eight shells.

Do not recentre the part in x-y. The roots are not distributed symmetrically about the axis,
so centring on the bounding box slides the whole lattice ~3.5 mm off its own axis and the
throat comes out eccentric. build() leaves the axis on x = y = 0 and asserts the radial
envelope matches the analytic one in both directions.

The file is not one closed surface. It is 2138 separately closed, mutually overlapping solids
-- which is what a slicer unions, and what trimesh reports as watertight with 2138 bodies.
There is no mesh-boolean backend in this environment, so that is the deliberate approach here
rather than a shortcut around a failed union.

Outputs
-------
e8_horn_torus_lattice_150mm.stl   one fused body, 150 mm across; the name states the class
                                  measured from MAJOR vs TUBE, so it tracks TUBE_FRAC
e8_torus_lattice.png              preview sheet (preview_e8_torus.py)
e8_torus_lattice_views.png        object gallery (preview_e8_torus_gallery.py)

Run:  python e8_torus_lattice.py
"""

import itertools
import math
import os

import numpy as np

from mesh_kit import Mesh, tube, write_stl

# ---------------------------------------------------------------- parameters

SIZE = 150.0        # bounding-box width of the finished lattice, mm
MAJOR = 37.76       # torus major radius, mm -- picked so normalising to SIZE is a no-op
                    # and the strut radii below are therefore literal
#: Tube radius as a fraction of MAJOR, which is what picks the torus out of the family:
#: below 1 a ring torus with an open hole, exactly 1 a horn torus whose hole closes to a
#: single point on the axis, above 1 a spindle. Horn is the default because it is where the
#: eight shells spread furthest apart -- the roots sit 6.7 mm apart at a/R = 1 against
#: 5.2 mm at 0.63 -- and because at a/R >= ~1.2 the shells crowd the axis and the roots
#: start to merge into each other (2.1 mm apart by 1.4), which no longer prints as a lattice.
TUBE_FRAC = 1.0
#: Braces are added a whole Coxeter edge-orbit at a time, shortest first, until every root
#: carries at least this many. Targeting the degree rather than a fixed orbit count keeps the
#: invariant that matters when TUBE_FRAC changes -- a different embedding reorders the orbits
#: by length, so a fixed count silently leaves some shell under-braced (58 orbits gives every
#: root 12 on a ring torus but only 9 on this one).
MIN_BRACE_DEGREE = 12

COIL_R = 1.40       # swept coil radius, mm
BRACE_R = 1.00      # brace strut radius, mm
NODE_R = 1.90       # root node radius, mm
#: How far short of a root's centre a brace stops. Running braces all the way in would put
#: two solids' cap-centre vertices at the same point, and once a slicer or trimesh merges
#: coincident vertices that is a non-manifold vertex. Stopping short buries the cap inside
#: the node sphere instead -- sqrt(BRACE_INSET^2 + BRACE_R^2) = 1.49 < NODE_R -- so no two
#: solids in the file share a vertex and each one stays separately closed.
BRACE_INSET = 1.10

SUBDIV = 24         # coil samples per node-to-node span (30 * SUBDIV per coil)
COIL_SIDES = 10
BRACE_SIDES = 8
NODE_LAT, NODE_LON = 10, 14

PHI = (1 + 5 ** 0.5) / 2
TUBE = MAJOR * TUBE_FRAC


# ---------------------------------------------------------------- E8

def e8_roots():
    """The 240 roots: +-e_i +- e_j (112) and (+-1/2)^8 with an even number of minus signs
    (128). All of norm^2 = 2."""
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(8)
                v[i], v[j] = si, sj
                out.append(v)
    for s in itertools.product((0.5, -0.5), repeat=8):
        if sum(1 for x in s if x < 0) % 2 == 0:
            out.append(np.array(s))
    return np.array(out)


def simple_roots():
    """Bourbaki's simple roots for E8, in the same coordinates as e8_roots()."""
    a = np.zeros((8, 8))
    a[0] = np.array([1, -1, -1, -1, -1, -1, -1, 1]) / 2
    a[1] = np.array([1, 1, 0, 0, 0, 0, 0, 0])
    a[2] = np.array([-1, 1, 0, 0, 0, 0, 0, 0])
    for k in range(3, 8):
        a[k, k - 2] = -1
        a[k, k - 1] = 1
    return a


def coxeter_element(simples):
    """c = s_1 s_2 ... s_8. Orthogonal, of order 30 -- the Coxeter number of E8."""
    c = np.eye(8)
    for s in simples:
        n = s / np.linalg.norm(s)
        c = c @ (np.eye(8) - 2 * np.outer(n, n))
    return c


def coxeter_orbits(roots, c):
    """Partition the roots into orbits of c, each listed in c-order."""
    idx = {tuple(np.round(r, 6)): i for i, r in enumerate(roots)}
    nxt = np.array([idx[tuple(np.round(c @ r, 6))] for r in roots])
    orbits, seen = [], set()
    for i in range(len(roots)):
        if i in seen:
            continue
        o, j = [], i
        while j not in seen:
            seen.add(j)
            o.append(j)
            j = int(nxt[j])
        orbits.append(o)
    return orbits, nxt


def edge_orbits(edges, nxt):
    """Group E8's 60-degree edges into orbits of the Coxeter element.

    c permutes the roots, so it permutes the edges too: 6720 edges fall into 224 orbits of
    30. Keeping *whole* orbits rather than everything under a length cutoff is what makes
    the braced subset Coxeter-equivariant -- and it forces the brace degree to be constant
    on each of the 8 root orbits, because an edge orbit hands every root of the orbits it
    touches the same number of ends.
    """
    seen, out = set(), []
    for e in edges:
        if e in seen:
            continue
        o, cur = [], e
        while cur not in seen:
            seen.add(cur)
            o.append(cur)
            a, b = int(nxt[cur[0]]), int(nxt[cur[1]])
            cur = (min(a, b), max(a, b))
        out.append(o)
    return out


def cluster(x, tol=1e-9):
    """Group near-equal values; returns (cluster means, per-element label). Rounding to a
    fixed decimal place instead would inject 1e-10 of noise into the golden-ratio test,
    and the true spread inside a ring is 3e-15."""
    order = np.argsort(x)
    groups, lab = [], np.zeros(len(x), int)
    for i in order:
        if not groups or x[i] - groups[-1][-1] > tol:
            groups.append([x[i]])
        else:
            groups[-1].append(x[i])
        lab[i] = len(groups) - 1
    return np.array([np.mean(g) for g in groups]), lab


def eigenplane(c, exponent, order=30):
    """Unit eigenvector of c for the eigenvalue exp(2*pi*i*exponent/order)."""
    w, V = np.linalg.eig(c)
    m = np.round(np.angle(w) / (2 * np.pi) * order).astype(int) % order
    k = int(np.argmin(np.abs(m - exponent)))
    v = V[:, k]
    return v / np.linalg.norm(v)


# ---------------------------------------------------------------- primitives

def sphere(centre, r, nlat=NODE_LAT, nlon=NODE_LON):
    """Lat-long sphere, outward-wound, poles fanned."""
    c = np.asarray(centre, float)
    verts = [c + np.array([0.0, 0.0, r])]
    for i in range(1, nlat):
        s, z = math.sin(math.pi * i / nlat), math.cos(math.pi * i / nlat)
        for j in range(nlon):
            a = 2 * math.pi * j / nlon
            verts.append(c + r * np.array([s * math.cos(a), s * math.sin(a), z]))
    verts.append(c - np.array([0.0, 0.0, r]))
    faces = []
    for j in range(nlon):
        faces.append((0, 1 + j, 1 + (j + 1) % nlon))
    for i in range(nlat - 2):
        b0, b1 = 1 + i * nlon, 1 + (i + 1) * nlon
        for j in range(nlon):
            j2 = (j + 1) % nlon
            faces.append((b0 + j, b1 + j, b1 + j2))
            faces.append((b0 + j, b1 + j2, b0 + j2))
    last, p = 1 + (nlat - 2) * nlon, len(verts) - 1
    for j in range(nlon):
        faces.append((p, last + (j + 1) % nlon, last + j))
    return verts, faces


def coil_path(t, th0, ph0, wind, n):
    """Sample points and sweep frame of one (1, wind) coil at tube radius t.

    er is perpendicular to the tangent identically: differentiating
    (MAJOR + t*cos phi) e_rho + t*sin phi * z gives a vector in span(e_theta, e_poloidal),
    and both of those are orthogonal to er. So the frame needs no parallel transport and
    carries no holonomy -- the swept tube closes on itself exactly.
    """
    u = 2 * np.pi * np.arange(n) / n
    th, ph = th0 + u, ph0 + wind * u
    ct, st, cp, sp = np.cos(th), np.sin(th), np.cos(ph), np.sin(ph)
    zero = np.zeros(n)
    up = np.zeros((n, 3))
    up[:, 2] = 1.0
    erho = np.stack([ct, st, zero], 1)
    etheta = np.stack([-st, ct, zero], 1)
    er = cp[:, None] * erho + sp[:, None] * up
    epol = -sp[:, None] * erho + cp[:, None] * up
    P = MAJOR * erho + t * er
    T = (MAJOR + t * cp)[:, None] * etheta + (wind * t) * epol
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    return P, er, T


def coil_solid(P, er, T, r, sides=COIL_SIDES):
    """Sweep a closed tube of radius r along closed polyline P with frame (er, T x er)."""
    n = len(P)
    Wv = np.cross(T, er)
    a = 2 * np.pi * np.arange(sides) / sides
    ring = (np.cos(a)[None, :, None] * er[:, None, :]
            + np.sin(a)[None, :, None] * Wv[:, None, :])
    verts = (P[:, None, :] + r * ring).reshape(-1, 3)
    faces = []
    for k in range(n):
        k2 = (k + 1) % n
        for s in range(sides):
            s2 = (s + 1) % sides
            faces.append((k * sides + s, k * sides + s2, k2 * sides + s2))
            faces.append((k * sides + s, k2 * sides + s2, k2 * sides + s))
    return verts, faces


# ---------------------------------------------------------------- build

def densify(pts, step):
    """Resample a polyline so no gap exceeds `step`."""
    if len(pts) < 2:
        return pts
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    out = [pts[:1]]
    for k, L in enumerate(seg):
        m = max(1, int(np.ceil(L / step)))
        u = np.arange(1, m + 1)[:, None] / m
        out.append(pts[k] + u * (pts[k + 1] - pts[k]))
    return np.vstack(out)


def occupied_volume(paths, scale, n=1_500_000, seed=0):
    """Monte-Carlo volume of the union of every capsule and sphere in the lattice.

    Summing primitive volumes would badly over-count -- the solids interpenetrate at all
    240 joints -- and there is no mesh-boolean backend here to do it exactly, so sample
    instead. Centrelines are densified to 0.1 of the local radius, which keeps the
    sphere-union scalloping error under 1% of volume.
    """
    from scipy.spatial import cKDTree
    pts = [densify(q, 0.1 * r) for q, r in paths]
    rad = np.concatenate([np.full(len(q), r) for q, (_, r) in zip(pts, paths)])
    pts = np.vstack(pts)
    tree = cKDTree(pts)
    pad = rad.max()
    lo, hi = pts.min(0) - pad, pts.max(0) + pad
    rng = np.random.default_rng(seed)
    q = rng.uniform(lo, hi, (n, 3))
    d, i = tree.query(q, k=8, workers=-1)
    frac = float((d <= rad[i]).any(1).mean())
    return frac * float(np.prod(hi - lo)) * scale ** 3


def overhangs(paths, scale, flat_deg=20.0):
    """Total centreline length, how much of it lies within flat_deg of the build plate, and
    the longest unbroken shallow run -- which is the bridge the slicer has to span or prop."""
    total = flat = longest = 0.0
    for q, _ in paths:
        if len(q) < 2:
            continue
        d = np.diff(q, axis=0)
        L = np.linalg.norm(d, axis=1)
        shallow = np.degrees(np.arcsin(np.abs(d[:, 2]) / L)) < flat_deg
        total += L.sum()
        flat += L[shallow].sum()
        run = 0.0
        for ok, ln in zip(shallow, L):                 # paths are open; coils were closed
            run = run + ln if ok else 0.0              # by repeating the first point
            longest = max(longest, run)
    return total * scale, flat * scale, longest * scale


def coil_clearance(paths, ncoils=8):
    """Smallest surface gap between two different coils. They sit on different shells and
    share no roots, so they must not touch. One tree per coil and cross-queries only -- a
    single tree over all of them returns same-coil neighbours for any sane k.
    """
    from scipy.spatial import cKDTree
    pts = [densify(q, 0.15) for q, _ in paths[:ncoils]]
    trees = [cKDTree(q) for q in pts]
    best = np.inf
    for a in range(ncoils):
        for b in range(ncoils):
            if a != b:
                best = min(best, float(trees[b].query(pts[a], workers=-1)[0].min()))
    return best


def components(n, edges):
    p = list(range(n))

    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    for a, b in edges:
        p[find(a)] = find(b)
    return len({find(i) for i in range(n)})


def scale_guess(P):
    """Normalising factor SIZE/extent, needed before the mesh exists so clearance asserts can
    be stated in finished millimetres."""
    return SIZE / max(np.ptp(P[:, 0]), np.ptp(P[:, 1]))


def build(verbose=True, metrics=None):
    """metrics=False skips the Monte-Carlo volume, coil-clearance and overhang passes, which
    are ~30 s of the run and of no use to the preview scripts."""
    if metrics is None:
        metrics = verbose
    roots = e8_roots()
    cox = coxeter_element(simple_roots())

    order = next(k for k in range(1, 61)
                 if np.allclose(np.linalg.matrix_power(cox, k), np.eye(8), atol=1e-9))
    w = np.linalg.eigvals(cox)
    exps = sorted((np.round(np.angle(w) / (2 * np.pi) * order).astype(int) % order).tolist())
    assert order == 30 and exps == [1, 7, 11, 13, 17, 19, 23, 29], (order, exps)

    orbits, nxt = coxeter_orbits(roots, cox)
    assert len(orbits) == 8 and all(len(o) == 30 for o in orbits)

    z1 = roots @ np.conj(eigenplane(cox, 1))
    z7 = roots @ np.conj(eigenplane(cox, 7))
    r1, theta, phi = np.abs(z1), np.angle(z1), np.angle(z7)

    rings, shell = cluster(r1)
    assert len(rings) == 8
    frac = rings / rings.max()
    pairs = [(i, j) for i in range(8) for j in range(8)
             if i != j and abs(frac[j] - PHI * frac[i]) < 1e-11]
    assert len(pairs) == 4 and len({k for p in pairs for k in p}) == 8, pairs
    resid = max(abs(frac[j] - PHI * frac[i]) for i, j in pairs)

    # coil winding read off the data: theta advances one station, phi advances `wind`
    station = 2 * np.pi / 30
    coils = []
    for o in orbits:
        dth = np.diff(np.unwrap(theta[o])) / station
        dph = np.diff(np.unwrap(phi[o])) / station
        assert np.allclose(dth, dth[0], atol=1e-6) and np.allclose(dph, dph[0], atol=1e-6)
        assert abs(abs(dth[0]) - 1.0) < 1e-6
        coils.append((float(r1[o[0]]), float(theta[o[0]]), float(phi[o[0]]),
                      float(round(dph[0] / dth[0])), o))
    assert sorted({c[3] for c in coils}) == [7.0], sorted({c[3] for c in coils})

    # ---- node positions on the torus
    t_of = TUBE * r1 / r1.max()
    rho = MAJOR + t_of * np.cos(phi)
    P = np.stack([rho * np.cos(theta), rho * np.sin(theta), t_of * np.sin(phi)], 1)

    D = np.linalg.norm(P[:, None] - P[None, :], axis=2)
    np.fill_diagonal(D, np.inf)

    # ---- braces: shortest Coxeter orbits of 60-degree edges, until every root has enough
    G = np.round(roots @ roots.T).astype(int)
    cand = [(i, j) for i, j in itertools.combinations(range(240), 2) if G[i, j] == 1]
    eorb = edge_orbits(cand, nxt)
    assert len(eorb) == 224 and all(len(o) == 30 for o in eorb)
    span = np.array([np.mean([D[a, b] for a, b in o]) for o in eorb])
    picked, deg = [], np.zeros(240, int)
    for k in np.argsort(span):
        if deg.min() >= MIN_BRACE_DEGREE:
            break
        picked.append(int(k))
        deg = np.bincount(np.array([e for j in picked for e in eorb[j]]).ravel(),
                          minlength=240)
    norb = len(picked)
    braces = np.array([e for j in picked for e in eorb[j]])
    bl = np.linalg.norm(P[braces[:, 0]] - P[braces[:, 1]], axis=1)

    shell_deg = [int(deg[shell == k][0]) for k in range(8)]
    assert all(len(set(deg[shell == k].tolist())) == 1 for k in range(8)),         "brace degree is not constant on a root orbit -- edge orbits were split"
    assert deg.min() >= MIN_BRACE_DEGREE, deg.min()

    # Two braces leaving the same root in nearly the same direction put their inset cap
    # centres nearly on top of each other, and coincident vertices between two solids are
    # exactly what BRACE_INSET exists to avoid. Test the separation itself rather than the
    # angle: what matters is that no vertex merge could ever weld them.
    caps = {}
    for a, b in braces:
        for u, x in ((a, b), (b, a)):
            d = P[x] - P[u]
            caps.setdefault(int(u), []).append(P[u] + BRACE_INSET * d / np.linalg.norm(d))
    gap = min(float(np.min(np.linalg.norm(np.asarray(c)[:, None] - np.asarray(c)[None, :],
                                         axis=2) + 9e9 * np.eye(len(c))))
              for c in caps.values() if len(c) > 1)
    assert gap * scale_guess(P) > 0.02, f"brace cap centres only {gap:.4f} mm apart"

    ribs = [(i, int(nxt[i])) for i in range(240)]
    ncomp = components(240, [(int(a), int(b)) for a, b in braces] + ribs)
    assert ncomp == 1, f"{ncomp} disconnected pieces"
    assert len(set(shell[braces.ravel()])) == 8, "a shell is braced only to itself"

    # ---- mesh
    mesh = Mesh()
    paths = []
    n = 30 * SUBDIV
    for rr, th0, ph0, wind, o in coils:
        cp, cer, ct = coil_path(TUBE * rr / r1.max(), th0, ph0, wind, n)
        assert np.abs(np.einsum("ij,ij->i", cer, ct)).max() < 1e-12
        assert np.abs(cp[::SUBDIV] - P[o]).max() < 1e-9   # coil hits its own 30 roots
        mesh.add_solid(*coil_solid(cp, cer, ct, COIL_R), tag="coil",
                       wall=int(shell[o[0]]))
        paths.append((np.vstack([cp, cp[:1]]), COIL_R))
    assert math.hypot(BRACE_INSET, BRACE_R) < NODE_R - 0.2
    for a, b in braces:
        d = P[b] - P[a]
        d /= np.linalg.norm(d)
        p0, p1 = P[a] + BRACE_INSET * d, P[b] - BRACE_INSET * d
        mesh.add_solid(*tube(p0, p1, BRACE_R, nseg=BRACE_SIDES), tag="brace", wall=8)
        k = max(2, int(np.linalg.norm(p1 - p0) / 0.25))
        paths.append((p0 + np.linspace(0, 1, k)[:, None] * (p1 - p0), BRACE_R))
    for i in range(240):
        mesh.add_solid(*sphere(P[i], NODE_R), tag="node", wall=int(shell[i]))
        paths.append((P[i][None, :], NODE_R))

    bad = mesh.validate()
    assert not bad, f"{len(bad)} non-manifold solids"

    # ---- normalise: SIZE across, sitting on z = 0, torus axis left on x = y = 0
    #
    # Do NOT recentre in x-y. The roots are not distributed symmetrically about the axis --
    # the point cloud's bounding-box centre sits 3.5 mm off it -- so centring the part on its
    # bounding box slides it sideways off its own axis and the hole comes out eccentric,
    # pinched to 26 mm on one side and open to 20 mm past centre on the other. The axis is
    # the feature worth keeping concentric; slicers centre the plate themselves.
    v = np.asarray(mesh.v, float)
    ext = v.max(0) - v.min(0)
    scale = SIZE / max(ext[0], ext[1])
    zfloor = v[:, 2].min()
    v[:, 2] -= zfloor
    mesh.v = [tuple(p) for p in v * scale]
    P = P.copy()
    P[:, 2] -= zfloor                                  # keep info["P"] in the mesh's frame

    # Concentricity check, and it only passes because the part was not slid off its axis
    # above: the material's radial extent has to match the analytic envelope in every
    # direction. One root sits exactly on the inner equator of the outer shell, so on a ring
    # torus that root's sphere is the narrowest point of the hole; on the horn torus the
    # inner equator has collapsed onto the axis, so the same root sits at the origin and the
    # hole radius is 0.
    skin = max(NODE_R, COIL_R)
    rho_v = np.hypot(*np.asarray(mesh.v, float)[:, :2].T)
    hole = float(rho_v.min())
    assert abs(rho_v.max() - (MAJOR + TUBE + skin) * scale) < 0.05, rho_v.max()
    assert abs(hole - max(0.0, (MAJOR - TUBE - skin) * scale)) < 0.05, hole
    kind = ("ring" if MAJOR > TUBE + 1e-9 else
            "horn" if abs(MAJOR - TUBE) < 1e-9 else "spindle")

    info = dict(order=order, exps=exps, shell=shell, pairs=pairs, scale=scale,
                overhang=overhangs(paths, scale) if metrics else None,
                coil_gap=(coil_clearance(paths) - 2 * COIL_R) * scale if metrics else None,
                resid=resid,
                volume=occupied_volume(paths, scale) if metrics else None,
                P=P * scale, D=D * scale, braces=braces, bl=bl * scale,
                cand=len(cand), t_of=t_of * scale, coils=coils, nxt=nxt, hole=hole,
                kind=kind, hub_roots=int((np.hypot(P[:, 0], P[:, 1]) < 1e-6).sum()),
                hub_passes=int(round(coils[0][3])),
                shell_deg=shell_deg, deg=deg, norb=norb)
    if verbose:
        report(mesh, info)
    return mesh, info


def report(mesh, info):
    v, f = mesh.arrays()
    ext = v.max(0) - v.min(0)
    s = info["scale"]
    rad = np.unique(np.round(info["t_of"], 6))
    print("E8 root system")
    print(f"  240 roots   Coxeter order {info['order']}   exponents {info['exps']}")
    print("  8 Coxeter orbits of 30, each a (1, 7) coil on the torus")
    print("\ntorus")
    print(f"  major radius {MAJOR * s:.2f} mm   tube radius {TUBE * s:.2f} mm"
          f"   (normalising scale {s:.4f})")
    print("  shell radii, in golden-ratio pairs:")
    for i, j in sorted(info["pairs"]):
        print(f"    {rad[i]:6.2f} mm : {rad[j]:6.2f} mm   ratio {rad[j] / rad[i]:.10f}"
              f"   phi = {PHI:.10f}")
    print(f"    worst deviation from exactly phi: {info['resid']:.2e}")
    print("\nstruts")
    print(f"  coils   8 closed swept tubes, r {COIL_R * s:.2f} mm")
    print(f"  braces  {len(info['braces'])} of {info['cand']} 60-degree root pairs "
          f"= {100 * len(info['braces']) / info['cand']:.1f}% of E8's edge set, "
          f"{info['bl'].min():.1f}-{info['bl'].max():.1f} mm, r {BRACE_R * s:.2f} mm")
    print(f"          {info['norb']} of 224 Coxeter edge-orbits, whole; brace degree per "
          f"shell {info['shell_deg']} (E8's own degree is 56)")
    print(f"  nodes   240 spheres, r {NODE_R * s:.2f} mm")
    print(f"  closest two roots {info['D'].min():.2f} mm apart, node surfaces clearing "
          f"{info['D'].min() - 2 * NODE_R * s:.2f} mm")
    print("\nmesh")
    print(f"  {len(mesh._solid_start)} solids, {len(f)} triangles, every solid edge-manifold")
    print(f"  bbox {ext[0]:.2f} x {ext[1]:.2f} x {ext[2]:.2f} mm   (bed is 320 x 320 x 320)")
    R, a = MAJOR * s, TUBE * s
    print(f"  {info['kind']} torus: R {R:.2f}, a {a:.2f}, a/R {a / R:.3f}; hole dia "
          f"{2 * info['hole']:.1f} mm, concentric with the axis")
    if info["kind"] == "horn":
        print(f"  the hole has closed to a point: the outer shell's inner equator collapses "
              f"onto the axis,")
        print(f"  so {info['hub_roots']} root and {info['hub_passes']} passes of the outer "
              f"coil meet at the origin -- a {info['hub_passes']}-fold hub, one pass per "
              f"poloidal turn")
    vol = info["volume"]
    print(f"  solid volume {vol / 1000:.1f} cm^3 -> {vol * 1.27e-3:.0f} g in PETG"
          f"   ({100 * vol / (ext[0] * ext[1] * ext[2]):.1f}% of the bounding box)")
    low = v[v[:, 2] < 0.4]
    patch, _ = cluster(np.arctan2(low[:, 1], low[:, 0]), tol=0.15)
    total, flatlen, bridge = info["overhang"]
    print()
    print("printability")
    print(f"  centreline length {total / 1000:.2f} m; {100 * flatlen / total:.0f}% of it lies "
          f"within 20 deg of the plate")
    print(f"  longest shallow run {bridge:.1f} mm -- supports or an organic-support profile "
          f"are needed, this will not bridge dry")
    print(f"  closest approach between two different coils: {info['coil_gap']:.2f} mm of air")
    print(f"  bed contact: {len(patch)} patches on a circle of diameter "
          f"{2 * np.hypot(low[:, 0], low[:, 1]).mean():.0f} mm -- the outer coil's "
          f"poloidal low points")


def main():
    mesh, info = build()
    out = f"e8_{info['kind']}_torus_lattice_150mm.stl"   # name states the measured class
    n = write_stl(out, mesh)
    print(f"\nwrote {out}  ({n} triangles, {os.path.getsize(out) / 1e6:.1f} MB)")

    import trimesh
    m = trimesh.load(out)
    print(f"  trimesh: extents {np.round(m.extents, 2)}   watertight {m.is_watertight}"
          f"   winding consistent {m.is_winding_consistent}   closed shells {m.body_count}")
    fits = np.all(m.extents <= 320.0)
    print(f"  fits the 320 mm bed: {fits}")
    assert m.is_watertight and m.is_winding_consistent and fits
    return mesh, info


if __name__ == "__main__":
    main()
