#!/usr/bin/env python3
"""Mesh primitives for the crystal-lattice models.

Everything here produces triangle soup for `mesh_kit.Mesh.add_solid`, following the same
contract as the rest of the repo: each solid is individually watertight and outward-wound,
and parts assembled from several overlapping solids rely on the slicer's union rather than
a CSG boolean. That is not laziness -- this environment has no working mesh-boolean backend
(no manifold3d, fcl or rtree), so overlapping watertight solids is the only route to an
embossed label or a sleeved hub. `mesh_kit.Mesh.validate()` is the real correctness check;
`trimesh.is_watertight` on a whole multi-solid part will read False by construction.
"""

import math
import struct

import numpy as np
from scipy.spatial import ConvexHull


# ---------------------------------------------------------------- output

def write_stl(path, mesh, header="PrintableSacredGeometry crystal lattices"):
    """Binary STL. Same as `mesh_kit.write_stl` but with a caller-supplied header, so
    each file says what it is instead of inheriting another design's label."""
    v, f = mesh.arrays()
    tri = v[f]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)
    with open(path, "wb") as fh:
        fh.write(header.encode()[:79].ljust(80, b" "))
        fh.write(struct.pack("<I", len(f)))
        for i in range(len(f)):
            fh.write(struct.pack("<12fH", *nrm[i], *tri[i, 0], *tri[i, 1], *tri[i, 2], 0))
    return len(f)


# ---------------------------------------------------------------- winding

def orient_outward(verts, faces, centre=None):
    """Flip any triangle whose normal points toward `centre`. Valid only for convex
    solids, which covers everything here except the sleeve (wound by hand, because its
    bore must face inward)."""
    v = np.asarray(verts, float)
    f = np.asarray(faces, np.int64)
    c = v.mean(axis=0) if centre is None else np.asarray(centre, float)
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    outward = np.einsum("ij,ij->i", n, tri.mean(axis=1) - c) > 0.0
    f = np.where(outward[:, None], f, f[:, ::-1])
    return v, f


# ---------------------------------------------------------------- convex bodies

def parallelepiped(A, B, C, origin=(0.0, 0.0, 0.0)):
    """The unit cell itself: eight corners spanned by the three cell vectors."""
    A, B, C = (np.asarray(x, float) for x in (A, B, C))
    o = np.asarray(origin, float)
    verts = [o + i * A + j * B + k * C
             for k in (0, 1) for j in (0, 1) for i in (0, 1)]
    quads = [(0, 1, 3, 2), (4, 5, 7, 6),        # k = 0, k = 1
             (0, 1, 5, 4), (2, 3, 7, 6),        # j = 0, j = 1
             (0, 2, 6, 4), (1, 3, 7, 5)]        # i = 0, i = 1
    faces = []
    for a, b, c, d in quads:
        faces += [(a, b, c), (a, c, d)]
    return orient_outward(verts, faces)


def hull_solid(points):
    """Convex hull of a point cloud, triangulated and outward-wound."""
    p = np.asarray(points, float)
    h = ConvexHull(p)
    return orient_outward(p, h.simplices)


def prism_on_plane(poly2d, origin, ex, ey, normal, depth):
    """Extrude a convex 2D polygon sitting on an arbitrary plane, `depth` along `normal`.
    Used for the embossed Pearson glyphs."""
    origin = np.asarray(origin, float)
    ex, ey = np.asarray(ex, float), np.asarray(ey, float)
    n = np.asarray(normal, float) * depth
    base = [origin + x * ex + y * ey for x, y in poly2d]
    verts = base + [p + n for p in base]
    m = len(base)
    faces = []
    for i in range(1, m - 1):                       # caps
        faces += [(0, i, i + 1), (m, m + i, m + i + 1)]
    for i in range(m):                              # walls
        j = (i + 1) % m
        faces += [(i, j, m + j), (i, m + j, m + i)]
    return orient_outward(verts, faces)


# ---------------------------------------------------------------- kit hardware

def hub(r, cut=-0.78, nu=28, nv=9):
    """A sphere truncated by a horizontal plane at z = cut*r and closed with a flat disc.

    A whole sphere touches the bed at a point and peels off; the truncation gives a
    pad of diameter 2*r*sqrt(1 - cut^2), about 8.7 mm at r = 7, which sticks."""
    phi_max = math.acos(cut)
    verts = [(0.0, 0.0, r)]
    for i in range(1, nv + 1):
        phi = phi_max * i / nv
        sz, sr = r * math.cos(phi), r * math.sin(phi)
        for k in range(nu):
            t = 2.0 * math.pi * k / nu
            verts.append((sr * math.cos(t), sr * math.sin(t), sz))
    verts.append((0.0, 0.0, r * cut))
    hubc = len(verts) - 1

    def ring(i, k):
        return 1 + (i - 1) * nu + (k % nu)

    faces = []
    for k in range(nu):
        faces.append((0, ring(1, k), ring(1, k + 1)))
    for i in range(1, nv):
        for k in range(nu):
            faces += [(ring(i, k), ring(i + 1, k), ring(i + 1, k + 1)),
                      (ring(i, k), ring(i + 1, k + 1), ring(i, k + 1))]
    for k in range(nu):
        faces.append((hubc, ring(nv, k + 1), ring(nv, k)))
    return orient_outward(verts, faces)


def sphere(r, nu=28, nv=16):
    """Plain UV sphere, for the fused one-piece frames where a hub is not a socket and
    needs no bed pad -- only the hubs actually touching the bed use `hub` instead."""
    verts = [(0.0, 0.0, r)]
    for i in range(1, nv):
        phi = math.pi * i / nv
        sz, sr = r * math.cos(phi), r * math.sin(phi)
        for k in range(nu):
            t = 2.0 * math.pi * k / nu
            verts.append((sr * math.cos(t), sr * math.sin(t), sz))
    verts.append((0.0, 0.0, -r))
    south = len(verts) - 1

    def ring(i, k):
        return 1 + (i - 1) * nu + (k % nu)

    faces = []
    for k in range(nu):
        faces.append((0, ring(1, k), ring(1, k + 1)))
        faces.append((south, ring(nv - 1, k + 1), ring(nv - 1, k)))
    for i in range(1, nv - 1):
        for k in range(nu):
            faces += [(ring(i, k), ring(i + 1, k), ring(i + 1, k + 1)),
                      (ring(i, k), ring(i + 1, k + 1), ring(i, k + 1))]
    return orient_outward(verts, faces)


def sleeve(direction, r0, r1, r_in, r_out, nseg=20):
    """A socket: an open-bore tube standing off along `direction` from r0 to r1.

    Wound by hand rather than through `orient_outward`, because the bore wall has to face
    inward. r0 sits inside the hub sphere, so the hub itself closes the bottom of the bore
    and the effective socket depth is r1 minus the hub radius -- no boolean subtraction
    needed anywhere."""
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    tmp = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(d, tmp)
    u /= np.linalg.norm(u)
    w = np.cross(d, u)

    rings = [u * math.cos(2 * math.pi * k / nseg) + w * math.sin(2 * math.pi * k / nseg)
             for k in range(nseg)]
    verts = ([r0 * d + r_out * e for e in rings] + [r1 * d + r_out * e for e in rings]
             + [r0 * d + r_in * e for e in rings] + [r1 * d + r_in * e for e in rings])
    n = nseg
    A, B, C, D = 0, n, 2 * n, 3 * n

    faces = []
    for k in range(n):
        k2 = (k + 1) % n
        faces += [(A + k, A + k2, B + k2), (A + k, B + k2, B + k)]        # outer wall
        faces += [(C + k, D + k2, C + k2), (C + k, D + k, D + k2)]        # bore, inward
        faces += [(B + k, B + k2, D + k2), (B + k, D + k2, D + k)]        # cap at r1
        faces += [(A + k, C + k, C + k2), (A + k, C + k2, A + k2)]        # cap at r0
    return verts, faces


# ---------------------------------------------------------------- orientation

def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def align(u, v):
    """Rotation matrix carrying unit vector u onto unit vector v (Rodrigues)."""
    u = np.asarray(u, float)
    v = np.asarray(v, float)
    c = float(u @ v)
    if c > 1.0 - 1e-12:
        return np.eye(3)
    if c < -1.0 + 1e-12:
        a = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        ax = np.cross(u, a)
        ax /= np.linalg.norm(ax)
        return 2.0 * np.outer(ax, ax) - np.eye(3)
    ax = np.cross(u, v)
    s = np.linalg.norm(ax)
    ax /= s
    K = np.array([[0.0, -ax[2], ax[1]], [ax[2], 0.0, -ax[0]], [-ax[1], ax[0], 0.0]])
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def best_up(dirs, samples=1200):
    """The 'up' that maximises the lowest socket elevation, i.e. the print orientation
    whose worst-drooping sleeve droops least.

    A body-centre node cannot be oriented well by pointing the mean socket direction up --
    for a tetrahedral spider the mean is exactly zero. Searching orientations instead finds
    the 2-fold axis, which leaves the worst sleeve at -35.3 degrees rather than -90."""
    d = np.asarray(dirs, float)
    cands = fibonacci_sphere(samples)
    elev = np.arcsin(np.clip(cands @ d.T, -1.0, 1.0))       # (samples, nsockets)
    return cands[int(np.argmax(elev.min(axis=1)))]


# ---------------------------------------------------------------- print analysis

def overhang_slope(verts, faces, bed_tol=0.25):
    """Shallowest downward-facing surface, in degrees from horizontal, ignoring faces
    resting on the bed.

    A vertical wall is 90 degrees and needs nothing; a horizontal ceiling is 0 and needs
    support; slicers give up somewhere around 45. Returns (min_slope_deg, area_below_45)."""
    v = np.asarray(verts, float)
    f = np.asarray(faces, np.int64)
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(n, axis=1)
    keep = area > 1e-12
    n, area, tri, f = n[keep], area[keep], tri[keep], f[keep]
    n = n / area[:, None]

    zmin = v[:, 2].min()
    down = (n[:, 2] < -1e-6) & (tri[:, :, 2].max(axis=1) > zmin + bed_tol)
    if not down.any():
        return 90.0, 0.0
    slope = np.degrees(np.arccos(np.clip(np.abs(n[down, 2]), 0.0, 1.0)))
    a = area[down] / 2.0
    return float(slope.min()), float(a[slope < 45.0].sum())


def face_frame(verts, faces, target=(0.0, 0.0, 1.0), tol=0.02):
    """The flat face whose normal is closest to `target`, as
    (centroid, ex, ey, normal, points in local 2D relative to the centroid).

    Returning the outline rather than a half-width lets the label fitter treat a rhombus,
    a hexagon and a square the same way -- it tests containment against the real face
    boundary instead of assuming the face is a rectangle."""
    v = np.asarray(verts, float)
    f = np.asarray(faces, np.int64)
    t = np.asarray(target, float)
    t = t / np.linalg.norm(t)
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1)
    n = n / np.maximum(ln, 1e-12)[:, None]
    dot = n @ t
    best = dot.max()
    sel = dot > best - tol
    normal = n[sel].mean(axis=0)
    normal /= np.linalg.norm(normal)
    pts = tri[sel].reshape(-1, 3)
    centre = pts.mean(axis=0)

    seed = np.array([1.0, 0.0, 0.0])
    if abs(float(seed @ normal)) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    ex = seed - float(seed @ normal) * normal
    ex /= np.linalg.norm(ex)
    ey = np.cross(normal, ex)
    local = pts - centre
    return centre, ex, ey, normal, np.stack([local @ ex, local @ ey], axis=1)
