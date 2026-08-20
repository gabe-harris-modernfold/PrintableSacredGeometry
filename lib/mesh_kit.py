#!/usr/bin/env python3
"""Mesh primitives shared by the geometry generators.

Extracted verbatim from the earlier per-design modules so the live design does
not import its mesh kit from a superseded one. Solids are individually watertight
and outward-wound; Mesh.validate() checks edge-manifoldness per solid."""

import math
import struct

import numpy as np


class Mesh:
    """Accumulates triangle soup, tagged by solid, with per-solid validation."""

    def __init__(self):
        self.v = []
        self.f = []
        self.tag = []
        self.wall = []                     # owning wall index, or -1
        self._solid_start = []

    def add_solid(self, verts, faces, tag, wall=-1):
        base = len(self.v)
        self.v.extend([tuple(map(float, p)) for p in verts])
        self._solid_start.append((base, len(self.f), len(faces)))
        for a, b, c in faces:
            self.f.append((base + a, base + b, base + c))
            self.tag.append(tag)
            self.wall.append(wall)

    def arrays(self):
        return np.asarray(self.v, float), np.asarray(self.f, np.int64)

    def pick_walls(self, walls):
        """Faces belonging to whole walls -- gives clean section views without
        the shards a centroid-plane cull leaves behind."""
        v, f = self.arrays()
        return v, f[np.isin(np.asarray(self.wall), list(walls))]

    def pick(self, walls=None, tags=None, drop_tags=None):
        v, f = self.arrays()
        keep = np.ones(len(f), bool)
        if walls is not None:
            keep &= np.isin(np.asarray(self.wall), list(walls))
        if tags is not None:
            keep &= np.isin(np.asarray(self.tag), list(tags))
        if drop_tags is not None:
            keep &= ~np.isin(np.asarray(self.tag), list(drop_tags))
        return v, f[keep]

    def validate(self):
        """Every solid must be edge-manifold: each undirected edge used twice,
        once in each direction."""
        bad = []
        for base, f0, nf in self._solid_start:
            edges = {}
            for a, b, c in self.f[f0:f0 + nf]:
                for e in ((a, b), (b, c), (c, a)):
                    edges[e] = edges.get(e, 0) + 1
            for (a, b), n in edges.items():
                if n != 1 or edges.get((b, a), 0) != 1:
                    bad.append((base, a, b, n))
                    break
        return bad

def ensure_ccw(poly):
    p = np.asarray(poly, float)
    area2 = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
    return p if area2 > 0 else p[::-1].copy()

def prism(poly_xy, z0, z1):
    """Extrude a CCW convex/simple-star polygon along +z. Fan-triangulated caps."""
    p = ensure_ccw(poly_xy)
    n = len(p)
    verts = [(x, y, z0) for x, y in p] + [(x, y, z1) for x, y in p]
    faces = []
    for i in range(1, n - 1):                       # bottom cap, normal -z
        faces.append((0, i + 1, i))
    for i in range(1, n - 1):                       # top cap, normal +z
        faces.append((n, n + i, n + i + 1))
    for i in range(n):                              # sides
        j = (i + 1) % n
        faces.append((i, j, n + j))
        faces.append((i, n + j, n + i))
    return verts, faces

def tube(p0, p1, r, nseg=8, cap=True):
    """Solid cylinder from p0 to p1."""
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    ax = p1 - p0
    L = np.linalg.norm(ax)
    ax = ax / L
    tmp = np.array([0.0, 0.0, 1.0]) if abs(ax[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(ax, tmp)
    u /= np.linalg.norm(u)
    w = np.cross(ax, u)
    ring = [u * math.cos(2 * math.pi * k / nseg) + w * math.sin(2 * math.pi * k / nseg)
            for k in range(nseg)]
    verts = [p0 + r * d for d in ring] + [p1 + r * d for d in ring]
    verts += [p0, p1]
    faces = []
    for k in range(nseg):
        k2 = (k + 1) % nseg
        faces.append((k, k2, nseg + k2))
        faces.append((k, nseg + k2, nseg + k))
        if cap:
            faces.append((2 * nseg, k2, k))
            faces.append((2 * nseg + 1, nseg + k, nseg + k2))
    return verts, faces

def write_stl(path, mesh):
    v, f = mesh.arrays()
    tri = v[f]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)
    with open(path, "wb") as fh:
        fh.write(b"hex mirror tube module - PrintableSacredGeometry".ljust(80, b" "))
        fh.write(struct.pack("<I", len(f)))
        for i in range(len(f)):
            fh.write(struct.pack("<12fH", *nrm[i], *tri[i, 0], *tri[i, 1],
                                 *tri[i, 2], 0))
    return len(f)

def _bisectors(poly):
    """Per-vertex (unit inward bisector, half interior angle)."""
    out = []
    n = len(poly)
    for i in range(n):
        vi = np.array(poly[i], float)
        vj = np.array(poly[(i + 1) % n], float)
        vk = np.array(poly[(i - 1) % n], float)
        e1 = (vj - vi) / np.linalg.norm(vj - vi)
        e2 = (vk - vi) / np.linalg.norm(vk - vi)
        half = math.acos(max(-1.0, min(1.0, float(e1 @ e2)))) / 2.0
        b = e1 + e2
        out.append((vi, b / np.linalg.norm(b), half))
    return out

def inset_poly(poly, d):
    """Shrink a polygon by perpendicular distance d on every edge."""
    return [vi + (d / math.sin(half)) * b for vi, b, half in _bisectors(poly)]

def normals(verts, faces):
    """Inward facet normals, and inward vertex normals (area-weighted mean)."""
    t = verts[faces]
    fn = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    area = np.linalg.norm(fn, axis=1, keepdims=True)
    fn = fn / area
    c = t.mean(axis=1)
    fn[np.einsum("ij,ij->i", fn[:, :2], c[:, :2]) > 0] *= -1.0   # radial test only
    vn = np.zeros_like(verts)
    for f, n, a in zip(faces, fn, area[:, 0]):
        vn[list(f)] += n * a
    vn /= np.linalg.norm(vn, axis=1, keepdims=True)
    return fn, vn

def pad_sites(tri, n):
    """Three support points on one facet, inset along the in-plane bisectors far
    enough that the pad edge clears every glass edge. Pad radius scales with the
    facet, because these run from ~52 mm down to ~12 mm a side."""
    side = float(np.mean([np.linalg.norm(tri[(i + 1) % 3] - tri[i])
                          for i in range(3)]))
    pr = float(np.clip(side / 9.0, 0.9, 2.5))
    out = []
    for i in range(3):
        a, b, c = tri[i], tri[(i + 1) % 3], tri[(i - 1) % 3]
        e1 = (b - a) / np.linalg.norm(b - a)
        e2 = (c - a) / np.linalg.norm(c - a)
        half = math.acos(float(np.clip(e1 @ e2, -1, 1))) / 2.0
        bis = e1 + e2
        bis /= np.linalg.norm(bis)
        d = max(side / 6.0, (pr + 0.8) / math.sin(half))
        out.append((a + d * bis, pr))
    return out
