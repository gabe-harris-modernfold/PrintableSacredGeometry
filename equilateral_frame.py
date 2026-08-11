#!/usr/bin/env python3
"""
Printable space frame for the (3,7) phyllotactic equilateral vortex tube.

The topology is HELICAL, not ring-based, so the ring-stack frames in
vortex_frame.py do not transfer. Vertices lie on one helix v_k = (R cos k*alpha,
R sin k*alpha, k*beta); each joins k+-3, k+-7, k+-10, giving two triangle families
(k, k+3, k+10) and (k, k+7, k+10). Euler: V - E + F = 56 - 148 + 92 = 0, an open
tube.

Consequences of the helix that the frame has to absorb:
  * Both ends are RAGGED -- the boundary is a helical staircase over 10 vertices,
    not a flat ring. So the base and rim are separate flat collars joined to the
    end nodes by posts of individually computed length.
  * There is no discrete rotational symmetry, only a screw symmetry, so nothing
    in the frame can be built by rotating one sector.

Everything structural stays OUTSIDE the mirror plane, so the bore is 100% mirror
and no strut is ever in the light path -- verified at the end against the facet
planes directly, not by a radial probe (which under-reads on twisted geometry).

Run:  python equilateral_frame.py
"""

import json
import math

import numpy as np

import mesh_kit as MK

GLASS = 1.5
RELIEF = 0.6
STRUT_R = 2.6
NODE_R = 3.4
PAD_LEN = 12.0          # must span from the strut ring in to the glass
BASE_T = 8.0            # flat collar thickness
RIM_T = 7.0
# Collars must start OUTBOARD of the largest wall radius R, not inboard of it: an
# inner rim at R-6 punched 11.9 mm through the bore where the wall dips toward its
# 52.4 mm inradius. Growing them outward also widens the footprint, which a 287 mm
# tall object wants anyway.
COLLAR_R0 = 14.0         # collar inner radius = R + this
COLLAR_R1 = 30.0        # collar outer radius = R + this
COLLAR_GAP = 2.0        # axial gap. Counter-intuitively this must stay SMALL:
                        # the end facets' outward normals tilt upward, so axial
                        # separation SUBTRACTS from clearance while radial
                        # separation adds. Raising the gap 3 -> 6.5 mm made the
                        # rim worse (0.542 -> 0.312). Radial is the only lever.
POST_R = 2.8
NODE_OFF = 8.0          # RADIAL node offset; see build()
NCOLLAR = 24            # polygon segments approximating the round collars


def load(path="equilateral_tube_solution.json"):
    return json.load(open(path))


def topology(sol):
    p, r, K = sol["p"], sol["r"], sol["K"]
    R, al, be = sol["R"], sol["alpha"], sol["beta"]
    q = p + r
    k = np.arange(K)
    verts = np.stack([R * np.cos(k * al), R * np.sin(k * al), k * be], axis=1)

    faces, edges = [], set()
    for i in range(K - q):
        faces.append((i, i + p, i + q))
        faces.append((i, i + r, i + q))
        edges |= {frozenset((i, i + p)), frozenset((i + p, i + q)),
                  frozenset((i, i + q)), frozenset((i, i + r)),
                  frozenset((i + r, i + q))}
    faces = np.asarray(faces, int)
    edges = np.asarray([sorted(e) for e in edges], int)
    assert len(verts) - len(edges) + len(faces) == 0, \
        f"not an open tube: V={len(verts)} E={len(edges)} F={len(faces)}"
    return verts, faces, edges, q


def annulus(z0, z1, r_in, r_out, n=NCOLLAR):
    """Flat round collar as n mitred quad prisms."""
    out = []
    for m in range(n):
        a0 = 2 * np.pi * m / n
        a1 = 2 * np.pi * (m + 1) / n
        poly = [(r_in * math.cos(a0), r_in * math.sin(a0)),
                (r_in * math.cos(a1), r_in * math.sin(a1)),
                (r_out * math.cos(a1), r_out * math.sin(a1)),
                (r_out * math.cos(a0), r_out * math.sin(a0))]
        out.append(MK.prism(poly, z0, z1))
    return out


def build(sol):
    verts, faces, edges, q = topology(sol)
    fn, vn = MK.normals(verts, faces)
    # Offset nodes RADIALLY, not along the vertex normal. Near the ragged helical
    # ends a vertex touches fewer facets, so its area-weighted normal stops
    # pointing outward and pushed 179 frame vertices up to 11.9 mm into the bore.
    # Every vertex sits at radius R, so a radial push guarantees the node clears
    # the largest wall radius there is. The (k,k+3) strut is the binding one: its
    # ends span ~38.8 deg of azimuth, so it dips to cos(19.4 deg) = 0.943 of the
    # node radius, which is why NODE_OFF has to be well above STRUT_R.
    rad = verts.copy()
    rad[:, 2] = 0.0
    rad /= np.linalg.norm(rad, axis=1, keepdims=True)
    nodes = verts + rad * NODE_OFF

    m = MK.Mesh()
    for a, b in edges:
        m.add_solid(*MK.tube(nodes[a], nodes[b], STRUT_R, nseg=7), tag="strut")
    for i, pnt in enumerate(nodes):
        m.add_solid(*MK.tube(pnt - vn[i] * NODE_R * 0.5,
                             pnt + vn[i] * NODE_R * 0.5, NODE_R, nseg=8),
                    tag="node")

    for f, n in zip(faces, fn):
        for site, pr in MK.pad_sites(verts[f], n):
            m.add_solid(*MK.tube(site - n * (GLASS + PAD_LEN), site - n * GLASS,
                                 pr, nseg=8), tag="pad")

    R = sol["R"]
    ri, ro = R + COLLAR_R0, R + COLLAR_R1
    zlo, zhi = float(verts[:, 2].min()), float(verts[:, 2].max())
    zb, zt = zlo - COLLAR_GAP, zhi + COLLAR_GAP
    for prism in annulus(zb - BASE_T, zb, ri, ro):
        m.add_solid(*prism, tag="base")
    for prism in annulus(zt, zt + RIM_T, ri, ro):
        m.add_solid(*prism, tag="rim")

    # Ragged ends: post each of the lowest/highest q nodes to its collar face.
    low = np.argsort(verts[:, 2])[:q]
    high = np.argsort(verts[:, 2])[-q:]
    # Posts SLANT outward to the collar mid-radius rather than dropping vertically.
    # Vertical posts at the node radius came within 1.225 mm of the mirror near the
    # ragged ends -- inside the 1.5 mm glass, so they would have fouled the mirror
    # rather than the frame. Slanting also triangulates the joint.
    rmid = R + 0.5 * (COLLAR_R0 + COLLAR_R1)
    for i in low:
        a = nodes[i]
        u = np.array([a[0], a[1], 0.0])
        u /= np.linalg.norm(u)
        m.add_solid(*MK.tube(u * rmid + np.array([0, 0, zb - 0.01]), a, POST_R,
                             nseg=7), tag="post")
    for i in high:
        a = nodes[i]
        u = np.array([a[0], a[1], 0.0])
        u /= np.linalg.norm(u)
        m.add_solid(*MK.tube(a, u * rmid + np.array([0, 0, zt + RIM_T + 0.01]),
                             POST_R, nseg=7), tag="post")
    return m, verts, faces, fn


def clearance(m, verts, faces, fn):
    """Signed clearance of every frame vertex to the mirror SURFACE.

    Must be true point-to-triangle distance, not point-to-plane. Two earlier
    versions were wrong: a radial probe from the axis under-reads on twisted
    geometry, and taking the nearest facet PLANE falls through to the opposite
    wall (reading -147 mm, ~one tube diameter) whenever a node's perpendicular
    foot lands outside its local triangle. For a point outside a triangle the
    closest point is always on the boundary, so: inside -> use the projection,
    outside -> use the nearest of the three edge segments. Exact either way.
    """
    v, _ = m.arrays()
    tri = verts[faces]
    A, B, C = tri[:, 0], tri[:, 1], tri[:, 2]
    e1, e2 = B - A, C - A
    d11 = np.einsum("ij,ij->i", e1, e1)
    d22 = np.einsum("ij,ij->i", e2, e2)
    d12 = np.einsum("ij,ij->i", e1, e2)
    den = d11 * d22 - d12 * d12
    segs = [(A, B), (B, C), (C, A)]
    seg_d = [(q - pp) for pp, q in segs]
    seg_l2 = [np.einsum("ij,ij->i", d, d) for d in seg_d]

    tagv = _vertex_tags(m, len(v))
    per = {}
    worst, bad = 1e9, 0
    for vi, pnt in enumerate(v):
        w = pnt - A
        sd = -np.einsum("ij,ij->i", w, fn)           # + = outboard of the mirror
        proj = w + fn * sd[:, None]
        b1 = np.einsum("ij,ij->i", proj, e1)
        b2 = np.einsum("ij,ij->i", proj, e2)
        u = (b1 * d22 - b2 * d12) / den
        vv = (b2 * d11 - b1 * d12) / den
        inside = (u >= 0) & (vv >= 0) & (u + vv <= 1)
        cp = A + e1 * u[:, None] + e2 * vv[:, None]
        best = np.where(inside, np.linalg.norm(pnt - cp, axis=1), np.inf)
        for (pp, _), d, l2 in zip(segs, seg_d, seg_l2):
            t = np.clip(np.einsum("ij,ij->i", pnt - pp, d) / l2, 0.0, 1.0)
            q = pp + d * t[:, None]
            dist = np.linalg.norm(pnt - q, axis=1)
            take = (~inside) & (dist < best)
            best = np.where(take, dist, best)
            cp = np.where(take[:, None], q, cp)
        k = int(np.argmin(best))
        c = float(np.dot(pnt - cp[k], -fn[k]))       # sign from the closest facet
        t = tagv[vi]
        if t not in per or c < per[t]:
            per[t] = c
        if c < worst:
            worst = c
        bad += int(c < -1e-6)
    return worst, bad, len(v), per


def _vertex_tags(m, nv):
    """Tag of the first solid using each mesh vertex, so clearance can be broken
    down per feature. Diagnosing this by hand three times was the tell."""
    tag = np.asarray(m.tag)
    _, f = m.arrays()
    out = np.empty(nv, object)
    out[:] = ""
    for i, face in enumerate(f):
        for j in face:
            if out[j] == "":
                out[j] = tag[i]
    return out


def report(sol):
    m, verts, faces, fn = build(sol)
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    lo, hi = v.min(axis=0), v.max(axis=0)
    tri = v[f]
    vol = float(np.einsum("ij,ij->i", tri[:, 0],
                          np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)
    ed = np.array([[np.linalg.norm(verts[fc[(i + 1) % 3]] - verts[fc[i]])
                    for i in range(3)] for fc in faces])

    print(f"({sol['p']},{sol['r']}) phyllotactic tube: R={sol['R']:.2f} mm, "
          f"divergence {math.degrees(sol['alpha']):.3f} deg")
    print(f"  vertices {len(verts)}  edges {len(m.tag) and ''}"
          f"{int((tag=='strut').sum()//28)}  mirrors {len(faces)}")
    print(f"  MIRRORS: all edges {ed.min():.6f} .. {ed.max():.6f} mm "
          f"-> 1 shape, equilateral {sol['s']}")
    print("  triangles by feature:", "  ".join(
        f"{t}={int((tag == t).sum())}"
        for t in ("strut", "node", "pad", "base", "rim", "post")))
    print(f"  bbox x {lo[0]:7.1f}..{hi[0]:6.1f}  y {lo[1]:7.1f}..{hi[1]:6.1f}"
          f"  z {lo[2]:7.1f}..{hi[2]:6.1f}")
    print(f"  envelope {hi[0]-lo[0]:.1f} x {hi[1]-lo[1]:.1f} x {hi[2]-lo[2]:.1f} mm"
          f"   (bed 320 cubed)")
    print(f"  solids {len(m._solid_start)}  triangles {len(f)}")
    print(f"  volume {vol/1000:.1f} cm3 -> {vol*1.27/1000:.0f} g PETG solid")
    bad = m.validate()
    print(f"  manifold check: {'PASS' if not bad else f'FAIL {bad[:2]}'}")
    worst, n_bad, n_tested, per = clearance(m, verts, faces, fn)
    print(f"  clearance to the mirror surface: {n_tested} frame vertices, "
          f"{n_bad} intruding; tightest {worst:+.4f} mm")
    print(f"    pads must read exactly {GLASS} (they carry the glass); everything"
          f" else must be >= {GLASS} or it fouls the mirror:")
    for t, c in sorted(per.items(), key=lambda kv: kv[1]):
        flag = "OK" if (c >= GLASS - 1e-6) else "FOULS GLASS"
        print(f"      {t:8s} {c:+8.3f} mm   {flag}")
    return m


if __name__ == "__main__":
    sol = load()
    m = report(sol)
    n = MK.write_stl("equilateral_tube_frame.stl", m)
    print(f"\nwrote equilateral_tube_frame.stl ({n} triangles)")
