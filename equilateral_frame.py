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

# Glazing lips. Measured first: every facet is within 3.8 deg of vertical, so
# n.zhat spans only +-0.066 and gravity pulls a mirror off its pads with at most
# 6.6% of 4.2 g = 2.7 mN. The remaining 99.8% is in-plane shear, which silicone
# carries at ~0.7 kPa against ~1-2 MPa capability. So these lips are NOT needed
# for gravity -- they are insurance against a silicone-to-PETG bond releasing,
# which is the genuinely weak link.
#
# A lip must sit in FRONT of the glass to retain it, so it necessarily enters the
# bore. Corners are the only place with room: beam hits clear facet VERTICES by
# 5.54 mm but clear facet EDGES by only 0.99 mm, so edge-mid lips would be struck.
# Two lips per facet, not three -- three makes the mirror impossible to install,
# since 1.5 mm glass cannot flex under them. Two lips plus an open third corner
# lets the glass slide in.
# FIRST ATTEMPT FAILED, and geometrically so: per-facet corner lips need a stem
# reaching outboard, and with a gap-free tiling every route crosses either the
# mirror being retained or its neighbour. Measured 1544 lip vertices buried in
# neighbouring glass. The only gap-free route out of the bore is at the tube ends.
#
# So retain at the LATTICE VERTICES instead. Six mirror corners meet at each one;
# clip VOID_R off every corner and a single peg passes through the resulting hole
# and caps all six. 56 pegs instead of 184 lips, and clipped corners are wanted
# anyway -- sharp corners on 1.5 mm glass chip and the chips scatter.
#
# Sizing is set by the beam: hits clear lattice VERTICES by 5.54 mm, so the head
# radius plus the beam radius must stay inside that.
VOID_R = 3.0            # corner clipped off each mirror, measured from the vertex
PEG_R = 2.1             # stem stays inside the void: sqrt(1.5^2+2.1^2)=2.58 < 3.0
HEAD_T = 1.2            # how far a tab protrudes into the bore
TAB_R0 = 2.0            # tab runs from here to TAB_R1, measured from the vertex
TAB_R1 = 5.0            # so it overlaps the clipped mirror corner by TAB_R1-VOID_R
#
# A single ROUND head per vertex does not work: it is aligned to the AVERAGED
# vertex normal, so its rim dips into the more steeply tilted neighbours by
# HEAD_R*sin(theta) -- measured 1.33 mm into 1.5 mm glass. Shrinking it enough to
# clear leaves nothing to cap with. So each corner gets its OWN tab, lying on its
# own facet's plane, branching laterally off the central stem in the void.
#
# Tabs go only on facets the beam never touches. Beam hits clear lattice vertices
# by just 5.54 mm, and a tab reaching TAB_R1 + HEAD_T/2 = 5.6 mm would be struck.
# Only 16 of 92 facets are ever lit, so skipping those costs almost no retention
# and removes the optical risk entirely.
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
    lit = lit_facets(sol)
    for i in range(len(verts)):
        m.add_solid(*MK.tube(nodes[i], verts[i], PEG_R, nseg=7), tag="peg")
    for fi, (f, n) in enumerate(zip(faces, fn)):
        if fi in lit:
            continue
        tri = verts[f]
        for j in range(3):
            a = tri[j]
            e1 = tri[(j + 1) % 3] - a
            e2 = tri[(j + 2) % 3] - a
            b = e1 / np.linalg.norm(e1) + e2 / np.linalg.norm(e2)
            b /= np.linalg.norm(b)
            c = n * (HEAD_T / 2.0)          # sit on THIS facet's plane, not the mean
            m.add_solid(*MK.tube(a + b * TAB_R0 + c, a + b * TAB_R1 + c,
                                 HEAD_T / 2.0, nseg=8), tag="tab")
    return m, verts, faces, fn


def lit_facets(sol):
    """Facet indices the beam ever strikes, over ALL arms. Screw symmetry only, so
    arm 0 is not representative."""
    import equilateral_tube as ET
    import ray_optics as RO
    tris, nrm, _ = ET.build_helix(sol["R"], sol["alpha"], sol["beta"],
                                 sol["p"], sol["r"], sol["K"])
    out = set()
    for k in range(6):
        p0, d0 = RO.launch_from_focus(sol["theta"], az=360.0 * k / 6)
        pts, log = RO.shoot(p0, d0, tris, nrm, 6, 90)
        if len(log) < 3:
            continue
        out |= {e["facet"] for e in log}
    return out


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
    per, buried = {}, {}
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
        # A vertex with clearance strictly inside (0, GLASS) is BURIED IN A MIRROR.
        # Reporting only the per-tag minimum hides this for lips, whose hook is
        # legitimately negative -- so count it separately. The lip stems sit just
        # outside their own triangle, which is a neighbouring facet's territory.
        if 1e-6 < c < GLASS - 1e-6:
            # No glass exists within VOID_R of a lattice vertex -- that corner is
            # clipped off, which is the whole point of the peg scheme. Without this
            # exemption the check reports the stems as buried in mirrors that have
            # been cut away.
            if np.min(np.linalg.norm(verts - pnt, axis=1)) > VOID_R:
                n_b, d_b = buried.get(t, (0, 0.0))
                buried[t] = (n_b + 1, max(d_b, c))
        if c < worst:
            worst = c
        bad += int(c < -1e-6)
    return worst, bad, len(v), per, buried


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


def beam_vs_pegs(m, sol):
    """The lips are the only frame feature inside the bore, so they are the only
    ones that can be hit. Check every arm, not just arm 0 -- this tube has screw
    symmetry only, so the arms are not congruent and arm 0 is not the worst."""
    import equilateral_tube as ET
    import ray_optics as RO
    tris, nrm, _ = ET.build_helix(sol["R"], sol["alpha"], sol["beta"],
                                 sol["p"], sol["r"], sol["K"])
    hits = []
    for k in range(6):
        p0, d0 = RO.launch_from_focus(sol["theta"], az=360.0 * k / 6)
        pts, log = RO.shoot(p0, d0, tris, nrm, 6, 90)
        if len(log) < 3:
            continue
        hits.extend(pts[1:len(log) + 1])
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    lipv = np.unique(f[np.isin(tag, ["peg", "tab"])])
    if not len(lipv) or not hits:
        return
    d = np.linalg.norm(np.asarray(hits)[:, None, :] - v[lipv][None, :, :], axis=2)
    print(f"    beam vs pegs: {len(hits)} hits, {len(lipv)} peg vertices, "
          f"closest approach {d.min():.2f} mm "
          f"({'CLEAR' if d.min() > 1.5 else 'BEAM STRIKES A PEG'})")


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
        for t in ("strut", "node", "pad", "base", "rim", "post", "peg", "tab")))
    print(f"  bbox x {lo[0]:7.1f}..{hi[0]:6.1f}  y {lo[1]:7.1f}..{hi[1]:6.1f}"
          f"  z {lo[2]:7.1f}..{hi[2]:6.1f}")
    print(f"  envelope {hi[0]-lo[0]:.1f} x {hi[1]-lo[1]:.1f} x {hi[2]-lo[2]:.1f} mm"
          f"   (bed 320 cubed)")
    print(f"  solids {len(m._solid_start)}  triangles {len(f)}")
    print(f"  volume {vol/1000:.1f} cm3 -> {vol*1.27/1000:.0f} g PETG solid")
    bad = m.validate()
    print(f"  manifold check: {'PASS' if not bad else f'FAIL {bad[:2]}'}")
    worst, n_bad, n_tested, per, buried = clearance(m, verts, faces, fn)
    print(f"  clearance to the mirror surface: {n_tested} frame vertices, "
          f"{n_bad} intruding; tightest {worst:+.4f} mm")
    print(f"    pads must read exactly {GLASS} (they carry the glass); everything"
          f" else must be >= {GLASS} or it fouls the mirror:")
    for t, c in sorted(per.items(), key=lambda kv: kv[1]):
        if t in ("peg", "tab"):
            flag = "by design (tab caps a clipped mirror corner)"
        else:
            flag = "OK" if (c >= GLASS - 1e-6) else "FOULS GLASS"
        print(f"      {t:8s} {c:+8.3f} mm   {flag}")
    if buried:
        for t, (n_b, d_b) in sorted(buried.items()):
            kind = "modelling overlap" if d_b < 0.05 else "REAL INTERFERENCE"
            print(f"    buried in glass: {t} x{n_b}, deepest {d_b:.4f} mm  ({kind})")
    else:
        print(f"    no frame vertex is buried inside a mirror (0..{GLASS} band)")
    beam_vs_pegs(m, sol)
    return m


if __name__ == "__main__":
    sol = load()
    m = report(sol)
    n = MK.write_stl("equilateral_tube_frame.stl", m)
    print(f"\nwrote equilateral_tube_frame.stl ({n} triangles)")
