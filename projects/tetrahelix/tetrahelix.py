#!/usr/bin/env python3
"""
Fuller's ten-cycle tetrahelix (Synergetics 933.00) as an open 3-sided-beam
scaffold, in two versions: the natural one that misses closure by 5.69 deg, and
a closed one that comes back exactly every ten units.

THE TEN-CYCLE
-------------
Stacking REGULAR tetrahedra face-to-face forces a twist of exactly
theta = arccos(-2/3) = 131.8103149 deg per tetrahedron. So the "ten" is not a
count of tetrahedra -- ten of those sweep 3.66 turns. It counts TRIPLE-BONDED
tetrahedra (Fuller's wording): steps along one of the three helical edge-strands,
each joining k to k+3. One strand step advances 3*theta mod 360 = 35.4309447 deg,
and ten of them make 354.309447 deg -- one turn less the published 5.690553 deg.

CLOSING IT (--closed, the default)
----------------------------------
To close exactly, ten strand steps must make 360, i.e. 3*theta = 36 mod 360, so
theta = 12, 132 or 252 deg. 132 is 0.1897 deg off the natural twist, so it is the
one that barely disturbs the solid. Fixing theta then leaves only R and H free
against three step lengths, so the cells cannot stay regular -- but they come
much closer than expected. Choosing H^2/R^2 = (cos3t - cos1t)/4 makes d1 and d3
EXACTLY equal, leaving d2 short by 0.2781%:

    d1 = d3 = 30.000 mm      d2 = 29.917 mm      (0.083 mm apart)

Two beam lengths, and the difference is under one layer height. R and H shift by
0.06% and 0.17% from regular. In exchange:

  * ten strand steps = 360.000000 deg, dead closed;
  * 30 vertex steps = exactly 11 turns, so v[k+30] sits exactly above v[k] --
    the column becomes vertically periodic, i.e. genuinely stackable, which the
    natural helix never is at any length.

--regular gives back the true Boerdijk-Coxeter solid: one beam length, regular
cells, and the 5.69 deg deficit that Fuller hung the DNA argument on. Use that
one if the irrationality is the point; use --closed if the closure is.

WHY A SCAFFOLD
--------------
Consecutive tetrahedra share a whole FACE, so the union of all their edges
collapses to three step families off one vertex helix -- k to k+1, k+2, k+3, and
nothing else. 3V-6 beams total, and the column is ~89% air.

Run:  python tetrahelix.py [--closed | --regular] [--tencycles 1] [--edge 30]
"""

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import mesh_kit as MK

THETA_REG = math.acos(-2.0 / 3.0)       # 131.8103149 deg -- regular cells
THETA_CLS = math.radians(132.0)         # exact ten-unit closure
TEN_CYCLE = 30          # vertex steps in one ten-cycle = 10 strand steps of 3

EDGE = 30.0             # nominal beam length (d1 = d3). 30 steps x ~0.316 x edge
                        # is the whole height, so this is the only knob that can
                        # bust the bed: 32.8 mm is the ceiling for one ten-cycle.
BEAM_R = 3.3            # circumradius of the triangular section (side 5.7 mm,
                        # wall-to-axis 1.65 mm) -- thin enough to see through,
                        # thick enough that a 0.4 mm nozzle lays 4 perimeters.
NODE_R = 4.6            # octahedral ball at each vertex. Must exceed BEAM_R or
                        # the six beams meeting there only kiss instead of fusing.

# The base is OFF by default: a 46 mm-radius tripod under a 31 mm-wide helix is a
# bigger object than the sculpture it carries, and it reads as a giant tetrahedron
# stuck on the bottom rather than as part of the stack. Bare, the helix lands on
# its single lowest node -- fine for a slicer with a brim, and the piece is only
# 52 g, but it is a point contact and will want holding while layer one goes down.
BASE_H = 24.0
BASE_R = 46.0
BASE_LEGS = 2


def geometry(edge=EDGE, closed=True):
    """(theta, R, H) for the requested helix, scaled so the d1 beam is `edge`.

    Regular: closed forms R = 3*sqrt3/10, H = 1/sqrt10 make all three step
    lengths equal -- that is what pins theta to arccos(-2/3) in the first place.
    Closed: theta is pinned to 132 deg instead, which spends the equal-length
    freedom, so H is set to the one value that keeps d1 == d3."""
    if not closed:
        return THETA_REG, 3.0 * math.sqrt(3.0) / 10.0 * edge, edge / math.sqrt(10.0)
    t = THETA_CLS
    c1, c3 = math.cos(t), math.cos(3 * t)
    x = (c3 - c1) / 4.0                         # H^2/R^2 that gives d1 == d3
    R = edge / math.sqrt(2.0 * (1.0 - c1) + x)
    return t, R, R * math.sqrt(x)


def helix(n_vert, edge=EDGE, closed=True):
    t, R, H = geometry(edge, closed)
    k = np.arange(n_vert)
    return np.stack([R * np.cos(k * t), R * np.sin(k * t), k * H], axis=1)


def topology(n_vert):
    """Edges and tetrahedra of the BC helix.

    Consecutive tetrahedra share a whole face, so {k,k+1,k+2,k+3} over all k
    yields only the three step families. The k->k+3 family splits into three
    interleaved strands by k mod 3; strand 0 is the one the ten-cycle is counted
    along, and it is tagged separately so it can be seen in the print."""
    edges, strand = [], []
    for m in (1, 2, 3):
        for k in range(n_vert - m):
            (strand if (m == 3 and k % 3 == 0) else edges).append((k, k + m))
    tets = [(k, k + 1, k + 2, k + 3) for k in range(n_vert - 3)]
    assert len(edges) + len(strand) == 3 * n_vert - 6, "edge count must be 3V-6"
    return edges, strand, tets


def beam(p0, p1, r, ref):
    """Triangular prism from p0 to p1, one apex aimed along `ref`.

    Written out rather than calling MK.tube(nseg=3) because tube() picks its
    cross-section basis from a fixed world vector, which leaves the three flats
    rolling at random along the helix -- the lattice then reads as pipe, not as
    faceted stock. Winding matches tube(): sides then two caps."""
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    ax = p1 - p0
    ax /= np.linalg.norm(ax)
    u = np.asarray(ref, float)
    u = u - ax * float(u @ ax)                  # project the aim off the axis
    if np.linalg.norm(u) < 1e-9:                # ref parallel to the beam
        u = np.cross(ax, [0.0, 0.0, 1.0])
        if np.linalg.norm(u) < 1e-9:
            u = np.cross(ax, [1.0, 0.0, 0.0])
    u /= np.linalg.norm(u)
    w = np.cross(ax, u)
    ring = [u * math.cos(2 * math.pi * k / 3) + w * math.sin(2 * math.pi * k / 3)
            for k in range(3)]
    verts = [p0 + r * d for d in ring] + [p1 + r * d for d in ring] + [p0, p1]
    faces = []
    for k in range(3):
        j = (k + 1) % 3
        faces += [(k, j, 3 + j), (k, 3 + j, 3 + k), (6, j, k), (7, 3 + k, 3 + j)]
    return verts, faces


def octa(c, r):
    """Octahedral node ball -- 8 triangles, watertight, and the octahedron is the
    tetrahelix's own dual cell, so the joints look native to the stack."""
    c = np.asarray(c, float)
    verts = [c + np.array(d) * r for d in
             ((1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))]
    faces = [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4),
             (1, 0, 5), (2, 1, 5), (3, 2, 5), (0, 3, 5)]
    return verts, faces


def radial(p):
    """Outward horizontal direction at p, the aim used for every beam apex."""
    u = np.array([p[0], p[1], 0.0])
    n = np.linalg.norm(u)
    return u / n if n > 1e-9 else np.array([1.0, 0.0, 0.0])


def base_frame(verts, base_r=BASE_R, n_legs=BASE_LEGS):
    """Tripod: a ground triangle plus legs up to the lowest helix vertices.

    The three lowest vertices sit ~132 deg apart in azimuth -- close to a tripod
    but not one, so the feet go on a true 120 deg triangle and each foot reaches
    for its nearest low vertices instead. That keeps the footprint stable without
    pretending the helix has a symmetry it does not have."""
    a0 = math.atan2(verts[0][1], verts[0][0])
    feet = [np.array([base_r * math.cos(a0 + 2 * math.pi * i / 3),
                      base_r * math.sin(a0 + 2 * math.pi * i / 3), 0.0])
            for i in range(3)]
    beams = [(feet[i], feet[(i + 1) % 3]) for i in range(3)]
    low = verts[:6]
    for f in feet:
        for j in np.argsort(np.linalg.norm(low - f, axis=1))[:n_legs]:
            beams.append((f, low[j]))
    return feet, beams


def steps_for(tencycles=1.0, tets=None, turns=None, closed=True):
    """Vertex steps. tencycles and tets are exact; turns must round."""
    if tets is not None:
        return max(3, int(tets) + 2)
    if turns is not None:
        t = THETA_CLS if closed else THETA_REG
        return max(3, int(round(turns * 2.0 * math.pi / t)))
    return max(3, int(round(tencycles * TEN_CYCLE)))


def build(tencycles=1.0, edge=EDGE, with_base=False, tets=None, turns=None,
          closed=True):
    steps = steps_for(tencycles, tets, turns, closed)
    n_vert = steps + 1
    verts = helix(n_vert, edge, closed)
    edges, strand, cells = topology(n_vert)

    if with_base:
        verts = verts + np.array([0.0, 0.0, BASE_H])

    m = MK.Mesh()
    for tag, group in (("strut", edges), ("strand", strand)):
        for a, b in group:
            mid = 0.5 * (verts[a] + verts[b])
            m.add_solid(*beam(verts[a], verts[b], BEAM_R, radial(mid)), tag=tag)
    for p in verts:
        m.add_solid(*octa(p, NODE_R), tag="node")

    feet = []
    if with_base:
        feet, bb = base_frame(verts)
        for p, q in bb:
            m.add_solid(*beam(p, q, BEAM_R, radial(0.5 * (p + q))), tag="pad")
        for f in feet:
            m.add_solid(*octa(f, NODE_R), tag="pad")
    return m, verts, edges + strand, strand, cells, feet, steps


def _verify(edge=EDGE):
    """Every claim in the docstring is asserted rather than trusted."""
    # Regular: all three step lengths equal, cells exactly regular.
    v = helix(8, edge, closed=False)
    for m in (1, 2, 3):
        assert np.allclose(np.linalg.norm(v[m:] - v[:-m], axis=1), edge, atol=1e-9)
    t = v[:4]
    vol = abs(np.dot(t[1] - t[0], np.cross(t[2] - t[0], t[3] - t[0]))) / 6.0
    assert abs(vol - edge ** 3 / (6 * math.sqrt(2))) < 1e-6, vol
    step = math.degrees(3 * THETA_REG) % 360.0
    assert abs(360.0 - 10 * step - 5.690553) < 1e-6, step

    # Closed: strand step exactly 36 deg, ten of them exactly 360, d1 == d3,
    # and the cells non-degenerate.
    assert abs((math.degrees(3 * THETA_CLS) % 360.0) - 36.0) < 1e-12
    v = helix(TEN_CYCLE + 1, edge, closed=True)
    d = [np.linalg.norm(v[m:] - v[:-m], axis=1) for m in (1, 2, 3)]
    for q in d:
        assert np.ptp(q) < 1e-9, "each step family must still be one length"
    assert abs(d[0][0] - edge) < 1e-9 and abs(d[2][0] - edge) < 1e-9, "d1 == d3"
    assert d[1][0] < d[0][0], "d2 is the short one"
    az = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
    assert abs((az[TEN_CYCLE] - az[0] + 180) % 360 - 180) < 1e-9, "must close"
    assert np.allclose(v[TEN_CYCLE, :2], v[0, :2], atol=1e-9), "v30 above v0"
    t = v[:4]
    assert abs(np.dot(t[1] - t[0], np.cross(t[2] - t[0], t[3] - t[0]))) > 1e-6


def report(edge, with_base, tencycles=1.0, tets=None, turns=None, closed=True):
    _verify(edge)
    m, verts, edges, strand, cells, feet, steps = build(
        tencycles, edge, with_base, tets, turns, closed)
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    lo, hi = v.min(axis=0), v.max(axis=0)
    tri = v[f]
    vol = float(np.einsum("ij,ij->i", tri[:, 0],
                          np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)
    t, R, H = geometry(edge, closed)
    step_deg = math.degrees(3 * t) % 360.0

    print(f"Fuller ten-cycle tetrahelix (Synergetics 933.00) -- "
          f"{'CLOSED (theta forced to 132 deg)' if closed else 'REGULAR (true BC helix)'}")
    print(f"  twist {math.degrees(t):.7f} deg per tetrahedron"
          + ("" if closed else "  (= arccos(-2/3))"))
    if closed:
        print(f"    natural twist is {math.degrees(THETA_REG):.7f}; forced by "
              f"{math.degrees(t - THETA_REG):+.4f} deg to buy exact closure")
    print(f"  one k->k+3 edge strand advances {step_deg:.7f} deg per step")
    print(f"    10 strand steps = {10*step_deg:.6f} deg = 1 turn "
          + (f"EXACTLY  <- THE TEN-CYCLE, CLOSED" if closed
             else f"less {360-10*step_deg:.6f} deg  <- THE TEN-CYCLE"))
    print(f"    {360/step_deg:.4f} strand units per turn  vs  B-DNA 10.4-10.5 bp")

    # Measured back off the built mesh, not restated from the constants above.
    az = np.degrees(np.arctan2(verts[:, 1], verts[:, 0]))
    clo = [(j, (az[j] - az[0] + 180) % 360 - 180) for j in range(3, len(verts), 3)]
    print("  strand closure measured on the model (azimuth vs vertex 0):")
    for i in range(0, min(len(clo), 20), 10):
        print("   ", "  ".join(f"v{j}:{d:+7.2f}" for j, d in clo[i:i + 10]))
    if len(verts) > TEN_CYCLE:
        dxy = float(np.linalg.norm(verts[TEN_CYCLE, :2] - verts[0, :2]))
        print(f"  v{TEN_CYCLE} sits {dxy:.6f} mm off the axis-line through v0 "
              f"-> {'vertically periodic, modules stack' if dxy < 1e-6 else 'no repeat'}")

    L = np.array([np.linalg.norm(verts[b] - verts[a]) for a, b in edges])
    uniq, cnt = np.unique(np.round(L, 6), return_counts=True)
    print(f"  BEAMS: {len(uniq)} distinct length(s), 3-sided, section side "
          f"{BEAM_R*math.sqrt(3):.2f} mm")
    for q, c in zip(uniq, cnt):
        print(f"      {q:9.4f} mm  x{c:3d}   ({100*(q-uniq.max())/uniq.max():+.4f}%"
              f" vs longest)")
    print(f"  helix radius {R:.4f} mm, rise {H:.4f} mm per step")
    print(f"  {steps} vertex steps = {steps/TEN_CYCLE:.3f} ten-cycles"
          f" = {steps*math.degrees(t)/360:.4f} raw turns")
    print(f"  vertices {len(verts)}  face-to-face additions {steps}  "
          f"closed tetrahedral cells {len(cells)}")
    print(f"  beams {len(edges)} (= 3V-6), of which {len(strand)} are the marked "
          f"strand")
    print("  triangles by feature:", "  ".join(
        f"{q}={int((tag == q).sum())}" for q in ("strut", "strand", "node", "pad")
        if (tag == q).any()))
    print(f"  bbox x {lo[0]:7.1f}..{hi[0]:6.1f}  y {lo[1]:7.1f}..{hi[1]:6.1f}"
          f"  z {lo[2]:7.1f}..{hi[2]:6.1f}")
    ext = hi - lo
    fits = "FITS" if (ext <= 320.0 + 1e-9).all() else "OVER BED"
    print(f"  envelope {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm   "
          f"(bed 320 cubed: {fits})")
    print(f"  solids {len(m._solid_start)}  triangles {len(f)}")
    rad = float(np.linalg.norm(v[:, :2], axis=1).max())
    hull = math.pi * rad ** 2 * ext[2]
    print(f"  volume {vol/1000:.1f} cm3 -> {vol*1.27/1000:.0f} g PETG "
          f"({100.0*vol/hull:.1f}% of the enclosing cylinder -- "
          f"{100.0-100.0*vol/hull:.1f}% open)")
    if not with_base:
        foot = verts[int(np.argmin(verts[:, 2]))]
        print(f"  no base: stands on the single lowest node at "
              f"r={np.linalg.norm(foot[:2]):.1f} mm -- point contact, wants a brim")
    bad = m.validate()
    print(f"  manifold check: {'PASS' if not bad else f'FAIL {bad[:2]}'} "
          f"({len(m._solid_start)} solids, each edge-manifold; they overlap at "
          f"the nodes by design, so the union self-intersects)")
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regular", action="store_true",
                    help="true Boerdijk-Coxeter solid: regular cells, one beam "
                         "length, 5.69 deg deficit (default is the closed helix)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--tencycles", type=float, default=1.0,
                   help="ten-cycles, 30 vertex steps each (default 1)")
    g.add_argument("--tets", type=int,
                   help="closed tetrahedral cells instead -- exact")
    g.add_argument("--turns", type=float,
                   help="raw 360 deg turns of the VERTEX helix; rounds")
    ap.add_argument("--edge", type=float, default=EDGE,
                    help="nominal beam length d1 = d3, mm (default 30)")
    ap.add_argument("--base", action="store_true",
                    help="add the tripod base (off by default)")
    ap.add_argument("--out", default="tetrahelix.stl")
    a = ap.parse_args()
    m = report(a.edge, a.base, a.tencycles, a.tets, a.turns, not a.regular)
    n = MK.write_stl(a.out, m)
    print(f"\nwrote {a.out} ({n} triangles)")


if __name__ == "__main__":
    main()
