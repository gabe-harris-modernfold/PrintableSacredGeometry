#!/usr/bin/env python3
"""Independent wall-thickness check on the stellated Platonic shells.

stellated_platonics.hollow() argues the wall is exactly WALL because the inner
surface is a uniform scale of the outer one and every flank plane sits at the same
distance from the centre. That argument is about the *flanks*. It says nothing
directly about the concave valleys along the core solid's edges, where two flanks
meet at a reflex angle and a naive offset would normally thin out.

So this measures it the dumb way instead: scatter points over the outer surface and
brute-force the distance from each to every inner triangle. Small meshes, so exact
point-triangle distance over all pairs is affordable and there is nothing to trust.

Run:  python verify_stellated_platonics.py
"""

import sys

import numpy as np
import trimesh

import stellated_platonics as SP

SAMPLES_PER_TRI = 256


def sample_triangles(tri, n):
    """Roughly uniform barycentric samples, plus the corners and edge midpoints --
    the corners are where a thin spot would hide."""
    u = np.random.default_rng(7).random((len(tri), n, 2))
    flip = u.sum(axis=2) > 1.0
    u[flip] = 1.0 - u[flip]
    w = np.stack([1.0 - u[..., 0] - u[..., 1], u[..., 0], u[..., 1]], axis=-1)
    pts = np.einsum("tsk,tkj->tsj", w, tri).reshape(-1, 3)
    corners = tri.reshape(-1, 3)
    mids = np.concatenate([0.5 * (tri[:, i] + tri[:, (i + 1) % 3]) for i in range(3)])
    return np.vstack([pts, corners, mids])


def point_tri_distance(P, tri, chunk=4096):
    """Min distance from each point in P to the triangle soup `tri` (Ericson)."""
    A, B, C = tri[:, 0], tri[:, 1], tri[:, 2]
    ab, ac = B - A, C - A
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    den = d00 * d11 - d01 * d01

    out = np.empty(len(P))
    for lo in range(0, len(P), chunk):
        p = P[lo:lo + chunk][:, None, :]
        ap = p - A[None]
        d20 = np.einsum("pij,ij->pi", ap, ab)
        d21 = np.einsum("pij,ij->pi", ap, ac)
        v = (d11 * d20 - d01 * d21) / den
        w = (d00 * d21 - d01 * d20) / den
        u = 1.0 - v - w

        # clamp the barycentric coordinates back onto the triangle
        inside = (u >= 0) & (v >= 0) & (w >= 0)
        v = np.clip(v, 0.0, 1.0)
        w = np.clip(w, 0.0, 1.0)
        s = v + w
        over = s > 1.0
        v = np.where(over, v / np.where(over, s, 1.0), v)
        w = np.where(over, w / np.where(over, s, 1.0), w)
        q = A[None] + v[..., None] * ab[None] + w[..., None] * ac[None]
        dist = np.linalg.norm(p - q, axis=2)

        # a point projecting inside the triangle still needs the plane distance
        n = np.cross(ab, ac)
        n = n / np.linalg.norm(n, axis=1, keepdims=True)
        plane = np.abs(np.einsum("pij,ij->pi", ap, n))
        dist = np.where(inside, np.minimum(dist, plane), dist)
        out[lo:lo + chunk] = dist.min(axis=1)
    return out


def main():
    ok = True
    print(f"target wall = {SP.WALL:.2f} mm   "
          f"({SAMPLES_PER_TRI} samples/triangle + corners + edge midpoints)")
    print("-" * 84)
    for name in SP.SOLIDS:
        rec = SP.build(name)
        m, sh = SP.as_mesh(rec, SP.SIZE)

        half = len(rec["tris"])
        V, F = np.asarray(m.vertices), np.asarray(m.faces)
        outer, inner = V[F[:half]], V[F[half:]]

        d = point_tri_distance(sample_triangles(outer, SAMPLES_PER_TRI), inner)
        lo, hi = float(d.min()), float(d.max())
        # thin spots must not exist; thick spots are fine and expected at the tips
        good = lo >= SP.WALL - 1e-3
        ok &= good and bool(m.is_watertight)
        print(f"{name:<13} min {lo:7.3f} mm   max {hi:7.3f} mm   "
              f"predicted {sh['wall_min']:.2f}-{sh['wall_max']:.2f}   "
              f"watertight {str(m.is_watertight):<5} "
              f"[{'PASS' if good else 'FAIL'}]")
    print("-" * 84)
    print("RESULT:", "walls hold" if ok else "THIN SPOT FOUND")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
