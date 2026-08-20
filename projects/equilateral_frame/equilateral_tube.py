#!/usr/bin/env python3
"""
Equilateral 50.8 mm mirrors only. Two questions:

  1. Does the obvious answer work?  An exact equilateral tiling of a constant-
     radius tube forces e1 = e2, hence twist = 180/N: the UNIFORM ANTIPRISM. But
     an antiprism's up-triangle puts its apex at theta = pi/N, exactly on the
     perpendicular bisector of its base, so every facet is mirror-symmetric about
     a meridional plane. Its normal lies in that plane, n.phihat = 0, L = r*v_phi
     is conserved at every bounce, and the ray can never reach the axis.
     -> checked in antiprism_check().

  2. Is there any other equilateral tiling of a tube?  Yes: drop the planar-ring
     topology. Put the vertices on a single HELIX, v_k = (R cos k*alpha,
     R sin k*alpha, k*beta), and let each vertex join k+-p, k+-r, k+-(p+r). Every
     triangle is (k, k+p, k+p+r) or (k, k+r, k+p+r), whose edges are the chords
     D(p), D(r), D(p+r). Demand all three equal s and you get an equilateral
     triangulated tube that is CHIRAL -- no meridional mirror plane, so
     n.phihat != 0 and the de-swirl is back on.

     (p, r) are the parastichy numbers. Fibonacci pairs give the phyllotactic
     tubes; the antiprism is the degenerate achiral member of the same family.

Run:  python equilateral_tube.py
"""

import math
import os
import sys

import numpy as np
from scipy.optimize import fsolve

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import ray_optics as V

S = 50.8            # required mirror edge, all three of them


# ------------------------------------------------- 1. the antiprism dead end --


def antiprism(N, s=S, M=6):
    """Uniform antiprism stack: the ONLY equilateral tiling with planar rings."""
    R = s / (2.0 * math.sin(math.pi / N))
    h = s * math.sqrt(1.0 - 1.0 / (4.0 * math.cos(math.pi / (2 * N)) ** 2))
    rings = []
    for j in range(M + 1):
        a = math.pi / N * j + 2.0 * np.pi * np.arange(N) / N
        rings.append(np.stack([R * np.cos(a), R * np.sin(a),
                               np.full(N, (M - j) * h)], axis=1))
    return rings, R, h


def antiprism_check(N=6, M=6):
    rings, R, h = antiprism(N, M=M)
    tris, nrm = V.facets(rings)
    ed = [sorted(float(np.linalg.norm(t[(i + 1) % 3] - t[i])) for i in range(3))
          for t in tris]
    c = tris.mean(axis=1)
    phat = np.stack([-c[:, 1], c[:, 0], np.zeros(len(c))], 1)
    phat /= np.linalg.norm(phat, axis=1, keepdims=True)
    azim = np.abs(np.einsum("ij,ij->i", nrm, phat))
    print(f"UNIFORM ANTIPRISM  N={N}  R={R:.3f}  band h={h:.3f}  "
          f"across corners {2*R:.1f} mm")
    print(f"  edge lengths: min {np.min(ed):.6f}  max {np.max(ed):.6f} mm "
          f"(target {S})")
    print(f"  |n.phihat| over all {len(nrm)} facets: max {azim.max():.3e}")
    print(f"  -> azimuthal normal component is zero, so L is conserved exactly")

    # trace and watch L
    p0, d0 = V.launch_from_focus(70.0)
    pts, log = V.shoot(p0, d0, tris, nrm, N, 60)
    Ls = [s2["L"] for s2 in log]
    print(f"  trace from the axis: {len(log)} bounces, "
          f"L ranges {min(Ls):+.4f} .. {max(Ls):+.4f}")
    print(f"  a ray launched ON the axis has L=0 and STAYS meridional -> it can "
          f"reach the axis, but it never acquires swirl: no vortex.")
    # and a swirling ray is locked out of the middle
    p1 = np.array([R * 0.9, 0.0, h * M * 0.5])
    d1 = np.array([-0.30, 0.80, -0.52])
    d1 /= np.linalg.norm(d1)
    pts1, log1 = V.shoot(p1, d1, tris, nrm, N, 60)
    L1 = np.array([s2["L"] for s2 in log1])
    r1 = np.array([s2["r"] for s2 in log1])
    print(f"  a SWIRLING ray: |L| = {abs(L1).min():.3f}..{abs(L1).max():.3f} "
          f"(constant), r never below {r1.min():.2f} mm -> caustic locks it out")
    return tris, nrm


# --------------------------------------- 2. chiral helical equilateral tubes --


def chord(R, alpha, beta, j):
    return math.hypot(2.0 * R * math.sin(j * alpha / 2.0), j * beta)


def solve_helix(p, r, s=S, guess=None):
    """Find (R, alpha, beta) making chords D(p) = D(r) = D(p+r) = s."""
    q = p + r

    def F(x):
        R, alpha, beta = x
        return [chord(R, alpha, beta, p) - s,
                chord(R, alpha, beta, r) - s,
                chord(R, alpha, beta, q) - s]

    best = None
    for a0 in np.linspace(0.15, 3.0, 40):
        for R0 in np.linspace(0.4 * s, 4.0 * s, 20):
            b0 = s / max(p, 1) * 0.3
            try:
                x, info, ok, _ = fsolve(F, [R0, a0, b0], full_output=True)
            except Exception:
                continue
            if ok != 1:
                continue
            R, alpha, beta = x
            if R < 0.2 * s or beta <= 1e-6 or not (0.02 < alpha < 2 * math.pi):
                continue
            if max(abs(v) for v in F(x)) > 1e-7:
                continue
            if best is None or R > best[0]:
                best = (R, alpha, beta)
    return best


def build_helix(R, alpha, beta, p, r, K):
    """Vertices on one helix; two triangle families, all equilateral."""
    k = np.arange(K)
    v = np.stack([R * np.cos(k * alpha), R * np.sin(k * alpha), k * beta], 1)
    q = p + r
    tris = []
    for i in range(K - q):
        tris.append([v[i], v[i + p], v[i + q]])
        tris.append([v[i], v[i + r], v[i + q]])
    t = np.asarray(tris, float)
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    c = t.mean(axis=1)
    n[np.einsum("ij,ij->i", n[:, :2], c[:, :2]) > 0] *= -1.0
    return t, n, v


if __name__ == "__main__":
    antiprism_check(6)
    print()
    print("=" * 72)
    print("CHIRAL HELICAL EQUILATERAL TUBES  (vertices on one helix)")
    print(f"{'p,r':>7} {'R mm':>8} {'alpha deg':>10} {'beta mm':>9} "
          f"{'edge err':>10} {'|n.phi| max':>12} {'bore r':>8}")
    for (p, r) in [(1, 2), (2, 3), (1, 3), (3, 5), (2, 5), (5, 8), (3, 7),
                   (4, 7), (5, 7)]:
        sol = solve_helix(p, r)
        if sol is None:
            print(f"{p},{r:<5} {'no solution':>8}")
            continue
        R, alpha, beta = sol
        t, n, v = build_helix(R, alpha, beta, p, r, 60)
        ed = np.array([[np.linalg.norm(tt[(i + 1) % 3] - tt[i]) for i in range(3)]
                       for tt in t])
        c = t.mean(axis=1)
        ph = np.stack([-c[:, 1], c[:, 0], np.zeros(len(c))], 1)
        ph /= np.linalg.norm(ph, axis=1, keepdims=True)
        azim = np.abs(np.einsum("ij,ij->i", n, ph)).max()
        # inscribed radius of the faceted bore = min distance from axis to a facet
        bore = min(float(np.linalg.norm(np.cross(tt[1] - tt[0], tt[2] - tt[0])) and
                         abs(np.dot(n[i], tt[0]))) for i, tt in enumerate(t))
        print(f"{p},{r:<5} {R:8.2f} {math.degrees(alpha):10.3f} {beta:9.3f} "
              f"{abs(ed - S).max():10.2e} {azim:12.4f} {bore:8.2f}")
