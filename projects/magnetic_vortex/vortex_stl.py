#!/usr/bin/env python3
"""
Parametric magnetic-vortex relief -> watertight binary STL (3D printable).

This version is a THIN DOUBLE-WALLED SHELL (default 1 mm):
  - outer surface : the vortex relief  z = lift + core(r) + ridge(r,theta)
  - inner surface : the outer surface offset inward by WALL along its normal
  - the two are joined by a thin lip at the rim -> a hollow dome, open underneath

Texture encoded statically:
  core(r)  = polarity * peak * exp(-(r/rc)^sharp)   (out-of-plane m_z, the peak)
  ridge    = spiral grooves whose handedness = chirality (the in-plane curl)

Edit PARAMETERS and re-run:  python3 vortex_stl.py
"""

import math, struct
import numpy as np

# ---------------- PARAMETERS (all mm) ----------------
R          = 30.0     # disk radius  (-> 60 mm diameter)
WALL       = 0.6      # shell wall thickness (the "double wall" gap)
PEAK       = 10.8     # core out-of-plane amplitude
RC_FRAC    = 0.16     # core width as fraction of R (smaller = sharper core)
SHARP      = 1.0      # core profile exponent: 1=pointy cusp, 2=rounded gaussian
RIDGE_AMP  = 1.2      # spiral ridge height
ARMS       = 6        # number of spiral arms
TWIST      = 1.3      # radial spiral rate
CHIRALITY  = +1       # +1 = counter-clockwise curl, -1 = clockwise
POLARITY   = +1       # +1 = core up (+z), -1 = core down (-z)
NR         = 160      # radial subdivisions  (resolution)
NA         = 280      # angular subdivisions (resolution)
OUT_PATH   = "magnetic_vortex.stl"
# -----------------------------------------------------

RC = RC_FRAC * R
LIFT = WALL                      # raise dome so inner rim lands near z=0

def z_core(r):
    return POLARITY * PEAK * math.exp(-((r / RC) ** SHARP))

def z_ridge(r, th):
    w = math.sin(math.pi * min(r / R, 1.0))
    return RIDGE_AMP * w * math.sin(ARMS * th + CHIRALITY * TWIST * 2.0 * math.pi * (r / R))

def z_out(r, th):
    return LIFT + z_core(r) + z_ridge(r, th)

# ---- outer surface grid (i = 0 center .. NR rim) ----
ang = [2.0 * math.pi * j / NA for j in range(NA)]
P = np.zeros((NR + 1, NA, 3))            # outer vertex positions
for i in range(NR + 1):
    r = R * i / NR
    for j in range(NA):
        th = ang[j]
        P[i, j] = (r * math.cos(th), r * math.sin(th), z_out(r, th))
# collapse center ring to a single apex point (average)
P[0, :] = P[0, :].mean(axis=0)

# ---- smooth vertex normals on the outer surface ----
N = np.zeros_like(P)
def addN(i0, j0, i1, j1, i2, j2):
    a, b, c = P[i0, j0], P[i1, j1], P[i2, j2]
    n = np.cross(b - a, c - a)
    for (ii, jj) in ((i0, j0), (i1, j1), (i2, j2)):
        N[ii, jj] += n
for i in range(NR):
    for j in range(NA):
        jn = (j + 1) % NA
        addN(i, j, i + 1, j, i, jn)
        addN(i, jn, i + 1, j, i + 1, jn)
ln = np.linalg.norm(N, axis=2, keepdims=True); ln[ln < 1e-12] = 1.0
N /= ln
if N[:, :, 2].mean() < 0:             # ensure outward (up) normals
    N = -N

# inner surface = outer offset inward by WALL
Q = P - WALL * N

# A sharp tip cannot carry a 1 mm wall: the inner flanks would cross.
# Clamp inner radius to be non-increasing toward the apex along each spoke,
# which makes the top of the point solid while keeping a true wall elsewhere.
for j in range(NA):
    maxr = float("inf")
    for i in range(NR, -1, -1):
        x, y = Q[i, j, 0], Q[i, j, 1]
        r = math.hypot(x, y)
        if r > maxr and r > 1e-9:
            s = maxr / r
            Q[i, j, 0] = x * s; Q[i, j, 1] = y * s
        else:
            maxr = r
Q[0, :] = Q[0, :].mean(axis=0)

# ---- assemble triangles ----
faces = []
verts = []
vid = {}
def v(p):
    k = (round(p[0], 5), round(p[1], 5), round(p[2], 5))
    if k not in vid:
        vid[k] = len(verts); verts.append(p)
    return vid[k]

# outer skin: apex fan (i=0->1) then quads (i>=1)
oap = v(P[0, 0])
for j in range(NA):
    jn = (j + 1) % NA
    faces.append((oap, v(P[1, j]), v(P[1, jn])))
for i in range(1, NR):
    for j in range(NA):
        jn = (j + 1) % NA
        a, b, c, d = v(P[i, j]), v(P[i + 1, j]), v(P[i, jn]), v(P[i + 1, jn])
        faces.append((a, b, c)); faces.append((c, b, d))
# inner skin (reversed winding -> normals point into cavity)
iap = v(Q[0, 0])
for j in range(NA):
    jn = (j + 1) % NA
    faces.append((iap, v(Q[1, jn]), v(Q[1, j])))
for i in range(1, NR):
    for j in range(NA):
        jn = (j + 1) % NA
        a, b, c, d = v(Q[i, j]), v(Q[i + 1, j]), v(Q[i, jn]), v(Q[i + 1, jn])
        faces.append((a, c, b)); faces.append((c, d, b))
# rim lip joining outer edge to inner edge (closes the shell)
for j in range(NA):
    jn = (j + 1) % NA
    o0, o1 = v(P[NR, j]), v(P[NR, jn])
    i0, i1 = v(Q[NR, j]), v(Q[NR, jn])
    faces.append((o0, i0, o1)); faces.append((o1, i0, i1))

V = np.asarray(verts, dtype=np.float64)

# ---- write binary STL ----
with open(OUT_PATH, "wb") as f:
    f.write(b"magnetic vortex shell".ljust(80, b" "))
    f.write(struct.pack("<I", len(faces)))
    for (ia, ib, ic) in faces:
        a, b, c = V[ia], V[ib], V[ic]
        n = np.cross(b - a, c - a); m = np.linalg.norm(n)
        n = n / m if m > 1e-12 else np.array([0.0, 0.0, 1.0])
        f.write(struct.pack("<3f", *n))
        f.write(struct.pack("<3f", *a)); f.write(struct.pack("<3f", *b)); f.write(struct.pack("<3f", *c))
        f.write(struct.pack("<H", 0))

# ---- validate watertight + report min vertical clearance ----
from collections import defaultdict
edges = defaultdict(int)
for (a, b, c) in faces:
    for e in ((a, b), (b, c), (c, a)):
        edges[tuple(sorted(e))] += 1
bad = sum(1 for x in edges.values() if x != 2)
gap = (P - Q)[:, :, 2]
print(f"wrote {OUT_PATH}: {len(verts)} verts, {len(faces)} tris")
print(f"watertight: euler V-E+F = {len(verts)-len(edges)+len(faces)}, non-manifold edges = {bad}")
print(f"wall = {WALL} mm, dia {2*R:.0f} mm, height ~{PEAK+LIFT:.1f} mm, "
      f"chirality={CHIRALITY:+d}, polarity={POLARITY:+d}")
