#!/usr/bin/env python3
"""
Parametric magnetic-vortex relief -> watertight binary STL (3D printable).

Geometry encodes the vortex texture statically:
  - core peak  : z_core(r)  = polarity * peak * exp(-(r/rc)^2)   (out-of-plane m_z)
  - curl       : spiral ridges whose handedness = chirality      (in-plane swirl)
The result is a closed manifold solid: relief top + flat bottom + cylindrical wall.

Edit the PARAMETERS block and re-run:  python3 vortex_stl.py
"""

import math, struct
import numpy as np

# ---------------- PARAMETERS (all mm) ----------------
R          = 30.0     # disk radius  (-> 60 mm diameter)
BASE_Z     = 12.0     # solid base plate height under the relief
PEAK       = 9.0      # core out-of-plane amplitude
RC_FRAC    = 0.16     # core width as fraction of R (smaller = sharper core)
SHARP      = 1.0      # core profile exponent: 1=pointy cusp, 2=rounded gaussian
RIDGE_AMP  = 1.3      # spiral ridge height
ARMS       = 6        # number of spiral arms
TWIST      = 1.3      # radial spiral rate (turns across the radius)
CHIRALITY  = +1       # +1 = counter-clockwise curl, -1 = clockwise
POLARITY   = +1       # +1 = core up (+z), -1 = core down (-z)
NR         = 140      # radial subdivisions  (resolution)
NA         = 260      # angular subdivisions (resolution)
OUT_PATH   = "magnetic_vortex.stl"
# -----------------------------------------------------

RC = RC_FRAC * R

def z_core(r):
    return POLARITY * PEAK * math.exp(-((r / RC) ** SHARP))

def z_ridge(r, th):
    # window -> 0 at center and rim so core stays clean and rim is a flat circle
    w = math.sin(math.pi * min(r / R, 1.0))
    return RIDGE_AMP * w * math.sin(ARMS * th + CHIRALITY * TWIST * 2.0 * math.pi * (r / R))

def z_top(r, th):
    return BASE_Z + z_core(r) + z_ridge(r, th)

# ---- build vertices ----
verts = []
def add(x, y, z):
    verts.append((x, y, z)); return len(verts) - 1

top_center = add(0.0, 0.0, z_top(0.0, 0.0))          # apex of relief
top_rings = []                                        # top_rings[i-1][j]
for i in range(1, NR + 1):
    r = R * i / NR
    ring = []
    for j in range(NA):
        th = 2.0 * math.pi * j / NA
        ring.append(add(r * math.cos(th), r * math.sin(th), z_top(r, th)))
    top_rings.append(ring)

bot_center = add(0.0, 0.0, 0.0)                        # flat bottom center
bot_rim = []
for j in range(NA):
    th = 2.0 * math.pi * j / NA
    bot_rim.append(add(R * math.cos(th), R * math.sin(th), 0.0))

# ---- build faces (consistent outward winding) ----
faces = []
# top center fan (normals up/out)
r1 = top_rings[0]
for j in range(NA):
    faces.append((top_center, r1[j], r1[(j + 1) % NA]))
# top quads between consecutive rings
for i in range(NR - 1):
    a, b = top_rings[i], top_rings[i + 1]
    for j in range(NA):
        jn = (j + 1) % NA
        faces.append((a[j], b[j], a[jn]))
        faces.append((a[jn], b[j], b[jn]))
# cylindrical wall: top outer ring -> bottom rim (normals point outward)
outer = top_rings[NR - 1]
for j in range(NA):
    jn = (j + 1) % NA
    faces.append((outer[j], bot_rim[j], outer[jn]))
    faces.append((outer[jn], bot_rim[j], bot_rim[jn]))
# bottom fan (normals down/out)
for j in range(NA):
    jn = (j + 1) % NA
    faces.append((bot_center, bot_rim[jn], bot_rim[j]))

V = np.asarray(verts, dtype=np.float64)

# ---- write binary STL with computed facet normals ----
with open(OUT_PATH, "wb") as f:
    f.write(b"magnetic vortex relief".ljust(80, b" "))
    f.write(struct.pack("<I", len(faces)))
    for (ia, ib, ic) in faces:
        a, b, c = V[ia], V[ib], V[ic]
        n = np.cross(b - a, c - a)
        ln = np.linalg.norm(n)
        n = n / ln if ln > 1e-12 else np.array([0.0, 0.0, 1.0])
        f.write(struct.pack("<3f", *n))
        f.write(struct.pack("<3f", *a))
        f.write(struct.pack("<3f", *b))
        f.write(struct.pack("<3f", *c))
        f.write(struct.pack("<H", 0))

print(f"wrote {OUT_PATH}: {len(verts)} verts, {len(faces)} triangles")
print(f"footprint {2*R:.0f} mm dia, height ~{BASE_Z + PEAK + RIDGE_AMP:.1f} mm, "
      f"chirality={CHIRALITY:+d}, polarity={POLARITY:+d}")
