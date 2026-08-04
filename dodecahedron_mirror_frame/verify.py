"""
Pin-fit test.

Seat a real pin in tile 0's groove and intersect it against every tile in the
assembled ball. If two mating tiles' half-grooves form one true bore the pin
touches no solid material anywhere. If the groove were not self-mating, the
half sticking across the mitre plane would collide with the neighbour.

Also reports how much of the pin is actually captured by the neighbour, so a
groove that lines up but is too shallow to grip would still be caught.
"""
import sys, os, math
import numpy as np
import trimesh
from trimesh.creation import cylinder

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dodeca_frame import (build_oriented, build_tile, face_rotations,
                          edge_frame, prism, r_i, mirror_inr, MIRROR_THK,
                          PIN_D, PIN_LEN, GROOVE_Z, ENGINE)

out = sys.argv[1]
solid, normals, fverts = build_oriented()
rot_deg = math.degrees(math.atan2(fverts[0][1], fverts[0][0]))
tile = build_tile(solid, fverts, rot_deg)

body = tile.copy(); body.apply_translation([0, 0, r_i])
rots = face_rotations(solid, normals)
frames = [body.copy() for _ in rots]
for f, (_, T) in zip(frames, rots):
    f.apply_transform(T)

mid, e = edge_frame(fverts, 0)
P = mid * (1.0 - GROOVE_Z / r_i)
pin = cylinder(radius=PIN_D / 2.0 - 0.05, sections=64,
               segment=[P - e * PIN_LEN / 2, P + e * PIN_LEN / 2])
print(f"pin        : d{PIN_D} x {PIN_LEN:.0f} mm, volume {pin.volume:.0f} mm3")

total = 0.0
for i, f in enumerate(frames):
    v = abs(trimesh.boolean.intersection([pin, f], engine=ENGINE).volume)
    if v > 1.0:
        print(f"  COLLISION with tile {i}: {v:.1f} mm3")
    total += v
print(f"interference: {total:.3f} mm3 -> "
      f"{'PASS, grooves are collinear' if total < 5 else 'FAIL'}")

# capture: how much pin volume sits inside each tile's groove void?
# grow the pin radially and see which tiles now hit it -- the two that share
# this edge must both be involved, and roughly equally.
fat = cylinder(radius=PIN_D / 2.0 + 0.8, sections=64,
               segment=[P - e * PIN_LEN / 2, P + e * PIN_LEN / 2])
caps = []
for i, f in enumerate(frames):
    v = abs(trimesh.boolean.intersection([fat, f], engine=ENGINE).volume)
    if v > 5.0:
        caps.append((i, v))
print(f"tiles gripping the pin: {[c[0] for c in caps]} "
      f"(volumes {[round(c[1]) for c in caps]} mm3)")
print("verdict    :",
      "PASS - exactly 2 tiles, balanced" if len(caps) == 2
      and abs(caps[0][1] - caps[1][1]) / max(caps[0][1], 1) < 0.1
      else "CHECK")


# ---- renders -------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def draw(ax, mesh, azim, elev, base=(0.56, 0.61, 0.69)):
    A, E = math.radians(azim), math.radians(elev)
    fwd = np.array([math.cos(E)*math.cos(A), math.cos(E)*math.sin(A), math.sin(E)])
    up = np.array([0, 0, 1.0]); rt = np.cross(up, fwd); rt /= np.linalg.norm(rt)
    up = np.cross(fwd, rt); M = np.vstack([rt, up, fwd])
    T = mesh.vertices[mesh.faces] @ M.T
    N = mesh.face_normals @ M.T
    keep = N[:, 2] < 0
    T, N = T[keep], N[keep]
    sh = np.clip(np.abs(N @ (M @ np.array([.35, .5, .9]))), 0, 1) * .72 + .28
    o = np.argsort(T[:, :, 2].mean(axis=1))[::-1]
    ax.add_collection(PolyCollection(
        T[o][:, :, :2],
        facecolors=np.clip(np.array(base)[None, :] * sh[o][:, None], 0, 1),
        edgecolors="none", antialiased=False))
    ax.set_xlim(T[:, :, 0].min()-4, T[:, :, 0].max()+4)
    ax.set_ylim(T[:, :, 1].min()-4, T[:, :, 1].max()+4)
    ax.set_aspect("equal"); ax.axis("off")


flat = trimesh.load(os.path.join(out, "dodeca-mirror-tile.stl"))
asm = trimesh.load(os.path.join(out, "dodeca-mirror-assembly.stl"))
frames_only = trimesh.load(os.path.join(out, "_frames_only.stl"))
pair = trimesh.util.concatenate([frames[0], frames[1], pin])

fig, ax = plt.subplots(2, 3, figsize=(16.5, 10.6), facecolor="white")
for a, (m, az, el, t, c) in zip(ax.ravel(), [
    (flat, 40, 38, "one tile as printed - pocket up, no supports needed", (.55, .60, .68)),
    (flat, 40, -34, "underside - flat on the bed, groove visible in the mitre", (.55, .60, .68)),
    (pair, 8, 6,  "two tiles + pin: the half-grooves form one bore", (.58, .63, .71)),
    (frames_only, 25, 20, "frame only, all 12 tiles", (.55, .60, .68)),
    (asm, 25, 20, "finished: 12 tiles + 12 mirrors", (.62, .68, .76)),
    (asm, 90, 4, "straight-on: mirror flush with the 8 mm rib", (.62, .68, .76)),
]):
    draw(a, m, az, el, c)
    a.set_title(t, fontsize=10.5, color="#333")
fig.tight_layout()
fig.savefig(os.path.join(out, "preview.png"), dpi=118, facecolor="white")
print("rendered preview.png")
