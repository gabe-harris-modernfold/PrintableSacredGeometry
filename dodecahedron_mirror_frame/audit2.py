"""
Follow-up on the two open questions from audit.py.

  G  what IS a single tile topologically? (to interpret the 24-body union)
  H  does a PAIRWISE union merge cleanly? (isolates iterated-boolean damage)
  I  tighten the bound on any gap between mating mitre faces
  J  can a pinned pair be joined by pure translation along the mitre normal?
     (this is the motion the pins DO permit; audit E only ruled out the
      radial wedge motion)
"""
import sys, os, math
import numpy as np
import trimesh
from trimesh.creation import cylinder

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dodeca_frame import (build_oriented, build_tile, face_rotations, prism,
                          edge_frame, r_i, GROOVE_Z, PIN_D, PIN_LEN, ENGINE)

V = lambda m: abs(m.volume) if m is not None and len(m.faces) else 0.0

solid, normals, fverts = build_oriented()
rot = math.degrees(math.atan2(fverts[0][1], fverts[0][0]))
tile = build_tile(solid, fverts, rot)
body = tile.copy(); body.apply_translation([0, 0, r_i])
rots = face_rotations(solid, normals)
frames = []
for _, T in rots:
    f = body.copy(); f.apply_transform(T); frames.append(f)

print("=" * 66)
print("G  topology of ONE tile")
parts = tile.split(only_watertight=False)
print(f"   connected bodies : {len(parts)}")
print(f"   watertight       : {tile.is_watertight}   euler: {tile.euler_number}")
print(f"   genus            : {(2 - tile.euler_number)//2} "
      f"(1 = a ring, as expected for a bezel with an open centre)")
print(f"   -> 12 rings, each euler 0, should union to euler 0 if merged.")
print(f"      audit.py reported euler 24 / 24 bodies for the 12-way union,")
print(f"      which is not a valid closed-surface value -> boolean artifact.")

print("=" * 66)
print("H  pairwise union of two mating tiles")
pair = trimesh.boolean.union([frames[0], frames[7]], engine=ENGINE)
pp = pair.split(only_watertight=False)
print(f"   bodies  : {len(pp)}  (1 = the two tiles fused, so faces do touch)")
print(f"   volume  : {V(pair)/1000:.3f} cm3  vs 2x tile "
      f"{2*V(tile)/1000:.3f} cm3   diff {abs(V(pair)-2*V(tile)):.3f} mm3")
print(f"   watertight: {pair.is_watertight}  euler: {pair.euler_number}")

print("=" * 66)
print("I  gap between mating mitre faces (inward push, fine steps)")
for delta in (0.005, 0.01, 0.02, 0.05, 0.10, 0.25):
    t = frames[0].copy(); t.apply_translation([0, 0, -delta])
    v = max(V(trimesh.boolean.intersection([t, frames[i]], engine=ENGINE))
            for i in range(1, 12))
    print(f"   -{delta:5.3f} mm inward: interference {v:9.2f} mm3"
          f"{'   <-- faces already touching' if v > 0.5 else ''}")

print("=" * 66)
print("J  pinned pair, separated along the mitre normal")
a = np.array(fverts[0], float)
b = np.array(fverts[1], float)
mid, e = edge_frame(fverts, 0)
radial = np.array([mid[0], mid[1], 0.0]); radial /= np.linalg.norm(radial)
n = np.cross(e, -a); n /= np.linalg.norm(n)
if np.dot(n, radial) < 0:
    n = -n
p = mid * (1.0 - GROOVE_Z / r_i)
pin = cylinder(radius=PIN_D / 2.0 - 0.05, sections=48,
               segment=[p - e * PIN_LEN / 2, p + e * PIN_LEN / 2])
print(f"   mitre normal n : {np.round(n, 4)}  (oriented TOWARD the neighbour,")
print(f"                    so separating tile 0 means translating along -n)")
for sign, lbl in ((-1.0, "-n (separate)"), (+1.0, "+n (into neighbour)")):
    print(f"   --- {lbl} ---")
    for delta in (0.0, 0.5, 2.0, 5.0, 15.0):
        t = frames[0].copy(); t.apply_translation(n * sign * delta)
        vt = V(trimesh.boolean.intersection([t, frames[7]], engine=ENGINE))
        vp = V(trimesh.boolean.intersection([t, pin], engine=ENGINE))
        print(f"     {delta:5.2f} mm: vs neighbour {vt:8.2f}   vs pin {vp:8.2f} mm3"
              f"{'   clean' if vt < 1 and vp < 1 else '   BLOCKED'}")
