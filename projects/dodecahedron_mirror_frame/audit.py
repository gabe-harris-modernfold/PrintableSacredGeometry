"""
Design audit. Checks things the build-time asserts do not cover:

  A  do the 12 tiles union into ONE closed shell (no gaps, not just no overlap)?
  B  do adjacent mirrors clear each other at the 116.565 deg fold?
  C  does any groove break through into the mirror pocket?
  D  can a tile be withdrawn/inserted radially (is the socket a true taper)?
  E  with pins fitted, can the closing tile still be inserted?
  F  bed fit including brim
"""
import sys, os, math
import numpy as np
import trimesh
from trimesh.creation import cylinder

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dodeca_frame import (build_oriented, build_tile, face_rotations, prism,
                          edge_frame, r_i, mirror_inr, face_inr, pocket_inr,
                          open_inr, MIRROR_THK, THK, GROOVE_Z, GROOVE_R,
                          GROOVE_LEN, PIN_D, PIN_LEN, RIB, CLEAR, ENGINE)

V = lambda m: abs(m.volume) if m is not None and len(m.faces) else 0.0

solid, normals, fverts = build_oriented()
rot = math.degrees(math.atan2(fverts[0][1], fverts[0][0]))
tile = build_tile(solid, fverts, rot)                 # local frame
body = tile.copy(); body.apply_translation([0, 0, r_i])   # body frame
mirror = prism(mirror_inr, r_i - MIRROR_THK, r_i, rot)
rots = face_rotations(solid, normals)

frames, mirrors = [], []
for _, T in rots:
    f = body.copy();   f.apply_transform(T); frames.append(f)
    m = mirror.copy(); m.apply_transform(T); mirrors.append(m)

print("=" * 66)
print("A  union of the 12 tiles")
uni = frames[0]
for f in frames[1:]:
    uni = trimesh.boolean.union([uni, f], engine=ENGINE)
bodies = uni.split(only_watertight=False)
print(f"   watertight      : {uni.is_watertight}")
print(f"   separate bodies : {len(bodies)}  (1 = ball closes with no gaps)")
print(f"   volume          : {V(uni)/1000:.2f} cm3 vs 12x tile "
      f"{12*V(tile)/1000:.2f} cm3   diff {abs(V(uni)-12*V(tile)):.2f} mm3")
print(f"   genus/euler     : euler_number {uni.euler_number}")

print("=" * 66)
print("B  adjacent mirror clearance")
clash = max(V(trimesh.boolean.intersection([mirrors[0], mirrors[i]], engine=ENGINE))
            for i in range(1, 12))
print(f"   mirror-mirror overlap : {clash:.4f} mm3 -> "
      f"{'PASS' if clash < 1e-3 else 'FAIL'}")
d = []
for i in range(1, 12):
    a, _, _ = trimesh.proximity.closest_point(mirrors[i], mirrors[0].vertices)
    b, _, _ = trimesh.proximity.closest_point(mirrors[0], mirrors[i].vertices)
    d.append(min(np.linalg.norm(mirrors[0].vertices - a, axis=1).min(),
                 np.linalg.norm(mirrors[i].vertices - b, axis=1).min()))
print(f"   closest neighbour gap : {min(d):.2f} mm (vertex-to-surface, "
      f"upper bound on true min)")

print("=" * 66)
print("C  groove vs mirror pocket")
mloc = mirror.copy(); mloc.apply_translation([0, 0, -r_i])
worst = 0.0
for k in range(5):
    mid, e = edge_frame(fverts, k)
    p = mid * (1.0 - GROOVE_Z / r_i) - np.array([0, 0, r_i])
    g = cylinder(radius=GROOVE_R, sections=48,
                 segment=[p - e * GROOVE_LEN / 2, p + e * GROOVE_LEN / 2])
    worst = max(worst, V(trimesh.boolean.intersection([g, mloc], engine=ENGINE)))
print(f"   groove into pocket : {worst:.4f} mm3 -> "
      f"{'PASS' if worst < 1e-3 else 'FAIL'}")
print(f"   groove crown to pocket floor : "
      f"{(GROOVE_Z - GROOVE_R) - MIRROR_THK:.2f} mm of material")
print(f"   groove crown to inner face   : "
      f"{THK - (GROOVE_Z + GROOVE_R):.2f} mm of material")

print("=" * 66)
print("D  radial withdrawal / insertion of one tile (no pins)")
for delta in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0):
    t = frames[0].copy(); t.apply_translation([0, 0, delta])   # outward
    out = max(V(trimesh.boolean.intersection([t, frames[i]], engine=ENGINE))
              for i in range(1, 12))
    t2 = frames[0].copy(); t2.apply_translation([0, 0, -delta])  # inward
    inn = max(V(trimesh.boolean.intersection([t2, frames[i]], engine=ENGINE))
              for i in range(1, 12))
    print(f"   +{delta:5.2f} mm outward: {out:9.1f} mm3   "
          f"-{delta:5.2f} mm inward: {inn:9.1f} mm3")

print("=" * 66)
print("E  same withdrawal, but with the 5 pins fitted")
pins = []
for k in range(5):
    mid, e = edge_frame(fverts, k)
    p = mid * (1.0 - GROOVE_Z / r_i)
    pins.append(cylinder(radius=PIN_D / 2.0 - 0.05, sections=48,
                         segment=[p - e * PIN_LEN / 2, p + e * PIN_LEN / 2]))
for delta in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0):
    t = frames[0].copy(); t.apply_translation([0, 0, delta])
    hit = sum(V(trimesh.boolean.intersection([t, p], engine=ENGINE)) for p in pins)
    print(f"   +{delta:5.2f} mm: pin interference {hit:8.2f} mm3"
          f"{'   <-- blocked' if hit > 1 else ''}")

print("=" * 66)
print("F  geometry / fit summary")
e = tile.extents
print(f"   tile bbox          : {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} mm")
print(f"   + 5 mm brim        : {e[0]+10:.1f} x {e[1]+10:.1f} mm")
print(f"   rim wall at edge   : {face_inr - pocket_inr:.2f} mm")
print(f"   mirror-to-mirror   : {2*RIB:.1f} mm along the faces "
      f"({2*(RIB+CLEAR):.1f} mm incl. clearance)")
print(f"   ledge at pocket flr: "
      f"{pocket_inr - open_inr*(1 - MIRROR_THK/r_i):.2f} mm")
print(f"   mitre face area    : ~{V(tile)/1000:.1f} cm3 tile, "
      f"5 joints per tile")
