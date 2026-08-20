"""
Dimensional review: are the chosen sizes appropriate, as opposed to merely
geometrically valid?

  K  first-layer footprint and cross-section through the print height
  L  actual wall widths at the places that matter
  M  what size mirror does the pocket really accept? (empirical, not nominal)
  N  filament / mass estimates
  O  bed fit against real printer sizes
"""
import sys, os, math
import numpy as np
import trimesh
from shapely.geometry import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dodeca_frame import (build_oriented, build_tile, prism, pent_pts,
                          r_i, face_inr, pocket_inr, open_inr, mirror_inr,
                          MIRROR_SIDE, MIRROR_THK, THK, CLEAR, RIB, LEDGE,
                          GROOVE_Z, GROOVE_R, P_INR, ENGINE)

V = lambda m: abs(m.volume) if m is not None and len(m.faces) else 0.0
solid, normals, fverts = build_oriented()
rot = math.degrees(math.atan2(fverts[0][1], fverts[0][0]))
tile = build_tile(solid, fverts, rot)                       # local: z 0 .. -14
flat = tile.copy(); flat.apply_translation([0, 0, THK])     # print: z 0 .. 14


def xsec_area(mesh, z):
    s = mesh.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
    if s is None:
        return 0.0, 0
    p, _ = s.to_2D()
    polys = p.polygons_full
    return sum(q.area for q in polys), len(polys)


print("=" * 70)
print("K  cross-section through the print height (z=0 is the bed)")
print("      z (mm)   area (mm2)   loops   note")
notes = {0.2: "first layer -> bed adhesion",
         5.0: "groove height",
         9.8: "just below pocket floor",
         10.2: "just above pocket floor (rim only)",
         13.8: "top / mirror-front plane"}
for z in (0.2, 2.0, 5.0, 9.8, 10.2, 13.8):
    a, n = xsec_area(flat, z)
    print(f"    {z:6.1f}   {a:9.1f}   {n:5d}   {notes.get(z,'')}")

print("=" * 70)
print("L  wall widths (radial, at the edge midpoints)")
mid_inr = lambda z: face_inr * (1.0 - abs(z) / r_i)      # mitre, local z
hol_inr = lambda z: open_inr * (1.0 - abs(z) / r_i)      # tapered hollow
print(f"    rim wall at face plane      : {face_inr - pocket_inr:5.2f} mm  "
      f"(visible rib per side)")
print(f"    ledge at pocket floor       : {mid_inr(-MIRROR_THK) - hol_inr(-MIRROR_THK):5.2f} mm")
print(f"    bottom ring (bed contact)   : {mid_inr(-THK) - hol_inr(-THK):5.2f} mm")
print(f"    groove crown -> pocket floor: {(GROOVE_Z - GROOVE_R) - MIRROR_THK:5.2f} mm")
print(f"    groove crown -> inner face  : {THK - (GROOVE_Z + GROOVE_R):5.2f} mm")
print(f"    groove wall radially inward : "
      f"{(mid_inr(-GROOVE_Z) - GROOVE_R*0.8507) - hol_inr(-GROOVE_Z):5.2f} mm")
print(f"    -> 0.4 mm nozzle: divide by ~0.42 for perimeter count")

print("=" * 70)
print("M  what mirror size does the pocket ACTUALLY accept?")
lo, hi = MIRROR_SIDE - 1.0, MIRROR_SIDE + 4.0
for _ in range(22):                      # bisect on largest fitting pentagon
    mid = (lo + hi) / 2.0
    test = prism(mid * P_INR, r_i - MIRROR_THK, r_i, rot)
    test.apply_translation([0, 0, -r_i])
    if V(trimesh.boolean.intersection([tile, test], engine=ENGINE)) < 0.01:
        lo = mid
    else:
        hi = mid
print(f"    nominal cut          : {MIRROR_SIDE:.2f} mm side")
print(f"    largest that seats   : {lo:.2f} mm side")
print(f"    usable oversize      : +{lo - MIRROR_SIDE:.2f} mm  "
      f"(radial clearance CLEAR = {CLEAR} mm)")
print(f"    hand-cut glass on a traced template is realistically +/-0.5 to 1 mm")
print(f"    per edge, and all five edges must fit -> this is the tight spot.")
for c in (0.6, 1.0, 1.5, 2.0):
    side = (mirror_inr + c) / P_INR
    print(f"      CLEAR {c:.1f} mm -> accepts up to {side:.2f} mm  "
          f"(+{side - MIRROR_SIDE:.2f}), rim wall {face_inr - (mirror_inr + c):.2f} mm")

print("=" * 70)
print("N  material")
solid_cm3 = 12 * V(tile) / 1000.0
print(f"    12 tiles, solid      : {solid_cm3:.0f} cm3")
xs, _ = xsec_area(flat, 5.0)
print(f"    ring cross-section   : {xs:.0f} mm2 at mid-height")
for frac, lbl in ((1.0, "100% solid"), (0.46, "3 perim + 15% infill")):
    g = solid_cm3 * frac * 1.24
    print(f"    {lbl:22s}: {solid_cm3*frac:5.0f} cm3 -> {g:5.0f} g PLA")
glass = 12 * V(prism(mirror_inr, 0, MIRROR_THK, rot)) / 1000.0
print(f"    glass at {MIRROR_THK} mm       : {glass:.0f} cm3 -> {glass*2.5/1000:.2f} kg")
print(f"    glass at 3 mm         : {glass*3/MIRROR_THK:.0f} cm3 -> "
      f"{glass*3/MIRROR_THK*2.5/1000:.2f} kg")

print("=" * 70)
print("O  bed fit")
e = tile.extents
print(f"    tile         : {e[0]:.1f} x {e[1]:.1f} mm  (pentagon diameter is")
print(f"                   1.902 x circumradius, so no rotation beats {e[0]:.1f})")
for name, bed in (("Prusa Mini", 180), ("Ender 3 / MK4", 220),
                  ("Bambu P1S/X1", 256), ("Prusa XL", 360)):
    bare = "fits" if e[0] <= bed else "NO"
    brim = "fits" if e[0] + 10 <= bed else "no brim"
    print(f"    {name:15s} {bed:3d} mm : bare {bare:4s} | +5mm brim {brim}")
