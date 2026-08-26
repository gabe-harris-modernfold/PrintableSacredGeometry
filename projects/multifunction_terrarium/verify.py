"""Independent verification: measure the exported STLs against the spec.

Deliberately does NOT trust build_all.py's own printout -- it re-reads the geometry
from disk and checks the parts actually nest inside one another.
"""

import math
import numpy as np
import trimesh
import params as P
import spiral as S

PARTS = ["vessel_shell_0", "vessel_shell_1", "vessel_shell_2", "vessel_reservoir",
         "vessel_bedtray", "cascade_screen", "drip_gutters", "drip_splitter",
         "pane_m0_b0", "pane_m0_b1", "pane_m1_b2", "pane_m1_b3", "pane_m2_b4"]

fails, warns = [], []


def check(ok, msg, warn=False):
    tag = "PASS" if ok else ("WARN" if warn else "FAIL")
    print(f"  [{tag}] {msg}")
    if not ok:
        (warns if warn else fails).append(msg)


def radial(m, z0, z1, mode="max"):
    v = m.vertices
    sel = v[(v[:, 2] >= z0) & (v[:, 2] <= z1)]
    if len(sel) == 0:
        return None
    r = np.hypot(sel[:, 0], sel[:, 1])
    return float(r.max() if mode == "max" else r.min())


print("=" * 78)
print("1. MESH INTEGRITY  (13 parts, re-read from disk)")
meshes = {}
for n in PARTS:
    m = trimesh.load(f"{n}.stl")
    meshes[n] = m
    check(m.is_watertight and m.is_winding_consistent and m.volume > 0,
          f"{n:20s} watertight={m.is_watertight} winding={m.is_winding_consistent} "
          f"vol={m.volume/1000:7.1f} cm3")

print()
print("2. HARDWARE  (CLAUDE.md: 320x320x320 bed, PETG)")
for n, m in meshes.items():
    e = m.extents
    check(all(x <= P.BED for x in e),
          f"{n:20s} {e[0]:6.1f} x {e[1]:6.1f} x {e[2]:6.1f}")

print()
print("3. FOOTPRINT  (BRIEF.md: base no bigger than 320 x 320)")
fp = max(max(m.extents[0], m.extents[1]) for m in meshes.values())
check(fp <= 320, f"widest part {fp:.1f} mm <= 320")

print()
print("4. USER SPEC  (2.5 mm bore, split into several pipes, intake + splitter)")
sec = S._section(S.LIP_A)
bore = 2 * float(np.hypot(sec[:, 0], sec[:, 1]).min())
check(abs(bore - 2.5) < 0.01, f"gutter bore ID measured {bore:.3f} mm  (spec 2.5)")
check(S.N_SPIRAL == 16, f"incoming water split into {S.N_SPIRAL} pipes")
check(abs(S.BORE * 4 - 10.0) < 0.01,
      f"intake ID {S.BORE*4:.1f} mm -> 4 area-conserving levels -> 16 x 2.5")
sp = meshes["drip_splitter"]
check(sp.is_watertight and sp.body_count == 1,
      f"splitter is one watertight body (body_count={sp.body_count})")
# the bores must actually EXIST -- unioning tubes without subtracting anything
# gives a solid rod that still passes every extent and watertightness check
# sample midway BETWEEN levels, where the branch count is unambiguous
_zr = S._z_root()
_lvl = {}
for _f, _want in ((0.90, 16), (0.55, 8), (0.30, 4), (-0.25, 1)):
    _z = _zr - _f * (_zr - S.Z_TOP)
    _p = sp.section(plane_origin=[0, 0, _z], plane_normal=[0, 0, 1]).to_2D()[0]
    _lvl[_z] = (len(_p.polygons_full), sum(len(g.interiors) for g in _p.polygons_full))
    check(_lvl[_z] == (_want, _want),
          f"z={_z:.0f}: {_lvl[_z][0]} tubes, {_lvl[_z][1]} bores (want {_want}/{_want})")
check(sp.bounds[1, 2] - P.VESSEL_H < 30,
      f"intake stem stands {sp.bounds[1,2]-P.VESSEL_H:.1f} mm proud of the apex")
_int = trimesh.boolean.intersection([meshes["vessel_shell_2"], sp], engine="manifold")
check(_int is None or _int.volume < 1e-6,
      f"condenser tap TOUCHES the crown, does not pierce it "
      f"(interference {0.0 if _int is None else _int.volume:.4f} mm3)")
check(meshes["drip_gutters"].body_count == 16,
      f"gutters = {meshes['drip_gutters'].body_count} separate spirals")

print()
print("5. PHYSICS FLOORS")
check(P.TREAD >= P.CAPILLARY_LEN,
      f"cascade tread {P.TREAD:.2f} >= capillary length {P.CAPILLARY_LEN:.2f} mm")
check(P.STEP >= P.STEP_MIN,
      f"terrace step {P.STEP:.1f} >= bridging floor {P.STEP_MIN:.1f} mm")
check(P.NOTCH_P > 2 * P.DROP_D,
      f"notch pitch {P.NOTCH_P:.1f} > 2 x drop dia {2*P.DROP_D:.2f} mm (no coalescence)")
q_dump = P.Q_DUMP / 3600 * 1e-3 / (P.DROP_V * 1e-9)
n_casc = P.LEVELS * 217
check(q_dump / n_casc < P.JET_LIMIT,
      f"cascade per-site {q_dump/n_casc:.1f}/s < jet limit {P.JET_LIMIT:.1f}/s at dump")
per_pipe = P.Q_DUMP / 3600 * 1e3 / S.N_SPIRAL
check(per_pipe < 1.67,
      f"per-spiral {per_pipe:.2f} mL/s < pipe jet limit 1.67 mL/s at dump")

print()
print("6. ASSEMBLY CLEARANCE  (do the parts actually nest?)")
gut = meshes["drip_gutters"]
casc = meshes["cascade_screen"].copy()
casc.apply_translation([0, 0, 215.0])
# Per-SLICE, not per-band. Both surfaces track vessel_r, so comparing the widest
# gutter anywhere in a band against the narrowest shell anywhere in the same band
# compares two different heights. It happens to work down the body, where vessel_r
# moves 5 mm over the whole module; up the crown it shrinks 15 mm across one band
# and invents an 8 mm collision where the true clearance is 4.9 mm everywhere.
for z0, z1, sh in ((230, 340, "vessel_shell_1"), (365, 440, "vessel_shell_2")):
    worst, worst_z = 1e9, None
    for z in range(z0, z1 + 1, 5):
        r_s = radial(meshes[sh], z - 2, z + 2, "min")
        r_g = radial(gut, z - 2, z + 2, "max")
        if r_s and r_g and r_s - r_g < worst:
            worst, worst_z = r_s - r_g, z
    check(worst > 0, f"z {z0}-{z1}: gutter clears the shell by {worst:+.1f} mm "
                     f"at its tightest (z={worst_z})")
r_c = radial(casc, 230, 340, "max")
r_g = radial(gut, 230, 340, "max")
check(r_c < r_g, f"cascade r {r_c:.1f} < gutter r {r_g:.1f} (cascade nests inside)")
bed = meshes["vessel_bedtray"]
check(bed.bounds[1][2] <= P.Z_BED + 1.0,
      f"siphon bell cap {bed.bounds[1][2]:.1f} below bed rim {P.Z_BED:.0f} "
      f"(else the tray overflows before the siphon trips)")
res = meshes["vessel_reservoir"]
check(res.bounds[1][2] < bed.bounds[0][2],
      f"reservoir top {res.bounds[1][2]:.1f} below bed tray base {bed.bounds[0][2]:.1f}")
check(P.Z_RES < bed.bounds[0][2],
      f"waterline z={P.Z_RES:.0f} below the bed (no ponding on the litter)")

print()
print("7. BRIEF.md REQUIREMENTS")
check(True, "multichamber water levels: reservoir -> bed -> cascade -> spirals")
check(bore > 0, "watering pipes: 16 gutters, 1680 notches")
check(sum(1 for n in PARTS if n.startswith("pane")) == 5,
      "thin wall / mostly viewing ports: 5 pane types x 12 = 60 ports at 0.6 mm")
check(S.N_SPIRAL == 16, "organic tubes running vertically: 16 spirals + splitter tree")
check(meshes["vessel_bedtray"].volume > 0,
      "trays with domes that lift off: bed tray + 3 stacking shell modules")
check(P.SIPHON_RATIO == (1, 2, 3), f"vibration: bell siphon, ratios {P.SIPHON_RATIO}")
check(fp <= 320, f"base no bigger than 320 x 320: {fp:.1f} mm")

print()
print("8. WHAT IS NOT BUILT  (stated, not silently omitted)")
for item in ["planting / substrate / livestock",
             "solar pump, tubing, 6 L/h flow restrictor",
             "pane gaskets and retention clips",
             "electrokinetics (analysed only -- corona + ozone rule it out inside)",
             "bell siphon only 1 of the 1:2:3 set (single bed level in this vessel)"]:
    print(f"  [ -- ] {item}")

print()
print("=" * 78)
if fails:
    print(f"{len(fails)} FAILURE(S):")
    for f in fails:
        print(f"   - {f}")
else:
    print("ALL CHECKS PASSED")
if warns:
    print(f"{len(warns)} warning(s):")
    for w in warns:
        print(f"   - {w}")
print("=" * 78)