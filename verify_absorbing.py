#!/usr/bin/env python3
"""Checks on the absorbing-edge dish.

The one that matters most: the rim must be a CONTINUOUS solid wall all the way
round, at every height. The slots open into the chamber and stop there. An
earlier version ran the slot cut out to the dish's outer radius and perforated
the rim, which would have drained the dish onto the speaker -- so that check is
first, and it probes the wall itself rather than inferring it from the outline.
"""

import numpy as np
import trimesh

import cymatics_dish as C
import cymatics_dish_absorbing as A

mesh = trimesh.load(A.OUT)
base = trimesh.load(C.OUT_DISH)
fails = []


def check(ok, msg):
    print(f"  [{'ok ' if ok else 'FAIL'}] {msg}")
    if not ok:
        fails.append(msg)


def ring(r, z, n):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return mesh.contains(np.column_stack(
        [r * np.cos(th), r * np.sin(th), np.full_like(th, z)]))


N = 12 * A.N_SLOT                       # ~12 samples per slot pitch

print("watertight boundary — nothing may perforate the rim or the floor")
r_wall = C.R_OUT - C.RIM_T / 2
for z in np.arange(0.2, C.RIM_H - C.RIM_ROUND, 0.4):
    if not ring(r_wall, z, N).all():
        check(False, f"rim wall is PERFORATED at z={z:.1f} — the dish would leak")
        break
else:
    check(True, f"rim wall solid all round at r={r_wall:.1f}, "
                f"z=0.2 to {C.RIM_H - C.RIM_ROUND:.1f} ({N} probes per height)")
r_floor = (A.R_FI + A.R_FO) / 2
check(ring(r_floor, -C.FLOOR_T / 2, N).all(),
      "floor solid under every slot")
check(ring(C.R_OUT - C.RIM_T / 2, -C.FLOOR_T / 2, N).all(),
      "floor solid under the rim")
check(mesh.is_watertight and mesh.is_winding_consistent,
      "mesh watertight and consistently wound")
check(len(mesh.split(only_watertight=False)) == 1, "single shell")

print("fence")
mid = ring(r_floor, A.SLOT_H / 2, N)
opens = int(np.sum(np.diff((~mid).astype(int)) > 0) + ((~mid)[0] and not (~mid)[-1]))
check(opens == A.N_SLOT,
      f"{opens} slots open through the fence, {(1 - mid.mean()) * 100:.0f}% open "
      f"(built {A.POROSITY:.0%}; probe resolution ±{100 / 12:.0f}%)")
check(ring(r_floor, A.SLOT_H + (A.FENCE_H - A.SLOT_H) / 2, N).all(),
      f"lintel continuous above the slots ({A.FENCE_H - A.SLOT_H:.1f} mm band)")
check(A.FENCE_T <= 1.5, f"fence {A.FENCE_T:.1f} mm, within the 1.5 mm limit")

print("chamber")
r_ch = C.R_IN - A.CHAMBER / 2
for z in (0.5, 3.0, 5.0):
    check(not ring(r_ch, z, 720).any(),
          f"chamber clear all round at z={z:.1f} (r={r_ch:.2f})")
check(A.SLOT_R_OUT < C.R_IN,
      f"slot cut stops {C.R_IN - A.SLOT_R_OUT:.2f} mm short of the rim's inner face")

print("unchanged from the base dish")
check(abs(2 * np.hypot(*mesh.vertices[:, :2].T).max() - C.OD) < 0.05,
      f"OD {2 * np.hypot(*mesh.vertices[:, :2].T).max():.2f} mm")
check(abs(mesh.bounds[0][2] + C.FLOOR_T) < 1e-3
      and abs(mesh.bounds[1][2] - C.RIM_H) < 0.01,   # rim apex, arc faceting
      f"height {mesh.bounds[0][2]:+.2f} to {mesh.bounds[1][2]:+.2f} mm")
dm = (mesh.volume - base.volume) * C.DENSITY
# a leaky build weighs LESS than this, because its slots also ate rim material
check(2.0 < dm < 6.0, f"fence adds {dm:.1f} g of plastic and removes no rim")
check(max(mesh.extents) <= 320.0, f"fits the bed: {mesh.extents[0]:.1f} mm")

print("capacity")
print(f"        water surface Ø{2 * A.R_FI:.1f} (was Ø{2 * C.R_IN:.0f})")
for d in (2.0, 3.0, 5.0):
    print(f"        {d:.0f} mm water = {np.pi * A.R_FI**2 * d / 1000:6.1f} ml")

print()
print("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURES: " + "; ".join(fails))
raise SystemExit(1 if fails else 0)
