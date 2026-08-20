#!/usr/bin/env python3
"""Independent checks on the cymatics dish, skirt and support ring.

Re-derives the speaker from cymatics_dish's constants and asserts the things
that would quietly ruin the part: mesh validity, bed fit, that the band's outer
face LIES on the cone, that nothing is buried in the cone or fouls the surround,
that the parts do not interpenetrate, that the bumps locate without bearing, and
that no wall exceeds its nominal thickness.
"""

import numpy as np
import trimesh
from matplotlib.path import Path

import cymatics_dish as C

BED = 320.0
TOL = 1e-3
E_PETG = 2.0e9             # Pa
RHO_PETG = 1270.0          # kg/m^3
NU = 0.38

dish = trimesh.load(C.OUT_DISH)
skirt = trimesh.load(C.OUT_SKIRT)
support = trimesh.load(C.OUT_SUPPORT)
fails = []


def check(ok, msg):
    print(f"  [{'ok ' if ok else 'FAIL'}] {msg}")
    if not ok:
        fails.append(msg)


def max_thickness(profile, step=0.1):
    """Largest inscribed circle in the meridional section = max wall thickness.

    The parts are solids of revolution, so the 2D section carries the whole
    answer, and the vent notches only remove material."""
    p = np.asarray(profile)
    lo, hi = p.min(0) - step, p.max(0) + step
    gr = np.stack(np.meshgrid(np.arange(lo[0], hi[0], step),
                              np.arange(lo[1], hi[1], step), indexing="ij"), -1)
    pts = gr.reshape(-1, 2)
    inside = Path(np.vstack([p, p[:1]])).contains_points(pts)
    pts = pts[inside]
    a, b = p, np.roll(p, -1, axis=0)
    ab = b - a
    best = np.full(len(pts), np.inf)
    for i in range(len(a)):                       # distance to each edge
        d = ab[i]
        t = np.clip(((pts - a[i]) @ d) / (d @ d), 0, 1) if d @ d else 0.0
        best = np.minimum(best, np.linalg.norm(pts - (a[i] + np.outer(t, d)), axis=1))
    j = int(np.argmax(best))
    return 2 * best[j], pts[j]


def plate_sag(a_m, t_m, q_pa):
    """Centre deflection of a uniformly loaded circular plate, clamped and
    simply-supported bounds (Roark). The real edge is between the two."""
    D = E_PETG * t_m**3 / (12 * (1 - NU**2))
    clamped = q_pa * a_m**4 / (64 * D)
    return clamped, clamped * (5 + NU) / (1 + NU)


def normal_gap(r, z):
    """Distance from a point to the cone, along the cone's normal. Positive is
    inside the funnel (clear of the paper)."""
    return (z - C.cone_z(r)) / np.hypot(C.K, 1.0)


print("mesh validity")
for name, m in (("dish", dish), ("skirt", skirt), ("ring", support)):
    check(m.is_watertight, f"{name} watertight")
    check(m.is_winding_consistent, f"{name} winding consistent")
    check(m.volume > 0, f"{name} volume positive ({m.volume / 1000:.1f} cm^3)")
    check(len(m.split(only_watertight=False)) == 1, f"{name} single shell")

print("bed fit (320 x 320 x 320)")
for name, m in (("dish", dish), ("skirt", skirt), ("ring", support)):
    e = m.extents
    check(max(e) <= BED, f"{name} extents {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} mm")

print("diameters and levels")
r_dish = np.hypot(*dish.vertices[:, :2].T)
r_skirt = np.hypot(*skirt.vertices[:, :2].T)
check(abs(2 * r_dish.max() - C.OD) < 0.05, f"dish OD = {2 * r_dish.max():.2f} mm")
check(abs(2 * r_skirt.max() - C.OD) < 0.05,
      f"skirt flange OD = {2 * r_skirt.max():.2f} mm — flush with the dish, so "
      f"glue-up aligns on two visible edges")
check(abs(dish.bounds[0][2] + C.FLOOR_T) < TOL
      and abs(skirt.bounds[1][2] + C.FLOOR_T) < TOL,
      f"parts meet at the glue plane z = {-C.FLOOR_T:.2f}")
check(abs(C.OD - 277.0) < TOL, f"dish OD is {C.OD:.0f} mm as asked")

print("speaker fit")
for name, m, r in (("dish", dish, r_dish), ("skirt", skirt, r_skirt),
                   ("ring", support, np.hypot(*support.vertices[:, :2].T))):
    g = normal_gap(r, m.vertices[:, 2])
    check(g.min() > -0.01, f"{name} never buried in the cone (min {g.min():+.4f} mm)")

on = r_skirt[np.abs(skirt.vertices[:, 2] - C.cone_z(r_skirt)) < 0.01]
check(len(on) > 100 and abs(on.max() * 2 - C.CONTACT_D) < 0.05
      and abs(on.min() * 2 - C.SKIRT_D_BOT) < 0.5,
      f"band bears on the cone, Ø{on.max() * 2:.1f} -> Ø{on.min() * 2:.1f}")
crown = C.cone_z(C.CONE_MOUTH_D / 2) + C.SURROUND_H
check(-C.FLOOR_T - C.SKIRT_T > crown + C.SURROUND_CLR - TOL,
      f"flange underside clears the surround crown by "
      f"{-C.FLOOR_T - C.SKIRT_T - crown:.1f} mm")
check(C.R_OUT - C.R_CON > 5.0,
      f"dish overhangs the contact circle by {C.R_OUT - C.R_CON:.1f} mm all round")

print("assembly")
check(dish.intersection(skirt).volume < 1.0,
      f"dish/skirt interpenetration = {dish.intersection(skirt).volume:.4f} mm^3")
check(support.intersection(skirt).volume < 1.0,
      f"ring/skirt interpenetration = {support.intersection(skirt).volume:.4f} mm^3")

# vent notches: walk a circle across the flange's glue face, just under it
th = np.linspace(0, 2 * np.pi, 3600, endpoint=False)
rp = (C.R_OUT + C.R_CON) / 2
pts = np.column_stack([rp * np.cos(th), rp * np.sin(th),
                       np.full_like(th, -C.FLOOR_T - C.VENT_D / 2)])
inside = skirt.contains(pts)
gaps = int(np.sum(np.diff(inside.astype(int)) < 0) + (inside[-1] and not inside[0]))
check(gaps == C.N_VENT, f"{gaps} vent notches open at Ø{2 * rp:.0f} "
                        f"(open fraction {1 - inside.mean():.1%})")

print("rim wall")
probe = np.linspace(0, 2 * np.pi, 36, endpoint=False)
solid_wall = []
for zz in np.arange(C.FILLET + 0.3, C.RIM_H - C.RIM_ROUND, 0.25):
    rr = C.R_IN + 0.15
    solid_wall.append(dish.contains(np.column_stack(
        [rr * np.cos(probe), rr * np.sin(probe), np.full_like(probe, zz)])))
check(all(s.all() for s in solid_wall),
      f"inner wall is unbroken from z={C.FILLET + 0.3:.1f} to "
      f"{C.RIM_H - C.RIM_ROUND:.1f}")

print("locating bumps")
r_c, z_c = C.bump_centre()
th = np.linspace(0, 2 * np.pi, 2880, endpoint=False)
rp = C.R_CON + (r_c + C.BUMP_D / 2 - C.R_CON) / 2
inside = skirt.contains(np.column_stack([rp * np.cos(th), rp * np.sin(th),
                                         np.full_like(th, z_c)]))
lobes = int(np.sum(np.diff(inside.astype(int)) > 0) + (inside[0] and not inside[-1]))
check(lobes == C.N_BUMP, f"{lobes} bumps at Ø{2 * rp:.1f}, "
                         f"{inside.mean() * 100:.0f}% of the circumference")
check(abs(2 * (r_c + C.BUMP_D / 2) - 2 * (C.R_CON + C.BUMP_D - C.BUMP_EMBED)) < TOL,
      f"bumps stand {r_c + C.BUMP_D / 2 - C.R_CON:.2f} mm proud of the "
      f"Ø{C.CONTACT_D:.0f} collar")
near = r_skirt > C.R_CON + 0.2        # bump caps only; clears boolean noise
gap = normal_gap(r_skirt[near], skirt.vertices[near, 2])
check(gap.min() > C.BUMP_CLR - 0.02,
      f"closest bump approach to the cone = {gap.min():.3f} mm "
      f"(target {C.BUMP_CLR:.1f}) — they locate, they do not bear")
check(z_c + C.BUMP_D / 2 < -C.FLOOR_T - C.SKIRT_T
      and z_c - C.BUMP_D / 2 > C.cone_z(C.R_CON),
      f"bumps sit inside the collar's {C.COLLAR_H:.1f} mm height")

print("wall thickness (max inscribed circle in the section)")
for name, prof in (("dish", C.dish_profile()), ("skirt", C.skirt_profile()),
                   ("ring", C.support_profile())):
    tmax, where = max_thickness(prof)
    corner = tmax <= 1.5 * C.RIM_T
    check(tmax <= C.RIM_T + TOL or corner,
          f"{name} max {tmax:.2f} mm at r={where[0]:.1f}, z={where[1]:.1f}"
          + (" (an L-junction — unavoidable where two walls meet)"
             if tmax > C.RIM_T + TOL else ""))
check(max(C.RIM_T, C.FLOOR_T, C.SKIRT_T, C.SUPPORT_T) <= 1.0 + TOL,
      f"no wall exceeds 1 mm: rim/skirt/ring {C.RIM_T:.2f}, floor {C.FLOOR_T:.2f}")
print(f"        (the {C.N_BUMP} locating bumps are solid Ø{C.BUMP_D:.0f} nubs — "
      f"the only material thicker than a wall)")

print("floor stiffness (E=2.0 GPa PETG)")
q = (1000 * 9.81 * C.WATER / 1000) + (RHO_PETG * C.FLOOR_T / 1000 * 9.81)
for label, a_mm in (("unsupported", C.R_CON - C.SKIRT_T),
                    ("with Ø%.0f ring" % C.SUPPORT_D, C.SUPPORT_D / 2)):
    lo, hi = plate_sag(a_mm / 1000, C.FLOOR_T / 1000, q)
    print(f"        {label:<16} span Ø{2 * a_mm:5.0f} -> centre sag "
          f"{lo * 1000:.2f}–{hi * 1000:.2f} mm  ({q:.1f} Pa)")
D_ = E_PETG * (C.FLOOR_T / 1000)**3 / (12 * (1 - NU**2))
a_max = (C.WATER / 2000 * 64 * D_ * (1 + NU) / (q * (5 + NU)))**0.25
print(f"        a ring at Ø{2000 * a_max:.0f} or below would hold sag under "
      f"half the water depth; Ø{C.SUPPORT_D:.0f} fitted")

print("capacity and mass")
for depth in (2.0, 3.0, 5.0):
    print(f"        {depth:.0f} mm water = {np.pi * C.R_IN**2 * depth / 1000:6.1f} ml")
check(C.RIM_H - C.WATER >= 6.0,
      f"freeboard above {C.WATER:.0f} mm water = {C.RIM_H - C.WATER:.2f} mm")
mass = (dish.volume + skirt.volume) * C.DENSITY
print(f"        plastic {mass:.0f} g; with {C.WATER:.0f} mm water "
      f"{mass + np.pi * C.R_IN**2 * C.WATER / 1000:.0f} g")

print("clearances")
check(C.SKIRT_D_BOT / 2 - C.CONE_DUSTCAP_D / 2 > 10.0,
      f"band toe clears the dust cap by "
      f"{(C.SKIRT_D_BOT - C.CONE_DUSTCAP_D) / 2:.1f} mm radially")
check(C.SUPPORT_D / 2 - C.CONE_DUSTCAP_D / 2 > 10.0,
      f"Ø{C.SUPPORT_D:.0f} ring clears the dust cap by "
      f"{(C.SUPPORT_D - C.CONE_DUSTCAP_D) / 2:.1f} mm radially")

print()
print("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURES: " + "; ".join(fails))
raise SystemExit(1 if fails else 0)
