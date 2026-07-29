"""
Add a domed pedestal STAND beneath the scrying pool.

The stand is a third energy collector: a convex hemispherical cup turned upside
down — dome pointing UP, hollow opening DOWN. Its rim rests on the ground and
its truncated apex forms a flat seat fused to the centre of the pool's base, so
the dome's axis points straight at the centre of the pool. Substrate vibration
entering the wide rim is carried up the shell and concentrated at the apex,
driving the pool's central antinode from below — the counterpart to the horns
(which gather from around) and the cups (which drive the surface).

The STAND HEIGHT is golden: H_stand = u / phi = 61.80 mm  (u = 100 mm seed).
"""

import numpy as np
import trimesh
import scrying_pool as sp
import add_collectors as ac

# ---- parameters ----------------------------------------------------------
STAND_H = sp.SEED / sp.PHI            # golden stand height (mm) = u/phi
R_SEAT  = sp.SEED / sp.PHI**2         # flat apex-seat radius (mm) = u/phi^2 (widened)
MOUTH_R = 2.0 * np.hypot(STAND_H, 0.5*sp.SEED/sp.PHI**2)  # mouth radius — frozen at 2x
WALL    = sp.WALL                     # dome wall = 2 phi (matches the pool)
NARC    = 220                         # samples along each dome arc
SECT    = 360                         # revolve sections
OUT     = "scrying_pool_stand.stl"
# --------------------------------------------------------------------------


def dome():
    """Spherical-cap shell: wide mouth (rim) at z=0, narrow apex seat at
    z=STAND_H, hollow opening downward (an inverted cup, convex up)."""
    Hs, r_s, r_m, wall = STAND_H, R_SEAT, MOUTH_R, WALL
    z_c = (r_s**2 + Hs**2 - r_m**2) / (2*Hs)   # sphere centre on the axis
    Rs = np.hypot(r_m, z_c)                    # outer sphere radius
    Ri = Rs - wall                             # inner sphere radius
    z_in = Hs - wall                           # underside of the seat
    zo = np.linspace(Hs, 0.0, NARC)
    outer = np.column_stack([np.sqrt(np.clip(Rs**2-(zo-z_c)**2, 0, None)), zo])
    zi = np.linspace(0.0, z_in, NARC)
    inner = np.column_stack([np.sqrt(np.clip(Ri**2-(zi-z_c)**2, 0, None)), zi])
    prof = np.vstack([[0.0, Hs], outer, inner, [0.0, z_in]])
    m = trimesh.creation.revolve(prof, sections=SECT)
    m.merge_vertices(); m.fix_normals()
    return m, r_m


def build():
    pool, info = ac.build()
    H = info["height"]                       # pool height (base->rim)
    pool.apply_translation([0, 0, STAND_H])  # lift pool onto the stand
    dm, Rd = dome()

    res = trimesh.boolean.union([pool, dm])
    res.merge_vertices()
    res = max(res.split(only_watertight=False), key=lambda c: c.area)
    res.remove_unreferenced_vertices(); res.fix_normals()

    water_z = STAND_H + H
    below = res.slice_plane([0, 0, water_z - 0.5], [0, 0, -1], cap=True)
    info.update(stand_h=STAND_H, dome_rad=Rd, dome_dia=2*Rd, r_seat=R_SEAT,
                water_z=water_z, total_h=water_z,
                holds_water=(below is not None and below.is_watertight),
                watertight=res.is_watertight,
                components=len(res.split(only_watertight=False)),
                faces=len(res.faces))
    return res, info


def main():
    m, i = build()
    m.export(OUT)
    print(f"stand height  H=u/phi   : {i['stand_h']:.3f} mm")
    print(f"dome diameter (footprint): {i['dome_dia']:.1f} mm")
    print(f"apex seat radius         : {i['r_seat']:.2f} mm")
    print(f"water surface z          : {i['water_z']:.2f} mm")
    print(f"overall height           : {i['total_h']:.2f} mm  ({i['total_h']/10:.2f} cm)")
    print(f"HOLDS WATER (no leak)    : {i['holds_water']}")
    print(f"watertight={i['watertight']}  components={i['components']}  faces={i['faces']}")
    print(f"exported {OUT}")


if __name__ == "__main__":
    main()
