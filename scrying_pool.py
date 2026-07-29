"""
Golden-ratio scrying pool.

A perfect circle of revolution: a small flat base disc that flares up a smooth
half-cosine wall to a wide, flat water surface -- a shallow, rigid, radially
symmetric vessel meant to show standing-wave / cymatic patterns on the water.

GOLDEN LADDER
-------------
The 10 cm depth cap is used as the golden SEED (u). Every principal dimension
is a power of phi off that seed:

        depth  :  base diameter  :  surface diameter   =   1 : phi^2 : phi^4

    depth               = u            = 100.000 mm
    inner base    dia   = u * phi^2    = 261.803 mm
    inner surface dia   = u * phi^4    = 685.410 mm   (area ~ 3690 cm^2)

CONSTRUCTION
------------
Two solids of revolution -- an outer envelope and the inner water cavity --
are subtracted (outer - cavity) to form a watertight shell:

  * wall / floor thickness  t = 2*phi mm  (a multiple of phi, < 5 mm, ~3 mm),
    applied as a true surface-NORMAL offset so the shallow flare keeps uniform
    thickness instead of thinning in the middle,
  * a rounded fillet blends the flat base into the wall (smooth base->surface),
  * a flat golden rim lip (width RIM_LIP) stiffens the perimeter.
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
PHI      = (1.0 + np.sqrt(5.0)) / 2.0     # golden ratio
SEED     = 100.0            # u: golden seed sizing base/surface (mm)
DEPTH    = SEED / PHI**3    # water depth (mm): u/phi^3 -> shallow, responsive
SCALE    = 1.0             # multiply every golden dimension (1.0 = full size)
WALL     = 2.0 * PHI       # wall / floor thickness t (mm) = a multiple of phi
FILLET   = True            # round the base<->wall inner corner (rf = t*phi)
RIM_LIP  = 2.0 * PHI**3    # rim annulus width at the top (mm), golden
SECTIONS = 720             # angular facets (~3 mm chords on the rim circle)
NWALL    = 320             # samples along the flared wall profile
NFLOOR   = 90              # samples across the flat base disc
ARC      = 28              # (unused) kept for backward compatibility
OUT      = "scrying_pool.stl"
# --------------------------------------------------------------------------


def smoothstep_wall(r0, r1, z0, z1, n):
    """Half-cosine flare from (r0,z0) to (r1,z1): zero slope at both ends."""
    s = np.linspace(0.0, 1.0, n)
    r = r0 + (r1 - r0) * (1.0 - np.cos(np.pi * s)) / 2.0
    z = z0 + (z1 - z0) * s
    return np.column_stack([r, z])


def outward_normals(pts):
    """Unit normals of a polyline, rotated to point toward +r / -z (out of
    the water, into the wall material)."""
    d = np.gradient(pts, axis=0)                 # tangents (dr, dz)
    n = np.column_stack([d[:, 1], -d[:, 0]])     # rotate -90deg -> (dz, -dr)
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    return n


def smoothstep(x):
    """Smoothstep 3x^2 - 2x^3 on [0,1]; zero slope at both ends."""
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def inner_profile(rb, rt, zf, H):
    """Water-contact curve, centre -> rim: a flat base disc that blends
    TANGENTIALLY into a smooth z(r) flared wall — no fillet, no vertical
    riser, so the bowl bottom is one continuous smooth surface (no ring)."""
    fr = np.linspace(0.0, rb, NFLOOR)
    floor = np.column_stack([fr, np.full_like(fr, zf)])
    wr = np.linspace(rb, rt, NWALL)
    wz = zf + (H - zf) * smoothstep((wr - rb) / (rt - rb))
    wall = np.column_stack([wr, wz])
    return np.vstack([floor, wall[1:]])


def build():
    u   = SEED * SCALE
    t   = WALL
    d   = DEPTH * SCALE            # water depth (decoupled from the seed)
    rb  = 0.5 * u * PHI**2         # flat base radius
    rt  = 0.5 * u * PHI**4         # inner surface (water) radius
    zf  = t                        # inner floor level (floor slab thickness t)
    H   = t + d                    # rim / water-surface level

    inner = inner_profile(rb, rt, zf, H)               # water-contact surface
    outer = inner + t * outward_normals(inner)         # uniform normal-offset shell

    # closed cross-section: inner (centre->rim) + outer (rim->centre); the
    # (rt,H)->(rt,H-t) step is the thin rim edge, revolved to a solid vessel
    prof = np.vstack([inner, outer[::-1]])
    vessel = trimesh.creation.revolve(prof, sections=SECTIONS)
    vessel.merge_vertices()
    vessel = max(vessel.split(only_watertight=False), key=lambda c: c.area)
    vessel.remove_unreferenced_vertices()
    vessel.fix_normals()

    # water capacity: revolve the cavity (inner surface + flat top cap)
    cavity = trimesh.creation.revolve(np.vstack([inner, [[0.0, H]]]),
                                      sections=SECTIONS)

    info = dict(
        u=u, t=t, depth=d, rf=0.0, rim=t,
        base_dia=2 * rb, surf_dia=2 * rt,
        surf_area_cm2=np.pi * rt**2 / 100.0,
        outer_dia=2 * outer[:, 0].max(), height=H,
        capacity_L=cavity.volume / 1e6,
        material_L=vessel.volume / 1e6,
        watertight=vessel.is_watertight,
        winding=vessel.is_winding_consistent,
        verts=len(vessel.vertices), faces=len(vessel.faces),
    )
    return vessel, info


def main():
    mesh, i = build()
    mesh.export(OUT)
    print(f"golden seed u          : {i['u']:.3f} mm")
    print(f"wall thickness t=2*phi : {i['t']:.3f} mm")
    print(f"floor->wall blend      : tangential smooth (no ring)")
    print(f"water depth            : {i['depth']:.3f} mm  ({i['depth']/10:.2f} cm)")
    print(f"inner base   diameter  : {i['base_dia']:.3f} mm")
    print(f"inner surface diameter : {i['surf_dia']:.3f} mm  "
          f"(area {i['surf_area_cm2']:.0f} cm^2)")
    print(f"outer diameter         : {i['outer_dia']:.3f} mm  "
          f"({i['outer_dia']/10:.2f} cm)")
    print(f"overall height         : {i['height']:.3f} mm  ({i['height']/10:.2f} cm)")
    print(f"water capacity         : {i['capacity_L']:.2f} L")
    print(f"material volume        : {i['material_L']:.3f} L")
    print(f"watertight={i['watertight']}  winding_ok={i['winding']}  "
          f"verts={i['verts']}  faces={i['faces']}")
    print(f"exported {OUT}")


if __name__ == "__main__":
    main()
