"""
Add concave "sound collector" cups around the rim of the scrying pool.

A raised torus COLLAR is fused onto the rim, and N concave spherical scoops are
carved into its inner face, each aimed down-and-inward at the water-surface
center (0,0,H). Sound striking a scoop is redirected/concentrated toward the
pool surface. Everything is above the water line, so the vessel stays
watertight, and the collar keeps every cup connected as one printable body.

Focal note: each scoop is a spherical mirror of focal length Rs/2 (near-field
with little cups). Raise CARVE_R for a flatter, longer-focus scoop; lower it
for a deeper, shorter-focus cup. AIM_Z sets how far below the collar the cups
point (the surface sits at z = H).
"""

import numpy as np
import trimesh

import scrying_pool as sp

# ---- parameters ----------------------------------------------------------
N_CUPS    = 21         # Fibonacci: cups all along the edge
COLLAR_A  = 30.0       # collar tube (minor) radius (mm)
COLLAR_DR = 8.0        # collar centerline beyond the water radius (mm)
COLLAR_Z  = 4.0        # collar lift above the rim (tube dips below the water
                       # line only as SOLID material -- it is never pierced there)
CARVE_R   = 34.0       # scoop sphere radius -> focal length CARVE_R/2
CARVE_D   = 16.0       # scoop bite depth into the collar (mm)
HORNS      = True      # add radial gather-to-center horns through the collar
HORN_Z     = 20.0      # horn-axis height above the water line (mouths stay dry)
HORN_MOUTH = 16.0      # outer mouth radius (gathers ambient sound)
HORN_THROAT = 6.0      # inner throat radius (emits toward center)
HORN_EXT   = 10.0      # mouth flare beyond the collar outer face (mm)
HORN_SEC   = 64        # horn revolve sections
POOL_SEC  = 360        # pool angular resolution for this build
TOR_MAJ   = 240        # torus major sections
TOR_MIN   = 48         # torus minor sections
SUB       = 3          # scoop icosphere subdivisions
OUT       = "scrying_pool_collectors.stl"
# --------------------------------------------------------------------------


def build():
    sp.SECTIONS = POOL_SEC
    pool, info = sp.build()
    rt = 0.5 * info["surf_dia"]
    H = info["height"]

    Rm = rt + COLLAR_DR                         # collar centerline radius
    collar = trimesh.creation.torus(Rm, COLLAR_A,
                                    major_sections=TOR_MAJ,
                                    minor_sections=TOR_MIN)
    collar.apply_translation([0, 0, H + COLLAR_Z])

    # every carving tool is clipped to ABOVE the water line (z >= H) so nothing
    # is ever removed from the retaining wall below it -> the pool cannot leak,
    # and each cup bottoms out exactly AT the water surface (surface touches cup)
    def clip_above_water(mesh):
        m = mesh.slice_plane([0, 0, H], [0, 0, 1], cap=True)
        return m if (m is not None and len(m.faces) > 0) else None

    off = COLLAR_A + CARVE_R - CARVE_D          # scoop-center offset inward
    carves = []
    for k in range(N_CUPS):
        th = 2 * np.pi * k / N_CUPS
        pt = np.array([Rm * np.cos(th), Rm * np.sin(th), H + COLLAR_Z])
        # aim down-and-inward at the water-surface center
        d = np.array([0.0, 0.0, H]) - pt
        d /= np.linalg.norm(d)
        c = trimesh.creation.icosphere(subdivisions=SUB, radius=CARVE_R)
        c.apply_translation(pt + d * off)
        c = clip_above_water(c)
        if c is not None:
            carves.append(c)

    # ---- gather-to-center horns (staggered, lifted clear of the water) ----
    if HORNS:
        r_in = Rm - COLLAR_A - 2.0                     # throat just inside face
        L = 2.0 * COLLAR_A + 4.0 + HORN_EXT            # mouth past outer face
        # solid frustum revolved about +z: throat (z=0) -> mouth (z=L)
        prof = np.array([[0, 0], [HORN_THROAT, 0], [HORN_MOUTH, L], [0, L]])
        for k in range(N_CUPS):
            tho = 2 * np.pi * (k + 0.5) / N_CUPS
            u = np.array([np.cos(tho), np.sin(tho), 0.0])   # radial (thru center)
            horn = trimesh.creation.revolve(prof, sections=HORN_SEC)
            horn.merge_vertices(); horn.fix_normals()
            horn.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], u))
            horn.apply_translation([r_in * u[0], r_in * u[1], H + HORN_Z])
            horn = clip_above_water(horn)
            if horn is not None:
                carves.append(horn)

    solid = trimesh.boolean.union([pool, collar])
    result = trimesh.boolean.difference([solid] + carves)

    # boolean can leave zero-volume sliver shells at the horn/cup seams;
    # weld first, then keep the real vessel body and drop the fragments
    result.merge_vertices()
    result = max(result.split(only_watertight=False), key=lambda c: c.area)
    result.remove_unreferenced_vertices()
    result.fix_normals()

    comps = result.split(only_watertight=False)

    # leak check: the basin below the water line must be a closed watertight cup
    below = result.slice_plane([0, 0, H - 0.5], [0, 0, -1], cap=True)
    holds = below is not None and below.is_watertight

    info.update(n_cups=N_CUPS, n_horns=N_CUPS if HORNS else 0,
                cup_focal=CARVE_R / 2.0,
                collar_dia=2 * (Rm + COLLAR_A), water_line=H, holds_water=holds,
                watertight=result.is_watertight, components=len(comps),
                verts=len(result.vertices), faces=len(result.faces))
    return result, info


def main():
    m, i = build()
    m.export(OUT)
    print(f"inner cups             : {i['n_cups']}  (focal ~{i['cup_focal']:.1f} mm)")
    print(f"gather-to-center horns : {i['n_horns']}  (mouths above water line)")
    print(f"water line z           : {i['water_line']:.2f} mm")
    print(f"HOLDS WATER (no leak)  : {i['holds_water']}  "
          f"<- basin below water line is watertight")
    print(f"collar outer diameter  : {i['collar_dia']/10:.2f} cm")
    print(f"water capacity         : {i['capacity_L']:.2f} L")
    print(f"watertight={i['watertight']}  components={i['components']}  "
          f"verts={i['verts']}  faces={i['faces']}")
    print(f"exported {OUT}")


if __name__ == "__main__":
    main()
