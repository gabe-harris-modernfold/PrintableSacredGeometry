#!/usr/bin/env python3
"""Assembly renders: each POC put together from its printed parts.

The STLs are laid out for the PLATE (every sphere at angle 0, every segment on
its own axis); this script applies the assembly transforms instead, which is
also a check that the parts do not interfere where they are supposed to mate.

    python preview_set.py            -> poc-set.png plus one image per POC
"""
import math
import os
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "lib"))
import render_mesh as R

STL = "stl/"


def load(name, rot_z=0.0, shift=(0, 0, 0)):
    m = trimesh.load(STL + name)
    if rot_z:
        m.apply_transform(trimesh.transformations.rotation_matrix(rot_z, (0, 0, 1)))
    m.apply_translation(shift)
    return m


def clarion_projector(shift=(0, 0, 0)):
    return [(load("clarion_p_throat.stl", shift=shift), "frame"),
            (load("clarion_p_bell.stl", shift=shift), "pad")]


def clarion_squillo(shift=(0, 0, 0)):
    return [(load("clarion_s_throat.stl", shift=shift), "frame"),
            (load("clarion_s_bell.stl", shift=shift), "pad"),
            (load("clarion_s_qrs_ring.stl", shift=shift), "frame")]


def halo(shift=(0, 0, 0)):
    """Spheres and plugs are stored in their PLATE positions; this is where they
    actually go -- which doubles as the chamber/sphere interference check."""
    import halo as HL
    out = [(load("halo_chamber.stl", shift=shift), "frame")]
    for i in range(4):
        a = 2 * math.pi * i / 4
        for h in "ab":
            m = load(f"halo_sphere{i+1}_{h}.stl", shift=(0, 0, HL.Z_NECK))
            m.apply_transform(trimesh.transformations.rotation_matrix(a, (0, 0, 1)))
            m.apply_translation(shift)
            out.append((m, "pad" if h == "a" else "frame"))
        r_in = HL.sphere_radius(HL.V_TGT[i])
        ctr = HL.R_CH + HL.L_N + r_in
        _, travel = HL.plug_dims(i, r_in)
        p = trimesh.load(STL + f"halo_plug{i+1}.stl")
        p.apply_transform(trimesh.transformations.rotation_matrix(-math.pi / 2,
                                                                  (1, 0, 0)))
        p.apply_translation((ctr, r_in - (travel + 24.0), HL.Z_NECK))
        p.apply_transform(trimesh.transformations.rotation_matrix(a, (0, 0, 1)))
        p.apply_translation(shift)
        out.append((p, "glass"))
    return out


def volute(shift=(0, 0, 0)):
    return [(load("volute_body.stl", shift=shift), "frame"),
            (load("volute_bell.stl", shift=shift), "pad"),
            (load("volute_mouthpiece.stl", shift=shift), "glass")]


def shoot(groups, out, az=38, el=20, w=1500, h=950):
    layers = [(np.asarray(m.vertices), np.asarray(m.faces), mat)
              for m, mat in groups]
    R.render(layers, az, el, w=w, h=h, out=out)


if __name__ == "__main__":
    shoot(clarion_projector(), "poc-a-projector.png")
    shoot(clarion_squillo(), "poc-a-squillo.png")
    shoot(halo(), "poc-c-halo.png", az=30, el=16)
    shoot(volute(), "poc-b-volute.png", az=200, el=26)
    # spread along the SCREEN-horizontal axis for the camera azimuth, or the
    # world-space offsets just stack the instruments on top of each other
    az = 28.0
    u = np.array([-math.sin(math.radians(az)), math.cos(math.radians(az)), 0.0])
    at = lambda d: tuple(u * d)
    shoot(clarion_projector(at(-980)) + clarion_squillo(at(-520))
          + halo(at(0)) + volute(at(620)),
          "poc-set.png", az=az, el=14, w=2000, h=800)
    print("wrote poc-a-projector.png poc-a-squillo.png poc-b-volute.png "
          "poc-c-halo.png poc-set.png")
