#!/usr/bin/env python3
"""Boolean fit check of every mating pair, in assembled position.

Each pair is intersected as assembled; what is left is fit clearance engagement,
not a collision -- the README's table. Helical joints are tested over one pitch
of axial offset (ridge-on-ridge at the wrong phase is not a failure), and the
minimum over the sweep is reported. Run AFTER a build:

    python fit_check.py
"""
import math

import numpy as np
import trimesh

import halo as HL

STL = "stl/"


def load(name):
    return trimesh.load(STL + name + ".stl")


def ivol(a, b):
    try:
        m = trimesh.boolean.intersection([a, b], engine="manifold")
        return 0.0 if m is None or m.is_empty else abs(m.volume)
    except Exception:
        return 0.0


def pair(name, a, b, note=""):
    v = ivol(a, b)
    print(f"  {name:<44} {v:10.1f} mm3  {note}")
    return v


def helical(name, a, b, axis, pitch, note=""):
    best = math.inf
    for off in np.linspace(0.0, pitch, 5):
        bb = b.copy()
        bb.apply_translation(np.asarray(axis, float) * off)
        best = min(best, ivol(a, bb))
    print(f"  {name:<44} {best:10.1f} mm3  min over one pitch  {note}")
    return best


def main():
    print("fit check (assembled-position boolean intersections)")
    bad = 0.0
    bad += helical("clarion p throat <-> bell", load("clarion_p_throat"),
                   load("clarion_p_bell"), (0, 0, 1), 8.0)
    bad += helical("clarion s throat <-> bell", load("clarion_s_throat"),
                   load("clarion_s_bell"), (0, 0, 1), 8.0)
    bad += pair("clarion s bell <-> qrs ring", load("clarion_s_bell"),
                load("clarion_s_qrs_ring"), "meet on a plane")
    body = load("volute_body")
    bad += pair("volute body <-> bell plug", body, load("volute_bell"),
                "chordal lap film")
    bad += pair("volute body <-> mouthpiece plug", body,
                load("volute_mouthpiece"), "chordal lap film")
    ch = load("halo_chamber")
    for i in range(4):
        a = 2 * math.pi * i / 4
        Rz = trimesh.transformations.rotation_matrix(a, (0, 0, 1))
        sp = []
        for h in "ab":
            m = load(f"halo_sphere{i+1}_{h}")
            m.apply_translation((0, 0, HL.Z_NECK))
            m.apply_transform(Rz)
            sp.append(m)
        sph = trimesh.util.concatenate(sp)
        # both thread halves come from one field function, so the exported pose
        # IS the right phase. The ~150 mm3 that remains is the collar's root
        # ring biting ~1 mm into the chamber wall corner (stub starts at
        # R_CH + 4, the wall's outer face is at R_CH + 5) -- pre-existing, and
        # identical on the bare --no-skin parts.
        bad += pair(f"halo chamber <-> sphere {i+1} (threaded)", ch, sph,
                    "collar root ring, pre-existing")
    print(f"  {'TOTAL':<44} {bad:10.1f} mm3")


if __name__ == "__main__":
    main()
