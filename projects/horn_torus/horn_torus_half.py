"""
Cut the horn torus in half on a VERTICAL plane through the central axis and
cap the cut face, so you can look straight at the cross-section.

For a horn torus that cross-section is two disks tangent at the origin -- the
"connecting center point" where the tube meets itself. (A ring torus would
show two separated disks with a gap; a horn torus shows them kissing.)
"""

import numpy as np
import trimesh

SRC = "horn_torus.stl"
OUT = "horn_torus_half.stl"


def main():
    m = trimesh.load(SRC)

    # keep the half with y >= 0 by intersecting with a big half-space box.
    # the box's y range starts at 0, so the cut face lands on the vertical
    # plane through the central axis.
    span = float(np.abs(m.bounds).max()) * 2.0 + 10.0
    box = trimesh.creation.box(extents=[2 * span, span, 2 * span])
    box.apply_translation([0, span / 2.0, 0])     # box covers y in [0, span]

    half = trimesh.boolean.intersection([m, box])

    if half is None or len(half.faces) == 0:
        raise RuntimeError("intersection produced empty mesh")

    half.fix_normals()
    print(f"half: verts={len(half.vertices)}, faces={len(half.faces)}, "
          f"watertight={half.is_watertight}")
    half.export(OUT)
    sz = np.round(half.bounds[1] - half.bounds[0], 1)
    print(f"exported {OUT}: size={sz} mm")


if __name__ == "__main__":
    main()
