"""
Split the scrying pool into printable pieces, each within a bed limit.

A Ø702 mm pool cannot be printed whole on a ~310 mm bed, and it cannot be done
in 3 pieces: the radius (351 mm) already exceeds the bed, so every piece that
spans center-to-rim is too big. Minimum feasible split:

    1 CENTER DISC   cut on the golden base circle  (Ø = u*phi^2 = 261.8 mm)
    N RIM WEDGES    equal arcs of the flaring wall, chord <= bed limit

The center is isolated with a cylinder boolean; each wedge is carved from the
rim ring by two radial plane slices (capped -> watertight solids). Every part
is checked so its bounding box fits the bed.
"""

import os
import numpy as np
import trimesh

import scrying_pool as sp

# ---- parameters ----------------------------------------------------------
BED     = 310.0     # max piece footprint (mm)
NW      = 7         # number of rim wedges
OUTDIR  = "scrying_pool_segments"
# --------------------------------------------------------------------------


def foot(mesh):
    """Footprint (x,y bounding box) in mm."""
    lo, hi = mesh.bounds
    return hi[0] - lo[0], hi[1] - lo[1]


def main():
    vessel, info = sp.build()
    rc = 0.5 * info["base_dia"]          # golden base radius = cut circle
    H = info["height"]

    os.makedirs(OUTDIR, exist_ok=True)
    cyl = trimesh.creation.cylinder(radius=rc, height=4 * H, sections=sp.SECTIONS)

    center = trimesh.boolean.intersection([vessel, cyl])
    ring = trimesh.boolean.difference([vessel, cyl])

    pieces = [("center_disc", center)]
    for k in range(NW):
        a, b = 2 * np.pi * k / NW, 2 * np.pi * (k + 1) / NW
        w = ring.slice_plane([0, 0, 0], [-np.sin(a), np.cos(a), 0], cap=True)
        w = w.slice_plane([0, 0, 0], [np.sin(b), -np.cos(b), 0], cap=True)
        w.merge_vertices(); w.fix_normals()
        pieces.append((f"wedge_{k+1}", w))

    print(f"cut circle (golden base) : Ø{2*rc:.2f} mm")
    print(f"{'piece':14s} {'footprint(mm)':>18s}  {'fits '+str(int(BED)):>8s}  "
          f"{'watertight':>10s}")
    allfit = True
    for name, m in pieces:
        fx, fy = foot(m)
        fit = max(fx, fy) <= BED
        allfit &= fit
        m.export(os.path.join(OUTDIR, name + ".stl"))
        print(f"{name:14s} {fx:8.1f} x {fy:6.1f}  {str(fit):>8s}  "
              f"{str(m.is_watertight):>10s}")
    print(f"\n{len(pieces)} pieces -> {OUTDIR}/   all within {int(BED)} mm: {allfit}")


if __name__ == "__main__":
    main()
