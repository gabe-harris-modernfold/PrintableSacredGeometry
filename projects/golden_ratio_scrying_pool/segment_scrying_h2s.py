"""
Bambu Lab H2S version of the full scrying-pool assembly (pool + collector
collar/horns + stand), split into exactly 2 printable halves.

The H2S build volume is 340 x 320 x 340 mm. The full-size assembly is
Ø761 x 123mm -- far too big for any single piece, and even a half of it
still spans the FULL diameter along its flat cut face (bisecting a disc
does not shrink its rim-to-rim distance), so a 2-piece split only fits the
H2S if the whole model is scaled down first.

H2S_SCALE = 0.33 shrinks the assembly to Ø302 x 107mm (comfortable margin
under the 340mm bed limit; wall thickness is NOT scaled -- it's already an
absolute print-strength constant, ~3.24mm -- so it does not get thinner).

CUT: a single vertical plane through the center, at the angle that falls
exactly between two collar cups on one side (avoiding a random asymmetric
slice). Because N_CUPS=21 is odd, cup and horn positions together occupy
every pi/21 angular step around the collar, so no diameter can dodge every
feature on both sides -- the opposite ray of this same cut line runs exactly
through the axis of one horn tunnel, splitting it into a lengthwise groove
in each half (rather than a lopsided partial cut), and through the center
of one cup, splitting it cleanly in half.

The cut itself is done as a BOOLEAN INTERSECTION with a half-space box
(not trimesh's slice_plane cap): the horn tunnel's lengthwise cut leaves an
annular (hole-in-hole) cross-section that slice_plane's simple planar
capping cannot close (verified: 3 broken faces, euler off by ~1), while the
CSG boolean engine resolves it correctly (verified watertight both halves).

Each half is then rotated 90 degrees so its flat cut face sits on the print
bed (print orientation), rather than standing the un-cut height on edge.
"""

import numpy as np
import trimesh

import scrying_pool as sp
import add_collectors as ac
import add_stand as st

# ---- parameters ------------------------------------------------------------
H2S_SCALE  = 0.33                 # golden-ladder scale factor (see docstring)
BED_X      = 340.0                # H2S build volume (mm)
BED_Y      = 320.0
BED_Z      = 340.0
CUT_ANGLE  = np.pi / ac.N_CUPS     # between cup0 and cup1 (see docstring)
PAD        = 50.0                 # half-space box overhang beyond the model (mm)
OUT_A      = "scrying_pool_h2s_half_A.stl"
OUT_B      = "scrying_pool_h2s_half_B.stl"
# -----------------------------------------------------------------------------


def half_space_box(lo, hi, side):
    """Axis-aligned box covering y>=0 (side=+1) or y<=0 (side=-1) of the
    model's bounds, padded well beyond every other extent."""
    ext = [hi[0] - lo[0] + 2 * PAD, (hi[1] - lo[1]) / 2 + PAD, hi[2] - lo[2] + 2 * PAD]
    box = trimesh.creation.box(extents=ext)
    cx, cz = (lo[0] + hi[0]) / 2, (lo[2] + hi[2]) / 2
    cy = side * ((hi[1] - lo[1]) / 4 + PAD / 2)
    box.apply_translation([cx, cy, cz])
    return box


def cut_half(assembly, side):
    lo, hi = assembly.bounds
    box = half_space_box(lo, hi, side)
    half = trimesh.boolean.intersection([assembly, box])
    half.merge_vertices()
    half.fix_normals()
    return half


def reorient_flat(mesh):
    """Rotate the cut face (currently at y=0, facing -y) onto the bed
    (facing -z), then drop the piece onto z=0."""
    R = trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0])
    mesh.apply_transform(R)
    mesh.apply_translation([0, 0, -mesh.bounds[0][2]])
    return mesh


def footprint(mesh):
    lo, hi = mesh.bounds
    return hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]


def main():
    sp.SCALE = H2S_SCALE
    assembly, info = st.build()

    # rotate so the chosen cut direction aligns with +X; the cut plane
    # becomes exactly y=0, so the half-space boxes are axis-aligned
    Rz = trimesh.transformations.rotation_matrix(-CUT_ANGLE, [0, 0, 1])
    assembly.apply_transform(Rz)

    half_a = cut_half(assembly, +1)
    half_b = cut_half(assembly, -1)

    half_a = reorient_flat(half_a)
    half_b = reorient_flat(half_b)

    half_a.export(OUT_A)
    half_b.export(OUT_B)

    print(f"golden scale             : {H2S_SCALE}")
    print(f"full assembly diameter   : {info['collar_dia']:.1f} mm  "
          f"(unscaled would be 761.4 mm)")
    print(f"H2S bed                  : {BED_X:.0f} x {BED_Y:.0f} x {BED_Z:.0f} mm")
    print(f"cut angle                : {np.degrees(CUT_ANGLE):.3f} deg "
          f"(between cup 0/1; splits cup 11 & horn 0 lengthwise on the far side)")
    print()
    for name, half in (("half_A", half_a), ("half_B", half_b)):
        fx, fy, fz = footprint(half)
        fits = fx <= BED_X and fy <= BED_Y and fz <= BED_Z
        print(f"{name:8s} footprint {fx:6.1f} x {fy:6.1f} mm   "
              f"print height {fz:6.1f} mm   fits_H2S={fits}   "
              f"watertight={half.is_watertight}")


if __name__ == "__main__":
    main()
