#!/usr/bin/env python3
"""
Absorbing-edge variant of the cymatics dish.

Same Ø277 dish, same skirt, same support ring -- the only change is a perforated
fence standing just inside the rim, with a narrow chamber behind it. Waves that
reach the edge push flow through the slots; the contraction and re-expansion
sheds vortices and burns the energy instead of reflecting it.

WHY THIS AND NOT A BEACH
------------------------
A beach (depth tapering to zero) is the textbook absorber, but it needs 2-3
wavelengths of run -- 30 mm here -- and it cannot be printed integral to a
floor-down dish anyway (a 5.7 degree ramp cantilevers at 84 degrees from
vertical). The fence does the same job in 1.8 mm of radius.

SIZING
------
Every number below came out of cymatics_fence_design.py, which builds the
impedance of the fence (slot inertance + linearised jet loss) and the chamber
(transmission line closed by the rim) against the characteristic impedance of a
capillary-gravity wave at this depth, then sweeps for the best worst-case
absorption over 10-120 Hz. Do not tune these by eye -- rerun that script.

    fence     0.8 mm thick, 1.2 mm slots, porosity 0.45
    chamber   1.0 mm
    predicted 54% absorbed at the weakest frequency, 75% mean

CAVEAT
------
The chamber is 1.0 mm wide and the capillary length in water is 2.7 mm, so the
chamber will wick full and its free surface is a strongly curved meniscus rather
than the flat channel the model assumes. That makes it STIFFER than modelled,
which is equivalent to a smaller chamber -- and smaller still works (0.75 mm
scores 46%/83%). The error lands on the safe side of the cliff at B = 1.75 mm,
where absorption collapses. If it underperforms, move CHAMBER, not the fence.
"""

import numpy as np
import trimesh

import cymatics_dish as C

# ---- fence, from cymatics_fence_design.py --------------------------------
FENCE_T  = 0.8             # fence wall thickness (mm) -- 2 perimeters at 0.4
SLOT_W   = 1.2             # slot width (mm)
POROSITY = 0.45            # open fraction of the fence
CHAMBER  = 1.0             # gap between the fence and the rim's inner face (mm)
FENCE_H  = 6.5             # fence height (mm)
SLOT_H   = 5.0             # slot height (mm) -- leaves a lintel tying the posts
SLOT_SINK = 0.05           # cut this far into the floor so the slot floor is
                           # flush and the boolean has no coplanar faces

OUT = "cymatics_dish_absorbing.stl"

R_FO = C.R_IN - CHAMBER            # fence outer face
R_FI = R_FO - FENCE_T              # fence inner face -- the water surface edge
PITCH = SLOT_W / POROSITY
N_SLOT = int(round(2 * np.pi * (R_FO + R_FI) / 2 / PITCH))

# The slots must open into the chamber and STOP there. Run them to the dish's
# outer radius and they perforate the rim itself, which drains the dish onto the
# speaker -- so the cut is bounded well short of the rim's inner face.
SLOT_R_IN  = R_FI - 1.0                        # inboard of the fence, into water
SLOT_R_OUT = R_FO + CHAMBER / 2                # into the chamber, short of the rim

assert FENCE_T <= 1.5, "fence thicker than asked"
assert SLOT_H < FENCE_H, "slots would cut the fence into loose posts"
assert R_FI > C.SUPPORT_D / 2 + 10, "fence would crowd the support ring"
assert SLOT_R_OUT < C.R_IN - 0.2, "slot cut would reach the rim wall and leak"
assert SLOT_R_IN < R_FI, "slot cut would not reach through the fence"


def fence_ring():
    """Plain tube standing on the floor, before the slots are cut."""
    return C.solid(np.array([[R_FI, 0.0], [R_FO, 0.0],
                             [R_FO, FENCE_H], [R_FI, FENCE_H]]))


def slot_cutter():
    """N_SLOT radial boxes through the fence, as one solid."""
    boxes = []
    for i in range(N_SLOT):
        b = trimesh.creation.box(extents=[SLOT_R_OUT - SLOT_R_IN, SLOT_W,
                                          SLOT_H + SLOT_SINK])
        b.apply_translation([(SLOT_R_IN + SLOT_R_OUT) / 2, 0.0,
                             (SLOT_H - SLOT_SINK) / 2])
        b.apply_transform(trimesh.transformations.rotation_matrix(
            2 * np.pi * i / N_SLOT, [0, 0, 1]))
        boxes.append(b)
    return trimesh.util.concatenate(boxes)


def main():
    dish = C.solid(C.dish_profile()).union(fence_ring())
    dish = dish.difference(slot_cutter())
    dish.fix_normals()
    dish.export(OUT)

    water = np.pi * R_FI**2 * C.WATER * 1e-3
    post = PITCH - SLOT_W
    print(f"base dish    : Ø{C.OD:.0f}, unchanged rim and floor")
    print(f"fence        : Ø{2 * R_FI:.1f} inner / Ø{2 * R_FO:.1f} outer, "
          f"{FENCE_T:.1f} mm thick, {FENCE_H:.1f} mm tall")
    print(f"slots        : {N_SLOT} × {SLOT_W:.1f} mm on a {PITCH:.2f} mm pitch "
          f"({POROSITY:.0%} open), {SLOT_H:.1f} mm tall, {post:.2f} mm posts")
    print(f"chamber      : {CHAMBER:.1f} mm to the rim's inner face")
    print(f"slot cut     : r {SLOT_R_IN:.1f} -> {SLOT_R_OUT:.1f}, stopping "
          f"{C.R_IN - SLOT_R_OUT:.1f} mm short of the rim at r={C.R_IN:.1f}")
    print(f"radial cost  : {CHAMBER + FENCE_T:.1f} mm per side - water surface "
          f"Ø{2 * C.R_IN:.0f} -> Ø{2 * R_FI:.1f}")
    print(f"water        : {water:.0f} ml at {C.WATER:.0f} mm "
          f"(was {np.pi * C.R_IN**2 * C.WATER / 1000:.0f} ml)")
    e = dish.extents
    print(f"part         : Ø{max(e[0], e[1]):.1f} x {e[2]:.1f} mm   "
          f"{dish.volume * C.DENSITY:.1f} g   "
          f"watertight={dish.is_watertight}  faces={len(dish.faces)}")
    print(f"exported {OUT}")


if __name__ == "__main__":
    main()
