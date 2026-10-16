#!/usr/bin/env python3
"""Hanging eye for the Negative Space sculpture.

A flat plate printed lying down, so every layer runs around the loop and the
load never pulls layers apart. Its tongue drops into a slot sawn across the top
of the 20 mm axis tube and is held by one M4 bolt through tube and tongue. Load
path: ribbon -> spokes -> tube -> M4 bolt (double shear) -> tongue -> eye.

Writes hanger_ring.stl. Part coordinates: plate in XY, tongue pointing -Y, bolt
hole centre at the origin, thickness along +Z."""

import math

import trimesh
from shapely.geometry import Point, box

THICK = 8.0                     # print height; saw the tube slot 8.4 mm wide
RING_OD, RING_ID = 56.0, 36.0   # 36 mm eye takes a ceiling hook or a cord loop
TUBE_OD = 20.0                  # tongue is as wide as the tube, flush in the slot
TONGUE_LEN = 32.0               # slot depth in the tube
BOLT_D = 4.3                    # M4 clearance
BOLT_EDGE = 10.0                # bolt centre to the tongue's lower end (tear-out edge)
NECK_W, NECK_LEN = 16.0, 6.0
ROUND = 2.0                     # corner radius in plan
PETG_G_MM3 = 1.27e-3
DESIGN_LOAD_N = 5.0 * 9.81      # whole sculpture, rounded up


def outline():
    tongue_top = TONGUE_LEN - BOLT_EDGE
    yc = tongue_top + NECK_LEN + 0.5 * RING_OD
    ring = Point(0, yc).buffer(0.5 * RING_OD, quad_segs=64)
    plate = (box(-0.5 * TUBE_OD, -BOLT_EDGE, 0.5 * TUBE_OD, tongue_top)
             .union(box(-0.5 * NECK_W, tongue_top - 1.0, 0.5 * NECK_W, yc))
             .union(ring))
    plate = plate.buffer(ROUND).buffer(-2 * ROUND).buffer(ROUND)   # round every corner
    return (plate.difference(Point(0, yc).buffer(0.5 * RING_ID, quad_segs=64))
            .difference(Point(0, 0).buffer(0.5 * BOLT_D, quad_segs=32))), yc


def main():
    poly, yc = outline()
    mesh = trimesh.creation.extrude_polygon(poly, THICK)
    mesh.export("hanger_ring.stl")

    # Eye as a thin ring under a diametral pull: M_max = P R / pi at the load point.
    rm = 0.25 * (RING_OD + RING_ID)
    b = 0.5 * (RING_OD - RING_ID)
    bend = DESIGN_LOAD_N * rm / math.pi / (THICK * b ** 2 / 6.0)
    tension = DESIGN_LOAD_N / (2.0 * THICK * b)
    tear = DESIGN_LOAD_N / (2.0 * THICK * (BOLT_EDGE - 0.5 * BOLT_D))
    print(f"hanger_ring.stl: extents {mesh.extents.round(1)} mm, watertight "
          f"{mesh.is_watertight}, {mesh.volume * PETG_G_MM3:.0f} g PETG")
    print(f"  eye centre {yc:.0f} mm above the bolt; at {DESIGN_LOAD_N:.0f} N the ring sees "
          f"{bend + tension:.1f} MPa, bolt tear-out shear {tear:.2f} MPa")


if __name__ == "__main__":
    main()
