#!/usr/bin/env python3
"""Incense vortex-ring base -- a hand-pumped smoke-ring generator for the 100 mm acrylic tube.

How it works
------------
An incense cone smoulders in a well sunk into the floor of a closed smoke chamber. The
chamber floor's centre is a latex membrane (a cut party balloon) stretched over a spigot.
Pressing the piano-key lever pushes the membrane up; the displaced slug of smoky air leaves
through a sharp-edged orifice under the tube, rolls up into a vortex ring and climbs the
tube. Releasing the key lets the membrane spring back, and the chamber refills through a
nitrile-film flap lying on a grille in the chamber floor, or back through the orifice.
No electronics.

Design numbers (see README.md)
------------------------------
* Ring diameter ~1.2x orifice; keep it under ~half the 96 mm bore -> orifice 25-40 mm,
  30 mm default. Wall image vortex slows a 36 mm ring only ~11 %.
* Clean rings want a slug length L = 2-4 orifice diameters. The Ø110 membrane pushed by an
  Ø80 disc displaces ~7150 mm^3 per mm of travel, so the 12 mm full stroke gives L/D = 4
  through the 30 mm orifice. The inserts have a bellmouth inlet, so the bore runs full, and a
  square exit edge, so the ring separates at the bore diameter: both keep that L/D honest.
* A soft ~0.1 s press gives rings rising ~0.4 m/s (~1 s up the tube); a sharp tap fires
  them at ~3 m/s.
* The cone's heat makes chamber + tube a chimney: the Hearth tornado runs on exactly this
  draft, ~0.1-0.3 Pa. Any opening low in the chamber would feed a steady smoke stream up the
  tube, so the well opens only into the chamber (four ducts bring chamber air down to the
  cone's base), and the refill valve is a 0.1 mm nitrile film lying flat on a floor grille:
  its own weight (~1.1 Pa) holds it shut against the draft, and a brisk key release lifts it.
  What remains is a faint buoyant exchange through the orifice itself (~4-6 mL/s).
* The key stops on a pin between two front posts rather than on the table, so the 12 mm
  stroke holds at a table edge or on a cloth.

Run from this directory:  python vortex_ring_base.py
Outputs: one STL per part in print orientation, plus preview PNGs.
"""

import math
import os
import sys

import numpy as np
import trimesh
from build123d import (Align, Axis, Box, Circle, Cylinder, Plane, Polygon, Pos, Rectangle,
                       RectangleRounded, Rot, Text, export_stl, extrude, revolve)

HERE = os.path.dirname(os.path.abspath(__file__))
C, MIN = Align.CENTER, Align.MIN

# ---------------------------------------------------------------- parameters

# tube interface -- identical to the Hearth / Wind bases, so one tube serves all of them
TUBE_OD, TUBE_WALL, TUBE_LEN = 100.0, 2.0, 393.7
FIT = 0.4                                   # radial slip clearance
GROOVE_IR = TUBE_OD / 2 - TUBE_WALL - FIT   # 47.6
GROOVE_OR = TUBE_OD / 2 + FIT               # 50.4
COLLAR_W, COLLAR_H = 2.4, 10.0

# smoke chamber
CH_R = 114.0            # chamber inner radius
WALL_T = 2.4            # chamber wall thickness
CH_H = 90.0             # floor top -> lid underside (keeps the lid ~90 mm above the ember)
PLATE_T = 5.0           # floor plate and lid thickness
FLOOR_GROOVE_D = 2.5
LID_GROOVE_D = 3.0
GROOVE_IN = CH_R - FIT                  # 113.6  wall-ring groove, inner radius
GROOVE_OUT = CH_R + WALL_T + FIT        # 116.8  wall-ring groove, outer radius
OUTER_R = GROOVE_OUT + 2.4              # 119.2  floor plate / lid radius

# membrane and pusher
MEM_R = 55.0            # spigot bore = free membrane radius (Ø110)
SPIG_T = 3.5            # spigot wall: leaves 2 mm of wall under the collar groove
SPIG_H = 14.0
GRV_D = 1.5             # spigot groove that the collar's bead presses the latex into
COL_GAP = 0.6           # collar bore clearance over the spigot, for the rolled-up latex
COL_T = 3.4             # collar wall
COL_SPLIT = 4.0         # gap the zip tie closes
DISC_R = 40.0           # pusher disc (Ø80)
DISC_T = 4.0
PRELOAD = 2.0           # disc lifts the membrane this far at rest -> keeps the key returned
STROKE = 12.0           # disc travel when the lever lands on the stop pin

# piano-key lever (runs along -Y)
LEVER_HW = 8.0          # half width
Y_PIV = -76.0           # pivot pin
Y_KEY = -184.0          # key end
ROD_R = 3.0             # cross rod under the disc
ROD_FIT = 0.1           # rod hole clearance in the lever; the rod is glued
PIN_R = 2.5             # pivot pin and stop pin
PIN_FIT = 0.15          # slip fit in the cheeks and posts: the pins come out for a balloon change
SADDLE_DROP = 6.5       # cross-rod centre -> disc underside
CHEEK_X0, CHEEK_X1 = LEVER_HW + 0.4, LEVER_HW + 4.4

# front posts: the lever passes between them, and a pin across them under its arm ends the stroke
POST_X0, POST_X1 = LEVER_HW + 1.8, LEVER_HW + 13.8
POST_Y0, POST_Y1 = -150.0, -127.0
Y_STOP = -134.0

# feet, out at the rim (with the posts they set the tipping margin)
FEET = [(sx * 73.5, sy * 73.5) for sx in (-1, 1) for sy in (-1, 1)]
FOOT_R0, FOOT_R1 = 13.5, 11.5           # radius at the plate, on the table

# cone well (behind the membrane, on +Y)
WELL_Y = 88.0
WELL_R = 20.0           # inner radius (Ø40 takes any bottle cap)
WELL_T = 2.4
CAP_R = 15.0            # the cap seats where the 45° funnel narrows to this radius
CAP_AND_CONE = 28.0     # cap (~6) + cone (~22): puts the cone tip at floor level
APEX_R = 4.0            # flat at the bottom of the funnel
DUCT_R = 2.0            # air ducts from the chamber floor down to the cone's base
DUCT_OFF = WELL_R + 4.6
DUCT_ANGS = (45, 135, 225, 315)

# lid hatch over the well (the two nylon hinges go here)
HATCH_HOLE_R = 21.0
HATCH_W, HATCH_L = 50.0, 48.0           # tangential (x) x radial (y)
HINGE_GAP = 4.0                         # knuckle clearance between hatch and hinge land

# refill flap valve: nitrile film lying on a grille in the floor, hinge edge pinned by a clip
FLAP_X = -89.0                          # grille centre, on -X: clear of the well, lever and feet
FLAP_SLOTS = (-12.5, 0.0, 12.5)         # slot centres along X, relative to FLAP_X
FLAP_SLOT_W, FLAP_SLOT_L = 9.0, 30.0    # slots run along Y
CLIP_X = -64.0                          # clip strip along the film's inner edge
CLIP_PEGS = (-18.0, 18.0)

# orifice inserts
ORIFICES = (25, 30, 35, 40)
INS_FLANGE_R, INS_BODY_R = 35.0, 28.0
REENTRANT = 12.0        # insert collar hanging below the lid
BELL_R = 6.0            # bellmouth radius at the insert's inlet

# ---------------------------------------------------------------- derived heights

# Heights are set so the key's far bottom corner would just reach the table at full stroke. The
# key end is chamfered to clear the table by ~1.3 mm there; the stop pin takes the load instead.
TH_FULL = math.asin(STROKE / (0.0 - Y_PIV))
Z_PIV = (Y_PIV - Y_KEY) * math.sin(TH_FULL) + 8.0 * math.cos(TH_FULL)
Z_ROD = Z_PIV
Z_DISC_BOT = Z_ROD + SADDLE_DROP
Z_DISC_TOP = Z_DISC_BOT + DISC_T
Z_MEM = Z_DISC_TOP - PRELOAD            # spigot rim = membrane clamp plane
Z_F = Z_MEM + SPIG_H                    # floor plate underside
Z_FT = Z_F + PLATE_T                    # floor plate top
Z_LB = Z_FT + CH_H                      # lid underside
Z_LT = Z_LB + PLATE_T                   # lid top (tube rests here)
Z_CAP = Z_FT - CAP_AND_CONE             # cap seat height in the well
Z_FUN_TOP = Z_CAP + (WELL_R - CAP_R)    # funnel starts (45°)
Z_BLEED = Z_CAP + 10.5                  # duct outlets into the well, just above the cap
Z_GRV = Z_MEM + 5.0                     # spigot groove starts
# The stop pin rests on the bottom of its hole; the lever's flat underside (8 mm below the pivot
# axis) is tangent to it exactly at full stroke.
Z_STOP = Z_PIV + (-8.0 - PIN_R + math.sin(TH_FULL) * (Y_STOP - Y_PIV)) / math.cos(TH_FULL)

# ---------------------------------------------------------------- helpers


def cyl(r, z0, z1, x=0.0, y=0.0):
    return Pos(x, y, z0) * Cylinder(r, z1 - z0, align=(C, C, MIN))


def ring(ri, ro, z0, z1, x=0.0, y=0.0):
    return cyl(ro, z0, z1, x, y) - cyl(ri, z0 - 1, z1 + 1, x, y)


def box(x0, x1, y0, y1, z0, z1):
    return Pos((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2) * Box(x1 - x0, y1 - y0, z1 - z0)


def rev(pts, arc=360.0):
    """Revolve an (r, z) profile about Z."""
    return revolve(Plane.XZ * Polygon(*pts, align=None), Axis.Z, revolution_arc=arc)


def xrod(r, x0, x1, y, z):
    """Cylinder along X."""
    return Pos((x0 + x1) / 2, y, z) * Rot(0, 90, 0) * Cylinder(r, x1 - x0)


def radial_rod(r, r0, r1, ang, z):
    """Cylinder along a radial direction at azimuth `ang` (deg), spanning radius r0..r1."""
    a = math.radians(ang)
    rm = (r0 + r1) / 2
    return Pos(rm * math.cos(a), rm * math.sin(a), z) * Rot(0, 0, ang) * Rot(0, 90, 0) * Cylinder(r, r1 - r0)


def arc_pts(cx, cz, r, a0, a1, n=10):
    return [(cx + r * math.cos(math.radians(a)), cz + r * math.sin(math.radians(a)))
            for a in np.linspace(a0, a1, n)]

# ---------------------------------------------------------------- parts (assembled coordinates)


def floor_plate():
    outline = Circle(OUTER_R) + Pos(0, -125) * RectangleRounded(84, 50, 8)   # tongue over the key
    part = extrude(Plane.XY.offset(Z_F) * outline, PLATE_T)
    part -= ring(GROOVE_IN, GROOVE_OUT, Z_FT - FLOOR_GROOVE_D, Z_FT + 1)
    part -= cyl(MEM_R, Z_F - 1, Z_FT + 1)

    # membrane spigot: rounded rim so the latex can't tear; the collar's bead presses the
    # rolled-up latex into the groove, so the membrane can't creep down off the rim
    ro = MEM_R + SPIG_T
    spig = ([(MEM_R, Z_F + 0.5), (MEM_R, Z_MEM + SPIG_T / 2)]
            + arc_pts(MEM_R + SPIG_T / 2, Z_MEM + SPIG_T / 2, SPIG_T / 2, 180, 360)
            + [(ro, Z_GRV), (ro - GRV_D, Z_GRV + GRV_D), (ro - GRV_D, Z_GRV + GRV_D + 1.0),
               (ro, Z_GRV + 2 * GRV_D + 1.0), (ro, Z_F + 0.5)])
    part += rev(spig)

    # cone well: straight bore and a closed 45° funnel the cap self-centres in. It opens only
    # into the chamber -- a vent to the room under the floor would turn chamber + tube into a
    # chimney. Four ducts, carried in ribs on the well wall, bring chamber air from the floor
    # down to the cone's base instead.
    v = 2 * math.sqrt(2.0)                                  # 2 mm wall measured normal to 45°
    zo = Z_FUN_TOP - WELL_R - v                             # outer funnel: z = r + zo
    z_floor = Z_FUN_TOP - (WELL_R - APEX_R)                 # flat well floor (no cone apex)
    part += Pos(0, WELL_Y, 0) * rev([(0, Z_F + 0.5), (WELL_R + WELL_T, Z_F + 0.5),
                                     (WELL_R + WELL_T, WELL_R + WELL_T + zo),
                                     (z_floor - 2.0 - zo, z_floor - 2.0), (0, z_floor - 2.0)])
    ducts = [(DUCT_OFF * math.cos(math.radians(a)), WELL_Y + DUCT_OFF * math.sin(math.radians(a)))
             for a in DUCT_ANGS]
    for x, y in ducts:
        part += cyl(DUCT_R + 2.0, Z_BLEED - DUCT_R - 2.0, Z_F + 0.5, x, y)
    part -= Pos(0, WELL_Y, 0) * rev([(0, Z_FT + 1), (WELL_R, Z_FT + 1), (WELL_R, Z_FUN_TOP),
                                     (APEX_R, z_floor), (0, z_floor)])
    for a, (x, y) in zip(DUCT_ANGS, ducts):
        part -= cyl(DUCT_R, Z_BLEED, Z_FT + 1, x, y)
        part -= Pos(0, WELL_Y, 0) * radial_rod(DUCT_R, WELL_R - 1.0, DUCT_OFF + 1.0, a, Z_BLEED)

    # refill flap valve: grille through the floor, blind peg holes for the film's clip strip
    for dx in FLAP_SLOTS:
        part -= box(FLAP_X + dx - FLAP_SLOT_W / 2, FLAP_X + dx + FLAP_SLOT_W / 2,
                    -FLAP_SLOT_L / 2, FLAP_SLOT_L / 2, Z_F - 1, Z_FT + 1)
    for y in CLIP_PEGS:
        part -= cyl(1.6, Z_FT - 3.5, Z_FT + 1, CLIP_X, y)

    # lever fulcrum: two cheeks, rest stop bar outboard of the pin
    for s in (-1, 1):
        part += box(*sorted((s * CHEEK_X0, s * CHEEK_X1)), -90, -65, Z_PIV - 7, Z_F + 0.5)
    part += box(-CHEEK_X1, CHEEK_X1, -90, -84, Z_PIV + 8, Z_F + 0.5)
    part -= xrod(PIN_R + PIN_FIT, -20, 20, Y_PIV, Z_PIV)

    # front posts under the tongue: the lever runs between them, the stop pin across them
    for s in (-1, 1):
        part += box(*sorted((s * POST_X0, s * POST_X1)), POST_Y0, POST_Y1, 0, Z_F + 0.5)
    part -= xrod(PIN_R + PIN_FIT, -(POST_X1 + 1), POST_X1 + 1, Y_STOP, Z_STOP + PIN_FIT)

    # feet out at the rim, tapered for a stiff root
    for x, y in FEET:
        part += Pos(x, y, 0) * rev([(0, 0), (FOOT_R1, 0), (FOOT_R0, Z_F + 0.5), (0, Z_F + 0.5)])
    return part


def chamber_wall():
    return ring(CH_R, CH_R + WALL_T, Z_FT - FLOOR_GROOVE_D, Z_LB + LID_GROOVE_D)


def valve_clip():
    """Strip that pins the hinge edge of the nitrile film flat on the chamber floor."""
    part = extrude(Plane.XY.offset(Z_FT) * Pos(CLIP_X, 0) * RectangleRounded(7.0, 44.0, 2.0), 2.5)
    for y in CLIP_PEGS:
        part += cyl(1.45, Z_FT - 3.0, Z_FT + 0.5, CLIP_X, y)
    return part


def lid():
    part = extrude(Plane.XY.offset(Z_LB) * Circle(OUTER_R), PLATE_T)
    part -= ring(GROOVE_IN, GROOVE_OUT, Z_LB - 1, Z_LB + LID_GROOVE_D)
    part += ring(GROOVE_IR - COLLAR_W, GROOVE_IR, Z_LT - 0.5, Z_LT + COLLAR_H)
    part += ring(GROOVE_OR, GROOVE_OR + COLLAR_W, Z_LT - 0.5, Z_LT + COLLAR_H)
    part -= cyl(INS_FLANGE_R + 0.2, Z_LT - 3.0, Z_LT + 1)         # insert counterbore
    part -= cyl(INS_BODY_R + 0.2, Z_LB - 1, Z_LT)
    part -= cyl(HATCH_HOLE_R, Z_LB - 1, Z_LT + 1, y=WELL_Y)
    # raised hinge land beside the hatch: its top is flush with the closed hatch, so a plain
    # butt hinge screws flat across both
    land = box(HATCH_W / 2 + HINGE_GAP, HATCH_W / 2 + HINGE_GAP + 20,
               WELL_Y - HATCH_L / 2, WELL_Y + HATCH_L / 2, Z_LT - 0.5, Z_LT + 3.0)
    part += land & cyl(OUTER_R - 2, Z_LT - 1, Z_LT + 4)
    return part


def hatch():
    part = extrude(Plane.XY.offset(Z_LT) * Pos(0, WELL_Y) * RectangleRounded(HATCH_W, HATCH_L, 5), 3.0)
    part += ring(HATCH_HOLE_R - 3.0, HATCH_HOLE_R - 0.4, Z_LT - 2.9, Z_LT + 0.5, y=WELL_Y)
    part += box(-HATCH_W / 2 - 9, -HATCH_W / 2 + 1, WELL_Y - 7, WELL_Y + 7, Z_LT, Z_LT + 3.0)  # lift tab
    return part


def orifice_insert(d):
    """Straight bore with a square exit edge flush with the top face, so the slug separates at
    the bore diameter, and a bellmouth inlet, so the bore runs full instead of necking down."""
    r = d / 2
    z_bot = Z_LB - REENTRANT
    prof = ([(INS_BODY_R, z_bot), (INS_BODY_R, Z_LT - 3.0), (INS_FLANGE_R, Z_LT - 3.0),
             (INS_FLANGE_R, Z_LT), (r, Z_LT), (r, z_bot + BELL_R)]
            + arc_pts(r + BELL_R, z_bot + BELL_R, BELL_R, 180, 270, 12)[1:])
    part = rev(prof)
    label = Plane.XY.offset(Z_LT) * Pos(0, 29.0) * Text(str(d), font_size=5.5, align=(C, C))
    part -= extrude(label, amount=-0.6)
    return part


def pusher_disc(dy=0.0, dz=0.0):
    prof = ([(0, Z_DISC_BOT), (DISC_R, Z_DISC_BOT), (DISC_R, Z_DISC_TOP - 2.0)]
            + arc_pts(DISC_R - 2.0, Z_DISC_TOP - 2.0, 2.0, 0, 90) + [(0, Z_DISC_TOP)])
    part = rev(prof)
    for s in (-1, 1):                                        # saddles that ride the cross rod
        part += box(*sorted((s * (LEVER_HW + 0.5), s * (LEVER_HW + 12.4))), -6, 6, Z_ROD, Z_DISC_BOT + 0.5)
    part -= xrod(ROD_R + 0.3, -(LEVER_HW + 10.4), LEVER_HW + 10.4, 0, Z_ROD)
    return Pos(0, dy, dz) * part


def membrane_clamp():
    """Split collar: slides up over the latex-covered spigot, its bead seats the latex in the
    spigot groove, and a zip tie (or rubber band) in its outer groove closes the split. The zip
    tie keeps pulling, so the grip doesn't depend on print accuracy."""
    ro = MEM_R + SPIG_T
    ri, rb, rout = ro + COL_GAP, ro - COL_GAP, ro + COL_GAP + COL_T
    z0, z1 = Z_MEM + 2.0, Z_MEM + 12.0
    part = rev([(ri, z0), (rout, z0), (rout, Z_MEM + 4.0), (rout - 1.0, Z_MEM + 5.0),
                (rout - 1.0, Z_MEM + 9.0), (rout, Z_MEM + 10.0), (rout, z1), (ri + 0.8, z1),
                (ri, z1 - 0.8), (ri, Z_GRV + 3.7), (rb, Z_GRV + 2.5), (rb, Z_GRV + 1.5),
                (ri, Z_GRV + 0.3)])
    part -= box(rb - 1, rout + 1, -COL_SPLIT / 2, COL_SPLIT / 2, z0 - 1, z1 + 1)
    return part


def lever():
    zp = Z_PIV
    prof = [(6, 5), (6, zp + 5), (-8, zp + 5), (-11, zp + 2), (-62, zp + 2), (-66, zp + 8),
            (-150, zp + 8), (-154, zp + 11), (-181, zp + 11), (Y_KEY, zp + 8), (Y_KEY, zp - 5),
            (Y_KEY + 8, zp - 8), (-72, zp - 8), (-58, 5)]
    part = extrude(Plane.YZ * Polygon(*prof, align=None), amount=LEVER_HW, both=True)
    part -= xrod(PIN_R + 0.3, -20, 20, Y_PIV, zp)
    part -= xrod(ROD_R + ROD_FIT, -20, 20, 0, Z_ROD)
    return part


def cross_rod():
    return xrod(ROD_R, -(LEVER_HW + 10.0), LEVER_HW + 10.0, 0, Z_ROD)


def pin(x_face):
    """Pin through two supports whose outer faces are at ±x_face, centred on the origin: head on
    -X, and on +X a groove for the hairpin clip."""
    tip = x_face + 4.0
    prof = [(0, -(x_face + 2.6)), (4.0, -(x_face + 2.6)), (4.0, -(x_face + 0.2)),
            (PIN_R, -(x_face + 0.2)), (PIN_R, x_face + 0.2), (PIN_R - 0.8, x_face + 0.2),
            (PIN_R - 0.8, x_face + 2.0), (PIN_R, x_face + 2.0), (PIN_R, tip - 0.5),
            (PIN_R - 0.5, tip), (0, tip)]
    return Rot(0, 90, 0) * rev(prof)          # built along Z, laid along X


def pivot_pin():
    return Pos(0, Y_PIV, Z_PIV) * pin(CHEEK_X1)


def stop_pin():
    return Pos(0, Y_STOP, Z_STOP) * pin(POST_X1)


def pin_clip(x_face=0.0):
    """Hairpin clip: slides sideways onto the pin groove (Ø3.4), its 3.0 mm slot snaps over and
    the pocket holds it. Lies against the support's outer face at x_face, legs pointing -Z."""
    s = Pos(0, -3.25) * RectangleRounded(7.2, 15.5, 1.0)
    s -= Pos(0, -5.75) * Rectangle(3.0, 11.6)
    s -= Circle(1.75)
    s -= Polygon((-2.7, -11.1), (2.7, -11.1), (1.5, -9.4), (-1.5, -9.4), align=None)
    return extrude(Plane.YZ.offset(x_face + 0.3) * s, 1.6)


def pivot_clip():
    return Pos(0, Y_PIV, Z_PIV) * pin_clip(CHEEK_X1)


def stop_clip():
    return Pos(0, Y_STOP, Z_STOP) * Rot(180, 0, 0) * pin_clip(POST_X1)   # legs up: the table is close

# ---------------------------------------------------------------- print orientation + export


def to_bed(part, flip=False, side=False):
    if flip:
        part = Rot(180, 0, 0) * part
    if side:
        part = Rot(0, -90, 0) * part          # X axis -> up; the pin head (-X) lands on the bed
    part = part.clean()                       # merge coplanar faces left by the unions
    bb = part.bounding_box()
    return Pos(-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z) * part


def export(part, name, **orient):
    path = os.path.join(HERE, name)
    export_stl(to_bed(part, **orient), path, tolerance=0.01, angular_tolerance=0.15)
    # Re-save through trimesh: merges the per-facet vertices and turns -0.0 into 0.0, which
    # MeshLab-family slicers otherwise treat as distinct vertices (phantom open seams).
    m = trimesh.load(path)
    m.vertices = m.vertices + 0.0
    m.export(path)
    m = trimesh.load(path)
    bodies = len(m.split(only_watertight=False))
    ext = m.extents
    ok = m.is_watertight and m.volume > 0 and bodies == 1 and max(ext) <= 320
    print(f"  {'OK ' if ok else 'BAD'} {name:24s} {ext[0]:6.1f} x {ext[1]:6.1f} x {ext[2]:6.1f} mm  "
          f"watertight={m.is_watertight} bodies={bodies} vol={m.volume / 1000:6.1f} cm3")
    return ok

# ---------------------------------------------------------------- checks


def rotate_lever(part, deg):
    return part.rotate(Axis((0, Y_PIV, Z_PIV), (1, 0, 0)), deg)


def mechanism_check(parts):
    th = math.degrees(TH_FULL)
    ok = True
    print("\nMechanism clearances (mm) through the stroke:")
    for frac in (0.0, 0.5, 1.0):
        a = th * frac
        lev = rotate_lever(parts["lever"], a)
        rod = rotate_lever(parts["rod"], a)
        rise = -Y_PIV * math.sin(math.radians(a))
        dy = -Y_PIV * (math.cos(math.radians(a)) - 1)
        disc = pusher_disc(dy, rise)
        d_floor = lev.distance_to(parts["floor"])
        d_clamp = lev.distance_to(parts["clamp"])
        d_disc = lev.distance_to(disc)
        d_rod_disc = rod.distance_to(disc)
        d_stop = lev.distance_to(parts["stop"])
        zmin = lev.bounding_box().min.Z
        print(f"  {frac * 100:3.0f}% press (lever {a:4.1f} deg, disc +{rise:4.1f}): "
              f"floor {d_floor:5.2f}  clamp {d_clamp:5.2f}  disc {d_disc:5.2f}  "
              f"rod-disc {d_rod_disc:4.2f}  stop pin {d_stop:5.2f}  key above table {zmin:5.2f}")
        # the stop bar is meant to touch the lever at rest, the stop pin at full stroke
        if frac > 0 and d_floor < 0.3:
            ok = False
        if d_clamp < 0.5 or d_disc < 0.3:
            ok = False
        if (frac < 1 and d_stop < 0.3) or (frac == 1 and d_stop > 0.05):
            ok = False
        if zmin < (1.0 if frac == 1 else 0.0):
            ok = False
    stroke_vol = math.pi * STROKE * (MEM_R ** 2 + MEM_R * DISC_R + DISC_R ** 2) / 3
    print(f"\nFull stroke displaces {stroke_vol / 1000:.0f} mL "
          f"({stroke_vol / STROKE / 1000:.2f} mL per mm of disc travel)")
    for d in ORIFICES:
        a = math.pi * d * d / 4
        print(f"  orifice Ø{d}: full-press slug L/D = {stroke_vol / a / d:4.1f}   "
              f"ring ~{1.2 * d:.0f} mm = {1.2 * d / (TUBE_OD - 2 * TUBE_WALL):.2f} of bore")
    return ok


def assembly_check(parts):
    """Every pair of parts in the assembly: report any that overlap by more than a sliver."""
    keys = list(parts)
    bad = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ba, bb = parts[a].bounding_box(), parts[b].bounding_box()
            if (ba.min.X > bb.max.X or bb.min.X > ba.max.X or ba.min.Y > bb.max.Y or
                    bb.min.Y > ba.max.Y or ba.min.Z > bb.max.Z or bb.min.Z > ba.max.Z):
                continue
            v = (parts[a] & parts[b]).volume
            if v > 0.5:
                bad.append(f"{a}/{b} {v:.1f} mm3")
    print("\nAssembly overlaps: " + (", ".join(bad) if bad else f"none among {len(keys)} parts"))
    return not bad


def stability_check(floor):
    """Support polygon of everything touching the table vs the tube axis."""
    from scipy.spatial import ConvexHull
    verts, _ = floor.tessellate(0.2, 0.3)
    pts = np.array([[p.X, p.Y] for p in verts if p.Z < 0.05])
    margin = -ConvexHull(pts).equations[:, 2].max()
    print(f"\nTipping: nearest edge of the support polygon is {margin:.1f} mm from the tube axis")
    return margin >= 80.0

# ---------------------------------------------------------------- previews


def tess(part, tol=0.2):
    import pyvista as pv
    v, t = part.tessellate(tol, 0.3)
    v = np.array([[p.X, p.Y, p.Z] for p in v])
    f = np.hstack([np.full((len(t), 1), 3), np.asarray(t)]).ravel()
    return pv.PolyData(v, f)


def previews(parts):
    import pyvista as pv
    pv.OFF_SCREEN = True
    tube = ring(TUBE_OD / 2 - TUBE_WALL, TUBE_OD / 2, Z_LT, Z_LT + TUBE_LEN)
    cap = Pos(0, WELL_Y, Z_CAP) * Cylinder(15.5, 6, align=(C, C, MIN))
    cone = Pos(0, WELL_Y, Z_CAP + 6) * (Rot(0, 0, 0) * revolve(
        Plane.XZ * Polygon((0, 0), (7.5, 0), (0, 22), align=None), Axis.Z))
    mem_pts = [(DISC_R, Z_DISC_TOP + 0.6), (MEM_R + SPIG_T / 2, Z_MEM - 0.6),
               (MEM_R + SPIG_T / 2, Z_MEM - 1.2), (DISC_R, Z_DISC_TOP)]
    membrane = rev(mem_pts) + cyl(DISC_R, Z_DISC_TOP, Z_DISC_TOP + 0.6)
    film = box(-110, CLIP_X + 3, -23, 23, Z_FT, Z_FT + 0.4)
    rings = []
    for z, d in ((Z_LT + 70, 36), (Z_LT + 190, 40), (Z_LT + 320, 44)):
        rr = d / 2
        rings.append(Pos(0, 0, z) * revolve(
            Plane.XZ * Polygon(*arc_pts(rr, 0, 4.5, 0, 360, 24)[:-1], align=None), Axis.Z))

    style = [
        ("floor", "#8a8f99", 1.0), ("wall", "#a7adb8", 1.0), ("lid", "#8a8f99", 1.0),
        ("hatch", "#c9a35a", 1.0), ("insert", "#c9a35a", 1.0), ("lever", "#3d6fb6", 1.0),
        ("rod", "#3d6fb6", 1.0), ("pin", "#3d6fb6", 1.0), ("stop", "#3d6fb6", 1.0),
        ("disc", "#3d6fb6", 1.0), ("clamp", "#c9a35a", 1.0), ("clip", "#c9a35a", 1.0),
        ("pclip", "#c9a35a", 1.0), ("sclip", "#c9a35a", 1.0),
    ]
    extras = [(membrane, "#d4537e", 1.0), (film, "#d4537e", 1.0), (cap, "#b8b8b0", 1.0),
              (cone, "#8c5a2b", 1.0)]
    extras += [(r, "#9a9a92", 0.85) for r in rings]

    def scene(section, out, cam, zoom, ortho=None):
        pl = pv.Plotter(off_screen=True, window_size=(1400, 1800 if section else 1400))
        pl.set_background("#f2f2f0")
        for key, col, op in style:
            m = tess(parts[key])
            if section:
                m = m.clip(normal=(1, 0, 0), origin=(0, 0, 0), invert=True)
            if m.n_points:                      # parts wholly on the cut-away side vanish
                pl.add_mesh(m, color=col, opacity=op, smooth_shading=False)
        for shape, col, op in extras:
            m = tess(shape)
            if section:
                m = m.clip(normal=(1, 0, 0), origin=(0, 0, 0), invert=True)
            if m.n_points:
                pl.add_mesh(m, color=col, opacity=op)
        t = tess(tube, 0.3)
        if section:
            t = t.clip(normal=(1, 0, 0), origin=(0, 0, 0), invert=True)
        pl.add_mesh(t, color="#b0c6d6", opacity=0.25)
        pl.camera_position = cam
        if ortho:
            pl.enable_parallel_projection()
            pl.camera.parallel_scale = ortho
        else:
            pl.camera.zoom(zoom)
        pl.screenshot(os.path.join(HERE, out))
        pl.close()
        print("  wrote", out)

    scene(True, "preview_section.png", [(1000, -35, 170), (0, -35, 170), (0, 0, 1)], 1.0, ortho=232)
    scene(False, "preview_base.png", [(420, -520, 380), (0, -20, 70), (0, 0, 1)], 1.15)
    scene(False, "preview_underside.png", [(330, -470, -300), (0, -40, 35), (0, 0, 1)], 1.25)

    # folio card for web/spiral-instruments.html: same cream ground, yellow parts and grey
    # tube as the OpenSCAD renders of instruments I-III
    pl = pv.Plotter(off_screen=True, window_size=(760, 1098))
    pl.set_background("#FFFFE6")
    for key in style:
        pl.add_mesh(tess(parts[key[0]]), color="#FAE14A" if key[0] != "lever" else "#E8C42A")
    for r in rings:
        pl.add_mesh(tess(r), color="#8C8C84", opacity=0.9)
    pl.add_mesh(tess(tube, 0.3), color="#DCDCD2", opacity=0.5)
    pl.camera_position = [(560, -1250, 900), (0, -30, 262), (0, 0, 1)]
    pl.camera.zoom(1.2)
    pl.screenshot(os.path.join(HERE, "preview_card.png"))
    pl.close()
    print("  wrote preview_card.png")

# ---------------------------------------------------------------- main


def main():
    print(f"Heights: pivot {Z_PIV:.2f}  membrane {Z_MEM:.2f}  floor {Z_F:.2f}-{Z_FT:.2f}  "
          f"lid {Z_LB:.2f}-{Z_LT:.2f}  stop pin {Z_STOP:.2f}  tube top {Z_LT + TUBE_LEN:.1f} mm")
    parts = {
        "floor": floor_plate(), "wall": chamber_wall(), "lid": lid(), "hatch": hatch(),
        "insert": orifice_insert(30), "lever": lever(), "rod": cross_rod(), "pin": pivot_pin(),
        "stop": stop_pin(), "disc": pusher_disc(), "clamp": membrane_clamp(), "clip": valve_clip(),
        "pclip": pivot_clip(), "sclip": stop_clip(),
    }
    ok = mechanism_check(parts)
    ok &= assembly_check(parts)
    ok &= stability_check(parts["floor"])

    print("\nExporting STLs (print orientation):")
    ok &= export(parts["floor"], "floor_plate.stl", flip=True)
    ok &= export(parts["wall"], "chamber_wall.stl")
    ok &= export(parts["lid"], "lid.stl")
    ok &= export(parts["hatch"], "hatch.stl", flip=True)
    for d in ORIFICES:
        ok &= export(orifice_insert(d), f"orifice_{d}mm.stl", flip=True)
    ok &= export(parts["lever"], "lever.stl", side=True)
    ok &= export(parts["disc"], "pusher_disc.stl", flip=True)
    ok &= export(parts["rod"], "cross_rod.stl", side=True)
    ok &= export(parts["pin"], "pivot_pin.stl", side=True)
    ok &= export(parts["stop"], "stop_pin.stl", side=True)
    ok &= export(parts["pclip"], "pin_clip.stl", side=True)
    ok &= export(parts["clamp"], "membrane_clamp.stl")
    ok &= export(parts["clip"], "valve_clip.stl", flip=True)

    if "--no-preview" not in sys.argv:
        print("\nRendering previews:")
        previews(parts)
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
