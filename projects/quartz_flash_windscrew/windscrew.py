#!/usr/bin/env python3
"""Quartz Flash Windscrew -- parametric build, PETG on a 320 mm bed.

Run from this directory:  python windscrew.py
Writes one STL per printable part (in print orientation) into stl/, plus
assembly.npz for preview.py, and prints a part table and a clearance table.

Styled as a 1926 instrument: Paris-Exposition Art Deco (ziggurat setbacks,
sunbursts, fluting, rondels, chevrons, Egyptian-revival eye) over a spirit-board
dial.  Motifs and the period typeface live in deco_motifs.py.

World frame: axle along +x (wind blows +x; the rotor nose is upwind at -x),
z up, yaw axis at x = y = 0.

  clear PETG   rotor          Liam-type conical Archimedes screw, 3 blades, one piece
               drum_shell     round tumbler drum, 12 dovetail rails; its inboard plate is
                              the planet carrier and hangs on two 608s (no through-axle)
  dark brown   tile x12       liner tiles on the rails; quartz chips epoxied on their faces
               rib_* x6       interchangeable lifters on the tiles (4 height/face profiles)
               gear_housing   stationary governor chamber + 57-tooth ring, on block_mouth
               drive_hub      pinch-clamped to the rod; the governor shoes grip it
               sun_governor   15-tooth sun on two 688s + cage disc; shoe_cage, shoe x3
               planet x3      21 teeth on 4 mm dowels in the carrier
               spacer_a/_b    inner-race spacers on the shaft
               rotor_clamp    keys into the rotor's nose hex and pinch-clamps the rod
               base           three-tier ziggurat, sunburst, MCMXXVI plaque
               post           fluted column with a stepped capital
               yaw_hub        stepped octagon, ziggurat corbels under the pad
               yoke_nose/_mouth  setback risers, knee rondels, banding; sunburst fan fin
               block_nose/_mouth octagonal 608 blocks; the nose block has stepped caps, the
                              mouth block two floating 608s, the gear housing and the pointer seat
               comm_disc      spirit-board dial (words + Major Arcana 0-XXI + four suits), the drum's lid
                              (swap to brass/gold filament at Z = DISC_T for the relief)
               pointer        stepped column, pierced fan bracket, arrow tip
               finial x2      stepped spires: over the nose collar and the dial post (brass)
               race_spacer    nose bearing's inner race to the nose collar
Hardware: plain 8 mm aluminium rod, 6061-T6 or 6082-T6 (5/16" in the US): ~480 mm axle +
~370 mm post; 5x 8 mm shaft collars (nose, drum hub, under the base, post top, yaw-hub
top); 7x 608-2RS, 2x 688-2RS, 3x 4x16 mm dowels, M3x12 x2 + nuts (pinch clamps), M5x30 x4
+ nuts (yoke feet), M5x35 x4 self-tapping (block feet), M3x10 countersunk x3 (housing),
M3x5 x3 (shoe cage), M3x10 countersunk x6 (lid), M3x20 x2 (pointer), up to 12x M5 nuts
(governor tuning mass), a 1.5 mm O-ring, cork sheet 1.5 mm, Loctite 641 for the inner races.
"""
import math, os, sys
import numpy as np
import trimesh
from trimesh.transformations import rotation_matrix
from shapely.geometry import Polygon, LineString, Point, box as sbox
from shapely.ops import unary_union

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import deco_motifs as DM
import gears as GR
import tumbler as TB
from shapely import affinity

BED = 320.0
RHO_PETG = 1.27e-3          # g/mm^3

# ------------------------------------------------------------------ hardware
ROD_CLEAR = 8.6             # M8 clearance
NUT_AF = 13.3               # M8 nut trap across flats (the rotor hub keeps these; the clamp keys into one)
ROD_D = 8.0                 # plain aluminium rod, 6061-T6 / 6082-T6 (5/16" in the US)
COLLAR_D, COLLAR_W = 16.0, 9.5   # 8 mm shaft collar with a set screw
B608_SLIP = 22.25           # mouth-block pockets: the outer races float axially
FLOAT = 1.5                 # axial float room in those pockets
NUT_H = 7.0
B608_D, B608_W = 22.15, 7.0 # press fit for PETG
B608_SHOULDER = 18.0        # bore behind the pocket: outer race only

# ------------------------------------------------------------------ rotor
R_HUB = 13.0
R_MAX = 150.0
ROTOR_L = 300.0             # nose to mouth
RULE_DEG = 50.0             # blade rulings vs the rotor's cross-plane: >= 45 -> prints unsupported nose-down
K_RULE = math.tan(math.radians(RULE_DEG))
Z_ROOT = ROTOR_L - (R_MAX - R_HUB) * K_RULE      # root height reached at the trailing edge
N_BLADE = 3
WRAP = 2 * math.pi          # one full turn per blade
W0 = 0.12                   # blade width at the nose, fraction of full
BLADE_T = 2.0
HUB_LEN = 150.0

# ------------------------------------------------------------------ frame
X_RISER = 190.0             # axle bearings at +-X_RISER
BLOCK_LEN = 16.0
BLOCK_RF = 18.0             # octagonal block, half across flats
FOOT_Z = (-30.0, 8.0)       # block foot / riser pad extent along the axle, block-local (+ = outboard)
BLOCK_BOLTS = (-24.0, -13.0)
SEAT_Z = (-16.0, 8.0)       # pointer seat on the mouth block, block-local z
SEAT_SCREWS = (-12.5, 4.5)
POINTER_DX = -4.0           # pointer sits 4 mm inboard of the block centre, clear of the gear housing
Z_AX = 200.0                # axle above the yaw-hub pad top
YOKE_T = 20.0               # yoke thickness (y)
ARM_W = 22.0
PAD_TOP = Z_AX - 24.0       # riser pad top = block foot bottom
RISER_C = X_RISER + sum(FOOT_Z) / 2              # riser centred under its pad
KNEE_R = 15.0
SETBACKS = ((26.0, PAD_TOP - 34), (30.0, PAD_TOP - 26), (34.0, PAD_TOP - 18))   # (width, z0), 8 tall
FIN_T = 4.0                 # fin is flush with the yoke's bed face (y = -YOKE_T/2)
FIN_C = (76.0, -4.0)        # fan centre, pad-top frame
FIN_R = 120.0
FIN_A = (math.radians(-90), math.radians(-12))
FIN_RAYS = 9
POST_L = 280.0
POST_D = 40.0
BASE = 200.0
BASE_T = 8.0
POST_Z0 = 12.0              # socket floor: 10 mm counterbore for the rod collar under the base
POST_GAP = 10.5             # post top to yaw hub: the post-top collar
YAW_H = 50.0
PAD_T = 10.0               # deep enough for the collar above the top yaw bearing
SKIRT_H = POST_GAP - 1.0   # yaw-hub skirt hiding the post-top collar (stops 1 mm above the post)
FOOT_BOLTS = (30.0, 47.0)   # yoke-foot bolt x positions, clear of the yaw-hub barrel
Z_PAD = POST_Z0 + POST_L + POST_GAP + YAW_H + PAD_T  # world z of the pad top
RELIEF = 1.0                # raised ornament height
ENGRAVE = 0.8               # the same ornament cut into bed-side faces

# ------------------------------------------------------------------ drum tumbler (physics in tumbler.py)
R_IN = TB.R_IN * 1000       # shell inner radius
DRUM_L = TB.DRUM_L * 1000   # inner length
R_W = TB.R_W * 1000         # quartz liner surface
SHELL_T = 2.4
N_RAILS = TB.N_RAILS
RAIL = (2.5, 4.0, 2.5)      # dovetail: half-width at its root, at its crown, height
TILE_FACE = R_IN - 4.2      # tile face; ~3 mm of epoxied quartz chips brings it to R_W
FIT = 0.25
RIB_BASE = 6.0              # rib leading-edge offset on the tile
PLATE_T = 7.0               # drum inboard plate = planet carrier, holds 608 #1
LUG_R = 66.0                # lid screws, on the dial's sector dividers

# ------------------------------------------------------------------ gearbox + governor
GEAR_M = 1.5
ZS, ZP, ZR = 15, 21, 57     # sun, planet, fixed ring: the carrier (drum) turns at 1/4.8 of the shaft
XS, XP = 0.15, -0.15        # profile shift: no undercut on the 15-tooth sun
FACE = 10.0
N_PLANETS = 3
A_GEAR = GEAR_M * (ZS + ZP) / 2
RATIO = 1 + ZR / ZS
B688_D, B688_W = 16.1, 5.0
PIN_D = 4.0                 # 4 x 16 mm steel dowels pressed into the carrier
HOUSING_R = 51.0
CHAMBER_R = 40.0
HUB_R = 16.0                # governor drive hub on the shaft
SHOE_W = 20.0
SHOE_AX = 11.0
SHOE_R = (HUB_R + 1.8, 36.0)                   # 1.5 mm cork pad on the inner face
CAGE_R = 38.0
X_CAP = X_RISER + BLOCK_LEN / 2                # mouth block's outboard face (no caps: the housing sits there)
X_NOSE_CAP = X_RISER + BLOCK_LEN / 2 + 4       # nose block's stepped-cap face
X_HB = (X_CAP, X_CAP + 4.0)                    # housing back plate, flush on the block face
X_HUB = (X_CAP + 5.7, X_CAP + 24.7)            # drive hub: 6 mm pinch-clamp collar, then the friction band
X_SHOE = (X_HUB[0] + 6.5, X_HUB[0] + 6.5 + SHOE_AX)
X_CAGE = (X_CAP + 25.5, X_CAP + 28.5)          # cage disc (the shoe cage base sits inboard of it)
X_RING = (X_CAP + 29.5, X_CAP + 40.5)          # ring gear
X_SUN = (X_CAGE[1], X_RING[1] - 0.5)
X_PLANET = (X_RING[0], X_RING[0] + FACE)
X_DRUM0 = X_RING[1] + 1.0                      # carrier plate face
X_DISC = X_DRUM0 + PLATE_T + DRUM_L            # the dial is the drum's outboard lid

# ------------------------------------------------------------------ disc
DISC_R = 110.0
DISC_T = 4.0
BEZEL_H = 2.0
WORDS = ["YES", "ABOVE", "MAYBE", "NO", "BELOW", "NOT HUMAN"]   # opposite pairs face each other
SYMBOLS = {"YES": lambda: DM.sun(), "NO": lambda: DM.crescent(),
           "ABOVE": lambda: DM.chevrons(up=True), "BELOW": lambda: DM.chevrons(up=False),
           "MAYBE": lambda: DM.sunrise(), "NOT HUMAN": lambda: DM.eye()}
NUMBERS = DM.ROMAN          # tarot spinner notation: the 22 Major Arcana, 0 (The Fool) to XXI (The World)
SUIT_R = 36.0               # the four suit signs sit in medallions in the sunburst
POINTER_APEX = 105.8        # arrow tip radius: just outside the number track

# ================================================================== helpers
def clean(m):
    """Collapse sub-micron slivers left by the booleans so the mesh still reloads
    watertight after STL's float32 rounding."""
    m = m.copy()
    m.merge_vertices(digits_vertex=4)
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.fix_normals()
    return m

def U(*ms):
    return trimesh.boolean.union([m for m in ms if m is not None], engine="manifold")

def D(a, *bs):
    return trimesh.boolean.difference([a] + list(bs), engine="manifold")

def I(a, b):
    return trimesh.boolean.intersection([a, b], engine="manifold")

def cyl(r, z0, z1, n=96):
    m = trimesh.creation.cylinder(radius=r, height=z1 - z0, sections=n)
    m.apply_translation([0, 0, (z0 + z1) / 2])
    return m

def hexp(af, z0, z1, rot=0.0):
    m = cyl(af / math.sqrt(3), z0, z1, n=6)
    m.apply_transform(rotation_matrix(math.radians(30 + rot), [0, 0, 1]))
    return m

def box(x0, x1, y0, y1, z0, z1):
    m = trimesh.creation.box(extents=[x1 - x0, y1 - y0, z1 - z0])
    m.apply_translation([(x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2])
    return m

def at(m, x=0.0, y=0.0, z=0.0):
    m = m.copy()
    m.apply_translation([x, y, z])
    return m

def along_x(m):
    """Map a z-axis part to the x axis (z -> x)."""
    m = m.copy()
    m.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))
    return m

def extrude(poly, z0, z1):
    """Extrude shapely (multi)polygon(s) along z."""
    ms = []
    for p in DM.as_polys(poly):
        e = trimesh.creation.extrude_polygon(p, z1 - z0)
        e.apply_translation([0, 0, z0])
        ms.append(e)
    return ms

def extrude_xz(poly, y0, y1):
    """Extrude shapely (multi)polygon(s) drawn in world (x, z) through y0..y1."""
    out = []
    for p in DM.as_polys(poly):
        m = trimesh.creation.extrude_polygon(p, y1 - y0)
        m.apply_translation([0, 0, y0])
        m.apply_transform(rotation_matrix(math.pi / 2, [1, 0, 0]))  # (x, y, z) -> (x, -z, y)
        m.apply_transform(np.diag([1, -1, 1, 1]))                    # flip y back
        m.fix_normals()
        out.append(m)
    return out[0] if len(out) == 1 else trimesh.util.concatenate(out) if out else None

def x_cyl(r, x0, x1, y=0.0, z=0.0, n=24):
    m = cyl(r, x0, x1, n)
    m.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))      # +z -> +x
    m.apply_translation([0, y, z])
    return m

def x_hex(af, x0, x1, y=0.0, z=0.0):
    m = hexp(af, x0, x1)
    m.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))
    m.apply_translation([0, y, z])
    return m

def pinch_cutters(r_out, z0, z1, y_screw, z_screw):
    """Split clamp on a bore along z: a slit on +y and an M3 screw across it (along x),
    head counterbore on -x, nut pocket on +x.  Grips a plain rod without denting it."""
    half = math.sqrt(max(r_out**2 - y_screw**2, 1.0))
    return [box(-0.6, 0.6, ROD_D / 2 - 0.2, r_out + 1, z0 - 1, z1 + 1),
            x_cyl(1.7, -half - 1, half + 1, y_screw, z_screw, 16),
            x_cyl(3.1, -half - 1, -3.4, y_screw, z_screw, 24),
            x_hex(5.7, 3.4, half + 1, y_screw, z_screw)]

def oct_prism(r_flat, z0, z1):
    return extrude(DM.octagon(r_flat), z0, z1)[0]

def oct_frustum(r0, z0, r1, z1):
    pts = [(*v, z0) for v in DM.octagon(r0).exterior.coords[:-1]] + \
          [(*v, z1) for v in DM.octagon(r1).exterior.coords[:-1]]
    return trimesh.convex.convex_hull(np.array(pts))

def revolve(profile, n=128):
    m = trimesh.creation.revolve(np.asarray(profile, float), sections=n)
    if m.volume < 0:
        m.invert()
    return m

def to_bed(m):
    m = m.copy()
    m.apply_translation(-m.bounds[0] * [0, 0, 1])
    return m

def face_relief(pattern, face_y, up, height=RELIEF):
    """Raised ornament on a y face (up=+1 grows toward +y) as a list of solids."""
    y0, y1 = (face_y - 0.2, face_y + height) if up > 0 else (face_y - height, face_y + 0.2)
    m = extrude_xz(pattern, y0, y1)
    return [m] if m is not None else []

def face_engrave(pattern, face_y, up, depth=ENGRAVE):
    """Cutter for the same ornament sunk into a face whose outward normal is `up`."""
    y0, y1 = (face_y - depth, face_y + 0.2) if up > 0 else (face_y - 0.2, face_y + depth)
    m = extrude_xz(pattern, y0, y1)
    return [m] if m is not None else []

# ================================================================== rotor
def blade(theta0, nt=200, ns=22):
    t = np.linspace(0, 1, nt)
    s = np.linspace(0, 1, ns)
    T, S = np.meshgrid(t, s, indexing="ij")
    width = (R_MAX - R_HUB) * (W0 + (1 - W0) * T**0.7)
    r0 = R_HUB - 1.5
    r = r0 + S * (width + (R_HUB - r0))
    th = theta0 + T * WRAP
    z = T * Z_ROOT + (r - r0) * K_RULE
    P = np.stack([r * np.cos(th), r * np.sin(th), z], -1)
    dt = np.gradient(P, axis=0)
    ds = np.gradient(P, axis=1)
    n = np.cross(dt, ds)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    A, B = P + n * BLADE_T / 2, P - n * BLADE_T / 2
    idx = lambda i, j, k: k * nt * ns + i * ns + j
    V = np.concatenate([A.reshape(-1, 3), B.reshape(-1, 3)])
    F = []
    for i in range(nt - 1):
        for j in range(ns - 1):
            a, b, c, d = idx(i, j, 0), idx(i + 1, j, 0), idx(i + 1, j + 1, 0), idx(i, j + 1, 0)
            F += [(a, b, c), (a, c, d)]
            a, b, c, d = idx(i, j, 1), idx(i + 1, j, 1), idx(i + 1, j + 1, 1), idx(i, j + 1, 1)
            F += [(a, c, b), (a, d, c)]
    def rim(pairs):
        for (i0, j0), (i1, j1) in pairs:
            a, b = idx(i0, j0, 0), idx(i1, j1, 0)
            c, d = idx(i1, j1, 1), idx(i0, j0, 1)
            F.extend([(a, d, c), (a, c, b)])
    rim([((i, 0), (i + 1, 0)) for i in range(nt - 1)])                         # root
    rim([((i + 1, ns - 1), (i, ns - 1)) for i in range(nt - 1)])               # tip edge
    rim([((0, j + 1), (0, j)) for j in range(ns - 1)])                         # leading
    rim([((nt - 1, j), (nt - 1, j + 1)) for j in range(ns - 1)])               # trailing
    m = trimesh.Trimesh(V, np.array(F), process=True)
    if m.volume < 0:
        m.invert()
    return m

def rotor():
    hub = cyl(R_HUB, 0, HUB_LEN)
    blades = [blade(2 * math.pi * k / N_BLADE) for k in range(N_BLADE)]
    m = U(hub, *blades)
    m = D(m, cyl(ROD_CLEAR / 2, -1, HUB_LEN + 1, 48),
          hexp(NUT_AF, -1, NUT_H), hexp(NUT_AF, HUB_LEN - NUT_H, HUB_LEN + 1))
    return clean(m)            # already in print orientation: nose on the bed

def overhang_report(m, limit_deg=45.0):
    n = m.face_normals
    a = m.area_faces
    down = n[:, 2] < -math.cos(math.radians(limit_deg)) + 1e-9
    on_bed = m.triangles_center[:, 2] < 0.5
    bad = down & ~on_bed
    return a[bad].sum() / a.sum()

# ================================================================== drum tumbler
def ann_sector(r0, r1, a0, a1, n=48):
    a = np.linspace(a0, a1, n)
    return Polygon([(r1 * math.cos(t), r1 * math.sin(t)) for t in a] +
                   [(r0 * math.cos(t), r0 * math.sin(t)) for t in a[::-1]])

def rot2(geom, deg):
    return affinity.rotate(geom, deg, origin=(0, 0))

def dovetail(w_wall, r_wall, w_top, r_top):
    """Dovetail centred on +y: half-width w_wall at radius r_wall, w_top at r_top."""
    return Polygon([(-w_wall, r_wall), (w_wall, r_wall), (w_top, r_top), (-w_top, r_top)])

def lid_angles():
    """Dial screw angles (clockwise from 12 o'clock in the dial frame): on the dividers."""
    return [(k + 0.5) * 2 * math.pi / len(WORDS) for k in range(len(WORDS))]

def drum_shell():
    """Local z along the axle from the inboard (carrier) face.  Clear PETG, printed on
    the carrier plate.  The dial screws on as the outboard lid."""
    z_top = PLATE_T + DRUM_L
    parts = [cyl(R_IN + SHELL_T, 0, PLATE_T, 180), cyl(R_IN + 4.0, 0, 3, 180),
             D(cyl(R_IN + SHELL_T, PLATE_T - 0.1, z_top, 180), cyl(R_IN, PLATE_T - 0.2, z_top + 1, 180)),
             D(cyl(R_IN + 4.4, z_top - 5, z_top, 180), cyl(R_IN + 1, z_top - 6, z_top + 1, 180)),
             cyl(15.0, PLATE_T - 0.1, PLATE_T + 8.0, 64)]                  # holds 608 #2; the collar sits beyond
    rail = dovetail(RAIL[0], R_IN + 0.6, RAIL[1], R_IN - RAIL[2])
    parts += extrude(unary_union([rot2(rail, 360 * k / N_RAILS) for k in range(N_RAILS)]), PLATE_T - 0.1, z_top)
    for a in lid_angles():
        lx, ly = -math.cos(a), math.sin(a)                       # dial frame -> drum-local frame
        tx, ty = -ly, lx
        pts = []
        for rr in (R_IN + 1.0, LUG_R + 4.0):
            for s in (-4.5, 4.5):
                for zz in (z_top - 6, z_top):
                    pts.append((rr * lx + s * tx, rr * ly + s * ty, zz))
        for s in (-4.5, 4.5):                                    # 45-degree gusset under the lug
            pts.append(((R_IN + 1.0) * lx + s * tx, (R_IN + 1.0) * ly + s * ty, z_top - 15))
        parts.append(trimesh.convex.convex_hull(np.array(pts)))
    m = U(*parts)
    cut = [cyl(B608_D / 2, -1, 2 * B608_W, 64), cyl(9.0, 2 * B608_W - 0.01, PLATE_T + 18, 48)]
    for k in range(N_PLANETS):
        q = math.pi / 2 + 2 * math.pi * k / N_PLANETS
        cut.append(at(cyl(PIN_D / 2 - 0.05, -1, PLATE_T + 0.5, 24), A_GEAR * math.cos(q), A_GEAR * math.sin(q)))
    for a in lid_angles():
        cut.append(at(cyl(1.25, z_top - 10, z_top + 1, 16), -LUG_R * math.cos(a), LUG_R * math.sin(a)))
    return clean(D(m, *cut))

def tile_profile():
    """Liner tile centred on +y: rides a shell rail, carries a rail for a rib, and takes
    epoxied quartz chips on its face."""
    half = math.pi / N_RAILS - math.radians(0.35)
    body = ann_sector(TILE_FACE, R_IN - 0.2, math.pi / 2 - half, math.pi / 2 + half)
    socket = dovetail(RAIL[0] + FIT, R_IN + 1, RAIL[1] + FIT, R_IN - RAIL[2] - FIT)
    rail = dovetail(RAIL[0], TILE_FACE + 0.4, RAIL[1], TILE_FACE - RAIL[2])
    return unary_union([body.difference(socket), rail])

def rib_profile(rib):
    """Lifter rib centred on +y, from (u tangential, v radial-inward) mapped onto the tile
    face.  Vertical through the chip layer, then the face leans back by beta over its
    height above the quartz.  Leading face at +u; slide it on end-for-end to reverse."""
    h_m, beta_deg = TB.RIBS[rib]
    h = h_m * 1000
    r0 = TILE_FACE - 0.1
    v0 = r0 - R_W
    lean = h * math.tan(math.radians(beta_deg))
    top_lead, top_trail = RIB_BASE - lean, RIB_BASE - lean - 3.0
    trail = min(-8.5, top_trail - 1.5)
    uv = [(trail, 0), (RIB_BASE, 0), (RIB_BASE, v0), (top_lead, v0 + h), (top_trail, v0 + h), (trail, 4.0)]
    sock = [(-(RAIL[0] + FIT), -0.6), (RAIL[0] + FIT, -0.6), (RAIL[1] + FIT, RAIL[2] + FIT),
            (-(RAIL[1] + FIT), RAIL[2] + FIT)]
    def to_xy(pts):
        dense = []
        for (u0, w0), (u1, w1) in zip(pts, pts[1:] + pts[:1]):
            for t in np.linspace(0, 1, 12, endpoint=False):
                dense.append((u0 + (u1 - u0) * t, w0 + (w1 - w0) * t))
        return Polygon([((r0 - w) * math.sin(u / r0), (r0 - w) * math.cos(u / r0)) for u, w in dense])
    # simplify: collinear points on the straight faces triangulate into slivers that clean() would drop
    return to_xy(uv).buffer(0).difference(to_xy(sock).buffer(0)).simplify(0.01)

def tile():
    return clean(extrude(tile_profile(), 0, DRUM_L - 1)[0])

def rib(name):
    return clean(extrude(rib_profile(name), 0, DRUM_L - 1)[0])

def chip_layer(k, with_rib):
    """Render-only: the epoxied quartz chips on tile k (not a printed part)."""
    half = math.pi / N_RAILS - math.radians(0.35)
    face = ann_sector(R_W, TILE_FACE + 0.05, math.pi / 2 - half, math.pi / 2 + half)
    keep_out = sbox(-(RAIL[1] + 0.8), 40, RAIL[1] + 0.8, 70)
    if with_rib:
        keep_out = rib_profile(TB.DEFAULT_RIB).buffer(0.6).union(keep_out)
    return rot2(face.difference(keep_out), 360 * k / N_RAILS - 90)

# ================================================================== gearbox + governor
def gear_shapes():
    if not hasattr(gear_shapes, "cache"):
        sun = GR.external(ZS, GEAR_M, XS)
        planet = GR.external(ZP, GEAR_M, XP)
        ring = GR.internal(ZR, ZP, GEAR_M, planet, r_outer=ZR * GEAR_M / 2 + 7, x_p=XP)
        gear_shapes.cache = (sun, planet, ring)
    return gear_shapes.cache

def housing_screw_xy():
    return [(15 * math.cos(math.pi / 2 + 2 * math.pi * k / 3), 15 * math.sin(math.pi / 2 + 2 * math.pi * k / 3))
            for k in range(3)]

def gear_housing():
    """Stationary: governor chamber + fixed ring gear.  Local z from the back plate's
    block-side face; prints on the back plate.  Three countersunk M3 into the block cap."""
    _, _, ring = gear_shapes()
    z_r0, z_r1 = X_RING[0] - X_HB[0], X_RING[1] - X_HB[0]
    parts = [cyl(HOUSING_R + 2, 0, 2, 180), cyl(HOUSING_R, 1.9, 4.1, 180),
             D(cyl(HOUSING_R, 4, z_r0 + 0.1, 180), cyl(CHAMBER_R, 3.9, z_r0 + 0.2, 180)),
             D(cyl(HOUSING_R, z_r0, z_r1, 180), cyl(ZR * GEAR_M / 2 + 6.5, z_r0 - 1, z_r1 + 1, 180)),
             D(cyl(HOUSING_R + 1.5, z_r1 - 3, z_r1, 180), cyl(HOUSING_R - 0.5, z_r1 - 4, z_r1 + 1, 180))]
    parts += extrude(ring, z_r0, z_r1)
    m = U(*parts)
    cut = [cyl(10.0, -1, 5, 64)]
    for x, y in housing_screw_xy():
        cut.append(at(cyl(1.7, -1, 5, 16), x, y))
        cut.append(at(revolve([(0, 2.3), (1.7, 2.3), (3.4, 4.02), (0, 4.02)], 24), x, y))
    for k in range(24):                                          # fluting round the chamber
        a = 2 * math.pi * k / 24
        c = trimesh.creation.capsule(height=12, radius=1.0, count=[10, 10])
        c.apply_translation([(HOUSING_R + 0.45) * math.cos(a), (HOUSING_R + 0.45) * math.sin(a),
                             6.5 - c.bounds[0, 2]])
        cut.append(c)
    return clean(D(m, *cut))

def drive_hub():
    """The governor's friction drum, pinch-clamped to the plain rod at its inboard collar
    (the shoes ride on the outboard 12 mm, clear of the screw recesses)."""
    L = X_HUB[1] - X_HUB[0]
    return clean(D(cyl(HUB_R, 0, L, 96), cyl(ROD_D / 2 + 0.05, -1, L + 1, 32),
                   *pinch_cutters(HUB_R, 0, L, 10.0, 3.0)))

def rotor_clamp():
    """Fixes the rotor to the plain rod at its nose: a hex boss keys into the rotor hub's
    nose trap, the flange pinch-clamps the rod, and a small boss bears on the nose
    bearing's inner race.  Local z from that boss face toward the rotor."""
    m = U(cyl(5.5, 0, 1.1, 48), cyl(R_HUB, 1, 9, 96), hexp(NUT_AF - 0.15, 8.9, 15))
    return clean(D(m, cyl(ROD_D / 2 + 0.05, -1, 16, 32), *pinch_cutters(R_HUB, 0, 15, 8.0, 5.0)))

SHOE_ANGLES = [math.pi / 2 + 2 * math.pi * k / 3 for k in range(3)]
SHOE_CAGE_SCREWS = [a + math.pi / 3 for a in SHOE_ANGLES]

def sun_governor():
    """Sun gear on two 688s plus the cage disc the shoe cage screws to.  Local z from
    the cage disc's inboard face; prints on that face."""
    sun, _, _ = gear_shapes()
    z_s0, z_s1 = X_SUN[0] - X_CAGE[0], X_SUN[1] - X_CAGE[0]
    z_688 = X_RING[0] + 0.5 - X_CAGE[0]
    m = U(cyl(CAGE_R, 0, z_s0 + 0.1, 128), *extrude(sun, z_s0, z_s1))
    cut = [cyl(6.1, -1, z_688, 48), cyl(B688_D / 2, z_688, z_s1 + 1, 48)]
    for a in SHOE_CAGE_SCREWS:
        cut.append(at(cyl(1.25, -1, 2.6, 16), 28 * math.cos(a), 28 * math.sin(a)))
    return clean(D(m, *cut))

def shoe_cage():
    """Guide walls for the three shoes; local z as sun_governor (base at z -2..0)."""
    z0 = X_SHOE[0] - 0.2 - X_CAGE[0]
    base = D(cyl(CAGE_R, -2, 0, 128), cyl(HUB_R + 2.5, -3, 1, 64))
    walls = []
    for a in SHOE_ANGLES:
        for s in (-1, 1):
            c = s * (SHOE_W / 2 + 0.4 + 1.5)
            w = box(c - 1.5, c + 1.5, 18, 34, z0, -1.9)
            w.apply_transform(rotation_matrix(a - math.pi / 2, [0, 0, 1]))
            walls.append(w)
    m = U(base, *walls)
    z_or = (X_SHOE[0] + X_SHOE[1]) / 2 - X_CAGE[0]
    cut = [D(cyl(CAGE_R + 2, z_or - 1.8, z_or + 1.8, 96), cyl(31.0, z_or - 3, z_or + 3, 96))]  # O-ring path
    for a in SHOE_CAGE_SCREWS:
        cut.append(at(cyl(1.7, -3, 1, 16), 28 * math.cos(a), 28 * math.sin(a)))
    return clean(D(m, *cut))

def shoe():
    """Flyweight shoe centred on +y, local z axial (0..SHOE_AX).  Cork pad on the inner
    face, O-ring in the outer groove, two pockets for up to 2 x M5 nuts each (tuning
    mass), opening toward the cage base so the nuts cannot fall out."""
    w = SHOE_W / 2
    ri, ro = SHOE_R
    xs = np.linspace(-w, w, 24)
    pts = [(x, math.sqrt(ro * ro - x * x)) for x in xs[::-1]] + [(x, math.sqrt(ri * ri - x * x)) for x in xs]
    m = extrude(Polygon(pts), 0, SHOE_AX)[0]
    groove = trimesh.creation.torus(major_radius=ro - 0.4, minor_radius=1.0, major_sections=96, minor_sections=16)
    groove.apply_translation([0, 0, SHOE_AX / 2])
    cut = [groove] + [at(hexp(8.2, SHOE_AX - 8.3, SHOE_AX + 1), x, 27.0) for x in (-5.0, 5.0)]
    return clean(D(m, *cut))

def planet():
    _, p, _ = gear_shapes()
    return clean(D(extrude(p, 0, FACE)[0], cyl(PIN_D / 2 + 0.1, -1, FACE + 1, 24)))

def spacer(length):
    return clean(D(cyl(5.5, 0, length, 48), cyl(4.3, -1, length + 1, 32)))

# ================================================================== bearing blocks
def _x_hole(r, x0, x1, z, n=24):
    """Hole along block-local x (world vertical) at block-local z."""
    h = cyl(r, x0, x1, n)
    h.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))      # +z -> +x
    h.apply_translation([0, 0, z])
    return h

def bearing_block(mouth):
    """Local frame: axle along +z (the block's outboard direction), local +x = world
    up.  The nose block has stepped caps on its outboard face (the print's top).  The
    mouth block carries two 608s, one at each end, so the overhung gearbox and drum do
    not sag the sun out of mesh; its outer races float axially for the rod's and the
    frame's different thermal expansion.  The gear housing screws to its outboard face."""
    L = BLOCK_LEN
    z0 = FOOT_Z[0]             # the body runs the full foot length so the part stands flat on z0
    parts = [oct_prism(BLOCK_RF, z0, L / 2),
             box(-24, -12, -YOKE_T / 2, YOKE_T / 2, z0, FOOT_Z[1]),
             box(-24, 0, -YOKE_T / 2, YOKE_T / 2, z0, L / 2)]
    if mouth:
        parts += [box(14, 21, -6, 6, SEAT_Z[0], SEAT_Z[1])]            # pointer seat
    else:
        parts += [cyl(15.0, L / 2 - 0.5, L / 2 + 2), cyl(12.5, L / 2 + 1.9, L / 2 + 4)]
    m = U(*parts)
    if mouth:
        cut = [cyl(B608_SLIP / 2, z0 - 1, z0 + B608_W + FLOAT),            # inboard 608, floating
               cyl(B608_SLIP / 2, L / 2 - B608_W - FLOAT, L / 2 + 1),      # outboard 608, floating
               cyl(B608_SHOULDER / 2, z0, L / 2)]
    else:
        cut = [cyl(B608_D / 2, z0 - 1, z0 + B608_W),                    # pocket faces the rotor
               cyl(B608_SHOULDER / 2, z0, L / 2 + 5)]
    cut += [_x_hole(2.7, -30, -6, zf) for zf in BLOCK_BOLTS]
    if mouth:
        cut += [_x_hole(1.4, 12, 22, zz, 16) for zz in SEAT_SCREWS]     # pointer screws, self-tap M3
        for x, y in housing_screw_xy():                                  # gear housing, self-tap M3
            cut.append(at(cyl(1.25, L / 2 - 9, L / 2 + 1, 16), -x, -y))
    return clean(D(m, *cut))   # local frame; print_orient stands it on its rotor-side face

def block_to_world(m_local, mouth):
    m = m_local.copy()
    R = np.eye(4)
    if mouth:
        R[:3, :3] = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
        m.apply_transform(R)
        m.apply_translation([X_RISER, 0, Z_PAD + Z_AX])
    else:
        R[:3, :3] = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        m.apply_transform(R)
        m.apply_translation([-X_RISER, 0, Z_PAD + Z_AX])
    return m

def race_spacer():
    """Tube from the nose bearing's inner race to the outboard collar, inside the block bore."""
    length = (BLOCK_LEN / 2 + 4) - (FOOT_Z[0] + B608_W) + 1
    return clean(D(cyl(6.0, 0, length, 48), cyl(4.3, -1, length + 1, 32)))

# ================================================================== yoke
def yoke_profile():
    foot = sbox(0, 0, 55, 14)
    arm = LineString([(40, 7), (RISER_C, 22)]).buffer(ARM_W / 2, cap_style=2)
    riser = LineString([(RISER_C, 15), (RISER_C, SETBACKS[0][1])]).buffer(ARM_W / 2, cap_style=2)
    knee = Point(RISER_C, 20).buffer(KNEE_R, 64)
    steps = [sbox(RISER_C - w / 2, z0, RISER_C + w / 2, z0 + 8) for w, z0 in SETBACKS]
    pad = sbox(X_RISER + FOOT_Z[0], PAD_TOP - 10, X_RISER + FOOT_Z[1], PAD_TOP)
    return unary_union([foot, arm, riser, knee, pad] + steps).intersection(sbox(0, -500, 500, 500))

def yoke_pattern():
    """Triple banding on the arm and riser, a rondel at the knee, a chevron under
    the setbacks, lozenges in the steps and on the foot."""
    a0, a1 = np.array([40.0, 7.0]), np.array([RISER_C, 22.0])
    L = np.linalg.norm(a1 - a0)
    d = (a1 - a0) / L
    n = np.array([-d[1], d[0]])
    f = [DM.stroke([a0 + d * 18 + n * o, a0 + d * (L - KNEE_R - 3) + n * o]) for o in (-5, 0, 5)]
    c = (RISER_C, 20.0)
    f += [DM.ring(11.2, 12.4, c), DM.ring(7.4, 8.6, c), Point(c).buffer(4.0, 48)]
    top = SETBACKS[0][1]
    f += [DM.stroke([(RISER_C + o, 20 + KNEE_R + 3), (RISER_C + o, top - 10)]) for o in (-5, 0, 5)]
    f.append(DM.chevron((RISER_C, top - 5), 14, 4, 1.2, up=True))
    f += [DM.lozenge((RISER_C, z0 + 4), 9, 4.6) for _, z0 in SETBACKS]
    f.append(DM.lozenge((14, 7), 12, 7))
    return unary_union(f).intersection(yoke_profile().buffer(-1.6))

def fin_profile():
    a = np.linspace(*FIN_A, FIN_RAYS + 1)
    body = DM.sector(FIN_C, FIN_R - 6, *FIN_A)
    bumps = [Point(FIN_C[0] + (FIN_R - 6) * math.cos(t), FIN_C[1] + (FIN_R - 6) * math.sin(t)).buffer(7.5, 48)
             for t in a]
    arm_z = lambda x: 7 + 15 * (x - 40) / (RISER_C - 40)            # arm centreline: the fin runs up into it
    x0, x1 = FIN_C[0] - 6, RISER_C + 8
    top = Polygon([(x0, FIN_C[1] - 3), (x1, FIN_C[1] - 3), (x1, arm_z(x1)), (x0, arm_z(x0))])  # overlaps the fan apex
    return unary_union([body, top] + bumps).intersection(sbox(FIN_C[0] - 8, -400, RISER_C + 10, 40))

def fin_pattern():
    """A rising-sun fan: tapered rays, an outer arc, a quarter sun at the hinge."""
    step = (FIN_A[1] - FIN_A[0]) / FIN_RAYS
    f = [DM.ray(FIN_C, FIN_A[0] + (i + 0.5) * step, 22, FIN_R - 22, 2.8, 1.2) for i in range(FIN_RAYS)]
    wedge = DM.sector(FIN_C, FIN_R + 10, *FIN_A)
    f.append(DM.ring(FIN_R - 17, FIN_R - 15.8, FIN_C).intersection(wedge))
    f.append(DM.sector(FIN_C, 9.0, *FIN_A))
    f.append(DM.ring(13.0, 14.2, FIN_C).intersection(wedge))
    face = fin_profile().difference(yoke_profile()).buffer(-1.5)
    return unary_union(f).intersection(face)

def yoke(mouth):
    prof = yoke_profile()
    pat = yoke_pattern()
    m = extrude_xz(prof, -YOKE_T / 2, YOKE_T / 2)
    adds = face_relief(pat, YOKE_T / 2, +1)
    cuts = face_engrave(pat, -YOKE_T / 2, -1)
    if mouth:
        adds.append(extrude_xz(fin_profile(), -YOKE_T / 2, -YOKE_T / 2 + FIN_T))
        fp = fin_pattern()
        adds += face_relief(fp, -YOKE_T / 2 + FIN_T, +1)
        cuts += face_engrave(fp, -YOKE_T / 2, -1)
    m = U(m, *adds)
    for x in FOOT_BOLTS:                                               # feet to yaw-hub pad
        cuts.append(at(cyl(2.7, -1, 15, 24), x))
    for dx in BLOCK_BOLTS:                                             # block feet: M5 self-tapping pilots
        cuts.append(at(cyl(2.1, PAD_TOP - 25, PAD_TOP + 1, 24), X_RISER + dx))
    m = D(m, *cuts)
    if not mouth:
        m.apply_transform(np.diag([-1, 1, 1, 1]))
        m.fix_normals()
    return clean(m)            # world-local (z from pad top); print lying on its -y face

# ================================================================== yaw hub, post, base
def yaw_hub():
    body = U(oct_prism(25, 0, 6), oct_frustum(25, 6, 22, 9), oct_prism(22, 9, 37),
             oct_frustum(22, 37, 25, 40), oct_prism(25, 40, YAW_H))
    pad = extrude(DM.clipped_rect(-55, -YOKE_T / 2 - 1, 55, YOKE_T / 2 + 1, 5), YAW_H, YAW_H + PAD_T)[0]
    # ziggurat corbels carry the pad ends (4 x 4 mm steps, printable body-down)
    corbels = []
    for s in (-1, 1):
        pts, x, z = [(20, YAW_H), (55, YAW_H)], 55.0, YAW_H
        while x > 24:
            z -= 4
            pts += [(x, z), (x - 4, z)]
            x -= 4
        pts.append((20, z))
        poly = Polygon([(s * px, pz) for px, pz in pts]).buffer(0)
        corbels.append(extrude_xz(poly, -YOKE_T / 2 - 1, YOKE_T / 2 + 1))
    # skirt over the post-top nut, 45-degree inner cone so it prints body-down
    skirt = D(oct_prism(25, -SKIRT_H, 0.5), revolve([(0, -SKIRT_H - 0.1), (23, -SKIRT_H - 0.1), (16, 0.01), (0, 0.01)], 96))
    m = U(body, pad, skirt, *corbels)
    cut = [cyl(B608_D / 2, -SKIRT_H - 1, B608_W), cyl(B608_D / 2, YAW_H - B608_W, YAW_H + 0.01),
           cyl(B608_SHOULDER / 2, 0, YAW_H), cyl(9, YAW_H - 1, YAW_H + PAD_T + 1)]
    for x in [-b for b in FOOT_BOLTS] + list(FOOT_BOLTS):
        cut.append(at(cyl(2.7, 0, YAW_H + PAD_T + 1, 24), x))
        cut.append(at(hexp(8.4, 0, YAW_H), x))                         # nut channel, open below
    return clean(D(m, *cut))   # world orientation; print body-down

def post():
    prof = [(0, 0), (POST_D / 2, 0), (POST_D / 2, 258), (23, 261), (23, 265), (26, 268),
            (26, 276), (22, 276), (22, POST_L), (0, POST_L)]
    m = revolve(prof)
    flutes = []
    for k in range(16):
        a = 2 * math.pi * k / 16
        c = trimesh.creation.capsule(height=196, radius=2.0, count=[16, 16])
        c.apply_translation([20.7 * math.cos(a), 20.7 * math.sin(a), 50 - c.bounds[0, 2]])   # stopped flutes from z 50
        flutes.append(c)
    m = D(m, cyl(4.5, -1, POST_L + 1, 32), *flutes)
    return clean(m)

def base():
    tiers = [(BASE / 2, 22, 0, BASE_T), (82, 18, BASE_T, 15), (62, 14, 15, 22)]
    parts = []
    for h, ch, z0, z1 in tiers:
        parts += extrude(DM.clipped_rect(-h, -h, h, h, ch), z0, z1)
    parts += [oct_prism(37, 22, 29), oct_prism(32, 29, 36), cyl(26, 36, 45)]
    rays = DM.sunburst((0, 0), 38, 60, 52, 32, 3.0, 1.2, a0=math.pi / 32)
    rays = rays.intersection(DM.clipped_rect(-62, -62, 62, 62, 14).buffer(-2.5)).difference(DM.octagon(37.5))
    parts += extrude(rays, 22 - 0.2, 22 + RELIEF)
    plaque = unary_union([DM.line_text("MCMXXVI", 6.5, (0, -72)),
                          DM.lozenge((-33, -72), 8, 4.6), DM.lozenge((33, -72), 8, 4.6),
                          sbox(-44, -66.0, 44, -65.0), sbox(-44, -79.0, 44, -78.0)])
    parts += extrude(plaque, 15 - 0.2, 15 + 0.8)
    m = U(*parts)
    cut = [cyl(POST_D / 2 + 0.25, POST_Z0, 46), cyl(ROD_D / 2 + 0.3, -1, POST_Z0 + 1, 32),
           cyl(COLLAR_D / 2 + 0.4, -1, COLLAR_W + 0.5, 64)]           # rod collar, recessed under the base
    for sx in (-1, 1):
        for sy in (-1, 1):
            cut.append(at(cyl(3.2, -1, BASE_T + 1, 24), sx * 78, sy * 78))
    return clean(D(m, *cut))

# ================================================================== communication disc
def comm_disc():
    """Spirit-board dial, also the drum's outboard lid.  Print frame: X = viewer's right,
    Y = up, relief on +Z."""
    base_m = cyl(DISC_R, 0, DISC_T, 180)
    hub = U(cyl(17, DISC_T - 0.5, DISC_T + 3, 64), cyl(14, DISC_T + 2.9, DISC_T + 6, 64),
            cyl(11.5, DISC_T + 5.9, DISC_T + 8, 64))
    n = len(WORDS)
    sector = 2 * math.pi / n
    f = [DM.ring(106.0, 107.2), DM.ring(104.5, 105.4)]                   # bezel step, number track
    nn = len(NUMBERS)
    for k in range(nn):                                                  # a tick between neighbouring trumps
        f.append(DM.place(sbox(-0.45, 0, 0.45, 3.6), 100.6, (k + 0.5) * 2 * math.pi / nn))
    for k, s in enumerate(NUMBERS):
        f += DM.arc_text(s, 96.4, k * 2 * math.pi / nn, 6.5, tracking=0.14)
    f += [DM.ring(91.0, 92.0), DM.ring(88.6, 89.5)]
    for k in range(n):                                                   # sector dividers
        a = (k + 0.5) * sector
        f.append(DM.place(sbox(-0.5, 57.0, 0.5, 88.0), 0.0, a))
        f.append(DM.place(DM.lozenge((0, 0), 3.8, 6.4), 74.0, a))
    for k, w in enumerate(WORDS):
        f += DM.arc_text(w, 71.0, k * sector, 8.2, max_arc=0.78 * 71.0 * sector)
        f.append(DM.place(SYMBOLS[w](), 62.5, k * sector))
    f += [DM.ring(54.6, 55.6), DM.ring(52.4, 53.3)]
    suit_at = [math.pi / 4 + j * math.pi / 2 for j in range(4)]         # wands, cups, swords, pentacles
    medallions = unary_union([Point(SUIT_R * math.sin(a), SUIT_R * math.cos(a)).buffer(9.4, 64) for a in suit_at])
    f.append(DM.sunburst((0, 0), 20.2, 50.5, 43.0, 36, 2.2, 0.95).difference(medallions))
    for a, glyph in zip(suit_at, DM.SUITS):
        f.append(DM.place(glyph(), SUIT_R, a))
        f.append(DM.place(DM.ring(7.9, 8.8), SUIT_R, a))
    f.append(DM.ring(18.4, 19.4))
    relief = unary_union(f).intersection(Point(0, 0).buffer(105.95, 180)).difference(Point(0, 0).buffer(17.4, 64))
    bezel = DM.ring(107.2 - 0.05, DISC_R, res=180)
    m = U(base_m, hub, *extrude(relief, DISC_T - 0.2, DISC_T + RELIEF),
          *extrude(bezel, DISC_T - 0.2, DISC_T + BEZEL_H))
    m = U(m, cyl(COLLAR_D / 2, DISC_T + 7.9, DISC_T + 8 + COLLAR_W + 0.5, 64))   # finial post
    cut = []
    for a in lid_angles():                                               # countersunk M3 into the drum lugs
        x, y = LUG_R * math.sin(a), LUG_R * math.cos(a)
        cut.append(at(cyl(1.7, -1, DISC_T + 2, 16), x, y))
        cut.append(at(revolve([(0, DISC_T + RELIEF - 1.6), (1.7, DISC_T + RELIEF - 1.6),
                               (3.3, DISC_T + RELIEF + 0.02), (0, DISC_T + RELIEF + 0.02)], 24), x, y))
    return clean(D(m, *cut))   # print orientation: dial up; world: print z -> +x, print Y -> +z

def disc_to_world(m_print):
    m = m_print.copy()
    R = np.eye(4)
    R[:3, :3] = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])          # X->y, Y->z, Z->x
    m.apply_transform(R)
    m.apply_translation([X_DISC, 0, Z_PAD + Z_AX])
    return m

# ================================================================== pointer + finials
def pointer():
    """World-local (x, z - axle), 10 mm thick in y.  Stepped plinth and crown, a pierced
    quarter-sun bracket in the corner, and a drop that ends in a 90-degree arrow
    (shaped in y so the reader, looking along -x, sees the point)."""
    xf = X_DISC + DISC_T + BEZEL_H + 3.0
    zt = DISC_R + 14.0
    xc = X_RISER + POINTER_DX
    body = unary_union([
        sbox(xc - 12, 21, xc + 11.5, 34),               # plinth (screws at +-8.5 from its centre)
        sbox(xc - 8, 34, xc + 8, 38),                   # step
        sbox(xc - 5, 34, xc + 5, zt + 8),               # shaft
        sbox(xc - 5, zt, xf + 8, zt + 8),               # beam
        sbox(xc - 7, zt + 8, xc + 7, zt + 11), sbox(xc - 5, zt + 11, xc + 5, zt + 14),
        sbox(xc - 3, zt + 14, xc + 3, zt + 17),         # stepped crown
        sbox(xf - 2, zt + 8, xf + 10, zt + 11), sbox(xf, zt + 11, xf + 8, zt + 13),
    ])
    fc = (xc + 5, zt)
    fan = DM.sector(fc, 34, -math.pi / 2, 0)
    for a in np.radians([-18, -36, -54, -72]):
        fan = fan.difference(DM.ray(fc, a, 11, 29, 2.4, 2.4))
    fan = fan.difference(DM.ring(6.0, 7.2, fc))
    m = extrude_xz(unary_union([body, fan]), -5, 5)
    drop = extrude_xz(sbox(xf, POINTER_APEX, xf + 8, zt + 1), -5, 5)
    r = POINTER_APEX
    arrow_yz = Polygon([(-5, zt + 9), (5, zt + 9), (5, r + 7), (3.5, r + 7), (3.5, r + 3.5),
                        (0, r), (-3.5, r + 3.5), (-3.5, r + 7), (-5, r + 7)])
    ap = trimesh.creation.extrude_polygon(arrow_yz, 12)                # (y, z) extruded along +x
    R = np.eye(4)
    R[:3, :3] = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])          # local (y, z, x) -> world
    ap.apply_transform(R)
    ap.apply_translation([xf - 2, 0, 0])
    m = U(m, I(drop, ap))
    cut = [at(cyl(1.7, 20, 36, 16), X_RISER + dz) for dz in SEAT_SCREWS]
    return clean(D(m, *cut))   # world-local (z from axle); print lying on a y face

def finial():
    """Stepped spire pressed over the nose shaft collar (it clears the stationary cap by
    1 mm), and over the dial's centre post."""
    d = COLLAR_W + 1.0
    m = U(cyl(12, 0, d + 2, 64), oct_prism(10, d + 1.9, d + 6), oct_prism(8, d + 5.9, d + 10),
          oct_prism(6, d + 9.9, d + 14), revolve([(0, d + 13.9), (4.6, d + 13.9), (0.9, d + 24), (0, d + 24)], 48))
    return clean(D(m, cyl(COLLAR_D / 2 + 0.15, -1, d, 64), cyl(ROD_D / 2 + 0.4, -1, d + 2, 32)))

# ================================================================== assembly + report
def report(name, m, material, count=1):
    ext = m.extents
    fits = all(e <= BED + 1e-6 for e in ext)
    print(f"{name:14s} x{count}  {material:10s} {ext[0]:6.1f} x {ext[1]:6.1f} x {ext[2]:6.1f} mm  "
          f"watertight={str(m.is_watertight):5s} {m.volume/1000:7.1f} cm3  ~{m.volume*RHO_PETG:5.0f} g solid"
          f"{'' if fits else '  ** EXCEEDS BED **'}")
    assert fits, f"{name} exceeds the bed"

FLIPPED = ("shoe_cage",)

def print_orient(name, m):
    if name in ("yoke_nose", "yoke_mouth", "pointer"):
        m = m.copy()
        m.apply_transform(rotation_matrix(math.pi / 2, [1, 0, 0]))    # -y face on the bed
        return to_bed(m)
    if name in FLIPPED:
        m = m.copy()
        m.apply_transform(rotation_matrix(math.pi, [1, 0, 0]))
        return to_bed(m)
    return to_bed(m)

def on_axle(m, x0, deg=0.0):
    """Place a part built along local z (local z = 0 at world x0) onto the axle."""
    m = m.copy()
    if deg:
        m.apply_transform(rotation_matrix(math.radians(deg), [0, 0, 1]))
    return at(along_x(m), x0, 0, Z_PAD + Z_AX)

def assembly(P):
    zax = Z_PAD + Z_AX
    W = [("base", P["base"], "brown"),
         ("post", at(P["post"], z=POST_Z0), "brown"),
         ("yaw_hub", at(P["yaw_hub"], z=Z_PAD - YAW_H - PAD_T), "brown")]
    for n in ("yoke_nose", "yoke_mouth"):
        W.append((n, at(P[n], z=Z_PAD), "brown"))
    W.append(("block_nose", block_to_world(P["block_nose"], False), "brown"))
    W.append(("block_mouth", block_to_world(P["block_mouth"], True), "brown"))
    r = P["rotor"].copy()
    r.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))     # print z -> world +x
    W.append(("rotor", at(r, -ROTOR_L / 2, 0, zax), "glass"))
    fc = P["rotor_clamp"].copy()
    W.append(("rotor_clamp", on_axle(fc, -X_RISER - FOOT_Z[0]), "brass"))
    # gearbox + governor
    W.append(("gear_housing", on_axle(P["gear_housing"], X_HB[0]), "brown"))
    W.append(("drive_hub", on_axle(P["drive_hub"], X_HUB[0]), "brown"))
    W.append(("sun_governor", on_axle(P["sun_governor"], X_CAGE[0]), "brass"))
    W.append(("shoe_cage", on_axle(P["shoe_cage"], X_CAGE[0]), "brown"))
    for k, a in enumerate(SHOE_ANGLES):
        W.append((f"shoe_{k}", on_axle(P["shoe"], X_SHOE[0], math.degrees(a) - 90), "brass"))
    for k in range(N_PLANETS):
        q = math.pi / 2 + 2 * math.pi * k / N_PLANETS
        pl = at(P["planet"], A_GEAR * math.cos(q), A_GEAR * math.sin(q))
        W.append((f"planet_{k}", on_axle(pl, X_PLANET[0]), "brass"))
    # drum: shell, liner tiles with quartz chips, ribs on alternate tiles, dial lid
    W.append(("drum_shell", on_axle(P["drum_shell"], X_DRUM0), "glass"))
    z_liner = X_DRUM0 + PLATE_T + 0.5
    for k in range(N_RAILS):
        deg = 360 * k / N_RAILS - 90
        W.append((f"tile_{k}", on_axle(P["tile"], z_liner, deg), "brown"))
        chips = extrude(chip_layer(k, k % 2 == 0), 0, DRUM_L - 1)
        if chips:
            W.append((f"chips_{k}", on_axle(trimesh.util.concatenate(chips), z_liner), "quartz"))
        if k % 2 == 0:
            W.append((f"rib_{k}", on_axle(P[f"rib_{TB.DEFAULT_RIB}"], z_liner, deg), "brown"))
    W.append(("comm_disc", disc_to_world(P["comm_disc"]), "brown"))
    W.append(("pointer", at(P["pointer"], z=zax), "brown"))
    fn = P["finial"].copy()
    fn.apply_transform(rotation_matrix(-math.pi / 2, [0, 1, 0]))    # local z -> world -x
    W.append(("finial_nose", at(fn, -(X_NOSE_CAP + 1), 0, zax), "brass"))
    W.append(("finial_disc", at(along_x(P["finial"]), X_DISC + DISC_T + 8, 0, zax), "brass"))
    collar = D(cyl(COLLAR_D / 2, 0, COLLAR_W, 48), cyl(ROD_D / 2, -1, COLLAR_W + 1, 24))
    nose_c = collar.copy(); nose_c.apply_transform(rotation_matrix(-math.pi / 2, [0, 1, 0]))
    W.append(("collar_nose", at(nose_c, -(X_NOSE_CAP + 1), 0, zax), "steel"))
    W.append(("collar_drum", on_axle(collar, X_DRUM0 + 2 * B608_W + 1.6), "steel"))
    rod_end = X_DRUM0 + 2 * B608_W + 1.6 + COLLAR_W
    W.append(("rod", at(along_x(cyl(ROD_D / 2, -(X_NOSE_CAP + 1 + COLLAR_W), rod_end, 24)), 0, 0, zax), "alu"))
    return W

def clearance_checks(W):
    """Rotating parts sweep surfaces of revolution about the axle, so a stationary
    part clears them iff its (x, r) image stays off theirs.  Surface samples of both
    in the (x, r) half-plane, nearest-neighbour distance = clearance (~0.5 mm res.).
    Parts turning at different speeds (rotor/shaft, governor, drum) are checked
    against each other the same way.  Meshing gears are checked in 2D by GR.mesh_check."""
    from scipy.spatial import cKDTree
    zax = Z_PAD + Z_AX
    parts = {n: m for n, m, _ in W}
    def xr(m, n):
        p, _ = trimesh.sample.sample_surface_even(m, n)
        return np.column_stack([p[:, 0], np.hypot(p[:, 1], p[:, 2] - zax)])
    shoes = trimesh.util.concatenate([parts[f"shoe_{k}"] for k in range(3)])
    drum = trimesh.util.concatenate([parts["drum_shell"], parts["comm_disc"]])
    parts["shoes"], parts["drum+lid"] = shoes, drum
    pairs = [("rotor", ["yoke_nose", "yoke_mouth", "block_nose", "block_mouth", "pointer", "gear_housing"]),
             ("drum+lid", ["gear_housing", "pointer", "block_mouth", "yoke_mouth", "sun_governor"]),
             ("drive_hub", ["gear_housing", "block_mouth", "sun_governor", "shoe_cage"]),
             ("sun_governor", ["gear_housing"]),
             ("shoe_cage", ["gear_housing"]),
             ("shoes", ["gear_housing", "sun_governor"]),
             ("finial_nose", ["block_nose"]),
             ("rotor_clamp", ["block_nose"]),
             ("finial_disc", ["pointer"])]
    worst = {}
    for rot, statics in pairs:
        tree = cKDTree(xr(parts[rot], 80000))
        for st in statics:
            d, _ = tree.query(xr(parts[st], 30000))
            worst[(rot, st)] = d.min()
            print(f"  clearance {rot:13s} <-> {st:13s} {d.min():6.1f} mm")
    return worst

if __name__ == "__main__":
    os.makedirs("stl", exist_ok=True)
    for stale in ("ico_inboard", "ico_outboard", "paddle"):
        if os.path.exists(os.path.join("stl", f"{stale}.stl")):
            os.remove(os.path.join("stl", f"{stale}.stl"))
    P = {}
    P["rotor"] = rotor()
    P["block_nose"] = bearing_block(False)
    P["block_mouth"] = bearing_block(True)
    P["race_spacer"] = race_spacer()
    P["rotor_clamp"] = rotor_clamp()
    P["gear_housing"] = gear_housing()
    P["drive_hub"] = drive_hub()
    P["sun_governor"] = sun_governor()
    P["shoe_cage"] = shoe_cage()
    P["shoe"] = shoe()
    P["planet"] = planet()
    P["spacer_a"] = spacer(X_RING[0] + 0.5 - X_HUB[1])       # drive hub -> inboard 688 inner race
    P["spacer_b"] = spacer(X_DRUM0 - (X_SUN[1]))              # outboard 688 -> drum 608 inner race
    P["drum_shell"] = drum_shell()
    P["tile"] = tile()
    for name in TB.RIBS:
        P[f"rib_{name}"] = rib(name)
    P["yoke_nose"] = yoke(False)
    P["yoke_mouth"] = yoke(True)
    P["yaw_hub"] = yaw_hub()
    P["post"] = post()
    P["base"] = base()
    P["comm_disc"] = comm_disc()
    P["pointer"] = pointer()
    P["finial"] = finial()

    mats = {"rotor": "clear", "drum_shell": "clear", "comm_disc": "brown+brass", "finial": "brass",
            "sun_governor": "brown", "shoe": "brown"}
    counts = {"finial": 2, "tile": N_RAILS, "shoe": 3, "planet": N_PLANETS}
    counts.update({f"rib_{n}": TB.N_RIBS for n in TB.RIBS})
    print(f"axle height {Z_PAD + Z_AX:.0f} mm, rotor top {Z_PAD + Z_AX + R_MAX:.0f} mm, "
          f"gearbox {X_HB[0]:.1f}-{X_RING[1]:.1f}, drum {X_DRUM0:.1f}-{X_DISC:.1f}, ratio {RATIO:.2f}:1")
    for name, m in P.items():
        pm = print_orient(name, m)
        path = os.path.join("stl", f"{name}.stl")
        pm.export(path)
        back = trimesh.load(path)
        report(name, back, mats.get(name, "dark brown"), counts.get(name, 1))
        assert back.is_watertight, f"{name} not watertight after reload"
    print(f"rotor unsupported overhang (>45 deg from vertical), area fraction: "
          f"{overhang_report(P['rotor']):.3f}")
    sun, pl, ring = gear_shapes()
    ov = GR.mesh_check(sun, pl, ring, ZS, ZP, ZR, GEAR_M)
    print(f"gear mesh: worst overlap over a mesh pitch {ov:.3f} mm^2")
    assert ov < 0.01, "gears interfere"
    W = assembly(P)
    worst = clearance_checks(W)
    assert all(v > 0.5 for v in worst.values()), "a stationary part touches a rotating one"
    np.savez_compressed("assembly.npz", **{f"{i:02d}_{n}_{mat}_v": m.vertices for i, (n, m, mat) in enumerate(W)},
                        **{f"{i:02d}_{n}_{mat}_f": m.faces for i, (n, m, mat) in enumerate(W)})
    print("wrote stl/*.stl and assembly.npz")
