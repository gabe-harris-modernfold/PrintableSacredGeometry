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
               ico_inboard    icosahedron half facing the riser (sleeve hole), beaded edges
               ico_outboard   icosahedron half facing the disc (nut-trap hub), beaded edges
  dark brown   base           three-tier ziggurat, sunburst, MCMXXVI plaque
               post           fluted column with a stepped capital
               yaw_hub        stepped octagon, ziggurat corbels under the pad
               yoke_nose/_mouth  setback risers, knee rondels, banding; sunburst fan fin
               block_nose/_mouth octagonal 608 blocks with stepped caps; mouth carries the
                              fixed sleeve (608 in its tip) and the pointer seat
               paddle         stationary striker inside the icosahedron
               comm_disc      spirit-board dial: words, symbols, number ring, sunburst
                              (swap to brass/gold filament at Z = DISC_T for the relief)
               pointer        stepped column, pierced fan bracket, arrow tip
               finial x2      stepped spires over the two axle-end nuts (brass filament)
               race_spacer    pushes the nose bearing's inner race
Hardware: M8 rod 550 mm (axle) + 360 mm (yaw/post), 5x 608-2RS, M8 nuts/washers,
M5x30 x4 + nuts (yoke feet), M5x35 x4 self-tapping (block feet into the risers),
M3x16 x6 (icosahedron flange), M3x20 x2 (pointer), M3x10 x1 (paddle).
"""
import math, os, sys
import numpy as np
import trimesh
from trimesh.transformations import rotation_matrix
from shapely.geometry import Polygon, LineString, Point, box as sbox
from shapely.ops import unary_union

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import deco_motifs as DM

BED = 320.0
RHO_PETG = 1.27e-3          # g/mm^3

# ------------------------------------------------------------------ hardware
ROD_CLEAR = 8.6             # M8 clearance
NUT_AF = 13.3               # M8 nut trap across flats
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
SEAT_Z = (-12.0, 12.0)      # pointer seat on the mouth block, block-local z
SEAT_SCREWS = (-8.5, 8.5)
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
YAW_H = 50.0
PAD_T = 8.0
SKIRT_H = 7.0              # yaw-hub skirt hiding the post-top nut (stops 1 mm above the post)
FOOT_BOLTS = (30.0, 47.0)   # yoke-foot bolt x positions, clear of the yaw-hub barrel
Z_PAD = BASE_T + POST_L + 8.0 + YAW_H + PAD_T  # world z of the pad top (post-top nut + washer = 8)
RELIEF = 1.0                # raised ornament height
ENGRAVE = 0.8               # the same ornament cut into bed-side faces

# ------------------------------------------------------------------ icosahedron tumbler
ICO_RC = 70.0               # circumradius
ICO_WALL = 2.4
ICO_EDGE = ICO_RC / math.sin(2 * math.pi / 5)
ICO_RI = ICO_EDGE * math.sqrt(3) * (3 + math.sqrt(5)) / 12      # inradius
RAIL_R = 1.8                # beaded edges -- a cut crystal in a leaded frame
STUD_R = 3.2
SLEEVE_D = 30.0
SLEEVE_HEX = 28.0           # tip section, across flats, carries the paddle
SLEEVE_GAP = 12.0           # block face to icosahedron face
X_ICO = X_RISER + BLOCK_LEN / 2 + SLEEVE_GAP + ICO_RI
SLEEVE_HOLE = 34.0
TAB_N = 6
STONE_GAP = 6.0             # paddle tip stays this far inside the walls
PADDLE_T = 10.0
PADDLE_HALF = 30.0
PADDLE_STONE_D = 16.5       # through-hole: one quartz stone epoxied through, both faces

# ------------------------------------------------------------------ disc
DISC_R = 110.0
DISC_T = 4.0
BEZEL_H = 2.0
WORDS = ["YES", "ABOVE", "MAYBE", "NO", "BELOW", "NOT HUMAN"]   # opposite pairs face each other
SYMBOLS = {"YES": lambda: DM.sun(), "NO": lambda: DM.crescent(),
           "ABOVE": lambda: DM.chevrons(up=True), "BELOW": lambda: DM.chevrons(up=False),
           "MAYBE": lambda: DM.sunrise(), "NOT HUMAN": lambda: DM.eye()}
NUMBERS = [str(i) for i in range(10)]
X_DISC = X_ICO + ICO_RI + 0.6                                     # disc back face
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

# ================================================================== icosahedron
def ico_solid(scale):
    ico = trimesh.creation.icosahedron()
    ico.apply_scale(ICO_RC / np.linalg.norm(ico.vertices[0]))
    n = ico.face_normals[0]
    ico.apply_transform(trimesh.geometry.align_vectors(n, [0, 0, -1]))
    # one bottom-face vertex onto +x so the tab ring is mirror-symmetric
    bot = ico.vertices[ico.faces[0]]
    ang = math.atan2(bot[0, 1], bot[0, 0])
    ico.apply_transform(rotation_matrix(-ang, [0, 0, 1]))
    ico.apply_scale(scale)
    return ico

def ico_inner_min_radius(zl):
    """Inscribed radius of the inner cross-section at local height zl (exact for a
    convex polytope: min over face planes of the distance to the slice line)."""
    inner = ico_solid((ICO_RI - ICO_WALL) / ICO_RI)
    d = np.einsum("ij,ij->i", inner.face_normals, inner.triangles[:, 0])
    nxy = np.linalg.norm(inner.face_normals[:, :2], axis=1)
    ok = nxy > 1e-6
    return np.min((d[ok] - inner.face_normals[ok, 2] * zl) / nxy[ok])

def ico_halves():
    outer = ico_solid(1.0)
    inner = ico_solid((ICO_RI - ICO_WALL) / ICO_RI)
    big = 2 * ICO_RC
    # beaded edges and vertex studs, clipped to the two face planes so each half
    # still sits flat on its outer face
    rails = [trimesh.creation.cylinder(radius=RAIL_R, segment=outer.vertices[e], sections=16)
             for e in outer.edges_unique]
    studs = []
    for v in outer.vertices:
        s = trimesh.creation.icosphere(subdivisions=2, radius=STUD_R)
        s.apply_translation(v)
        studs.append(s)
    dressed = I(U(outer, *rails, *studs), box(-big, big, -big, big, -ICO_RI, ICO_RI))
    shell = D(dressed, inner)
    lower = I(shell, box(-big, big, -big, big, -big, 0))
    upper = I(shell, box(-big, big, -big, big, 0, big))
    tabs_l, tabs_u, holes = [], [], []
    for k in range(TAB_N):
        a = 2 * math.pi * k / TAB_N
        loc, _, _ = outer.ray.intersects_location([[0, 0, -0.2]], [[math.cos(a), math.sin(a), 0]])
        rw = np.linalg.norm(loc[:, :2], axis=1).max()
        R = rotation_matrix(a, [0, 0, 1])
        pts_l = [(rw - 3, -6, -4), (rw - 3, 6, -4), (rw + 9, -6, -4), (rw + 9, 6, -4),
                 (rw - 3, -6, 0), (rw - 3, 6, 0), (rw + 9, -6, 0), (rw + 9, 6, 0),
                 (rw - 3, -6, -15), (rw - 3, 6, -15)]
        t = trimesh.convex.convex_hull(np.array(pts_l))
        t.apply_transform(R)
        tabs_l.append(t)
        tu = t.copy()
        tu.apply_transform(np.diag([1, 1, -1, 1]))
        tu.fix_normals()
        tabs_u.append(tu)
        holes.append(at(cyl(1.7, -20, 20, 24), (rw + 4.5) * math.cos(a), (rw + 4.5) * math.sin(a)))
    lower = D(U(lower, *tabs_l), inner, *holes)
    upper = D(U(upper, *tabs_u), inner, *holes)
    lower = D(lower, cyl(SLEEVE_HOLE / 2, -ICO_RC, 0))          # inboard: sleeve hole
    boss = cyl(16, ICO_RI - 12, ICO_RI - 1)                     # outboard: hub boss inside
    upper = U(upper, I(boss, ico_solid(1.0)))
    upper = D(upper, cyl(ROD_CLEAR / 2, 0, ICO_RC, 48), hexp(NUT_AF, ICO_RI - NUT_H, ICO_RC))
    return clean(lower), clean(upper)

# ================================================================== bearing blocks
def _x_hole(r, x0, x1, z, n=24):
    """Hole along block-local x (world vertical) at block-local z."""
    h = cyl(r, x0, x1, n)
    h.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))      # +z -> +x
    h.apply_translation([0, 0, z])
    return h

def bearing_block(mouth):
    """Local frame: axle along +z (the block's outboard direction), local +x = world
    up.  Octagonal body, stepped caps on the outboard face (the print's top)."""
    L = BLOCK_LEN
    z0 = FOOT_Z[0]             # the body runs the full foot length so the part stands flat on z0
    parts = [oct_prism(BLOCK_RF, z0, L / 2),
             box(-24, -12, -YOKE_T / 2, YOKE_T / 2, z0, FOOT_Z[1]),
             box(-24, 0, -YOKE_T / 2, YOKE_T / 2, z0, L / 2)]
    sl = X_ICO - (X_RISER + L / 2)
    if mouth:
        parts += [box(14, 21, -6, 6, SEAT_Z[0], SEAT_Z[1]),            # pointer seat
                  cyl(17.2, L / 2 - 0.5, L / 2 + 2), cyl(15.9, L / 2 + 1.9, L / 2 + 4),
                  cyl(SLEEVE_D / 2, L / 2 - 1, L / 2 + sl - 12),
                  hexp(SLEEVE_HEX, L / 2 + sl - 12, L / 2 + sl)]
    else:
        parts += [cyl(15.0, L / 2 - 0.5, L / 2 + 2), cyl(12.5, L / 2 + 1.9, L / 2 + 4)]
    m = U(*parts)
    cut = [cyl(B608_D / 2, z0 - 1, z0 + B608_W),                        # pocket faces the rotor
           cyl(B608_SHOULDER / 2, z0, L / 2 + (0 if mouth else 5))]
    cut += [_x_hole(2.7, -30, -6, zf) for zf in BLOCK_BOLTS]
    if mouth:
        tip = L / 2 + sl
        cut += [cyl(10.5 / 2, L / 2, tip), cyl(B608_D / 2, tip - B608_W, tip + 1)]
        cut += [_x_hole(1.4, 12, 22, zz, 16) for zz in SEAT_SCREWS]     # pointer screws, self-tap M3
        h = cyl(1.4, 0, 20, 16)                                          # paddle set screw, from world below
        h.apply_transform(rotation_matrix(-math.pi / 2, [0, 1, 0]))
        h.apply_translation([0, 0, tip - 6])
        cut.append(h)
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
    """Tube from the nose bearing's inner race to the outboard nut, inside the block bore."""
    length = (BLOCK_LEN / 2 + 4) - (FOOT_Z[0] + B608_W) + 1
    return clean(D(cyl(6.0, 0, length, 48), cyl(4.3, -1, length + 1, 32)))

# ================================================================== paddle
def paddle_profile():
    """Paddle outline in world (x - X_ICO, z - axle)."""
    xs = np.linspace(-PADDLE_HALF, PADDLE_HALF, 41)
    tips = [ico_inner_min_radius(x) - STONE_GAP for x in xs]
    pts = [(x, -t) for x, t in zip(xs, tips)]
    pts = [(xs[-1], -17.0)] + pts[::-1] + [(xs[0], -17.0)]
    return Polygon(pts).buffer(0)

def paddle():
    blade_poly = unary_union([paddle_profile(), sbox(-12, -24, 0, -15)])
    z_st = -(ico_inner_min_radius(0) - STONE_GAP - PADDLE_STONE_D / 2 - 4)
    blade_poly = blade_poly.difference(Point(0.0, z_st).buffer(PADDLE_STONE_D / 2, 64))
    m = extrude_xz(blade_poly, -PADDLE_T / 2, PADDLE_T / 2)
    c = at(along_x(hexp(36, 0, 12)), -12)
    m = U(m, c)
    bore = at(along_x(hexp(SLEEVE_HEX + 0.4, -1, 14)), -13)
    m = D(m, bore, at(cyl(1.7, -30, 0, 16), -6))
    return clean(m)            # world-local; print lying on a y face

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
    cut = [cyl(POST_D / 2 + 0.25, BASE_T, 46), cyl(4.4, -1, BASE_T + 1, 32), hexp(NUT_AF, -1, NUT_H - 1)]
    for sx in (-1, 1):
        for sy in (-1, 1):
            cut.append(at(cyl(3.2, -1, BASE_T + 1, 24), sx * 78, sy * 78))
    return clean(D(m, *cut))

# ================================================================== communication disc
def comm_disc():
    """Spirit-board dial in print frame: X = viewer's right, Y = up, relief on +Z."""
    base_m = cyl(DISC_R, 0, DISC_T, 180)
    hub = U(cyl(17, DISC_T - 0.5, DISC_T + 3, 64), cyl(14, DISC_T + 2.9, DISC_T + 6, 64),
            cyl(11.5, DISC_T + 5.9, DISC_T + 8, 64))
    n = len(WORDS)
    sector = 2 * math.pi / n
    f = [DM.ring(106.0, 107.2), DM.ring(104.5, 105.4)]                   # bezel step, number track
    for k in range(50):                                                  # railroad ticks
        if k % 5:
            f.append(DM.place(sbox(-0.45, 0, 0.45, 3.6), 100.6, k * 2 * math.pi / 50))
    for k, s in enumerate(NUMBERS):
        f += DM.arc_text(s, 95.6, k * 2 * math.pi / len(NUMBERS), 7.2)
    f += [DM.ring(91.0, 92.0), DM.ring(88.6, 89.5)]
    for k in range(n):                                                   # sector dividers
        a = (k + 0.5) * sector
        f.append(DM.place(sbox(-0.5, 57.0, 0.5, 88.0), 0.0, a))
        f.append(DM.place(DM.lozenge((0, 0), 3.8, 6.4), 74.0, a))
    for k, w in enumerate(WORDS):
        f += DM.arc_text(w, 71.0, k * sector, 8.2, max_arc=0.78 * 71.0 * sector)
        f.append(DM.place(SYMBOLS[w](), 62.5, k * sector))
    f += [DM.ring(54.6, 55.6), DM.ring(52.4, 53.3)]
    f.append(DM.sunburst((0, 0), 20.2, 50.5, 43.0, 36, 2.2, 0.95))
    f.append(DM.ring(18.4, 19.4))
    relief = unary_union(f).intersection(Point(0, 0).buffer(105.95, 180)).difference(Point(0, 0).buffer(17.4, 64))
    bezel = DM.ring(107.2 - 0.05, DISC_R, res=180)
    m = U(base_m, hub, *extrude(relief, DISC_T - 0.2, DISC_T + RELIEF),
          *extrude(bezel, DISC_T - 0.2, DISC_T + BEZEL_H))
    m = D(m, cyl(ROD_CLEAR / 2, -1, DISC_T + 9, 48), hexp(NUT_AF, DISC_T + 8 - NUT_H, DISC_T + 9))
    return clean(m)            # print orientation: dial up; world: print z -> +x, print Y -> +z

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
    xc = X_RISER
    body = unary_union([
        sbox(xc - 12, 21, xc + 12, 34),                 # plinth (screws at +-8.5)
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
    cut = [at(cyl(1.7, 20, 36, 16), xc + dx) for dx in SEAT_SCREWS]
    return clean(D(m, *cut))   # world-local (z from axle); print lying on a y face

def finial():
    """Stepped spire pressed over an axle-end M8 nut (collar clears the stationary face by 1 mm)."""
    m = U(cyl(12, 0, 4, 64), oct_prism(10, 3.9, 8), oct_prism(8, 7.9, 12), oct_prism(6, 11.9, 16),
          revolve([(0, 15.9), (4.6, 15.9), (0.9, 26), (0, 26)], 48))
    return clean(D(m, hexp(NUT_AF, -1, 5.5), cyl(4.4, -1, 11, 32)))

# ================================================================== assembly + report
def report(name, m, material, count=1):
    ext = m.extents
    fits = all(e <= BED + 1e-6 for e in ext)
    print(f"{name:14s} x{count}  {material:10s} {ext[0]:6.1f} x {ext[1]:6.1f} x {ext[2]:6.1f} mm  "
          f"watertight={str(m.is_watertight):5s} {m.volume/1000:7.1f} cm3  ~{m.volume*RHO_PETG:5.0f} g solid"
          f"{'' if fits else '  ** EXCEEDS BED **'}")
    assert fits, f"{name} exceeds the bed"

def print_orient(name, m):
    if name in ("yoke_nose", "yoke_mouth", "pointer", "paddle"):
        m = m.copy()
        m.apply_transform(rotation_matrix(math.pi / 2, [1, 0, 0]))    # -y face on the bed
        return to_bed(m)
    if name == "ico_outboard":
        m = m.copy()
        m.apply_transform(rotation_matrix(math.pi, [1, 0, 0]))
        return to_bed(m)
    return to_bed(m)

def assembly(P):
    zax = Z_PAD + Z_AX
    W = [("base", P["base"], "brown"),
         ("post", at(P["post"], z=BASE_T), "brown"),
         ("yaw_hub", at(P["yaw_hub"], z=Z_PAD - YAW_H - PAD_T), "brown")]
    for n in ("yoke_nose", "yoke_mouth"):
        W.append((n, at(P[n], z=Z_PAD), "brown"))
    W.append(("block_nose", block_to_world(P["block_nose"], False), "brown"))
    W.append(("block_mouth", block_to_world(P["block_mouth"], True), "brown"))
    r = P["rotor"].copy()
    r.apply_transform(rotation_matrix(math.pi / 2, [0, 1, 0]))     # print z -> world +x
    W.append(("rotor", at(r, -ROTOR_L / 2, 0, zax), "glass"))
    for n in ("ico_inboard", "ico_outboard"):
        W.append((n, at(along_x(P[n]), X_ICO, 0, zax), "glass"))
    W.append(("paddle", at(P["paddle"], X_ICO, 0, zax), "brown"))
    W.append(("comm_disc", disc_to_world(P["comm_disc"]), "brown"))
    W.append(("pointer", at(P["pointer"], z=zax), "brown"))
    x_cap = X_RISER + BLOCK_LEN / 2 + 4                              # stepped cap face
    fn = P["finial"].copy()
    fn.apply_transform(rotation_matrix(-math.pi / 2, [0, 1, 0]))    # local z -> world -x
    W.append(("finial_nose", at(fn, -(x_cap + 1), 0, zax), "brass"))
    W.append(("finial_disc", at(along_x(P["finial"]), X_DISC + DISC_T + 8 + 1, 0, zax), "brass"))
    x_end = X_DISC + DISC_T + 8 + 6.5
    W.append(("rod", at(along_x(cyl(4, -(x_cap + 6.5), x_end, 24)), 0, 0, zax), "steel"))
    return W

def clearance_checks(W):
    """Rotating parts sweep surfaces of revolution about the axle, so a stationary
    part clears them iff its (x, r) image stays off theirs.  Surface samples of both
    in the (x, r) half-plane, nearest-neighbour distance = clearance (~0.5 mm res.)."""
    from scipy.spatial import cKDTree
    zax = Z_PAD + Z_AX
    parts = {n: m for n, m, _ in W}
    def xr(m, n):
        p, _ = trimesh.sample.sample_surface_even(m, n)
        return np.column_stack([p[:, 0], np.hypot(p[:, 1], p[:, 2] - zax)])
    pairs = [("rotor", ["yoke_nose", "yoke_mouth", "block_nose", "block_mouth", "pointer"]),
             ("ico_inboard", ["block_mouth", "paddle", "pointer", "yoke_mouth"]),
             ("ico_outboard", ["paddle", "pointer"]),
             ("comm_disc", ["pointer", "yoke_mouth"]),
             ("finial_nose", ["block_nose"]),
             ("finial_disc", ["pointer"])]
    worst = {}
    for rot, statics in pairs:
        tree = cKDTree(xr(parts[rot], 80000))
        for st in statics:
            d, _ = tree.query(xr(parts[st], 30000))
            worst[(rot, st)] = d.min()
            print(f"  clearance {rot:13s} <-> {st:12s} {d.min():6.1f} mm")
    return worst

if __name__ == "__main__":
    os.makedirs("stl", exist_ok=True)
    P = {}
    P["rotor"] = rotor()
    P["ico_inboard"], P["ico_outboard"] = ico_halves()
    P["block_nose"] = bearing_block(False)
    P["block_mouth"] = bearing_block(True)
    P["race_spacer"] = race_spacer()
    P["paddle"] = paddle()
    P["yoke_nose"] = yoke(False)
    P["yoke_mouth"] = yoke(True)
    P["yaw_hub"] = yaw_hub()
    P["post"] = post()
    P["base"] = base()
    P["comm_disc"] = comm_disc()
    P["pointer"] = pointer()
    P["finial"] = finial()

    mats = {"rotor": "clear", "ico_inboard": "clear", "ico_outboard": "clear",
            "comm_disc": "brown+brass", "finial": "brass"}
    counts = {"finial": 2}
    print(f"axle height {Z_PAD + Z_AX:.0f} mm, rotor top {Z_PAD + Z_AX + R_MAX:.0f} mm, "
          f"icosahedron centre x={X_ICO:.1f}, disc back x={X_DISC:.1f}")
    for name, m in P.items():
        pm = print_orient(name, m)
        path = os.path.join("stl", f"{name}.stl")
        pm.export(path)
        back = trimesh.load(path)
        report(name, back, mats.get(name, "dark brown"), counts.get(name, 1))
        assert back.is_watertight, f"{name} not watertight after reload"
    print(f"rotor unsupported overhang (>45 deg from vertical), area fraction: "
          f"{overhang_report(P['rotor']):.3f}")
    W = assembly(P)
    worst = clearance_checks(W)
    assert all(v > 0.5 for v in worst.values()), "a stationary part touches a rotating one"
    np.savez_compressed("assembly.npz", **{f"{i:02d}_{n}_{mat}_v": m.vertices for i, (n, m, mat) in enumerate(W)},
                        **{f"{i:02d}_{n}_{mat}_f": m.faces for i, (n, m, mat) in enumerate(W)})
    print("wrote stl/*.stl and assembly.npz")
