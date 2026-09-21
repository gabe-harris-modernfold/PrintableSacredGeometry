#!/usr/bin/env python3
"""Weihnachtspyramide -- candle-driven carousel, parametric, PETG on a 320 mm bed.

Run from this directory:  python pyramid.py
Writes one STL per printable part plus assembly.stl (preview only) and prints a
part table (extents, watertight, volume, mass).

Architecture (see README.md):
  base_plate      square plate, four tea-light rings, four leaning post bosses
  mast            tapered central post (FRAME='mast', default)
  post x4 + apex_collar   leaning posts and crown (FRAME='posts')
  pin_carrier     stub pressed into the collar; holds a 1.5 mm steel pin (or the
                  all-printed pin_printed alternative)
  rotor_hub       vee cup down over the pin, eight angled blade slots
  blade x8        flat plates, printed flat, slid into the hub slots
"""
import math, os, sys
import numpy as np
import trimesh
from trimesh.transformations import rotation_matrix

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from mesh_kit import prism

# ----------------------------------------------------------------- parameters
BED = 320.0                 # printable cube edge, mm
RHO_PETG = 1.27e-3          # g/mm^3

PLATE = 250.0               # square base, edge
PLATE_T = 6.0
PLATE_CORNER = 22.0         # corner chamfer

TEALIGHT_D = 39.0           # standard tea light cup
CUP_ID = TEALIGHT_D + 1.2
CUP_WALL = 2.0
CUP_H = 10.0                # ring height above plate
CUP_FLOOR_HOLE = 22.0       # push-out hole under the cup, also an air gap
R_CANDLE = 88.0             # candle centres, on the axes
FLAME_ABOVE_PLATE = PLATE_T + 16.0 + 30.0   # cup + tea light + flame tip

FRAME = "mast"              # "mast": one central post. "posts": four leaning posts + collar
MAST_D_BASE, MAST_D_TOP = 18.0, 18.0   # constant-section mast (set MAST_D_TOP smaller to taper)
MAST_BOSS_D, MAST_BOSS_TOP_D = 52.0, 30.0   # flared boss, foot to crown
MAST_BOSS_H = 24.0          # boss height above the plate; the mast socket runs the full depth
MAST_SOCKET_DEPTH = PLATE_T + MAST_BOSS_H   # through the plate and the boss

N_POST = 4
POST_D = 9.0
POST_FIT = 0.35             # diametral clearance in sockets
R_POST_BASE = 112.0         # post axis at z=0 (on the diagonals)
R_POST_TOP = 20.0           # post axis at z=Z_APEX
Z_APEX = 300.0
BOSS_D = 18.0
POST_Z0, POST_Z1 = 2.0, 296.0
SPLIT_POSTS = False         # two halves + printed dowel if the post will not fit
DOWEL_D = 5.0

COLLAR_D = 40.0
COLLAR_Z0, COLLAR_Z1 = 268.0, 300.0
CARRIER_D = 6.0             # pin-carrier stub
CARRIER_FIT = 0.25
STEEL_PIN_D = 1.5           # 1.5 mm steel rod / large sewing needle, dome-tipped
PIN_ABOVE_COLLAR = 34.0     # pin tip above collar top; hub bottom clears the carrier flange by ~5 mm

HUB_D = 36.0
HUB_H = 30.0
HUB_BORE = 7.0              # clearance bore for the pin, from below
JEWEL_D = 4.0               # vee-jewel OD if you have one; else the printed vee is the cup
JEWEL_SOCKET_Z = 22.0       # socket floor height in the hub
SOCKET_H = 3.5
N_BLADE = 8
BLADE_T = 1.0
SLOT_FIT = 0.12
PITCH_DEG = 40.0            # blade plane vs. horizontal
ROTOR_R = 120.0             # blade tip radius (240 mm disc, prints flat)
ROOT_CHORD = 28.0
TIP_CHORD = 52.0
TAB_LEN = 10.0              # blade length inside the hub beyond the rim
TIP_SWEEP = 22.0            # aft sweep of the tip centreline, mm
SPAR_W, SPAR_H = 5.0, 0.8   # leading-edge spar strip: width, height below the plate
FENCE_T, FENCE_H, FENCE_LEN = 1.0, 6.0, 30.0   # tip fence wall: thickness, height, span covered
SLOT_INNER_R = 7.5          # slot stops here; wall to bore = 4 mm
SLOT_Z = HUB_H / 2 - 2.0    # slot centre height in the hub

OUT = os.path.dirname(os.path.abspath(__file__))
ENGINE = "manifold"

# ------------------------------------------------------------------ helpers
def cyl(d, z0, z1, sections=64, at=(0, 0)):
    c = trimesh.creation.cylinder(radius=d / 2, height=z1 - z0, sections=sections)
    c.apply_translation([at[0], at[1], (z0 + z1) / 2])
    return c

def cyl_between(p0, p1, d, sections=48):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    return trimesh.creation.cylinder(radius=d / 2, segment=[p0, p1], sections=sections)

def box(sx, sy, sz, at=(0, 0, 0)):
    b = trimesh.creation.box(extents=[sx, sy, sz])
    b.apply_translation(at)
    return b

def union(*ms):
    return trimesh.boolean.union(list(ms), engine=ENGINE)

def diff(a, *bs):
    return trimesh.boolean.difference([a] + list(bs), engine=ENGINE)

def poly_solid(poly_xy, z0, z1):
    v, f = prism(poly_xy, z0, z1)
    return trimesh.Trimesh(np.asarray(v), np.asarray(f), process=True)

def post_axis(k):
    """(point at z=0, direction per mm of z) of post k."""
    a = math.radians(45 + 90 * k)
    u = np.array([math.cos(a), math.sin(a), 0.0])
    p_base = u * R_POST_BASE
    p_top = u * R_POST_TOP + np.array([0, 0, Z_APEX])
    d = p_top - p_base
    return p_base, d / d[2]

def post_point(k, z):
    p0, dz = post_axis(k)
    return p0 + dz * z

# -------------------------------------------------------------------- parts
def base_plate():
    h = PLATE / 2
    c = PLATE_CORNER
    outline = [(-h + c, -h), (h - c, -h), (h, -h + c), (h, h - c),
               (h - c, h), (-h + c, h), (-h, h - c), (-h, -h + c)]
    parts = [poly_solid(outline, 0, PLATE_T)]
    holes = []
    for i in range(4):
        a = math.radians(90 * i)
        at = (R_CANDLE * math.cos(a), R_CANDLE * math.sin(a))
        parts.append(cyl(CUP_ID + 2 * CUP_WALL, PLATE_T - 1, PLATE_T + CUP_H, at=at))
        holes.append(cyl(CUP_ID, PLATE_T, PLATE_T + CUP_H + 1, at=at))
        holes.append(cyl(CUP_FLOOR_HOLE, -1, PLATE_T + 1, at=at))
    if FRAME == "posts":
        for k in range(N_POST):
            parts.append(cyl_between(post_point(k, -1.0), post_point(k, 22.0), BOSS_D))
            holes.append(cyl_between(post_point(k, POST_Z0), post_point(k, 40.0), POST_D + POST_FIT))
    else:
        prof = np.array([(0.0, PLATE_T - 1), (MAST_BOSS_D / 2, PLATE_T - 1),
                         (MAST_BOSS_TOP_D / 2, PLATE_T + MAST_BOSS_H), (0.0, PLATE_T + MAST_BOSS_H)])
        parts.append(trimesh.creation.revolve(prof, sections=96))
        holes.append(cyl(MAST_D_BASE + POST_FIT, -1.0, MAST_SOCKET_DEPTH + 1))   # through-hole
    body = union(*parts)
    body = diff(body, box(PLATE + 40, PLATE + 40, 20, at=(0, 0, -10)))   # flat underside
    return diff(body, *holes)

def post():
    p0, p1 = post_point(0, POST_Z0), post_point(0, POST_Z1)
    L = float(np.linalg.norm(p1 - p0))
    if not SPLIT_POSTS:
        m = trimesh.creation.cylinder(radius=POST_D / 2, height=L, sections=48)
        m.apply_translation([0, 0, L / 2])
        return {"post": m}, L
    half = L / 2
    lower = trimesh.creation.cylinder(radius=POST_D / 2, height=half, sections=48)
    lower.apply_translation([0, 0, half / 2])
    lower = diff(lower, cyl(DOWEL_D + 0.3, half - 12, half + 1))
    return {"post_half": lower, "post_dowel": cyl(DOWEL_D, 0, 24)}, L

def mast():
    """Central post, sits on the table through the base (z=0) up to z=COLLAR_Z1; pin-carrier socket on top."""
    L = COLLAR_Z1
    rb, rt = MAST_D_BASE / 2, MAST_D_TOP / 2
    prof = np.array([(0.0, 0.0), (rb, 0.0), (rb, MAST_SOCKET_DEPTH), (rt, L), (0.0, L)])
    m = trimesh.creation.revolve(prof, sections=64)
    return diff(m, cyl(CARRIER_D + CARRIER_FIT, L - 12, L + 1))

def apex_collar():
    parts = [cyl(COLLAR_D, COLLAR_Z0, COLLAR_Z1)]
    holes = [cyl(CARRIER_D + CARRIER_FIT, COLLAR_Z1 - 12, COLLAR_Z1 + 1)]
    for k in range(N_POST):
        parts.append(cyl_between(post_point(k, COLLAR_Z0 - 2), post_point(k, COLLAR_Z1 - 2), BOSS_D))
        holes.append(cyl_between(post_point(k, COLLAR_Z0 - 20), post_point(k, POST_Z1), POST_D + POST_FIT))
    body = union(*parts)
    body = diff(body, box(200, 200, 40, at=(0, 0, COLLAR_Z1 + 20)))       # flat crown
    body = diff(body, box(200, 200, 40, at=(0, 0, COLLAR_Z0 - 20)))       # flat underside
    body = diff(body, *holes)
    body.apply_translation([0, 0, -COLLAR_Z0])
    return body

def pin_carrier():
    """Stub into the collar + flange; bore for the steel pin. Prints flange-down."""
    body = union(cyl(CARRIER_D, 0, 12), cyl(10.0, 12, 15))
    return diff(body, cyl(STEEL_PIN_D + 0.1, 3, 16))

def pin_printed():
    """All-printed fallback: stub + cone, 0.6 mm tip. Print point-up, slow."""
    cone_h = PIN_ABOVE_COLLAR - 3.0 - 0.3
    cone = trimesh.creation.cone(radius=3.0, height=cone_h, sections=48)
    cone.apply_translation([0, 0, 15])
    tip = trimesh.creation.icosphere(subdivisions=2, radius=0.3)
    tip.apply_translation([0, 0, 15 + cone_h])
    return union(cyl(CARRIER_D, 0, 12), cyl(10.0, 12, 15), cone, tip)

def rotor_hub():
    hub = cyl(HUB_D, 0, HUB_H)
    holes = [cyl(HUB_BORE, -1, JEWEL_SOCKET_Z)]
    # jewel socket whose roof is a 90-degree printed vee that doubles as the cup
    holes.append(cyl(JEWEL_D + 0.2, JEWEL_SOCKET_Z - 0.01, JEWEL_SOCKET_Z + SOCKET_H))
    vee = trimesh.creation.cone(radius=(JEWEL_D + 0.2) / 2, height=(JEWEL_D + 0.2) / 2, sections=48)
    vee.apply_translation([0, 0, JEWEL_SOCKET_Z + SOCKET_H - 0.01])
    holes.append(vee)
    slot_w = BLADE_T + SLOT_FIT
    for i in range(N_BLADE):
        a = math.radians(360 * i / N_BLADE)
        s = box(HUB_D, slot_w, ROOT_CHORD, at=(HUB_D / 2 + SLOT_INNER_R, 0, 0))
        s.apply_transform(rotation_matrix(math.radians(PITCH_DEG), [1, 0, 0]))
        s.apply_transform(rotation_matrix(a, [0, 0, 1]))
        s.apply_translation([0, 0, SLOT_Z])
        holes.append(s)
    return diff(hub, *holes)

def blade_outline():
    """Scimitar planform in blade XY: x radial (x=0 sits at SLOT_INNER_R), y chord.
    Leading edge at -y (the blade moves toward -y), trailing edge at +y. The
    centreline sweeps back toward the trailing edge over the outer half and the
    chord closes to a rounded point, like a falcon primary."""
    L = ROTOR_R - SLOT_INNER_R
    x_tab = HUB_D / 2 - SLOT_INNER_R                # straight tab inside the hub
    xs = np.linspace(x_tab, L, 80)
    t = (xs - x_tab) / (L - x_tab)                  # 0 at hub rim, 1 at tip
    chord = ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * np.sin(np.pi * np.clip(t / 0.85, 0, 1) / 2)
    chord = np.where(t > 0.55, chord * np.sqrt(np.clip(1 - ((t - 0.55) / 0.45) ** 2, 0, 1)), chord)
    chord = np.maximum(chord, 0.6)
    sweep = TIP_SWEEP * np.clip((t - 0.35) / 0.65, 0, 1) ** 2.2   # centreline bends aft
    le = sweep - chord / 2
    te = sweep + chord / 2
    pts = [(0.0, -ROOT_CHORD / 2)] + list(zip(xs, le)) + list(zip(xs[::-1], te[::-1])) + [(0.0, ROOT_CHORD / 2)]
    return pts

def blade():
    """Flat 1 mm plate + leading-edge spar + tip fence, all on the pressure (-z) face so
    the part prints flat, fence up, with no supports. Exported flipped for printing."""
    from shapely.geometry import Polygon, box as sbox
    from shapely.ops import unary_union
    outline = Polygon(blade_outline()).buffer(0)
    plate = trimesh.creation.extrude_polygon(outline, BLADE_T)
    x_tab = HUB_D / 2 - SLOT_INNER_R
    L = ROTOR_R - SLOT_INNER_R
    # leading-edge spar: a SPAR_W strip inside the LE, from just outside the hub to 80 % span
    inner = outline.buffer(-SPAR_W)
    spar_poly = outline.difference(inner).intersection(sbox(x_tab + 1.5, -200, 0.80 * L, 0))
    spar_poly = spar_poly.intersection(outline.buffer(-0.01))
    # tip fence: FENCE_T wall following the outer FENCE_LEN of the outline
    fence_poly = outline.difference(outline.buffer(-FENCE_T)).intersection(sbox(L - FENCE_LEN, -200, L + 1, 200))
    solids = [plate]
    for poly, h in ((spar_poly, SPAR_H), (fence_poly, FENCE_H)):
        geoms = getattr(poly, "geoms", [poly])
        for g in geoms:
            if g.area < 0.5:
                continue
            m = trimesh.creation.extrude_polygon(g.buffer(0.02), h)
            m.apply_translation([0, 0, -h + 0.05])      # hangs below the plate's -z face
            solids.append(m)
    return union(*solids)

def blade_print_orientation(m):
    """Flip so spar and fence point up: prints flat on the plate's smooth face."""
    m = m.copy()
    m.apply_transform(rotation_matrix(math.pi, [1, 0, 0]))
    m.apply_translation([0, 0, -m.bounds[0][2]])
    return m

# ---------------------------------------------------------------- assembly
def rotor_hub_z():
    """Hub bottom in world z: the pin dome seats ~0.8 mm into the vee."""
    return COLLAR_Z1 + PIN_ABOVE_COLLAR - (JEWEL_SOCKET_Z + SOCKET_H + 0.8)

def assembled(parts, post_len):
    world = [("base_plate", parts["base_plate"].copy())]
    if FRAME == "mast":
        world.append(("mast", parts["mast"].copy()))
    for k in range(N_POST if FRAME == "posts" else 0):
        p0, p1 = post_point(k, POST_Z0), post_point(k, POST_Z1)
        m = trimesh.creation.cylinder(radius=POST_D / 2, height=post_len, sections=48)
        m.apply_translation([0, 0, post_len / 2])
        d = (p1 - p0) / np.linalg.norm(p1 - p0)
        m.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], d))
        m.apply_translation(p0)
        world.append((f"post_{k}", m))
    if FRAME == "posts":
        c = parts["apex_collar"].copy(); c.apply_translation([0, 0, COLLAR_Z0]); world.append(("apex_collar", c))
    pc = parts["pin_carrier"].copy(); pc.apply_translation([0, 0, COLLAR_Z1 - 12]); world.append(("pin_carrier", pc))
    world.append(("steel_pin", cyl(STEEL_PIN_D, COLLAR_Z1 - 9, COLLAR_Z1 + PIN_ABOVE_COLLAR, sections=24)))
    z_hub = rotor_hub_z()
    h = parts["rotor_hub"].copy(); h.apply_translation([0, 0, z_hub]); world.append(("rotor_hub", h))
    for i in range(N_BLADE):
        b = parts["blade"].copy()
        b.apply_translation([SLOT_INNER_R, 0, -BLADE_T / 2])
        b.apply_transform(rotation_matrix(math.radians(PITCH_DEG), [1, 0, 0]))
        b.apply_transform(rotation_matrix(math.radians(360 * i / N_BLADE), [0, 0, 1]))
        b.apply_translation([0, 0, z_hub + SLOT_Z])
        world.append((f"blade_{i}", b))
    for i in range(4):   # tea lights, preview only
        a = math.radians(90 * i)
        world.append((f"tealight_{i}", cyl(TEALIGHT_D, PLATE_T, PLATE_T + 16,
                                           at=(R_CANDLE * math.cos(a), R_CANDLE * math.sin(a)))))
    return world, z_hub

def clean(m):
    """Collapse sub-micron slivers left by the booleans so the mesh still reloads
    watertight after STL's float32 rounding."""
    m = m.copy()
    m.merge_vertices(digits_vertex=4)
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    return m

def report(name, m):
    ext = m.extents
    fits = all(e <= BED + 1e-6 for e in ext)
    vol = m.volume
    print(f"{name:14s} {ext[0]:7.1f} x {ext[1]:6.1f} x {ext[2]:6.1f} mm  "
          f"watertight={str(m.is_watertight):5s} vol={vol/1000:7.1f} cm3  mass~{vol*RHO_PETG:6.1f} g"
          f"{'' if fits else '  ** EXCEEDS BED **'}")
    assert m.is_watertight, f"{name} not watertight"
    assert m.volume > 0, f"{name} inverted"

if __name__ == "__main__":
    parts = {"base_plate": base_plate()}
    post_parts, post_len = post()
    if FRAME == "posts":
        parts.update(post_parts)
        parts["apex_collar"] = apex_collar()
    else:
        parts["mast"] = mast()
    parts["pin_carrier"] = pin_carrier()
    parts["pin_printed"] = pin_printed()
    parts["rotor_hub"] = rotor_hub()
    parts["blade"] = blade()
    lean = math.degrees(math.atan((R_POST_BASE - R_POST_TOP) / Z_APEX))
    print(f"frame: {FRAME}" + (f", post length {post_len:.1f} mm, lean {lean:.1f} deg" if FRAME == "posts" else ""))
    for name, m in parts.items():
        parts[name] = m = clean(m)
        report(name, m)
        path = os.path.join(OUT, f"{name}.stl")
        (blade_print_orientation(m) if name == "blade" else m).export(path)
        assert trimesh.load(path).is_watertight, f"{name} not watertight after reload"
    world, z_hub = assembled(parts, post_len)
    scene = trimesh.util.concatenate([m for _, m in world])
    scene.export(os.path.join(OUT, "assembly.stl"))
    z_blade = z_hub + SLOT_Z
    print(f"rotor hub bottom z={z_hub:.1f}, blade plane centre z={z_blade:.1f}, "
          f"flame tip z={FLAME_ABOVE_PLATE:.1f}, blades {z_blade - FLAME_ABOVE_PLATE:.0f} mm above flames")
    print(f"total height {scene.bounds[1][2]:.1f} mm")
    rotor_mass = (parts['rotor_hub'].volume + N_BLADE * parts['blade'].volume) * RHO_PETG
    print(f"rotor mass (solid hub) ~ {rotor_mass:.0f} g; blade {parts['blade'].volume*RHO_PETG:.1f} g each")
