"""
Dodecahedron mirror-ball frame generator.

12 identical printable bezel tiles that 12 pentagonal mirrors drop into. Each
mirror is cut from one 10 cm-side regular hexagonal mirror tile (the 20 x 17.32
x 10 cm spec).

Tiles mate on mitre faces cut at half the 116.565 deg dihedral. Alignment is by
a half-round groove running along each shared edge: two mating tiles form a
full 4 mm bore. That works because a dodecahedron has a 2-fold rotation axis
through every edge midpoint, and a groove whose axis is parallel to the edge and
centred on the mitre plane maps onto itself under that rotation.

Output
  dodeca-mirror-tile.stl      one tile, oriented for support-free printing (x12)
  dodeca-mirror-pin.stl       one alignment pin (x30, or use 4 mm dowel/rod)
  dodeca-mirror-assembly.stl  preview of the finished ball (do not print)
  pentagon-cut-template.svg   1:1 cutting template for the hexagons
"""

import math
import numpy as np
import trimesh
from trimesh.creation import extrude_polygon, cylinder
from shapely.geometry import Polygon

# ---------------------------------------------------------------- parameters

MIRROR_SIDE = 100.0   # mm, pentagon side you cut from each hexagon
MIRROR_THK  = 4.0     # mm, mirror glass thickness
CLEAR       = 1.2     # mm, radial pocket clearance (accepts 100.0 - 101.7).
                      # 0.6 was too tight: it left only +0.87 mm of oversize
                      # for glass cut by hand to a traced template, where
                      # +/-1 mm per edge is normal and all five must fit.
                      # Costs 0.6 mm of rim wall (3.40 -> 2.80 mm, still ~7
                      # perimeters at 0.42 mm) and changes nothing else.
RIB         = 5.0     # mm, visible frame rib per side (inradius inset).
                      # Must exceed CLEAR + MIRROR_THK*face_inr/r_i + ~1.2.
                      # The mitre leans in 0.618 mm per mm of depth, so the
                      # pocket wall is a WEDGE: it loses 2.47 mm of width over
                      # the 4 mm pocket. At RIB 4.0 / CLEAR 1.2 the wall bottom
                      # was 0.33 mm -- under one extrusion width, so the
                      # slicer would drop or break the first traces of the rim.
LEDGE       = 6.0     # mm, ledge reach under the mirror
THK         = 14.0    # mm, frame thickness (radial)

PIN_D       = 4.0     # mm, pin diameter
PIN_LEN     = 30.0    # mm, pin length
GROOVE_R    = 2.15    # mm, groove radius (0.3 mm total clearance on the pin)
GROOVE_LEN  = 34.0    # mm, groove length along the shared edge
GROOVE_Z    = 9.0     # mm, groove axis depth below the face plane

ENGINE = "manifold"

# ------------------------------------------------------ pentagon / dodeca math

P_INR = 1.0 / (2.0 * math.tan(math.radians(36.0)))   # side -> inradius
P_CIR = 1.0 / (2.0 * math.sin(math.radians(36.0)))   # side -> circumradius
D_INR = 1.1135163644             # dodeca edge -> inradius
D_CIR = 1.4012585384             # dodeca edge -> circumradius
DIHEDRAL = 116.5650512           # degrees

mirror_inr = MIRROR_SIDE * P_INR
face_inr   = mirror_inr + RIB
edge_a     = face_inr / P_INR
r_i        = edge_a * D_INR
r_c        = edge_a * D_CIR
pocket_inr = mirror_inr + CLEAR
open_inr   = pocket_inr - LEDGE


def max_pentagon_in_hexagon(hex_side):
    """Largest regular pentagon centred in a regular hexagon, over all rotations."""
    hex_inr = hex_side * math.sqrt(3) / 2.0
    normals = [math.radians(30 + 60 * j) for j in range(6)]
    best = (0.0, 0.0)
    for deg in np.linspace(0, 72, 72001):
        th = math.radians(deg)
        worst = max(math.cos(th + math.radians(72 * k) - n)
                    for k in range(5) for n in normals)
        r = hex_inr / worst
        if r > best[0]:
            best = (r, deg)
    return best[0] / P_CIR, best[1]


def pent_pts(inradius, z, rot_deg):
    """5 body-coordinate points of a regular pentagon at height z."""
    rc = inradius / P_INR * P_CIR
    return np.array([
        [rc * math.cos(math.radians(rot_deg + 72 * k)),
         rc * math.sin(math.radians(rot_deg + 72 * k)), z] for k in range(5)])


def prism(inradius, z0, z1, rot_deg=0.0):
    """Straight pentagonal prism between two heights (body coordinates)."""
    p = pent_pts(inradius, 0.0, rot_deg)
    m = extrude_polygon(Polygon(p[:, :2]), height=z1 - z0)
    m.apply_translation([0, 0, z0])
    return m


def cone(inradius, rot_deg, grow=1.45):
    """Pyramid from the body centre through a face-plane pentagon.

    Because the apex is exactly at the body centre, this tapers at the same
    rate as the shell, so a ledge cut with it keeps a constant width at every
    depth. A straight prism would instead be overtaken by the mitre taper.
    """
    base = pent_pts(inradius, r_i, rot_deg)
    return trimesh.convex.convex_hull(np.vstack([np.zeros(3), base * grow]))


# ------------------------------------------------------------ base dodecahedron

def dodecahedron(edge):
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    inv = 1.0 / phi
    v = [(sx, sy, sz) for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]
    for s1 in (1, -1):
        for s2 in (1, -1):
            v += [(0, s1 * inv, s2 * phi), (s1 * inv, s2 * phi, 0),
                  (s1 * phi, 0, s2 * inv)]
    return trimesh.convex.convex_hull(np.array(v, float) * (edge / (2.0 / phi)))


def face_normals_unique(mesh, tol=1e-6):
    out = []
    for n in mesh.face_normals:
        if not any(np.linalg.norm(n - m) < tol for m in out):
            out.append(n)
    return np.array(out)


def build_oriented():
    """Dodecahedron centred at the origin, face 0's outward normal on +Z."""
    d = dodecahedron(edge_a)
    d.apply_translation(-d.centroid)
    d.apply_transform(trimesh.geometry.align_vectors(
        face_normals_unique(d)[0], [0, 0, 1.0]))
    top = d.vertices[d.vertices[:, 2] > r_i - 1e-4]
    ang = math.degrees(math.atan2(top[0][1], top[0][0]))
    d.apply_transform(trimesh.transformations.rotation_matrix(
        math.radians(90.0 - ang), [0, 0, 1]))
    top = d.vertices[d.vertices[:, 2] > r_i - 1e-4]
    order = np.argsort([math.atan2(p[1], p[0]) for p in top])
    return d, face_normals_unique(d), top[order]


# ------------------------------------------------------------------- tile build

def edge_frame(fverts, k):
    """(midpoint, edge unit vector) for edge k of the face, body coordinates."""
    a = np.array(fverts[k], float)
    b = np.array(fverts[(k + 1) % 5], float)
    e = (b - a) / np.linalg.norm(b - a)
    return (a + b) / 2.0, e


def build_tile(solid, fverts, rot_deg):
    inner = solid.copy()
    inner.apply_scale((r_i - THK) / r_i)
    shell = trimesh.boolean.difference([solid, inner], engine=ENGINE)

    # the wedge's lateral faces are exactly the 5 mitre planes, because its
    # apex sits at the body centre
    wedge = trimesh.convex.convex_hull(
        np.vstack([np.zeros(3), np.array(fverts, float) * 1.45]))
    tile = trimesh.boolean.intersection([shell, wedge], engine=ENGINE)

    pocket = prism(pocket_inr, r_i - MIRROR_THK, r_i + 12.0, rot_deg)
    tile = trimesh.boolean.difference(
        [tile, pocket, cone(open_inr, rot_deg)], engine=ENGINE)

    grooves = []
    for k in range(5):
        mid, e = edge_frame(fverts, k)
        p = mid * (1.0 - GROOVE_Z / r_i)      # on the mitre plane, GROOVE_Z deep
        grooves.append(cylinder(radius=GROOVE_R, sections=48,
                                segment=[p - e * GROOVE_LEN / 2,
                                         p + e * GROOVE_LEN / 2]))
    tile = trimesh.boolean.difference([tile] + grooves, engine=ENGINE)

    tile.apply_translation([0, 0, -r_i])       # local: face plane at z = 0
    return tile


def build_pin():
    r, L, ch = PIN_D / 2.0, PIN_LEN, 0.6
    prof = [(0, 0), (r - ch, 0), (r, ch), (r, L - ch), (r - ch, L), (0, L)]
    return trimesh.creation.revolve(prof, sections=48)


# --------------------------------------------------------- symmetry / assembly

def ordered_face_verts(mesh, n):
    sel = mesh.vertices[mesh.vertices @ n > r_i - 1e-3]
    c = sel.mean(axis=0)
    u = sel[0] - c
    u = u / np.linalg.norm(u)
    w = np.cross(n, u)
    return sel[np.argsort([math.atan2(np.dot(p - c, w), np.dot(p - c, u))
                           for p in sel])]


def kabsch(A, B):
    U, _, Vt = np.linalg.svd(A.T @ B)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    return Vt.T @ D @ U.T


def face_rotations(mesh, normals):
    """Exact dodecahedral rotations carrying face 0 onto every face."""
    V0 = ordered_face_verts(mesh, normals[0])
    out = []
    for n in normals:
        Vk = ordered_face_verts(mesh, n)
        best = (1e9, None)
        for s in range(5):
            Vs = np.roll(Vk, s, axis=0)
            R = kabsch(V0, Vs)
            err = float(np.linalg.norm((R @ V0.T).T - Vs))
            if err < best[0]:
                best = (err, R)
        T = np.eye(4)
        T[:3, :3] = best[1]
        out.append((best[0], T))
    return out


# ------------------------------------------------------------------ svg template

def write_svg(path, hex_side=100.0, pent_side=MIRROR_SIDE):
    W, H = 297.0, 210.0
    cx, cy = W / 2.0, H / 2.0
    hexpts = [(cx + hex_side * math.cos(math.radians(60 * k)),
               cy - hex_side * math.sin(math.radians(60 * k))) for k in range(6)]
    rc = pent_side * P_CIR
    pentpts = [(cx + rc * math.cos(math.radians(72 * k)),
                cy - rc * math.sin(math.radians(72 * k))) for k in range(5)]
    f = lambda p: " ".join(f"{x:.3f},{y:.3f}" for x, y in p)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm"
     viewBox="0 0 {W} {H}">
  <rect x="0" y="0" width="{W}" height="{H}" fill="#fff"/>
  <polygon points="{f(hexpts)}" fill="none" stroke="#999"
           stroke-width="0.4" stroke-dasharray="3,2"/>
  <polygon points="{f(pentpts)}" fill="none" stroke="#000" stroke-width="0.6"/>
  <circle cx="{cx:.3f}" cy="{cy:.3f}" r="1" fill="#000"/>
  <g font-family="sans-serif" font-size="4.2" fill="#000">
    <text x="8" y="14">PENTAGON CUT TEMPLATE - print at 100% / actual size, A4 or Letter LANDSCAPE</text>
    <text x="8" y="21">Solid line = cut. Dashed = your 10 cm hexagon, for registration.</text>
    <text x="8" y="28">Align one pentagon corner to one hexagon corner, as drawn. Cut all 5 lines.</text>
    <text x="8" y="{H - 17:.1f}">Pentagon side {pent_side:.1f} mm  |  hexagon side {hex_side:.1f} mm  |  12 needed</text>
  </g>
  <g stroke="#000" stroke-width="0.4" font-family="sans-serif" font-size="4.2">
    <line x1="8" y1="{H - 8:.1f}" x2="108" y2="{H - 8:.1f}"/>
    <line x1="8" y1="{H - 11:.1f}" x2="8" y2="{H - 5:.1f}"/>
    <line x1="108" y1="{H - 11:.1f}" x2="108" y2="{H - 5:.1f}"/>
    <text x="112" y="{H - 6.5:.1f}" stroke="none">100 mm - check with a ruler before cutting glass</text>
  </g>
</svg>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(svg)


# ------------------------------------------------------------------------- main

if __name__ == "__main__":
    import os, sys
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)

    mx, mxrot = max_pentagon_in_hexagon(100.0)
    print(f"max pentagon in 100 mm hex : {mx:.2f} mm side @ {mxrot:.0f} deg")
    print(f"designing for              : {MIRROR_SIDE:.1f} mm  "
          f"({mx - MIRROR_SIDE:.1f} mm of cutting margin)")
    print(f"dodecahedron edge          : {edge_a:.2f} mm")
    print(f"ball dia face-face         : {2*r_i:.1f} mm")
    print(f"ball dia point-point       : {2*r_c:.1f} mm")
    print(f"dihedral / mitre           : {DIHEDRAL:.3f} / {DIHEDRAL/2:.3f} deg")

    solid, normals, fverts = build_oriented()
    rot_deg = math.degrees(math.atan2(fverts[0][1], fverts[0][0]))
    tile = build_tile(solid, fverts, rot_deg)

    ext = tile.extents
    print(f"\ntile watertight            : {tile.is_watertight}")
    print(f"tile volume                : {tile.volume/1000:.1f} cm3  "
          f"(x12 = {12*tile.volume/1000:.0f} cm3 solid)")
    print(f"tile bbox                  : {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm")

    # The pocket wall is a wedge, not a parallel wall: the mitre plane leans
    # inward at face_inr/r_i per mm of depth, so the wall is widest at the
    # mirror face and narrowest where it meets the ledge. The narrow end is
    # what the slicer has to print first, so that is the number that matters.
    taper = face_inr / r_i
    rim_top = face_inr - pocket_inr
    rim_base = rim_top - MIRROR_THK * taper
    EXTRUSION = 0.42
    print(f"\nmitre taper                : {taper:.4f} mm/mm "
          f"({math.degrees(math.atan(taper)):.2f} deg from vertical)")
    print(f"pocket wall at mirror face : {rim_top:.2f} mm")
    print(f"pocket wall at the ledge   : {rim_base:.2f} mm "
          f"({rim_base/EXTRUSION:.1f} traces at {EXTRUSION} mm) -> "
          f"{'PASS' if rim_base >= 3*EXTRUSION else 'FAIL, too thin to print'}")
    print(f"glass the pocket accepts   : {MIRROR_SIDE:.1f} to "
          f"{(mirror_inr + CLEAR)/P_INR:.2f} mm side (+{CLEAR/P_INR:.2f})")
    print(f"bed contact ring           : "
          f"{face_inr*(1-THK/r_i) - open_inr*(1-THK/r_i):.2f} mm wide")

    mirror = prism(mirror_inr, r_i - MIRROR_THK, r_i, rot_deg)
    mloc = mirror.copy(); mloc.apply_translation([0, 0, -r_i])
    clash = abs(trimesh.boolean.intersection([tile, mloc], engine=ENGINE).volume)
    print(f"frame/mirror clash         : {clash:.2f} mm3  "
          f"-> {'PASS' if clash < 1 else 'FAIL, mirror will not seat'}")

    # support-free orientation: inner surface flat on the bed, pocket opening up.
    # Both walls lean out 31.7 deg from vertical and the ledge grows off the bed.
    flat = tile.copy()
    flat.apply_translation([0, 0, -flat.bounds[0][2]])
    flat.export(os.path.join(out, "dodeca-mirror-tile.stl"))
    build_pin().export(os.path.join(out, "dodeca-mirror-pin.stl"))

    rots = face_rotations(solid, normals)
    print(f"symmetry residual          : {max(r[0] for r in rots):.1e} mm")
    body = tile.copy(); body.apply_translation([0, 0, r_i])
    frames, mirrors = [], []
    for _, T in rots:
        a = body.copy();   a.apply_transform(T); frames.append(a)
        b = mirror.copy(); b.apply_transform(T); mirrors.append(b)

    worst = max(abs(trimesh.boolean.intersection(
        [frames[0], frames[i]], engine=ENGINE).volume)
        for i in range(1, len(frames)))
    print(f"tile-tile interference     : {worst:.3f} mm3  "
          f"-> {'PASS' if worst < 1 else 'FAIL'}")

    asm = trimesh.util.concatenate(frames + mirrors)
    asm.export(os.path.join(out, "dodeca-mirror-assembly.stl"))
    trimesh.util.concatenate(frames).export(os.path.join(out, "_frames_only.stl"))
    print(f"assembly bbox              : {np.round(asm.extents, 1)}")

    glass = 12 * mirror.volume / 1000.0
    print(f"glass volume               : {glass:.0f} cm3 -> {glass*2.5/1000:.2f} kg at 2.5 g/cm3")

    write_svg(os.path.join(out, "pentagon-cut-template.svg"))
    print("\nwrote:", os.path.abspath(out))
