#!/usr/bin/env python3
"""Concept massings for the Liber Ventorum hanging wind harp: four proposals.

Brief (2026-11-05): centre-radiating, hangs from a hook, a counterweight steadies it against the wind,
clear PETG tubes, steel music strings, Helmholtz chambers, and the octahedron present in the form.

These are massing models for looking at, not printable parts. Run from this directory:

    python concepts.py

Writes concept_a.png .. concept_d.png (hero renders), concepts_lineup.png (all four at one scale,
hung from one beam) and prints the numbers quoted in the research doc.
"""
import itertools
import math

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

pv.OFF_SCREEN = True

C_SOUND = 343.0
G = 9.81
RHO_STEEL = 7850.0
NU_AIR = 1.5e-5
SQ2, SQ3 = math.sqrt(2), math.sqrt(3)

BG = "#f6f5f1"
CLEAR = dict(color="#bfdde9", opacity=0.42, specular=0.9, specular_power=50, smooth_shading=True)
CHAMBER = dict(color="#d9a441", opacity=0.72, specular=0.5, specular_power=25)
OCTA = dict(color="#d9a441", opacity=0.72, specular=0.5, specular_power=25,
            show_edges=True, edge_color="#7a5714", line_width=2)
STEEL = dict(color="#3d434b", specular=0.9, specular_power=60, smooth_shading=True)
BRONZE = dict(color="#5e4630", specular=0.6, specular_power=30, show_edges=True,
              edge_color="#3a2a1b", line_width=1)
NODE = dict(color="#e9eef0", opacity=0.85, specular=0.6, smooth_shading=True)
CORD = dict(color="#262626", smooth_shading=True)


# ---------------------------------------------------------------- primitives (mm, z up)
def octahedron(center=(0, 0, 0), r=100.0, rot=None):
    v = np.array([[r, 0, 0], [-r, 0, 0], [0, r, 0], [0, -r, 0], [0, 0, r], [0, 0, -r]], float)
    if rot is not None:
        v = v @ rot.T
    v += np.asarray(center, float)
    tris = [(4, 0, 2), (4, 2, 1), (4, 1, 3), (4, 3, 0), (5, 2, 0), (5, 1, 2), (5, 3, 1), (5, 0, 3)]
    faces = np.array([[3, *t] for t in tris]).ravel()
    return pv.PolyData(v, faces)


def tube(p0, p1, r):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    length = float(np.linalg.norm(d))
    return pv.Cylinder(center=(p0 + p1) / 2, direction=d / length, radius=r, height=length,
                       resolution=40, capping=True)


def string(p0, p1, r=1.3):
    return pv.Line(np.asarray(p0, float), np.asarray(p1, float)).tube(radius=r, n_sides=10)


def ball(c, r):
    return pv.Sphere(radius=r, center=np.asarray(c, float), theta_resolution=24, phi_resolution=24)


def hook(z_attach, z_ring):
    """Cord from the piece up to a ring hanging on a hook (the hook is drawn as the ring)."""
    ring = pv.ParametricTorus(ringradius=26, crosssectionradius=4.5)
    ring = ring.rotate_x(90, inplace=False).translate((0, 0, z_ring), inplace=False)
    return [(string((0, 0, z_attach), (0, 0, z_ring - 26), 2.0), CORD), (ring, STEEL)]


def align_to_z(v):
    """Rotation matrix taking unit vector v onto +z (Rodrigues)."""
    v = np.asarray(v, float) / np.linalg.norm(v)
    z = np.array([0.0, 0.0, 1.0])
    k = np.cross(v, z)
    s, c = np.linalg.norm(k), float(v @ z)
    if s < 1e-12:
        return np.eye(3)
    k /= s
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * kx + (1 - c) * (kx @ kx)


def counterweight(z_top, r=120):
    """Octahedral bob: printed PETG shell filled with steel shot, hung by its top vertex."""
    return [(octahedron((0, 0, z_top - r), r), BRONZE)]


# ---------------------------------------------------------------- the four concepts
def concept_a():
    """Plumb Octahedron: the 12 strings ARE the octahedron's edges; the bob tensions them."""
    parts = []
    rr = 609.5 / SQ2                      # 24 in edges, the patch's string length
    top, bot = np.array([0, 0, rr]), np.array([0, 0, -rr])
    tips = [np.array(p, float) * rr for p in [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]]
    hub = 125.0
    parts.append((octahedron(r=hub), OCTA))
    for sx, sy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:          # necks on the lower faces
        n = np.array([sx, sy, -1.0]) / SQ3
        c0 = np.array([sx, sy, -1.0]) * hub / 3
        parts.append((tube(c0, c0 + n * 26, 15), CHAMBER))
    for p in tips:
        u = p / np.linalg.norm(p)
        parts.append((tube(u * hub, p - u * 18, 15), CLEAR))
        parts.append((ball(p, 20), NODE))
    parts.append((tube((0, 0, hub), top - [0, 0, 18], 15), CLEAR))
    parts.append((ball(top, 20), NODE))
    parts.append((ball(bot, 10), STEEL))
    for i, p in enumerate(tips):
        q = tips[(i + 1) % 4]
        parts += [(string(top, p), STEEL), (string(p, q), STEEL), (string(p, bot), STEEL)]
    parts.append((string(bot, bot - [0, 0, 40], 2.0), CORD))
    parts += counterweight(bot[2] - 40)
    parts += hook(top[2] + 20, top[2] + 110)
    return parts, top[2] + 110


def concept_b():
    """Axis Monochord: one long vertical string, hook to bob, threaded through the body."""
    parts = []
    r = 130.0
    parts.append((octahedron(r=r), OCTA))
    for s in itertools.product((1, -1), repeat=3):                # 8 face tubes = 8 winds
        n = np.array(s, float) / SQ3
        c0 = np.array(s, float) * r / 3
        parts.append((tube(c0, c0 + n * 300, 15), CLEAR))
    for a in [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]:      # necks at the waist
        a = np.array(a, float)
        parts.append((tube(a * (r - 25), a * (r + 22), 13), CHAMBER))
    z_hook, z_bob = 1000.0, -800.0
    parts.append((string((0, 0, z_hook - 26), (0, 0, r + 8), 1.4), STEEL))
    parts.append((string((0, 0, -r - 8), (0, 0, z_bob), 1.4), STEEL))
    parts.append((tube((0, 0, r), (0, 0, r + 8), 20), NODE))       # bridges at the two vertices
    parts.append((tube((0, 0, -r), (0, 0, -r - 8), 20), NODE))
    ring = pv.ParametricTorus(ringradius=26, crosssectionradius=4.5)
    parts.append((ring.rotate_x(90, inplace=False).translate((0, 0, z_hook), inplace=False), STEEL))
    parts += counterweight(z_bob)
    return parts, z_hook


def tensegrity_nodes(b=150.0, c=300.0):
    """Expanded octahedron at equilibrium: nodes on Jessen's icosahedron (0, ±b, ±2b) cyclic."""
    nodes = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            nodes += [(0, s1 * b, s2 * c), (s2 * c, 0, s1 * b), (s1 * b, s2 * c, 0)]
    nodes = np.array(nodes, float)
    struts, cables = [], []
    cable_len = math.sqrt(6) * b
    for i, j in itertools.combinations(range(len(nodes)), 2):
        d = float(np.linalg.norm(nodes[i] - nodes[j]))
        if abs(d - 2 * c) < 1e-6 and np.count_nonzero(np.abs(nodes[i] - nodes[j]) > 1e-6) == 1:
            struts.append((i, j))
        elif abs(d - cable_len) < 1e-6:
            cables.append((i, j))
    return nodes, struts, cables


def concept_c():
    """Tonos Tensegrity: 6 clear struts, 24 steel strings at one shared tension."""
    parts = []
    b, c = 150.0, 300.0
    nodes, struts, cables = tensegrity_nodes(b, c)
    assert len(struts) == 6 and len(cables) == 24, (len(struts), len(cables))
    rot = align_to_z((0, b, c))                                   # hang from node (0, b, c)
    nodes = nodes @ rot.T
    for i, j in struts:
        u = (nodes[j] - nodes[i]) / np.linalg.norm(nodes[j] - nodes[i])
        parts.append((tube(nodes[i] + u * 14, nodes[j] - u * 14, 20), CLEAR))
    for i, j in cables:
        parts.append((string(nodes[i], nodes[j], 1.2), STEEL))
    for p in nodes:
        parts.append((ball(p, 16), NODE))
    cham = 90.0                                                    # chamber in the still centre
    parts.append((octahedron(r=cham, rot=rot), OCTA))
    for a in np.vstack([np.eye(3), -np.eye(3)]):
        parts.append((string((a * cham) @ rot.T, (a * b) @ rot.T, 0.9), STEEL))
    ztop, zbot = nodes[:, 2].max(), nodes[:, 2].min()
    parts.append((string((0, 0, zbot), (0, 0, zbot - 60), 2.0), CORD))
    parts += counterweight(zbot - 60)
    parts += hook(ztop, ztop + 100)
    return parts, ztop + 100


def concept_d():
    """Star of Valihas: 12 tube-zither spokes radiating from an octahedral chamber."""
    parts = []
    r = 110.0
    parts.append((octahedron(r=r), OCTA))
    spokes = [(np.array(a, float) * r, np.array(a, float) * (r + 300))
              for a in [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]]
    for s in itertools.product((1, -1), repeat=3):
        u = np.array(s, float) / SQ3
        c0 = np.array(s, float) * r / 3
        spokes.append((c0, c0 + u * 300))
    for p0, p1 in spokes:
        parts.append((tube(p0, p1, 15), CLEAR))
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        side = np.cross(u, [0, 0, 1.0])
        side /= np.linalg.norm(side)
        perp = np.cross(u, side)
        for sgn in (1, -1):
            off = perp * sgn * 23
            a, z = p0 + u * 40 + off, p1 - u * 14 + off
            parts.append((string(a, z, 1.1), STEEL))
            parts.append((tube(p0 + u * 40 + perp * sgn * 14, a, 3.5), NODE))   # bridges
            parts.append((tube(p1 - u * 14 + perp * sgn * 14, z, 3.5), NODE))
        parts.append((ball(p1, 17), NODE))
    parts.append((tube((0, 0, r), (0, 0, r + 150), 12), CLEAR))
    parts.append((string((0, 0, -r), (0, 0, -r - 320), 2.0), CORD))
    parts += counterweight(-r - 320)
    parts += hook(r + 150, r + 230)
    return parts, r + 230


CONCEPTS = [("a", "A  Plumb Octahedron", concept_a), ("b", "B  Axis Monochord", concept_b),
            ("c", "C  Tonos Tensegrity", concept_c), ("d", "D  Star of Valihas", concept_d)]


# ---------------------------------------------------------------- rendering
def view(parts, out, size, focal, half_height, az=32.0, el=14.0):
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.set_background(BG)
    pl.enable_depth_peeling(number_of_peels=10)
    for mesh, style in parts:
        pl.add_mesh(mesh, **style)
    a, e = math.radians(az), math.radians(el)
    d = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    focal = np.asarray(focal, float)
    pl.enable_parallel_projection()
    pl.camera.focal_point = focal
    pl.camera.position = focal + d * 6000
    pl.camera.up = (0, 0, 1)
    pl.camera.parallel_scale = half_height
    pl.enable_anti_aliasing("ssaa")
    pl.screenshot(out)
    pl.close()
    print("wrote", out)


def bounds_of(parts):
    b = np.array([m.bounds for m, _ in parts])
    return np.array([b[:, 0].min(), b[:, 1].max(), b[:, 2].min(), b[:, 3].max(),
                     b[:, 4].min(), b[:, 5].max()])


def font(px):
    for name in ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, px)
        except OSError:
            continue
    return ImageFont.load_default()


def render_all():
    for key, title, fn in CONCEPTS:
        parts, _ = fn()
        bb = bounds_of(parts)
        focal = [(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2, (bb[4] + bb[5]) / 2]
        half = 0.56 * max(bb[5] - bb[4], 1.15 * (bb[1] - bb[0]))
        if key != "b":
            view(parts, f"concept_{key}.png", (1200, 1400), focal, half)
            continue
        # the monochord is 2 m tall: whole piece on the left, the body close up on the right
        view(parts, "_b_full.png", (520, 1400), focal, half)
        view(parts, "_b_body.png", (1080, 1400), (0, 0, 0), 470)
        sheet = Image.new("RGB", (1600, 1400), BG)
        sheet.paste(Image.open("_b_full.png").convert("RGB"), (0, 0))
        sheet.paste(Image.open("_b_body.png").convert("RGB"), (520, 0))
        ImageDraw.Draw(sheet).line([(520, 80), (520, 1320)], fill="#c9c7c0", width=2)
        sheet.save("concept_b.png")
        print("wrote concept_b.png")

    # lineup: same scale, every hook ring at the same height, one beam above
    beam_z, half, size = 1100.0, 1250.0, (620, 1700)
    tiles = []
    for key, title, fn in CONCEPTS:
        parts, z_ring = fn()
        shift = beam_z - z_ring
        parts = [(m.translate((0, 0, shift), inplace=False), s) for m, s in parts]
        view(parts, f"_lineup_{key}.png", size, (0, 0, beam_z - half + 120), half, az=32, el=8)
        tiles.append((Image.open(f"_lineup_{key}.png").convert("RGB"), title))
    w, h = size
    sheet = Image.new("RGB", (w * 4, h + 90), BG)
    draw = ImageDraw.Draw(sheet)
    px_per_mm = h / (2 * half)
    beam_y = int(h / 2 - (beam_z - (beam_z - half + 120)) * px_per_mm)
    for i, (im, title) in enumerate(tiles):
        sheet.paste(im, (i * w, 0))
        draw.text((i * w + 40, h + 22), title, fill="#222222", font=font(34))
    draw.rectangle([20, beam_y - 14, w * 4 - 20, beam_y - 4], fill="#8c8c8c")
    x0, y1 = w * 4 - 70, h - 60                                   # 1 m scale bar
    y0 = int(y1 - 1000 * px_per_mm)
    draw.rectangle([x0, y0, x0 + 8, y1], fill="#222222")
    draw.text((x0 - 70, (y0 + y1) // 2 - 18), "1 m", fill="#222222", font=font(32))
    sheet.save("concepts_lineup.png")
    print("wrote concepts_lineup.png")


# ---------------------------------------------------------------- numbers for the doc
def mu(d_mm):
    return RHO_STEEL * math.pi * (d_mm / 1000) ** 2 / 4


def f1(L, T, d_mm):
    return math.sqrt(T / mu(d_mm)) / (2 * L)


def helmholtz(vol_l, neck_d_mm=30.0, neck_len_mm=20.0, n=1):
    a = n * math.pi * (neck_d_mm / 2000) ** 2
    leff = (neck_len_mm + 1.5 * neck_d_mm / 2) / 1000
    return C_SOUND / (2 * math.pi) * math.sqrt(a / (vol_l / 1000 * leff))


def octa_inner_litres(r_vertex_mm, wall_mm=2.4):
    ri = r_vertex_mm - wall_mm * SQ3
    return 4 / 3 * ri ** 3 / 1e6


def strouhal(re):
    return -3.3265 / re + 0.1816 + 1.6e-4 * re if re < 178 else 0.19


def cents(r):
    return 1200 * math.log2(r)


def numbers():
    cw = 10.0                                   # kg of counterweight
    print("\n== A Plumb Octahedron")
    rr = 609.5 / SQ2
    t_low = cw * G / 4 / math.cos(math.radians(45))
    print(f"vertex radius {rr:.0f} mm, height vertex-to-vertex {2 * rr:.0f} mm, edge 609.5 mm")
    print(f"tube hub->tip {rr - 125:.0f} mm (one print)")
    print(f"lower-string tension per 10 kg bob: {t_low:.1f} N")
    for d in (0.3, 0.4, 0.5):
        print(f"  steel {d} mm, 0.6095 m, {t_low:.1f} N -> f1 {f1(0.6095, t_low, d):.0f} Hz")
    print(f"bob for 33-55 N: {33 / (G / 4 / math.cos(math.radians(45))):.1f}-"
          f"{55 / (G / 4 / math.cos(math.radians(45))):.1f} kg; octave needs 4x the mass")
    v = octa_inner_litres(125)
    print(f"hub inner volume {v:.2f} L -> Helmholtz 1 neck {helmholtz(v):.0f} Hz, "
          f"4 necks {helmholtz(v, n=4):.0f} Hz (lumped)")

    print("\n== B Axis Monochord")
    body = 2.5
    for name, L, T in (("upper", 0.87, (cw + body) * G), ("lower", 0.67, cw * G)):
        for d in (0.5, 1.0):
            f = f1(L, T, d)
            re = 3.0 * d / 1000 / NU_AIR
            fw = strouhal(re) * 3.0 / (d / 1000)
            print(f"  {name} {L} m, {T:.0f} N, {d} mm: f1 {f:.0f} Hz; wind tone at 3 m/s {fw:.0f} Hz "
                  f"-> harmonic {round(fw / f)}")
    print(f"  interval upper/lower (0.5 mm): {cents(f1(0.67, cw * G, 0.5) / f1(0.87, (cw + body) * G, 0.5)):.0f} cents")
    v = octa_inner_litres(130)
    print(f"body inner volume {v:.2f} L -> 1 neck {helmholtz(v):.0f} Hz, 4 necks {helmholtz(v, n=4):.0f} Hz")

    print("\n== C Tonos Tensegrity")
    b = 150.0
    L = math.sqrt(6) * b / 1000
    print(f"struts {4 * b:.0f} mm, strings {L * 1000:.0f} mm, strut force = sqrt(6) x string tension")
    for d in (0.6, 0.5, 0.4):
        print(f"  steel {d} mm at 50 N: f1 {f1(L, 50, d):.0f} Hz")
    print("  ratio 0.6:0.5:0.4 mm -> 10:12:15, a just minor triad")
    print(f"  strut compression at 50 N strings: {math.sqrt(6) * 50:.0f} N")
    for od, idd in ((30, 26), (40, 36)):
        i = math.pi / 64 * (od ** 4 - idd ** 4)
        print(f"  Euler buckling, {od} mm PETG tube 600 mm pinned: {math.pi ** 2 * 2000 * i / 600 ** 2:.0f} N")
    v = octa_inner_litres(90)
    print(f"centre chamber {v:.2f} L -> 1 neck {helmholtz(v):.0f} Hz")

    print("\n== D Star of Valihas")
    tube_f = C_SOUND / (4 * (0.300 + 0.61 * 0.013))
    print(f"spoke tube 300 mm stopped: {tube_f:.0f} Hz")
    for d in (0.3, 0.4):
        t = (2 * 0.27 * tube_f) ** 2 * mu(d)
        print(f"  steel {d} mm x 0.27 m tuned to the tube: {t:.1f} N")
    v = octa_inner_litres(110)
    print(f"core inner volume {v:.2f} L -> 1 neck {helmholtz(v):.0f} Hz")

    print("\n== bob and wind")
    shot = 4.6                                   # kg/L bulk steel shot
    vol = cw / shot
    r_oct = (3 * vol * 1e6 / 4) ** (1 / 3)
    print(f"10 kg of steel shot = {vol:.1f} L -> octahedral bob vertex radius {r_oct:.0f} mm")
    for name, area, mass in (("A", 0.09, 11.2), ("B", 0.12, 11.3), ("C", 0.15, 11.7), ("D", 0.16, 11.5)):
        tilts = [math.degrees(math.atan(0.5 * 1.2 * u * u * area / (mass * G))) for u in (10, 20)]
        print(f"  {name}: tilt {tilts[0]:.1f} deg at 10 m/s, {tilts[1]:.1f} deg at 20 m/s")
    print(f"  pendulum period, 1 m hook-to-mass: {2 * math.pi * math.sqrt(1 / G):.1f} s")


if __name__ == "__main__":
    numbers()
    render_all()
