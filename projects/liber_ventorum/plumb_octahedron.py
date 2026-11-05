#!/usr/bin/env python3
"""Concept A developed, rev 1: the Plumb Octahedron.

Design rule for rev 1: every string that can take its tension from hanging mass does, so heat and PETG
creep cannot detune it. The four waist strings are left temperature-sensitive on purpose.

    python plumb_octahedron.py

Prints the rev-1 numbers quoted in the research doc and writes plumb_octahedron_rev1.png.
"""
import math

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

from concepts import (BG, BRONZE, CHAMBER, CLEAR, CORD, OCTA, STEEL, ball, font, octahedron, string,
                      strouhal, tube)

pv.OFF_SCREEN = True

G, C_SOUND, NU_AIR = 9.81, 343.0, 1.5e-5
RHO_STEEL, BREAK_MPA = 7850.0, 2600.0          # plain steel guitar string, approximate
RHO_PETG, CTE_PETG, CTE_STEEL = 1.27e-3, 65e-6, 12e-6   # g/mm^3, 1/K, 1/K
SHOT_KG_PER_L = 4.6
COS45 = math.cos(math.radians(45))

EDGE = 609.5                                   # mm, saddle to saddle (24 in)
RV = EDGE / math.sqrt(2)                       # vertex radius, 431 mm
HUB_R, TIP_R, END_R = 125.0, 55.0, 40.0        # octahedral hub, tip nodes, top/bottom nodes
TUBE_OD, TUBE_ID = 30.0, 26.0
TUBE_LEN = (RV - TIP_R) - HUB_R                # hub vertex to tip-node vertex

GAUGES_IN = [0.018, 0.015, 0.012, 0.010]       # 6:5:4:3.33 -> 1 : 6/5 : 3/2 : 9/5
EARTH = ["D3", "F3", "A3", "C4"]
HEAVEN = ["E3", "G3", "B3", "D4"]
ROOT_HZ = 146.832                              # D3
WHOLE_TONE = 9 / 8
WAIST_OFFSET_HZ = 0.8                          # waist strings sit this far above their earth partner
SPRING_N_PER_MM = 2.0

# frame mass estimates, kg
M_HUB, M_TIPNODE, M_TOP, M_BOTTOM = 0.45, 0.12, 0.06, 0.08
M_TUBE = math.pi / 4 * (TUBE_OD ** 2 - TUBE_ID ** 2) * TUBE_LEN * RHO_PETG / 1000


def d_mm(g_in):
    return g_in * 25.4


def mu(g_in):
    return RHO_STEEL * math.pi * (d_mm(g_in) / 1000) ** 2 / 4


def freq(T, g_in, L=EDGE / 1000):
    return math.sqrt(T / mu(g_in)) / (2 * L)


def harmonic_at(U, g_in, f1):
    d = d_mm(g_in) / 1000
    re = U * d / NU_AIR
    if re < 47:
        return 0
    return max(1, round(strouhal(re) * U / d / f1))


def numbers():
    out = {}
    L = EDGE / 1000
    t_low = (2 * L * ROOT_HZ) ** 2 * mu(GAUGES_IN[0])
    m_bob = 4 * COS45 * t_low / G - M_BOTTOM
    out.update(t_low=t_low, m_bob=m_bob)
    print(f"vertex radius {RV:.1f} mm; tube {TUBE_LEN:.0f} mm; tube mass {M_TUBE * 1000:.0f} g")
    print(f"earth tension {t_low:.2f} N -> bob {m_bob:.2f} kg "
          f"({m_bob / SHOT_KG_PER_L:.2f} L of shot, octahedral vertex radius "
          f"{(3 * m_bob / SHOT_KG_PER_L * 1e6 / 4) ** (1 / 3):.0f} mm)")

    t_up = t_low * WHOLE_TONE ** 2
    m_eff_needed = (t_up - t_low) * COS45 / G
    m_frame = M_TIPNODE + M_TUBE + M_HUB / 4
    trim = m_eff_needed - m_frame
    out.update(t_up=t_up, trim=trim, m_frame=m_frame)
    no_trim_cents = 1200 * math.log2(math.sqrt((t_low + m_frame * G / COS45) / t_low))
    print(f"heaven tension {t_up:.2f} N; load per tip needed {m_eff_needed:.3f} kg; frame gives "
          f"{m_frame:.3f} kg (alone: +{no_trim_cents:.0f} cents); trim per tip {trim:.3f} kg")

    print("\nstring  earth            heaven           break-load share (heaven)")
    for g_in, e, h in zip(GAUGES_IN, EARTH, HEAVEN):
        area = math.pi * d_mm(g_in) ** 2 / 4
        brk = BREAK_MPA * area
        print(f"  .{round(g_in * 1000):03d}  {e} {freq(t_low, g_in):6.1f} Hz   {h} {freq(t_up, g_in):6.1f} Hz"
              f"   {t_up / brk * 100:4.0f}% of {brk:.0f} N")

    total = m_bob + M_BOTTOM + 4 * (M_TIPNODE + M_TUBE + trim) + M_HUB + M_TOP
    out["total"] = total
    print(f"\ntotal hanging mass {total:.1f} kg")
    for u in (10, 20):
        f = 0.5 * 1.2 * u * u * 1.0 * 0.09
        print(f"  tilt at {u} m/s: {math.degrees(math.atan(f / (total * G))):.1f} deg")
    bob_drag = 0.5 * 1.2 * 10 ** 2 * 0.8 * 0.05
    dt = bob_drag / (2 * math.sin(math.radians(45)))
    print(f"  bob drag at 10 m/s {bob_drag:.1f} N -> lower tensions +/-{dt:.1f} N "
          f"(+/-{1200 * math.log2(math.sqrt((t_low + dt) / t_low)):.0f} cents)")
    per_9c = m_bob * ((2 ** (9 / 1200)) ** 2 - 1)
    print(f"  shot to move the earth chord 9 cents: {per_9c * 1000:.0f} g")

    print("\nwind -> harmonic (earth strings), onset wind per gauge")
    for g_in in GAUGES_IN:
        f1 = freq(t_low, g_in)
        onset = 47 * NU_AIR / (d_mm(g_in) / 1000)
        row = [harmonic_at(u, g_in, f1) for u in (1.5, 2, 3, 4, 5)]
        print(f"  .{round(g_in * 1000):03d} onset {onset:.2f} m/s; harmonics at 1.5/2/3/4/5 m/s: {row}")

    print(f"\ntheta: waist offset {WAIST_OFFSET_HZ} Hz -> beat {WAIST_OFFSET_HZ} x n; 4-8 Hz needs n = "
          f"{4 / WAIST_OFFSET_HZ:.0f}-{8 / WAIST_OFFSET_HZ:.0f}")
    for u in (1.5, 2, 2.5, 3, 3.5, 4, 5):
        n = harmonic_at(u, GAUGES_IN[0], freq(t_low, GAUGES_IN[0]))
        print(f"  D string at {u} m/s: harmonic {n}, beat {n * WAIST_OFFSET_HZ:.1f} Hz")

    print("\ntemperature")
    grow = CTE_PETG * RV - CTE_STEEL * EDGE / math.sqrt(2)     # mm/K, tip radius vs steel
    span = grow * math.sqrt(2)
    for g_in in GAUGES_IN:
        ea = 200e3 * math.pi * d_mm(g_in) ** 2 / 4 / EDGE        # N/mm, string axial stiffness
        k = 1 / (1 / SPRING_N_PER_MM + 1 / ea)
        f0 = freq(t_low, g_in) + WAIST_OFFSET_HZ
        dT = k * span
        print(f"  waist .{round(g_in * 1000):03d}: {dT:.3f} N/K -> {f0 * dT / (2 * t_low):.3f} Hz/K "
              f"(no spring: {f0 * ea * span / (2 * t_low):.2f} Hz/K)")
    ang = math.degrees(math.asin(min(1, (RV + 20 * CTE_PETG * RV) / EDGE)))
    print(f"  gravity strings: 20 K moves the string angle 45.00 -> {ang:.3f} deg, "
          f"{1200 * math.log2(math.sqrt(COS45 / math.cos(math.radians(ang)))):.2f} cents")

    print("\nstructure")
    t_w = t_low * (1 + WAIST_OFFSET_HZ / ROOT_HZ) ** 2
    comp = COS45 * (t_up + t_low) + 2 * COS45 * t_w
    area = math.pi / 4 * (TUBE_OD ** 2 - TUBE_ID ** 2)
    inertia = math.pi / 64 * (TUBE_OD ** 4 - TUBE_ID ** 4)
    pcr = math.pi ** 2 * 2000 * inertia / TUBE_LEN ** 2
    print(f"  tube compression {comp:.0f} N, stress {comp / area:.2f} MPa, Euler {pcr:.0f} N "
          f"(x{pcr / comp:.1f})")
    print(f"  top node carries {total * G:.0f} N; each heaven string {t_up:.1f} N")

    v_in = 4 / 3 * (HUB_R - 3 * math.sqrt(3)) ** 3 / 1e6
    a = math.pi * 0.015 ** 2
    leff = a * C_SOUND ** 2 / ((2 * math.pi * ROOT_HZ) ** 2 * v_in / 1000)
    print(f"  hub inner volume {v_in:.2f} L; one 30 mm neck for D3 needs L_eff {leff * 1000:.1f} mm "
          f"-> neck length {leff * 1000 - 1.5 * 15:.1f} mm (lumped)")
    return out


# ---------------------------------------------------------------- render
NODE_OCTA = dict(color="#e8eef0", opacity=0.82, specular=0.6, show_edges=True,
                 edge_color="#8a969b", line_width=1)


def model():
    parts, anchors = [], {}
    top, bot = np.array([0, 0, RV]), np.array([0, 0, -RV])
    tips = [np.array(p, float) * RV for p in [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]]
    parts.append((octahedron(r=HUB_R), OCTA))
    for sx, sy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
        n = np.array([sx, sy, -1.0]) / math.sqrt(3)
        c0 = np.array([sx, sy, -1.0]) * HUB_R / 3
        parts.append((tube(c0, c0 + n * 22, 15), CHAMBER))
    for p in tips:
        u = p / RV
        parts.append((tube(u * HUB_R, u * (RV - TIP_R), TUBE_OD / 2), CLEAR))
        parts.append((octahedron(p, TIP_R), NODE_OCTA))
    parts.append((octahedron(top, END_R), NODE_OCTA))
    parts.append((octahedron(bot, END_R), NODE_OCTA))
    for i, p in enumerate(tips):
        q = tips[(i + 1) % 4]
        parts += [(string(top, p, 1.2), STEEL), (string(p, q, 1.2), STEEL), (string(p, bot, 1.2), STEEL)]
    z_bob_top = bot[2] - END_R - 45
    parts.append((string(bot - [0, 0, END_R], (0, 0, z_bob_top), 2.0), CORD))
    bob_r = 125.0
    parts.append((octahedron((0, 0, z_bob_top - bob_r), bob_r), BRONZE))
    z_ring = top[2] + END_R + 75
    parts.append((string(top + [0, 0, END_R], (0, 0, z_ring - 26), 2.0), CORD))
    parts.append((tube((0, 0, top[2] + END_R + 18), (0, 0, top[2] + END_R + 34), 9), STEEL))
    ring = pv.ParametricTorus(ringradius=26, crosssectionradius=4.5)
    parts.append((ring.rotate_x(90, inplace=False).translate((0, 0, z_ring), inplace=False), STEEL))

    t1, t2, t3, t4 = tips
    anchors = {
        "hook": ((0, 0, z_ring + 10), "Hook and ball-bearing swivel"),
        "top": (top, "Top node: the four heaven strings meet"),
        "heaven": ((top + t1) / 2, "Heaven strings (upper edges): Em7, tensioned by bob + tip trim"),
        "tip": (t1, "Tip node: trim weight, spring and peg"),
        "waist": ((t1 + t4) / 2, "Waist strings on springs: slow beats that follow the sun"),
        "tube": ((t3 * 0.55), "Clear PETG tube, 30 x 26 x 251 mm"),
        "hub": (np.array([0, 0, -HUB_R * 0.45]), "Hub: Helmholtz chamber tuned to D3, necks below"),
        "earth": ((t3 + bot) / 2, "Earth strings (lower edges): Dm7, tensioned by the bob alone"),
        "bottom": (bot, "Bottom node and shackle"),
        "bob": (np.array([0, 0, z_bob_top - bob_r]), "Bob: 11.8 kg of steel shot in a PETG shell"),
    }
    return parts, anchors


def render(out="plumb_octahedron_rev1.png"):
    parts, anchors = model()
    w, h = 1320, 1900
    b = np.array([m.bounds for m, _ in parts])
    zlo, zhi = b[:, 4].min(), b[:, 5].max()
    focal = np.array([0.0, 0.0, (zlo + zhi) / 2])
    half = 0.53 * (zhi - zlo)
    az, el = math.radians(32), math.radians(14)
    d = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    pl = pv.Plotter(off_screen=True, window_size=(w, h))
    pl.set_background(BG)
    pl.enable_depth_peeling(number_of_peels=10)
    for mesh, style in parts:
        pl.add_mesh(mesh, **style)
    pl.enable_parallel_projection()
    pl.camera.focal_point = focal
    pl.camera.position = focal + d * 6000
    pl.camera.up = (0, 0, 1)
    pl.camera.parallel_scale = half
    pl.enable_anti_aliasing("ssaa")
    pl.screenshot("_rev1_raw.png")
    pl.close()

    fwd = -d
    right = np.cross(fwd, [0, 0, 1.0])
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    s = (h / 2) / half
    margin = 430
    canvas = Image.new("RGB", (w + 2 * margin, h), BG)
    canvas.paste(Image.open("_rev1_raw.png").convert("RGB"), (margin, 0))
    draw = ImageDraw.Draw(canvas)
    fnt = font(25)
    pts = []
    for key, (p, text) in anchors.items():
        q = np.asarray(p, float) - focal
        x, y = margin + w / 2 + (q @ right) * s, h / 2 - (q @ up) * s
        pts.append([x, y, text])
    left = sorted([p for p in pts if p[0] < margin + w / 2], key=lambda p: p[1])
    rgt = sorted([p for p in pts if p[0] >= margin + w / 2], key=lambda p: p[1])
    for side, group in (("l", left), ("r", rgt)):
        last = -1e9
        for x, y, text in group:
            ty = max(y, last + 92)
            last = ty
            words, lines, cur = text.split(), [], ""
            for wd in words:
                trial = (cur + " " + wd).strip()
                if draw.textlength(trial, font=fnt) > margin - 60:
                    lines.append(cur)
                    cur = wd
                else:
                    cur = trial
            lines.append(cur)
            tx = margin - 20 if side == "l" else margin + w + 20
            draw.line([(x, y), (tx, ty)], fill="#6f6f6f", width=2)
            draw.ellipse([x - 5, y - 5, x + 5, y + 5], fill="#2a2a2a")
            for i, line in enumerate(lines):
                ly = ty - 14 + i * 31
                if side == "l":
                    draw.text((tx - 8 - draw.textlength(line, font=fnt), ly), line, fill="#222222", font=fnt)
                else:
                    draw.text((tx + 8, ly), line, fill="#222222", font=fnt)
    canvas.save(out)
    print("wrote", out)


if __name__ == "__main__":
    numbers()
    render()
