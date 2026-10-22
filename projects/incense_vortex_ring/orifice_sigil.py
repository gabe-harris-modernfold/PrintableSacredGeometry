#!/usr/bin/env python3
"""Non-circular orifice inserts for the Incense Vortex Ring -- the sigil generator.

ORIFICE_GRIMOIRE.md section 3 asks for a tool that takes any outline, checks the six
constraints a mouth must satisfy for the air to carry it, scales it to a chosen equal-area
diameter, and emits a drop-in insert with the same Ø70 flange / Ø56 body / 12 mm re-entrant
collar as the round inserts in vortex_ring_base.py.  This is that tool.

Differences from the round inserts (orifice_insert in vortex_ring_base.py):
* the inlet is a 45° chamfer about 6 mm deep instead of a bellmouth (a bellmouth needs a
  revolve; a chamfer works on any outline).  It is built as a ruled loft from the outline
  to its 6 mm round-cornered offset, i.e. a rolling-ball chamfer, so sharp corners do not
  spike through the body wall;
* the exit edge is square and flush with the top face, as before;
* the label is a short engraved text on the flange.

Outlines
--------
Built-in names:  triangle square pentagon hexagon ellipse slot twin star5 star6 vesica
  (ellipse and slot take --ar; twin is two Ø20 circles 26 mm apart; star5 is the pentagram
  and star6 the hexagram: both fail constraint 2 by design, grimoire Tier III, and need
  --force; star5 defaults to Ø25-equivalent area so its points clear the Ø46 fit limit).
--points "x,y x,y x,y ..."       a closed polygon, any units (it is rescaled to area)
--svg "M 0 0 L 30 0 ... Z"       a simple SVG path: M L H V C S Q T A Z, absolute or relative.
                                 SVG y points down; it is flipped so the mouth, seen from
                                 above (down the tube), reads as drawn.
--points-file / --svg-file       the same from a file (an .svg file: its first path's d="").

The six constraints (grimoire section 3), as implemented
--------------------------------------------------------
1. One closed, non-self-intersecting outline.  Two or more disjoint outlines are accepted
   only for the built-in twin or with --multi (the air treats them as Twin Circles).
2. At most 7 corners, at most 2 of them reflex.  A vertex is a corner when the outline turns
   by more than CORNER_DEG there; sampled curves turn less per step and do not count.
3. Smallest limb or neck width >= 8 mm.  Measured by morphological opening: the outline is
   eroded by r and dilated back; a neck narrower than 2r splits the eroded shape, a limb
   narrower than 2r is removed outright.  Corner caps are removed too, so a removed piece
   only counts as a limb when it reaches more than 3r beyond the opened shape, i.e. when
   it is at least ~1.5x longer than it is wide (a corner sharper than ~29° also fails).
   The reported width is twice the largest r that survives this test.
4. Area between the Ø25 and Ø35 circles (491-962 mm²), circumscribed circle about the
   insert axis <= Ø46 (leaves a 5 mm wall in the Ø56 body).
5. No interior islands.  An outline inside another is an island; cut spokes into the
   drawing so it becomes one outline (the spokes are corners and count under 2).
6. Corners print sharp enough: informational (0.4 mm nozzle -> ~0.3 mm corner radius).

Usage (run from this directory)
-------------------------------
  python orifice_sigil.py --all                 the five inserts of grimoire section 6
  python orifice_sigil.py triangle              one built-in at its grimoire size
  python orifice_sigil.py square --area-dia 25  rescaled to the Ø25-equivalent area
  python orifice_sigil.py --points "0,0 40,0 20,35" --area-dia 30 --name wedge --label WDG
  python orifice_sigil.py --svg "M0 0 L30 0 L30 20 L10 20 Z" --name key --check-only
  python orifice_sigil.py star5 --force         emit it even though a check fails
Options: --rotate DEG, --no-pymeshlab (skip lib/validate_mesh.py), --check-only.
"""

import argparse
import math
import os
import re
import subprocess
import sys

import numpy as np
import shapely
import trimesh
from shapely.geometry import Point, Polygon as SPolygon
from shapely.ops import unary_union
from shapely.validation import explain_validity

from build123d import Align, Kind, Plane, Polygon, Pos, Text, extrude, loft, offset

from vortex_ring_base import (DISC_R, INS_BODY_R, INS_FLANGE_R, MEM_R, PLATE_T, REENTRANT,
                              STROKE, TUBE_OD, TUBE_WALL, Z_LB, Z_LT, cyl, export)

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "lib")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # Ø, °, ² on a Windows console
C = Align.CENTER

# ---------------------------------------------------------------- parameters

CHAMFER = 6.0           # 45° inlet chamfer depth (replaces BELL_R of the round inserts)
LABEL_Y = 29.0          # label centre radius on the flange top, as in orifice_insert
LABEL_FONT = 5.0
LABEL_DEPTH = 0.6
LABEL_MAX_W = 28.0      # text wider than this is shrunk to fit the flange
CORNER_DEG = 20.0       # a vertex turning more than this is a corner
MIN_WIDTH = 8.0         # constraint 3
D_MIN, D_MAX = 25.0, 35.0           # constraint 4, equal-area diameter window
CIRCUM_MAX = 46.0                   # constraint 4, circumscribed circle
MAX_CORNERS, MAX_REFLEX = 7, 2      # constraint 2
BORE = TUBE_OD - 2 * TUBE_WALL      # 96
STROKE_VOL = math.pi * STROKE * (MEM_R ** 2 + MEM_R * DISC_R + DISC_R ** 2) / 3   # 85.8 mL
CURVE_STEP_DEG = 2.0    # sampling of built-in curves (turn per step, well under CORNER_DEG)

# the five inserts of grimoire section 6
GRIMOIRE_SET = (
    ("triangle", dict(area_dia=25.0)),
    ("square", dict(area_dia=30.0)),
    ("ellipse", dict(area_dia=30.0, ar=2.0)),
    ("slot", dict(area_dia=25.0, ar=4.0)),
    ("twin", dict()),
)

# ---------------------------------------------------------------- outlines


def eq_area(d):
    return math.pi * d * d / 4


def regular_polygon(n, circum_r=10.0):
    """Regular n-gon.  Odd n: a vertex on +Y (triangle point-up, seen from above).  Even n: an
    edge on top, so the square is born a square and the T/2 flip makes it the diamond."""
    a0 = math.pi / 2 + (math.pi / n if n % 2 == 0 else 0.0)
    return [(circum_r * math.cos(a0 + 2 * math.pi * k / n),
             circum_r * math.sin(a0 + 2 * math.pi * k / n)) for k in range(n)]


def arc(cx, cy, r, a0, a1, step_deg=CURVE_STEP_DEG):
    n = max(2, int(math.ceil(abs(a1 - a0) / step_deg)) + 1)
    return [(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a)))
            for a in np.linspace(a0, a1, n)]


def ellipse_pts(a, b):
    n = int(360 / CURVE_STEP_DEG)
    return [(a * math.cos(t), b * math.sin(t)) for t in np.linspace(0, 2 * math.pi, n, endpoint=False)]


def slot_pts(length, width):
    """Stadium: overall length x width, long axis on X."""
    r, c = width / 2, (length - width) / 2
    return arc(c, 0, r, -90, 90)[:-1] + arc(-c, 0, r, 90, 270)[:-1]


def circle_pts(cx, cy, r):
    return arc(cx, cy, r, 0, 360)[:-1]


def star_pts(n=5, circum_r=15.0):
    """Regular {n/2} star-polygon outline (n points, n notches), point on +Y: the pentagram for
    n = 5, the hexagram (two overlaid triangles) for n = 6.  The notch radius is where the
    chords skipping one vertex cross: R·cos(2π/n)/cos(π/n)."""
    inner = circum_r * math.cos(math.pi * 2 / n) / math.cos(math.pi / n)
    pts = []
    for k in range(2 * n):
        r = circum_r if k % 2 == 0 else inner
        a = math.pi / 2 + k * math.pi / n
        pts.append((r * math.cos(a), r * math.sin(a)))
    return pts


def vesica_pts(r=15.0):
    """Vesica piscis: two circles of radius r through each other's centres; cusps on ±Y."""
    half = r / 2        # centres at (±r/2, 0); cusps at (0, ±r·√3/2)
    return (arc(-half, 0, r, -60, 60)[:-1] + arc(half, 0, r, 120, 240)[:-1])


def builtin(name, ar=None):
    """Return (outlines, default equal-area diameter, descriptor) for a built-in name."""
    if name in ("triangle", "square", "pentagon", "hexagon"):
        n = {"triangle": 3, "square": 4, "pentagon": 5, "hexagon": 6}[name]
        return [regular_polygon(n)], (25.0 if name == "triangle" else 30.0), "polygon"
    if name == "ellipse":
        ar = ar or 2.0
        return [ellipse_pts(ar * 10.0, 10.0)], 30.0, "ellipse"
    if name == "slot":
        ar = ar or 4.0
        return [slot_pts(ar * 10.0, 10.0)], 25.0, "slot"
    if name == "twin":
        return [circle_pts(-13.0, 0, 10.0), circle_pts(13.0, 0, 10.0)], None, "twin"
    if name in ("star", "star5"):
        return [star_pts(5)], 25.0, "star"
    if name == "star6":
        return [star_pts(6)], 30.0, "star"       # its points fit the Ø46 limit at Ø30-eq
    if name == "vesica":
        return [vesica_pts()], 30.0, "vesica"
    raise SystemExit(f"unknown outline '{name}'")


BUILTINS = ("triangle", "square", "pentagon", "hexagon", "ellipse", "slot", "twin", "star5", "star6", "vesica", "star")

# ---------------------------------------------------------------- input parsing

_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def parse_points(text):
    nums = [float(v) for v in _NUM.findall(text)]
    if len(nums) < 6 or len(nums) % 2:
        raise SystemExit("--points needs at least three x,y pairs")
    return [[(nums[i], nums[i + 1]) for i in range(0, len(nums), 2)]]


def _bezier(ctrl, n=16):
    ctrl = np.asarray(ctrl, dtype=float)
    k = len(ctrl) - 1
    out = []
    for t in np.linspace(0, 1, n + 1)[1:]:
        p = np.zeros(2)
        for i, c in enumerate(ctrl):
            p += math.comb(k, i) * (1 - t) ** (k - i) * t ** i * c
        out.append((p[0], p[1]))
    return out


def _svg_arc(p0, rx, ry, phi_deg, large, sweep, p1):
    """SVG elliptical arc from p0 to p1, sampled (implementation notes F.6.5)."""
    x0, y0 = p0
    x1, y1 = p1
    if (x0, y0) == (x1, y1):
        return []
    rx, ry = abs(rx), abs(ry)
    if rx == 0 or ry == 0:
        return [p1]
    phi = math.radians(phi_deg)
    cp, sp = math.cos(phi), math.sin(phi)
    dx, dy = (x0 - x1) / 2, (y0 - y1) / 2
    x1p = cp * dx + sp * dy
    y1p = -sp * dx + cp * dy
    lam = (x1p / rx) ** 2 + (y1p / ry) ** 2
    if lam > 1:
        rx, ry = rx * math.sqrt(lam), ry * math.sqrt(lam)
    num = rx * rx * ry * ry - rx * rx * y1p * y1p - ry * ry * x1p * x1p
    den = rx * rx * y1p * y1p + ry * ry * x1p * x1p
    coef = math.sqrt(max(0.0, num / den)) if den else 0.0
    if large == sweep:
        coef = -coef
    cxp, cyp = coef * rx * y1p / ry, -coef * ry * x1p / rx
    cx = cp * cxp - sp * cyp + (x0 + x1) / 2
    cy = sp * cxp + cp * cyp + (y0 + y1) / 2

    def ang(ux, uy, vx, vy):
        a = math.atan2(ux * vy - uy * vx, ux * vx + uy * vy)
        return a
    t1 = ang(1, 0, (x1p - cxp) / rx, (y1p - cyp) / ry)
    dt = ang((x1p - cxp) / rx, (y1p - cyp) / ry, (-x1p - cxp) / rx, (-y1p - cyp) / ry)
    if not sweep and dt > 0:
        dt -= 2 * math.pi
    elif sweep and dt < 0:
        dt += 2 * math.pi
    n = max(2, int(math.ceil(abs(math.degrees(dt)) / CURVE_STEP_DEG)))
    out = []
    for t in np.linspace(t1, t1 + dt, n + 1)[1:]:
        ex, ey = rx * math.cos(t), ry * math.sin(t)
        out.append((cp * ex - sp * ey + cx, sp * ex + cp * ey + cy))
    out[-1] = p1
    return out


def parse_svg_path(d):
    """Minimal SVG path parser: M L H V C S Q T A Z (and relative forms).  Returns a list of
    closed outlines, y flipped so the mouth reads as drawn when seen from above."""
    tokens = re.findall(r"[MmLlHhVvCcSsQqTtAaZz]|" + _NUM.pattern, d)
    outlines, cur, start = [], [], (0.0, 0.0)
    pos, last_c = (0.0, 0.0), None
    i, cmd = 0, None

    def take(n):
        nonlocal i
        vals = [float(v) for v in tokens[i:i + n]]
        if len(vals) < n:
            raise SystemExit(f"SVG path: '{cmd}' needs {n} numbers")
        i += n
        return vals

    def close():
        nonlocal cur
        if len(cur) >= 3:
            outlines.append(cur)
        cur = []

    while i < len(tokens):
        if re.match(r"[A-Za-z]", tokens[i]):
            cmd = tokens[i]
            i += 1
            if cmd in "Zz":
                close()
                pos = start
                last_c = None
                continue
        rel = cmd.islower()
        c = cmd.upper()
        ox, oy = pos if rel else (0.0, 0.0)
        if c == "M":
            x, y = take(2)
            close()
            pos = start = (ox + x, oy + y)
            cur = [pos]
            cmd = "l" if rel else "L"          # subsequent pairs are implicit linetos
            last_c = None
        elif c == "L":
            x, y = take(2)
            pos = (ox + x, oy + y)
            cur.append(pos)
            last_c = None
        elif c == "H":
            (x,) = take(1)
            pos = (ox + x, pos[1])
            cur.append(pos)
            last_c = None
        elif c == "V":
            (y,) = take(1)
            pos = (pos[0], oy + y)
            cur.append(pos)
            last_c = None
        elif c in "CS":
            if c == "C":
                x1, y1, x2, y2, x, y = take(6)
                c1 = (ox + x1, oy + y1)
            else:
                x2, y2, x, y = take(4)
                c1 = (2 * pos[0] - last_c[0], 2 * pos[1] - last_c[1]) if last_c else pos
            c2, end = (ox + x2, oy + y2), (ox + x, oy + y)
            cur += _bezier([pos, c1, c2, end])
            pos, last_c = end, c2
        elif c in "QT":
            if c == "Q":
                x1, y1, x, y = take(4)
                c1 = (ox + x1, oy + y1)
            else:
                x, y = take(2)
                c1 = (2 * pos[0] - last_c[0], 2 * pos[1] - last_c[1]) if last_c else pos
            end = (ox + x, oy + y)
            cur += _bezier([pos, c1, end])
            pos, last_c = end, c1
        elif c == "A":
            rx, ry, rot, large, sweep, x, y = take(7)
            end = (ox + x, oy + y)
            cur += _svg_arc(pos, rx, ry, rot, bool(large), bool(sweep), end)
            pos, last_c = end, None
        else:
            raise SystemExit(f"SVG path: unsupported command '{cmd}'")
    close()
    if not outlines:
        raise SystemExit("SVG path: no closed outline found")
    return [[(x, -y) for x, y in o] for o in outlines]


def load_svg_file(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"\sd\s*=\s*\"([^\"]+)\"", text) or re.search(r"\sd\s*=\s*'([^']+)'", text)
    return parse_svg_path(m.group(1) if m else text)


def clean_outline(pts):
    """Drop repeated / near-duplicate points and a closing point equal to the first."""
    out = []
    for p in pts:
        if not out or math.dist(out[-1], p) > 1e-6:
            out.append(tuple(float(v) for v in p))
    if len(out) > 1 and math.dist(out[0], out[-1]) < 1e-6:
        out.pop()
    return out

# ---------------------------------------------------------------- geometry + checks


def to_shapely(outlines):
    return [SPolygon(o) for o in outlines]


def corners(pts):
    """(corner count, reflex count, sharpest interior angle in degrees) of a closed polyline."""
    p = np.asarray(pts, dtype=float)
    n = len(p)
    a = p - np.roll(p, 1, axis=0)             # incoming edge vectors
    b = np.roll(p, -1, axis=0) - p            # outgoing
    cross = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    dot = (a * b).sum(axis=1)
    turn = np.degrees(np.arctan2(cross, dot))         # signed turning angle at each vertex
    area2 = (p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]).sum()
    orient = 1.0 if area2 > 0 else -1.0
    is_corner = np.abs(turn) > CORNER_DEG
    reflex = is_corner & (np.sign(turn) != orient)
    sharpest = 180.0 - np.abs(turn[is_corner]).max() if is_corner.any() else 180.0
    return int(is_corner.sum()), int(reflex.sum()), float(sharpest), n


def _opening_ok(poly, r):
    """True when eroding by r and dilating back removes nothing but corner caps."""
    er = poly.buffer(-r)
    if er.is_empty or er.geom_type != "Polygon":
        return False                                # vanished, or a neck split it
    opened = er.buffer(r)
    removed = poly.difference(opened)
    parts = list(removed.geoms) if hasattr(removed, "geoms") else [removed]
    reach = 0.0
    for comp in parts:
        if comp.is_empty or comp.area < 0.05:
            continue
        reach = max(reach, max(opened.distance(Point(c)) for c in comp.exterior.coords))
    return reach <= 3 * r                           # beyond 3r it is a limb, not a cap


def limb_width(poly):
    """Twice the largest erosion radius r that passes _opening_ok (module docstring,
    constraint 3): a 0.25 mm scan up to the inscribed radius, then bisected to 0.01 mm."""
    r_max = shapely.maximum_inscribed_circle(poly, 0.01).length
    lo, hi = 0.0, None
    r = 0.25
    while r <= r_max + 1e-9:
        if not _opening_ok(poly, r):
            hi = r
            break
        lo = r
        r += 0.25
    if hi is None:
        hi = r_max + 0.25
    while hi - lo > 0.01:
        mid = (lo + hi) / 2
        if _opening_ok(poly, mid):
            lo = mid
        else:
            hi = mid
    return 2 * lo, 2 * r_max


def analyse(outlines, multi_ok=False):
    """Run the six constraints and the hydraulics.  Returns a dict with 'checks' (list of
    (number, title, passed, detail)) and the metrics."""
    polys = to_shapely(outlines)
    checks = []

    # 1. one closed, non-self-intersecting outline
    bad = [explain_validity(p) for p in polys if not p.is_valid]
    valid = [p for p in polys if p.is_valid]
    disjoint = len(polys) > 1 and all(polys[i].disjoint(polys[j]) for i in range(len(polys))
                                      for j in range(i + 1, len(polys)))
    inside = [(i, j) for i, a in enumerate(valid) for j, b in enumerate(valid) if i != j and b.within(a)]
    if bad:
        checks.append((1, "one closed, non-self-intersecting outline", False, "; ".join(bad)))
    elif inside:
        checks.append((1, "one closed, non-self-intersecting outline", False,
                       f"{len(polys)} outlines, {len(inside)} inside another: an island (see 5)"))
    elif len(polys) == 1:
        checks.append((1, "one closed, non-self-intersecting outline", True, "1 outline, simple"))
    elif disjoint and multi_ok:
        checks.append((1, "one closed, non-self-intersecting outline", True,
                       f"{len(polys)} disjoint outlines, accepted as a deliberate pair (twin)"))
    elif disjoint:
        checks.append((1, "one closed, non-self-intersecting outline", False,
                       f"{len(polys)} separate openings: the air treats them as Twin Circles; pass --multi to accept"))
    else:
        checks.append((1, "one closed, non-self-intersecting outline", False,
                       f"{len(polys)} outlines that touch or overlap each other"))

    # 2. corners
    n_c = n_r = 0
    sharpest = 180.0
    for o in outlines:
        c, r, s, _ = corners(o)
        n_c += c
        n_r += r
        sharpest = min(sharpest, s)
    ok2 = n_c <= MAX_CORNERS and n_r <= MAX_REFLEX
    checks.append((2, f"<= {MAX_CORNERS} corners, <= {MAX_REFLEX} reflex", ok2,
                   f"{n_c} corners (turn > {CORNER_DEG:.0f}°), {n_r} reflex"))

    # 3. limb / neck width, measured on the open region (outer outlines minus any islands)
    region = unary_union([p for k, p in enumerate(valid) if not any(j == k for _, j in inside)])
    for _, j in inside:
        region = region.difference(valid[j])
    widths = [limb_width(p) for p in (region.geoms if hasattr(region, "geoms") else [region]) if not p.is_empty]
    w_min = min(w for w, _ in widths) if widths else 0.0
    insc = min(d for _, d in widths) if widths else 0.0
    checks.append((3, f"smallest limb/neck width >= {MIN_WIDTH:.0f} mm", w_min >= MIN_WIDTH - 1e-6,
                   f"{w_min:.1f} mm (largest inscribed circle Ø{insc:.1f})"))

    # 4. area window and circumscribed circle
    area = sum(p.area for p in polys)
    d_eq = math.sqrt(4 * area / math.pi)
    circum = 2 * max(math.hypot(x, y) for o in outlines for x, y in o)
    ok4 = (eq_area(D_MIN) * (1 - 5e-3) <= area <= eq_area(D_MAX) * (1 + 5e-3)) and circum <= CIRCUM_MAX + 1e-6
    checks.append((4, f"area Ø{D_MIN:.0f}-Ø{D_MAX:.0f} eq ({eq_area(D_MIN):.0f}-{eq_area(D_MAX):.0f} mm²), circumscribed <= Ø{CIRCUM_MAX:.0f}",
                   ok4, f"{area:.0f} mm² = Ø{d_eq:.1f} eq; circumscribed Ø{circum:.1f} about the axis"))

    # 5. islands
    islands = len(inside)
    checks.append((5, "no interior islands (unless spoked)", islands == 0,
                   "none" if islands == 0 else f"{islands}: cut three spokes so it becomes one outline"))

    # 6. corners print sharp enough (informational)
    checks.append((6, "corners print sharp enough", True,
                   f"sharpest corner {sharpest:.0f}°; 0.4 mm nozzle leaves ~0.3 mm radius, below the core radius"
                   if sharpest < 180 else "no corners"))

    perim = sum(p.length for p in polys)
    if area > 1e-9:
        d_h = 4 * area / perim
        ld = (STROKE_VOL / area) / d_h          # slug length / hydraulic diameter
        press = min(100.0, 100.0 * 4.0 / ld)
    else:
        d_h = ld = press = float("nan")
    bb = unary_union(polys).bounds
    return dict(checks=checks, area=area, d_eq=d_eq, perim=perim, d_h=d_h, ld=ld, press=press,
                circum=circum, width=w_min, corners=n_c, reflex=n_r, bbox=(bb[2] - bb[0], bb[3] - bb[1]),
                ring=1.2 * circum, passed=all(c[2] for c in checks))


def describe(kind, outlines, a):
    """Human dimensions for the report / README."""
    w, h = a["bbox"]
    if kind == "polygon":
        return f"side {a['perim'] / a['corners']:.1f}"
    if kind == "ellipse":
        return f"{w:.1f} x {h:.1f}"
    if kind == "slot":
        return f"{w:.1f} x {h:.1f}, rounded ends"
    if kind == "star":
        return f"{a['corners'] // 2}-point star, Ø{a['circum']:.1f} over the points"
    if kind == "twin":
        d = math.sqrt(4 * SPolygon(outlines[0]).area / math.pi)
        c0, c1 = SPolygon(outlines[0]).centroid, SPolygon(outlines[1]).centroid
        return f"2 x Ø{d:.1f}, {c0.distance(c1):.1f} apart"
    return f"{w:.1f} x {h:.1f} bbox"


def report(name, kind, outlines, a):
    print(f"\n=== {name}: {describe(kind, outlines, a)} ===")
    for num, title, ok, detail in a["checks"]:
        print(f"  [{'PASS' if ok else 'FAIL'}] {num}. {title}: {detail}")
    print(f"  area {a['area']:.0f} mm² (Ø{a['d_eq']:.1f} eq)  perimeter {a['perim']:.1f}  "
          f"D_h {a['d_h']:.1f}  full-press slug L/D_h {a['ld']:.1f}  -> press to ~{a['press']:.0f} %")
    print(f"  ring at birth ~{a['ring']:.0f} mm across = {a['ring'] / BORE:.2f} of the {BORE:.0f} mm bore")
    print("  RESULT:", "PASS" if a["passed"] else "FAIL")

# ---------------------------------------------------------------- transform


def centre_scale(outlines, area_dia=None, rotate=0.0):
    polys = to_shapely(outlines)
    area = sum(p.area for p in polys)
    if area > 1e-9:
        cen = unary_union(polys).centroid
    else:                                       # self-intersecting: centre on the bounding box
        b = unary_union([p.envelope for p in polys]).bounds
        cen = Point((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
    k = math.sqrt(eq_area(area_dia) / area) if (area_dia and area > 1e-9) else 1.0
    th = math.radians(rotate)
    out = []
    for o in outlines:
        pts = []
        for x, y in o:
            x, y = (x - cen.x) * k, (y - cen.y) * k
            pts.append((x * math.cos(th) - y * math.sin(th), x * math.sin(th) + y * math.cos(th)))
        out.append(pts)
    return out

# ---------------------------------------------------------------- CAD


def sigil_insert(outlines, label):
    """Ø70 flange / Ø56 body / 12 mm re-entrant collar, as orifice_insert(); the bore is the
    outline(s), with a 45° rolling-ball chamfer at the inlet and a square edge at the exit."""
    z_bot = Z_LB - REENTRANT
    part = cyl(INS_BODY_R, z_bot, Z_LT - 3.0) + cyl(INS_FLANGE_R, Z_LT - 3.0, Z_LT)
    for pts in outlines:
        face = Polygon(*pts, align=None)
        part -= extrude(Plane.XY.offset(z_bot - 1.0) * face, amount=REENTRANT + PLATE_T + 2.0)
        # inlet chamfer: loft from the 7 mm round-cornered offset 1 mm below the inlet face to
        # the outline itself CHAMFER above it (45°); the body's wall clips the rest
        try:
            lo = Plane.XY.offset(z_bot - 1.0) * offset(face, amount=CHAMFER + 1.0, kind=Kind.ARC)
            hi = Plane.XY.offset(z_bot + CHAMFER) * face
            cham = loft([lo, hi], ruled=True)
        except Exception as e:                       # fallback: sharp-cornered taper, clipped
            print(f"  (round-cornered chamfer failed: {e}; using a sharp taper)")
            cham = extrude(Plane.XY.offset(z_bot + CHAMFER) * face, amount=-(CHAMFER + 1.0), taper=-45)
            cham &= cyl(INS_BODY_R - 2.0, z_bot - 2.0, z_bot + CHAMFER + 1.0)
        part -= cham
    if label:
        txt = Text(label, font_size=LABEL_FONT, align=(C, C))
        w = txt.bounding_box().size.X
        if w > LABEL_MAX_W:
            txt = Text(label, font_size=LABEL_FONT * LABEL_MAX_W / w, align=(C, C))
        part -= extrude(Plane.XY.offset(Z_LT) * Pos(0, LABEL_Y) * txt, amount=-LABEL_DEPTH)
    return part


def validate_with_pymeshlab(path):
    r = subprocess.run([sys.executable, os.path.join(LIB, "validate_mesh.py"), path],
                       capture_output=True, text=True)
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT:")), r.stderr.strip()[-200:])
    print(f"  validate_mesh.py: {line}")
    return r.returncode == 0

# ---------------------------------------------------------------- driver


def run_one(name, outlines, kind, label, area_dia, rotate, multi_ok, force, check_only, pymeshlab):
    outlines = [clean_outline(o) for o in outlines]
    outlines = centre_scale(outlines, area_dia, rotate)
    a = analyse(outlines, multi_ok)
    report(name, kind, outlines, a)
    row = dict(name=name, dims=describe(kind, outlines, a), **{k: a[k] for k in ("area", "d_eq", "d_h", "ld", "press", "circum")},
               passed=a["passed"], file=None, mesh_ok=None, pml_ok=None, extents=None)
    if check_only:
        return row
    if not a["passed"] and not force:
        print("  not generated (pass --force to emit it anyway)")
        return row
    fname = f"orifice_{name}_{a['d_eq']:.0f}eq.stl"
    print(f"  building {fname} ...")
    part = sigil_insert(outlines, label)
    ok = export(part, fname, flip=True)
    m = trimesh.load(os.path.join(HERE, fname))
    row.update(file=fname, mesh_ok=ok, extents=tuple(round(float(v), 1) for v in m.extents))
    if pymeshlab:
        row["pml_ok"] = validate_with_pymeshlab(os.path.join(HERE, fname))
    return row


def summary(rows):
    print("\nSummary")
    print(f"  {'insert':22s} {'dimensions (mm)':28s} {'area':>5s} {'D_eq':>5s} {'D_h':>5s} {'L/D_h':>5s} {'press':>5s}  {'circ':>5s}  checks  file")
    for r in rows:
        f = r["file"] or "-"
        mesh = "" if r["mesh_ok"] is None else ("  mesh OK" if r["mesh_ok"] else "  mesh BAD")
        pml = "" if r["pml_ok"] is None else ("  pymeshlab OK" if r["pml_ok"] else "  pymeshlab FAIL")
        ext = "" if r["extents"] is None else f"  {r['extents'][0]:.0f}x{r['extents'][1]:.0f}x{r['extents'][2]:.0f}"
        print(f"  {r['name']:22s} {r['dims']:28s} {r['area']:5.0f} {r['d_eq']:5.1f} {r['d_h']:5.1f} {r['ld']:5.1f} {r['press']:4.0f}%  {r['circum']:5.1f}  "
              f"{'PASS' if r['passed'] else 'FAIL'}    {f}{ext}{mesh}{pml}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("shape", nargs="?", choices=BUILTINS, help="built-in outline")
    ap.add_argument("--all", action="store_true", help="the five inserts of grimoire section 6")
    ap.add_argument("--points", help='closed polygon "x,y x,y x,y ..."')
    ap.add_argument("--points-file", help="file of x,y pairs")
    ap.add_argument("--svg", help="SVG path data (M L H V C S Q T A Z)")
    ap.add_argument("--svg-file", help=".svg file (first path) or a file of path data")
    ap.add_argument("--area-dia", type=float, help="scale to the area of this circle diameter (mm)")
    ap.add_argument("--ar", type=float, help="aspect ratio for ellipse (default 2) / slot (default 4)")
    ap.add_argument("--rotate", type=float, default=0.0, help="rotate the mouth (deg, CCW from above)")
    ap.add_argument("--name", help="output name (orifice_<name>_<D>eq.stl); default from the shape")
    ap.add_argument("--label", help="engraved label (about 8 characters fit); default from the name")
    ap.add_argument("--multi", action="store_true", help="accept several disjoint outlines")
    ap.add_argument("--force", action="store_true", help="emit the STL even if a check fails")
    ap.add_argument("--check-only", action="store_true", help="report only, no CAD")
    ap.add_argument("--no-pymeshlab", action="store_true", help="skip lib/validate_mesh.py")
    args = ap.parse_args()
    pml = not args.no_pymeshlab
    rows = []

    if args.all:
        for shape, kw in GRIMOIRE_SET:
            outlines, _, kind = builtin(shape, kw.get("ar"))
            name = f"{shape}{kw['ar']:.0f}" if "ar" in kw else shape
            label = default_label(shape, kw.get("area_dia"), outlines)
            rows.append(run_one(name, outlines, kind, label, kw.get("area_dia"), 0.0, shape == "twin",
                                False, args.check_only, pml))
    else:
        if args.shape:
            outlines, dia, kind = builtin(args.shape, args.ar)
            area_dia = args.area_dia or dia
            name = args.name or (f"{args.shape}{(args.ar or (2 if args.shape == 'ellipse' else 4)):.0f}"
                                 if args.shape in ("ellipse", "slot") else
                                 "star5" if args.shape == "star" else args.shape)
            label = args.label or default_label(args.shape, area_dia, outlines)
            multi_ok = args.multi or args.shape == "twin"
        else:
            if args.points:
                outlines = parse_points(args.points)
            elif args.points_file:
                outlines = parse_points(open(args.points_file, encoding="utf-8").read())
            elif args.svg:
                outlines = parse_svg_path(args.svg)
            elif args.svg_file:
                outlines = load_svg_file(args.svg_file)
            else:
                ap.error("give a built-in shape, --points/--points-file, --svg/--svg-file, or --all")
            kind = "custom"
            area_dia = args.area_dia or 30.0
            name = args.name or "sigil"
            label = args.label or name[:6].upper()
            multi_ok = args.multi
        rows.append(run_one(name, outlines, kind, label, area_dia, args.rotate, multi_ok,
                            args.force, args.check_only, pml))

    summary(rows)
    ok = all(r["passed"] and r["mesh_ok"] is not False and r["pml_ok"] is not False for r in rows)
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


def default_label(shape, area_dia, outlines):
    short = {"triangle": "TRI", "square": "SQ", "pentagon": "PENT", "hexagon": "HEX", "ellipse": "ELL",
             "slot": "SLOT", "twin": "TWIN", "star": "STAR5", "star5": "STAR5", "star6": "STAR6",
             "vesica": "VES"}[shape]
    if shape == "twin":
        d = math.sqrt(4 * SPolygon(outlines[0]).area / math.pi)
        return f"2x{d:.0f}"
    return f"{short} {area_dia:.0f}" if area_dia else short


if __name__ == "__main__":
    sys.exit(main())
