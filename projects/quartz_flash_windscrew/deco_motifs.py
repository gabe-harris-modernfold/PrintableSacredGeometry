#!/usr/bin/env python3
"""Art Deco motif library (Paris Exposition 1925, the look of a 1926 instrument).

2D shapely geometry in mm, in a local (u, v) plane; callers place and extrude it.
Every stroke is >= 0.9 mm so it survives a 0.4 mm nozzle (two extrusion lines).

Typeface: Century Gothic Bold -- a geometric sans of the Futura (1927) /
Twentieth Century family.  Broadway, the archetypal Deco face, was tested and
rejected: its hairlines lose ~12 % of the glyph area to a 0.8 mm minimum stroke
even at 14 mm caps, i.e. they would not print.
"""
import math, os
import numpy as np
from shapely.geometry import Polygon, LineString, Point, box as sbox
from shapely.ops import unary_union, transform
from shapely import affinity

FONT_FILE = r"C:\Windows\Fonts\GOTHICB.TTF"
TRACKING = 0.26             # wide Deco letter-spacing, fraction of cap height
MIN_STROKE = 0.9

# ------------------------------------------------------------------ primitives
def ring(r0, r1, c=(0.0, 0.0), res=128):
    p = Point(c)
    return p.buffer(r1, res).difference(p.buffer(r0, res))

def octagon(r_flat, c=(0.0, 0.0)):
    rv = r_flat / math.cos(math.pi / 8)
    return Polygon([(c[0] + rv * math.cos(math.pi / 8 + k * math.pi / 4),
                     c[1] + rv * math.sin(math.pi / 8 + k * math.pi / 4)) for k in range(8)])

def clipped_rect(x0, y0, x1, y1, ch):
    """Rectangle with 45-degree clipped corners -- the Deco plan form."""
    return Polygon([(x0 + ch, y0), (x1 - ch, y0), (x1, y0 + ch), (x1, y1 - ch),
                    (x1 - ch, y1), (x0 + ch, y1), (x0, y1 - ch), (x0, y0 + ch)])

def lozenge(c, w, h):
    return Polygon([(c[0] - w / 2, c[1]), (c[0], c[1] + h / 2), (c[0] + w / 2, c[1]), (c[0], c[1] - h / 2)])

def stroke(pts, w=1.2):
    return LineString(pts).buffer(w / 2, cap_style=2, join_style=2)

def ray(c, a, r0, r1, w0, w1):
    """Tapered ray at math angle a (radians, CCW from +u)."""
    d = np.array([math.cos(a), math.sin(a)])
    n = np.array([-d[1], d[0]])
    c = np.asarray(c, float)
    return Polygon([c + d * r0 + n * w0 / 2, c + d * r1 + n * w1 / 2,
                    c + d * r1 - n * w1 / 2, c + d * r0 - n * w0 / 2])

def sunburst(c, r0, r_long, r_short, n, w0=2.0, w1=MIN_STROKE, a0=0.0, a1=2 * math.pi, closed=True):
    """Alternating long/short tapered rays -- the signature motif of 1925."""
    k = np.arange(n)
    ang = a0 + (a1 - a0) * (k / n if closed else k / max(n - 1, 1))
    return unary_union([ray(c, a, r0, r_long if i % 2 == 0 else r_short, w0, w1)
                        for i, a in enumerate(ang)])

def sector(c, r, a0, a1, res=96):
    a = np.linspace(a0, a1, res)
    return Polygon([tuple(c)] + [(c[0] + r * math.cos(t), c[1] + r * math.sin(t)) for t in a])

def chevron(c, w, h, t=1.2, up=True):
    s = 1 if up else -1
    return stroke([(c[0] - w / 2, c[1] - s * h / 2), (c[0], c[1] + s * h / 2),
                   (c[0] + w / 2, c[1] - s * h / 2)], t)

# ------------------------------------------------------------------ disc symbols (v = outward/up)
def sun(c=(0, 0), r=2.6, n=8, r_ray=(3.7, 5.6)):
    return unary_union([Point(c).buffer(r, 48),
                        sunburst(c, r_ray[0], r_ray[1], r_ray[1] - 1.0, n, 1.6, MIN_STROKE, math.pi / 8)])

def crescent(c=(0, 0), r=5.0):
    outer = Point(c).buffer(r, 64)
    bite = Point((c[0] + 0.42 * r, c[1] + 0.18 * r)).buffer(0.84 * r, 64)
    return outer.difference(bite)

def sunrise(c=(0, 0), r=3.6):
    half = sector((c[0], c[1] - 1.2), r, 0, math.pi)
    rays = sunburst((c[0], c[1] - 1.2), r + 1.2, r + 3.0, r + 2.2, 5, 1.4, MIN_STROKE,
                    math.radians(20), math.radians(160), closed=False)
    horizon = sbox(c[0] - r - 3.4, c[1] - 2.5, c[0] + r + 3.4, c[1] - 1.5)
    return unary_union([half, rays, horizon])

def eye(c=(0, 0), w=11.0, h=5.2):
    """Egyptian-revival eye (the Tutankhamun craze of the 1920s): vesica + iris."""
    d = (w * w / 4 - h * h / 4) / h
    R = d + h / 2
    ves = Point((c[0], c[1] - d)).buffer(R, 96).intersection(Point((c[0], c[1] + d)).buffer(R, 96))
    inner = affinity.scale(ves, (w - 2 * 1.0) / w, (h - 2 * 1.0) / h, origin=c)
    return unary_union([ves.difference(inner), Point(c).buffer(1.7, 32)])

def chevrons(c=(0, 0), up=True, w=7.0, h=3.0, gap=2.6):
    s = 1 if up else -1
    return unary_union([chevron((c[0], c[1] + s * gap / 2), w, h, 1.2, up),
                        chevron((c[0], c[1] - s * gap / 2), w, h, 1.2, up)])

# ------------------------------------------------------------------ type
def glyph(ch, cap):
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(fname=FONT_FILE) if os.path.exists(FONT_FILE) else \
        FontProperties(family="DejaVu Sans", weight="bold")
    hcap = TextPath((0, 0), "H", size=100, prop=fp).get_extents().height
    tp = TextPath((0, 0), ch, size=100 * cap / hcap, prop=fp)
    rings = sorted([r for r in tp.to_polygons() if len(r) >= 3], key=lambda r: -abs(Polygon(r).area))
    shape = Polygon()
    for r in rings:
        p = Polygon(r).buffer(0)
        if p.is_empty:
            continue
        shape = shape.difference(p) if shape.contains(p.representative_point()) else shape.union(p)
    ext = tp.get_extents()
    return shape, ext.x0, ext.x1

def _layout(text, cap):
    gl, wd = [], []
    for ch in text:
        if ch == " ":
            gl.append(None); wd.append(0.42 * cap)
            continue
        g, x0, x1 = glyph(ch, cap)
        gl.append((g, x0)); wd.append(x1 - x0)
    gap = TRACKING * cap
    return gl, wd, sum(wd) + gap * (len(text) - 1), gap

def line_text(text, cap, c=(0, 0)):
    gl, wd, total, gap = _layout(text, cap)
    out, s = [], c[0] - total / 2
    for g, w in zip(gl, wd):
        if g is not None:
            out.append(affinity.translate(g[0], s - g[1], c[1] - cap / 2))
        s += w + gap
    return unary_union(out)

def place(geom, r, alpha):
    """Map local (u, v) -- v pointing outward -- to the point at radius r, angle alpha
    (radians clockwise from 12 o'clock, viewer frame: X right, Y up)."""
    up = np.array([math.sin(alpha), math.cos(alpha)])
    right = np.array([math.cos(alpha), -math.sin(alpha)])
    def tf(x, y, z=None):
        x, y = np.asarray(x), np.asarray(y)
        return (right[0] * x + up[0] * (r + y), right[1] * x + up[1] * (r + y))
    return transform(tf, geom)

def arc_text(text, r_base, alpha_c, cap, max_arc=None):
    """Curved text, letter tops outward, reading clockwise, centred on alpha_c."""
    while True:
        gl, wd, total, gap = _layout(text, cap)
        if max_arc is None or total <= max_arc:
            break
        cap *= 0.95
    out, s = [], -total / 2
    for g, w in zip(gl, wd):
        if g is not None:
            local = affinity.translate(g[0], -(g[1] + w / 2), 0)
            out.append(place(local, r_base, alpha_c + (s + w / 2) / r_base))
        s += w + gap
    return out

def as_polys(geom):
    if geom.is_empty:
        return []
    return [geom] if geom.geom_type == "Polygon" else [g for g in geom.geoms if g.geom_type == "Polygon"]
