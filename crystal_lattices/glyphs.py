#!/usr/bin/env python3
"""An eleven-glyph stroke font, just enough to emboss Pearson symbols.

The set of Pearson symbols for the 14 Bravais lattices needs exactly eleven characters --
lowercase `a m o t h c` for the crystal system and uppercase `P S I F R` for the centring
(R for the rhombohedral lattice hR) -- so a full font is unnecessary. `check_coverage`
enforces that: it is what caught R going missing on the first pass.

Each glyph is a list of open polylines in an em box with
the baseline at y = 0 and the cap height at y = 1; `text_quads` sweeps every segment into
a rectangle of the requested stroke width and returns them ready to extrude.

Strokes are butt-capped rectangles plus an octagonal dot at every polyline vertex, which
gives round joins and round caps without a boolean union -- overlapping polygons extrude
into overlapping prisms, which every slicer unions on import. Extending the segments
instead (a projecting cap) is cheaper but serrates the outside of every curve, since each
cap overshoots the vertex on a convex polyline.

Curves are polygonal; at the 10 mm cap height used for the models a 16-gon ring is
indistinguishable from a circle.
"""

import math

X_HEIGHT = 0.62      # lowercase body height, fraction of cap height
STROKE = 0.16        # default stroke width in em
ADVANCE_GAP = 0.10   # inter-letter space in em


def _ring(cx, cy, r, n=16, start=0.0, sweep=360.0, close=True):
    """Polygonal arc. `close` joins the last point back to the first."""
    pts = []
    steps = n if close else max(2, int(round(n * sweep / 360.0)))
    for k in range(steps + (0 if close else 1)):
        t = math.radians(start + sweep * k / steps)
        pts.append((cx + r * math.cos(t), cy + r * math.sin(t)))
    if close:
        pts.append(pts[0])
    return pts


_XR = X_HEIGHT / 2.0     # lowercase bowl radius


# Each entry: (advance width in em, [polyline, ...])
GLYPHS = {
    "o": (2 * _XR + 0.10, [_ring(_XR + 0.05, _XR, _XR)]),
    "c": (2 * _XR + 0.04, [_ring(_XR + 0.05, _XR, _XR, start=40, sweep=280, close=False)]),
    "a": (2 * _XR + 0.10, [_ring(_XR + 0.05, _XR, _XR),
                           [(2 * _XR + 0.05, 0.0), (2 * _XR + 0.05, X_HEIGHT)]]),
    "m": (0.92, [[(0.05, 0.0), (0.05, X_HEIGHT)],
                 [(0.43, 0.0), (0.43, X_HEIGHT)],
                 [(0.81, 0.0), (0.81, X_HEIGHT)],
                 [(0.05, X_HEIGHT), (0.81, X_HEIGHT)]]),
    "h": (0.62, [[(0.05, 0.0), (0.05, 1.0)],
                 [(0.51, 0.0), (0.51, X_HEIGHT)],
                 [(0.05, X_HEIGHT), (0.51, X_HEIGHT)]]),
    "t": (0.56, [[(0.28, 0.0), (0.28, 0.86)],
                 [(0.05, X_HEIGHT), (0.51, X_HEIGHT)]]),
    "P": (0.70, [[(0.08, 0.0), (0.08, 1.0)],
                 [(0.08, 1.0), (0.44, 1.0), (0.60, 0.90),
                  (0.60, 0.64), (0.44, 0.54), (0.08, 0.54)]]),
    "R": (0.70, [[(0.08, 0.0), (0.08, 1.0)],
                 [(0.08, 1.0), (0.44, 1.0), (0.60, 0.90),
                  (0.60, 0.64), (0.44, 0.54), (0.08, 0.54)],
                 [(0.38, 0.54), (0.62, 0.0)]]),
    "S": (0.68, [[(0.60, 0.86), (0.46, 1.0), (0.19, 1.0), (0.05, 0.86),
                  (0.05, 0.68), (0.20, 0.56), (0.45, 0.44), (0.60, 0.32),
                  (0.60, 0.14), (0.46, 0.0), (0.19, 0.0), (0.05, 0.14)]]),
    # Serifed, so it cannot be read as a lowercase l or a 1.
    "I": (0.46, [[(0.23, 0.0), (0.23, 1.0)],
                 [(0.05, 1.0), (0.41, 1.0)],
                 [(0.05, 0.0), (0.41, 0.0)]]),
    "F": (0.64, [[(0.08, 0.0), (0.08, 1.0)],
                 [(0.08, 1.0), (0.58, 1.0)],
                 [(0.08, 0.56), (0.46, 0.56)]]),
}


def text_width(s, cap=1.0, stroke=STROKE):
    """Width of `s` in the same units as `cap`."""
    if not s:
        return 0.0
    w = sum(GLYPHS[ch][0] for ch in s) + ADVANCE_GAP * (len(s) - 1)
    return (w + stroke) * cap


def text_polys(s, cap=1.0, stroke=STROKE, centre=True):
    """Every stroke of `s` as a convex CCW polygon in 2D, scaled to cap height `cap`.

    Butt-capped rectangles for the segments, octagonal dots at the vertices for the joins.
    With `centre` the run is centred on the origin, which is what the embossing code
    wants -- it places polygons relative to a face centroid.
    """
    polys, pen = [], 0.0
    for ch in s:
        adv, strokes = GLYPHS[ch]
        # One dot per distinct vertex. Several glyphs share a vertex between strokes --
        # P's stem meets its bowl at (0.08, 1.0), a closed ring repeats its first point --
        # and emitting a dot per occurrence puts two *exactly coincident* octagons in the
        # file. Overlapping solids are fine; coincident ones are not: merging duplicate
        # vertices on import welds them into a doubled shell with zero volume, which is
        # what verify.py reported for all 14 labelled cells.
        seen = set()
        for poly in strokes:
            for (x0, y0), (x1, y1) in zip(poly, poly[1:]):
                polys.extend(_segment_quad(pen + x0, y0, pen + x1, y1, stroke))
            for x, y in poly:
                key = (round(x, 9), round(y, 9))
                if key in seen:
                    continue
                seen.add(key)
                polys.append(_dot(pen + x, y, stroke / 2.0))
        pen += adv + ADVANCE_GAP

    dx = -text_width(s, 1.0, stroke) / 2.0 + stroke / 2.0 if centre else 0.0
    dy = -0.5 if centre else 0.0
    return [[((x + dx) * cap, (y + dy) * cap) for x, y in p] for p in polys]


def _segment_quad(x0, y0, x1, y1, w):
    """One stroke segment as a butt-capped rectangle."""
    dx, dy = x1 - x0, y1 - y0
    ln = math.hypot(dx, dy)
    if ln < 1e-12:
        return []
    px, py = -dy / ln * w / 2.0, dx / ln * w / 2.0
    return [[(x0 + px, y0 + py), (x0 - px, y0 - py), (x1 - px, y1 - py), (x1 + px, y1 + py)]]


def _dot(cx, cy, r, n=8):
    """Octagonal join/cap, circumscribing the stroke so it does not undercut the width."""
    rc = r / math.cos(math.pi / n)
    return [(cx + rc * math.cos(2 * math.pi * (k + 0.5) / n),
             cy + rc * math.sin(2 * math.pi * (k + 0.5) / n)) for k in range(n)]


PEARSON_CHARS = set("amothcPSIFR")


def check_coverage(symbols):
    """Every character of every Pearson symbol must have a glyph."""
    need = set("".join(symbols))
    missing = need - set(GLYPHS)
    if missing:
        raise KeyError(f"no glyph for {sorted(missing)}")
    return sorted(need)


if __name__ == "__main__":
    from lattices import LATTICES

    syms = [l.pearson for l in LATTICES]
    print("characters used:", "".join(check_coverage(syms)))
    print(f"widest symbol: {max(text_width(s, 10.0) for s in syms):.2f} mm at 10 mm cap")
    for s in syms:
        print(f"  {s}  width {text_width(s, 10.0):5.2f} mm  "
              f"{len(text_polys(s, 10.0))} polys")
