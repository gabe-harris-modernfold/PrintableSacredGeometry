#!/usr/bin/env python3
"""Web assets for the site page web/indirect-perception.html.

For each Negative Space sculpture writes one compact preview mesh (.nsm) to
site/assets/negative-space/ for the in-page 3D viewer: the ribbon, the rod in place, and
the hidden solid it implies, each indexed and quantised to 16 bits over a shared box.
The ~300k-triangle icosahedron ribbon is decimated for the preview only. The printable
STLs are copied unchanged to site/stl/ for download, and the page's hero image is
converted to WebP.

.nsm layout (little-endian): b"NSM1", uint32 chunk count, float32[3] box min,
float32[3] box max, then per chunk: uint32 vertex count, uint32 triangle count,
uint16[3 n] positions (padded to 4 bytes), uint32[3 m] indices.

Run from this directory: python web_assets.py"""

import math
import os
import shutil
import struct

import numpy as np
import pyvista as pv
import trimesh
from PIL import Image
from scipy.spatial import ConvexHull

import negative_space_cube as ncube
import negative_space_desk as nd
import negative_space_ico as nico
import negative_space_pyramid as npyr

SITE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "site")
ASSETS = os.path.join(SITE, "assets", "negative-space")
STL_DIR = os.path.join(SITE, "stl")
HERO = os.path.join(SITE, "..", "indirect-perception.png")
PREVIEW_TRIS = 80_000          # decimate any ribbon above this for the viewer


def hull(verts):
    v = np.asarray(verts, float)
    return trimesh.Trimesh(v, ConvexHull(v).simplices).convex_hull


def decimate(mesh, target):
    if len(mesh.faces) <= target:
        return mesh
    faces = np.hstack([np.full((len(mesh.faces), 1), 3), mesh.faces]).ravel()
    pd = pv.PolyData(mesh.vertices, faces).decimate(1.0 - target / len(mesh.faces))
    return trimesh.Trimesh(pd.points, pd.faces.reshape(-1, 4)[:, 1:])


def write_nsm(path, meshes):
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    hi = np.max([m.bounds[1] for m in meshes], axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    out = [b"NSM1", struct.pack("<I", len(meshes)), struct.pack("<6f", *lo, *hi)]
    for m in meshes:
        q = np.round((m.vertices - lo) / span * 65535).astype("<u2")
        pos = q.tobytes()
        pos += b"\0" * (-len(pos) % 4)
        out += [struct.pack("<II", len(m.vertices), len(m.faces)), pos,
                m.faces.astype("<u4").tobytes()]
    data = b"".join(out)
    with open(path, "wb") as f:
        f.write(data)
    return len(data)


def kanizsa_svg(path, side=330.0, support=0.44):
    """The classic Kanizsa triangle, drawn exactly: three black discs at the corners of an
    upright triangle, an outlined triangle pointing down behind them, and a white triangle
    painted over both. The white triangle has no edge of its own; its corners bite the discs
    and its sides break the outline. support = disc diameter / side (Shipley & Kellman)."""
    rc = side / math.sqrt(3.0)                    # circumradius of the illusory triangle
    r = support * side / 2.0
    up = [(rc * math.cos(a), -rc * math.sin(a)) for a in np.radians([90, 210, 330])]
    down = [(0.85 * rc * math.cos(a), -0.85 * rc * math.sin(a)) for a in np.radians([270, 30, 150])]
    pts = lambda p: " ".join(f"{x:.2f},{y:.2f}" for x, y in p)
    w, h = 600, 560
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="{-w / 2} {-h / 2 - 48} {w} {h}" width="{w}" height="{h}">
<title>Kanizsa triangle</title>
<rect x="{-w / 2}" y="{-h / 2 - 48}" width="{w}" height="{h}" fill="#ffffff"/>
<polygon points="{pts(down)}" fill="none" stroke="#111111" stroke-width="7" stroke-linejoin="miter"/>
{"".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="#111111"/>' for x, y in up)}
<polygon points="{pts(up)}" fill="#ffffff"/>
</svg>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)


def catmull_rom(pts, n=12):
    """Smooth open curve through pts (closed ends clamped)."""
    p = np.vstack([pts[0], pts, pts[-1]])
    out = []
    for i in range(1, len(p) - 2):
        for t in np.linspace(0.0, 1.0, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * (2 * p[i] + (-p[i - 1] + p[i + 1]) * t
                              + (2 * p[i - 1] - 5 * p[i] + 4 * p[i + 1] - p[i + 2]) * t2
                              + (-p[i - 1] + 3 * p[i] - 3 * p[i + 1] + p[i + 2]) * t3))
    out.append(p[-2])
    return np.array(out)


def rubin_svg(path):
    """Rubin's vase: a black vase whose two sides are the profiles of two faces looking at
    each other. Profile given as (height, half-width of the vase) from rim to foot."""
    prof = np.array([(0, 150), (14, 150), (26, 116), (62, 104), (112, 96), (140, 86),
                     (150, 88), (172, 112), (200, 78), (226, 36), (238, 44), (250, 74),
                     (270, 78), (284, 56), (296, 76), (310, 58), (334, 90), (366, 62),
                     (380, 70), (404, 112), (440, 120), (484, 126), (506, 162),
                     (520, 162)], float)
    c = catmull_rom(prof[:, ::-1])               # (half-width, height)
    right = [(x, y - 260) for x, y in c]
    left = [(-x, y) for x, y in right[::-1]]
    d = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in right + left) + " Z"
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="-300 -280 600 560" width="600" height="560">
<title>Rubin's vase</title>
<rect x="-300" y="-280" width="600" height="560" fill="#ffffff"/>
<path d="{d}" fill="#111111"/>
</svg>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)


def necker_svg(path, s=230.0, dx=110.0, dy=-80.0):
    """Necker's cube: twelve equal edges in oblique projection, no hidden lines, no shading,
    so either square can be the front face."""
    front = [(-s / 2, -s / 2), (s / 2, -s / 2), (s / 2, s / 2), (-s / 2, s / 2)]
    front = [(x - dx / 2, y - dy / 2) for x, y in front]
    back = [(x + dx, y + dy) for x, y in front]
    lines = [(front[i], front[(i + 1) % 4]) for i in range(4)]
    lines += [(back[i], back[(i + 1) % 4]) for i in range(4)]
    lines += [(front[i], back[i]) for i in range(4)]
    body = "".join(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}"/>'
                   for a, b in lines)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="-300 -280 600 560" width="600" height="560">
<title>Necker cube</title>
<rect x="-300" y="-280" width="600" height="560" fill="#ffffff"/>
<g stroke="#111111" stroke-width="6" stroke-linecap="round">{body}</g>
</svg>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)


def mooney(out_twotone, out_reveal):
    """A Dalmatian-style two-tone image made from the round window of the hero picture:
    blurred, then thresholded against its own local mean, so only patches remain."""
    from scipy import ndimage as ndi
    img = Image.open(HERO).convert("RGB")
    w, h = img.size
    cx, cy, r = 0.4917 * w, 0.31 * h, 0.175 * w
    box = tuple(int(round(v)) for v in (cx - r, max(cy - r, 0), cx + r, cy + r))
    crop = img.crop(box).resize((840, round(840 * (box[3] - box[1]) / (box[2] - box[0]))),
                                Image.BICUBIC)
    a = np.asarray(crop.convert("L"), float)
    b = ndi.gaussian_filter(a, 8.0)
    dark = b < ndi.gaussian_filter(a, 80.0) - 2.0
    Image.fromarray(np.where(dark, 17, 255).astype(np.uint8)).convert("1").save(
        out_twotone, optimize=True)
    crop.save(out_reveal, "WEBP", quality=82, method=6)


def main():
    os.makedirs(ASSETS, exist_ok=True)
    kanizsa_svg(os.path.join(ASSETS, "kanizsa-triangle.svg"))
    rubin_svg(os.path.join(ASSETS, "rubin-vase.svg"))
    necker_svg(os.path.join(ASSETS, "necker-cube.svg"))
    mooney(os.path.join(ASSETS, "two-tone.png"), os.path.join(ASSETS, "two-tone-reveal.webp"))
    print("kanizsa-triangle.svg, rubin-vase.svg, necker-cube.svg, two-tone.png")
    rod = nd.to_trimesh(nd.rod_part())
    models = {
        "sphere": ("negative_space_desk.stl", rod, trimesh.creation.icosphere(4, nd.R)),
        "icosahedron": ("negative_space_ico.stl", rod, hull(nico.VERTS)),
        "pyramid": ("negative_space_pyramid.stl",
                    nd.to_trimesh(nd.rod_part(npyr.TUBE_BOTTOM, npyr.TUBE_TOP)),
                    hull(npyr.VERTS)),
        "cube": ("negative_space_cube.stl", rod, hull(ncube.VERTS)),
    }
    for name, (stl, rod_mesh, ghost) in models.items():
        ribbon = decimate(trimesh.load(stl), PREVIEW_TRIS)
        size = write_nsm(os.path.join(ASSETS, f"{name}.nsm"), [ribbon, rod_mesh, ghost])
        print(f"{name}.nsm  {len(ribbon.faces)} ribbon triangles, {size // 1024} KB")

    for stl in ["negative_space_desk.stl", "negative_space_desk_rod.stl",
                "negative_space_ico.stl", "negative_space_pyramid.stl",
                "negative_space_pyramid_rod.stl", "negative_space_cube.stl"]:
        shutil.copy2(stl, STL_DIR)
    print(f"copied 6 STLs to {os.path.normpath(STL_DIR)}")

    for src, name, width in [(HERO, "indirect-perception", 1400)] + [
            (f"negative_space_{k}.png", f"render-{k}", 1800)
            for k in ("desk", "ico", "pyramid", "cube")]:
        img = Image.open(src).convert("RGB")
        img.thumbnail((width, width), Image.LANCZOS)
        out = os.path.join(ASSETS, f"{name}.webp")
        img.save(out, "WEBP", quality=84, method=6)
        print(f"{name}.webp  {img.size}, {os.path.getsize(out) // 1024} KB")


if __name__ == "__main__":
    main()
