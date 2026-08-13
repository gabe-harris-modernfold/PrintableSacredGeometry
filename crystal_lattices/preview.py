#!/usr/bin/env python3
"""Preview renders for the crystal-lattice set.

previews/bravais_frames.png   the 14 lattices as ball-and-stick cells -- the picture the
                              kit assembles into, and the quickest way to see that the
                              centring points and their spiders are where they should be
previews/solid_cells.png      the 14 solid comparison cells with their embossed labels
previews/assembled_frames.png the 14 one-piece fused frames, each shown in the print
                              orientation `build.py` chose for it, standing on the bed
previews/supercells.png       the 14 lattices as repeating blocks of cells, with the unit
                              cell at the origin drawn heavier
previews/wigner_seitz.png     cube, truncated octahedron, rhombic dodecahedron

Matplotlib only, Agg backend, so this runs headless.

Three things had to be worked around, all of them matplotlib's 3D axes rather than the
geometry -- see `draw_parts`, `sharp_edges` and the two `best_view_*` functions.
"""

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402
import numpy as np                                                     # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection  # noqa: E402
from scipy.spatial import ConvexHull                                   # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from build import (KIT_LONGEST, MONO_LONGEST, SOLID_LONGEST, SUPER_N,   # noqa: E402
                   SUPER_LONGEST, WS_A, build_assembled, build_solid,
                   build_supercell, build_ws, wigner_seitz)
from lattices import LATTICES, frame                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "previews")

CORNER = "#3b6fb6"
CENTRE = "#d64545"
STRUT = "#6b7684"
FACE = "#c2ccd9"
LABEL = "#d08a1e"
EDGE = "#39424e"
LIGHT = np.array([0.35, -0.55, 0.76])


# ---------------------------------------------------------------- drawing

def shade(base, normals):
    """Cheap lambert shading so the facets read as a solid rather than a silhouette."""
    c = np.clip(np.asarray(normals) @ LIGHT, 0.0, 1.0) * 0.55 + 0.45
    rgb = np.array(matplotlib.colors.to_rgb(base))
    return np.clip(c[:, None] * rgb[None, :], 0.0, 1.0)


def _normals(tri):
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    return np.divide(n, ln, out=np.zeros_like(n), where=ln > 0)


def draw_parts(ax, parts):
    """Draw several (verts, faces, colour) groups as ONE Poly3DCollection.

    Matplotlib depth-sorts collections as wholes, not triangle by triangle across
    collections, so a cell added as one collection and its embossed label as another draws
    the cell on top and the label vanishes. Merging them into a single collection puts every
    triangle into the same sort and the raised glyphs show."""
    tris, cols = [], []
    for v, f, base in parts:
        if len(f) == 0:
            continue
        t = np.asarray(v, float)[np.asarray(f)]
        tris.append(t)
        cols.append(shade(base, _normals(t)))
    if not tris:
        return
    ax.add_collection3d(Poly3DCollection(np.concatenate(tris),
                                         facecolors=np.concatenate(cols),
                                         edgecolors="none", linewidths=0))


def sharp_edges(verts, faces, min_angle=12.0):
    """Edges where the two adjacent facets actually turn a corner.

    Flat-shaded triangle soup with no outlines is hard to read -- gypsum's oblique slab
    came out looking like a self-intersecting bowtie. Drawing the real polyhedron edges
    (and not the triangulation diagonals across each flat face) restores the shape."""
    v = np.asarray(verts, float)
    f = np.asarray(faces, np.int64)
    n = _normals(v[f])
    adj = {}
    for fi, (a, b, c) in enumerate(f):
        for p, q in ((a, b), (b, c), (c, a)):
            adj.setdefault((min(p, q), max(p, q)), []).append(fi)
    cos_lim = math.cos(math.radians(min_angle))
    out = []
    for (p, q), fs in adj.items():
        if len(fs) != 2 or float(n[fs[0]] @ n[fs[1]]) < cos_lim:
            out.append((v[p], v[q]))
    return out


def equalise(ax, pts, pad=1.04, radius=None):
    """Centre the view on `pts`. Pass `radius` to force a scale shared across panels.

    Per-panel autoscaling defeats the solid figure: every cell has a 60 mm longest edge by
    construction, but corundum's bounding box is 128 mm and pyrite's is 60, so autoscaling
    draws corundum's 60 mm edges at half the length of pyrite's and the set stops being
    comparable. One radius for all panels makes equal edges look equal."""
    pts = np.asarray(pts, float)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    c = (lo + hi) / 2.0
    r = (float((hi - lo).max()) / 2.0 * pad) if radius is None else radius
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


# ---------------------------------------------------------------- camera

def _project(pts, elev, azim):
    el, az = math.radians(elev), math.radians(azim)
    f = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    r = np.array([-math.sin(az), math.cos(az), 0.0])
    u = np.cross(f, r)
    return np.stack([np.asarray(pts, float) @ r, np.asarray(pts, float) @ u], axis=1)


def _scan(pts, score, elevs, step):
    best, best_score = (18, -58), -np.inf
    for elev in elevs:
        for azim in range(-180, 180, step):
            s = score(_project(pts, elev, azim))
            if s > best_score:
                best, best_score = (elev, azim), s
    return best


def best_view_separation(pts, elevs=(12, 18, 24, 30), step=5):
    """Camera that spreads the projected nodes out, for the ball-and-stick figure.

    A single shared camera does not survive this set: beryl's 120-degree base ends up
    viewed along its short diagonal and collapses to a sliver, and hemimorphite puts two
    corners exactly on top of its body centre."""
    def score(s):
        d = np.linalg.norm(s[:, None, :] - s[None, :, :], axis=2)
        span = d.max()
        return -np.inf if span < 1e-9 else d[~np.eye(len(s), dtype=bool)].min() / span
    return _scan(pts, score, elevs, step)


def best_view_bulk(pts, elevs=(14, 20, 26), step=6):
    """Camera that maximises the projected silhouette area, for the solid figure.

    Here the *shape* is the message, so the node-separation objective is wrong: for
    corundum it picks the view straight down the 3-fold axis, where a rhombohedron
    elongated 2.5 : 1 reads as a small cube. Silhouette area shows the elongation."""
    def score(s):
        try:
            return float(ConvexHull(s).volume)      # 2D hull: .volume is the area
        except Exception:
            return -np.inf
    return _scan(pts, score, elevs, step)


def common_radius(longest, pad=1.05):
    """Half-width of the largest cell bounding box in the set, so every panel can share
    one scale."""
    return max(float(np.ptp(cell_corners(l, longest), axis=0).max())
               for l in LATTICES) / 2.0 * pad


def cell_corners(latt, longest):
    m = latt.vectors(longest=longest)
    return np.array([i * m[0] + j * m[1] + k * m[2]
                     for i in (0, 1) for j in (0, 1) for k in (0, 1)])


# ---------------------------------------------------------------- figures

def frames_figure(path):
    fig = plt.figure(figsize=(19, 6.6))
    fig.patch.set_facecolor("white")
    rad = common_radius(KIT_LONGEST)
    for i, latt in enumerate(LATTICES):
        ax = fig.add_subplot(2, 7, i + 1, projection="3d")
        m = latt.vectors(longest=KIT_LONGEST)
        f = frame(latt)
        pos = np.array([np.asarray(fr) @ m for fr in f.frac])

        ax.add_collection3d(Line3DCollection([(pos[a], pos[b]) for a, b in f.bonds],
                                             colors=STRUT, linewidths=1.5))
        for kind, col, size in (("corner", CORNER, 46), ("centre", CENTRE, 66)):
            sel = [k for k, t in enumerate(f.kind) if t == kind]
            if sel:
                ax.scatter(pos[sel, 0], pos[sel, 1], pos[sel, 2], c=col, s=size,
                           depthshade=False, edgecolors="white", linewidths=0.7)
        equalise(ax, pos, radius=rad)
        elev, azim = best_view_separation(pos)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"{latt.pearson}  {latt.mineral}\n"
                     f"{len(f.frac)} nodes  {len(f.bonds)} struts", fontsize=9, pad=-2)
    fig.suptitle("The 14 Bravais lattices as printable ball-and-stick cells "
                 f"(longest edge {KIT_LONGEST:.0f} mm; red = centring point)",
                 fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=115)
    plt.close(fig)


def solids_figure(path):
    fig = plt.figure(figsize=(19, 6.6))
    fig.patch.set_facecolor("white")
    rad = common_radius(SOLID_LONGEST)
    for i, latt in enumerate(LATTICES):
        ax = fig.add_subplot(2, 7, i + 1, projection="3d")
        mesh, cap = build_solid(latt)
        body = mesh.pick(drop_tags=["label"])
        draw_parts(ax, [(*body, FACE), (*mesh.pick(tags=["label"]), LABEL)])
        ax.add_collection3d(Line3DCollection(sharp_edges(*body), colors=EDGE,
                                             linewidths=0.8))
        equalise(ax, mesh.arrays()[0], radius=rad)
        elev, azim = best_view_bulk(cell_corners(latt, SOLID_LONGEST))
        ax.view_init(elev=elev, azim=azim)
        a, b, c, al, be, ga = latt.cell()
        ax.set_title(f"{latt.pearson}  {latt.mineral}\n"
                     f"1 : {b / a:.2f} : {c / a:.2f}   {al:.0f}/{be:.0f}/{ga:.0f}",
                     fontsize=9, pad=-2)
    fig.suptitle(f"Solid comparison set, longest cell edge {SOLID_LONGEST:.0f} mm "
                 "(axial ratios a : b : c normalised to a, and the interaxial angles)",
                 fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=115)
    plt.close(fig)


def assembled_figure(path):
    """The fused frames as they sit on the bed, so the chosen orientation is checkable:
    the bed is drawn as a grey plane at z = 0 and the camera is low enough to see which
    hubs are feet."""
    fig = plt.figure(figsize=(19, 6.8))
    fig.patch.set_facecolor("white")
    rad = common_radius(MONO_LONGEST) + 12.0
    for i, latt in enumerate(LATTICES):
        ax = fig.add_subplot(2, 7, i + 1, projection="3d")
        mesh, rep = build_assembled(latt)
        draw_parts(ax, [(*mesh.pick(tags=["strut"]), STRUT),
                        (*mesh.pick(tags=["hub"]), CORNER)])
        v = mesh.arrays()[0]
        lo, hi = v.min(axis=0), v.max(axis=0)
        bed = np.array([[[lo[0] - 6, lo[1] - 6, 0.0], [hi[0] + 6, lo[1] - 6, 0.0],
                         [hi[0] + 6, hi[1] + 6, 0.0], [lo[0] - 6, hi[1] + 6, 0.0]]])
        ax.add_collection3d(Poly3DCollection(bed, facecolors="#e6e9ee",
                                             edgecolors="#b8c0cc", linewidths=0.6))
        equalise(ax, v, radius=rad)
        ax.view_init(elev=12, azim=-62)
        need = "NEEDS SUPPORT" if rep["islands"] else "no support"
        ax.set_title(f"{latt.pearson}  {latt.mineral}\n"
                     f"face {rep['face']} on the bed\n"
                     f"{rep['bridges']} bridges, max {rep['span']:.0f} mm\n{need}",
                     fontsize=8, pad=-4)
    fig.suptitle("One-piece fused frames in their chosen print orientation "
                 f"(longest cell edge {MONO_LONGEST:.0f} mm; grey plane is the bed)",
                 fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=115)
    plt.close(fig)


def supercell_figure(path):
    """The repeating blocks. The origin cell is drawn in a second colour on top of the
    heavier geometry, so the repeating unit is identifiable in the render as well as on the
    printed part."""
    fig = plt.figure(figsize=(19, 7.4))
    fig.patch.set_facecolor("white")
    for i, latt in enumerate(LATTICES):
        ax = fig.add_subplot(2, 7, i + 1, projection="3d")
        mesh, rep = build_supercell(latt)
        draw_parts(ax, [(*mesh.pick(tags=["strut"]), STRUT),
                        (*mesh.pick(tags=["hub"]), CORNER),
                        (*mesh.pick(tags=["cell"]), CENTRE)])
        v = mesh.arrays()[0]
        equalise(ax, v)
        ax.view_init(elev=16, azim=-60)
        need = f"{rep['islands']} islands" if rep["islands"] else "no support"
        ax.set_title(f"{latt.pearson}  {latt.mineral}\n"
                     f"{rep['cells']} cells: {rep['nodes']} hubs, {rep['struts']} struts\n"
                     f"{need}", fontsize=8.5, pad=-4)
    fig.suptitle(f"The 14 lattices as {SUPER_N}x{SUPER_N}x{SUPER_N} repeating blocks "
                 f"(one cell edge {SUPER_LONGEST:.0f} mm; red = the unit cell that repeats)",
                 fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=115)
    plt.close(fig)


def ws_figure(path):
    fig = plt.figure(figsize=(13.5, 5.2))
    fig.patch.set_facecolor("white")
    names = {"cP": "cube", "cI": "truncated octahedron", "cF": "rhombic dodecahedron"}
    bz = {"cP": "cube", "cI": "rhombic dodecahedron", "cF": "truncated octahedron"}
    for i, (pearson, name, pts, vol, rot) in enumerate(wigner_seitz()):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        mesh, cap = build_ws(pearson, pts, rot)
        body = mesh.pick(drop_tags=["label"])
        draw_parts(ax, [(*body, FACE), (*mesh.pick(tags=["label"]), LABEL)])
        ax.add_collection3d(Line3DCollection(sharp_edges(*body), colors=EDGE,
                                             linewidths=0.8))
        equalise(ax, mesh.arrays()[0])
        ax.view_init(elev=20, azim=-58)
        ax.set_title(f"{pearson} Wigner-Seitz cell -- {names[pearson]}\n"
                     f"volume a$^3$/{round(WS_A ** 3 / vol)}     "
                     f"Brillouin zone: {bz[pearson]}", fontsize=10, pad=-6)
    fig.suptitle(f"Cubic Wigner-Seitz cells at a common lattice constant a = {WS_A:.0f} mm"
                 " -- all three space-fill, and cI and cF swap solids in reciprocal space",
                 fontsize=12, y=0.97)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path, dpi=115)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, fn in (("bravais_frames.png", frames_figure),
                     ("solid_cells.png", solids_figure),
                     ("assembled_frames.png", assembled_figure),
                     ("supercells.png", supercell_figure),
                     ("wigner_seitz.png", ws_figure)):
        p = os.path.join(OUT, name)
        fn(p)
        print(f"wrote {os.path.relpath(p, HERE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
