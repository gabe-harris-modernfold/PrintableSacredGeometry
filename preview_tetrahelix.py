#!/usr/bin/env python3
"""
Renders the ten-cycle tetrahelix scaffold, with the counted edge-strand in green.

The green run is the whole argument made visible: ten beams, k to k+3 to k+6 ...
to k+30, spiralling once around the column. The last two tiles are the payoff --
the same axial view of the closed helix and of the true Boerdijk-Coxeter one.
Closed, the green ring meets itself. Regular, it misses by 5.69 deg, and that gap
is the irrationality Fuller hung the DNA argument on.

Counting single tetrahedra instead of strand steps gives 3.66 turns and no
closure at all, which is why the strand is what gets marked here.

No culling: a centroid-plane cut through a space frame does not open it up, it
shreds it into speckle. The lattice is ~89% air already.

Run:  python preview_tetrahelix.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import render_mesh as RM
import tetrahelix as TH

TW, TH_PX = 380, 900


def layers_for(closed):
    m, verts, edges, strand, cells, feet, steps = TH.build(
        1.0, TH.EDGE, False, closed=closed)
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    struts = (v, f[np.isin(tag, ["strut", "node"])], "frame")
    marked = (v, f[tag == "strand"], "beam", False)     # never culled
    base = (v, f[tag == "pad"], "pad")
    lay = [struts, marked] + ([base] if (tag == "pad").any() else [])
    return lay, struts, marked, v, verts, steps


def main():
    lay, struts, marked, v, verts, steps = layers_for(True)
    rlay, rstruts, rmarked, rv, rverts, _ = layers_for(False)
    mid = (0.0, 0.0, 0.5 * (v[:, 2].min() + v[:, 2].max()))
    rmid = (0.0, 0.0, 0.5 * (rv[:, 2].min() + rv[:, 2].max()))
    node = tuple(verts[len(verts) // 2])

    views = [
        (f"closed ten-cycle: {steps} additions", 20, 6, mid, None, lay),
        ("three-quarter", 55, 24, mid, None, lay),
        ("detail, 3 tetrahedra", 35, 12, node, 46.0, lay),
        ("axial, CLOSED: strand meets exactly", 90, 88, mid, None,
         [struts, marked]),
        ("axial, regular BC: misses by 5.69 deg", 90, 88, rmid, None,
         [rstruts, rmarked]),
    ]
    scratch = os.path.join(os.environ.get("TEMP", "."), "tetrahelix_views")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, target, span, lz) in enumerate(views):
        out = os.path.join(scratch, f"th{i}.png")
        RM.render(lz, az, el, w=TW, h=TH_PX, target=target, span=span, out=out)
        im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TH_PX - 20, TW, TH_PX], fill=(238, 238, 235))
        d.text((8, TH_PX - 15), label, fill=(60, 60, 60))
        tiles.append(im)
    sheet = Image.new("RGB", (len(tiles) * TW, TH_PX), (242, 242, 240))
    for i, im in enumerate(tiles):
        sheet.paste(im, (i * TW, 0))
    sheet.save("tetrahelix_views.png")
    print("wrote tetrahelix_views.png", sheet.size)


if __name__ == "__main__":
    main()
