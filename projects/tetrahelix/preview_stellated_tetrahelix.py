#!/usr/bin/env python3
"""
Renders the stellated three-ribbon tetrahelix beside the plain tube it grew from,
one colour per ribbon.

A glued cell takes the colour of the ribbon its base face belonged to, so the three
ribbons are still legible after stellation -- the spikes come in three helical bands.
Nothing was fitted to make them fit: each is an exactly regular tetrahedron on an
exactly equilateral face, and the whole layer is collision-free.

The last tile is the limit. A second layer would add 126 more cells and 77 of them
collide, so one layer is all there is: this doubles the sculpture, it does not
triple it.

Run:  python preview_stellated_tetrahelix.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import render_mesh as RM
import ribbon_tetrahelix as RT

TW, TH_PX = 340, 950

RM.MAT["spkA"] = {"col": (86, 226, 122), "amb": 0.38, "spec": 0.28, "shin": 22.0}
RM.MAT["spkB"] = {"col": (236, 132, 92), "amb": 0.38, "spec": 0.28, "shin": 22.0}
RM.MAT["spkC"] = {"col": (112, 168, 240), "amb": 0.38, "spec": 0.28, "shin": 22.0}


def stell_layers():
    m, v, faces, outer, spikes, home, strips, rf, _g = RT.build_stellated()
    V, F = m.arrays()
    tag = np.asarray(m.tag)
    lay = [(V, F[np.isin(tag, ["bore", "fillet"])], "frame")]
    for i, mat in enumerate(("spkA", "spkB", "spkC")):
        sel = F[tag == f"spike{i}"]
        if len(sel):
            lay.append((V, sel, mat, False))
    return lay, V


def plain_layers():
    m, v, faces, home, strips, _g = RT.build(3)
    V, F = m.arrays()
    tag = np.asarray(m.tag)
    lay = [(V, F[np.isin(tag, ["wall", "fillet"])], "frame")]
    for i, mat in enumerate(("spkA", "spkB", "spkC")):
        sel = F[tag == f"ribbon{i}"]
        if len(sel):
            lay.append((V, sel, mat, False))
    return lay, V


def main():
    st, vs = stell_layers()
    pl, vp = plain_layers()
    mids = (0.0, 0.0, 0.5 * (vs[:, 2].min() + vs[:, 2].max()))
    midp = (0.0, 0.0, 0.5 * (vp[:, 2].min() + vp[:, 2].max()))
    span = 2.25 * float(np.linalg.norm(vs[:, :2], axis=1).max())

    views = [
        ("plain tube, 43 mm across", 20, 6, midp, span, pl),
        ("STELLATED, 86 mm across", 20, 6, mids, span, st),
        ("three-quarter", 62, 22, mids, None, st),
        ("detail: 42 glued cells, 3 bands", 35, 10, (0.0, 0.0, mids[2]), 130.0, st),
        ("axial: bore still open", 90, 87, mids, None, st),
        ("the notch: 7.356 deg on k->k+1", 55, 4, (0.0, 0.0, mids[2]), 78.0, st),
    ]
    scratch = os.path.join(os.environ.get("TEMP", "."), "stellated_views")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, target, sp, lz) in enumerate(views):
        out = os.path.join(scratch, f"st{i}.png")
        RM.render(lz, az, el, w=TW, h=TH_PX, target=target, span=sp, out=out)
        im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TH_PX - 20, TW, TH_PX], fill=(238, 238, 235))
        d.text((8, TH_PX - 15), label, fill=(60, 60, 60))
        tiles.append(im)
    sheet = Image.new("RGB", (len(tiles) * TW, TH_PX), (242, 242, 240))
    for i, im in enumerate(tiles):
        sheet.paste(im, (i * TW, 0))
    sheet.save("stellated_tetrahelix_views.png")
    print("wrote stellated_tetrahelix_views.png", sheet.size)


if __name__ == "__main__":
    main()
