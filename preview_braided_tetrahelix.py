#!/usr/bin/env python3
"""
Renders the three-tetrahelix rope, one colour per WHOLE tetrahelix, and the thing
the module is for: the stack.

Each colour is a complete Boerdijk-Coxeter column -- its own 39 cells, its own 3V-6
beams -- bent onto the lay. Grey rods are the ties between plies. If a colour ever
appeared to hand off to another the construction would be wrong; each runs the full
length.

The stack tiles are the argument. Three modules, each rotated by the module's own
lay advance and raised by its height, and the braid runs on through both seams: the
cut faces meet flat, every cell is whole, every consecutive pair still shares its
face. The seam tile is deliberately at a hard zoom so you can look for a step and
not find one.

Last tile is the P=30 period, where the lay goes almost straight and psi lands on
132.000 deg -- Fuller's ten-cycle, recovered as a convergent rather than a fudge.

Run:  python preview_braided_tetrahelix.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import braided_tetrahelix as BT
import render_mesh as RM

TW, TH_PX = 330, 950

RM.MAT["strandA"] = {"col": (86, 226, 122), "amb": 0.40, "spec": 0.26, "shin": 20.0}
RM.MAT["strandB"] = {"col": (236, 132, 92), "amb": 0.40, "spec": 0.26, "shin": 20.0}
RM.MAT["strandC"] = {"col": (112, 168, 240), "amb": 0.40, "spec": 0.26, "shin": 20.0}


def layers_for(stack=1, period=BT.PERIOD, steps=BT.STEPS, base=False):
    m, rope, _cut = BT.build(BT.EDGE, steps, period, BT.CLEAR, base=base,
                             stack=stack)
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    root = np.array([t.split("#")[0] for t in tag])
    lay = [(v, f[root == "tie"], "frame", False),
           (v, f[root == "base"], "pad")]
    # one colour per tetrahelix, held across every module in the stack
    lay += [(v, f[np.isin(root, [f"strand{s}", f"node{s}"])], mat, False)
            for s, mat in enumerate(("strandA", "strandB", "strandC"))]
    return lay, v, rope


def main():
    one, v1, r = layers_for(1, base=True)
    three, v3, _ = layers_for(3)
    tencyc, v30, r30 = layers_for(1, period=30, steps=30)
    mid1 = (0.0, 0.0, 0.5 * (v1[:, 2].min() + v1[:, 2].max()))
    mid3 = (0.0, 0.0, 0.5 * (v3[:, 2].min() + v3[:, 2].max()))
    mid30 = (0.0, 0.0, 0.5 * (v30[:, 2].min() + v30[:, 2].max()))
    seam = (0.0, 0.0, r.height)

    views = [
        ("one module, 39 cells per ply", 20, 6, mid1, None, one),
        ("axial: 3-fold, plies clear", 90, 87, mid1, None, one),
        ("THREE STACKED, +90.8 deg each", 20, 6, mid3, None, three),
        ("seam close-up: no step, no gap", 35, 4, seam, 90.0, three),
        ("cells at a seam, harder zoom", 60, 8, seam, 52.0, three),
        ("--period 30: psi = 132, the ten-cycle", 20, 6, mid30, None, tencyc),
    ]
    scratch = os.path.join(os.environ.get("TEMP", "."), "braided_tetrahelix_views")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, target, span, lz) in enumerate(views):
        out = os.path.join(scratch, f"bt{i}.png")
        RM.render(lz, az, el, w=TW, h=TH_PX, target=target, span=span, out=out)
        im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TH_PX - 20, TW, TH_PX], fill=(238, 238, 235))
        d.text((8, TH_PX - 15), label, fill=(60, 60, 60))
        tiles.append(im)
    sheet = Image.new("RGB", (len(tiles) * TW, TH_PX), (242, 242, 240))
    for i, im in enumerate(tiles):
        sheet.paste(im, (i * TW, 0))
    sheet.save("braided_tetrahelix_views.png")
    print("wrote braided_tetrahelix_views.png", sheet.size)


if __name__ == "__main__":
    main()
