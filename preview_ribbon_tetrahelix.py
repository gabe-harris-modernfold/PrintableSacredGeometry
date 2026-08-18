#!/usr/bin/env python3
"""
Renders the one-, two- and three-ribbon tetrahelix side by side, one colour per
ribbon, to show they are the same tube cut three ways.

The grey in the grooves is the wall underneath. What changes between the three tiles
is only which edge family got cut: k->k+1 leaves a single ribbon that spirals the
whole length, k->k+2 leaves two, k->k+3 leaves three. Nothing about the tetrahedra
changes -- same 21 regular cells, same tube, same 45 mm across.

The axial tile is the bore: the interior faces are left out, so it is open end to
end, which is what makes it a folded-ribbon model rather than a solid chain.

Run:  python preview_ribbon_tetrahelix.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import render_mesh as RM
import ribbon_tetrahelix as RT

TW, TH_PX = 330, 950

RM.MAT["ribA"] = {"col": (86, 226, 122), "amb": 0.40, "spec": 0.26, "shin": 20.0}
RM.MAT["ribB"] = {"col": (236, 132, 92), "amb": 0.40, "spec": 0.26, "shin": 20.0}
RM.MAT["ribC"] = {"col": (112, 168, 240), "amb": 0.40, "spec": 0.26, "shin": 20.0}


def layers_for(ribbons):
    m, v, faces, home, strips, _g = RT.build(ribbons)
    V, F = m.arrays()
    tag = np.asarray(m.tag)
    lay = [(V, F[tag == "wall"], "frame")]
    for i, mat in enumerate(("ribA", "ribB", "ribC")):
        sel = F[tag == f"ribbon{i}"]
        if len(sel):
            lay.append((V, sel, mat, False))
    return lay, V


def main():
    one, v1 = layers_for(1)
    two, v2 = layers_for(2)
    three, v3 = layers_for(3)
    mid = (0.0, 0.0, 0.5 * (v3[:, 2].min() + v3[:, 2].max()))

    views = [
        ("ONE ribbon: cut k->k+1", 20, 6, mid, None, one),
        ("TWO ribbons: cut k->k+2", 20, 6, mid, None, two),
        ("THREE ribbons: cut k->k+3", 20, 6, mid, None, three),
        ("three-quarter", 62, 22, mid, None, three),
        ("detail: seam grooves = the 3 helices", 35, 8,
         (0.0, 0.0, mid[2]), 110.0, three),
        ("axial: open bore, 38 mm", 90, 87, mid, None, three),
    ]
    scratch = os.path.join(os.environ.get("TEMP", "."), "ribbon_tetrahelix_views")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, target, span, lz) in enumerate(views):
        out = os.path.join(scratch, f"rt{i}.png")
        RM.render(lz, az, el, w=TW, h=TH_PX, target=target, span=span, out=out)
        im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TH_PX - 20, TW, TH_PX], fill=(238, 238, 235))
        d.text((8, TH_PX - 15), label, fill=(60, 60, 60))
        tiles.append(im)
    sheet = Image.new("RGB", (len(tiles) * TW, TH_PX), (242, 242, 240))
    for i, im in enumerate(tiles):
        sheet.paste(im, (i * TW, 0))
    sheet.save("ribbon_tetrahelix_views.png")
    print("wrote ribbon_tetrahelix_views.png", sheet.size)


if __name__ == "__main__":
    main()
