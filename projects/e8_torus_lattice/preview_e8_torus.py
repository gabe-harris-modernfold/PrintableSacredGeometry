#!/usr/bin/env python3
"""Preview sheet for the E8 torus lattice.

Six panels. The last three carry the argument.

The Coxeter-plane panel is drawn, not rendered: it is the flat picture of what the tube's
cross-section is. Projecting all 240 roots onto the m = 1 eigenplane of the Coxeter element
puts them on 8 concentric rings of 30 -- the Petrie projection -- and those eight radii,
colour-paired here, are exactly the eight shell radii of the print. Each colour is a
golden-ratio pair.

Then one orbit alone, green. Three-quarters on, it is a coil that goes once around the torus
while winding 7 times around the tube. From above, its 30 roots sit at 30 azimuths exactly
12 degrees apart -- one thirtieth of a turn per Coxeter step, which is the same statement as
"the Coxeter element has order 30". One curve, both facts.

Run:  python preview_e8_torus.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import e8_torus_lattice as E8
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import render_mesh as RM

TW, TH = 460, 620
BG = (242, 242, 240)
PAIR_COLS = [(199, 84, 80), (74, 132, 176), (86, 152, 96), (198, 146, 62)]


def coxeter_plane_panel(info):
    """The Petrie projection: 240 roots, 8 rings of 30, one colour per golden-ratio pair."""
    im = Image.new("RGB", (TW, TH), BG)
    d = ImageDraw.Draw(im)
    cx, cy = TW / 2, (TH - 20) / 2
    rad, _ = E8.cluster(info["t_of"], tol=1e-6)
    s = 0.40 * TW / rad.max()
    col_of = {}
    for c, (i, j) in zip(PAIR_COLS, sorted(info["pairs"])):
        col_of[i] = col_of[j] = c

    for k, r in enumerate(rad):
        rr = r * s
        d.ellipse([cx - rr, cy - rr, cx + rr, cy + rr], outline=(214, 214, 210))

    # the Petrie polygon itself: the outer ring's 30 vertices in azimuth order
    P, shell = info["P"], info["shell"]
    th = np.arctan2(P[:, 1], P[:, 0])
    for k, r in enumerate(rad):
        a = np.sort(th[shell == k])
        pts = [(cx + r * s * np.cos(t), cy - r * s * np.sin(t)) for t in a]
        d.line(pts + pts[:1], fill=(196, 196, 192) if k < 7 else (150, 150, 146),
               width=2 if k == 7 else 1)
    for k, r in enumerate(rad):
        for t in th[shell == k]:
            x, y = cx + r * s * np.cos(t), cy - r * s * np.sin(t)
            d.ellipse([x - 3.4, y - 3.4, x + 3.4, y + 3.4], fill=col_of[k],
                      outline=(250, 250, 248))

    d.text((10, 10), "240 roots -> 8 rings of 30", fill=(70, 70, 70))
    for n, (c, (i, j)) in enumerate(zip(PAIR_COLS, sorted(info["pairs"]))):
        y = 30 + 16 * n
        d.rectangle([10, y + 3, 20, y + 11], fill=c)
        d.text((26, y), f"{rad[i]:5.2f} : {rad[j]:5.2f} mm   x phi", fill=(70, 70, 70))
    return im


def main():
    mesh, info = E8.build(verbose=False, metrics=False)
    v = np.asarray(mesh.v, float)
    mid = (0.0, 0.0, 0.5 * (v[:, 2].min() + v[:, 2].max()))

    whole = [(*mesh.pick(tags=["coil", "node"]), "frame"),
             (*mesh.pick(tags=["brace"]), "pad")]

    # one Coxeter orbit: its coil and its own 30 roots share a wall index (the shell)
    pick = 6
    orbit = [(*mesh.pick(walls=[w for w in range(9) if w != pick]), "frame"),
             (*mesh.pick(walls=[pick]), "beam", False)]

    nb = len(info["braces"])
    views = [                                       # (label, az, el, layers) -- None = drawn
        (f"the lattice: 8 coils, 240 roots, {nb} braces", 52, 26, whole),
        ("down the torus axis", 90, 89.5, whole),
        ("edge on: 8 nested shells inside the tube", 0, 3, whole),
        ("the Coxeter plane, which is the tube cross-section", 0, 0, None),
        ("one Coxeter orbit: 30 roots on a (1,7) coil", 52, 26, orbit),
        ("the same orbit from above: 30 azimuths, 12 deg apart", 90, 89.5, orbit),
    ]

    scratch = os.path.join(os.environ.get("TEMP", "."), "e8_torus_views")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, layers) in enumerate(views):
        if layers is None:
            im = coxeter_plane_panel(info)
        else:
            out = os.path.join(scratch, f"e8_{i}.png")
            RM.render(layers, az, el, w=TW, h=TH, target=mid, out=out)
            im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TH - 20, TW, TH], fill=(238, 238, 235))
        d.text((8, TH - 15), label, fill=(60, 60, 60))
        tiles.append(im)

    cols = 3
    rows = (len(tiles) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * TW, rows * TH), BG)
    for i, im in enumerate(tiles):
        sheet.paste(im, ((i % cols) * TW, (i // cols) * TH))
    sheet.save("e8_torus_lattice.png")
    print("wrote e8_torus_lattice.png", sheet.size)


if __name__ == "__main__":
    main()
