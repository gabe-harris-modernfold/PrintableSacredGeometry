#!/usr/bin/env python3
"""The 8 Coxeter coils of the E8 horn torus lattice, alone, coloured by golden-ratio pair.

The braces are stripped out so the coils can be seen. Each of the 8 Coxeter orbits of 30
roots is one closed (1, 7) coil, and the 8 tube radii they sit at form four exact
golden-ratio pairs, one colour each:

  7.89 : 12.76    15.69 : 25.39    18.97 : 30.69    23.32 : 37.73 mm

That pairing is the E8 -> H4 folding: E8's roots project to two 600-cells, one phi times the
size of the other, so the shells pair off the same way. Four of the panels show one pair on
its own, the smaller radius in the lighter tone, which is the only way to actually see the
ratio rather than be told it. The flat panel is the same four colours on the Coxeter-plane
projection, where those radii come from in the first place.

The outermost coil, gold, is the one that reaches the axis -- on a horn torus a = R, so its
inner equator collapses to the origin and it passes through that single point 7 times.

Run:  python preview_e8_coils.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import e8_torus_lattice as E8
import preview_e8_torus as PT
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import render_mesh as RM

TW = 700
BG = (242, 242, 240)


def lighten(c, f=0.52):
    return tuple(int(round(x + (245 - x) * f)) for x in c)


for n, c in enumerate(PT.PAIR_COLS):
    RM.MAT[f"pair{n}"] = {"col": c, "amb": 0.42, "spec": 0.22, "shin": 18.0}
    RM.MAT[f"pair{n}lo"] = {"col": lighten(c), "amb": 0.50, "spec": 0.20, "shin": 18.0}


def main():
    mesh, info = E8.build(verbose=False, metrics=False)
    v = np.asarray(mesh.v, float)
    mid = (0.0, 0.0, 0.5 * (v[:, 2].min() + v[:, 2].max()))
    rad, _ = E8.cluster(info["t_of"], tol=1e-6)
    pairs = sorted(info["pairs"])                  # (smaller shell, phi * it)

    pair_of = {}
    for n, (i, j) in enumerate(pairs):
        pair_of[i] = pair_of[j] = n

    # coil + its own 30 roots share a wall index, which is the shell index
    def coils(walls, mats):
        return [(*mesh.pick(tags=["coil", "node"], walls=[k]), m)
                for k, m in zip(walls, mats)]

    every = coils(range(8), [f"pair{pair_of[k]}" for k in range(8)])

    # the Coxeter-plane diagram is written for a 460x620 tile; borrow it at this tile size
    PT.TW, PT.TH = TW, TW

    views = [("all 8 coils, three-quarter", 52, 24, every),
             ("all 8, down the axis", 90, 89.5, every),
             ("all 8, edge on", 0, 2, every),
             (None, 0, 0, None)]
    for n, (i, j) in enumerate(pairs):
        views.append((f"{rad[i]:.2f} : {rad[j]:.2f} mm   ratio {rad[j] / rad[i]:.6f}",
                      52, 24, coils([i, j], [f"pair{n}lo", f"pair{n}"])))

    scratch = os.path.join(os.environ.get("TEMP", "."), "e8_coils")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, layers) in enumerate(views):
        if layers is None:
            im = PT.coxeter_plane_panel(info)
            label = "where the radii come from: the Coxeter plane"
        else:
            out = os.path.join(scratch, f"c{i}.png")
            # every panel shares one span so the four pairs are directly comparable
            RM.render(layers, az, el, w=TW, h=TW, target=mid, span=170.0, out=out)
            im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TW - 22, TW, TW], fill=(238, 238, 235))
        d.text((9, TW - 16), label, fill=(60, 60, 60))
        tiles.append(im)

    cols = 4
    rows = (len(tiles) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * TW, rows * TW), BG)
    for i, im in enumerate(tiles):
        sheet.paste(im, ((i % cols) * TW, (i // cols) * TW))
    sheet.save("e8_coils_golden_pairs.png")
    print("wrote e8_coils_golden_pairs.png", sheet.size)

    RM.render(every, 50, 23, w=1600, h=1120, target=mid,
              out="e8_coils_golden_pairs_hero.png")
    print("wrote e8_coils_golden_pairs_hero.png")


if __name__ == "__main__":
    main()
