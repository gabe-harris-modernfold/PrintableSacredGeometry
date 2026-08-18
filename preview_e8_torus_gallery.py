#!/usr/bin/env python3
"""Object gallery for the E8 horn torus lattice -- the thing itself rather than the argument.

Eight panels: a hero three-quarter, down the axis, edge on, a wedge seen down the tube, the
outer coil isolated, a joint out on the rim, the eight coils alone coloured by golden-ratio
pair, and the braces alone with the coils stripped out.

The wedge is cut by whole solid, not by triangle. A centroid-plane face cull through a space
frame does not open it up, it shreds it into speckle -- so each coil, brace and node is kept
or dropped entire, on the side of the plane its own centroid falls. Cutting it in half and
looking at the cut face does not work either: the far half of what is kept stands in the way.

Run:  python preview_e8_torus_gallery.py
"""

import os

import numpy as np
from PIL import Image, ImageDraw

import e8_torus_lattice as E8
import render_mesh as RM
from preview_e8_torus import PAIR_COLS

TW = 620
BG = (242, 242, 240)

for n, c in enumerate(PAIR_COLS):                      # runtime-only materials
    RM.MAT[f"pair{n}"] = {"col": c, "amb": 0.42, "spec": 0.22, "shin": 18.0}
RM.MAT["edge"] = {"col": (150, 122, 96), "amb": 0.34, "spec": 0.14, "shin": 16.0}


def solid_side(mesh, keep):
    """Face indices of the solids whose own centroid satisfies keep(centroid)."""
    v = np.asarray(mesh.v, float)
    starts = mesh._solid_start
    out = np.zeros(len(mesh.f), bool)
    for n, (base, f0, nf) in enumerate(starts):
        end = starts[n + 1][0] if n + 1 < len(starts) else len(v)
        if keep(v[base:end].mean(0)):
            out[f0:f0 + nf] = True
    return out


def main():
    mesh, info = E8.build(verbose=False, metrics=False)
    v, f = mesh.arrays()
    tagc = np.asarray(mesh.tag)
    mid = (0.0, 0.0, 0.5 * (v[:, 2].min() + v[:, 2].max()))
    nb = len(info["braces"])

    whole = [(*mesh.pick(tags=["coil", "node"]), "frame"),
             (*mesh.pick(tags=["brace"]), "pad")]

    wedge = solid_side(mesh, lambda c: abs(np.degrees(np.arctan2(c[1], c[0]))) <= 55.0)
    cut = [(v, f[wedge & np.isin(tagc, ["coil", "node"])], "frame", False),
           (v, f[wedge & (tagc == "brace")], "pad", False)]

    pair_of = {}
    for n, (i, j) in enumerate(sorted(info["pairs"])):
        pair_of[i] = pair_of[j] = n
    coils_only = [(*mesh.pick(tags=["coil", "node"], walls=[k]), f"pair{pair_of[k]}")
                  for k in range(8)]
    braces_only = [(*mesh.pick(tags=["brace", "node"]), "edge")]
    # The 7-fold hub cannot be photographed in situ -- it sits at the dead centre of a dense
    # lattice, so a close-up of it is just a wall of struts. Isolate the curve that actually
    # goes through the axis instead: only the outermost shell reaches rho = 0, and it does so
    # once per poloidal turn, so this one closed coil crosses itself 7 times at the origin.
    outer = [(*mesh.pick(tags=["coil", "node"], walls=list(range(7))), "frame"),
             (*mesh.pick(tags=["coil", "node"], walls=[7]), "beam", False)]

    P = info["P"]
    rho = np.hypot(P[:, 0], P[:, 1])
    rim = tuple(P[int(np.argmax(rho))])
    assert rho.min() < 1e-6, "no root on the axis -- not a horn torus?"

    views = [
        (f"150 mm across, 240 roots, 8 coils, {nb} braces", 52, 24, mid, None, whole),
        ("down the axis: the throat has closed", 90, 89.5, mid, None, whole),
        ("edge on", 0, 2, mid, None, whole),
        ("a 55 deg wedge, seen down the tube", 0, 6, mid, 120.0, cut),
        ("the outer coil alone: through the axis 7 times", 52, 24, mid, None, outer),
        ("a joint out on the rim, 62 mm across", 46, 20, rim, 62.0, whole),
        ("the 8 coils alone, by golden-ratio pair", 52, 24, mid, None, coils_only),
        (f"the {nb} E8 braces alone, coils stripped", 52, 24, mid, None, braces_only),
    ]

    scratch = os.path.join(os.environ.get("TEMP", "."), "e8_torus_gallery")
    os.makedirs(scratch, exist_ok=True)
    tiles = []
    for i, (label, az, el, target, span, layers) in enumerate(views):
        out = os.path.join(scratch, f"g{i}.png")
        RM.render(layers, az, el, w=TW, h=TW, target=target, span=span, out=out)
        im = Image.open(out).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([0, TW - 20, TW, TW], fill=(238, 238, 235))
        d.text((8, TW - 15), label, fill=(60, 60, 60))
        tiles.append(im)

    cols = 4
    rows = (len(tiles) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * TW, rows * TW), BG)
    for i, im in enumerate(tiles):
        sheet.paste(im, ((i % cols) * TW, (i // cols) * TW))
    sheet.save("e8_torus_lattice_views.png")
    print("wrote e8_torus_lattice_views.png", sheet.size)

    RM.render(whole, 48, 22, w=1500, h=1050, target=mid, out="e8_torus_lattice_hero.png")
    print("wrote e8_torus_lattice_hero.png")


if __name__ == "__main__":
    main()
