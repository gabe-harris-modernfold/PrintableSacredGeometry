#!/usr/bin/env python3
"""
Renders the (3,7) phyllotactic equilateral vortex tube: frame, glass and beam.

The beam shown is the TIME-REVERSED outbound trace, so every arm terminates
exactly on the focus at the base centre. Orbit views cull the glass facets on the
camera side, otherwise the mirrors hide the vortex; the frame is never culled by
centroid because that just shreds a space frame into speckle.

Run:  python preview_equilateral.py
"""

import math
import os

import numpy as np
from PIL import Image, ImageDraw

import equilateral_frame as EF
import equilateral_tube as ET
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import mesh_kit as MK
import ray_optics as RO
import render_mesh as RM

GLASS = 1.5
GAP = 0.4               # butt gap between mirrors, so the tiling reads
TW, TH = 440, 700


def glass_mesh(tris, nrm, keep=None):
    m = MK.Mesh()
    for i, (t, n) in enumerate(zip(tris, nrm)):
        if keep is not None and not keep(i):
            continue
        c = t.mean(axis=0)
        front = [c + (q - c) * max(0.0, (np.linalg.norm(q - c) - GAP)
                                   / np.linalg.norm(q - c)) for q in t]
        back = [q - n * GLASS for q in front]
        faces = [(0, 2, 1), (3, 4, 5)]
        for a in range(3):
            b = (a + 1) % 3
            faces += [(a, b, 3 + b), (a, 3 + b, 3 + a)]
        m.add_solid(list(front) + list(back), faces, tag="glass", wall=i)
    return m


def beam_mesh(tris, nrm, sol, r=1.0):
    """All arms, each an outbound trace from the focus reversed."""
    m = MK.Mesh()
    n_arms = 0
    for k in range(6):
        p0, d0 = RO.launch_from_focus(sol["theta"], az=360.0 * k / 6)
        pts, log = RO.shoot(p0, d0, tris, nrm, 6, 90)
        if len(log) < 3:
            continue
        n_arms += 1
        for a, b in zip(pts[::-1][:-1], pts[::-1][1:]):
            if np.linalg.norm(b - a) > 1e-6:
                m.add_solid(*MK.tube(a, b, r, nseg=8), tag="beam")
    return m, n_arms


def main():
    sol = EF.load()
    frame, verts, faces, fn = EF.build(sol)
    tris, nrm, _ = ET.build_helix(sol["R"], sol["alpha"], sol["beta"],
                                  sol["p"], sol["r"], sol["K"])
    H = (sol["K"] - 1) * sol["beta"]
    g = glass_mesh(tris, nrm)
    b, n_arms = beam_mesh(tris, nrm, sol)
    print(f"({sol['p']},{sol['r']}) tube: {len(tris)} equilateral {sol['s']} mm "
          f"mirrors, {n_arms} arms, H = {H:.0f} mm")

    F = (*frame.arrays(), "frame")
    G = (*g.arrays(), "glass")
    B = (*b.arrays(), "beam", False)
    mid = (0.0, 0.0, 0.45 * H)
    near = lambda a: (lambda c: (c[:, :2] @ np.array(
        [math.cos(math.radians(a)), math.sin(math.radians(a))])) > 0.0)

    sheets = {
        "equilateral_frame_views.png": [
            ("printed frame only", [F], 28, 15, None),
            ("assembled", [F, G], 28, 15, None),
            ("cut open, beam", [G, B], 10, 14, near(10)),
            ("down the axis", [G, B], 90, 86, None),
        ],
        "equilateral_vortex.png": [
            ("cut open, az 0deg", [G, B], 0, 15, near(0)),
            ("cut open, az 60deg", [G, B], 60, 15, near(60)),
            ("down the axis", [G, B], 90, 86, None),
            ("assembled, no cut", [G], 30, 16, None),
        ],
    }
    scratch = os.path.join(os.environ.get("TEMP", "."), "eq_views")
    os.makedirs(scratch, exist_ok=True)
    for name, views in sheets.items():
        tiles = []
        for i, (label, layers, az, el, cull) in enumerate(views):
            out = os.path.join(scratch, f"{name[:6]}{i}.png")
            RM.render(layers, az, el, w=TW, h=TH, target=mid, cull=cull, out=out)
            im = Image.open(out).convert("RGB")
            d = ImageDraw.Draw(im)
            d.rectangle([0, TH - 20, TW, TH], fill=(238, 238, 235))
            d.text((8, TH - 15), label, fill=(60, 60, 60))
            tiles.append(im)
        sheet = Image.new("RGB", (len(tiles) * TW, TH), (242, 242, 240))
        for i, im in enumerate(tiles):
            sheet.paste(im, (i * TW, 0))
        sheet.save(name)
        print("wrote", name, sheet.size)


if __name__ == "__main__":
    main()
