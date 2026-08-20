#!/usr/bin/env python3
"""Contact sheet for stellated_platonics.py.

Row 1  the Platonic solid.
Row 2  the one-level stellation grown on it, in the printed orientation.
Row 3  the same shell cut through the centre, so the 4 mm wall is visible.

All three rows share one scale and one camera, so the five are honestly comparable.
The two solids that cannot stellate are tinted blue-grey and carry their dual
compound instead.

Run:  python preview_stellated_platonics.py
"""

import os

from PIL import Image, ImageDraw

import render_mesh as RM
import stellated_platonics as SP

TW, TH = 420, 460
BAR = 22
AZ, EL = 28, 18

RM.MAT["core"] = {"col": (120, 126, 136), "amb": 0.34, "spec": 0.16, "shin": 18.0}
RM.MAT["star"] = {"col": (214, 168, 74), "amb": 0.34, "spec": 0.14, "shin": 30.0}
RM.MAT["subs"] = {"col": (128, 150, 168), "amb": 0.34, "spec": 0.12, "shin": 26.0}
RM.MAT["void"] = {"col": (74, 70, 66), "amb": 0.30, "spec": 0.06, "shin": 12.0}

CAPTIONS = {
    "tetrahedron": "tetrahedron  D=70.53",
    "cube": "cube  D=90.00",
    "octahedron": "octahedron  D=109.47",
    "dodecahedron": "dodecahedron  D=116.57",
    "icosahedron": "icosahedron  D=138.19",
}
RESULTS = {
    "tetrahedron": "no stellation - dual compound",
    "cube": "no stellation - dual compound",
    "octahedron": "stella octangula",
    "dodecahedron": "small stellated dodecahedron",
    "icosahedron": "small triambic icosahedron",
}


def tile(layers, span, target, label, out, cull=None):
    RM.render(layers, az=AZ, el=EL, w=TW, h=TH - BAR, target=target, span=span,
              cull=cull, out=out)
    im = Image.new("RGB", (TW, TH), (242, 242, 240))
    im.paste(Image.open(out).convert("RGB"), (0, 0))
    ImageDraw.Draw(im).text((8, TH - 16), label, fill=(60, 60, 60))
    return im


def main():
    scratch = os.path.join(os.environ.get("TEMP", "."), "stellated_platonics_views")
    os.makedirs(scratch, exist_ok=True)

    # the stella octangula's tetrahedra have edge 150*sqrt(2), so the frame has to
    # hold well over the 150 mm print size
    span = 2.25 * SP.SIZE
    target = (0.0, 0.0, 0.5 * SP.SIZE)
    # cut on x, which the camera faces almost square-on, so the section rim
    # reads as a band rather than a sliver
    half = lambda c: c[:, 0] > 0.0

    rows = [[], [], []]
    for name in SP.SOLIDS:
        r = SP.build(name)
        mat = "star" if r["true"] else "subs"
        s = SP.scale_factor(r, SP.SIZE)
        outerV, outerF = r["verts"] * s, r["tris"]
        lift = -outerV[:, 2].min()                   # the same drop() the STL gets
        coreV, coreF = SP.core_triangles(r)
        shellV, shellF, sh = SP.hollow(outerV, outerF, SP.WALL)

        nout = len(outerF)
        plan = [([(coreV * s, coreF, "core")], CAPTIONS[name], None),
                ([(outerV, outerF, mat)], RESULTS[name], None),
                # the cavity lining goes dark, so the lit rim between the two
                # surfaces IS the wall
                ([(shellV, shellF[nout:], "void"), (shellV, shellF[:nout], mat)],
                 f"cut: {sh['wall_min']:.1f} mm wall, hollow", half)]
        for row, (layers, cap, cull) in enumerate(plan):
            lifted = []
            for V, F, mat_ in layers:
                V = V.copy()
                V[:, 2] += lift
                lifted.append((V, F, mat_))
            out = os.path.join(scratch, f"{name}_{row}.png")
            rows[row].append(tile(lifted, span, target, cap, out, cull))

    sheet = Image.new("RGB", (5 * TW, 3 * TH), (242, 242, 240))
    for ri, tiles in enumerate(rows):
        for ci, im in enumerate(tiles):
            sheet.paste(im, (ci * TW, ri * TH))
    d = ImageDraw.Draw(sheet)
    for ci in range(1, 5):
        d.line([(ci * TW, 0), (ci * TW, 3 * TH)], fill=(214, 214, 210))
    for ri in range(1, 3):
        d.line([(0, ri * TH), (5 * TW, ri * TH)], fill=(214, 214, 210))
    sheet.save("stellated_platonics_views.png")
    print("wrote stellated_platonics_views.png", sheet.size)


if __name__ == "__main__":
    main()
