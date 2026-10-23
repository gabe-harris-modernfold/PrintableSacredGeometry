#!/usr/bin/env python3
"""The eight stackable plates of the matched set, and the pin that registers them.

Every figure is drawn with circles of radius R = 15 mm on one triangular lattice, inside a
common 170 mm ring frame with three pin holes. Stack any of them on the pins and their
lines fall on one another. Geometry lives in lib/sacred_grid.py.

Run from this directory: python plates.py
  01_circle.stl ... 08_tree_of_life.stl   the plates, 3 mm thick
  plate_pin.stl                           D-section pin, prints lying on its flat (make 3)
  preview_plates.png                      all eight from above
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

from PIL import Image, ImageDraw

import sacred_grid as sg


def draw(outline, size=420, pad=8):
    """Top view of one plate outline, filled, holes knocked out."""
    img = Image.new("RGB", (size, size), (246, 245, 241))
    d = ImageDraw.Draw(img)
    s = (size - 2 * pad) / (2 * sg.FRAME_OUT)

    def px(coords):
        return [(size / 2 + x * s, size / 2 - y * s) for x, y in coords]
    for poly in getattr(outline, "geoms", [outline]):
        d.polygon(px(poly.exterior.coords), fill=(52, 58, 70))
        for hole in poly.interiors:
            d.polygon(px(hole.coords), fill=(246, 245, 241))
    return img


def main():
    ok = True
    tiles = []
    for k, name in enumerate(sg.PLATES, 1):
        path = f"{k:02d}_{name}.stl"
        mesh = sg.to_trimesh(sg.plate(name))
        ok &= sg.check(mesh, path)
        mesh.export(path)
        tiles.append((name.replace("_", " "), draw(sg.plate_outline(name))))
    pin = sg.to_trimesh(sg.pin())
    ok &= sg.check(pin, "plate_pin.stl")
    pin.export("plate_pin.stl")

    sheet = Image.new("RGB", (4 * 420, 2 * 448), (246, 245, 241))
    d = ImageDraw.Draw(sheet)
    for k, (label, img) in enumerate(tiles):
        x, y = (k % 4) * 420, (k // 4) * 448
        sheet.paste(img, (x, y))
        d.text((x + 14, y + 424), f"{k + 1:02d}  {label}", fill=(60, 64, 72))
    sheet.save("preview_plates.png")
    print("wrote preview_plates.png")
    print(f"\nR = {sg.R} mm; frame {2 * sg.FRAME_OUT:.0f} mm; plates {sg.PLATE_T} mm; "
          f"pins {sg.PIN_D} mm in {sg.PIN_HOLE} mm holes at {sg.PIN_ANGLES} deg")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
