"""Render preview_waveform_solid.png from run_pure.mjs's waveform STLs: the solid lying whole, and as the two
print halves on their cut faces.

    node run_pure.mjs && python render_waveform_preview.py"""
import os
import sys

import numpy as np
import trimesh
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
import render_mesh as rm  # noqa: E402

rm.MAT["cream"] = {"col": (232, 222, 200), "amb": 0.3, "spec": 0.2, "shin": 18.0}
rm.BG = (11, 13, 20)

PANELS = [("js_wave_fins.stl", -62, 14, "waveform solid: the spoken-vowel demo, 72 fins, thickness by pitch"),
          ("js_wave_fins_halves.stl", -60, 38, "print as halves: cut faces down, every half-disc upright")]


def main():
    tiles = []
    for name, az, el, label in PANELS:
        m = trimesh.load(os.path.join(HERE, name), force="mesh")
        out = os.path.join(HERE, "_tile.png")
        rm.render([(np.asarray(m.vertices, float), np.asarray(m.faces), "cream")], az, el, w=1000, h=620, out=out)
        im = Image.open(out).convert("RGB")
        ImageDraw.Draw(im).text((14, 10), label, fill=(200, 200, 205))
        tiles.append(im)
        os.remove(out)
    sheet = Image.new("RGB", (sum(t.width for t in tiles), max(t.height for t in tiles)), rm.BG)
    x = 0
    for t in tiles:
        sheet.paste(t, (x, 0))
        x += t.width
    path = os.path.join(HERE, "preview_waveform_solid.png")
    sheet.save(path)
    print("wrote", path, sheet.size)


if __name__ == "__main__":
    main()
