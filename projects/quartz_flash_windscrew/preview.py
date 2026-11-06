#!/usr/bin/env python3
"""Preview renders.  Run after windscrew.py, from this directory:  python preview.py"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import render_mesh
from render_mesh import render

render_mesh.MAT["brown"] = {'col': (82, 54, 38), 'amb': 0.36, 'spec': 0.22, 'shin': 18.0}     # walnut / Bakelite
render_mesh.MAT["clear"] = {'col': (196, 214, 222), 'amb': 0.40, 'spec': 0.9, 'shin': 60.0}
render_mesh.MAT["steel"] = {'col': (200, 204, 210), 'amb': 0.3, 'spec': 0.9, 'shin': 60.0}
render_mesh.MAT["brass"] = {'col': (201, 162, 82), 'amb': 0.42, 'spec': 0.75, 'shin': 40.0}

def load():
    d = np.load("assembly.npz")
    keys = sorted({k[:-2] for k in d.files})
    out = []
    for k in keys:
        _, *name, mat = k.split("_")
        out.append(("_".join(name), d[k + "_v"], d[k + "_f"], "clear" if mat == "glass" else mat))
    return out

def layers(parts, only=None, skip=()):
    L = []
    for name, v, f, mat in parts:
        if (only and name not in only) or name in skip:
            continue
        if name == "comm_disc":
            # everything above the dial face prints after the filament swap
            fc = v[f].mean(axis=1)
            top = fc[:, 0] > v[:, 0].min() + 4.05
            L.append((v, f[~top], "brown"))
            L.append((v, f[top], "brass"))
            continue
        L.append((v, f, mat))
    return L

def sheet(P):
    """One-page general-arrangement sheet: four views, captioned."""
    from PIL import Image, ImageDraw, ImageFont
    views = [("Three-quarter, downwind side", -35, 16), ("Upwind, from the nose", 215, 14),
             ("Side elevation", -90, 2), ("Dial end", 0, 6)]
    tiles = []
    for i, (cap, az, el) in enumerate(views):
        tmp = f"_sheet_{i}.png"
        render(layers(P), az=az, el=el, w=1000, h=1150, out=tmp)
        tiles.append((cap, Image.open(tmp).convert("RGB")))
        os.remove(tmp)
    W, H, pad, head = 1000, 1150, 24, 70
    img = Image.new("RGB", (2 * W + 3 * pad, head + 2 * (H + 46) + 2 * pad), (242, 242, 240))
    d = ImageDraw.Draw(img)
    try:
        f_title = ImageFont.truetype(r"C:\Windows\Fonts\GOTHICB.TTF", 40)
        f_cap = ImageFont.truetype(r"C:\Windows\Fonts\GOTHIC.TTF", 26)
    except OSError:
        f_title = f_cap = ImageFont.load_default()
    d.text((pad, 16), "QUARTZ FLASH WINDSCREW  \u00b7  MCMXXVI", fill=(82, 54, 38), font=f_title)
    for i, (cap, t) in enumerate(tiles):
        x = pad + (i % 2) * (W + pad)
        y = head + pad + (i // 2) * (H + 46)
        img.paste(t, (x, y))
        d.text((x + 8, y + H + 8), cap, fill=(60, 60, 60), font=f_cap)
    img.save("preview_sheet.png")

if __name__ == "__main__":
    P = load()
    zax = 554.0
    render(layers(P), az=-35, el=16, w=1200, h=1400, out="preview_assembly.png")
    render(layers(P), az=-90, el=2, w=1400, h=1100, out="preview_side.png")
    render(layers(P, only=("comm_disc", "pointer", "finial_disc")), az=0, el=0, w=1100, h=1100,
           out="preview_disc.png")
    render(layers(P, skip=("comm_disc", "pointer", "finial_disc")), az=-60, el=-10, w=1100, h=900,
           target=(265, 0, zax), span=260, out="preview_tumbler.png")
    render(layers(P, only=("base", "post", "yaw_hub", "yoke_nose", "yoke_mouth")), az=-90, el=8,
           w=1300, h=1100, target=(0, 0, 300), span=560, out="preview_frame.png")
    render(layers(P, only=("base", "post", "yaw_hub")), az=-70, el=24, w=900, h=900,
           target=(0, 0, 60), span=260, out="preview_base.png")
    sheet(P)
    print("wrote preview_*.png")
