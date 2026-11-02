#!/usr/bin/env python3
"""Headless renders of the uterus candle holder, with a trimmed birthday candle in the internal
orifice and a flame.

Run from this directory, after cervix_candle.py and dish.py:
    python preview_cervix_candle.py
        cervix_candle_views.png      front, three-quarter, back, side and on its dish
        cervix_candle_vs_photo.png   the reference photo beside a candle-lit front render
        cervix_candle_canal_vs_photo.png  the canal and lips of both, enlarged
        cervix_candle_closeup.png    the canal and candle up close
"""

import numpy as np
import pyvista as pv
from PIL import Image

import cervix_candle as cc
import dish as dish_mod

ALABASTER = "#efe9df"

#          title           camera direction     distance
VIEWS = [("front",         (0.0, -1.0, 0.10),   300),
         ("three-quarter", (-0.7, -0.75, 0.35), 300),
         ("back",          (0.35, 1.0, 0.25),   300),
         ("side",          (1.0, 0.05, 0.15),   300),
         ("on its dish",   (-0.6, -0.8, 0.55),  420)]


def scene():
    u = cc.Uterus()
    dx, dy = dish_mod.centre_offset()
    body = pv.read("cervix_candle_upright.stl").translate((dx, dy, 0))
    foot = u.floor_contact()
    cy = dy + u.y_axis
    candle = pv.Cylinder(center=(dx, cy, 0.5 * (foot + u.z_wick)), direction=(0, 0, 1),
                         radius=cc.CANDLE_D / 2, height=u.z_wick - foot)
    flame_c = np.array([dx, cy, u.z_wick + 6.0])
    flame = pv.ParametricEllipsoid(2.2, 2.2, 6.0).translate(tuple(flame_c))
    dish = pv.read("dish.stl").translate((0, 0, -dish_mod.SOCKET_FLOOR_Z))
    return body, candle, flame, flame_c, dish


def views(out="cervix_candle_views.png"):
    body, candle, flame, _, dish = scene()
    pl = pv.Plotter(off_screen=True, shape=(1, len(VIEWS)), window_size=(520 * len(VIEWS), 640))
    c = np.array([0, 0, 50.0])
    for i, (title, d, dist) in enumerate(VIEWS):
        pl.subplot(0, i)
        pl.set_background("white")
        pl.add_mesh(body, color=ALABASTER, smooth_shading=True, specular=0.2, ambient=0.2)
        pl.add_mesh(candle, color="#ffffff", smooth_shading=True, ambient=0.35)
        pl.add_mesh(flame, color="#ffb030", smooth_shading=True, ambient=1.0)
        pl.add_mesh(dish, color="#c9c3b8", smooth_shading=True, ambient=0.2)
        pl.add_text(title, font_size=10, color="#333333")
        d = np.array(d, float) / np.linalg.norm(d)
        pl.camera.position = (c + d * dist).tolist()
        pl.camera.focal_point = c.tolist()
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30
    pl.screenshot(out)
    print("wrote", out)


def versus_photo(out="cervix_candle_vs_photo.png", canal_out="cervix_candle_canal_vs_photo.png"):
    """Front render lit like the photo (soft front fill plus a warm point light at the flame),
    framed to the photo's crop, beside the photo; and the canal region of both, enlarged."""
    body, candle, flame, flame_c, _ = scene()
    photo = Image.open("reference_photo.png").convert("RGB")
    H = 680
    W = round(photo.width * H / photo.height)
    pl = pv.Plotter(off_screen=True, window_size=(W, H), lighting="none")
    pl.set_background("#2a2520")
    pl.add_mesh(body, color=ALABASTER, smooth_shading=True, specular=0.1, ambient=0.18,
                diffuse=0.85)
    pl.add_mesh(candle, color="#fff6e8", smooth_shading=True, ambient=0.45)
    pl.add_mesh(flame, color="#ffd27a", ambient=1.0, diffuse=0.0)
    key = pv.Light(position=(260, -330, 260), focal_point=(0, 0, 50), color="#fff8ef",
                   intensity=0.7)
    fill = pv.Light(position=(-300, -300, 60), focal_point=(0, 0, 50), color="#f4efe8",
                    intensity=0.3)
    pl.add_light(fill)
    warm = pv.Light(position=tuple(flame_c), color="#ffb060", intensity=0.9, positional=True,
                    cone_angle=180, attenuation_values=(0.2, 0.02, 0.0))
    pl.add_light(key)
    pl.add_light(warm)
    # The photo spans about 110 mm top to bottom at the sculpture's scale (100 mm from the top
    # of the dome at y = 28 to the base at y = 342), so frame the same height.
    span = photo.height * HEIGHT_PER_PX()
    cz = 100.0 - (photo.height / 2 - 28) * HEIGHT_PER_PX()
    pl.camera.position = (0, -900, cz)
    pl.camera.focal_point = (0, 0, cz)
    pl.camera.up = (0, 0, 1)
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = span / 2
    render = Image.fromarray(pl.screenshot(None, return_img=True))
    sheet = Image.new("RGB", (2 * W + 20, H), (20, 20, 20))
    sheet.paste(photo.resize((W, H), Image.LANCZOS), (0, 0))
    sheet.paste(render, (W + 20, 0))
    sheet.save(out)
    print("wrote", out)
    # the canal and lips enlarged, photo beside render
    k = H / photo.height
    box = (120, 200, 240, photo.height)
    crops = [im.crop(tuple(round(v * k) for v in box)).resize(
             (round((box[2] - box[0]) * k * 3), round((box[3] - box[1]) * k * 3)), Image.LANCZOS)
             for im in (photo.resize((W, H), Image.LANCZOS), render)]
    pair = Image.new("RGB", (2 * crops[0].width + 10, crops[0].height), (20, 20, 20))
    pair.paste(crops[0], (0, 0))
    pair.paste(crops[1], (crops[0].width + 10, 0))
    pair.save(canal_out)
    print("wrote", canal_out)


def HEIGHT_PER_PX():
    import photo_trace as pt
    return cc.HEIGHT / (pt.BASE - 28)


def closeup(out="cervix_candle_closeup.png"):
    """The canal and candle from the front and from below-left, to check for ledges."""
    body, candle, flame, _, _ = scene()
    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1300, 800))
    for i, (d, dist) in enumerate([((0.0, -1.0, 0.0), 120), ((-0.55, -0.8, -0.25), 120)]):
        pl.subplot(0, i)
        pl.set_background("#2b2b2b")
        pl.add_mesh(body, color=ALABASTER, smooth_shading=True, specular=0.25, ambient=0.2)
        pl.add_mesh(candle, color="#ffffff", smooth_shading=True, ambient=0.4)
        pl.add_mesh(flame, color="#ffb030", ambient=1.0)
        c = np.array([0, 0, 26.0])
        d = np.array(d) / np.linalg.norm(d)
        pl.camera.position = (c + d * dist).tolist()
        pl.camera.focal_point = c.tolist()
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30
    pl.screenshot(out)
    print("wrote", out)


if __name__ == "__main__":
    views()
    versus_photo()
    closeup()
