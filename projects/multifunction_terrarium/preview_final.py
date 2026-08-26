"""Assembled previews of the whole terrarium."""

import numpy as np, trimesh, pyvista as pv
import params as P

pv.global_theme.allow_empty_mesh = True
GLASS = "#d8ecf6"; WATER = "#3f9ad0"; DEEP = "#22648f"
MOSS = "#5d9a67"; STONE = "#8d99a1"; POOL = "#54b4e0"


def pvm(t):
    f = np.hstack([np.full((len(t.faces), 1), 3), t.faces]).ravel()
    return pv.PolyData(t.vertices.copy(), f)


def load(n, dz=0.0):
    m = trimesh.load(n)
    if dz:
        m.apply_translation([0, 0, dz])
    return m


shells = [load(f"vessel_shell_{k}.stl") for k in range(3)]
res = load("vessel_reservoir.stl")
bed = load("vessel_bedtray.stl")
casc = load("cascade_screen.stl", 215.0)
gut = load("drip_gutters.stl")
spl = load("drip_splitter.stl")

water = trimesh.creation.cylinder(radius=143, height=P.Z_RES, sections=64)
water.apply_translation([0, 0, P.Z_RES / 2])

BODY = [(res, STONE), (water, POOL), (bed, MOSS),
        (casc, WATER), (gut, WATER), (spl, DEEP)]


def scene(pl, cut=None, glass_op=0.15):
    for s in shells:
        g = s.slice_plane([0, 0, 0], cut) if cut is not None else s
        if g is None or len(g.faces) == 0:
            continue
        pl.add_mesh(pvm(g), color=GLASS,
                    opacity=0.6 if cut is not None else glass_op,
                    smooth_shading=True, specular=0.9, specular_power=32)
    for m, c in BODY:
        g = m.slice_plane([0, 0, 0], cut) if cut is not None else m
        if g is None or len(g.faces) == 0:
            continue
        pl.add_mesh(pvm(g), color=c, opacity=0.5 if c == POOL else 1.0,
                    smooth_shading=True, specular=0.6, specular_power=22)


CEN = (0, 0, 280)

pl = pv.Plotter(off_screen=True, window_size=(1900, 1250), shape=(2, 2))

pl.subplot(0, 0)
scene(pl)
pl.add_text("ASSEMBLED   304 across corners x 540 tall", font_size=12)
pl.camera_position = [(1500, -1350, 1150), CEN, (0, 0, 1)]
pl.camera.zoom(1.45)

pl.subplot(0, 1)
scene(pl, cut=[0, -1, 0])
pl.add_text("SECTION   reservoir / bed + siphon / cascade / drip spirals", font_size=12)
pl.camera.SetParallelProjection(True)
pl.camera_position = [(0, -2400, 300), CEN, (0, 0, 1)]
pl.camera.parallel_scale = 320

pl.subplot(1, 0)
scene(pl, glass_op=0.09)
pl.add_text("STANDING EYE   1550 mm, 700 mm back, object on a 750 table", font_size=12)
pl.camera_position = [(700, 0, 800), (0, 0, 270), (0, 0, 1)]
pl.camera.zoom(1.05)

pl.subplot(1, 1)
scene(pl, glass_op=0.12)
pl.add_text("CROWN   condenser tap, 16 spirals, ovoid head", font_size=12)
pl.camera_position = [(620, -560, 700), (0, 0, 440), (0, 0, 1)]
pl.camera.zoom(1.5)

pl.screenshot("terrarium_final.png")
print("wrote terrarium_final.png")
