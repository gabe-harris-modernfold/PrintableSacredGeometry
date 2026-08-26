"""Orthographic-ish view sheet + cutaways of what is built so far."""
import numpy as np, trimesh, pyvista as pv
import params as P

pv.global_theme.allow_empty_mesh = True
SHELL = "#e3f1f9"; SCREEN = "#5cb4de"; CUT = "#f2b46a"

def pvm(t):
    f = np.hstack([np.full((len(t.faces),1),3), t.faces]).ravel()
    return pv.PolyData(t.vertices.copy(), f)

scr0 = trimesh.load("cascade_screen.stl")
shells, screens = [], []
for k in range(P.N_MODULE):
    s = trimesh.load(f"shell_module_{k}.stl"); s.apply_translation([0,0,k*P.MOD_H]); shells.append(s)
    c = scr0.copy(); c.apply_translation([0,0,k*P.MOD_H]); screens.append(c)
ALL = shells + screens

def add(pl, meshes, plane=None, keep=None, shell_op=1.0):
    for m, col in [(s, SHELL) for s in shells] + [(c, SCREEN) for c in screens]:
        if meshes == "screens" and col == SHELL: continue
        g = m
        if plane is not None:
            g = m.slice_plane(plane, keep)
            if g is None or len(g.faces) == 0: continue
        pl.add_mesh(pvm(g), color=col, smooth_shading=True, specular=.55,
                    specular_power=22, opacity=shell_op if col == SHELL else 1.0)

CEN = (0,0,P.TOTAL_H/2)

# ---------------------------------------------------------------- sheet A
pl = pv.Plotter(off_screen=True, window_size=(1800,1250), shape=(2,2))

pl.subplot(0,0)                                   # elevation, face on
add(pl, "all")
pl.add_text("ELEVATION - face on", font_size=12)
pl.camera.SetParallelProjection(True)
pl.camera_position=[(2600,0,315),CEN,(0,0,1)]; pl.camera.parallel_scale=360

pl.subplot(0,1)                                   # plan
add(pl, "all")
pl.add_text("PLAN - from above (310 across corners)", font_size=12)
pl.camera.SetParallelProjection(True)
pl.camera_position=[(0,0,2600),(0,0,315),(0,1,0)]; pl.camera.parallel_scale=190

pl.subplot(1,0)                                   # half section
add(pl, "all", plane=[0,0,0], keep=[0,-1,0])
pl.add_text("HALF SECTION - shell / 30mm theatre / cascade", font_size=12)
pl.camera.SetParallelProjection(True)
pl.camera_position=[(0,-2600,315),CEN,(0,0,1)]; pl.camera.parallel_scale=360

pl.subplot(1,1)                                   # iso cutaway
add(pl, "all", plane=[0,0,0], keep=[-1,0,0])
pl.add_text("ISO CUTAWAY", font_size=12)
pl.camera_position=[(1500,-1300,1150),CEN,(0,0,1)]; pl.camera.zoom(1.3)
pl.screenshot("views_A.png")

# ---------------------------------------------------------------- sheet B
pl = pv.Plotter(off_screen=True, window_size=(1800,1250), shape=(2,2))

pl.subplot(0,0)                                   # what a standing viewer sees
add(pl, "all", shell_op=0.85)
pl.add_text("STANDING EYE  (1550mm, 700mm back, object on a 750 table)", font_size=12)
pl.camera_position=[(700,0,800),(0,0,315),(0,0,1)]

pl.subplot(0,1)                                   # port band detail
add(pl, "all")
pl.add_text("PORT BAND - aimed aperture + theatre + terraces behind", font_size=12)
pl.camera_position=[(560,-120,470),(140,0,430),(0,0,1)]

pl.subplot(1,0)                                   # plan section through a band
add(pl, "all", plane=[0,0,455], keep=[0,0,-1])
pl.add_text("PLAN SECTION at z=455 - hex shell, gap, hex cascade", font_size=12)
pl.camera.SetParallelProjection(True)
pl.camera_position=[(0,0,-2600),(0,0,455),(0,1,0)]; pl.camera.parallel_scale=175

pl.subplot(1,1)                                   # terrace detail
add(pl, "screens")
pl.add_text("TERRACE DETAIL - 10mm risers, 2.86mm treads, notched lips", font_size=12)
pl.camera_position=[(300,-90,330),(95,0,300),(0,0,1)]; pl.camera.zoom(2.4)
pl.screenshot("views_B.png")
print("wrote views_A.png views_B.png")
