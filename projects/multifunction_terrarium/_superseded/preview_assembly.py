"""Stack what exists so far: 3 shells + 3 cascade screens."""
import numpy as np, trimesh, pyvista as pv
import params as P

def pvm(t):
    f = np.hstack([np.full((len(t.faces),1),3), t.faces]).ravel()
    return pv.PolyData(t.vertices.copy(), f)

shells, screens = [], []
scr = trimesh.load("cascade_screen.stl")
for k in range(P.N_MODULE):
    s = trimesh.load(f"shell_module_{k}.stl"); s.apply_translation([0,0,k*P.MOD_H])
    shells.append(s)
    c = scr.copy(); c.apply_translation([0,0,k*P.MOD_H]); screens.append(c)

pl = pv.Plotter(off_screen=True, window_size=(1500,1000), shape=(1,2))
pl.subplot(0,0)
for s in shells:  pl.add_mesh(pvm(s), color="#eaf6fc", opacity=0.30, smooth_shading=True)
for c in screens: pl.add_mesh(pvm(c), color="#7fc6e8", smooth_shading=True, specular=.6)
pl.add_text("assembled: shell (ghosted) + cascade", font_size=11)
pl.camera_position=[(1500,-1350,1250),(0,0,315),(0,0,1)]
pl.camera.zoom(1.25)

pl.subplot(0,1)
for s in shells:
    h = s.slice_plane([0,0,0],[0,-1,0])
    pl.add_mesh(pvm(h), color="#d8eef8", smooth_shading=True, specular=.5)
for c in screens:
    h = c.slice_plane([0,0,0],[0,-1,0])
    pl.add_mesh(pvm(h), color="#3d9fd0", smooth_shading=True, specular=.6)
pl.add_text("cut away - cascade inside, ports aimed", font_size=11)
pl.camera_position=[(1450,-1100,900),(0,0,315),(0,0,1)]
pl.camera.zoom(1.25)
pl.screenshot("assembly.png")
print("wrote assembly.png")
