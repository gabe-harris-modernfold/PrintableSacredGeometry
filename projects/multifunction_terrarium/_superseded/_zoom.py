import numpy as np, trimesh, pyvista as pv
m = trimesh.load("cascade_screen.stl")
def pvm(t):
    f = np.hstack([np.full((len(t.faces),1),3), t.faces]).ravel()
    return pv.PolyData(t.vertices, f)
pl = pv.Plotter(off_screen=True, window_size=(1600,800), shape=(1,2))
pl.subplot(0,0)
pl.add_mesh(pvm(m), color="#cfe9f7", smooth_shading=False, show_edges=False, specular=.7)
pl.add_text("lip close-up: are the notch teeth there?", font_size=11)
pl.camera_position=[(150,-40,120),(95,0,112),(0,0,1)]
pl.camera.zoom(3.2)
pl.subplot(0,1)
pl.add_mesh(pvm(m), color="#cfe9f7", show_edges=True, edge_color="#2a6d92",
            line_width=1, specular=.4)
pl.add_text("wireframe - one lip from above", font_size=11)
pl.camera_position=[(120,-25,175),(100,0,150),(0,0,1)]
pl.camera.zoom(4.5)
pl.screenshot("cascade_teeth.png")
print("ok")
