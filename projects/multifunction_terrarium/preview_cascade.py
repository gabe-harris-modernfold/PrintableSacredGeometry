"""Iso + radial section of one cascade module, so the staircase can be checked by eye."""
import numpy as np, trimesh, pyvista as pv
import params as P

pv.global_theme.allow_empty_mesh = True
m = trimesh.load("cascade_screen.stl")

def pvm(tm):
    f = np.hstack([np.full((len(tm.faces), 1), 3), tm.faces]).ravel()
    return pv.PolyData(tm.vertices, f)

# --- radial wedge, to read the profile
wedge = m.slice_plane([0, 0, 0], [0, 1, 0]).slice_plane([0, 0, 0], [-1, 0, 0])

pl = pv.Plotter(off_screen=True, window_size=(1700, 850), shape=(1, 2))
pl.subplot(0, 0)
pl.add_mesh(pvm(m), color="#bfe6f5", smooth_shading=True, specular=0.6,
            specular_power=25)
pl.add_text("cascade screen - one module", font_size=11)
pl.camera_position = [(340, -300, 300), (0, 0, 105), (0, 0, 1)]

pl.subplot(0, 1)
pl.add_mesh(pvm(wedge), color="#bfe6f5", smooth_shading=True, specular=0.6)
pl.add_text("quarter cut - 21 terraces, lip over next tread", font_size=11)
pl.camera_position = [(330, -330, 230), (55, -55, 105), (0, 0, 1)]
pl.screenshot("cascade_overview.png")

# --- true 2D profile, the thing that actually has to be right
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sec = m.section(plane_origin=[0, 0, 0], plane_normal=[0, 1, 0])
sl, _ = sec.to_2D()
fig, ax = plt.subplots(figsize=(6, 9))
for ent in sl.entities:
    p = sl.vertices[ent.points]
    ax.plot(p[:, 0], p[:, 1], lw=0.7, color="#0b6ea8")
ax.set_aspect("equal"); ax.set_title("cascade profile (section)  step 10, tread 2.86")
ax.set_xlabel("mm"); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig("cascade_profile.png", dpi=130)
print("wrote cascade_overview.png, cascade_profile.png")
