"""Preview the pool on its domed collector stand."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml
import scrying_pool as sp

BG, PANEL, LINE, INK, CY, TEXT, DIM, WARN = ("#0c2c4d", "#0a2540", "#9cc2e6",
    "#eaf2ff", "#6fd0e8", "#d4e2f4", "#7f9bbd", "#e0a24b")
plt.rcParams["font.family"] = "sans-serif"

Hstand = sp.SEED/sp.PHI
Hpool = sp.WALL + sp.DEPTH
water_z = Hstand + Hpool
rt = 0.5*sp.SEED*sp.PHI**4

mesh = trimesh.load("scrying_pool_stand.stl")
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_stand.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse", targetfacenum=140000, preservenormal=True)
ms.save_current_mesh("_stand_dec.stl")
dec = trimesh.load("_stand_dec.stl")

fig = plt.figure(figsize=(16, 7.4), facecolor=BG)

# ---- center cross-section ------------------------------------------------
ax = fig.add_axes([0.05, 0.08, 0.42, 0.84]); ax.set_aspect("equal")
ax.set_facecolor(PANEL)
lines = trimesh.intersections.mesh_plane(mesh, [0, 1, 0], [0, 0, 0])
ax.fill_between([-rt, rt], [water_z, water_z], [Hstand, Hstand], color="#12466e", alpha=0.55, lw=0)
ax.axhline(water_z, color=CY, lw=1.1, ls="--")
for seg in lines:
    ax.plot(seg[:, 0], seg[:, 2], color=INK, lw=1.2)
ax.set_xlim(-210, 210); ax.set_ylim(-8, water_z+40)
for s in ax.spines.values(): s.set_color(LINE)
ax.tick_params(colors=TEXT, labelsize=8)
ax.annotate("water surface", (150, water_z), (60, water_z+28), color=CY, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=CY))
ax.annotate("hollow dome — convex up,\nopening down (inverted cup)", (44, Hstand*0.5),
            (-205, Hstand*0.72), color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=INK))
ax.annotate("apex seat → pool centre", (0, Hstand), (-200, Hstand+34), color=WARN, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=WARN))
ax.annotate("", (200, 0), (200, Hstand), arrowprops=dict(arrowstyle="<->", color=DIM))
ax.text(186, Hstand/2, f"stand H = u/φ = {Hstand:.1f} mm", rotation=90, va="center",
        ha="right", color=DIM, fontsize=8.5)
ax.set_title("SECTION ON y=0 — pool on the domed collector stand",
             color=CY, fontsize=11, weight="bold")
ax.set_xlabel("radius (mm)", color=TEXT); ax.set_ylabel("z (mm)", color=TEXT)

# ---- iso half-cutaway ----------------------------------------------------
ax3 = fig.add_axes([0.50, 0.04, 0.48, 0.90], projection="3d"); ax3.patch.set_alpha(0)
half = dec.slice_plane([0, 0, 0], [0, 1, 0], cap=True)
L = np.array([0.35, 0.55, 0.85]); L /= np.linalg.norm(L)
it = np.clip(0.3 + 0.7*np.clip(half.face_normals@L, 0, 1), 0, 1)
cols = np.column_stack([0.62*it, 0.78*it, 0.92*it, np.ones(len(it))])
ax3.add_collection3d(Poly3DCollection(half.vertices[half.faces], facecolor=cols,
                     lw=0.05, edgecolor=(0.4, 0.55, 0.72, 0.25)))
lo, hi = half.bounds
ax3.set_xlim(lo[0], hi[0]); ax3.set_ylim(-(hi[0]-lo[0])/2, (hi[0]-lo[0])/2)
ax3.set_zlim(0, hi[0]-lo[0])
ax3.set_box_aspect((hi[0]-lo[0], hi[0]-lo[0], hi[0]-lo[0]))
ax3.view_init(elev=14, azim=-62); ax3.set_axis_off()
ax3.set_title("Half cutaway — pool + dome pedestal", color=CY, fontsize=11, weight="bold")

fig.savefig("scrying_pool_stand.png", dpi=125, facecolor=BG, bbox_inches="tight")
print("wrote scrying_pool_stand.png")
