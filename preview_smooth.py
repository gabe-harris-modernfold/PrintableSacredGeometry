"""Confirm the smooth bowl floor (no ring) — section + interior cutaway."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml
import scrying_pool as sp

BG, PANEL, LINE, INK, CY, TEXT, DIM = ("#0c2c4d", "#0a2540", "#9cc2e6",
    "#eaf2ff", "#6fd0e8", "#d4e2f4", "#7f9bbd")
plt.rcParams["font.family"] = "sans-serif"
H = sp.WALL + sp.DEPTH
rt = 0.5*sp.SEED*sp.PHI**4

mesh = trimesh.load("scrying_pool_collectors.stl")
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse", targetfacenum=150000, preservenormal=True)
ms.save_current_mesh("_smooth_dec.stl")
dec = trimesh.load("_smooth_dec.stl")

fig = plt.figure(figsize=(16, 6.6), facecolor=BG)

ax = fig.add_axes([0.05, 0.12, 0.44, 0.78]); ax.set_facecolor(PANEL)
seg = trimesh.intersections.mesh_plane(mesh, [0, 1, 0], [0, 0, 0])
ax.fill_between([-rt, rt], [H, H], [0, 0], color="#12466e", alpha=0.55, lw=0)
ax.axhline(H, color=CY, lw=1.1, ls="--")
for s in seg: ax.plot(s[:, 0], s[:, 2], color=INK, lw=1.3)
ax.set_xlim(-360, 360); ax.set_ylim(-3, 62); ax.set_aspect("equal")
for sp_ in ax.spines.values(): sp_.set_color(LINE)
ax.tick_params(colors=TEXT, labelsize=8)
ax.set_title("SECTION — floor blends tangentially into the wall (no ring)",
             color=CY, fontsize=11, weight="bold")
ax.set_xlabel("radius (mm)", color=TEXT); ax.set_ylabel("z (mm)", color=TEXT)
ax.annotate("flat centre\n(uniform depth)", (0, sp.WALL), (-120, 40), color=DIM, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=DIM))
ax.annotate("smooth tangential blend", (175, 10), (40, 45), color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=INK))

ax3 = fig.add_axes([0.50, 0.05, 0.48, 0.9], projection="3d"); ax3.patch.set_alpha(0)
half = dec.slice_plane([0, 0, 0], [0, 1, 0], cap=True)
L = np.array([0.3, 0.55, 0.85]); L /= np.linalg.norm(L)
it = np.clip(0.32 + 0.68*np.clip(half.face_normals@L, 0, 1), 0, 1)
cols = np.column_stack([0.6*it, 0.76*it, 0.9*it, np.ones(len(it))])
ax3.add_collection3d(Poly3DCollection(half.vertices[half.faces], facecolor=cols, lw=0))
lo, hi = half.bounds
ax3.set_xlim(lo[0], hi[0]); ax3.set_ylim(-(hi[0]-lo[0])/2, (hi[0]-lo[0])/2); ax3.set_zlim(0, hi[0]-lo[0])
ax3.set_box_aspect((hi[0]-lo[0], hi[0]-lo[0], hi[0]-lo[0]))
ax3.view_init(elev=26, azim=-58); ax3.set_axis_off()
ax3.set_title("Interior cutaway — continuous smooth bowl", color=CY, fontsize=11, weight="bold")

fig.savefig("scrying_pool_smoothfloor.png", dpi=125, facecolor=BG, bbox_inches="tight")
print("wrote scrying_pool_smoothfloor.png")
