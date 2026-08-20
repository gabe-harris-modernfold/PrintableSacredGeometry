"""Final hero render of the completed golden scrying pool with collectors."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml

import scrying_pool as sp
import add_collectors as ac

# geometry facts for the water disc + captions
_, info = ac.build()
rt = 0.5 * info["surf_dia"]; H = info["height"]

# decimate for a light but detailed render
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse",
                targetfacenum=160000, preservenormal=True)
ms.save_current_mesh("_hero_dec.stl")
m = trimesh.load("_hero_dec.stl")
print("render faces:", len(m.faces))

# two-light shading
L1 = np.array([0.3, 0.5, 0.9]); L1 /= np.linalg.norm(L1)
L2 = np.array([-0.6, -0.3, 0.4]); L2 /= np.linalg.norm(L2)
n = m.face_normals
inten = 0.30 + 0.62*np.clip(n@L1, 0, 1) + 0.20*np.clip(n@L2, 0, 1)
inten = np.clip(inten, 0, 1)
base = np.array([0.42, 0.60, 0.80])
cols = np.column_stack([np.clip(base[0]*inten,0,1), np.clip(base[1]*inten,0,1),
                        np.clip(base[2]*inten,0,1), np.ones(len(inten))])


def water(ax):
    a = np.linspace(0, 2*np.pi, 160)
    ring = np.column_stack([rt*np.cos(a), rt*np.sin(a), np.full_like(a, H)])
    disc = Poly3DCollection([ring], facecolor=(0.20, 0.68, 0.80, 0.45),
                            edgecolor=(0.12, 0.5, 0.62, 0.7), linewidths=1.2)
    ax.add_collection3d(disc)


def draw(ax, elev, azim):
    pc = Poly3DCollection(m.vertices[m.faces], linewidths=0); pc.set_facecolor(cols)
    ax.add_collection3d(pc)
    water(ax)
    lo, hi = m.bounds
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(0, hi[0]-lo[0])
    ax.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[0]-lo[0]))
    ax.view_init(elev=elev, azim=azim); ax.set_axis_off()


fig = plt.figure(figsize=(16, 9), facecolor="white")
axh = fig.add_axes([0.0, 0.0, 0.72, 1.0], projection="3d"); draw(axh, 24, 38)
axt = fig.add_axes([0.70, 0.30, 0.30, 0.40], projection="3d"); draw(axt, 84, 0)

fig.text(0.36, 0.95, "Golden Scrying Pool — φ⁴ : φ² : 1  ·  Ø70 cm  ·  18.6 L",
         ha="center", fontsize=16, weight="bold", color="#1c3b52")
fig.text(0.36, 0.075,
         "perfect circle · smooth φ flare · 21 concave surface cups + 21 "
         "gather-to-center horns · watertight, single body",
         ha="center", fontsize=11.5, color="#33566e")
fig.savefig("scrying_pool_final.png", dpi=130, bbox_inches="tight",
            facecolor="white")
print("wrote scrying_pool_final.png")
