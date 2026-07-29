"""Preview the collector pool: inner cups + gather-to-center horns."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml

SRC = "scrying_pool_collectors.stl"
full = trimesh.load(SRC)

ms = ml.MeshSet(); ms.load_new_mesh(SRC)
ms.apply_filter("meshing_decimation_quadric_edge_collapse",
                targetfacenum=130000, preservenormal=True)
ms.save_current_mesh("_preview_dec.stl")
dec = trimesh.load("_preview_dec.stl")
print("decimated faces:", len(dec.faces))

LIGHT = np.array([0.3, 0.4, 0.9]); LIGHT /= np.linalg.norm(LIGHT)


def shade(mesh, ax, elev, azim, bounds=None):
    inten = 0.4 + 0.6 * np.clip(mesh.face_normals @ LIGHT, 0, 1)
    cols = np.column_stack([0.44*inten, 0.62*inten, 0.82*inten, np.ones(len(inten))])
    pc = Poly3DCollection(mesh.vertices[mesh.faces], linewidths=0); pc.set_facecolor(cols)
    ax.add_collection3d(pc)
    lo, hi = mesh.bounds if bounds is None else bounds
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]))
    ax.view_init(elev=elev, azim=azim); ax.set_axis_off()


fig = plt.figure(figsize=(16, 12), facecolor="white")
ax1 = fig.add_subplot(2, 2, 1, projection="3d"); shade(dec, ax1, 82, 0)
ax1.set_title("Top — 21 inner cups + 21 horn throats", fontsize=11)
ax2 = fig.add_subplot(2, 2, 2, projection="3d"); shade(dec, ax2, 18, -52)
ax2.set_title("Iso — inner cups redirect to the near surface", fontsize=11)
ax3 = fig.add_subplot(2, 2, 3, projection="3d"); shade(dec, ax3, 3, -52)
ax3.set_title("Low side — horn mouths gather ambient sound", fontsize=11)

# close-up wedge viewed from outside
c = full.triangles_center
ang = np.arctan2(c[:, 1], c[:, 0])
keep = (np.abs(ang) < np.radians(34)) & (np.linalg.norm(c[:, :2], axis=1) > 250)
sub = full.submesh([keep], append=True)
ax4 = fig.add_subplot(2, 2, 4, projection="3d")
shade(sub, ax4, 6, 178, bounds=sub.bounds)
ax4.set_title("Close-up — horn mouth (outer) between two inner cups", fontsize=11)

fig.tight_layout()
fig.savefig("scrying_pool_collectors.png", dpi=118, bbox_inches="tight")
print("wrote scrying_pool_collectors.png")
