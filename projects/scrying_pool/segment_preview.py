"""MeshLab QA on every segment + a top-view print-layout preview."""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import trimesh
import pymeshlab as ml

BED = 310.0
files = sorted(glob.glob("scrying_pool_segments/*.stl"))

# ---- MeshLab QA ----------------------------------------------------------
print(f"{'piece':16s} {'holes':>5s} {'nonmanifE':>9s} {'comp':>4s} {'genus':>5s} {'wt':>3s}")
for f in files:
    ms = ml.MeshSet(); ms.load_new_mesh(f)
    tm = ms.get_topological_measures()
    print(f"{os.path.basename(f):16s} {tm.get('number_holes','?'):>5} "
          f"{tm.get('non_two_manifold_edges','?'):>9} "
          f"{tm.get('connected_components_number','?'):>4} "
          f"{tm.get('genus','?'):>5} "
          f"{str(tm.get('is_mesh_watertight','?'))[0]:>3}")

# ---- layout preview ------------------------------------------------------
n = len(files)
cols = 4
rows = int(np.ceil(n / cols))
fig, axes = plt.subplots(rows, cols, figsize=(3.1 * cols, 3.1 * rows),
                         facecolor="white")
axes = np.atleast_1d(axes).ravel()
for ax, f in zip(axes, files):
    m = trimesh.load(f)
    v = m.vertices
    # top-view outline: project onto xy, draw the convex-ish footprint via edges
    tri = v[m.faces][:, :, :2]
    ax.add_collection(plt.matplotlib.collections.PolyCollection(
        tri, facecolors="#7fb5de", edgecolors="none", alpha=0.15))
    lo, hi = m.bounds
    cx, cy = (lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2
    ax.add_patch(Rectangle((cx - BED/2, cy - BED/2), BED, BED, fill=False,
                           ec="#cc5555", ls="--", lw=1))
    fx, fy = hi[0]-lo[0], hi[1]-lo[1]
    ax.set_title(f"{os.path.basename(f)[:-4]}\n{fx:.0f} × {fy:.0f} mm",
                 fontsize=9)
    ax.set_aspect("equal"); ax.set_xlim(cx-BED/2-15, cx+BED/2+15)
    ax.set_ylim(cy-BED/2-15, cy+BED/2+15); ax.set_xticks([]); ax.set_yticks([])
for ax in axes[n:]:
    ax.axis("off")
fig.suptitle(f"Print layout — {n} pieces, each in a {int(BED)} mm bed (dashed)",
             fontsize=12)
fig.tight_layout()
fig.savefig("scrying_pool_segments.png", dpi=110, bbox_inches="tight")
print("\nwrote scrying_pool_segments.png")
