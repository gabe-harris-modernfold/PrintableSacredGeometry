"""Render preview PNGs from nickel_box_handle.stl (OpenSCAD's own --render of the
Minkowski-rounded part takes ~10 min per image; this takes seconds).

Run from this directory after regenerating the STL:
    python render_previews.py
"""
import numpy as np, trimesh, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# print-orientation STL: print (x, y, z) = (-Z_cap, Y_cap, X_cap); X_cap runs along the box
part = trimesh.load("nickel_box_handle.stl")
to_cap = trimesh.transformations.rotation_matrix(np.radians(90), [0, 1, 0])   # print -> cap frame
cap = part.copy(); cap.apply_transform(to_cap)

LIGHT = np.array([0.3, -0.5, 0.8]); LIGHT /= np.linalg.norm(LIGHT)

def colours(mesh, base):
    shade = 0.35 + 0.65 * np.clip(mesh.face_normals @ LIGHT, 0, 1)
    return np.c_[np.outer(shade, base), np.ones_like(shade)]

def iso(meshes, bases, elev, azim, fname, title):
    fig = plt.figure(figsize=(10, 8)); ax = fig.add_subplot(111, projection="3d")
    allv = np.vstack([m.vertices for m in meshes])
    for m, b in zip(meshes, bases):
        ax.add_collection3d(Poly3DCollection(m.vertices[m.faces], facecolors=colours(m, b), edgecolors="none"))
    c = (allv.max(0) + allv.min(0)) / 2; r = (allv.max(0) - allv.min(0)).max() / 2
    ax.set_xlim(c[0]-r, c[0]+r); ax.set_ylim(c[1]-r, c[1]+r); ax.set_zlim(c[2]-r, c[2]+r)
    ax.set_box_aspect((1, 1, 1)); ax.set_proj_type("ortho"); ax.view_init(elev, azim); ax.set_axis_off()
    ax.set_title(title); fig.tight_layout(); fig.savefig(fname, dpi=110); plt.close(fig)

def end_view(mesh, fname, title):
    """Orthographic view along -X_cap (looking at the outside end face). Painter's sort on X."""
    tri = mesh.vertices[mesh.faces]
    order = np.argsort(-tri[:, :, 0].mean(1))          # far (large x) first
    cols = colours(mesh, [0.85, 0.75, 0.35])[order]
    poly = tri[order][:, :, [1, 2]]                    # (Y, Z) plane
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.add_collection(PolyCollection(poly, facecolors=cols, edgecolors="none"))
    ax.set_aspect("equal"); ax.autoscale(); ax.margins(0.05); ax.set_axis_off(); ax.set_title(title)
    fig.tight_layout(); fig.savefig(fname, dpi=110); plt.close(fig)

if __name__ == "__main__":
    end_view(cap, "preview_end.png", "End view - grip with finger grooves, paracord ears at the lower corners")
    iso([part], [[0.85, 0.75, 0.35]], 35, -55, "preview_part.png", "Print orientation (end plate on the bed)")
    BOX_L, BOX_W, BOX_H = 240, 88, 100
    box = trimesh.creation.box([BOX_L, BOX_W, BOX_H]); box.apply_translation([BOX_L/2, 0, 0])
    cap_a = cap.copy(); cap_a.apply_translation([-4, 0, 0])
    cap_b = cap.copy(); cap_b.apply_transform(trimesh.transformations.reflection_matrix([0, 0, 0], [1, 0, 0]))
    cap_b.apply_translation([BOX_L + 4, 0, 0])
    iso([box, cap_a, cap_b], [[0.87, 0.72, 0.53], [0.27, 0.51, 0.71], [0.27, 0.51, 0.71]], 28, -60,
        "preview_assembly.png", "Box of nickels with a handle cap on each end")
    print("previews written")
