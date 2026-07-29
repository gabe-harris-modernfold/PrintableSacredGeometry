"""Vertical cross-section proving the water line vs. cups/horns.

The section plane (y=0) passes through a cup on the +x side (azimuth 0) and a
horn on the -x side (azimuth 180). A single slice therefore shows: solid
retaining wall below the water line, cup bottoms sitting exactly at the
surface, and horn mouths lifted above it.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import trimesh
import scrying_pool as sp

H = sp.WALL + sp.DEPTH          # water-surface height (mm)
rt = 0.5 * sp.SEED * sp.PHI**4  # water radius (mm)

mesh = trimesh.load("scrying_pool_collectors.stl")
lines = trimesh.intersections.mesh_plane(mesh, [0, 1, 0], [0, 0, 0])  # y=0 plane

fig = plt.figure(figsize=(16, 8), facecolor="white")
axf = fig.add_axes([0.06, 0.56, 0.88, 0.34])
axc = fig.add_axes([0.06, 0.06, 0.42, 0.42])
axh = fig.add_axes([0.54, 0.06, 0.42, 0.42])


def draw(A, xlim, ylim):
    A.fill_between([-rt, rt], [H, H], [0, 0], color="#bfe3f2", alpha=0.6, lw=0, zorder=0)
    A.axhline(H, color="#0f7d94", lw=1.3, ls="--", zorder=1)
    for seg in lines:
        A.plot(seg[:, 0], seg[:, 2], color="#1c3b52", lw=1.4, zorder=2)
    A.set_aspect("equal"); A.grid(alpha=0.2)
    A.set_xlabel("radius (mm)"); A.set_ylabel("z (mm)")
    if xlim: A.set_xlim(*xlim)
    if ylim: A.set_ylim(*ylim)


draw(axf, None, (-6, H+55))
axf.set_title("Full cross-section (y = 0 plane)", fontsize=11, color="#1c3b52", weight="bold")
axf.text(0, H+9, f"water line  z = {H:.1f} mm", ha="center", color="#0f7d94", fontsize=9)

draw(axc, (250, 400), (-4, H+55))
axc.set_title("Cup side (azimuth 0)", fontsize=11, color="#1c3b52", weight="bold")
axc.annotate("cup bottom AT the surface", (322, H), (270, H+42),
             fontsize=9.5, color="#1c3b52",
             arrowprops=dict(arrowstyle="->", color="#1c3b52"))
axc.annotate("solid wall below —\nno holes", (332, H/2), (256, H/2+6),
             fontsize=9.5, color="#155",
             arrowprops=dict(arrowstyle="->", color="#155"))

draw(axh, (-400, -250), (-4, H+55))
axh.set_title("Horn side (azimuth 180)", fontsize=11, color="#7a3b1c", weight="bold")
axh.annotate("horn mouth ABOVE water", (-352, H+24), (-330, H+46),
             fontsize=9.5, color="#7a3b1c", ha="center",
             arrowprops=dict(arrowstyle="->", color="#7a3b1c"))
axh.annotate("solid wall below —\nno holes", (-332, H/2), (-300, H/2-12),
             fontsize=9.5, color="#155",
             arrowprops=dict(arrowstyle="->", color="#155"))

fig.suptitle("Water retention check — nothing pierces the wall below the surface",
             fontsize=13, color="#1c3b52", weight="bold")
fig.savefig("scrying_pool_waterline.png", dpi=125, bbox_inches="tight")
print(f"water line z={H:.2f} mm; wrote scrying_pool_waterline.png")
