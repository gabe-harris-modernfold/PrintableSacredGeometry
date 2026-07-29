"""Preview render for the scrying pool: annotated cross-section + 3D cutaway."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import scrying_pool as sp

mesh, info = sp.build()
PHI, u, t = sp.PHI, info["u"], info["t"]
rb, rt = info["base_dia"] / 2, info["surf_dia"] / 2
zf, H = t, info["height"]

# full cross-section polygon (same construction as the solid)
rf, rim = info["rf"], info["rim"]
wall = sp.smoothstep_wall(rb, rt, zf + rf, H, 240)
outer_wall = wall + t * sp.outward_normals(wall)
a = np.linspace(-np.pi / 2, 0.0, sp.ARC)
fil = np.column_stack([(rb - rf) + rf * np.cos(a), (zf + rf) + rf * np.sin(a)])
inner = np.vstack([[0, zf], [rb - rf, zf], fil, wall])            # axis -> rim
outer = np.vstack([outer_wall[::-1], [rb + t, 0], [0, 0]])        # rim -> axis
px = np.r_[inner[:, 0], rt + rim, outer[:, 0]]
py = np.r_[inner[:, 1], H, outer[:, 1]]

fig = plt.figure(figsize=(15, 6.6), facecolor="white")

# ---- cross-section -------------------------------------------------------
ax = fig.add_subplot(1, 2, 1)
for sgn in (1, -1):
    ax.fill(sgn * px, py, color="#6ea8d8", alpha=0.35, lw=0)
    ax.plot(sgn * px, py, color="#25506f", lw=1.8)
ax.fill_between([-rt, rt], [H, H], [zf, zf], color="#bfe3f2", alpha=0.5, lw=0)  # water
ax.plot([-rt, rt], [H, H], color="#1a9ac0", lw=1.2, ls="--")

ax.annotate(f"water surface  Ø{2*rt:.1f} mm  =  u·φ⁴",
            (0, H), (0, H + 26), ha="center", color="#1a9ac0", fontsize=10,
            arrowprops=dict(arrowstyle="-", color="#1a9ac0", lw=.6))
ax.annotate(f"base  Ø{2*rb:.1f} mm  =  u·φ²",
            (0, zf), (0, -26), ha="center", color="#25506f", fontsize=10,
            arrowprops=dict(arrowstyle="-", color="#25506f", lw=.6))
ax.annotate("", (rt + 22, zf), (rt + 22, H), arrowprops=dict(arrowstyle="<->", color="#555"))
ax.text(rt + 30, (zf + H) / 2, "depth 100 mm = u", color="#555", fontsize=9,
        va="center", rotation=90)
ax.text(0, H / 2, f"wall t = 2φ = {t:.2f} mm", ha="center", color="#25506f", fontsize=9)
ax.set_title("Cross-section — golden ladder  depth : base : surface = 1 : φ² : φ⁴",
             fontsize=11, pad=14)
ax.set_aspect("equal"); ax.set_xlabel("mm"); ax.grid(alpha=0.25)
ax.set_ylim(-46, H + 46)

# ---- 3D half cutaway -----------------------------------------------------
half = mesh.slice_plane([0, 0, 0], [0, 1, 0], cap=True)
ax3 = fig.add_subplot(1, 2, 2, projection="3d")
pc = Poly3DCollection(half.vertices[half.faces], facecolor="#7fb5de",
                      edgecolor="#4c7fa6", alpha=0.95, linewidths=0.05)
ax3.add_collection3d(pc)
lo, hi = half.bounds
w = hi[0] - lo[0]
ax3.set_xlim(lo[0], hi[0]); ax3.set_ylim(-w/2, w/2); ax3.set_zlim(0, w)
ax3.set_box_aspect((1, 1, (hi[2]-lo[2]) / w))
ax3.view_init(elev=18, azim=-60)
ax3.set_title(f"Half cutaway — Ø{info['outer_dia']/10:.1f} cm × {info['height']/10:.2f} cm"
              f"  ·  holds {info['capacity_L']:.1f} L", fontsize=11)
ax3.set_axis_off()

fig.tight_layout()
fig.savefig("scrying_pool_preview.png", dpi=115, bbox_inches="tight")
print("wrote scrying_pool_preview.png")
