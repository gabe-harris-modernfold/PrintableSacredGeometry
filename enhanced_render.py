"""Annotated info-graphic of the completed golden scrying pool."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
import trimesh, pymeshlab as ml

import scrying_pool as sp
import add_collectors as ac

INK, MID, ACC = "#1c3b52", "#33566e", "#0f7d94"
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

# ---- geometry facts ------------------------------------------------------
_, info = ac.build()
PHI = sp.PHI
rt = 0.5 * info["surf_dia"]; rb = 0.5 * info["base_dia"]; H = info["height"]
t = info["t"]; rf = info["rf"]; rim = info["rim"]; d = info["depth"]
cap = info["capacity_L"]
Rm = rt + ac.COLLAR_DR; A = ac.COLLAR_A; Z = ac.COLLAR_Z; HZ = ac.HORN_Z


def ex(x):
    return int(round(np.log(x) / np.log(PHI)))


nb, ns, nd = ex(info["base_dia"]/d), ex(info["surf_dia"]/d), ex(sp.SEED/d)
ladder = f"1 : φ{str(nb).translate(SUP)} : φ{str(ns).translate(SUP)}"
depth_lbl = f"u/φ{str(nd).translate(SUP)}"

# ---- render mesh (decimated) ---------------------------------------------
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse",
                targetfacenum=150000, preservenormal=True)
ms.save_current_mesh("_enh_dec.stl")
m = trimesh.load("_enh_dec.stl")

L1 = np.array([0.3, 0.5, 0.9]); L1 /= np.linalg.norm(L1)
L2 = np.array([-0.6, -0.3, 0.4]); L2 /= np.linalg.norm(L2)
nrm = m.face_normals
it = np.clip(0.30 + 0.62*np.clip(nrm@L1, 0, 1) + 0.20*np.clip(nrm@L2, 0, 1), 0, 1)
cols = np.column_stack([0.42*it, 0.60*it, 0.80*it, np.ones(len(it))])

fig = plt.figure(figsize=(17.5, 10), facecolor="white")
axh = fig.add_axes([-0.02, 0.06, 0.62, 0.86], projection="3d")
ELEV, AZIM = 22, 38
pc = Poly3DCollection(m.vertices[m.faces], linewidths=0); pc.set_facecolor(cols)
axh.add_collection3d(pc)
a = np.linspace(0, 2*np.pi, 200)
axh.add_collection3d(Poly3DCollection(
    [np.column_stack([rt*np.cos(a), rt*np.sin(a), np.full_like(a, H)])],
    facecolor=(0.20, 0.68, 0.80, 0.42), edgecolor=(0.1, 0.5, 0.6, 0.7), lw=1.2))
lo, hi = m.bounds
axh.set_xlim(lo[0], hi[0]); axh.set_ylim(lo[1], hi[1]); axh.set_zlim(0, hi[0]-lo[0])
axh.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[0]-lo[0]))
axh.view_init(elev=ELEV, azim=AZIM); axh.set_axis_off()
fig.canvas.draw()


def call(feat, frac, text, ha="left"):
    x2, y2, _ = proj3d.proj_transform(*feat, axh.get_proj())
    axh.annotate(text, xy=(x2, y2), xycoords=axh.transData,
                 xytext=frac, textcoords="axes fraction",
                 fontsize=10.5, color=INK, ha=ha, va="center", weight="bold",
                 arrowprops=dict(arrowstyle="-", color=ACC, lw=1.1,
                                 connectionstyle="arc3,rad=0.15"),
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=ACC, lw=0.8))


def pol(r, deg, z):
    return np.array([r*np.cos(np.radians(deg)), r*np.sin(np.radians(deg)), z])


call(pol(0, 0, H), (0.52, 0.98),
     f"Water plane  Ø685.4 mm  ·  {d:.1f} mm deep\nsurface meets the cup "
     f"bottoms  ·  highly responsive", ha="center")
call(pol(Rm+A, 12, H+HZ), (0.99, 0.60),
     "Gather-to-center horn ×21\nmouth above the water line\n"
     "axes converge at center")
call(pol(rt-8, 74, H+Z+6), (0.02, 0.80),
     f"Concave surface cup ×21\nbottoms sit AT the water surface")
call(pol(235, 45, 12), (0.02, 0.50),
     "Smooth φ flare\nwall = 2φ = 3.24 mm")
call(pol(0, 0, 3), (0.26, 0.02),
     "Base Ø261.8 mm\nflat footprint", ha="center")

# ---- cross-section (vertical exaggeration for the shallow profile) --------
EXAG = 5.0
axc = fig.add_axes([0.60, 0.53, 0.40, 0.40])
Y = lambda v: v * EXAG
wall = sp.smoothstep_wall(rb, rt, t+rf, H, 220)
ow = wall + t*sp.outward_normals(wall)
ang = np.linspace(-np.pi/2, 0, sp.ARC)
fil = np.column_stack([(rb-rf)+rf*np.cos(ang), (t+rf)+rf*np.sin(ang)])
inner = np.vstack([[0, t], [rb-rf, t], fil, wall])
outer = np.vstack([ow[::-1], [rb+t, 0], [0, 0]])
px = np.r_[inner[:, 0], rt+rim, outer[:, 0]]
py = np.r_[inner[:, 1], H, outer[:, 1]]
for s in (1, -1):
    axc.fill(s*px, Y(py), color="#6ea8d8", alpha=0.32, lw=0)
    axc.plot(s*px, Y(py), color=INK, lw=1.6)
axc.fill_between([-rt, rt], [Y(H), Y(H)], [Y(t), Y(t)], color="#bfe3f2", alpha=0.5, lw=0)
axc.plot([-rt, rt], [Y(H), Y(H)], color=ACC, lw=1, ls="--")
axc.annotate("", (rt+20, Y(t)), (rt+20, Y(H)), arrowprops=dict(arrowstyle="<->", color=MID))
axc.text(rt+30, Y((t+H)/2), f"depth {d:.1f} = {depth_lbl}", rotation=90,
         va="center", fontsize=9, color=MID)
axc.annotate("", (-rb, -40), (rb, -40), arrowprops=dict(arrowstyle="<->", color=MID))
axc.text(0, -66, "base Ø261.8 = u·φ²", ha="center", fontsize=9, color=MID)
axc.annotate("", (-rt, Y(H)+40), (rt, Y(H)+40), arrowprops=dict(arrowstyle="<->", color=MID))
axc.text(0, Y(H)+66, "surface Ø685.4 = u·φ⁴", ha="center", fontsize=9, color=MID)
axc.set_title(f"Cross-section — golden ladder  {ladder}   (height ×{EXAG:.0f})",
              fontsize=11, color=INK, weight="bold")
axc.set_aspect("equal"); axc.axis("off"); axc.set_ylim(-95, Y(H)+95)

# ---- spec panel ----------------------------------------------------------
axs = fig.add_axes([0.60, 0.05, 0.40, 0.44]); axs.axis("off")
rows = [
    ("Golden seed  u", "100 mm  (sizes base & surface)"),
    ("Ratio  depth:base:surface", ladder),
    ("φ (golden ratio)", "1.618033989…"),
    ("Water surface Ø", "685.41 mm   (u·φ⁴)"),
    ("Base Ø", "261.80 mm   (u·φ²)"),
    ("Water depth", f"{d:.2f} mm   ({depth_lbl})"),
    ("Wall / floor", "3.236 mm   (2φ)"),
    ("Base fillet / rim lip", "5.24 mm / 8.47 mm"),
    ("Pool Ø / height", f"691.9 mm / {H:.1f} mm"),
    ("Collectors", "21 cups (at surface) + 21 horns (above)"),
    ("Collar outer Ø", f"{info['collar_dia']:.1f} mm"),
    ("Water capacity", f"{cap:.2f} L"),
    ("Holds water", "yes — wall solid below the surface"),
    ("Mesh", "watertight · 1 body · genus 21"),
]
axs.text(0.0, 1.0, "SPECIFICATION", fontsize=12, weight="bold", color=INK,
         transform=axs.transAxes)
y = 0.92
for k, v in rows:
    axs.text(0.0, y, k, fontsize=9.5, color=MID, transform=axs.transAxes)
    axs.text(0.46, y, v, fontsize=9.5, color=INK, weight="bold",
             transform=axs.transAxes)
    y -= 0.070

fig.text(0.30, 0.965, "GOLDEN SCRYING POOL", ha="center", fontsize=20,
         weight="bold", color=INK)
fig.text(0.30, 0.935,
         "a shallow φ-proportioned energy-collection plane — vibration "
         "transmuted to the water surface", ha="center", fontsize=11.5, color=MID)
fig.savefig("scrying_pool_annotated.png", dpi=135, bbox_inches="tight",
            facecolor="white")
print(f"depth={d:.2f} mm  ladder=1:phi^{nb}:phi^{ns}  capacity={cap:.2f} L")
print("wrote scrying_pool_annotated.png")
