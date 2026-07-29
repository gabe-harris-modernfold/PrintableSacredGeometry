"""Anatomy of the collector ring — named / discussed (blueprint style)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
import trimesh, pymeshlab as ml
import scrying_pool as sp
import add_collectors as ac

BG, PANEL, HEAD = "#0c2c4d", "#0a2540", "#123457"
LINE, INK, CY = "#9cc2e6", "#eaf2ff", "#6fd0e8"
TEXT, DIM, GRID, WARN = "#d4e2f4", "#7f9bbd", "#3d5f86", "#e0a24b"
FIGW, FIGH = 18.0, 11.5
plt.rcParams["font.family"] = "sans-serif"

rt = 0.5*sp.SEED*sp.PHI**4
H = sp.WALL + sp.DEPTH
Rm = rt + ac.COLLAR_DR
A = ac.COLLAR_A
Rout, Rin = Rm + A, Rm - A
N = ac.N_CUPS
foc = ac.CARVE_R/2

fig = plt.figure(figsize=(FIGW, FIGH), facecolor=BG)
dec = fig.add_axes([0, 0, 1, 1]); dec.set_xlim(0, FIGW); dec.set_ylim(0, FIGH)
dec.set_aspect("equal"); dec.axis("off")
for gx in np.arange(0.5, FIGW, 0.5): dec.axvline(gx, color=GRID, lw=0.4, alpha=0.12)
for gy in np.arange(0.5, FIGH, 0.5): dec.axhline(gy, color=GRID, lw=0.4, alpha=0.12)
dec.add_patch(Rectangle((0.3, 0.3), FIGW-0.6, FIGH-0.6, fill=False, ec=LINE, lw=1.6))
for cx, cy in [(0.8, FIGH-0.8), (FIGW-0.8, FIGH-0.8), (0.8, 0.8), (FIGW-0.8, 0.8)]:
    dec.add_patch(Circle((cx, cy), 0.14, fill=False, ec=LINE, lw=1.0))
    dec.plot([cx-0.24, cx+0.24], [cy, cy], color=LINE, lw=0.9)
    dec.plot([cx, cx], [cy-0.24, cy+0.24], color=LINE, lw=0.9)
fig.text(0.5, 0.955, "ANATOMY OF THE COLLECTOR RING", ha="center", color=INK, fontsize=22, weight="bold")
fig.text(0.5, 0.928, "THE GOLDEN SCRYING POOL  ·  named & described", ha="center", color=CY, fontsize=11.5, weight="bold")
fig.text(0.5, 0.018, "DWG. SGP-φ-001  ·  RING DETAIL  ·  UNITS: mm", ha="center", color=DIM, fontsize=9, family="monospace")

mesh = trimesh.load("scrying_pool_collectors.stl")
seg = trimesh.intersections.mesh_plane(mesh, [0, 1, 0], [0, 0, 0])


def secpanel(rect, title):
    ax = fig.add_axes(rect); ax.set_facecolor(PANEL); ax.set_aspect("equal")
    for s in ax.spines.values(): s.set_color(LINE)
    ax.tick_params(colors=TEXT, labelsize=7.5)
    ax.set_title(title, color=CY, fontsize=10.5, weight="bold")
    return ax


# ---- 3D detail -----------------------------------------------------------
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse", targetfacenum=200000, preservenormal=True)
ms.save_current_mesh("_anat_dec.stl")
dm = trimesh.load("_anat_dec.stl")
c = dm.triangles_center
ang = np.degrees(np.arctan2(c[:, 1], c[:, 0]))
keep = (np.abs(ang) < 32) & (np.linalg.norm(c[:, :2], axis=1) > 250)
sub = dm.submesh([keep], append=True)
ax3 = fig.add_axes([0.005, 0.36, 0.39, 0.52], projection="3d"); ax3.patch.set_alpha(0)
Ll = np.array([0.4, 0.5, 0.85]); Ll /= np.linalg.norm(Ll)
it = np.clip(0.32 + 0.68*np.clip(sub.face_normals@Ll, 0, 1), 0, 1)
cols = np.column_stack([0.6*it, 0.76*it, 0.9*it, np.ones(len(it))])
ax3.add_collection3d(Poly3DCollection(sub.vertices[sub.faces], facecolor=cols, lw=0))
lo, hi = sub.bounds
ax3.set_xlim(lo[0], hi[0]); ax3.set_ylim(lo[1], hi[1]); ax3.set_zlim(lo[2], lo[2]+(hi[0]-lo[0]))
ax3.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[0]-lo[0]))
ax3.view_init(elev=20, azim=6); ax3.set_axis_off()
ax3.set_title("FIG. A  —  RING DETAIL (isometric)", color=CY, fontsize=10.5, weight="bold", y=0.98)
fig.canvas.draw()


def call3(feat, frac, txt, ha="left"):
    x2, y2, _ = proj3d.proj_transform(*feat, ax3.get_proj())
    ax3.annotate(txt, xy=(x2, y2), xycoords=ax3.transData, xytext=frac,
                 textcoords="axes fraction", fontsize=9, color=INK, ha=ha, va="center",
                 weight="bold", arrowprops=dict(arrowstyle="-", color=CY, lw=1.0),
                 bbox=dict(boxstyle="round,pad=0.25", fc=PANEL, ec=CY, lw=0.7))


def pol(r, deg, z): return (r*np.cos(np.radians(deg)), r*np.sin(np.radians(deg)), z)
call3(pol(325, 3, H+3), (0.02, 0.86), "reflecting cup")
call3(pol(Rout, 8.6, H+ac.HORN_Z), (0.98, 0.72), "horn mouth", ha="right")
call3(pol(Rm, -14, H+ac.COLLAR_Z), (0.02, 0.30), "collar body")
call3(pol(rt, -2, H), (0.98, 0.16), "rim / water line", ha="right")

# ---- top plan schematic --------------------------------------------------
axp = fig.add_axes([0.40, 0.37, 0.27, 0.53]); axp.set_aspect("equal"); axp.axis("off")
axp.add_patch(Circle((0, 0), rt, fc="#123f66", ec=CY, lw=1.0))                 # water
axp.add_patch(Circle((0, 0), Rout, fc="none", ec=LINE, lw=1.4))               # collar outer
axp.add_patch(Circle((0, 0), Rin, fc="none", ec=LINE, lw=1.0, alpha=0.7))     # collar inner
for k in range(N):
    a1 = 2*np.pi*k/N
    axp.add_patch(Circle((325*np.cos(a1), 325*np.sin(a1)), 15, fc=PANEL, ec=INK, lw=1.0))  # cup
    a2 = 2*np.pi*(k+0.5)/N
    u = np.array([np.cos(a2), np.sin(a2)])
    axp.plot([Rin*u[0], (Rout+6)*u[0]], [Rin*u[1], (Rout+6)*u[1]], color=WARN, lw=1.4)       # horn
    axp.plot([0, Rin*u[0]], [0, Rin*u[1]], color=WARN, lw=0.4, alpha=0.35)                   # axis
axp.plot(0, 0, marker="o", ms=5, color=INK)
axp.set_xlim(-Rout*1.35, Rout*1.35); axp.set_ylim(-Rout*1.28, Rout*1.38)
axp.set_title("FIG. B  —  TOP PLAN  (21 cups + 21 horns, half-pitch offset)",
              color=CY, fontsize=10.5, weight="bold")
axp.annotate("reflecting cup ×21", (325, 0), (250, 250), color=INK, fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color=INK))
ah = np.radians(360/N*0.5)
axp.annotate("gather-to-center horn ×21", ((Rout)*np.cos(ah), (Rout)*np.sin(ah)),
             (-120, 330), color=WARN, fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color=WARN))
axp.text(0, 70, "axes converge\nat centre", color=WARN, fontsize=8, ha="center", style="italic")
axp.text(0, -60, "water plane", color=CY, fontsize=9, ha="center")

# ---- cup radial section --------------------------------------------------
axc = secpanel([0.02, 0.055, 0.30, 0.275], "FIG. C  —  SECTION THROUGH A CUP")
axc.fill_between([250, 395], [H, H], [0, 0], color="#12466e", alpha=0.5, lw=0)
axc.axhline(H, color=CY, lw=1.0, ls="--")
for s in seg:
    if s[:, 0].mean() > 0: axc.plot(s[:, 0], s[:, 2], color=INK, lw=1.3)
axc.plot(325, H+foc, marker="+", ms=9, color=WARN, mew=1.6)
axc.set_xlim(295, 392); axc.set_ylim(-3, 66)
axc.set_xlabel("radius (mm)", color=TEXT, fontsize=8); axc.set_ylabel("z", color=TEXT, fontsize=8)
axc.annotate("reflecting cup\n(concave scoop)", (322, H+6), (300, 55), color=INK, fontsize=8.3,
             arrowprops=dict(arrowstyle="->", color=INK))
axc.annotate("cup mouth @ water line", (330, H), (298, 20), color=CY, fontsize=8.3,
             arrowprops=dict(arrowstyle="->", color=CY))
axc.annotate(f"focal pt  f≈{foc:.0f} mm", (325, H+foc), (345, 55), color=WARN, fontsize=8.3,
             arrowprops=dict(arrowstyle="->", color=WARN))
axc.text(370, 8, "collar body", color=DIM, fontsize=8)

# ---- horn radial section -------------------------------------------------
axh = secpanel([0.335, 0.055, 0.30, 0.275], "FIG. D  —  SECTION THROUGH A HORN")
axh.fill_between([250, 395], [H, H], [0, 0], color="#12466e", alpha=0.5, lw=0)
axh.axhline(H, color=CY, lw=1.0, ls="--")
for s in seg:
    if s[:, 0].mean() < 0: axh.plot(-s[:, 0], s[:, 2], color=INK, lw=1.3)
axh.set_xlim(295, 392); axh.set_ylim(-3, 66)
axh.set_xlabel("radius (mm)", color=TEXT, fontsize=8); axh.set_ylabel("z", color=TEXT, fontsize=8)
axh.annotate("", (301, H+ac.HORN_Z), (321, H+ac.HORN_Z),          # beam: throat -> centre
             arrowprops=dict(arrowstyle="-|>", color=WARN, lw=1.8))
axh.text(299, H+ac.HORN_Z+5, "axis → centre", color=WARN, fontsize=8, style="italic")
axh.annotate(f"horn mouth\n(outer Ø{2*ac.HORN_MOUTH:.0f})", (376, H+ac.HORN_Z), (344, 57),
             color=INK, fontsize=8.3, arrowprops=dict(arrowstyle="->", color=INK))
axh.annotate(f"throat (inner Ø{2*ac.HORN_THROAT:.0f})", (322, H+ac.HORN_Z-2), (300, 20),
             color=INK, fontsize=8.3, arrowprops=dict(arrowstyle="->", color=INK))

# ---- nomenclature panel --------------------------------------------------
import textwrap
axn = fig.add_axes([0.685, 0.055, 0.29, 0.83]); axn.axis("off")
axn.add_patch(Rectangle((0, 0), 1, 1, transform=axn.transAxes, fc=PANEL, ec=LINE, lw=1.0))
axn.add_patch(Rectangle((0, 0.955), 1, 0.045, transform=axn.transAxes, fc=HEAD, ec=LINE, lw=1.0))
axn.text(0.03, 0.977, "NOMENCLATURE", transform=axn.transAxes, color=INK, fontsize=12, weight="bold", va="center")
items = [
    ("COLLAR", "raised golden ring fused to the pool rim; carries both collector sets."),
    ("SEAM", "continuous watertight fusion of collar to the pool's flared wall."),
    ("REFLECTING CUP  (×21)", f"concave scoop on the INNER face; its bottom sits at the water line, concentrating energy onto the surface (f≈{foc:.0f} mm)."),
    ("CUP MOUTH", "open rim of a cup, at the still water surface."),
    ("HORN  (×21)", "radial funnel through the collar, staggered half a pitch from the cups."),
    ("HORN MOUTH", f"wide OUTER opening (Ø{2*ac.HORN_MOUTH:.0f} mm); gathers ambient sound."),
    ("HORN THROAT", f"narrow INNER opening (Ø{2*ac.HORN_THROAT:.0f} mm); emits toward the centre."),
    ("HORN AXIS", "all 21 axes converge on the pool centre, above the surface."),
    ("INNER / OUTER FACE", "water-side (cups) vs world-side (horn mouths)."),
    ("RIM / WATER LINE", "perimeter of the still surface where cups meet water."),
    ("PITCH", f"{360/N:.2f}° between like features; cups/horns offset {360/N/2:.2f}°."),
]
y = 0.945
for name, desc in items:
    axn.text(0.035, y, name, transform=axn.transAxes, color=CY, fontsize=9.4, weight="bold", va="top")
    y -= 0.030
    for ln in textwrap.wrap(desc, 52):
        axn.text(0.05, y, ln, transform=axn.transAxes, color=TEXT, fontsize=8.5, va="top")
        y -= 0.025
    y -= 0.010

fig.savefig("scrying_pool_anatomy.png", dpi=125, facecolor=BG, bbox_inches="tight")
print("wrote scrying_pool_anatomy.png")
