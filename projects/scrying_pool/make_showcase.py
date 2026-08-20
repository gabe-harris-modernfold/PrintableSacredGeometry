"""Engineering / blueprint showcase sheet for the golden scrying pool."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml

import scrying_pool as sp

# ---- blueprint palette ---------------------------------------------------
BG    = "#0c2c4d"      # blueprint field
PANEL = "#0a2540"
LINE  = "#9cc2e6"      # pale line-blue (frames, rules)
INK   = "#eaf2ff"      # near-white
CY    = "#6fd0e8"      # cyan accent
TEXT  = "#d4e2f4"
DIM   = "#7f9bbd"
GRID  = "#3d5f86"

PHI = sp.PHI
depth = sp.DEPTH
rt = 0.5 * sp.SEED * PHI**4
c_wave = np.sqrt(9.81 * depth / 1000.0)
H = sp.WALL + sp.DEPTH

plt.rcParams["font.family"] = "sans-serif"

# ---- hero mesh -----------------------------------------------------------
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse",
                targetfacenum=150000, preservenormal=True)
ms.save_current_mesh("_show_dec.stl")
m = trimesh.load("_show_dec.stl")
L1 = np.array([0.25, 0.5, 0.9]); L1 /= np.linalg.norm(L1)
L2 = np.array([-0.7, -0.2, 0.5]); L2 /= np.linalg.norm(L2)
nrm = m.face_normals
it = np.clip(0.28 + 0.70*np.clip(nrm@L1, 0, 1) + 0.26*np.clip(nrm@L2, 0, 1), 0, 1)
steel = np.array([0.62, 0.78, 0.92])
cols = np.column_stack([np.clip(steel[0]*it, 0, 1), np.clip(steel[1]*it, 0, 1),
                        np.clip(steel[2]*it, 0, 1), np.ones(len(it))])

fig = plt.figure(figsize=(18, 11.5), facecolor=BG)

# ---- backdrop: grid, frame, registration marks ---------------------------
dec = fig.add_axes([0, 0, 1, 1]); dec.set_xlim(0, 18); dec.set_ylim(0, 11.5)
dec.set_aspect("equal"); dec.axis("off")
for gx in np.arange(0.5, 18, 0.5):
    dec.axvline(gx, color=GRID, lw=0.4, alpha=0.18)
for gy in np.arange(0.5, 11.5, 0.5):
    dec.axhline(gy, color=GRID, lw=0.4, alpha=0.18)
dec.add_patch(plt.Rectangle((0.35, 0.35), 17.30, 10.80, fill=False, ec=LINE, lw=1.6))
dec.add_patch(plt.Rectangle((0.5, 0.5), 17.0, 10.5, fill=False, ec=LINE, lw=0.6, alpha=0.6))
for cx, cy in [(0.9, 10.6), (17.1, 10.6), (0.9, 0.9), (17.1, 0.9)]:
    dec.add_patch(Circle((cx, cy), 0.16, fill=False, ec=LINE, lw=1.0))
    dec.plot([cx-0.28, cx+0.28], [cy, cy], color=LINE, lw=0.9)
    dec.plot([cx, cx], [cy-0.28, cy+0.28], color=LINE, lw=0.9)

# ---- title ---------------------------------------------------------------
fig.text(0.5, 0.955, "THE GOLDEN SCRYING POOL", ha="center", color=INK,
         fontsize=29, weight="bold")
fig.text(0.5, 0.919, "φ-PROPORTIONED VIBRATION PLANE  ·  AS ABOVE, SO BELOW",
         ha="center", color=CY, fontsize=12.5, weight="bold")
fig.text(0.5, 0.893,
         "a golden-ratio vessel for the transmutation of vibration into form upon the waters",
         ha="center", color=DIM, fontsize=11, style="italic")

# ---- hero render ---------------------------------------------------------
axh = fig.add_axes([0.005, 0.32, 0.585, 0.60], projection="3d")
axh.patch.set_alpha(0)
pc = Poly3DCollection(m.vertices[m.faces], linewidths=0); pc.set_facecolor(cols)
axh.add_collection3d(pc)
a = np.linspace(0, 2*np.pi, 200)
axh.add_collection3d(Poly3DCollection(
    [np.column_stack([rt*np.cos(a), rt*np.sin(a), np.full_like(a, H)])],
    facecolor=(0.08, 0.20, 0.34, 0.85), edgecolor=CY, lw=1.0))
lo, hi = m.bounds
axh.set_xlim(lo[0], hi[0]); axh.set_ylim(lo[1], hi[1]); axh.set_zlim(0, hi[0]-lo[0])
axh.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[0]-lo[0]))
axh.view_init(elev=24, azim=40); axh.set_axis_off()
fig.text(0.30, 0.345, "FIG. 1  —  ISOMETRIC  ·  COLLECTOR ASSEMBLY  (Ø 761 mm)",
         ha="center", color=DIM, fontsize=9.5)

# ---- cymatic mandala -----------------------------------------------------
axm = fig.add_axes([0.605, 0.40, 0.36, 0.44]); axm.set_aspect("equal"); axm.axis("off")
N = 700
xx, yy = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
r = np.hypot(xx, yy); th = np.arctan2(yy, xx)
patt = (np.cos(9*np.pi*r) * (0.72 + 0.28*np.cos(12*th))
        + 0.45*np.cos(5*np.pi*r) + 0.3*np.cos(21*th)*np.exp(-2*r))
patt = patt * np.clip(1.15 - r, 0, 1)
patt[r > 1] = np.nan
cmap = LinearSegmentedColormap.from_list(
    "blue", [BG, "#0f3f6b", "#2f89b0", CY, INK])
axm.imshow(patt, extent=(-1, 1, -1, 1), cmap=cmap, origin="lower", interpolation="bilinear")
axm.add_patch(Circle((0, 0), 1.0, fill=False, ec=LINE, lw=1.6))
axm.add_patch(Circle((0, 0), 1.06, fill=False, ec=LINE, lw=0.6, alpha=0.5))
for k in range(21):
    ang = 2*np.pi*k/21
    axm.plot(1.11*np.cos(ang), 1.11*np.sin(ang), marker="o", ms=3.4, color=INK)
axm.plot(0, 0, marker="o", ms=6, color=INK)
axm.set_xlim(-1.2, 1.2); axm.set_ylim(-1.2, 1.2)
axm.set_title("FIG. 2  —  FARADAY STANDING WAVE (form born of vibration)",
              color=CY, fontsize=10, weight="bold", pad=6)

# ---- panels --------------------------------------------------------------
def panel(x, y, w, h, title, lines):
    ax = fig.add_axes([x, y, w, h]); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                 fc=PANEL, ec=LINE, lw=1.0))
    ax.add_patch(plt.Rectangle((0, 0.86), 1, 0.14, transform=ax.transAxes,
                 fc="#123457", ec=LINE, lw=1.0))
    ax.text(0.04, 0.925, title, transform=ax.transAxes, color=INK,
            fontsize=13, weight="bold", va="center")
    yy = 0.80
    for s in lines:
        ax.text(0.045, yy, s, transform=ax.transAxes, color=TEXT,
                fontsize=9.4, va="top")
        yy -= 0.066 + 0.056*s.count("\n")

phys = [
    "◆  Cymatics: a vibrated surface forms Faraday standing waves —\n"
    "    nodal geometry appearing at half the drive frequency.",
    f"◆  Shallow by design (23.6 mm): gravity-wave speed\n"
    f"    c = √(gh) ≈ {c_wave:.2f} m/s. Less mass → less damping →\n"
    "    a livelier, faster-settling plane.",
    "◆  Rigid circular rim → clean Bessel modes Jₙ(kr)·cos nθ\n"
    "    (concentric rings and petals).",
    "◆  φ, the most irrational ratio: depth:base:surface =\n"
    "    1 : φ⁵ : φ⁷ spaces resonances quasiperiodically,\n"
    "    suppressing degenerate mode overlap.",
    "◆  21 cups (f = R/2 ≈ 17 mm) focus higher partials; 21 horns\n"
    "    gather ambient sound to a central radiating antinode.",
]
herm = [
    "◇  The still water is a speculum — a threshold mirror\n"
    "    between the seen and the unseen; the scryer reads the\n"
    "    forms that vibration calls forth.",
    "◇  Vibration made visible is the Logos: solve et coagula,\n"
    "    spirit condensing into pattern upon the waters.",
    "◇  φ, the divine proportion — the measure of shell and\n"
    "    flower — tunes the vessel to nature's own harmony.",
    "◇  21 = 3 × 7: the triad through the seven spheres; twin\n"
    "    rings gather the music of the spheres to one point.",
    "◇  Horns draw the outer world inward, cups return it to the\n"
    "    surface — the circuit of correspondence, above below.",
]
panel(0.045, 0.045, 0.44, 0.255, "PHYSICS  ·  THE RESPONSIVE PLANE", phys)
panel(0.515, 0.045, 0.44, 0.255, "HERMETIC  ·  CORRESPONDENCES", herm)

# ---- engineering title-block ribbon --------------------------------------
fig.text(0.30, 0.315,
         "φ = 1.618033989…    surface Ø 685.4 mm    depth 23.6 mm    3.84 L    "
         "1 : φ⁵ : φ⁷    watertight / holds water",
         ha="center", color=CY, fontsize=10.5, family="monospace")
fig.text(0.5, 0.02,
         "DWG. SGP-φ-001   ·   MATERIAL: PLA / RESIN (SEALED)   ·   UNITS: mm   "
         "·   MESH: WATERTIGHT · 1 BODY · GENUS 21   ·   SHEET 1/1",
         ha="center", color=DIM, fontsize=9, family="monospace")

fig.savefig("scrying_pool_showcase.png", dpi=130, facecolor=BG, bbox_inches="tight")
print(f"c_wave={c_wave:.3f} m/s; wrote scrying_pool_showcase.png (blueprint)")
