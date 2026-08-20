"""True 1:1 A4-landscape PDF of the pentagon cutting template."""
import sys, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dodeca_frame import MIRROR_SIDE, P_CIR, max_pentagon_in_hexagon

W, H = 297.0, 210.0          # A4 landscape, mm
HEX, PENT = 100.0, MIRROR_SIDE
cx, cy = W / 2.0, H / 2.0 - 4.0

fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.axis("off")

hexp = [(cx + HEX * math.cos(math.radians(60 * k)),
         cy + HEX * math.sin(math.radians(60 * k))) for k in range(6)]
rc = PENT * P_CIR
pentp = [(cx + rc * math.cos(math.radians(72 * k)),
          cy + rc * math.sin(math.radians(72 * k))) for k in range(5)]

ax.add_patch(plt.Polygon(hexp, closed=True, fill=False, ec="#999",
                         lw=0.5, ls=(0, (4, 2))))
ax.add_patch(plt.Polygon(pentp, closed=True, fill=False, ec="black", lw=0.8))
ax.plot([cx], [cy], "k+", ms=6, mew=0.6)

t = dict(family="DejaVu Sans", ha="left", va="baseline")
ax.text(8, H - 10, "PENTAGON CUT TEMPLATE", fontsize=9, weight="bold", **t)
ax.text(8, H - 16.5,
        "Print at 100% / Actual Size (NOT 'fit to page'). A4 or Letter, landscape.",
        fontsize=7, **t)
ax.text(8, H - 22,
        "Solid line = cut.  Dashed = your 10 cm hexagon, for registration.",
        fontsize=7, **t)
ax.text(8, H - 27.5,
        "Align one pentagon corner to one hexagon corner, as drawn. Cut all 5 lines. 12 needed.",
        fontsize=7, **t)
ax.text(8, 20, f"Pentagon side {PENT:.1f} mm   |   hexagon side {HEX:.1f} mm   |   "
               f"theoretical max {max_pentagon_in_hexagon(HEX)[0]:.1f} mm",
        fontsize=7, **t)

# 100 mm scale bar for verifying print scale
ax.plot([8, 108], [11, 11], "k-", lw=0.7)
for x in (8, 108):
    ax.plot([x, x], [8, 14], "k-", lw=0.7)
ax.text(112, 9.5, "100 mm - measure this before cutting glass", fontsize=7, **t)

out = sys.argv[1]
fig.savefig(os.path.join(out, "pentagon-cut-template.pdf"),
            format="pdf", facecolor="white")
print("wrote pentagon-cut-template.pdf (A4 landscape, 1:1)")
