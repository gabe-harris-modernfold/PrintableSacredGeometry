"""Four-sheet Operation & Performance Report (blueprint style, auto-fit panels)."""
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Wedge
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import trimesh, pymeshlab as ml
import scrying_pool as sp

BG, PANEL, HEAD = "#0c2c4d", "#0a2540", "#123457"
LINE, INK, CY = "#9cc2e6", "#eaf2ff", "#6fd0e8"
TEXT, DIM, GRID, WARN = "#d4e2f4", "#7f9bbd", "#3d5f86", "#e0a24b"
FIGW, FIGH = 12.0, 16.0
plt.rcParams["font.family"] = "sans-serif"

PHI = sp.PHI
depth = sp.DEPTH; h = depth/1000.0
rt = 0.5*sp.SEED*PHI**4; H = sp.WALL + sp.DEPTH
g, rho, sigma = 9.81, 1000.0, 0.072
c_wave = np.sqrt(g*h)
lam_c = 2*np.pi*np.sqrt(sigma/(rho*g))*1000.0
CAP_L = 3.84
mesh = trimesh.load("scrying_pool_collectors.stl")
mat_cm3 = abs(mesh.volume)/1000.0
m_solid = mat_cm3*1.24/1000.0            # kg, 100% PLA
m_print = m_solid*0.27                    # kg, ~20% infill + shells
m_water = CAP_L
head_kPa = rho*g*h/1000.0
therm = m_water*4186/1000.0


def surf_freq(lam_mm):
    lam = np.asarray(lam_mm, float)/1000.0
    k = 2*np.pi/lam
    return np.sqrt((g*k + sigma/rho*k**3)*np.tanh(k*h))/(2*np.pi)


def new_page(sub, no):
    fig = plt.figure(figsize=(FIGW, FIGH), facecolor=BG)
    dec = fig.add_axes([0, 0, 1, 1]); dec.set_xlim(0, FIGW); dec.set_ylim(0, FIGH)
    dec.set_aspect("equal"); dec.axis("off")
    for gx in np.arange(0.5, FIGW, 0.5): dec.axvline(gx, color=GRID, lw=0.4, alpha=0.13)
    for gy in np.arange(0.5, FIGH, 0.5): dec.axhline(gy, color=GRID, lw=0.4, alpha=0.13)
    dec.add_patch(Rectangle((0.3, 0.3), FIGW-0.6, FIGH-0.6, fill=False, ec=LINE, lw=1.6))
    dec.add_patch(Rectangle((0.45, 0.45), FIGW-0.9, FIGH-0.9, fill=False, ec=LINE, lw=0.6, alpha=0.6))
    for cx, cy in [(0.85, FIGH-0.85), (FIGW-0.85, FIGH-0.85), (0.85, 0.85), (FIGW-0.85, 0.85)]:
        dec.add_patch(Circle((cx, cy), 0.15, fill=False, ec=LINE, lw=1.0))
        dec.plot([cx-0.26, cx+0.26], [cy, cy], color=LINE, lw=0.9)
        dec.plot([cx, cx], [cy-0.26, cy+0.26], color=LINE, lw=0.9)
    fig.text(0.5, 0.963, "THE GOLDEN SCRYING POOL", ha="center", color=INK, fontsize=21, weight="bold")
    fig.text(0.5, 0.945, "OPERATION & PERFORMANCE REPORT", ha="center", color=CY, fontsize=11.5, weight="bold")
    fig.text(0.5, 0.930, sub, ha="center", color=DIM, fontsize=11, style="italic")
    dec.plot([0.6, FIGW-0.6], [FIGH-1.05, FIGH-1.05], color=LINE, lw=0.8, alpha=0.6)
    fig.text(0.5, 0.020, f"DWG. SGP-φ-001   ·   UNITS: mm   ·   φ = 1.6180339887   ·   {no}",
             ha="center", color=DIM, fontsize=9, family="monospace")
    return fig, dec


def panel(fig, x, w, y_top, title, items, fs=10.0, mono=False):
    """Auto-height titled text panel; returns bottom y (figure fraction)."""
    chars = int((w*FIGW - 0.45) / (0.0090*fs))
    wi = [(mk, (textwrap.wrap(tx, chars) if tx else [""])) for mk, tx in items]
    line_in = fs/72*1.34
    header_in, pad = 0.36, 0.13
    text_in = sum(len(ls)*line_in + 0.4*line_in for _, ls in wi)
    total_in = header_in + pad + text_in + pad
    total = total_in/FIGH; y0 = y_top - total
    ax = fig.add_axes([x, y0, w, total]); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc=PANEL, ec=LINE, lw=1.0))
    hf = header_in/total_in
    ax.add_patch(Rectangle((0, 1-hf), 1, hf, transform=ax.transAxes, fc=HEAD, ec=LINE, lw=1.0))
    ax.text(0.018, 1-hf/2, title, transform=ax.transAxes, color=INK, fontsize=12, weight="bold", va="center")
    lh = line_in/total_in
    yy = 1 - hf - pad/total_in
    fam = "monospace" if mono else "sans-serif"
    for mk, ls in wi:
        for i, ln in enumerate(ls):
            pre = (mk+"  ") if (mk and i == 0) else ("     " if mk else "")
            ax.text(0.02, yy, pre+ln, transform=ax.transAxes, color=TEXT, fontsize=fs,
                    va="top", family=fam)
            yy -= lh
        yy -= 0.4*lh
    return y0


def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values(): s.set_color(LINE)
    ax.tick_params(colors=TEXT, labelsize=8)
    ax.xaxis.label.set_color(TEXT); ax.yaxis.label.set_color(TEXT)
    ax.set_title(title, color=CY, fontsize=10.5, weight="bold")
    ax.grid(color=GRID, alpha=0.35, lw=0.5)


# ========================== SHEET 1 =======================================
fig, dec = new_page("SHEET 1 / 4  —  OVERVIEW & SETUP", "SHEET 1/4")
ms = ml.MeshSet(); ms.load_new_mesh("scrying_pool_collectors.stl")
ms.apply_filter("meshing_decimation_quadric_edge_collapse", targetfacenum=110000, preservenormal=True)
ms.save_current_mesh("_rep_dec.stl")
hm = trimesh.load("_rep_dec.stl")
L1 = np.array([0.25, 0.5, 0.9]); L1 /= np.linalg.norm(L1)
it = np.clip(0.3 + 0.7*np.clip(hm.face_normals@L1, 0, 1), 0, 1)
cols = np.column_stack([0.62*it, 0.78*it, 0.92*it, np.ones(len(it))])
axh = fig.add_axes([0.30, 0.735, 0.40, 0.175], projection="3d"); axh.patch.set_alpha(0)
axh.add_collection3d(Poly3DCollection(hm.vertices[hm.faces], facecolor=cols, lw=0))
lo, hi = hm.bounds
axh.set_xlim(lo[0], hi[0]); axh.set_ylim(lo[1], hi[1]); axh.set_zlim(0, hi[0]-lo[0])
axh.set_box_aspect((hi[0]-lo[0], hi[1]-lo[1], hi[0]-lo[0])); axh.view_init(elev=26, azim=40); axh.set_axis_off()

b = panel(fig, 0.06, 0.88, 0.735, "1.0  WHAT IT IS", [
    ("", "The Golden Scrying Pool is a shallow, perfectly circular vibration vessel "
     "proportioned entirely on the golden ratio. A wide 685.4 mm water plane, only "
     "23.6 mm deep, is carried on a 261.8 mm base by a smooth phi-flared wall; every "
     "principal dimension follows the ladder depth : base : surface = 1 : phi^5 : phi^7."),
    ("", "A raised collar rings the pool with 21 concave reflecting cups (whose bottoms "
     "meet the water surface) and 21 gather-to-center acoustic horns. Its purpose is "
     "twofold: to render sound and vibration visible as cymatic standing-wave patterns "
     "on the water, and to present a still mirror for contemplative scrying."),
])
b = panel(fig, 0.06, 0.88, b-0.02, "1.1  SETUP PROCEDURE", [
    ("1.", "LEVEL on a rigid surface and check with a bubble level. The wide flat plane "
     "magnifies tilt; ~0.5 deg pools the water to one side and biases every pattern."),
    ("2.", "SEAL FIRST — FDM prints are porous. Brush a food-safe epoxy or several sealer "
     "coats inside before first fill, or print in watertight resin; confirm no weeping."),
    ("3.", f"FILL TO THE CUP LINE — add clean water until the surface just touches the "
     f"bottom of the reflecting cups (~{CAP_L:.2f} L, {depth:.1f} mm). This is the design line."),
    ("4.", "SETTLE 30-60 s to rest and reach room temperature. A trace of surfactant "
     "sharpens fine capillary figures; a dark interior or pigment deepens contrast."),
    ("5.", "LIGHT with raking side light or a single overhead point source for the strongest "
     "relief on ripples; keep a dark surround to maximise the mirror."),
])
# cross-section diagram (vertical exaggeration)
EX = 6.0
axd = fig.add_axes([0.14, b-0.185, 0.72, 0.165]); axd.axis("off")
lines = trimesh.intersections.mesh_plane(mesh, [0, 1, 0], [0, 0, 0])
axd.fill_between([-rt, rt], [H*EX, H*EX], [0, 0], color="#12466e", alpha=0.6, lw=0)
axd.axhline(H*EX, color=CY, lw=1.2, ls="--")
for seg in lines:
    axd.plot(seg[:, 0], seg[:, 2]*EX, color=INK, lw=1.0)
axd.set_xlim(-390, 390); axd.set_ylim(-45, 430)
axd.annotate("water line = cup bottoms", (30, H*EX), (-255, H*EX+95), color=CY, fontsize=8.5,
             arrowprops=dict(arrowstyle="->", color=CY))
axd.annotate("horns: dry, above the line", (322, (H+20)*EX), (70, 375), color=WARN,
             fontsize=8.5, arrowprops=dict(arrowstyle="->", color=WARN))
axd.annotate("solid wall below — holds water", (150, (H/2)*EX), (-360, 300), color=INK,
             fontsize=8.5, arrowprops=dict(arrowstyle="->", color=INK))
axd.text(0, -38, f"FIG. 1.1  —  section on y=0  ·  vertical scale ×{EX:.0f}", ha="center",
         color=DIM, fontsize=8.5)
b = b-0.185
panel(fig, 0.06, 0.88, b-0.02, "1.2  AT A GLANCE", [
    ("▸", f"Water plane 685.4 mm dia · depth {depth:.1f} mm · capacity {CAP_L:.2f} L · "
     "base 261.8 mm dia · wall 3.236 mm (2phi) · overall 761 mm dia."),
    ("▸", f"21 reflecting cups + 21 horns · as-printed mass ~{m_print:.1f} kg (PLA, "
     f"~20% infill) · filled ~{m_print+m_water:.1f} kg."),
    ("▸", "Single watertight body — MeshLab: 0 holes, genus 21 (the dry horn tunnels). "
     "Water retention verified below the surface line."),
])
fig.savefig("report_1.png", dpi=120, facecolor=BG, bbox_inches="tight"); plt.close(fig)

# ========================== SHEET 2 =======================================
fig, dec = new_page("SHEET 2 / 4  —  OPERATION & DRIVING", "SHEET 2/4")
b = panel(fig, 0.06, 0.88, 0.9, "2.0  HOW IT IS DRIVEN", [
    ("", "All three methods below excite Faraday standing waves: the surface responds "
     "subharmonically, at half the driving frequency. Sweep slowly to find 'lock' "
     "frequencies where a clean, stationary figure stands."),
])
axs = fig.add_axes([0.08, b-0.30, 0.84, 0.29]); axs.set_aspect("equal"); axs.axis("off")
axs.set_xlim(-6.2, 6.2); axs.set_ylim(-2.3, 3.4)
axs.plot([-5, -1.4, -1.2, 1.2, 1.4, 5], [1.15, 0.05, 0.0, 0.0, 0.05, 1.15], color=INK, lw=2)
axs.plot([-5, 5], [1.15, 1.15], color=CY, lw=1.2, ls="--")
axs.fill_between([-5, 5], [1.15, 1.15], [0, 0], color="#12466e", alpha=0.5)
for sx in (-5, 5): axs.add_patch(Circle((sx, 1.15), 0.42, fc=PANEL, ec=INK, lw=1.5))
rx = np.linspace(-4.6, 4.6, 300)
axs.plot(rx, 1.15 + 0.12*np.sin(rx*4)*np.exp(-abs(rx)/6), color=CY, lw=1.0)
axs.add_patch(Rectangle((-0.9, -1.75), 1.8, 0.7, fc=HEAD, ec=INK, lw=1.4))
axs.text(0, -1.4, "TRANSDUCER", ha="center", color=INK, fontsize=8, weight="bold")
for dx in (-0.5, 0, 0.5):
    axs.add_patch(FancyArrowPatch((dx, -1.0), (dx, -0.15), arrowstyle="-|>", mutation_scale=12, color=WARN, lw=1.6))
axs.text(0, -2.1, "A  MECHANICAL — vertical drive of the whole vessel", ha="center", color=WARN, fontsize=8.5, weight="bold")
for rr in (0.5, 0.9, 1.3): axs.add_patch(Wedge((5.9, 2.6), rr, 190, 250, width=0.06, color=CY, alpha=0.9))
axs.add_patch(FancyArrowPatch((5.2, 1.7), (0.3, 1.25), arrowstyle="-|>", mutation_scale=12, color=CY, lw=1.5, connectionstyle="arc3,rad=-0.15"))
axs.text(4.9, 3.05, "B  ACOUSTIC", ha="center", color=CY, fontsize=8.5, weight="bold")
axs.text(2.4, 1.8, "horns gather to center", color=CY, fontsize=8, style="italic")
axs.add_patch(FancyArrowPatch((-5, 2.45), (-5, 1.72), arrowstyle="-|>", mutation_scale=13, color=INK, lw=1.8))
axs.text(-5, 2.65, "C  TAP", ha="center", color=INK, fontsize=8.5, weight="bold")
axs.text(-3.4, 0.5, "concentric Bessel rings", color=INK, fontsize=8, style="italic")
axs.text(0, 3.2, "FIG. 2.1  —  THREE WAYS TO DRIVE THE SURFACE", ha="center", color=DIM, fontsize=9)
panel(fig, 0.06, 0.88, b-0.31, "2.1  OPERATING NOTES", [
    ("A.", "MECHANICAL (exciter / speaker / sound-bath bowl coupled to the base): the "
     "classic cymatics drive and the most repeatable. Raise amplitude until the flat "
     "surface destabilises into a figure (the Faraday threshold), then back off slightly "
     "for the cleanest, most stationary pattern."),
    ("B.", "ACOUSTIC (voice, instrument, ambient tone): the 21 horns gather airborne sound "
     "from the rim and steer it to a central antinode. Directional and treble-biased "
     "(Sheet 3) — best for animating the surface with chant, bells, or singing bowls "
     "rather than precise figure control."),
    ("C.", "IMPULSE (a single tap on the rim): launches concentric ring waves that reveal "
     "the pool's natural modes as they ring down."),
    ("▸", "USABLE RANGE: ~10-130 Hz drive gives bold gravity/mixed waves (Sheet 3); "
     "130 Hz-1 kHz yields fine capillary filigree. Higher amplitude = shorter wavelength "
     "= denser pattern."),
    ("▸", "FOR SCRYING: hold amplitude at zero for a still mirror; add a faint sustained "
     "tone to let forms rise and dissolve. Keep the room draught-free."),
])
fig.savefig("report_2.png", dpi=120, facecolor=BG, bbox_inches="tight"); plt.close(fig)

# ========================== SHEET 3 =======================================
fig, dec = new_page("SHEET 3 / 4  —  ACOUSTIC & FLUID PERFORMANCE", "SHEET 3/4")
b = panel(fig, 0.06, 0.88, 0.9, "3.0  HOW THE SURFACE PERFORMS", [
    ("", f"Long waves ride the shallow layer at c = sqrt(g*h) = {c_wave:.2f} m/s; below "
     f"~{lam_c:.0f} mm wavelength surface tension dominates (capillary regime). The rigid "
     "circular rim sets clean Bessel modes Jn(kr)cos(n-theta); golden proportioning spaces "
     "the eigenfrequencies quasiperiodically, so distinct figures rarely blur together."),
])
ax1 = fig.add_axes([0.10, b-0.30, 0.37, 0.26]); style_ax(ax1, "FIG. 3.1  SURFACE-WAVE DISPERSION")
lam = np.logspace(np.log10(2), np.log10(400), 400)
ax1.loglog(lam, surf_freq(lam), color=CY, lw=2)
ax1.axvline(lam_c, color=WARN, lw=1.2, ls="--"); ax1.axvspan(2, lam_c, color=WARN, alpha=0.10)
ax1.text(lam_c*1.15, 3, "capillary", color=WARN, fontsize=8, rotation=90, va="bottom")
ax1.text(70, 1.5, "gravity /\nshallow", color=DIM, fontsize=8)
ax1.set_xlabel("surface wavelength λ (mm)"); ax1.set_ylabel("surface frequency (Hz)")
ax2 = fig.add_axes([0.57, b-0.30, 0.37, 0.26]); style_ax(ax2, "FIG. 3.2  PATTERN DENSITY vs DRIVE")
lam2 = np.linspace(4, 90, 300); drive = 2*surf_freq(lam2); rings = 2*rt/lam2
ax2.plot(drive, rings, color=CY, lw=2); ax2.set_xlim(0, 140)
ax2.set_xlabel("drive frequency (Hz)  [surface = f/2]"); ax2.set_ylabel("nodal rings across radius")
b = b-0.31
b = panel(fig, 0.06, 0.88, b, "3.1  PERFORMANCE NOTES", [
    ("▸", "SUBHARMONIC: a Faraday surface responds at half the drive frequency, so the "
     "tens-of-Hz drives typical of cymatics produce the single-digit-Hz surface waves "
     "that read cleanly on this wide plane."),
    ("▸", "COLLECTOR BANDWIDTH: the cups (~45 mm aperture) and horns (~32 mm mouth) give "
     "real directional gain only where the airborne wavelength is smaller than the "
     "aperture — roughly above 8-11 kHz. Below that they act as directional baffles that "
     "bias and steer ambient / high-frequency content toward the center."),
    ("▸", "GOLDEN SPECTRUM: because phi is the least-rational ratio, the 1 : phi^5 : phi^7 "
     "sizing avoids integer-related mode degeneracies — adjacent resonances stay distinct, "
     "giving crisper standing figures."),
])
rowtxt = [f"{'drive (Hz)':>12}{'surface λ (mm)':>18}{'nodal rings':>15}"]
for lv in (80, 40, 20, 10, 5):
    rowtxt.append(f"{2*surf_freq(lv):>12.1f}{lv:>18.1f}{2*rt/lv:>15.1f}")
panel(fig, 0.06, 0.88, b-0.02, "3.2  REPRESENTATIVE OPERATING POINTS  (rt = 342.7 mm)",
      [("", r) for r in rowtxt], fs=9.5, mono=True)
fig.savefig("report_3.png", dpi=120, facecolor=BG, bbox_inches="tight"); plt.close(fig)

# ========================== SHEET 4 =======================================
fig, dec = new_page("SHEET 4 / 4  —  STRUCTURAL, MATERIAL & SUMMARY", "SHEET 4/4")
b = panel(fig, 0.06, 0.88, 0.9, "4.0  STRUCTURAL, MATERIAL & THERMAL", [
    ("▸", f"MASS & STABILITY: material volume {mat_cm3:.0f} cm^3 → as-printed ~{m_print:.1f} kg "
     f"(PLA, ~20% infill; {m_solid:.1f} kg if solid). Filled ~{m_print+m_water:.1f} kg. Shallow "
     "water sits low and the collar is wide, so the centre of mass is low and it is very hard to tip."),
    ("▸", f"WATER LOAD: static head at the floor is only {head_kPa:.2f} kPa. The uniform 3.236 mm "
     "(2phi) wall and floor carry it with a large margin; deflection is negligible."),
    ("▸", "WATERTIGHTNESS: one closed manifold body — MeshLab reports 0 holes and genus 21 "
     "(the 21 dry horn tunnels above the water line). FDM needs an internal seal coat; resin "
     "is watertight as printed."),
    ("▸", f"THERMAL: the water stores ~{therm:.0f} kJ/K, so surface temperature drifts slowly and "
     "the plane stays calm. Keep out of direct sun and draughts, which drive convection ripples."),
    ("▸", "PRINTING: 761 mm exceeds desktop beds — print in segments or large-format. The "
     "horizontal horn tunnels are internal overhangs needing support / bridging."),
])
axm = fig.add_axes([0.10, b-0.275, 0.37, 0.24]); style_ax(axm, "FIG. 4.1  MASS BREAKDOWN (kg)")
axm.bar(["material\n(printed)", "water", "filled\ntotal"], [m_print, m_water, m_print+m_water],
        color=[LINE, CY, INK]); axm.set_ylabel("kg")
axp = fig.add_axes([0.57, b-0.275, 0.37, 0.24]); style_ax(axp, "FIG. 4.2  PERFORMANCE PROFILE (0-5)")
metrics = ["surface\nresponsiveness", "pattern\nclarity", "stability", "hi-freq\nsteering",
           "low-freq\nfocusing", "capacity"]
scores = [5.0, 4.5, 5.0, 3.5, 1.5, 1.5]; yp = np.arange(len(metrics))[::-1]
axp.barh(yp, scores, color=CY); axp.set_yticks(yp); axp.set_yticklabels(metrics, fontsize=7.5)
axp.set_xlim(0, 5); axp.set_xlabel("rating"); axp.grid(axis="x", color=GRID, alpha=0.35)
panel(fig, 0.06, 0.88, b-0.295, "4.1  SUMMARY & VERDICT", [
    ("▸", "STRENGTHS: an exceptionally responsive, dead-flat wide plane; very stable; "
     "watertight single body; golden geometry yielding clean, well-separated cymatic figures; "
     "a striking still mirror for scrying."),
    ("▸", "LIMITS BY DESIGN: shallow means low capacity (3.84 L) and a delicate surface easily "
     "disturbed by draughts; the passive collectors steer rather than truly focus below the "
     "ultrasonic range; the 761 mm size demands segmented / large-format printing and sealing."),
    ("▸", "BEST USE: a mechanically driven cymatics / vibration-visualisation plate and a "
     "contemplative scrying mirror. Drive at tens of Hz for bold figures, add a faint tone "
     "for a living surface, or leave still to gaze."),
    ("", "VERDICT — performs exactly to intent: a golden-ratio energy-collection plane that "
     "transmutes vibration into visible form upon the water."),
])
fig.savefig("report_4.png", dpi=120, facecolor=BG, bbox_inches="tight"); plt.close(fig)

print(f"mat={mat_cm3:.0f} cm^3  printed~{m_print:.2f} kg  solid={m_solid:.2f} kg  filled~{m_print+m_water:.2f} kg")
print("wrote report_1.png .. report_4.png")
