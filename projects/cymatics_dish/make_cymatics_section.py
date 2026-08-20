#!/usr/bin/env python3
"""Cross-section of the cymatics dish bridging over the assumed 12" surround."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import cymatics_dish as C

OUT = "cymatics_dish_section.png"
INK, WATER_C, CONE_C, SKIRT_C = "#1b1b1f", "#3d8fd1", "#b04a2a", "#2e7d54"
RING_C = "#9a6b1f"


def fill_both(ax, prof, **kw):
    """Fill the section on both sides of the axis as two separate polygons.

    Mirroring an annular profile into ONE polygon closes it across the axis and
    floods the middle -- which is what a plain tube like the support ring does."""
    p = np.asarray(prof)
    ax.fill(p[:, 0], p[:, 1], **kw)
    kw.pop("label", None)
    ax.fill(-p[:, 0], p[:, 1], **kw)


def speaker(ax, lw=3):
    """Cone, dust cap dome, surround roll and frame lip, both sides."""
    r_mouth, r_cap = C.CONE_MOUTH_D / 2, C.CONE_DUSTCAP_D / 2
    for s in (1, -1):
        ax.plot([s * r_mouth, s * r_cap], [C.cone_z(r_mouth), C.cone_z(r_cap)],
                color=CONE_C, lw=lw, solid_capstyle="round")
        a = np.linspace(0, np.pi / 2, 40)
        ax.plot(s * r_cap * np.cos(a), C.cone_z(r_cap) + 12 * np.sin(a),
                color=CONE_C, lw=2, alpha=.55)
        a = np.linspace(np.pi, 0, 40)                    # half-roll surround
        ax.plot(s * (r_mouth + 9 + 9 * np.cos(a)),
                C.cone_z(r_mouth) + C.SURROUND_H * np.sin(a),
                color=CONE_C, lw=2, alpha=.55)
        ax.plot([s * (r_mouth + 18), s * (r_mouth + 34)],
                [C.cone_z(r_mouth)] * 2, color=CONE_C, lw=4, alpha=.55)


def main():
    r_mouth, r_cap = C.CONE_MOUTH_D / 2, C.CONE_DUSTCAP_D / 2
    r_bot = C.SKIRT_D_BOT / 2
    rc, zc = C.bump_centre()

    fig, ax = plt.subplots(figsize=(13, 5.6))
    speaker(ax)

    fill_both(ax, C.dish_profile(), color=INK, zorder=3,
              label=f"dish — {C.RIM_T:.0f} mm rim, {C.FLOOR_T:.2f} mm floor")
    fill_both(ax, C.skirt_profile(), color=SKIRT_C, zorder=3,
              label=f"skirt — flange + collar + band ({C.SKIRT_T:.0f} mm)")
    fill_both(ax, C.support_profile(), color=RING_C, alpha=.85, zorder=3,
              label=f"Ø{C.SUPPORT_D:.0f} floor support ring")

    ax.fill_between([-C.R_IN, C.R_IN], 0, C.WATER, color=WATER_C, alpha=.45, zorder=4)
    ax.plot([-C.R_IN, C.R_IN], [C.WATER] * 2, color=WATER_C, lw=1.6, zorder=5)

    def note(x, z, txt, dx, dz, col=INK):
        ax.annotate(txt, xy=(x, z), xytext=(x + dx, z + dz), fontsize=9, color=col,
                    ha="left" if dx > 0 else "right", va="center",
                    arrowprops=dict(arrowstyle="-", lw=.9, color=col, alpha=.7))

    note(0, C.WATER / 2, f"{C.WATER:.0f} mm water — "
                         f"{np.pi * C.R_IN**2 * C.WATER / 1000:.0f} ml", 14, 26, WATER_C)
    note(C.R_OUT, C.RIM_H / 2, f"dish: Ø{C.OD:.0f}, {C.RIM_T:.0f} mm rim,\n"
                               f"{C.RIM_H:.2f} mm tall, {C.FLOOR_T:.2f} mm flat floor",
         46, 22)
    note(-C.R_IN, C.RIM_H / 2, f"Ø{2 * C.R_IN:.0f} water surface,\n"
                               f"{C.RIM_H - C.WATER:.2f} mm freeboard", 30, 24)
    note(-(C.R_OUT + C.R_CON) / 2, -C.FLOOR_T - C.SKIRT_T,
         f"Ø{C.OD:.0f} flange — glue face + {C.N_VENT} vents,\n"
         f"carries the {C.R_OUT - C.R_CON:.1f} mm overhang",
         -30, -30, SKIRT_C)
    note(rc, zc, f"{C.N_BUMP}× Ø{C.BUMP_D:.0f} bumps on the collar,\n"
                 f"{rc + C.BUMP_D / 2 - C.R_CON:.1f} mm proud, "
                 f"{C.BUMP_CLR:.1f} mm off the cone", 40, 26)
    note(-C.R_CON, -C.FLOOR_T - C.COLLAR_H / 2,
         f"{C.COLLAR_H:.1f} mm collar —\nbridges the surround", -20, -8, SKIRT_C)
    note(-(C.R_CON + r_bot) / 2, (C.cone_z(C.R_CON) + C.cone_z(r_bot)) / 2,
         f"{C.SKIRT_T:.0f} mm band, Ø{C.CONTACT_D:.0f}→Ø{C.SKIRT_D_BOT:.0f}\n"
         f"(bond to the cone here)", -26, -14, SKIRT_C)
    note(-r_mouth, C.cone_z(r_mouth), '12" cone (assumed) —\nmouth / frame plane',
         -14, 20, CONE_C)
    note(-r_bot, C.cone_z(r_bot), f"Ø{C.SKIRT_D_BOT:.0f} toe, feathered to "
                                  f"{C.SKIRT_T_TOE:.1f} mm —\nclears the dust cap by "
                                  f"{r_bot - r_cap:.0f} mm", -20, -16, SKIRT_C)

    # --- inset: rim, overhang, flange, collar and a bump ---
    ins = ax.inset_axes([0.795, 0.055, 0.185, 0.66])
    speaker(ins, lw=2.5)
    fill_both(ins, C.dish_profile(), color=INK, zorder=3)
    fill_both(ins, C.skirt_profile(), color=SKIRT_C, zorder=3)
    ins.fill_between([C.R_IN - 12, C.R_IN], 0, C.WATER,
                     color=WATER_C, alpha=.45, zorder=4)
    for s in (1, -1):
        ins.add_patch(plt.Circle((s * rc, zc), C.BUMP_D / 2, color=SKIRT_C, zorder=4))
    ins.set_xlim(C.R_CON - 10, C.R_OUT + 12)
    ins.set_ylim(C.cone_z(C.R_CON) - 2, C.RIM_H + 1.5)
    ins.set_aspect("equal")
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_title("overhang, flange, collar, bump", fontsize=8.0, pad=3)
    for sp in ins.spines.values():
        sp.set_alpha(.35)
    ax.indicate_inset_zoom(ins, edgecolor=INK, alpha=.35)

    ax.set_title(f"Cymatics petri dish — Ø{C.OD:.0f} × {C.RIM_H:.0f} mm, bridging "
                 f"over a Ø{C.CONE_MOUTH_D:.0f}/{C.CONE_DEPTH:.0f} mm cone "
                 f"({np.degrees(np.arctan(C.K)):.1f}° slope)", fontsize=12, pad=12)
    ax.set_xlabel("radius (mm)")
    ax.set_ylabel("height above the water side of the floor (mm)")
    ax.set_aspect("equal")
    ax.grid(alpha=.15)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    ax.set_xlim(-262, 262)
    ax.set_ylim(-96, 40)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
