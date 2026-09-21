"""Render the study's figures from a completed run.

    python figures.py --run run_ra1e7

Writes PNGs beside the run directory. Everything is headless.
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import body_model as bm
import vortex_analysis as va

INK = "#1C2028"
GOLD = "#8F6F24"
EMBER = "#C25B1E"
SMOKE = "#3A4150"
MUTED = "#656B79"

PLUME_CMAP = LinearSegmentedColormap.from_list(
    "plume", ["#0B0E14", "#1D2740", "#3E4A6B", "#8A6A46", "#D4913F", "#F6D9A0"])

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.edgecolor": "#C9C8C0",
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.titlesize": 10, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.grid": False,
})


def load(run):
    d = np.load(os.path.join(run, "fields.npz"))
    with open(os.path.join(run, "meta.json")) as fh:
        meta = json.load(fh)
    return d, meta


def _extent(meta):
    dx = meta["dx"]
    return [-meta["body_x"] * dx, (meta["nx"] - meta["body_x"]) * dx,
            -meta["body_y"] * dx, (meta["ny"] - meta["body_y"]) * dx]


def _sponge_cells(meta):
    """Rows of damped fluid below the open top, which must not be plotted as flow.

    Older runs predate this being recorded, so fall back to the solver's default
    sponge depth of 0.60 m.
    """
    n = meta.get("sponge_cells")
    return int(n) if n else int(round(0.60 / meta["dx"]))


# ---------------------------------------------------------------------------

def fig_body_model(path):
    """The measured boundary condition: silhouette coloured by skin temperature."""
    dx = 0.005
    nx, ny = 260, 380
    solid, seg_id, _ = bm.rasterise(nx, ny, dx, nx // 2, 8)

    temp = np.full((ny, nx), np.nan)
    for idx, s in enumerate(bm.SEGMENTS):
        temp[seg_id == idx] = s["T_skin"]

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(9.6, 6.4), gridspec_kw={"width_ratios": [1, 1.35]})

    ext = [-(nx // 2) * dx, (nx - nx // 2) * dx, -8 * dx, (ny - 8) * dx]
    im = ax.imshow(temp, origin="lower", extent=ext, cmap="inferno",
                   vmin=29.5, vmax=35.2, interpolation="nearest")
    ax.set_xlim(-0.52, 0.52)
    ax.set_ylim(-0.05, 1.80)
    ax.set_aspect("equal")
    ax.set_xlabel("m")
    ax.set_ylabel("height above floor (m)")
    ax.set_title("Measured skin temperature as the wall condition", loc="left")
    ax.axhline(0.8, color=GOLD, lw=1.1, ls="--")
    ax.text(-0.50, 0.845, "0.80 m — laminar → turbulent\n(Ra$_x$ ≈ 3.8×10⁸)",
            color=GOLD, fontsize=7.5, va="bottom", ha="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("°C")

    rows = bm.segment_table()
    order = np.argsort([-r["T_skin"] for r in rows])
    names = [rows[i]["name"].replace("_", " ") for i in order]
    y = np.arange(len(order))
    axb.barh(y, [rows[i]["h_c"] for i in order], height=0.38,
             color=SMOKE, label="convective $h_c$")
    axb.barh(y + 0.40, [rows[i]["h_r"] for i in order], height=0.38,
             color=EMBER, label="radiative $h_r$")
    axb.set_yticks(y + 0.20)
    axb.set_yticklabels(names)
    axb.invert_yaxis()
    axb.set_xlabel("W m$^{-2}$ K$^{-1}$  (16-segment manikin, standing, still air)")
    axb.set_title("Radiation is the larger channel over the torso — and it heats\n"
                  "the walls, not the air", loc="left")
    axb.legend(frameon=False, loc="lower right")
    for i, k in enumerate(order):
        axb.text(rows[k]["h_r"] + 0.12, i + 0.40,
                 f"{100 * rows[k]['radiative_share']:.0f}% rad",
                 fontsize=6.8, va="center", color=MUTED)
    axb.set_xlim(0, 7.4)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_field(d, meta, path, frame=-1):
    """Temperature with streamlines, plus the silhouette."""
    theta, u, v = d["theta"][frame], d["u"][frame], d["v"][frame]
    solid = d["solid"]
    ext = _extent(meta)
    dx = meta["dx"]

    fig, ax = plt.subplots(figsize=(6.0, 8.2))
    sh = np.ma.masked_where(solid, theta)
    ax.imshow(sh, origin="lower", extent=ext, cmap=PLUME_CMAP,
              vmin=0.0, vmax=0.55, interpolation="bilinear")
    ax.imshow(np.ma.masked_where(~solid, solid.astype(float)), origin="lower",
              extent=ext, cmap="gray_r", vmin=0, vmax=1.6, interpolation="nearest")

    ny, nx = theta.shape
    ys = (np.arange(ny) - meta["body_y"]) * dx
    xs = (np.arange(nx) - meta["body_x"]) * dx
    um = np.where(solid, np.nan, u)
    vm = np.where(solid, np.nan, v)
    speed = np.hypot(np.nan_to_num(um), np.nan_to_num(vm))
    ax.streamplot(xs, ys, np.nan_to_num(um), np.nan_to_num(vm),
                  color="#FFFFFF", linewidth=0.5, density=1.5,
                  arrowsize=0.5, broken_streamlines=True)

    ax.axhline(ys[meta["probe_row"]], color=GOLD, lw=0.9, ls=":")
    ax.text(xs[2], ys[meta["probe_row"]] + 0.045, "probe, 0.50 m above the head",
            color=GOLD, fontsize=7.5, va="bottom",
            bbox=dict(fc="#0B0E14", ec="none", alpha=0.72, pad=1.6))
    ax.set_xlim(-0.85, 0.85)
    ax.set_ylim(-0.05, meta["ny"] * dx - _sponge_cells(meta) * dx - 0.05)
    ax.set_aspect("equal")
    ax.set_xlabel("m")
    ax.set_ylabel("height above floor (m)")
    ax.set_title(f"Temperature and streamlines, t = {d['snap_t'][frame]:.1f} s\n"
                 f"Ra$_{{sim}}$ = {meta['ra']:.0e}", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_criteria(d, meta, path, frame=-1):
    """Vorticity vs Q vs lambda_ci over the plume — why lambda_ci is used."""
    u, v, solid = d["u"][frame], d["v"][frame], d["solid"]
    dx = meta["dx"]
    fluid = ~solid
    j0 = meta["body_y"] + int(round(1.45 / dx))
    j1 = meta["ny"] - _sponge_cells(meta)
    i0 = meta["body_x"] - int(round(0.55 / dx))
    i1 = meta["body_x"] + int(round(0.55 / dx))

    sub = lambda A: A[j0:j1, i0:i1]
    fields = [("vorticity  ∂v/∂x − ∂u/∂y", va.vorticity(u, v, dx), "RdBu_r"),
              ("Q-criterion", va.q_criterion(u, v, dx), "PuOr_r"),
              ("swirling strength  λ$_{ci}$", va.swirling_strength(u, v, dx), "magma")]

    ext = [(i0 - meta["body_x"]) * dx, (i1 - meta["body_x"]) * dx,
           (j0 - meta["body_y"]) * dx, (j1 - meta["body_y"]) * dx]

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 5.6), sharey=True)
    for ax, (title, fld, cmap) in zip(axes, fields):
        s = np.ma.masked_where(~sub(fluid), sub(fld))
        lim = float(np.nanpercentile(np.abs(s.compressed()), 99.2)) or 1.0
        kw = dict(vmin=-lim, vmax=lim) if cmap != "magma" else dict(vmin=0, vmax=lim)
        im = ax.imshow(s, origin="lower", extent=ext, cmap=cmap,
                       interpolation="bilinear", **kw)
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_xlabel("m")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

    # Search the same plume region the fitting uses, or the domain's big corner
    # recirculations win the scan and nothing gets marked inside the crop.
    jj, ii = np.mgrid[0:meta["ny"], 0:meta["nx"]]
    region = (fluid & (jj > j0) & (jj < j1) & (ii > i0) & (ii < i1))
    cores = va.find_cores(u, v, dx=1.0, mask=region, threshold=0.20, min_cells=8)
    for c in cores[:6]:
        axes[2].plot((c["col"] - meta["body_x"]) * dx,
                     (c["row"] - meta["body_y"]) * dx,
                     "o", mfc="none", mec="#7CFFB2", mew=1.3, ms=11)
    axes[0].set_ylabel("height above floor (m)")
    fig.suptitle("Vorticity floods the whole shear layer; λ$_{ci}$ narrows it to where "
                 "streamlines actually spiral (circles: detected cores).",
                 x=0.012, ha="left", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _rankine_stats(pool, tol=1e-3):
    """Population comparison of the two idealisations over every fitted core."""
    o = np.array([f["oseen"]["r2"] for f in pool])
    r = np.array([f["rankine"]["r2"] for f in pool])
    return dict(
        n=len(pool),
        oseen_mean=float(o.mean()), rankine_mean=float(r.mean()),
        oseen_median=float(np.median(o)), rankine_median=float(np.median(r)),
        oseen_win=float(np.mean(o > r + tol)),
        rankine_win=float(np.mean(r > o + tol)),
        tie=float(np.mean(np.abs(o - r) <= tol)),
        best_overall=float(np.maximum(o, r).max()),
        frac_above_r2_09=float(np.mean(np.maximum(o, r) > 0.9)),
    )


def fig_rankine(fits, path):
    """Azimuthal profiles of detected cores against Rankine and Oseen.

    Cores are sampled across the quality distribution, not taken from the top of
    it. Showing the four best-fitting cores would be picking the four that most
    flatter whichever profile was used to rank them.
    """
    pool = [f for f in fits if f.get("rankine") and f.get("oseen")]
    if not pool:
        return False
    stats_all = _rankine_stats(pool)
    quality = np.array([max(f["oseen"]["r2"], f["rankine"]["r2"]) for f in pool])
    order = np.argsort(quality)
    picks = [order[int(q * (len(order) - 1))] for q in (0.90, 0.70, 0.45, 0.20)]
    fits = [pool[i] for i in dict.fromkeys(picks)]
    n = len(fits)
    fig, axes = plt.subplots(1, n, figsize=(3.3 * n, 3.5), squeeze=False)
    for ax, f in zip(axes[0], fits):
        r = np.array(f["profile"]["r"])
        vt = np.array(f["profile"]["v_theta"])
        sgn = 1.0 if f["oseen"]["gamma"] >= 0 else -1.0
        rr = np.linspace(r.min(), r.max(), 220)
        ax.plot(r, sgn * vt, "o", ms=4, color=INK, label="simulated", zorder=3)
        ax.plot(rr, sgn * va.rankine_profile(rr, f["rankine"]["gamma"], f["rankine"]["a"]),
                "-", color=GOLD, lw=1.6,
                label=f"Rankine  $R^2$={f['rankine']['r2']:.3f}")
        ax.plot(rr, sgn * va.oseen_profile(rr, f["oseen"]["gamma"], f["oseen"]["a"]),
                "--", color=EMBER, lw=1.6,
                label=f"Oseen  $R^2$={f['oseen']['r2']:.3f}")
        ax.axvline(f["rankine"]["a"], color=GOLD, lw=0.8, ls=":")
        ax.set_xlabel("r  (cells from core centre)")
        ax.legend(frameon=False, fontsize=7.2, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0][0].set_ylabel(r"$v_\theta$  (lattice units)")
    fig.subplots_adjust(top=0.80)
    # State the verdict over the whole population, not the four drawn. On a
    # synthetic Oseen vortex the gap is decisive (1.000 vs 0.866); on real shed
    # cores neither wins, because they are not steady axisymmetric vortices.
    s = stats_all
    fig.suptitle(
        f"Four cores spanning the fit-quality range, drawn from {s['n']} fitted.\n"
        f"Oseen has the better mean $R^2$ ({s['oseen_mean']:.2f} vs {s['rankine_mean']:.2f}) "
        f"because it fails less badly, but Rankine is the better fit on "
        f"{100 * s['rankine_win']:.0f}% of individual cores. Only "
        f"{100 * s['frac_above_r2_09']:.0f}% of cores reach $R^2$ = 0.9 under either.",
        x=0.012, ha="left", fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


def fig_shedding(d, meta, path, tmin=None):
    """Probe time series and its spectrum, with the measured 5 s cadence marked."""
    t = d["hist_t"]
    vp = d["hist_vprobe"] / meta["u_lat"]
    if tmin is None:
        tmin = 0.40 * float(t[-1])          # drop spin-up and the onset step
    res = va.shedding_frequency(t, vp, tmin=tmin)

    fig, (ax, axp) = plt.subplots(1, 2, figsize=(10.2, 3.6),
                                  gridspec_kw={"width_ratios": [1.5, 1]})
    ax.plot(t, vp, lw=0.9, color=SMOKE)
    ax.axvspan(t[0], tmin, color=MUTED, alpha=0.13, lw=0)
    ax.text(0.5 * (t[0] + tmin), ax.get_ylim()[1] * 0.06, "spin-up,\nexcluded",
            ha="center", fontsize=7.5, color=MUTED)
    ax.set_xlabel("t  (s)")
    ax.set_ylabel(r"$v / U_{ff}$ at the probe")
    ax.set_title("Vertical velocity 0.50 m above the head", loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    if res:
        axp.semilogx(res["freqs"][1:], res["psd"][1:], color=SMOKE, lw=1.1)
        axp.axvline(res["freq_hz"], color=EMBER, lw=1.4,
                    label=f"sim, shedding line  {res['period_s']:.1f} s")
        if np.isfinite(res.get("freq_peakcount_hz", np.nan)):
            axp.axvline(res["freq_peakcount_hz"], color=EMBER, lw=1.0, ls=":",
                        label=f"sim, peak count  {res['period_peakcount_s']:.1f} s")
        if np.isfinite(res.get("freq_slow_hz", np.nan)) and \
                abs(res["freq_slow_hz"] - res["freq_hz"]) > 1e-9:
            axp.axvline(res["freq_slow_hz"], color=MUTED, lw=1.0, ls="-.",
                        label=f"slow meander  {res['period_slow_s']:.0f} s")
        axp.axvline(1 / 5.0, color=GOLD, lw=1.3, ls="--",
                    label="measured  5 s (TR-PIV)")
        axp.set_xlabel("frequency (Hz)")
        axp.set_ylabel("PSD")
        axp.legend(frameon=False, fontsize=8)
        axp.set_title("Shedding cadence", loc="left")
        axp.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return res


def torso_layer_profile(d, meta, frames=None, threshold=0.10):
    """Layer half-thickness beside the torso, time-averaged.

    The solver's own profile marched outward from each row's rightmost solid cell,
    which between shoulder and wrist is the arm, not the trunk — so it reported the
    arm's layer and spiked wherever the arm began and ended. Measure from the run of
    solid cells containing the midline instead, and stop at the arm rather than
    crossing it.
    """
    dx = meta["dx"]
    solid = d["solid"]
    ny, nx = solid.shape
    bx, by = meta["body_x"], meta["body_y"]
    if frames is None:
        frames = range(int(0.6 * len(d["snap_t"])), len(d["snap_t"]))

    acc, cnt = np.zeros(ny), np.zeros(ny)
    for k in frames:
        v = d["v"][k]
        for j in range(by, min(ny - 1, by + int(1.70 / dx))):
            row = solid[j]
            if not row[bx]:
                continue
            edge = bx
            while edge + 1 < nx and row[edge + 1]:
                edge += 1                      # right edge of the central run
            i = edge + 1
            prof, cols = [], []
            while i < nx and not row[i]:        # stop at the arm
                prof.append(v[j, i]); cols.append(i); i += 1
            if len(prof) < 4:
                continue
            prof = np.asarray(prof)
            if prof.max() <= 0:
                continue
            below = np.nonzero(prof < threshold * prof.max())[0]
            width = (below[0] if below.size else len(prof)) * dx
            acc[j] += width; cnt[j] += 1

    ok = cnt > 0
    heights = (np.nonzero(ok)[0] - by) * dx
    widths = acc[ok] / cnt[ok]
    return heights, widths


def fig_boundary_layer(d, meta, path):
    """Simulated layer growth, with the Ra-scaling correction made explicit."""
    h, w = torso_layer_profile(d, meta)
    if h.size == 0:
        h, w = d["bl_height"], d["bl_width"]
    if h.size == 0:
        return
    ratio = (meta["ra_real"] / meta["ra"]) ** 0.25

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(100 * w, h, "-", color=SMOKE, lw=1.6,
            label=f"simulated, Ra = {meta['ra']:.0e}")
    ax.plot(100 * w / ratio, h, "--", color=GOLD, lw=1.6,
            label=f"scaled by (Ra$_{{real}}$/Ra$_{{sim}}$)$^{{-1/4}}$ = 1/{ratio:.1f}")
    ax.plot([0.4], [0.15], "o", color=EMBER, ms=7)
    ax.annotate("measured 2–4 mm\nat the ankle", (0.4, 0.15), (3.5, 0.20),
                color=EMBER, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=EMBER, lw=0.9))
    ax.plot([20], [1.25], "o", color=EMBER, ms=7)
    ax.annotate("measured ≈ 20 cm\nat chest height", (20, 1.25), (23, 0.95),
                color=EMBER, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=EMBER, lw=0.9))
    ax.axhline(0.80, color=GOLD, lw=0.9, ls=":")
    ax.text(0.5, 0.82, "0.80 m — measured transition", color=GOLD, fontsize=7.5)
    ax.set_xlabel("layer half-thickness beside the trunk (cm)")
    ax.set_ylabel("height above floor (m)")
    ax.set_xlim(0, 32)
    ax.set_ylim(0, 1.75)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Boundary-layer growth along the body\n"
                 "Shape and trend only — at this Ra the absolute thickness is not a "
                 "prediction", loc="left", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_segment_conductance(meta, path):
    """Simulated per-segment wall conductance against the measured h_c."""
    cond = meta.get("segment_conductance") or {}
    if not cond:
        return
    rows = [r for r in bm.segment_table() if r["name"] in cond]
    if len(rows) < 3:
        return
    sim = np.array([cond[r["name"]] for r in rows])
    meas = np.array([r["h_c"] for r in rows])
    # Compare shape, not magnitude: a 2D slice cannot reproduce 3D coefficients.
    sim_n = sim / sim.mean()
    meas_n = meas / meas.mean()
    order = np.argsort(-meas_n)
    names = [rows[i]["name"].replace("_", " ") for i in order]
    y = np.arange(len(order))

    fig, (ax, axs) = plt.subplots(1, 2, figsize=(10.2, 4.4),
                                  gridspec_kw={"width_ratios": [1.3, 1]})
    ax.barh(y, meas_n[order], height=0.38, color=EMBER, label="measured $h_c$ (normalised)")
    ax.barh(y + 0.40, sim_n[order], height=0.38, color=SMOKE,
            label="simulated wall conductance (normalised)")
    ax.set_yticks(y + 0.20)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.axvline(1.0, color=MUTED, lw=0.7, ls=":")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_xlabel("relative to the segment mean")
    ax.set_title("Per-segment convective conductance", loc="left")

    axs.plot(meas_n, sim_n, "o", color=INK, ms=6)
    for i, r_ in enumerate(rows):
        axs.annotate(r_["name"].replace("_", " "), (meas_n[i], sim_n[i]),
                     (4, 3), textcoords="offset points", fontsize=7, color=MUTED)
    lim = [min(meas_n.min(), sim_n.min()) - 0.1, max(meas_n.max(), sim_n.max()) + 0.1]
    axs.plot(lim, lim, "--", color=MUTED, lw=0.8)

    # One extreme point can carry a Pearson coefficient on its own, so report the
    # rank correlation and the value with the foot removed alongside it.
    stats_ = dict(
        pearson=float(stats.pearsonr(meas_n, sim_n)[0]),
        spearman=float(stats.spearmanr(meas_n, sim_n)[0]),
    )
    keep = [i for i, r_ in enumerate(rows) if r_["name"] != "foot"]
    if len(keep) > 2:
        stats_["pearson_no_foot"] = float(stats.pearsonr(meas_n[keep], sim_n[keep])[0])
        stats_["spearman_no_foot"] = float(stats.spearmanr(meas_n[keep], sim_n[keep])[0])
    fc = None
    cond_map = {r_["name"]: c for r_, c in zip(rows, sim)}
    if "foot" in cond_map and cond_map.get("chest"):
        fc = cond_map["foot"] / cond_map["chest"]
        stats_["foot_over_chest_sim"] = float(fc)
        stats_["foot_over_chest_measured"] = 5.1 / 3.0

    axs.set_xlabel("measured (normalised)")
    axs.set_ylabel("simulated (normalised)")
    axs.set_title(f"Pearson r = {stats_['pearson']:.2f}  ·  Spearman "
                  f"{stats_['spearman']:.2f}\nwithout the foot, Pearson falls to "
                  f"{stats_.get('pearson_no_foot', float('nan')):.2f}",
                  loc="left", fontsize=9)
    axs.spines[["top", "right"]].set_visible(False)
    if fc:
        ax.text(0.02, -0.24, f"The ordering survives at the extremes; the magnitude "
                             f"does not. Simulated foot/chest {fc:.0f}× vs measured 1.7×.",
                transform=ax.transAxes, fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return stats_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run_ra1e7")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or args.run
    os.makedirs(out, exist_ok=True)

    fig_body_model(os.path.join(out, "fig_body_model.png"))
    print("body model figure written")

    d, meta = load(args.run)

    # Illustrate a frame that actually contains shed cores. The last snapshot is
    # an arbitrary phase of the shedding cycle and is often the quietest one.
    dxm = meta["dx"]
    fluid0 = ~d["solid"]
    head0 = meta["body_y"] + 1.70 / dxm
    jtop0 = meta["ny"] - _sponge_cells(meta)
    j0g, i0g = np.mgrid[0:meta["ny"], 0:meta["nx"]]
    plume = (fluid0 & (j0g > head0) & (j0g < jtop0)
             & (np.abs(i0g - meta["body_x"]) < 0.55 / dxm))
    best, best_n = len(d["snap_t"]) - 1, -1
    for k in range(int(0.5 * len(d["snap_t"])), len(d["snap_t"]), 3):
        n = len(va.find_cores(d["u"][k], d["v"][k], mask=plume,
                              threshold=0.20, min_cells=8))
        if n > best_n:
            best, best_n = k, n
    print(f"illustrating frame {best} (t={d['snap_t'][best]:.1f} s, {best_n} cores)")

    fig_field(d, meta, os.path.join(out, "fig_field.png"), frame=best)
    fig_criteria(d, meta, os.path.join(out, "fig_criteria.png"), frame=best)
    res = fig_shedding(d, meta, os.path.join(out, "fig_shedding.png"))
    fig_boundary_layer(d, meta, os.path.join(out, "fig_boundary_layer.png"))
    r = fig_segment_conductance(meta, os.path.join(out, "fig_segment_conductance.png"))

    # Fit cores shed above the head, over the statistically steady part of the run.
    # The search is confined to the plume: the standing recirculation cells filling
    # the corners of the domain are large, slow and rotating, so they dominate a
    # whole-domain lambda_ci scan and crowd out the structures actually of interest.
    dx = meta["dx"]
    fluid = ~d["solid"]
    head = meta["body_y"] + 1.70 / dx
    j_top = meta["ny"] - _sponge_cells(meta)
    jj, ii = np.mgrid[0:meta["ny"], 0:meta["nx"]]
    region = (fluid & (jj > head) & (jj < j_top)
              & (np.abs(ii - meta["body_x"]) < 0.55 / dx))

    fits = []
    k0 = int(0.45 * len(d["snap_t"]))
    for k in range(k0, len(d["snap_t"])):
        u, v = d["u"][k], d["v"][k]
        for c in va.find_cores(u, v, mask=region, threshold=0.20, min_cells=8)[:4]:
            f = va.fit_core(u, v, c, radius=max(9, int(2.4 * np.sqrt(c["cells"]))),
                            mask=fluid)
            if f and f.get("oseen") and f["oseen"]["r2"] > 0.5:
                f["t"] = float(d["snap_t"][k])
                fits.append(f)
    # Refitting costs most of this script's runtime, so keep the results.
    with open(os.path.join(out, "core_fits.json"), "w") as fh:
        json.dump([{k: v for k, v in f.items() if k != "profile"} | {
            "profile": f.get("profile")} for f in fits], fh)
    rank_stats = _rankine_stats([f for f in fits
                                 if f.get("rankine") and f.get("oseen")]) if fits else {}
    ok = fig_rankine(fits, os.path.join(out, "fig_rankine.png"))

    summary = dict(
        n_cores_fitted=len(fits),
        rankine_r2_mean=float(np.mean([f["rankine"]["r2"] for f in fits])) if fits else None,
        oseen_r2_mean=float(np.mean([f["oseen"]["r2"] for f in fits])) if fits else None,
        core_radius_cells_mean=float(np.mean([f["oseen"]["a"] for f in fits])) if fits else None,
        core_radius_m_mean=float(np.mean([f["oseen"]["a"] for f in fits]) * meta["dx"]) if fits else None,
        rankine_vs_oseen=rank_stats,
        shedding_hz=res["freq_hz"] if res else None,
        shedding_period_s=res["period_s"] if res else None,
        shedding_period_peakcount_s=res.get("period_peakcount_s") if res else None,
        shedding_window_s=[res["t_start"], res["t_end"]] if res else None,
        strouhal=va.strouhal(res["freq_hz"], bm.BODY_HEIGHT, meta["u_ff"]) if res else None,
        strouhal_measured=va.strouhal(1 / 5.0, bm.BODY_HEIGHT, meta["u_ff"]),
        conductance=r,
        v_peak_over_uff=float(np.max(d["hist_vmax"]) / meta["u_lat"]),
        rankine_figure=bool(ok),
    )

    # The measured 0.2-0.3 m/s is the peak in the plume 0.5 m above the head, not the
    # global maximum, which sits in the near-wall jet and is a different quantity.
    if "hist_vrowmax" in d:
        rowmax = np.asarray(d["hist_vrowmax"])
        half = rowmax[len(rowmax) // 2:]
        summary.update(
            probe_peak_over_uff=float(half.mean() / meta["u_lat"]),
            probe_peak_ms=float(half.mean() / meta["u_lat"] * meta["u_ff"]),
            probe_peak_ms_max=float(half.max() / meta["u_lat"] * meta["u_ff"]),
            measured_probe_ms="0.2-0.3",
        )
    summary["diverged"] = meta.get("diverged")
    summary["seconds_completed"] = meta.get("seconds_completed")
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
