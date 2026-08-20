#!/usr/bin/env python3
"""What patterns this dish should produce, from the dish's own dimensions.

Standing waves on a shallow water layer in a circular basin with vertical
walls. No normal flow at the wall gives a Neumann boundary, so the surface
modes are

    eta(r, th) = J_m(k r) cos(m th),     J_m'(k a) = 0  ->  k = j'_{m,n} / a

with a = the dish's inner radius. Gravity-capillary dispersion at depth h:

    omega^2 = (g k + sigma k^3 / rho) * tanh(k h)

Vertical shaking drives this parametrically (Faraday), so the pattern responds
at HALF the drive frequency: a 60 Hz tone paints a 30 Hz mode.

The drive frequency fixes k, and therefore the pattern's SCALE, tightly. It
does not fix the angular order m -- a whole family of (m, n) share nearly the
same k, and which one you actually see is picked by nonlinear selection and by
how symmetric your rig is. The panels below show representative members.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import jn, jnp_zeros

import cymatics_dish as C

OUT = "cymatics_patterns.png"
G, SIGMA, RHO = 9.81, 0.0728, 1000.0        # SI
A = C.R_IN / 1000.0                          # basin radius, m
DRIVE = [15.0, 30.0, 55.0, 90.0]             # drive frequencies to render, Hz
INK = "#1b1b1f"
ACC = "#2a6ea8"


def omega2(k, h):
    return (G * k + SIGMA * k**3 / RHO) * np.tanh(k * h)


def wavenumber(f, h):
    """k for a surface mode at frequency f (Hz) on depth h (m)."""
    w2 = (2 * np.pi * f)**2
    return brentq(lambda k: omega2(k, h) - w2, 1e-3, 1e5)


def cellular_m(f, h):
    """Angular order giving near-isotropic cells rather than concentric rings.

    A mode J_m(kr)cos(m th) is locally a plane wave with radial wavenumber
    sqrt(k^2 - m^2/r^2) and azimuthal m/r. Those match -- square cells -- when
    m = k r / sqrt(2); taken at r = 0.7a that is m ~ 0.5 k a. Faraday selection
    in a real dish lands in this neighbourhood, which is why cymatics photos
    show cells and stars, not bullseyes."""
    return max(1, int(round(0.5 * wavenumber(f, h) * A)))


def nearest_mode(f, h, m):
    """(j', n) whose eigenfrequency is closest to f for angular order m."""
    target = wavenumber(f, h) * A
    zeros = jnp_zeros(m, max(2, int(target / np.pi) + 3))
    n = int(np.argmin(np.abs(zeros - target)))
    return zeros[n], n + 1


def main():
    h = C.WATER / 1000.0
    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.5, 1.0], hspace=.30, wspace=.14)

    # --- top row: representative patterns ---------------------------------
    g = np.linspace(-A, A, 900)
    X, Y = np.meshgrid(g, g)
    R, TH = np.hypot(X, Y), np.arctan2(Y, X)
    for i, fd in enumerate(DRIVE):
        fp = fd / 2.0                                   # Faraday subharmonic
        m = cellular_m(fp, h)
        jz, n = nearest_mode(fp, h, m)
        eta = jn(m, jz * R / A) * np.cos(m * TH)
        eta = np.where(R <= A, eta / np.abs(eta[R <= A]).max(), np.nan)

        ax = fig.add_subplot(gs[0, i])
        ax.imshow(eta.T, extent=[-A, A, -A, A], origin="lower",
                  cmap="RdBu_r", vmin=-1, vmax=1, interpolation="bilinear")
        ax.add_patch(plt.Circle((0, 0), A, fill=False, color=INK, lw=1.4))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{fd:.0f} Hz drive → {fp:.1f} Hz pattern\n"
                     f"m={m}, n={n}   λ ≈ {2000 * np.pi / (jz / A):.1f} mm",
                     fontsize=9.5, pad=6)

    # --- bottom: scale vs drive frequency, one curve per depth ------------
    ax = fig.add_subplot(gs[1, :])
    fr = np.linspace(8, 200, 500)
    for hh, ls in ((0.002, "--"), (0.003, "-"), (0.005, ":")):
        lam = np.array([2000 * np.pi / wavenumber(f / 2, hh) for f in fr])
        ax.plot(fr, lam, ls, color=INK, lw=1.7, label=f"{hh * 1000:.0f} mm water")
    for fd in DRIVE:
        ax.axvline(fd, color=ACC, lw=1.0, alpha=.5, ls="-")
        ax.text(fd, 62, f"{fd:.0f}", color=ACC, fontsize=8.5, ha="center")
    ax.set_xlabel("drive frequency (Hz)")
    ax.set_ylabel("pattern wavelength (mm)")
    ax.set_yscale("log")
    ax.set_xlim(8, 200)
    ax.set_ylim(3, 80)
    ax.grid(alpha=.2, which="both")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title("Drive frequency sets the pattern scale; depth barely matters "
                 "above ~40 Hz (capillary takes over from gravity)", fontsize=10)

    tw = ax.twinx()                                   # same curve, useful units
    tw.set_yscale("log")
    tw.set_ylim(2 * C.R_IN / 3.0, 2 * C.R_IN / 80.0)
    tw.set_ylabel("cells across the dish")
    tw.minorticks_off()
    tw.set_yticks([5, 10, 20, 40, 80])
    tw.set_yticklabels(["5", "10", "20", "40", "80"])

    fig.text(0.5, 0.455,
             "Single eigenmodes are shown. A real dish rings several degenerate "
             "modes at once, which fills in the quiet centre and rotates the "
             "cells into stars and polygons —\nthe cell SIZE is the robust "
             "prediction, the exact figure is not. Patterns are subharmonic: a "
             "60 Hz tone paints a 30 Hz pattern.",
             ha="center", va="center", fontsize=8.8, color="#555")

    fig.suptitle(f"Predicted cymatic patterns — Ø{2 * C.R_IN:.0f} mm water "
                 f"surface, {C.WATER:.0f} mm deep, vertical (Faraday) drive",
                 fontsize=12.5)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"wrote {OUT}")

    print(f"basin Ø{2 * C.R_IN:.0f} mm, depth {C.WATER:.0f} mm")
    for fd in DRIVE:
        m = cellular_m(fd / 2, h)
        jz, n = nearest_mode(fd / 2, h, m)
        k = jz / A
        print(f"  {fd:5.0f} Hz drive -> {fd / 2:5.1f} Hz mode  m={m:3d} n={n:3d}  "
              f"k={k:7.1f} 1/m  lambda={2000 * np.pi / k:5.1f} mm  "
              f"{2 * C.R_IN / (1000 * np.pi / k):.0f} cells across")


if __name__ == "__main__":
    main()
