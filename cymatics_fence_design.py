#!/usr/bin/env python3
"""Size the perforated wave absorber (Jarlan fence) from the wave physics.

A perforated fence standing in the water with a narrow chamber behind it, closed
by the dish's rim. Waves push flow through the slots; the sudden contraction and
re-expansion sheds vortices and burns energy. The chamber sets the phase of what
comes back out.

MODEL
-----
Locally 1D (the chamber is millimetres across a >130 mm radius, so curvature is
irrelevant). Work in (p, q): p = pressure at the surface, q = depth-integrated
flux. For a travelling capillary-gravity wave on depth h,

    omega^2 = (g + sigma k^2 / rho) k tanh(k h)          dispersion
    q       = (omega / k) * eta                          continuity
    p       = (rho g + sigma k^2) eta                    restoring pressure
    => Z_c  = p / q = rho * omega / tanh(k h)            characteristic impedance

The chamber is a transmission line of length B terminated by a rigid wall, so
looking outward from the fence it presents

    Z_ch = -i * Z_c * cot(k B)

The fence adds a series impedance: an inertance from the slug of water in each
slot, and a resistance from the jet loss. The jet loss is quadratic in velocity,
so it is Lorentz-linearised at a representative wave amplitude -- absorption is
genuinely amplitude-dependent and this is where the honest uncertainty lives.

    u_slot = q / (eps * h)
    Z_w    = rho * (8/(3*pi)) * K * U / (eps*h)  +  i*omega*rho*l_eff/(eps*h)
    K      = (1/(Cc*eps))^2 - 1                  sharp-edged slot, Cc ~ 0.6
    l_eff  = d + 2*0.4*w                         slot length plus end corrections

Then R = (Z_in - Z_c)/(Z_in + Z_c) with Z_in = Z_w + Z_ch, and the fraction of
incident energy absorbed is 1 - |R|^2.

The comparison case -- the plain vertical wall the current dish has -- is
R = 1 exactly: nothing absorbed at any frequency.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

OUT = "cymatics_fence_design.png"
G, SIGMA, RHO, NU = 9.81, 0.0728, 1000.0, 1.0e-6
CC = 0.6                      # contraction coefficient, sharp-edged slot
H = 3.0e-3                    # water depth (m)
AMP = 0.3e-3                  # representative wave amplitude (m)
BAND = (10.0, 120.0)          # wave frequencies to care about (Hz)
INK, ACC, WARM = "#1b1b1f", "#2a6ea8", "#b04a2a"


def wavenumber(f, h=H):
    w2 = (2 * np.pi * f) ** 2
    return brentq(lambda k: (G * k + SIGMA * k**3 / RHO) * np.tanh(k * h) - w2,
                  1e-3, 1e5)


def absorption(k, omega, eps, B, d, w, h=H, amp=AMP):
    """Fraction of incident wave energy absorbed. k, omega may be arrays."""
    Zc = RHO * omega / np.tanh(k * h)

    t = np.tan(k * B)
    Zch = -1j * Zc / np.where(np.abs(t) < 1e-9, 1e-9, t)

    q = (omega / k) * amp                      # flux amplitude of the wave
    U = q / (eps * h)                          # slot velocity amplitude
    K = (1.0 / (CC * eps)) ** 2 - 1.0
    l_eff = d + 0.8 * w                        # slug length plus end corrections
    Zw = RHO * (8 / (3 * np.pi)) * K * U / (eps * h) \
        + 1j * omega * RHO * l_eff / (eps * h)

    R = (Zw + Zch - Zc) / (Zw + Zch + Zc)
    return 1.0 - np.abs(R) ** 2


FS = np.linspace(*BAND, 48)
KS = np.array([wavenumber(f) for f in FS])
WS = 2 * np.pi * FS


def band_score(eps, B, d, w):
    """Worst-case absorption across the band -- optimise the weakest point, not
    the average, so the fence has no dead frequency."""
    a = absorption(KS, WS, eps, B, d, w)
    return a.min(), a.mean()


def main():
    print(f"depth {H * 1000:.0f} mm, design amplitude {AMP * 1000:.1f} mm, "
          f"band {BAND[0]:.0f}-{BAND[1]:.0f} Hz")
    print("sweeping porosity, chamber, fence thickness and slot width "
          "for the worst point in the band\n")

    epss = np.arange(0.15, 0.86, 0.025)
    Bs = np.arange(0.5, 9.1, 0.25) * 1e-3
    ds = np.array([0.6, 0.8, 1.0, 1.2, 1.5]) * 1e-3      # <= 1.5 mm as asked
    ws = np.array([0.8, 1.2, 1.6, 2.0, 2.5, 3.0]) * 1e-3
    best = None
    for d in ds:
        for w in ws:
            for e in epss:
                for B in Bs:
                    lo, avg = band_score(e, B, d, w)
                    if best is None or lo > best[0]:
                        best = (lo, avg, e, B, d, w)
    lo, avg, eps_o, B_o, d_o, w_o = best
    print(f"unconstrained best: porosity {eps_o:.3f}, chamber {B_o * 1000:.2f} mm, "
          f"fence {d_o * 1000:.1f} mm, slots {w_o * 1000:.1f} mm")
    print(f"      -> worst {lo * 100:.0f}%, mean {avg * 100:.0f}% absorbed")

    # what actually got built: the fence rounded up to 2 perimeters at a 0.4
    # nozzle, which costs almost nothing in performance
    eps_b, B_b, d, w = 0.45, 1.0e-3, 0.8e-3, 1.2e-3
    lo, avg = band_score(eps_b, B_b, d, w)
    print(f"AS BUILT: porosity {eps_b:.3f}, chamber {B_b * 1000:.2f} mm, "
          f"fence {d * 1000:.1f} mm, slots {w * 1000:.1f} mm")
    print(f"      -> worst {lo * 100:.0f}%, mean {avg * 100:.0f}% absorbed")
    print(f"      slot pitch = w/eps = {w / eps_b * 1000:.2f} mm")
    print(f"      slug length / (eps*h) = "
          f"{(d + 0.8 * w) / (eps_b * H):.2f}  (the inertance ratio — this is "
          f"what limits the fence)")
    print(f"      radial cost = fence + chamber = "
          f"{(d + B_b) * 1000:.1f} mm per side\n")

    # re-grid at the chosen d, w so the map is comparable
    grid = np.zeros((len(epss), len(Bs)))
    for i, e in enumerate(epss):
        for j, B in enumerate(Bs):
            grid[i, j] = band_score(e, B, d, w)[0]

    fs = np.linspace(*BAND, 240)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.0),
                                  gridspec_kw={"width_ratios": [1.25, 1]})

    for e, B, ls, lab in ((eps_b, B_b, "-",
                           f"AS BUILT: ε={eps_b:.2f}, B={B_b*1000:.1f} mm"),
                          (0.15, B_b, "--", f"ε=0.15 (too tight), B={B_b*1000:.1f} mm"),
                          (0.80, B_b, ":", f"ε=0.80 (too open), B={B_b*1000:.1f} mm"),
                          (eps_b, 2.0e-3, "-.", f"ε={eps_b:.2f}, B=2.0 mm — the cliff")):
        a = absorption(np.array([wavenumber(f) for f in fs]),
                       2 * np.pi * fs, e, B, d, w)
        ax.plot(fs, np.array(a) * 100, ls, color=INK if ls == "-" else "#7a7a80",
                lw=2.2 if ls == "-" else 1.4, label=lab)
    ax.axhline(0, color=WARM, lw=2, label="plain vertical wall (today)")
    ax.set_xlabel("wave frequency (Hz)")
    ax.set_ylabel("energy absorbed (%)")
    ax.set_ylim(-4, 100)
    ax.set_xlim(*BAND)
    ax.grid(alpha=.2)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.set_title(f"Absorption vs frequency — {d*1000:.1f} mm fence, "
                 f"{w*1000:.1f} mm slots, {H*1000:.0f} mm water", fontsize=10.5)

    im = ax2.pcolormesh(Bs * 1000, epss, grid * 100, cmap="viridis",
                        vmin=0, vmax=100, shading="auto")
    ax2.plot(B_b * 1000, eps_b, "o", ms=9, mfc="none", mec="w", mew=2)
    ax2.set_xlabel("chamber width B (mm)")
    ax2.set_ylabel("porosity ε")
    ax2.set_title("Worst-case absorption across the band", fontsize=10.5)
    fig.colorbar(im, ax=ax2, label="% absorbed at the weakest frequency")

    fig.suptitle("Perforated wave absorber — sizing from the capillary-gravity "
                 "dispersion", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")

    print("\nabsorption at the frequencies that matter:")
    print("  drive   pattern(f/2)   lambda    absorbed at f    absorbed at f/2")
    for fd in (20.0, 30.0, 55.0, 90.0, 120.0):
        lam = 2 * np.pi / wavenumber(fd) * 1000
        print(f"  {fd:5.0f} Hz  {fd/2:6.1f} Hz    {lam:5.1f} mm     "
              f"{absorption(wavenumber(fd), 2*np.pi*fd, eps_b, B_b, d, w)*100:4.0f}%            "
              f"{absorption(wavenumber(fd/2), np.pi*fd, eps_b, B_b, d, w)*100:4.0f}%")
    return eps_b, B_b, d, w


if __name__ == "__main__":
    main()
