#!/usr/bin/env python3
"""The eye's spectral window, for the "window" section of web/indirect-perception.html.

Measured inputs (data/, downloaded from the Colour & Vision Research Laboratory, cvrl.org):
  linss2_10e_1.csv     Stockman & Sharpe (2000) 2-deg cone fundamentals, linear energy, 390-830 nm
  lensss_1.csv         Stockman & Sharpe lens optical density, 390-830 nm (32-year-old standard)
  lensws.csv           Wyszecki & Stiles (1982) lens density, 360-710 nm
  macss_1.csv          Stockman & Sharpe 2-deg macular pigment optical density, 390-830 nm
  linCIE2008v2e_1.csv  CIE (2008) 2-deg photopic luminosity function V(lambda), 390-830 nm

Outside 390-830 nm there are no fundamentals, so the photopigments are modelled with the
Govardovskii et al. (2000) A1 template (alpha + beta bands), self-screened at the Stockman &
Sharpe peak densities. This script fits each template's lambda-max to the measured
fundamentals (after removing lens and macular filtering), reports how well it fits, and
writes window_data.js: the arrays the page draws, so the chart is computed from measured
spectra plus one fitted, checked model.

Run from this directory: python spectral_window.py"""

import json
import os

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "window_data.js")
PAGE = os.path.join(HERE, "..", "..", "web", "indirect-perception.html")
PEAK_OD = {"L": 0.50, "M": 0.50, "S": 0.40}      # Stockman & Sharpe (2000), 2-deg


def load(name, cols=2):
    a = np.genfromtxt(os.path.join(DATA, name), delimiter=",")
    a = a[~np.isnan(a[:, 0])]
    return a[:, :cols] if cols else a


def govardovskii(lam, lmax):
    """A1 photopigment absorbance template, normalised to 1 at lmax (Govardovskii et al. 2000)."""
    x = lmax / lam
    a = 0.8795 + 0.0459 * np.exp(-(lmax - 300.0) ** 2 / 11940.0)
    alpha = 1.0 / (np.exp(69.7 * (a - x)) + np.exp(28.0 * (0.922 - x))
                   + np.exp(-14.9 * (1.104 - x)) + 0.674)
    lmb = 189.0 + 0.315 * lmax
    bb = -40.5 + 0.195 * lmax
    beta = 0.26 * np.exp(-((lam - lmb) / bb) ** 2)
    return alpha + beta


def absorptance(lam, lmax, od):
    a = 1.0 - 10.0 ** (-od * govardovskii(lam, lmax))
    return a / a.max() if np.ndim(lam) and len(np.atleast_1d(lam)) > 1 else a


def lens_density(lam):
    """Stockman & Sharpe lens density, extended below 390 nm with the Wyszecki & Stiles
    template scaled to meet it at 390 nm; zero beyond the data."""
    ss = load("lensss_1.csv")
    ws = load("lensws.csv")
    k = np.interp(390.0, ss[:, 0], ss[:, 1]) / np.interp(390.0, ws[:, 0], ws[:, 1])
    d = np.interp(lam, ss[:, 0], ss[:, 1], right=0.0)
    lo = lam < 390.0
    d[lo] = k * np.interp(lam[lo], ws[:, 0], ws[:, 1])
    return d, k


def water_alpha():
    """Absorption coefficient of liquid water (1/m) from Segelstein (1981), via the CC0
    refractiveindex.info database: alpha = 4 pi k / lambda. Returns (lambda_um, alpha)."""
    import re
    text = open(os.path.join(DATA, "water_segelstein.yml"), encoding="utf-8").read()
    rows = [l.split() for l in text.split("data: |", 1)[1].splitlines()
            if re.match(r"^\s+[0-9.]+E?[-+0-9]*\s", l)]
    a = np.array([[float(x) for x in r[:3]] for r in rows if len(r) >= 3])
    return a[:, 0], 4.0 * np.pi * a[:, 2] / (a[:, 0] * 1e-6)


def main():
    fund = load("linss2_10e_1.csv", cols=0)
    lam = fund[:, 0]
    dlens, k = lens_density(lam)
    dmac = np.interp(lam, *load("macss_1.csv").T, right=0.0)
    media = 10.0 ** (-(dlens + dmac))
    report = {}
    fits = {}
    for i, cone in enumerate("LMS"):
        e = np.nan_to_num(fund[:, i + 1])
        ok = e > 0
        q = e / lam                                   # energy -> quantal sensitivity
        q = q / q.max()
        target = np.log10(q[ok])

        def err(lmax):
            a = 1.0 - 10.0 ** (-PEAK_OD[cone] * govardovskii(lam, lmax))
            p = a * media
            p = p / p.max()
            return np.sqrt(np.mean((np.log10(p[ok]) - target) ** 2))

        res = minimize_scalar(err, bounds=(400.0, 600.0), method="bounded")
        a = 1.0 - 10.0 ** (-PEAK_OD[cone] * govardovskii(lam, res.x))
        p = a * media
        p /= p.max()
        dev = np.abs(np.log10(p[ok]) - target)
        vis = ok & (lam <= 700)
        fits[cone] = float(res.x)
        report[cone] = dict(lmax=round(res.x, 1), rms_log=round(res.fun, 3),
                            max_log_390_700=round(float(np.max(np.abs(
                                np.log10(p[vis]) - np.log10(q[vis])))), 3),
                            max_log_all=round(float(dev.max()), 3),
                            reach=f"{lam[ok].min():.0f}-{lam[ok].max():.0f} nm")
    for c, r in report.items():
        print(c, r)
    print(f"lens: Wyszecki & Stiles scaled x{k:.2f} below 390 nm")

    # arrays for the page: 5 nm from 300 to 1100 nm
    grid = np.arange(300, 1101, 5, dtype=float)
    dl, _ = lens_density(grid)
    dl[grid < 360] = np.nan                           # no lens data: drawn as opaque
    dm = np.interp(grid, *load("macss_1.csv").T, left=np.nan, right=0.0)
    vl = load("linCIE2008v2e_1.csv")
    v = np.interp(grid, vl[:, 0], vl[:, 1] / vl[:, 0], left=np.nan, right=np.nan)
    v = v / np.nanmax(v)                              # quantal V(lambda)
    pig = {c: (1.0 - 10.0 ** (-PEAK_OD[c] * govardovskii(grid, fits[c]))) for c in "LMS"}
    for c in pig:
        pig[c] /= pig[c].max()
    r = lambda a, n=4: [None if not np.isfinite(x) else float(f"{x:.{n}g}") for x in a]
    # measured cone fundamentals (energy units, each scaled to its peak) for the metamer demo
    f5 = fund[(fund[:, 0] % 5 == 0)]
    fn = np.nan_to_num(f5[:, 1:4])
    fn = fn / fn.max(axis=0)
    # water absorption on a log-wavelength grid, 0.2-30 um, for the "why the window is here" chart
    wl, wa = water_alpha()
    wgrid = np.logspace(np.log10(0.2), np.log10(30.0), 300)
    wlog = np.interp(np.log10(wgrid), np.log10(wl), np.log10(wa))
    k = int(np.argmin(np.where((wgrid > 0.3) & (wgrid < 0.8), wlog, 99)))
    print(f"water: clearest at {1000 * wgrid[k]:.0f} nm ({10 ** wlog[k]:.4f} /m); "
          f"at 9.35 um alpha {10 ** np.interp(np.log10(9.35), np.log10(wl), np.log10(wa)):.3g} /m")
    data = dict(lam0=300, step=5, lens=r(dl), mac=r(dm), V=r(v),
                wl0=0.2, wl1=30.0, wa=r(wlog),
                f0=int(f5[0, 0]), fL=r(fn[:, 0]), fM=r(fn[:, 1]), fS=r(fn[:, 2]),
                L=r(pig["L"]), M=r(pig["M"]), S=r(pig["S"]),
                lmax={c: round(fits[c], 1) for c in "LMS"},
                fit={c: report[c]["rms_log"] for c in "LMS"})
    js = "window.EYE_WINDOW = " + json.dumps(data, separators=(",", ":")) + ";"
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("// generated by projects/indirect-perception/spectral_window.py; do not edit\n")
        f.write(js + "\n")
    print(f"wrote {OUT} ({os.path.getsize(OUT) // 1024} KB)")

    # inline the same data into the page, between its markers
    with open(PAGE, encoding="utf-8") as f:
        page = f.read()
    a, b = "/* EYE_WINDOW:BEGIN */", "/* EYE_WINDOW:END */"
    i, j = page.index(a) + len(a), page.index(b)
    page = page[:i] + "\n" + js + "\n" + page[j:]
    with open(PAGE, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"inlined into {os.path.relpath(PAGE, HERE)}")
    mac = load("macss_1.csv")
    k = int(np.argmax(mac[:, 1]))
    print(f"macular peak {mac[k, 1]:.3f} at {mac[k, 0]:.0f} nm -> passes {100 * 10 ** -mac[k, 1]:.0f}%")


if __name__ == "__main__":
    main()
