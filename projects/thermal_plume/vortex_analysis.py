"""Vortex identification and Rankine fitting on the simulated plume.

This mirrors the pipeline the PIV groups ran on real thermal plumes: compute the
velocity-gradient tensor, isolate rotating cores from the surrounding shear and
strain with the Q-criterion and the swirling-strength criterion (lambda_ci), then
fit an azimuthal velocity profile to each core.

Two profiles are fitted to every core, which is the point of the exercise:

    Rankine   v(r) = Gam/(2 pi a^2) r   for r < a,   Gam/(2 pi r)   for r > a
    Oseen     v(r) = Gam/(2 pi r) [1 - exp(-1.2564 r^2/a^2)]

Rankine has a corner at r = a. A real viscous core does not, so the comparison
between the two fits is a measurement of how far the idealisation is from the flow
rather than a decoration on it.
"""

import numpy as np
from scipy import ndimage, optimize, signal

OSEEN_C = 1.25643      # makes `a` the radius of peak azimuthal velocity


# ---------------------------------------------------------------------------
# Velocity-gradient invariants
# ---------------------------------------------------------------------------

def gradients(u, v, dx=1.0):
    dudy, dudx = np.gradient(u, dx)
    dvdy, dvdx = np.gradient(v, dx)
    return dudx, dudy, dvdx, dvdy


def vorticity(u, v, dx=1.0):
    dudx, dudy, dvdx, dvdy = gradients(u, v, dx)
    return dvdx - dudy


def q_criterion(u, v, dx=1.0):
    """Q = 0.5(|Omega|^2 - |S|^2). Positive where rotation beats strain."""
    dudx, dudy, dvdx, dvdy = gradients(u, v, dx)
    s12 = 0.5 * (dudy + dvdx)
    om12 = 0.5 * (dvdx - dudy)
    s_sq = dudx ** 2 + dvdy ** 2 + 2.0 * s12 ** 2
    om_sq = 2.0 * om12 ** 2
    return 0.5 * (om_sq - s_sq)


def swirling_strength(u, v, dx=1.0):
    """lambda_ci: imaginary part of the complex eigenvalue pair of grad(u).

    Non-zero only where a fluid particle traces a spiral, so unlike vorticity it
    does not light up in pure shear. For a 2x2 tensor the eigenvalues are complex
    when the discriminant of the characteristic polynomial is negative.
    """
    dudx, dudy, dvdx, dvdy = gradients(u, v, dx)
    tr = dudx + dvdy
    det = dudx * dvdy - dudy * dvdx
    disc = tr ** 2 - 4.0 * det
    return np.where(disc < 0.0, 0.5 * np.sqrt(np.maximum(-disc, 0.0)), 0.0)


# ---------------------------------------------------------------------------
# Core detection
# ---------------------------------------------------------------------------

def find_cores(u, v, dx=1.0, mask=None, threshold=0.25, min_cells=12):
    """Connected regions of strong swirling strength, one core each.

    `threshold` is a fraction of the field's peak lambda_ci. Returns a list of
    dicts with the core's centroid (row, col), circulation sign, and cell count.
    """
    lci = swirling_strength(u, v, dx)
    if mask is not None:
        lci = np.where(mask, lci, 0.0)
    peak = lci.max()
    if peak <= 0:
        return []
    labels, n = ndimage.label(lci > threshold * peak)
    vort = vorticity(u, v, dx)
    cores = []
    for lab in range(1, n + 1):
        sel = labels == lab
        cells = int(sel.sum())
        if cells < min_cells:
            continue
        weight = lci * sel
        total = weight.sum()
        rows, cols = np.nonzero(sel)
        rc = float((weight[sel] * rows).sum() / total)
        cc = float((weight[sel] * cols).sum() / total)
        cores.append(dict(row=rc, col=cc, cells=cells,
                          lci_peak=float(lci[sel].max()),
                          sign=float(np.sign(vort[sel].mean())),
                          mean_vorticity=float(vort[sel].mean())))
    cores.sort(key=lambda c: -c["lci_peak"])
    return cores


# ---------------------------------------------------------------------------
# Profile fitting
# ---------------------------------------------------------------------------

def rankine_profile(r, gamma, a):
    a = max(a, 1e-9)
    inner = gamma / (2.0 * np.pi * a ** 2) * r
    outer = gamma / (2.0 * np.pi * np.maximum(r, 1e-9))
    return np.where(r < a, inner, outer)


def oseen_profile(r, gamma, a):
    a = max(a, 1e-9)
    r = np.maximum(r, 1e-9)
    return gamma / (2.0 * np.pi * r) * (1.0 - np.exp(-OSEEN_C * r ** 2 / a ** 2))


def azimuthal_profile(u, v, row, col, radius, nbins=18, mask=None):
    """Azimuthally averaged tangential velocity around a core centre.

    The core's own translation is removed first: a vortex drifting upward in the
    plume would otherwise read as a spurious asymmetry in v_theta.
    """
    ny, nx = u.shape
    j0, j1 = max(0, int(row - radius) - 1), min(ny, int(row + radius) + 2)
    i0, i1 = max(0, int(col - radius) - 1), min(nx, int(col + radius) + 2)
    jj, ii = np.mgrid[j0:j1, i0:i1]
    dy, dxx = jj - row, ii - col
    r = np.hypot(dy, dxx)

    sub_u, sub_v = u[j0:j1, i0:i1], v[j0:j1, i0:i1]
    sub_mask = np.ones_like(r, dtype=bool) if mask is None else mask[j0:j1, i0:i1]

    core = sub_mask & (r <= max(1.5, 0.25 * radius))
    uc = sub_u[core].mean() if core.any() else 0.0
    vc = sub_v[core].mean() if core.any() else 0.0

    with np.errstate(invalid="ignore", divide="ignore"):
        tx, ty = -dy / r, dxx / r          # unit tangential vector, counter-clockwise
    v_theta = (sub_u - uc) * tx + (sub_v - vc) * ty

    good = sub_mask & (r > 0.5) & (r <= radius) & np.isfinite(v_theta)
    if good.sum() < nbins:
        return None
    edges = np.linspace(0.5, radius, nbins + 1)
    idx = np.digitize(r[good], edges) - 1
    rs, vs, ns = [], [], []
    for b in range(nbins):
        sel = idx == b
        if sel.sum() >= 3:
            rs.append(0.5 * (edges[b] + edges[b + 1]))
            vs.append(float(v_theta[good][sel].mean()))
            ns.append(int(sel.sum()))
    if len(rs) < 6:
        return None
    return np.array(rs), np.array(vs), np.array(ns)


def _r2(y, yhat):
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def fit_core(u, v, core, radius, mask=None):
    """Fit Rankine and Oseen profiles to one core. Returns None if unfittable."""
    prof = azimuthal_profile(u, v, core["row"], core["col"], radius, mask=mask)
    if prof is None:
        return None
    r, vt, _ = prof

    peak_i = int(np.argmax(np.abs(vt)))
    a0 = max(1.0, r[peak_i])
    gam0 = vt[peak_i] * 2.0 * np.pi * a0
    out = dict(core)

    for name, fn in (("rankine", rankine_profile), ("oseen", oseen_profile)):
        try:
            popt, _ = optimize.curve_fit(
                lambda rr, gam, a: fn(rr, gam, a), r, vt,
                p0=[gam0, a0], maxfev=20000,
                bounds=([-np.inf, 0.5], [np.inf, radius]))
            pred = fn(r, *popt)
            out[name] = dict(gamma=float(popt[0]), a=float(popt[1]),
                             r2=float(_r2(vt, pred)),
                             rmse=float(np.sqrt(((vt - pred) ** 2).mean())))
        except (RuntimeError, ValueError):
            out[name] = None

    out["profile"] = dict(r=r.tolist(), v_theta=vt.tolist())
    return out


# ---------------------------------------------------------------------------
# Shedding cadence
# ---------------------------------------------------------------------------

def shedding_frequency(t, sig, tmin=None, fmin=0.06, fmax=3.0):
    """Dominant frequency of a probe signal, via Welch PSD.

    `tmin` drops the spin-up transient. It matters more than it looks: run the
    whole record and the spectrum locks onto the one-off step where the plume
    stops growing and starts oscillating, reporting a period several times the
    real one. Restrict to the statistically steady window instead.

    A peak-count estimate is returned alongside as an independent cross-check,
    because the two disagreeing is the signal that the band is wrong.
    """
    t = np.asarray(t, dtype=float)
    sig = np.asarray(sig, dtype=float)
    if tmin is not None:
        sel = t >= tmin
        t, sig = t[sel], sig[sel]
    if t.size < 64:
        return None
    fs = 1.0 / np.median(np.diff(t))
    x = signal.detrend(sig, type="linear")
    # Remove any residual slow wander so it cannot outweigh the shedding line.
    win = max(3, int((1.0 / fmin) * fs) | 1)
    x = x - ndimage.uniform_filter1d(x, win)

    nper = min(len(x), max(128, len(x) // 3))
    freqs, psd = signal.welch(x, fs=fs, nperseg=nper)
    band = (freqs >= fmin) & (freqs <= min(fmax, 0.5 * fs))
    if not band.any():
        return None
    fb, pb = freqs[band], psd[band]

    # Take the strongest interior local maximum, not the global one. The plume also
    # meanders slowly as a whole, and that energy piles up against the low-frequency
    # edge of the band, where argmax lands on the boundary rather than on a line.
    loc, _ = signal.find_peaks(pb)
    if loc.size:
        fpk = float(fb[loc[np.argmax(pb[loc])]])
        f_slow = float(fb[np.argmax(pb)])
    else:
        fpk = float(fb[np.argmax(pb)])
        f_slow = fpk

    # Independent estimate: mean interval between prominent maxima.
    prom = 0.25 * (np.percentile(x, 95) - np.percentile(x, 5))
    peaks, _ = signal.find_peaks(x, prominence=max(prom, 1e-12),
                                 distance=max(1, int(0.25 * fs / fmax)))
    f_count = np.nan
    if peaks.size >= 3:
        f_count = 1.0 / float(np.median(np.diff(t[peaks])))

    return dict(freq_hz=fpk, period_s=1.0 / fpk if fpk > 0 else np.inf,
                freq_peakcount_hz=f_count,
                period_peakcount_s=1.0 / f_count if f_count > 0 else np.nan,
                freq_slow_hz=f_slow,
                period_slow_s=1.0 / f_slow if f_slow > 0 else np.nan,
                n_peaks=int(peaks.size), t_start=float(t[0]), t_end=float(t[-1]),
                freqs=freqs, psd=psd)


def strouhal(freq_hz, length, velocity):
    """St = f L / U -- the dimensionless form the cadence must be compared in."""
    return freq_hz * length / velocity
