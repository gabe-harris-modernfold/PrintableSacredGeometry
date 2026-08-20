#!/usr/bin/env python3
"""faraday_disc.py -- 200 mm tactile Faraday-instability teaching disc.

Three concentric zones, each an axisymmetric Bessel standing wave frozen at
maximum displacement, driven at 30 / 60 / 90 Hz.  Faraday response is
subharmonic, so the surface waves are at 15 / 30 / 45 Hz; wavelengths come from
the gravity-capillary dispersion relation for deep water.

Zone amplitudes are slope-matched: zone A is pinned at 10 mm peak-to-valley and
B and C are scaled so every zone has the same maximum surface slope, which keeps
the outer rings from becoming blade-thin ridges a fingertip cannot enter.

Rims carry raised text: drive frequency, wave frequency, actual peak-to-valley.

Output: faraday_disc.stl plus preview PNGs.  Run from the repo root.
"""

import gc
import math
import struct

import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq, minimize_scalar

import matplotlib
matplotlib.use("Agg")
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path
from matplotlib.textpath import TextPath

# --------------------------------------------------------------------------
# physics
# --------------------------------------------------------------------------

G_ACCEL = 9.81          # m/s^2
SIGMA = 0.072           # N/m, clean water/air
RHO = 1000.0            # kg/m^3


def wavelength_mm(f_wave_hz):
    """Gravity-capillary wavelength (mm) for a deep-water surface wave."""
    omega2 = (2.0 * math.pi * f_wave_hz) ** 2

    def residual(k):                      # k in rad/m
        return G_ACCEL * k + (SIGMA / RHO) * k ** 3 - omega2

    k = brentq(residual, 1e-6, 1e6, xtol=1e-12, rtol=1e-14)
    return 2.0 * math.pi / k * 1000.0


# --------------------------------------------------------------------------
# disc parameters
# --------------------------------------------------------------------------

R_DISC = 100.0          # mm radius
ZONE_A_PV = 10.0        # mm peak-to-valley, zone A (user spec)
BASE_SLAB = 4.0         # mm of solid material under the deepest trough
# The rims are wayfinding, not the subject.  Keep them low and narrow enough
# that the wave field stays the dominant feature; they only have to be tall
# enough to read as a boundary against the tallest crest either side.
RIM_RISE = 4.0          # mm rim top above the wave midplane
RIM_RAMP = 0.7          # mm radial ramp from midplane up to rim top
RIM_W = 4.6             # mm minimum rim width; the J0-zero snapping below can
                        # widen any given rim by up to half a wavelength
RIM3_MIN = 5.5          # mm floor for the outer rim, which carries two strings
EDGE_CHAMFER = 1.0      # mm radial chamfer at the outer edge
EDGE_DROP = 1.2         # mm the chamfer falls from the rim top

TEXT_EMBOSS = 0.5       # mm raised text
MASK_PX = 0.02          # mm raster pitch for the glyph masks

GRID = 0.22             # mm target triangle edge on the top surface

DRIVE_HZ = (30.0, 60.0, 90.0)


LAMBDA = [wavelength_mm(f / 2.0) for f in DRIVE_HZ]
KWAVE = [2.0 * math.pi / lam for lam in LAMBDA]
J0_ZEROS = jn_zeros(0, 80)

# Every zone is the same axisymmetric mode, A * J0(k r) -- only k differs.  Zone
# A runs from the centre to a zero of J0, so it ends exactly on the midplane;
# the annuli are cut at zeros of their own J0 for the same reason.  All six zone
# edges therefore sit at the midplane and the rims meet them without a step.
J0_MIN = j0(jn_zeros(1, 1)[0])                     # first J0 minimum, -0.4028
AMP_A = ZONE_A_PV / (1.0 - J0_MIN)                 # J0 spans +1 .. -0.4028
_J0_MAX_R = 0.0                                    # zone A's crest is at r = 0

# Slope match: max |dh/dr| of A*J0(k r) is A*k*max|J1|.  Zone A sets the target;
# B and C are scaled to the same steepest slope, so no zone turns into
# blade-thin ridges a fingertip cannot enter.
_J1_PEAK = -minimize_scalar(lambda x: -abs(j1(x)), bounds=(1.0, 3.0),
                            method="bounded").fun
MAX_SLOPE = AMP_A * KWAVE[0] * _J1_PEAK


def build_layout():
    """Radial boundaries, every zone cut at zeros of its own J0.

    Two things are being traded off.  Ring count: counting rings alone would
    hand the radius to zone C, whose wavelength is shortest, and leave zone B
    with barely a full period to feel -- so score the *weakest* zone.  And rim
    width: a rim runs from one zone's last zero to the next zone's first zero at
    least RIM_W away, so snapping can inflate it by up to half a wavelength.
    Zone A's terminating zero is a free parameter that shifts every downstream
    boundary, so search it too and prefer layouts whose widest rim is slimmest.
    """
    ra, rb, rc = (J0_ZEROS / k for k in KWAVE)
    best = None
    for i_a in range(2, 7):
        for n_b in range(4, 14):
            i_b = int(np.searchsorted(rb, ra[i_a] + RIM_W))
            if i_b + n_b >= len(rb):
                break
            i_c = int(np.searchsorted(rc, rb[i_b + n_b] + RIM_W))
            # -1: searchsorted lands on the first zero *past* the limit
            n_c = int(np.searchsorted(rc, R_DISC - RIM3_MIN)) - 1 - i_c
            if n_c < 4 or i_c + n_c >= len(rc):
                continue
            widest = max(rb[i_b] - ra[i_a], rc[i_c] - rb[i_b + n_b],
                         R_DISC - rc[i_c + n_c])
            key = (min(i_a + 1, n_b, n_c), -round(widest, 3), n_b + n_c)
            if best is None or key > best[0]:
                best = (key, i_a, i_b, n_b, i_c, n_c)
    _, i_a, i_b, n_b, i_c, n_c = best

    zones = [(0.0, ra[i_a]), (rb[i_b], rb[i_b + n_b]),
             (rc[i_c], rc[i_c + n_c])]
    rims = [(zones[0][1], zones[1][0]), (zones[1][1], zones[2][0]),
            (zones[2][1], R_DISC)]
    return zones, rims


ZONES, RIMS = build_layout()
R_ZONE_A = ZONES[0][1]


def _count_rings(i):
    """Crests plus troughs in a zone -- the features a finger actually counts."""
    r = np.linspace(*ZONES[i], 40001)
    d = np.sign(np.diff(j0(KWAVE[i] * r)))
    # zone A also has the central antinode sitting on its r = 0 boundary
    return int((np.diff(d) != 0).sum()) + (1 if ZONES[i][0] == 0.0 else 0)


RINGS = [_count_rings(i) for i in range(3)]


def _zone_amp(i):
    """Scale a J0 annulus so its steepest slope matches zone A's."""
    r = np.linspace(*ZONES[i], 40001)
    return MAX_SLOPE / np.abs(KWAVE[i] * j1(KWAVE[i] * r)).max()


AMP = [AMP_A, _zone_amp(1), _zone_amp(2)]


def _zone_extremes(i):
    r = np.linspace(*ZONES[i], 40001)
    h = AMP[i] * j0(KWAVE[i] * r)
    return h.max() - h.min(), h.min()


ZONE_PV = [_zone_extremes(i)[0] for i in range(3)]

H_MIN = min(_zone_extremes(i)[1] for i in range(3))  # deepest point
Z_MID = BASE_SLAB - H_MIN                 # midplane height above the bed
Z_EDGE = Z_MID + RIM_RISE - EDGE_DROP     # top surface height at r = R_DISC


def smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def surface_profile(r):
    """Height above the midplane of the axisymmetric surface, at radius r."""
    r = np.atleast_1d(np.asarray(r, float))
    h = np.zeros_like(r)

    for idx, (z0, z1) in enumerate(ZONES):
        m = (r >= z0) & (r <= z1) if idx == 0 else (r > z0) & (r <= z1)
        h[m] = AMP[idx] * j0(KWAVE[idx] * r[m])

    for i, (a, b) in enumerate(RIMS):
        outer = (i == len(RIMS) - 1)
        m = (r > a) & (r <= b)
        if not m.any():
            continue
        rr = r[m]
        up = RIM_RISE * smoothstep((rr - a) / RIM_RAMP)
        if outer:
            down = RIM_RISE - EDGE_DROP * (
                1.0 - smoothstep((b - rr) / EDGE_CHAMFER))
        else:
            down = RIM_RISE * smoothstep((b - rr) / RIM_RAMP)
        h[m] = np.minimum(up, down)
    return h


# --------------------------------------------------------------------------
# rim decoration: raised text
# --------------------------------------------------------------------------

_FONT = FontProperties(family="DejaVu Sans", weight="bold")
_BIG = 1000.0           # build glyphs large, then scale: keeps curve flattening
                        # fine relative to a 2-3 mm cap height
_CAP_UNITS = TextPath((0, 0), "H", size=_BIG, prop=_FONT).get_extents().height


def glyph_polys(s, cap_mm):
    """Glyph outlines in mm, left edge at u = 0 and the inked bounding box
    centred on v = 0 so descenders do not drift off the rim.

    Kept as separate contours on purpose.  Path.contains_points fills the
    counters of a compound TextPath -- '0' comes out a solid blob -- and its
    `radius` argument is not a Minkowski offset either: one sign erodes the
    stems, the other floods the counters.  Rasterising contour by contour and
    XOR-ing gives the even-odd fill the glyphs actually need."""
    tp = TextPath((0, 0), s, size=_BIG, prop=_FONT)
    scale, bb = cap_mm / _CAP_UNITS, tp.get_extents()
    out = []
    for p in tp.to_polygons():
        q = np.asarray(p, float) * scale
        q[:, 0] -= bb.x0 * scale
        q[:, 1] -= 0.5 * (bb.y0 + bb.y1) * scale
        out.append(q)
    return out


def text_metrics(s, cap_mm):
    """Arc width and half-height of the inked area, in mm."""
    a = np.vstack(glyph_polys(s, cap_mm))
    return a[:, 0].max() - a[:, 0].min(), 0.5 * (a[:, 1].max() - a[:, 1].min())


def text_mask(s, cap_mm):
    """Binary even-odd raster of the string, plus its origin in mm."""
    polys = glyph_polys(s, cap_mm)
    a = np.vstack(polys)
    x0, y0 = a[:, 0].min() - 0.1, a[:, 1].min() - 0.1
    nx = int(np.ceil((a[:, 0].max() + 0.1 - x0) / MASK_PX))
    ny = int(np.ceil((a[:, 1].max() + 0.1 - y0) / MASK_PX))
    m = np.zeros((ny, nx), bool)
    for q in polys:                       # bbox-limited: one glyph at a time
        i0 = max(0, int((q[:, 0].min() - x0) / MASK_PX) - 1)
        i1 = min(nx, int((q[:, 0].max() - x0) / MASK_PX) + 2)
        j0 = max(0, int((q[:, 1].min() - y0) / MASK_PX) - 1)
        j1 = min(ny, int((q[:, 1].max() - y0) / MASK_PX) + 2)
        gx, gy = np.meshgrid(x0 + (np.arange(i0, i1) + 0.5) * MASK_PX,
                             y0 + (np.arange(j0, j1) + 0.5) * MASK_PX)
        m[j0:j1, i0:i1] ^= Path(q).contains_points(
            np.column_stack([gx.ravel(), gy.ravel()])).reshape(gy.shape)
    return m, x0, y0


GAP = 6.0            # mm of blank rim between neighbouring strings
CAP_FLOOR = 2.2      # mm smallest capital height that still prints: measured,
                     # the label's thinnest stroke is 0.45 mm here, just over
                     # one 0.4 mm nozzle.  Below this the glyphs break up.
CAP_CEIL = 3.0       # mm; a rim widened by zero-snapping must not grow its type


def build_rim_decor():
    """Lay the label strings out around each rim.

    Prefers three repeats per rim so a label is always within reach of the hand,
    shrinking the type to fit; falls back to two, then one, rather than let the
    capital height drop below CAP_FLOOR."""
    decor = []
    for i, (a, b) in enumerate(RIMS):
        outer = (i == len(RIMS) - 1)
        r_c = 0.5 * (a + b)
        flat = (b - a) - RIM_RAMP - (EDGE_CHAMFER if outer else RIM_RAMP)
        cap0 = float(np.clip(flat * 0.72, CAP_FLOOR, CAP_CEIL))

        strings = [(f"{DRIVE_HZ[i]:.0f} Hz → {DRIVE_HZ[i] / 2:.0f} Hz "
                    f"· {ZONE_PV[i]:.1f} mm p-v", cap0)]
        if outer:
            # set in caps: at this size lowercase counters close up below the
            # nozzle width, so 'a' and 'e' would print as solid blobs
            strings.append(("DRIVE → WAVE · SUBHARMONIC", cap0))

        w0 = sum(text_metrics(s, c)[0] for s, c in strings)
        circ = 2.0 * math.pi * r_c

        for n in (3, 2, 1):
            s = min(1.0, (circ / n - GAP * len(strings)) / w0)
            if n == 1 or min(c * s for _, c in strings) >= CAP_FLOOR:
                break
        # never let a string (descenders included) overflow the rim's flat top
        s = min(s, *[min(1.0, (flat * 0.42) / text_metrics(t, c)[1])
                     for t, c in strings])

        sized = []
        for t, c in strings:
            cap = max(c * s, CAP_FLOOR)
            w, hh = text_metrics(t, cap)
            sized.append((t, cap, w, hh))
        group_w = sum(x[2] for x in sized) + GAP * len(sized)

        placed = []
        slot = 2.0 * math.pi / n
        for g in range(n):
            u = 0.5 * (circ / n - group_w) + GAP * 0.5
            for t, cap, w, hh in sized:
                placed.append((t, cap, hh, g * slot, u, w))
                u += w + GAP
        decor.append({"r_c": r_c, "flat_lo": a + RIM_RAMP,
                      "flat_hi": b - (EDGE_CHAMFER if outer else RIM_RAMP),
                      "items": placed, "n": n,
                      "caps": [x[1] for x in sized]})
    return decor


DECOR = build_rim_decor()
_PATH_CACHE = {}


def decor_height(r_arr, theta):
    """Emboss height for matching arrays of radius and angle."""
    r_arr = np.asarray(r_arr, float)
    theta = np.asarray(theta, float)
    add = np.zeros(theta.shape)
    for rim in DECOR:
        idx_rim = np.flatnonzero((r_arr >= rim["flat_lo"])
                                 & (r_arr <= rim["flat_hi"]))
        if idx_rim.size == 0:
            continue
        r_c = rim["r_c"]
        v = r_arr[idx_rim] - r_c
        th = theta[idx_rim]
        for payload, cap, half_h, th0, u0, w in rim["items"]:
            # Reading direction runs clockwise: with letter tops pointing
            # outward (+v), the glyph's own +x axis is the -theta tangent.
            du = np.mod(th0 - th, 2.0 * math.pi) * r_c - u0
            sel = (np.abs(v) <= half_h + 0.2) & (du >= -0.2) & (du <= w + 0.2)
            if not sel.any():
                continue
            key = (payload, cap)
            if key not in _PATH_CACHE:
                _PATH_CACHE[key] = text_mask(payload, cap)
            mask, x0, y0 = _PATH_CACHE[key]
            ny, nx = mask.shape
            ix = ((du[sel] - x0) / MASK_PX).astype(np.int64)
            iy = ((v[sel] - y0) / MASK_PX).astype(np.int64)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            inside = np.zeros(ok.shape, bool)
            inside[ok] = mask[iy[ok], ix[ok]]
            hit = idx_rim[np.flatnonzero(sel)[inside]]
            add[hit] = np.maximum(add[hit], TEXT_EMBOSS)
    return add


# --------------------------------------------------------------------------
# mesh
# --------------------------------------------------------------------------

def ring_counts():
    n_rings = int(round(R_DISC / GRID))
    radii = np.linspace(0.0, R_DISC, n_rings + 1)
    counts = [1] + [max(6, int(round(2.0 * math.pi * r / GRID)))
                    for r in radii[1:]]
    return radii, counts


def stitch(n0, base0, n1, base1):
    """Triangulate the strip between two concentric rings of unequal counts by
    merging their vertices in angle order.  Ring 0 is the inner ring."""
    if n0 == 1:
        j = np.arange(n1)
        return np.column_stack([np.full(n1, base0), base1 + j,
                                base1 + (j + 1) % n1])
    e0 = (np.arange(1, n0 + 1)) / n0
    e1 = (np.arange(1, n1 + 1)) / n1
    keys = np.concatenate([e0, e1])
    is0 = np.concatenate([np.ones(n0, bool), np.zeros(n1, bool)])
    order = np.argsort(keys, kind="stable")
    t = is0[order]
    i_at = np.concatenate([[0], np.cumsum(t)[:-1]])
    j_at = np.concatenate([[0], np.cumsum(~t)[:-1]])

    tri = np.empty((n0 + n1, 3), np.int64)
    i0 = i_at % n0
    j0_ = j_at % n1
    tri[t, 0] = base0 + i0[t]
    tri[t, 1] = base1 + j0_[t]
    tri[t, 2] = base0 + (i_at[t] + 1) % n0
    tri[~t, 0] = base0 + i0[~t]
    tri[~t, 1] = base1 + j0_[~t]
    tri[~t, 2] = base1 + (j_at[~t] + 1) % n1
    return tri


def build_mesh():
    radii, counts = ring_counts()
    verts = []
    bases = []
    cursor = 0
    for r, n in zip(radii, counts):
        bases.append(cursor)
        cursor += n
        if n == 1:
            theta = np.zeros(1)
        else:
            theta = np.arange(n) * (2.0 * math.pi / n)
        z = (Z_MID + surface_profile(r)[0]
             + decor_height(np.full(n, r), theta))
        verts.append(np.column_stack([r * np.cos(theta), r * np.sin(theta), z]))
    top = np.vstack(verts)

    faces = [stitch(counts[j], bases[j], counts[j + 1], bases[j + 1])
             for j in range(len(counts) - 1)]

    n_out = counts[-1]
    out_top = bases[-1] + np.arange(n_out)
    off = len(top)
    theta = np.arange(n_out) * (2.0 * math.pi / n_out)
    bot_ring = np.column_stack([R_DISC * np.cos(theta), R_DISC * np.sin(theta),
                                np.zeros(n_out)])
    bot_centre = np.array([[0.0, 0.0, 0.0]])
    verts_all = np.vstack([top, bot_ring, bot_centre])

    j = np.arange(n_out)
    jn = (j + 1) % n_out
    bot = off + j
    botn = off + jn
    wall = np.vstack([
        np.column_stack([out_top, bot, botn]),
        np.column_stack([out_top, botn, bases[-1] + jn]),
    ])
    c = off + n_out
    bottom = np.column_stack([np.full(n_out, c), botn, bot])

    faces = np.vstack(faces + [wall, bottom])
    return verts_all, faces


def write_binary_stl(path, verts, faces, header=b"faraday disc"):
    tri = verts[faces]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)
    rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    assert rec.itemsize == 50, rec.itemsize
    data = np.zeros(len(faces), rec)
    data["n"] = nrm
    data["v"] = tri
    with open(path, "wb") as fh:
        fh.write(header.ljust(80, b" ")[:80])
        fh.write(struct.pack("<I", len(faces)))
        fh.write(data.tobytes())


def render_previews(px=1500):
    """Top-down height raster (definitive check that the rim decoration landed)
    plus a shaded 3D view."""
    import matplotlib.pyplot as plt

    ax1 = np.linspace(-R_DISC, R_DISC, px)
    gx, gy = np.meshgrid(ax1, ax1)
    r = np.hypot(gx, gy).ravel()
    th = np.arctan2(gy, gx).ravel()
    inside = r <= R_DISC
    z = np.full(r.shape, np.nan)
    z[inside] = (Z_MID + surface_profile(r[inside])
                 + decor_height(r[inside], th[inside]))
    z = z.reshape(px, px)

    # relief shading: light from the upper left, plus the raw height as colour
    gy_, gx_ = np.gradient(np.nan_to_num(z), 2.0 * R_DISC / px)
    shade = (1.0 - 0.6 * gx_ + 0.6 * gy_) / np.sqrt(1 + gx_ ** 2 + gy_ ** 2)
    shade[np.isnan(z)] = np.nan

    fig, axes = plt.subplots(1, 2, figsize=(17, 8.6))
    axes[0].imshow(z, extent=[-R_DISC, R_DISC] * 2, cmap="viridis",
                   origin="lower")
    axes[0].set_title("height field (mm above bed)")
    axes[1].imshow(shade, extent=[-R_DISC, R_DISC] * 2, cmap="gray",
                   origin="lower")
    axes[1].set_title("relief shading")
    for a in axes:
        a.set_xlabel("mm")
        a.set_aspect("equal")
    fig.tight_layout()
    fig.savefig("faraday_disc_preview.png", dpi=110)
    plt.close(fig)

    # rim close-ups, unwrapped, so the lettering and braille can be read
    fig, axes = plt.subplots(len(RIMS), 1, figsize=(15, 7))
    for i, rim in enumerate(DECOR):
        a, b = RIMS[i]
        u = np.linspace(0, 2.0 * math.pi * rim["r_c"] / rim["n"], 2400)
        vv = np.linspace(a, b, 170)
        uu, rr = np.meshgrid(u, vv)
        tt = -uu / rim["r_c"]
        zz = (surface_profile(rr.ravel())
              + decor_height(rr.ravel(), tt.ravel())).reshape(rr.shape)
        axes[i].imshow(zz, origin="lower", cmap="gray", aspect="equal",
                       extent=[0, u[-1], a, b])
        axes[i].set_title(f"rim {i + 1} unwrapped "
                          f"({rim['n']} label repeats around the ring)",
                          fontsize=9)
        axes[i].set_ylabel("r (mm)", fontsize=8)
    axes[-1].set_xlabel("arc length (mm)", fontsize=8)
    fig.tight_layout()
    fig.savefig("faraday_disc_rims.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 4))
    rp = np.linspace(0, R_DISC, 6000)
    ax.plot(rp, Z_MID + surface_profile(rp), lw=1.0, color="#185FA5")
    ax.axhline(Z_MID, color="0.7", lw=0.7, ls="--")
    for i, (z0, z1) in enumerate(ZONES):
        ax.axvspan(z0, z1, color="#378ADD", alpha=0.10)
        ax.text(0.5 * (z0 + z1), 15.2, f"{DRIVE_HZ[i]:.0f} Hz",
                ha="center", fontsize=9)
    ax.set_xlabel("radius (mm)")
    ax.set_ylabel("height above bed (mm)")
    ax.set_ylim(0, 16.5)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig("faraday_disc_section.png", dpi=130)
    plt.close(fig)


def report():
    names = "A B C".split()
    print("Faraday disc -- 200 mm, three concentric zones\n")
    print(f"{'zone':<5}{'drive':>8}{'wave':>8}{'lambda':>10}{'p-v':>8}"
          f"{'r range':>18}{'rings':>8}")
    for i, nm in enumerate(names):
        z0, z1 = ZONES[i]
        print(f"{nm:<5}{DRIVE_HZ[i]:>6.0f}Hz{DRIVE_HZ[i]/2:>6.0f}Hz"
              f"{LAMBDA[i]:>9.2f}mm{ZONE_PV[i]:>6.2f}mm"
              f"{f'{z0:.2f} - {z1:.2f}':>18}{RINGS[i]:>9}")
    print(f"\nmax surface slope, all zones : {math.degrees(math.atan(MAX_SLOPE)):.1f} deg")
    for i, (a, b) in enumerate(RIMS):
        caps = ", ".join(f"{c:.1f}" for c in DECOR[i]["caps"])
        print(f"rim {i + 1:<26}: {a:.2f} - {b:.2f} mm  "
              f"({b - a:.2f} wide, {DECOR[i]['n']} repeats, cap {caps} mm)")
    print(f"midplane above bed           : {Z_MID:.2f} mm")
    print(f"base slab under deepest point: {BASE_SLAB:.2f} mm")


def main():
    report()
    print("\nrendering previews ...")
    render_previews()
    print("\nbuilding mesh ...")
    verts, faces = build_mesh()
    print(f"  vertices {len(verts):,}   triangles {len(faces):,}")
    write_binary_stl("faraday_disc.stl", verts, faces)

    # the glyph masks and the raw mesh arrays together are enough to push
    # trimesh's load-time hashing into a MemoryError; drop them first
    del verts, faces
    _PATH_CACHE.clear()
    gc.collect()

    import trimesh
    m = trimesh.load("faraday_disc.stl")
    print(f"\nwatertight : {m.is_watertight}")
    print(f"winding ok : {m.is_winding_consistent}")
    print(f"extents    : {np.round(m.extents, 2)} mm")
    print(f"volume     : {m.volume / 1000.0:.1f} cm^3 "
          f"({m.volume / 1000.0 * 1.27:.0f} g solid PETG)")


if __name__ == "__main__":
    main()
