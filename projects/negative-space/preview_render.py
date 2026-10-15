#!/usr/bin/env python3
"""Preview render of the Negative Space ribbon described in README.md.

Not the STL generator (no segment joints): the ribbon is sampled on a (path,
width) grid, every hole is cut by dropping grid cells, and the result is given
its radial thickness as a closed shell. The close-up panel swaps in a finer grid
around Jul 1. Temperature uses the README's monthly normals interpolated to
days, so the band width is approximate.

    python preview_render.py            # README numbers -> negative_space_preview.png
    python preview_render.py <variant>  # see VARIANTS   -> negative_space_preview_<variant>.png"""

import math
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import pyvista as pv

import hanger_ring as hr

pv.OFF_SCREEN = True

R = 500.0                       # sphere radius, mm
N_TURNS = 12.37                 # set per variant in main()
YEAR = 2027
DAYS = 365
SYNODIC = 29.530589
REF_NEW_MOON = datetime(2000, 1, 6, 18, 14, tzinfo=timezone.utc)
MONTHLY_F = [29, 33, 43, 54, 64, 73, 76, 75, 68, 56, 44, 34]
MONTH_LEN = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
T_MIN, T_MAX = 29.0, 76.0
# turns, band width at the Jan / Jul normal (mm), moon aperture / width, thickness (mm),
# day-of-month numerals (digit height, slot width, bridge, which days). "sparse" is one
# turn per two lunations and the rest below it one per three, so full moons stack on
# two / three meridians.
VARIANTS = {
    "readme": dict(turns=12.37, w=(40.0, 90.0), aperture=0.6),
    "light": dict(turns=12.37, w=(27.0, 60.0), aperture=0.72),
    "sparse": dict(turns=DAYS / (2 * SYNODIC), w=(20.0, 45.0), aperture=0.72),
    "minimal": dict(turns=DAYS / (3 * SYNODIC), w=(16.0, 36.0), aperture=0.72),
    "numbered": dict(turns=DAYS / (3 * SYNODIC), w=(32.0, 72.0), aperture=0.72, thick=6.0,
                     numbers=dict(h=10.0, stroke=1.0, bridge=0.75, days=range(1, 32))),
    "bold": dict(turns=DAYS / (3 * SYNODIC), w=(32.0, 72.0), aperture=0.72, thick=6.0,
                 numbers=dict(h=20.0, stroke=2.5, bridge=1.2,
                              days=(1, 5, 10, 15, 20, 25, 30)), hang=True),
}
TAPER_LAT = 60.0                # band narrows poleward of this latitude
NOTCH_W, NOTCH_D = 2.0, 6.0     # daily notch on the upper edge
MONTH_NOTCH_W = 6.0             # month boundary: a wider notch, never a cut across the band
NOTCH_MAX = 0.2                 # notch depth capped at this share of the band width
APERTURE_CELL = 0.8             # aperture capped at this fraction of the day cell
LIGAMENT = 4.0                  # min plastic between a cut and the band edge
END_CLEAR = 1.5                 # min plastic between a moon and a notch
MIN_FEATURE = 2.5               # thinnest crescent slot / gibbous sliver worth cutting
# Integrity: no single cut may span more than MAX_CUT of the band width, and every
# cross-section across the band must keep at least MIN_NET of it as plastic.
MAX_CUT, MIN_NET = 0.5, 0.4
# Stencil seven-segment numerals: every segment is its own slot, so no counter falls out.
# They sit in a lane along the lower edge, clear of the notches, and may run into the
# neighbouring days; moons move to a lane above them.
DIGIT_ASPECT, DIGIT_SPACE = 0.7, 0.15    # digit width, gap between digits, / height
NUM_GAP = 3.0                   # plastic between the numeral and moon lanes
MIN_MOON = 8.0                  # numerals only where the moon can still be this wide
SEG_TABLE = np.array([[c in segs for c in "abcdefg"] for segs in
                      ["abcdef", "bc", "abdeg", "abcdg", "bcfg", "acdfg", "acdefg",
                       "abc", "abcdefg", "abcdfg"]])
PETG_G_MM3, ALU_G_MM3, STEEL_G_MM3 = 1.27e-3, 2.70e-3, 7.90e-3
TUBE_R, TUBE_WALL, SPOKE_R = 10.0, 1.0, 1.0    # aluminium axis tube, stainless spokes
# Hanging: the tube stops just under the bottom pole and rises HANG_TUBE_TOP above the
# top pole, where hanger_ring.py's tongue sits in a slot, bolted at HANG_BOLT.
HANG_TUBE_TOP = 50.0
HANG_BOLT = HANG_TUBE_TOP - (hr.TONGUE_LEN - hr.BOLT_EDGE)
NT, NV = 24000, 60              # grid cells along / across the band, full views
FINE_STEP, CLOSE_SPAN = 0.25, 0.025   # close-up grid step (mm), half-span in t


def lat_lon(t):
    """Centerline in radians; Jan 1 at the top pole, the year descends."""
    return np.radians(90.0 - 180.0 * t), np.radians(360.0 * N_TURNS * t)


def daily_temp_f():
    starts = np.concatenate([[0], np.cumsum(MONTH_LEN)[:-1]])
    mids = starts + np.array(MONTH_LEN) / 2.0
    return np.interp(np.arange(DAYS) + 0.5, mids, MONTHLY_F, period=DAYS)


def moon_phase():
    """Phase angle per day at noon UTC: 0 new, pi full."""
    jan1 = datetime(YEAR, 1, 1, 12, tzinfo=timezone.utc)
    ages = np.array([((jan1 + timedelta(days=d)) - REF_NEW_MOON).total_seconds()
                     / 86400.0 for d in range(DAYS)])
    return 2.0 * math.pi * np.mod(ages / SYNODIC, 1.0)


def sphere_pt(lat, lon, r=None):
    r = R if r is None else r                  # R read at call time: other scales override it
    return np.stack([r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon),
                     r * np.sin(lat)], axis=-1)


def seven_segments(u, v, h, w, s, g):
    """Slot masks a..g for a digit box u in [0, w], v in [0, h]; slots s wide,
    bridges g wide."""
    mid_u = (u >= s + g) & (u <= w - s - g)
    left, right = u <= s, u >= w - s
    upper, lower = v >= h / 2 + g / 2, v <= h / 2 - g / 2
    return [mid_u & (v >= h - s), right & upper, right & lower, mid_u & (v <= s),
            left & lower, left & upper, mid_u & (np.abs(v - h / 2) <= s / 2)]


class Ribbon:
    """Band width and every cut, as functions of path parameter t and offset y
    (mm across the band, +y toward the January pole). Cuts read correctly from
    outside: +x runs east along the path, waxing moons are lit on the +x side."""

    def __init__(self, w_min, w_max, aperture, numbers=None, root=None, notches=True):
        """root = (r, length, width): within `length` mm of axis radius `r` the band
        flares back up to `width`, so its ends can be fused to a pole cap.
        notches=False drops the day and month notches, leaving moons evenly spaced."""
        self.w_min, self.w_max, self.numbers, self.root = w_min, w_max, numbers, root
        self.notches = notches
        self.temp, self.phase = daily_temp_f(), moon_phase()
        self.illum = 0.5 * (1.0 - np.cos(self.phase))

        # Centerline arc length, used for along-path positions of every cut.
        self.tv = np.linspace(0.0, 1.0, NT + 1)
        self.cl = sphere_pt(*lat_lon(self.tv))
        self.s_v = np.concatenate(
            [[0.0], np.cumsum(np.linalg.norm(np.diff(self.cl, axis=0), axis=1))])
        self.s_day = np.interp(np.arange(DAYS + 1) / DAYS, self.tv, self.s_v)

        starts = np.concatenate([[0], np.cumsum(MONTH_LEN)])
        self.month_starts = starts[1:-1]
        month = np.searchsorted(starts, np.arange(DAYS), side="right") - 1
        self.dom = np.arange(DAYS) - starts[month] + 1
        month_at = np.isin(np.arange(DAYS + 1), self.month_starts) & notches
        clear = np.where(month_at, MONTH_NOTCH_W / 2, NOTCH_W / 2) + END_CLEAR
        c0, c1 = clear[:-1], clear[1:]          # at each day's start / end

        cell = np.diff(self.s_day)
        w = self.width((np.arange(DAYS) + 0.5) / DAYS)
        lane_top, lane_bot = 0.5 * w - LIGAMENT, -0.5 * w + LIGAMENT
        r_cap = np.full(DAYS, np.inf)
        self.labels = np.zeros((0, 3))
        if numbers:
            # A cross-section through a digit's stacked vertical slots and a moon
            # must still keep MIN_NET of the width.
            digit_cut = numbers["h"] - numbers["bridge"]
            moon_room = (1.0 - MIN_NET) * w - digit_cut
            lane_ok = ((w >= 2 * LIGAMENT + numbers["h"] + NUM_GAP + MIN_MOON)
                       & (moon_room >= MIN_MOON))
            lane_bot = np.where(lane_ok, lane_bot + numbers["h"] + NUM_GAP, lane_bot)
            self.labels = self.place_labels(lane_ok)
            for left, nw, _ in self.labels:
                under = (self.s_day[:-1] < left + nw) & (self.s_day[1:] > left)
                r_cap[under] = 0.5 * moon_room[under]
        self.y_moon = 0.5 * (lane_top + lane_bot)
        self.r_moon = np.clip(np.minimum.reduce([
            0.5 * aperture * w, 0.5 * MAX_CUT * w, 0.5 * (lane_top - lane_bot),
            0.5 * APERTURE_CELL * cell, 0.5 * cell - np.maximum(c0, c1), r_cap]), 0.0, None)
        # Skip crescents too thin to cut; open gibbous moons whose sliver would snap.
        d_moon = 2.0 * self.r_moon
        self.moon_none = d_moon * self.illum < MIN_FEATURE
        self.moon_full = ~self.moon_none & (d_moon * (1.0 - self.illum) < MIN_FEATURE)

    def width(self, t):
        lat, _ = lat_lon(t)
        day = np.clip((t * DAYS).astype(int), 0, DAYS - 1)
        w_t = self.w_min + (self.w_max - self.w_min) * (self.temp[day] - T_MIN) / (T_MAX - T_MIN)
        taper = np.clip(np.cos(lat) / math.cos(math.radians(TAPER_LAT)), 0.0, 1.0)
        if self.root is None:
            return w_t * taper
        r_cap, length, w_root = self.root
        flare = np.clip(1.0 - (R * np.cos(lat) - r_cap) / length, 0.0, 1.0)
        return np.maximum(w_t * taper, w_root * flare)

    def holes(self, t, y):
        lat_c, _ = lat_lon(t)
        half_w = 0.5 * self.width(t)
        day = np.clip((t * DAYS).astype(int), 0, DAYS - 1)
        stretch = np.cos(lat_c + y / R) / np.maximum(np.cos(lat_c), 1e-6)
        x = (np.interp(t, self.tv, self.s_v) - self.s_day[day]) * stretch
        x_end = (self.s_day[day + 1] - self.s_day[day]) * stretch
        edge_dist = np.minimum(x, x_end - x)
        boundary = np.where(x < x_end - x, day, day + 1)

        notch_w = np.where(np.isin(boundary, self.month_starts), MONTH_NOTCH_W, NOTCH_W)
        notch = ((edge_dist < notch_w / 2)
                 & (y > half_w - np.minimum(NOTCH_D, 2 * NOTCH_MAX * half_w)))

        # Moon aperture: the illuminated shape, waxing lit on the +x side.
        r, ph = self.r_moon[day], self.phase[day]
        xx, yy = x - 0.5 * x_end, y - self.y_moon[day]
        in_disc = xx ** 2 + yy ** 2 < r ** 2
        term = np.cos(ph) * np.sqrt(np.clip(r ** 2 - yy ** 2, 0.0, None))
        lit = np.where(ph < math.pi, xx >= term, xx <= -term)
        cut = in_disc & (lit | self.moon_full[day]) & ~self.moon_none[day]
        if self.notches:
            cut |= notch
        if len(self.labels):
            cut |= self.numeral(t, y, half_w, stretch)
        return cut

    def place_labels(self, lane_ok):
        """Numeral boxes along the centerline as rows of (left s, width, day of month):
        centred on their day, dropped where the band is too narrow for both lanes
        or a box would crowd its neighbour."""
        num = self.numbers
        dig_w, sp = DIGIT_ASPECT * num["h"], DIGIT_SPACE * num["h"]
        boxes = []
        for d in range(DAYS):
            dom = int(self.dom[d])
            if dom not in num["days"] or not lane_ok[max(d - 1, 0):d + 2].all():
                continue
            nw = 2 * dig_w + sp if dom >= 10 else dig_w
            left = 0.5 * (self.s_day[d] + self.s_day[d + 1] - nw)
            if boxes and left < boxes[-1][0] + boxes[-1][1] + 2 * END_CLEAR:
                continue
            boxes.append((left, nw, dom))
        return np.array(boxes, float).reshape(-1, 3)

    def numeral(self, t, y, half_w, stretch):
        """Slots of the nearest numeral box, in the lane along the lower edge."""
        num, lab = self.numbers, self.labels
        dig_w, sp = DIGIT_ASPECT * num["h"], DIGIT_SPACE * num["h"]
        s = np.interp(t, self.tv, self.s_v)
        centre = lab[:, 0] + 0.5 * lab[:, 1]
        k = np.clip(np.searchsorted(centre, s), 0, len(lab) - 1)
        km = np.maximum(k - 1, 0)
        k = np.where(np.abs(s - centre[km]) < np.abs(s - centre[k]), km, k)
        dom = lab[k, 2].astype(int)
        u = (s - lab[k, 0]) * stretch
        second = (dom >= 10) & (u >= dig_w + 0.5 * sp)
        u = u - np.where(second, dig_w + sp, 0.0)
        v = y + half_w - LIGAMENT
        on = SEG_TABLE[np.where((dom >= 10) & ~second, dom // 10, dom % 10)]
        hit = np.zeros(np.broadcast(u, v).shape, bool)
        segs = seven_segments(u, v, num["h"], dig_w, num["stroke"], num["bridge"])
        for i, seg in enumerate(segs):
            hit |= seg & on[..., i]
        return hit & (u >= 0) & (u <= dig_w) & (v >= 0) & (v <= num["h"])

    def integrity(self, step=0.3, chunk=4000):
        """Scan cross-sections across the band about every `step` mm along the whole
        ribbon. Returns (smallest share of the width left as plastic, t there) and
        (longest single cut as a share of the width, in mm, t there), ignoring the
        polar tips where the band is under 10 mm wide."""
        n_t, n_v = int(self.s_v[-1] / step), int(self.w_max / step)
        t_all = (np.arange(n_t) + 0.5) / n_t
        v = (np.arange(n_v) + 0.5) / n_v - 0.5
        net, run = np.ones(n_t), np.zeros(n_t)
        for i in range(0, n_t, chunk):
            t = t_all[i:i + chunk]
            cut = self.holes(t[:, None], v[None, :] * self.width(t)[:, None])
            net[i:i + chunk] = 1.0 - cut.mean(axis=1)
            cur, best = np.zeros(len(t)), np.zeros(len(t))
            for j in range(n_v):
                cur = (cur + 1.0) * cut[:, j]
                best = np.maximum(best, cur)
            run[i:i + chunk] = best / n_v
        w = self.width(t_all)
        net[w < 10.0], run[w < 10.0] = 1.0, 0.0
        i_net, i_run = np.argmin(net), np.argmax(run)
        return (net[i_net], t_all[i_net]), (run[i_run], run[i_run] * w[i_run], t_all[i_run])

    def solid(self, t0, t1, nt, nv, thick):
        """Ribbon between t0 and t1 as a closed shell: inner face on the sphere,
        outer face `thick` mm out, walls on every edge and hole. Returns the mesh
        and the band and plastic areas of the inner face."""
        tv = np.linspace(t0, t1, nt + 1)
        vv = np.linspace(-0.5, 0.5, nv + 1)
        lat_v, lon_v = lat_lon(tv)
        lat_g = lat_v[:, None] + vv[None, :] * self.width(tv)[:, None] / R
        pts = sphere_pt(lat_g, np.broadcast_to(lon_v[:, None], lat_g.shape)).reshape(-1, 3)

        tc = 0.5 * (tv[:-1] + tv[1:])
        vc = 0.5 * (vv[:-1] + vv[1:])
        keep = ~self.holes(tc[:, None], vc[None, :] * self.width(tc)[:, None])
        g = pts.reshape(nt + 1, nv + 1, 3)
        cell_area = 0.5 * np.linalg.norm(np.cross(g[1:, 1:] - g[:-1, :-1],
                                                  g[:-1, 1:] - g[1:, :-1]), axis=2)

        ii, jj = np.nonzero(keep)
        a = ii * (nv + 1) + jj
        quads = np.stack([a, a + nv + 1, a + nv + 2, a + 1], axis=1)
        n_pts = len(pts)
        e0, e1 = quads.ravel(), np.roll(quads, -1, axis=1).ravel()
        edge = e0.astype(np.int64) * n_pts + e1
        bnd = ~np.isin(e1.astype(np.int64) * n_pts + e0, edge)   # no twin: boundary
        walls = np.stack([e1[bnd], e0[bnd], e0[bnd] + n_pts, e1[bnd] + n_pts], axis=1)
        faces = np.vstack([quads[:, ::-1], quads + n_pts, walls])
        cells = np.hstack([np.full((len(faces), 1), 4), faces]).ravel()
        mesh = pv.PolyData(np.vstack([pts, pts * (R + thick) / R]), cells)
        centre = 0.25 * (g[:-1, :-1] + g[1:, :-1] + g[:-1, 1:] + g[1:, 1:]) * (R + thick / 2) / R
        moment = (cell_area[keep][:, None] * centre[keep]).sum(axis=0)
        return mesh, cell_area.sum(), cell_area[keep].sum(), moment

    def full_moon_az(self):
        full = [d for d in range(1, DAYS - 1) if self.illum[d] >= self.illum[d - 1]
                and self.illum[d] >= self.illum[d + 1]]
        fm_lon = np.array([lat_lon((d + 0.5) / DAYS)[1] for d in full])
        k = round(DAYS / SYNODIC / N_TURNS)    # lunations per turn = full-moon meridians
        return math.atan2(np.sin(k * fm_lon).mean(), np.cos(k * fm_lon).mean()) / k


def build_skeleton(cl, hang):
    """Axis tube, spokes and either a floor base or the hanging eye. Returns the
    meshes and (mass g, centroid) of each metal part for the balance check."""
    if hang:
        base_z, top_z = -R - 20.0, R + HANG_TUBE_TOP
    else:
        base_z, top_z = -R - 320.0, R + 20.0
    parts = [pv.Cylinder(center=(0, 0, 0.5 * (base_z + top_z)), direction=(0, 0, 1),
                         radius=TUBE_R, height=top_z - base_z, resolution=32)]
    masses = [(math.pi * (TUBE_R ** 2 - (TUBE_R - TUBE_WALL) ** 2) * (top_z - base_z)
               * ALU_G_MM3, np.array([0.0, 0.0, 0.5 * (base_z + top_z)]))]
    if not hang:
        parts.append(pv.Cylinder(center=(0, 0, base_z - 12), direction=(0, 0, 1), radius=220,
                                 height=24, resolution=96))
    n_spokes = int(2 * N_TURNS)
    for k in range(n_spokes):
        t = (k + 0.5) / (2 * N_TURNS)
        p = cl[int(round(t * NT))]
        horiz = math.hypot(p[0], p[1])
        if horiz < 4 * TUBE_R:
            continue
        start = p * (TUBE_R / horiz)          # true radius, starting at the tube wall
        parts.append(pv.Tube(pointa=start, pointb=p, radius=SPOKE_R, n_sides=12))
        masses.append((math.pi * SPOKE_R ** 2 * np.linalg.norm(p - start) * STEEL_G_MM3,
                       0.5 * (start + p)))
    return parts, masses


def hanger_meshes():
    """hanger_ring.stl stood upright on the tube, bolt at HANG_BOLT above the top
    pole, plus a cord rising out of the eye. Also returns the height of the hang point."""
    eye = pv.read("hanger_ring.stl")
    eye.rotate_x(90.0, point=(0.0, 0.0, 0.0), inplace=True)      # plate into the XZ plane
    eye.translate((0.0, 0.5 * hr.THICK, R + HANG_BOLT), inplace=True)
    _, yc = hr.outline()
    cord_z = R + HANG_BOLT + yc + 0.5 * hr.RING_ID - 1.5
    cord = pv.Tube(pointa=(0, 0, cord_z), pointb=(0, 0, cord_z + 260.0), radius=1.5,
                   n_sides=12)
    return eye, cord, R + HANG_BOLT + yc + 0.5 * hr.RING_ID


def main():
    global N_TURNS
    variant = sys.argv[1] if len(sys.argv) > 1 else "readme"
    cfg = VARIANTS[variant]
    N_TURNS = cfg["turns"]
    (w_min, w_max), aperture = cfg["w"], cfg["aperture"]
    thick, numbers, hang = cfg.get("thick", 3.0), cfg.get("numbers"), cfg.get("hang", False)
    out = "negative_space_preview.png" if variant == "readme" else \
        f"negative_space_preview_{variant}.png"

    rib = Ribbon(w_min, w_max, aperture, numbers)
    jul1 = sum(MONTH_LEN[:6])
    ta, tb = jul1 / DAYS - CLOSE_SPAN, jul1 / DAYS + CLOSE_SPAN
    pieces = [rib.solid(t0, t1, max(2, round(NT * (t1 - t0))), NV, thick)
              for t0, t1 in [(0.0, ta), (ta, tb), (tb, 1.0)]]
    span_len = np.interp(tb, rib.tv, rib.s_v) - np.interp(ta, rib.tv, rib.s_v)
    fine = rib.solid(ta, tb, int(span_len / FINE_STEP), int(w_max / FINE_STEP), thick)[0]
    skeleton, metal = build_skeleton(rib.cl, hang)
    fm_az = rib.full_moon_az()

    band = sum(p[1] for p in pieces)
    plastic = sum(p[2] for p in pieces)
    sphere = 4.0 * math.pi * R ** 2
    print(f"{variant}: {N_TURNS:.2f} turns, band {w_min:.0f}-{w_max:.0f} mm x {thick:.0f} mm, "
          f"aperture {aperture} x width, numerals "
          f"{'%.0f mm' % numbers['h'] if numbers else 'off'}")
    print(f"  {sum(p[0].n_cells for p in pieces)} faces, full-moon meridian at "
          f"{math.degrees(fm_az):.1f} deg; numerals on {len(rib.labels)} of {DAYS} days")
    print(f"  sphere covered by band {100 * band / sphere:.1f} %, by plastic "
          f"{100 * plastic / sphere:.1f} %; holes remove {100 * (1 - plastic / band):.1f} % "
          f"of the band")
    print(f"  ribbon mass at {thick:.0f} mm PETG: {plastic * thick * PETG_G_MM3 / 1000:.2f} kg")
    (net, t_net), (run, run_mm, t_run) = rib.integrity()
    date = lambda t: (datetime(YEAR, 1, 1) + timedelta(days=int(t * DAYS))).strftime("%b %d")
    ok = net >= MIN_NET - 0.01 and run <= MAX_CUT + 0.01
    print(f"  integrity {'OK' if ok else 'FAIL'}: thinnest cross-section keeps {100 * net:.0f} % "
          f"of the band ({date(t_net)}); longest single cut spans {100 * run:.0f} % "
          f"({run_mm:.1f} mm, {date(t_run)}); limits {100 * MIN_NET:.0f} % / "
          f"{100 * MAX_CUT:.0f} %")

    printed = []                                   # PETG parts besides the ribbon
    if hang:
        eye, cord, hang_z = hanger_meshes()
        printed.append(eye)
        # Balance: a hung piece turns until its centre of mass is under the eye.
        m_rib = plastic * thick * PETG_G_MM3
        parts = [(m_rib, sum(p[3] for p in pieces) / plastic)] + metal
        m_tot = sum(m for m, _ in parts)
        com = sum(m * c for m, c in parts) / m_tot
        off = math.hypot(com[0], com[1])
        tilt = math.degrees(math.atan2(off, hang_z - com[2]))
        print(f"  hanging: {m_tot / 1000:.2f} kg (ribbon {m_rib / 1000:.2f}, tube and spokes "
              f"{(m_tot - m_rib) / 1000:.2f}); centre of mass {off:.1f} mm off the axis, "
              f"{hang_z - com[2]:.0f} mm below the eye -> hangs {tilt:.1f} deg off vertical")
        skeleton.append(cord)

    near = rib.cl[int(round(jul1 / DAYS * NT))]
    if numbers:                                    # aim at the numeral lane instead
        t1 = np.array((jul1 + 0.5) / DAYS)
        drop = 0.5 * rib.width(t1) - LIGAMENT - 0.5 * numbers["h"]
        lat1, lon1 = lat_lon(t1)
        near = sphere_pt(lat1 - drop / R, lon1)
    off_az = math.atan2(near[1], near[0]) + (0.22 if numbers else 0.35)  # tube off-centre
    close_dir = (math.cos(off_az), math.sin(off_az), 0.12 if numbers else 0.3)
    whole = [p[0] for p in pieces]
    close = [pieces[0][0], fine, pieces[2][0]]

    # (title, view direction, up, focal point, distance or None to fit, ribbon meshes)
    views = [
        ("Three-quarter view", (math.cos(fm_az + 0.6), math.sin(fm_az + 0.6), 0.42),
         (0, 0, 1), (0.0, 0.0, -90.0), None, whole),
        ("Full-moon meridian, face on", (math.cos(fm_az), math.sin(fm_az), 0.0),
         (0, 0, 1), (0.0, 0.0, -90.0), None, whole),
        ("From the January pole", (0.0, 0.0, 1.0),
         (math.cos(fm_az), math.sin(fm_az), 0), (0.0, 0.0, -90.0), None, whole),
        ("Close-up: Jul 1 month notch, day notches, phases"
         + (", day numbers" if numbers else ""),
         close_dir, (0, 0, 1), tuple(near), 300.0 + 22.0 * numbers["h"] if numbers
         else 480.0, close),
    ]
    pl = pv.Plotter(off_screen=True, shape=(1, len(views)),
                    window_size=(900 * len(views), 1100), border=False)
    for col, (title, direction, up, focal, dist, meshes) in enumerate(views):
        pl.subplot(0, col)
        pl.set_background("#1c1e23")
        for m in meshes + printed:
            pl.add_mesh(m, color="#e8dfcf", specular=0.15, ambient=0.18)
        for part in skeleton:
            pl.add_mesh(part, color="#9aa3ad", smooth_shading=True, specular=0.5)
        d = np.asarray(direction, float)
        d /= np.linalg.norm(d)
        pl.camera.focal_point = focal
        pl.camera.position = tuple(d * (dist or 4000.0) + np.asarray(focal))
        pl.camera.up = up
        pl.camera.view_angle = 30.0
        if dist is None:
            pl.reset_camera()
            pl.camera.zoom(1.15)
        else:
            pl.reset_camera_clipping_range()
        pl.add_text(title, position="upper_left", font_size=13, color="#d6d0c4")
    pl.add_text(f"Negative Space: 1 m ribbon, {N_TURNS:.2f} turns, band {w_min:.0f}-"
                f"{w_max:.0f} x {thick:.0f} mm, moon {YEAR}",
                position="lower_right", font_size=10, color="#8a8f98")
    pl.screenshot(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
