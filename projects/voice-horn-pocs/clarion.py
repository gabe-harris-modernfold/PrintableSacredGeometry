#!/usr/bin/env python3
"""POC-A CLARION -- printable STLs for both Clarion lanes.

Cycle 3 split this POC into two honest devices, and Part XVI calls for printing
both (docs/voice-horn-pocs.md):

  projector  the cycle-5 KNEE genome, straight off pareto_c5.csv:
             cup 40 -> tube d39.2 x 16.5 -> step 1.22 -> 501.3 mm bore,
             r(u) = r_t + (r_m-r_t) u^n + A(r_m-r_t) exp(-((u-x0)/0.12)^2) sin(pi u)
             r_t 21.695, r_m 143.575, n 1.7373, A -0.1411, x0 0.80
  squillo    the Rev-1 hand geometry, exempt from the ripple constraint:
             cup 40 -> tube d22 x 30 -> 1:5.8 step -> 280 mm exponential bore
             r(u) = 26.5 (125/26.5)^u, plus the 63-well QRS Motto Ring.

Walls are the F1 sandwich as GEOMETRY ONLY: 5 mm of solid STL wall, to be sliced
as two perimeters over 40% Grid infill (never Cubic -- F3).

Run from this directory:  python clarion.py [--no-skin] [--variant p|s|both]
"""
import argparse
import math

import numpy as np
import trimesh

import horn_lib as H

OUT = "stl/"

# ---- geometry constants ------------------------------------------------------
R_CUP, L_CUP = 22.5, 40.0            # fixed lip cup in every Forge cycle
STEP_RAMP = 1.0                      # chamfer standing in for the abrupt step

# cycle-5 KNEE genes (pareto_c5.csv, knee of G_speech vs Gamma)
P = dict(re=19.6005, le=16.4731, rt=21.6950, rm=143.5751, L=501.2632,
         n=1.7373, A=-0.1411, x0=0.80)
# Rev-1 hand geometry (the squillo lane; c5 script's `Xhand` + r_exp)
S = dict(re=11.0, le=30.0, rt=26.5, rm=125.0, L=280.0)

# joint hardware
SP_LEN, SP_WALL, TH_PITCH, TH_DEPTH, CLR, SOCK_WALL = 20.0, 3.0, 8.0, 1.5, 0.35, 3.5
BOLT_D, BOLT_N = 4.4, 12

# QRS Motto Ring (RULE S4 / Skin section): N = 7, 3.5 kHz, 9 periods
QRS_N, QRS_PERIODS, QRS_UNIT = 7, 9, 7.0
QRS_DIV = 1.6                        # divider wall thickness, 4 extrusions
QRS_DEPTHS = [(k * k % QRS_N) * QRS_UNIT for k in range(QRS_N)]   # 0,7,28,14,14,28,7
QRS_LEN = 40.0
FLANGE_R, FLANGE_T, BOLT_R = 143.0, 6.0, 137.5


def r_flare_p(x):
    u = np.asarray(x, float) / P["L"]
    d = P["rm"] - P["rt"]
    return (P["rt"] + d * u ** P["n"]
            + P["A"] * d * np.exp(-((u - P["x0"]) / 0.12) ** 2) * np.sin(math.pi * u))


def r_flare_s(x):
    return S["rt"] * (S["rm"] / S["rt"]) ** (np.asarray(x, float) / S["L"])


class Bore:
    """Piecewise bore radius of the whole signal path, lips (z=0) to mouth."""

    def __init__(self, g, flare):
        self.g, self.flare = g, flare
        self.z_tube = L_CUP
        self.z_ramp = L_CUP + g["le"] - STEP_RAMP
        self.z_throat = L_CUP + g["le"]
        self.z_end = self.z_throat + g["L"]

    def r(self, z):
        z = np.atleast_1d(np.asarray(z, float))
        g = self.g
        out = np.empty_like(z)
        cup = z < self.z_tube
        out[cup] = R_CUP + (g["re"] - R_CUP) * z[cup] / L_CUP
        tub = (z >= self.z_tube) & (z < self.z_ramp)
        out[tub] = g["re"]
        rmp = (z >= self.z_ramp) & (z < self.z_throat)
        out[rmp] = g["re"] + (g["rt"] - g["re"]) * (z[rmp] - self.z_ramp) / STEP_RAMP
        fl = z >= self.z_throat
        out[fl] = self.flare(np.minimum(z[fl], self.z_end) - self.z_throat)
        return out

    def stations(self, z0, z1, dz=1.5, extra=()):
        """Feature-aware station list: breaks pinned wherever the bore kinks."""
        brk = [self.z_tube, self.z_ramp, self.z_throat] + list(extra)
        cuts = sorted({z0, z1, *[b for b in brk if z0 < b < z1]})
        out = []
        for a, b in zip(cuts[:-1], cuts[1:]):
            n = max(2, int(math.ceil((b - a) / dz)) + 1)
            out.append(np.linspace(a, b, n, endpoint=False))
        return np.unique(np.concatenate(out + [np.array([z1])]))

    def acoustics(self, name):
        m = math.log((self.g["rm"] / self.g["rt"]) ** 2) / (self.g["L"] / 1000.0)
        fc = m * 343.0 / (4 * math.pi)
        fm = 343000.0 / (2 * math.pi * self.g["rm"])
        ft = 343000.0 / (4 * (self.g["le"] + 0.85 * self.g["re"]))
        print(f"{name}: total {self.z_end:.1f} mm | throat A {math.pi*self.g['rt']**2:.0f} "
              f"mouth A {math.pi*self.g['rm']**2:.0f} mm2 | step "
              f"{(self.g['rt']/self.g['re'])**2:.2f}:1")
        print(f"   f_c {fc:.0f} Hz [Law 14.7]   f_match {fm:.0f} Hz [Law 14.8]   "
              f"tube 1/4-wave {ft:.0f} Hz [15.11]")
        return fc, fm, ft


# ---- field helpers ----------------------------------------------------------

def lip_bulge(z, z0):
    """Rounded mouthpiece rim over the first 4 mm (sand before use)."""
    t = np.clip((z - z0) / 4.0, 0.0, 1.0)
    return -4.0 * (1.0 - np.sin(0.5 * math.pi * t))


def body(bore, z0, z1, nt=192, dz=1.5, theta=None, lip=False,
         male_at=None, female_at=None, flange_at=None, extra=(),
         skin_field=None):
    """Shell fields for one horn segment; returns (mesh, z, r_bore).

    `skin_field(z, th) -> (nz, nt)` is added to the OUTER radius last -- the
    exterior texture (RULE S3). Callers window it clear of spigots, sockets
    and flanges; the joints stay glass."""
    th = H.grid_theta(nt) if theta is None else np.asarray(theta, float)
    ex = list(extra)
    for j in (male_at, female_at):
        if j is not None:
            ex += [j, j + SP_LEN, j + SP_LEN + 0.5, j + SP_LEN + 10.0]
    if flange_at is not None:
        ex += [flange_at, flange_at + FLANGE_T]
    z = bore.stations(z0, z1, dz=dz, extra=[e for e in ex if z0 < e < z1])
    rb = bore.r(z)
    ri = np.repeat(rb[:, None], len(th), 1)
    ro = ri + H.WALL
    if lip:
        ro += lip_bulge(z, z0)[:, None]

    if male_at is not None:                       # spigot: thinner wall + ridge
        j = male_at
        w = np.where(z >= j, SP_WALL,
                     H.WALL + (SP_WALL - H.WALL) * np.clip((z - (j - 8.0)) / 8.0, 0, 1))
        ro = ri + w[:, None] + H.thread(z, th, j + 1.5, j + SP_LEN - 1.5,
                                        TH_PITCH, TH_DEPTH)
    if female_at is not None:                     # socket: male field + clearance
        j = female_at
        inside = z <= j + SP_LEN
        base = ri + SP_WALL + CLR
        rid = H.thread(z, th, j + 1.5, j + SP_LEN - 1.5, TH_PITCH, TH_DEPTH)
        ri = np.where(inside[:, None], base + rid, ri)
        collar = ri[:, :1] * 0 + (rb[:, None] + SP_WALL + TH_DEPTH + CLR + SOCK_WALL)
        blend = np.clip((z - (j + SP_LEN + 0.5)) / 10.0, 0.0, 1.0)[:, None]
        ro = np.where(inside[:, None], collar,
                      collar * (1 - blend) + (rb[:, None] + H.WALL) * blend)
    if flange_at is not None:                     # annular bolt pad
        f = (z >= flange_at) & (z <= flange_at + FLANGE_T)
        ro = np.where(f[:, None], np.maximum(ro, FLANGE_R), ro)
    if skin_field is not None:
        ro = ro + skin_field(z, th)
    return H.shell(z, ri, ro, theta=th), z, rb


class Grip:
    """The fixed right-hand print of a throat segment: a +2 mm bulge band with
    four canted finger pits and an opposed thumb pit carved back to the nominal
    surface -- the pits read deep, the 5 mm wall is never thinned (RULE S3's
    blind-hand test: the hand falls into the print, everywhere else is jagged).

    Clocking: braille rings start at theta 0 and run counter-clockwise, so the
    fingers land opposite them (around pi) and the thumb sits just clockwise of
    the braille -- braille under the left index while the right hand holds."""

    RISE, DEPTH, CANT = 2.0, 2.0, 0.12            # mm, mm, rad per finger

    def __init__(self, bore, z0, z1, fingers_z, thumb_z, ramp=8.0):
        self.z0, self.z1, self.ramp = z0, z1, ramp
        r_at = lambda zc: float(bore.r(zc)[0]) + H.WALL + self.RISE
        self.pits = [dict(z=zc, th=math.pi + 0.55 + self.CANT * k, sz=8.5,
                          sa=26.0, r=r_at(zc), depth=self.DEPTH)
                     for k, zc in enumerate(fingers_z)]
        self.pits.append(dict(z=thumb_z, th=-0.95, sz=11.0, sa=15.0,
                              r=r_at(thumb_z), depth=self.DEPTH))

    def bulge(self, z):
        return self.RISE * H.smooth_band(z, self.z0, self.z1, self.ramp)

    def field(self, z, th):
        return self.bulge(z)[:, None] - H.pit_field(z, th, self.pits)

    def sink_extra(self, zs, ths):
        return H.pit_depth_at(zs, ths, self.pits)


def qrs_theta():
    """Non-uniform theta: two tight samples on each divider edge, 6 per well."""
    nsec = QRS_N * QRS_PERIODS
    span = 2 * math.pi / nsec
    half = 0.5 * QRS_DIV / S["rm"]                # divider half-width in radians
    out = []
    for k in range(nsec):
        c = k * span
        out += [c - half - 1e-4, c - half, c + half, c + half + 1e-4]
        out += list(c + half + 1e-4 + (span - 2 * half - 2e-4)
                    * np.linspace(0.08, 0.92, 6))
    return np.asarray(out, float) % (2 * math.pi)


def qrs_depth_field(th):
    """Per-theta well depth (0 on the dividers) and the depth the OUTER wall must
    clear -- on a divider that is the deeper of its two neighbours, so the wall
    steps with the sequence instead of running at max depth all the way round.
    Depths are (n^2 mod 7) x 7 mm."""
    nsec = QRS_N * QRS_PERIODS
    span = 2 * math.pi / nsec
    half = 0.5 * QRS_DIV / S["rm"]
    k = np.floor((th + half) / span).astype(int) % nsec
    ctr = (k * span) % (2 * math.pi)
    dth = (th - ctr + math.pi) % (2 * math.pi) - math.pi
    own = np.array([QRS_DEPTHS[i % QRS_N] for i in k])
    prev = np.array([QRS_DEPTHS[(i - 1) % QRS_N] for i in k])
    on_div = np.abs(dth) <= half + 1e-9
    return (np.where(on_div, 0.0, own),
            np.where(on_div, np.maximum(own, prev), own))


# ---- parts ------------------------------------------------------------------

def _skin_profile(bore, zp, grip=None):
    """Outer-surface profile the ornament sits on: wall + any grip bulge."""
    rp = bore.r(zp) + H.WALL
    return rp + grip.bulge(zp) if grip is not None else rp


def thorn_field(bore, z, z_lo, z_hi, n, h_lo, h_hi, m_law, rake=1.0, tip=0.6,
                grip=None):
    """Phyllotaxis dome/stud field, height on the flare law (RULE S2/S3).

    With `grip`, sites ride the bulge -- and inside a finger pit they drop to
    the carved floor, staying their own height proud instead of floating over
    it or towering out of it."""
    zp = np.linspace(z[0], z[-1], 400)
    rp = _skin_profile(bore, zp, grip)
    ramp = lambda zz: (zz < bore.z_ramp - 2.0) | (zz > bore.z_throat + 1.0)
    P_, N_, frac, area = H.phyllotaxis(zp, rp, n,
                                       weight=lambda zz: (0.35 + 0.65 *
                                       np.clip((zz - z_lo) / max(z_hi - z_lo, 1), 0, 1))
                                       * ramp(zz),   # no sites on the step cliff:
                                       z_lo=z_lo, z_hi=z_hi)  # its normal is ~ -z
    x = (P_[:, 2] - z_lo) / max(z_hi - z_lo, 1e-9)
    h = h_lo + (h_hi - h_lo) * (np.exp(m_law * x) - 1) / (math.exp(m_law) - 1)
    ax = N_ + rake * np.array([0.0, 0.0, 1.0])
    rb = np.minimum(0.42 * math.sqrt(area), tip + h * math.tan(math.radians(30)))
    if grip is not None:
        d = grip.sink_extra(P_[:, 2], np.arctan2(P_[:, 1], P_[:, 0]))
        nu = N_ / np.linalg.norm(N_, axis=1, keepdims=True)
        P_ = P_ - nu * d[:, None]
    return H.frusta(P_, ax, h, np.maximum(rb, tip + 0.3), r_tip=tip, sink=1.5)


def claw_field(bore, z, z_lo, z_hi, n, h_lo, h_hi, m_law, tip=0.18):
    """The sharpened thorn crown: two-segment claws curling toward the mouth,
    4:1 slenderness, genuinely sharp tips (they will scratch a careless grab --
    that is RULE S3's 'painful to seize' made literal)."""
    zp = np.linspace(z[0], z[-1], 400)
    rp = bore.r(zp) + H.WALL
    P_, N_, frac, area = H.phyllotaxis(zp, rp, n,
                                       weight=lambda zz: 0.35 + 0.65 *
                                       np.clip((zz - z_lo) / max(z_hi - z_lo, 1), 0, 1),
                                       z_lo=z_lo, z_hi=z_hi)
    x = (P_[:, 2] - z_lo) / max(z_hi - z_lo, 1e-9)
    h = h_lo + (h_hi - h_lo) * (np.exp(m_law * x) - 1) / (math.exp(m_law) - 1)
    rb = np.clip(h / 7.0, 0.7, 0.42 * math.sqrt(area))
    return H.claws(P_, N_, h, rb, r_tip=tip, curl=np.array([0.0, 0.0, 1.0]),
                   curl_amt=0.55)                 # [roots, tips]: union both


def braille_rings(bore, lines, z0, pitch=12.0, grip=None):
    """Grip-zone braille: dots inside the <=1 mm pleasure-dome budget (S3/S4)."""
    zp = np.linspace(z0 - 5, z0 + pitch * len(lines) + 10, 200)
    rp = _skin_profile(bore, zp, grip)
    pts, nrm = [], []
    for i, line in enumerate(lines):
        d, _ = H.braille_dots(line)
        p, n = H.wrap_dots(d, z0 + i * pitch, zp, rp)
        pts.append(p); nrm.append(n)
    P_ = np.concatenate(pts); N_ = np.concatenate(nrm)
    if grip is not None:                          # a dot inside a finger pit
        d = grip.sink_extra(P_[:, 2], np.arctan2(P_[:, 1], P_[:, 0]))
        P_ = P_ - N_ * d[:, None]                 # rides the carved floor, so it
    return H.frusta(P_, N_, np.full(len(P_), 0.6), np.full(len(P_), 0.75),
                    r_tip=0.45, nseg=6, sink=0.35)  # stays 0.6 proud and legible


def pop_cage(r_bore, pitch=2.5, strand=0.4, gap=3.0, h=9.0, wall=1.6):
    """Two offset crossed grids across the bore, right after the area step.

    A plosive is a momentum jet and cannot turn corners; the voice is a pressure
    wave whose wavelength dwarfs the spacing. RULE W2 keeps it out of the
    constriction -- this sits one area-step downstream."""
    r_in = r_bore - 0.35                          # slip fit into the bore
    segs = []
    for lay, z in enumerate((2.0, 2.0 + gap)):
        off = 0.0 if lay == 0 else 0.5 * pitch
        for s in np.arange(-r_in + 1.0 + off, r_in - 1.0, pitch):
            y = math.sqrt(max(r_in ** 2 - s ** 2, 0.0)) - 0.2
            if y < 1.0:
                continue
            if lay == 0:
                segs.append([(s, -y, z), (s, y, z)])
            else:
                segs.append([(-y, s, z), (y, s, z)])
    grid = H.struts(np.asarray(segs, float), 0.5 * strand, nseg=5)
    ring = H.shell(np.linspace(0, h, 6),
                   np.full((6, 96), r_in - 0.3), np.full((6, 96), r_in + wall))
    return H.union(ring, grid)


def fleece_coupon(r_in=22.0, thick=6.0, h=30.0, cell_wall=1.8, cell_lumen=3.4,
                  phi_wall=0.85, phi_lumen=0.96, strand_min=0.4):
    """Graded open-lattice sleeve for the RULE W5 A/B test (pass 7).

    RULE W4 asks for porosity 0.85 against the wall grading to 0.96 at the lumen
    face, so the STRAND diameter is solved per layer from the strut length that
    layer actually carries -- picking one strand size and hoping lands at 0.96
    everywhere, which is transparent (W3) but does no flow conditioning. The
    registration band is on the WALL side only: a face sheet at the lumen face is
    exactly what W5 forbids (-11.6 dB in the FEM)."""
    nr = max(2, int(round(thick / 2.0)) + 1)
    radii = np.linspace(r_in + 0.6, r_in + thick, nr)
    t_layer = thick / nr
    parts, report = [], []
    for i, r in enumerate(radii):
        f = i / max(nr - 1, 1)                    # 0 at the lumen, 1 at the wall
        cell = cell_lumen + (cell_wall - cell_lumen) * f
        phi = phi_lumen + (phi_wall - phi_lumen) * f
        nth = max(6, int(round(2 * math.pi * r / cell)))
        nz = max(2, int(round(h / cell)))
        th = np.linspace(0, 2 * math.pi, nth, endpoint=False) + 0.3 * i
        zs = np.linspace(1.0, h - 1.0, nz)
        length = nz * 2 * math.pi * r + nth * (zs[-1] - zs[0])
        layer_vol = 2 * math.pi * r * t_layer * h
        strand = max(strand_min,
                     2.0 * math.sqrt((1 - phi) * layer_vol / (math.pi * length)))
        sr = 0.5 * strand
        phi_got = 1 - math.pi * sr ** 2 * length / layer_vol
        report.append((cell, strand, phi_got))
        for z in zs:                              # hoops
            p = np.stack([r * np.cos(th), r * np.sin(th), np.full(nth, z)], -1)
            parts.append(H.tube_path(p, sr, nseg=5, closed=True))
        for t in th:                              # rails
            parts.append(H.tube_path([(r * math.cos(t), r * math.sin(t), zs[0]),
                                      (r * math.cos(t), r * math.sin(t), zs[-1])],
                                     sr, nseg=5))
        if i:                                     # radial ties to the layer below
            rp = radii[i - 1]
            for t in th[::max(1, 4 - i)]:
                for zm in zs[1::max(1, len(zs) // 4)]:
                    parts.append(H.tube_path(
                        [(rp * math.cos(t), rp * math.sin(t), zm),
                         (r * math.cos(t), r * math.sin(t), zm)], sr, nseg=5))
    band = H.shell(np.linspace(0, h, 4),
                   np.full((4, 96), r_in + thick),
                   np.full((4, 96), r_in + thick + 1.2))
    print("fleece coupon (lumen -> wall):  " + "  ".join(
        f"cell {c:.1f} strand {d:.2f} phi {p:.2f}" for c, d, p in report))
    return H.union(band, *parts)


def build_projector(skin=True):
    bore = Bore(P, r_flare_p)
    fc, fm, ft = bore.acoustics("CLARION projector (c5 KNEE)")
    z_split = 280.0
    rows = []

    grip = Grip(bore, 44.0, 150.0, fingers_z=(58.0, 81.0, 104.0, 127.0),
                thumb_z=72.0) if skin else None
    m, z, _ = body(bore, 0.0, z_split + SP_LEN, lip=True, male_at=z_split,
                   nt=288 if skin else 192, dz=0.8 if skin else 1.5,
                   skin_field=grip.field if skin else None)
    if skin:
        m = H.union(m,
                    thorn_field(bore, z, 44.0, 150.0, 900, 0.5, 1.0, 1.0, rake=0.0,
                                tip=0.45, grip=grip),
                    braille_rings(bore, ["clarion p",
                                         f"{fc:.0f} {fm:.0f} {ft:.0f}",
                                         "psg 2026"], 62.0, grip=grip))
    rows.append(H.report("clarion_p_throat", m, OUT + "clarion_p_throat.stl"))

    m, z, _ = body(bore, z_split, bore.z_end, female_at=z_split, dz=2.0)
    if skin:
        m = H.union(m, *claw_field(bore, z, 330.0, bore.z_end - 2.0, 1500,
                                   2.0, 12.0, 2.4))
    rows.append(H.report("clarion_p_bell", m, OUT + "clarion_p_bell.stl"))

    rows.append(H.report("clarion_p_popcage", pop_cage(P["rt"]),
                         OUT + "clarion_p_popcage.stl"))
    return rows


def build_squillo(skin=True):
    bore = Bore(S, r_flare_s)
    fc, fm, ft = bore.acoustics("CLARION squillo (Rev-1 hand)")
    z_split = 145.0
    z_ring = bore.z_end - QRS_LEN - FLANGE_T      # ring = flange + 40 mm of wells
    rows = []

    grip = Grip(bore, 44.0, 124.0, fingers_z=(52.0, 74.0, 96.0, 118.0),
                thumb_z=64.0) if skin else None
    m, z, _ = body(bore, 0.0, z_split + SP_LEN, lip=True, male_at=z_split,
                   nt=288 if skin else 192, dz=0.8 if skin else 1.5,
                   skin_field=grip.field if skin else None)
    if skin:
        m = H.union(m,
                    thorn_field(bore, z, 44.0, 120.0, 700, 0.5, 1.0, 1.0, rake=0.0,
                                tip=0.45, grip=grip),
                    braille_rings(bore, ["clarion s", "303 437 2900", "psg 2026"],
                                  74.0, grip=grip))   # above the 1:5.8 step cliff
    rows.append(H.report("clarion_s_throat", m, OUT + "clarion_s_throat.stl"))

    m, z, _ = body(bore, z_split, z_ring, female_at=z_split,
                   flange_at=z_ring - FLANGE_T)
    m = H.difference(m, bolt_ring(z_ring - FLANGE_T - 1, FLANGE_T + 2))
    if skin:
        m = H.union(m, *claw_field(bore, z, 190.0, z_ring - FLANGE_T - 6.0, 700,
                                   1.6, 10.5, 2.0))
    rows.append(H.report("clarion_s_bell", m, OUT + "clarion_s_bell.stl"))

    rows.append(H.report("clarion_s_qrs_ring", qrs_ring(bore, z_ring),
                         OUT + "clarion_s_qrs_ring.stl"))
    rows.append(H.report("clarion_s_popcage", pop_cage(S["rt"]),
                         OUT + "clarion_s_popcage.stl"))
    return rows


def bolt_ring(z0, h):
    cuts = []
    for k in range(BOLT_N):
        a = 2 * math.pi * k / BOLT_N + math.pi / BOLT_N
        c = (BOLT_R * math.cos(a), BOLT_R * math.sin(a))
        cuts.append(H.cyl((c[0], c[1], z0), (c[0], c[1], z0 + h), 0.5 * BOLT_D, 20))
    return trimesh.util.concatenate(cuts)


def qrs_ring(bore, z_flange):
    """The Motto Ring: 63 quadratic-residue wells, 9 periods of N = 7.

    Depths (n^2 mod 7) x 7 mm = 0/7/28/14/14/28/7 -> design frequency 3.5 kHz;
    below it the wells are invisible (RULE S1), above it they scatter the
    2.5-5 kHz harshness band off-axis, and the 63 well floors are the 63 slots
    Schiller's motto is engraved into.

    Printed wells-up: every well floor is then a flat surface laid over solid
    material -- no bridging -- and the 0.16 mm layers the lettering wants apply
    to this part alone while the bodies print at 0.3 mm (F2). The outer wall
    steps with the sequence, so the ring's own silhouette IS the crown; at
    dia 316 in a 320 bed it has no room for added ornament."""
    th = qrs_theta()
    dep, wall_dep = qrs_depth_field(th)
    z_wells = z_flange + FLANGE_T
    z = np.unique(np.concatenate([
        np.linspace(z_flange, z_wells, 5), np.linspace(z_wells, z_wells + 0.6, 3),
        np.linspace(z_wells + 0.6, bore.z_end, int(QRS_LEN / 1.5))]))
    rb = bore.r(z)
    gate = np.clip((z - z_wells) / 0.6, 0.0, 1.0)
    ri = rb[:, None] + dep[None, :] * gate[:, None]
    ro = rb[:, None] + wall_dep[None, :] + H.WALL
    ro = np.where((z <= z_wells)[:, None], np.maximum(ro, FLANGE_R), ro)
    m = H.shell(z, ri, ro, theta=th)
    return H.difference(m, bolt_ring(z_flange - 1.0, FLANGE_T + 2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skin", action="store_true",
                    help="add the ornament pass (claws, grip print, braille); "
                         "default is the plain smooth bodies")
    ap.add_argument("--variant", default="both", choices=["p", "s", "both"])
    a = ap.parse_args()
    rows = []
    if a.variant in ("p", "both"):
        rows += build_projector(skin=a.skin)
    if a.variant in ("s", "both"):
        rows += build_squillo(skin=a.skin)
    rows.append(H.report("fleece_coupon", fleece_coupon(),
                         OUT + "fleece_coupon.stl"))
    print()
    H.print_table(rows)
