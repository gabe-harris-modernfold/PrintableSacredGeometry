#!/usr/bin/env python3
"""The third Orphic sculpture: the empty coil shedding its skin. It is the serpent of
`orphic_egg.py --serpent`, unchanged, with its old skin peeled back from the head and bunched
over the rear half of the body. Ahead of the peel front the new scales show. At the front the
old skin is torn into flared tongues and pushed up into concertina folds. Behind it the skin
lies loose and crumpled, with tears that show the new scales, a few of them peeled up as flaps.
Round the tail on the ground it has slumped into crumpled sheets.

The old skin is built in the serpent's own tube coordinates: s along the centreline, phi round
it from the belly, rho out from it. It has an outer radius and, where it lifts off the body,
an inner one too, so the skin is solid and fused to the body where it lies down and a 1.3 mm
shell where it lifts. Its surface is embossed with the net of scale outlines that a shed skin
keeps, so it reads differently from the overlapping scales of the new skin. The sheets on the
ground are height fields on the bed with the same net.

Run from this directory:  python orphic_shedding.py [--draft]
    orphic_shedding.stl   the print, flat at z = 0
"""

import math
import os
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import orphic_egg as oe
import sacred_grid as sg
from orphic_egg import BODY_R, FAR, ROWS, SCALE_LEN, Z_G, smax, smin, smoothstep

FRONT = 450.0               # peel front, mm along the body (about 1.6 turns below the neck)
JAG = ((1, 6.0, 0.4), (2, 3.0, 2.1), (3, 3.5, 1.1), (5, 2.0, 2.9), (8, 1.2, 0.3))
                            # tear line round the body: harmonic, amplitude, phase
LOFT = 1.1                  # old skin's outer surface over the base radius (scale tips 0.68)
TH = 1.3                    # skin thickness where it lifts off: three 0.42 mm walls
FRILL = 17.0                # longest tongue ahead of the front
TONGUES = 4                 # the cuff is torn into this many
FLARE = 7.5                 # lift of the tongue tips off the body
RUFFLE = (7, 2.0)           # waves round the cuff's edge, their amplitude at the tips
COLLAR = (2.6, 4.5, 4.5)    # the skin rolled back at the front: height, centre behind it, half-width
BUNCH = (2.4, 6.0, 16.0)    # concertina folds behind the collar: height, pitch, decay length
CRUMPLE = (1.0, 2.0)        # crumple height along the body, and on the tail lying on the ground
CRUMPLE_WAVES = ((23.0, 1, 0.3), (14.0, -2, 1.9), (31.0, 3, 4.1), (9.5, 2, 2.6), (17.0, -4, 5.2),
                 (11.0, 5, 0.9))     # wavelength along the body, harmonic round it, phase
CREASE_WAVES = ((19.0, -1, 2.2), (12.5, 3, 0.4), (27.0, 2, 3.3), (8.0, -3, 1.4))
FINE_CREASES = ((6.5, 4, 0.2), (5.0, -5, 2.7), (8.5, 6, 1.1), (7.2, -3, 4.4))
NET = (0.35, 0.78)          # scale-outline net: groove depth, and where on the diamond it starts
TEARS, PEELED = 18, 4       # holes in the old skin, how many of them are peeled up as flaps
RIM = 0.6                   # a tear's edge curls up this much
SEED = 7

CURTAINS = ((462.0, 512.0, "next"), (545.0, 585.0, "next"), (600.0, 665.0, "ground"))
                            # drapes of old skin: span of the coil they hang from, and what they
                            # fall to -- the coil beneath, or the ground
DRAPE = ((2.2, 11.0), (1.0, 6.5))   # vertical folds: amplitude, wavelength
BELLY_OUT = 2.5             # a drape bellies out this far between its ends
HEM = 5.0                   # a drape on the ground splays out this far at its foot
K_CURTAIN = 1.0

SHEET_H, PILE = 1.2, 3.4    # sheets on the ground: base thickness, height of their folds
SHEETS = ((802.0, 5.0, 15.0, 8.0, 0.3), (851.0, 6.0, 13.0, 7.5, -0.4), (900.0, 4.0, 11.0, 6.5, 0.5),
          (846.0, -3.0, 10.0, 5.0, 0.2))
                            # s on the tail, centre beyond the tube wall (+ outward), half-length,
                            # half-width, turn off the tail's tangent (radians)
K_SKIN, K_SPILL = 0.4, 1.5  # smooth-union radii onto the body


def waves(sv, phi, table):
    n = sum(np.sin(2 * math.pi * sv / lam + k * phi + c) for lam, k, c in table)
    return np.clip(n / (0.8 * math.sqrt(len(table))), -1.0, 1.0)


def lattice(u, v):
    """Distance in the scale lattice to the nearest diamond's centre: 0 there, 1 on its outline."""
    best = None
    for j in (np.floor(v), np.floor(v) + 1):
        cu = np.round(u - 0.5 * j) + 0.5 * j
        qd = 2 * np.abs(u - cu) + np.abs(v - j)
        best = qd if best is None else np.minimum(best, qd)
    return best


def net(u, v):
    """Groove depth of the scale-outline net on the old skin."""
    return NET[0] * smoothstep(NET[1], 0.98, lattice(u, v))


class Slough:
    """The old skin on the body, in the serpent's tube coordinates."""

    def __init__(self, body):
        self.b = body
        jag = sum(a for _, a, _ in JAG)
        self.j0 = int(np.argmax(body.s >= FRONT - jag - FRILL - 6.0))
        self.tree = cKDTree(body.C[self.j0:])
        self.reach = BODY_R + LOFT + FLARE + RUFFLE[1] + BUNCH[0] + CRUMPLE[1] + 2.0
        self.tears = self._tears()

    def front(self, phi):
        return FRONT + sum(a * np.sin(k * phi + c) for k, a, c in JAG)

    def tips(self, phi):
        """s of the cuff's torn edge: tongues of uneven length, longest on the back."""
        tongue = np.abs(np.sin(TONGUES * phi / 2 + 0.6)) ** 0.5 * (0.75 + 0.25 * np.sin(2 * phi + 1.3))
        back = 0.5 + 0.5 * smoothstep(0.2 * math.pi, 0.7 * math.pi, np.abs(phi))
        return self.front(phi) - FRILL * tongue * back

    def shape(self, sv, phi):
        """Outer and inner radius of the old skin, before its net, tears and rims."""
        rr = np.interp(sv, self.b.s, self.b.r)
        e = self.front(phi) - sv                            # > 0: ahead of the front, in the cuff
        ahead = np.clip(e, 0.0, None) / FRILL
        behind = np.clip(-e, 0.0, None)
        lift = (FLARE + RUFFLE[1] * np.sin(RUFFLE[0] * phi + 0.35 * e + 0.7)) * ahead ** 2
        h, at, half = COLLAR
        roll = h * np.exp(-((behind - at) / half) ** 2) * (0.8 + 0.2 * np.sin(3 * phi + 0.5))
        h, pitch, decay = BUNCH
        past = np.clip(behind - at, 0.0, None)
        fold = h * np.exp(-past / decay) * (0.5 - 0.5 * np.cos(2 * math.pi * past / pitch
                                                               + 0.5 * np.sin(phi + 0.3)))
        amp = CRUMPLE[0] + (CRUMPLE[1] - CRUMPLE[0]) * smoothstep(740.0, 820.0, sv)
        crumple = amp * (0.6 * ((waves(sv, phi, CRUMPLE_WAVES) + 1) / 2) ** 1.6
                         + 0.5 * (1 - np.abs(waves(sv, phi, CREASE_WAVES))) ** 6
                         + 0.35 * (1 - np.abs(waves(sv, phi, FINE_CREASES))) ** 5)
        ro = rr + LOFT + lift + roll + fold + crumple
        ri = rr - 3.0 + (ro - TH - rr + 3.0) * smoothstep(-1.5, 0.5, e)   # solid until it lifts
        return ro, ri

    def _tears(self):
        rng = np.random.default_rng(SEED)
        b, out = self.b, []
        while len(out) < TEARS:
            s = rng.uniform(FRONT + 28.0, b.s[-1] - 10.0)
            phi = rng.choice([-1.0, 1.0]) * rng.uniform(0.35 * math.pi, math.pi)
            rr = np.interp(s, b.s, b.r)
            a = rng.uniform(3.5, 8.0) * np.clip(rr / BODY_R, 0.45, 1.0)
            bb = a * rng.uniform(0.75, 1.0) if len(out) < PEELED else a * rng.uniform(0.45, 0.8)
            near = [t for t in out if math.hypot(s - t[0], ((phi - t[1] + math.pi) % (2 * math.pi) - math.pi)
                                                 * (rr + LOFT)) < a + t[2] + 3.0]
            if near:
                continue
            peel = math.radians(rng.uniform(22.0, 38.0)) if len(out) < PEELED else 0.0
            out.append((s, phi, a, bb, rng.uniform(0, math.pi), peel, rng.uniform(0, 2 * math.pi, 4)))
        return out

    def _tear_fields(self, sv, phi, rho, rr, tex):
        """Distance to the nearest hole (negative inside), and the field of the peeled flaps."""
        holes = np.full(len(sv), 1e3)
        flaps = np.full(len(sv), FAR)
        for s0, p0, a, bb, th, peel, ph in self.tears:
            m = np.nonzero(np.abs(sv - s0) < a + 4.0)[0]
            if not len(m):
                continue
            arc = ((phi[m] - p0 + math.pi) % (2 * math.pi) - math.pi) * (rr[m] + LOFT)
            ds = sv[m] - s0
            x = math.cos(th) * ds + math.sin(th) * arc
            y = -math.sin(th) * ds + math.cos(th) * arc
            ang = np.arctan2(y / bb, x / a)
            wob = sum(w * np.sin(k * ang + c) for k, w, c in zip((2, 3, 5, 7), (0.10, 0.08, 0.06, 0.04), ph))
            d = (np.hypot(x / a, y / bb) - 1 - wob) * bb
            holes[m] = np.minimum(holes[m], d)
            if peel:
                # the torn piece, hinged at one end of the hole and curling up off it
                xi = np.clip(x + a, 0.0, None)
                ro = rr[m] + LOFT + math.tan(peel) * xi + 0.05 * xi * xi
                g = np.sqrt(1 + (math.tan(peel) + 0.1 * xi) ** 2)
                shell = smax(rho[m] - ro + tex[m], ro - TH - rho[m], 0.2) / g
                side = d + 0.8 * (2 * smoothstep(0.0, 4.0, xi) - 1)   # past the hinge, clear elsewhere
                flaps[m] = np.minimum(flaps[m], smax(shell, side, 0.3))
        return holes, flaps

    def local(self, sv, phi, rho, fine):
        b = self.b
        rr = np.interp(sv, b.s, b.r)
        ro, ri = self.shape(sv, phi)
        hs, hp = 0.05, 0.01                                 # slope of each radius, for true distances
        s1, s0 = self.shape(sv + hs, phi), self.shape(sv - hs, phi)
        p1, p0 = self.shape(sv, phi + hp), self.shape(sv, phi - hp)
        go = np.sqrt(1 + ((s1[0] - s0[0]) / (2 * hs)) ** 2 + ((p1[0] - p0[0]) / (2 * hp * ro)) ** 2)
        gi = np.sqrt(1 + ((s1[1] - s0[1]) / (2 * hs)) ** 2
                     + ((p1[1] - p0[1]) / (2 * hp * np.maximum(ri, 1.0))) ** 2)
        tex = net(np.interp(sv, b.s, b.u), phi / (2 * math.pi) * ROWS) if fine else np.zeros_like(sv)
        holes, flaps = self._tear_fields(sv, phi, rho, rr, tex)
        rim = RIM * smoothstep(2.5, 0.0, holes)
        d = smax((rho - ro - rim + tex) / go, (ri - rho) / gi, 0.25)
        g, gp = self.tips(phi), (self.tips(phi + hp) - self.tips(phi - hp)) / (2 * hp)
        d = smax(d, (g - sv) / np.sqrt(1 + (gp / (rr + LOFT)) ** 2), 0.35)   # the cuff's torn edge
        d = smax(d, -holes, 0.35)
        return np.minimum(d, flaps)

    def sdf(self, P, fine=True):
        b = self.b
        out = np.full(len(P), FAR)
        dist, j = self.tree.query(P, distance_upper_bound=self.reach, workers=-1)
        hit = np.isfinite(dist)
        i, p = j[hit] + self.j0, P[hit]
        w = p - b.C[i]
        t = np.einsum("ij,ij->i", w, b.T[i])
        t = np.clip(t, np.where(i == self.j0, 0.0, -b.ds), np.where(i == len(b.C) - 1, 0.0, b.ds))
        w -= t[:, None] * b.T[i]
        phi = np.arctan2(np.einsum("ij,ij->i", w, b.B[i]), np.einsum("ij,ij->i", w, b.N[i]))
        out[hit] = self.local(b.s[i] + t, phi, np.linalg.norm(w, axis=1), fine)
        return out


class Curtain:
    """A drape of old skin hanging from the outer flank of a coil, down to the coil beneath
    or to the ground, in cylindrical coordinates about the coil's axis: azimuth psi, height z,
    radius q. Both ends start inside a tube, so the drape grows out of the coil it hangs from
    and into the one it lands on, and it prints as a near-vertical wall."""

    Q_REF = 40.0                # radius at which psi is turned into arc length

    def __init__(self, body, s0, s1, to, rng):
        az = np.unwrap(np.arctan2(body.C[:, 1], body.C[:, 0]))
        k = np.nonzero((body.s >= s0) & (body.s <= s1))[0][::-1]       # azimuth falls along s
        self.th = az[k]
        self.qu = np.hypot(body.C[k, 0], body.C[k, 1]) + 0.6 * body.r[k]
        self.zu = body.C[k, 2]
        if to == "ground":
            self.ql, self.zl = self.qu + HEM, np.full(len(k), -2.0)
        else:
            below = np.nonzero(body.s > s1)[0][::-1]
            sl = np.interp(self.th - 2 * math.pi, az[below], body.s[below])
            cl = np.c_[[np.interp(sl, body.s, body.C[:, j]) for j in range(3)]].T
            self.ql = np.hypot(cl[:, 0], cl[:, 1]) + 0.6 * np.interp(sl, body.s, body.r)
            self.zl = cl[:, 2]
        self.ground = to == "ground"
        self.ph = rng.uniform(0, 2 * math.pi, 8)
        span = self.th[-1] - self.th[0]
        self.holes = [(self.th[0] + span * rng.uniform(0.25, 0.75), rng.uniform(0.3, 0.6),
                       rng.uniform(3.5, 5.5), rng.uniform(2.5, 4.0)) for _ in range(2)]

    def ends(self, psi):
        return [np.interp(psi, self.th, a) for a in (self.qu, self.zu, self.ql, self.zl)]

    def radius(self, psi, z):
        qu, zu, ql, zl = self.ends(psi)
        f = np.clip((zu - z) / (zu - zl), 0.0, 1.0)
        arc = psi * self.Q_REF
        ph = self.ph
        folds = sum(a * np.sin(2 * math.pi * arc / lam + p + 0.9 * np.sin(2 * math.pi * z / w + q))
                    for (a, lam), p, q, w in zip(DRAPE, ph[:2], ph[2:4], (19.0, 11.0)))
        folds += 0.8 * (1 - np.abs(np.sin(2 * math.pi * arc / 15.0 + 0.7 * z / 6.0 + ph[4]))) ** 6
        q = qu + (ql - qu) * f + BELLY_OUT * np.sin(math.pi * f) + folds * np.sin(math.pi * f) ** 0.5
        if self.ground:
            q = q + HEM * smoothstep(8.0, 0.0, z) ** 2      # the foot splays onto the bed
        return q, f, zu, zl

    def sdf(self, P, fine=True):
        out = np.full(len(P), FAR)
        mid = 0.5 * (self.th[0] + self.th[-1])
        th = np.arctan2(P[:, 1], P[:, 0])
        psi = mid + (th - mid + math.pi) % (2 * math.pi) - math.pi
        m = np.nonzero((psi > self.th[0] - 0.15) & (psi < self.th[-1] + 0.15)
                       & (P[:, 2] < self.zu.max() + 2.0) & (P[:, 2] > self.zl.min() - 2.0))[0]
        if not len(m):
            return out
        psi, z, q = psi[m], P[m, 2], np.hypot(P[m, 0], P[m, 1])
        qc, f, zu, zl = self.radius(psi, z)
        e, ez = 0.004, 0.05
        dpsi = (self.radius(psi + e, z)[0] - self.radius(psi - e, z)[0]) / (2 * e * np.maximum(q, 1.0))
        dz = (self.radius(psi, z + ez)[0] - self.radius(psi, z - ez)[0]) / (2 * ez)
        g = np.sqrt(1 + dpsi ** 2 + dz ** 2)
        tex = net(psi * self.Q_REF / (SCALE_LEN * BODY_R), z * ROWS / (2 * math.pi * BODY_R)) if fine else 0.0
        d = smax((q - qc - TH / 2 + tex) / g, (qc - TH / 2 - q) / g, 0.2)
        # torn sides, wandering down the drape
        ph = self.ph
        lo = self.th[0] + (4.0 + 3.0 * np.sin(9 * f + ph[5]) + 1.5 * np.sin(23 * f + ph[6])) / self.Q_REF
        hi = self.th[-1] - (4.0 + 3.0 * np.sin(8 * f + ph[7]) + 1.5 * np.sin(19 * f + ph[5])) / self.Q_REF
        d = smax(d, np.maximum(lo - psi, psi - hi) * q, 0.4)
        for c, fh, a, b in self.holes:
            zc = np.interp(c, self.th, self.zu) * (1 - fh) + np.interp(c, self.th, self.zl) * fh
            k = np.hypot((psi - c) * q / a, (z - zc) / b)
            wob = 0.12 * np.sin(3 * np.arctan2(z - zc, (psi - c) * q) + ph[0])
            d = smax(d, -(k - 1 - wob) * min(a, b), 0.35)
        out[m] = np.maximum(d, np.maximum(z - zu, zl - z))
        return out


class Spill:
    """Crumpled sheets of the old skin slumped onto the ground beside the tail: height fields
    on the bed, each a ragged oval with folds and the scale-outline net."""

    def __init__(self, body):
        rng = np.random.default_rng(SEED + 1)
        self.sheets = []
        for s, off, a, bb, turn in SHEETS:
            i = int(np.argmin(np.abs(body.s - s)))
            c = body.C[i, :2]
            out = c / np.linalg.norm(c)
            t = body.T[i, :2] / np.linalg.norm(body.T[i, :2])
            ang = math.atan2(t[1], t[0]) + turn
            centre = c + out * math.copysign(body.r[i] + abs(off), off)
            self.sheets.append((centre, ang, a, bb, rng.uniform(0, 2 * math.pi, 4)))
        k = rng.uniform(0, 2 * math.pi, 6)
        self.folds = [(math.cos(q), math.sin(q), lam, rng.uniform(0, 2 * math.pi))
                      for q, lam in zip(k, (8.0, 11.0, 5.5, 14.0, 9.0, 6.5))]

    def height(self, x, y):
        n = sum(np.sin(2 * math.pi * (cx * x + cy * y) / lam + c) for cx, cy, lam, c in self.folds)
        n = np.clip(n / (0.8 * math.sqrt(len(self.folds))), -1, 1)
        return 0.5 * ((n + 1) / 2) ** 2 + (1 - np.abs(n)) ** 4

    def sdf(self, P, fine=True):
        out = np.full(len(P), FAR)
        m = np.nonzero(P[:, 2] < SHEET_H + PILE + 3.0)[0]
        x, y, z = P[m].T
        d2 = np.full(len(m), 1e3)
        for c, ang, a, bb, ph in self.sheets:
            lx = math.cos(ang) * (x - c[0]) + math.sin(ang) * (y - c[1])
            ly = -math.sin(ang) * (x - c[0]) + math.cos(ang) * (y - c[1])
            th = np.arctan2(ly / bb, lx / a)
            wob = sum(w * np.sin(k * th + p) for k, w, p in zip((2, 3, 5, 9), (0.12, 0.10, 0.07, 0.04), ph))
            d2 = np.minimum(d2, (np.hypot(lx / a, ly / bb) - 1 - wob) * bb)
        e = 0.05
        pile = PILE * smoothstep(0.0, 3.5, -d2)             # folds die away to a thin ragged edge
        h = SHEET_H + pile * self.height(x, y)
        gx = pile * (self.height(x + e, y) - self.height(x - e, y)) / (2 * e)
        gy = pile * (self.height(x, y + e) - self.height(x, y - e)) / (2 * e)
        if fine:
            h = h - net(x / (SCALE_LEN * BODY_R), y * ROWS / (2 * math.pi * BODY_R))
        out[m] = smax(d2, (z - h) / np.sqrt(1 + gx * gx + gy * gy), 0.5)
        return out


class Shedding(oe.Piece):
    """The empty coil with its old skin coming away."""

    def __init__(self):
        super().__init__(with_egg=False)
        self.skin = Slough(self.body)
        self.spill = Spill(self.body)
        rng = np.random.default_rng(SEED + 2)
        self.curtains = [Curtain(self.body, s0, s1, to, rng) for s0, s1, to in CURTAINS]

    def field(self, P, fine=True):
        d = smin(self.body.sdf(P, fine), self.skin.sdf(P, fine), K_SKIN)
        for c in self.curtains:
            d = smin(d, c.sdf(P, fine), K_CURTAIN)
        d = smin(d, self.spill.sdf(P, fine), K_SPILL)
        d = smin(d, self.head.sdf(P), oe.K_HEAD)
        d = smin(d, self.key.sdf(P), oe.K_KEY)
        return np.maximum(d, Z_G - P[:, 2])

    def bounds(self):
        lo, hi = super().bounds()
        g = self.skin.reach - BODY_R
        return lo - [g, g, 0.0], hi + [g, g, g]


def main():
    t0 = time.time()
    piece = Shedding()
    sk = piece.skin
    print(f"    old skin from s = {FRONT:.0f} mm to the tail tip ({piece.body.s[-1]:.0f} mm), "
          f"{TEARS} tears, {PEELED} peeled up, {len(SHEETS)} sheets on the ground")
    m = oe.mesh(piece)
    print(f"    meshed in {time.time() - t0:.0f} s")
    m = oe.decimate(m, oe.TARGET_FACES)
    V = m.vertices.astype(np.float32).astype(float)
    V[:, 2] = np.maximum(V[:, 2], Z_G)                 # decimation nudges a few microns below the bed
    m = trimesh.Trimesh(V, m.faces)
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    m = oe.carve_pupils(piece, m)
    path = "orphic_shedding_draft.stl" if oe.DRAFT else "orphic_shedding.stl"
    m.export(path)
    back = trimesh.load(path)
    ok = sg.check(back, path)
    tip, area, contact, c = oe.stability(back)
    print(f"    on the bed {contact:.0f} mm2 within a {area:.0f} mm2 footprint, "
          f"centre of mass {c.round(1)}, tips at {tip:.1f} deg")
    print(f"    done in {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
