#!/usr/bin/env python3
"""The Orphic Egg: the world-egg with a serpent wound three times round it, head raised
over the crown. Sculpted after the reference pendant (reference_*.webp) without its hanger
loop, and given a foot so it stands: the serpent's last turn slides off the egg onto the
ground and its tail lies round the egg's tip.

Sized from the shared grid (lib/sacred_grid.py, R = 15 mm):
    egg       8R tall x 5R wide -- Fibonacci 8:5, near phi -- blunt end up
    serpent   body diameter R, coils 2.25R apart along the egg's meridian
    head      about 2R long, carried 1.2R above the crown
    pupils    exact vesicae piscis, each holding a vortex: a spiral winding into its inscribed circle
    key       a skeleton key bitten crosswise, its bow the Tripod of Life in a ring

The whole piece is one signed-distance field -- egg, serpent tube with imbricate dorsal
scales and belly plates, sculpted head -- joined by smooth unions, so every coil is
filleted onto the egg and the neck grows out of the skull. It is meshed by marching cubes
in narrow-band slabs (exact field only within BAND of the surface, a trilinear coarse
field elsewhere), Taubin-smoothed, and quadric-decimated to stay under the repo's 50 MB STL
limit. The pupils are too fine for that voxel; they are meshed separately and cut out last.

Run from this directory:  python orphic_egg.py [--draft] [--serpent]
    orphic_egg.stl       the print, flat at z = 0 (stands on the egg tip and the tail)
    orphic_serpent.stl   --serpent: the same serpent with no egg inside its coils
"""

import math
import os
import sys
import time

import numpy as np
import trimesh
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes
from trimesh.smoothing import filter_taubin

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import sacred_grid as sg

R = sg.R

EGG_H = 8 * R
EGG_W = 5 * R
EGG_SHIFT = 0.05 * EGG_H    # Hugelschaffer offset: widest section moved toward the blunt end
EGG_SINK = 3.0              # tip pressed this far below the bed, for a contact disc

BODY_R = R / 2
EMBED = 0.9                 # coil centreline sits this much nearer the egg than BODY_R
PITCH = 2.25 * R            # coil spacing along the meridian
MU_TOP = 23.0               # top coil: meridian distance below the crown
ROOT_DEG = 140.0            # azimuth where the neck leaves the top coil (front is -90)
CLIMB, CLIMB_RAD = 4.0, 1.2 # top coil climbs the crown this much over its last radians
NECK_K = (1.0, 1.8)         # neck curve tangent lengths at the head and root, x chord;
                            # tightest visible bend then 1.6 body radii
TAIL_TURNS = 0.85           # turns of tail on the ground after the body lands
TAPER_TURNS = 0.3           # tail taper starts this many turns before landing
TAPER_POW = 0.9
TIP_R = 1.2
GROUND = 0.82               # tail centreline height / its radius: 1.14 r contact strip
CURL = 12.0                 # tail tip curls this far out from the egg

ROWS = 26                   # dorsal scale rows round the body (even, so the stagger closes)
SCALE_LEN = 0.56            # scale length / body radius
SCALE_H = 0.09              # scale relief / body radius
BELLY = math.radians(40)    # half-width of the ventral plate band
SEAM = math.radians(5)      # groove between belly plates and dorsal scales

HEAD_AT = (-3.0, 3.0)       # back of the skull, x y
HEAD_LIFT = 20.0            # back of the skull above the crown
HEAD_SCALE = 1.2            # Head's sculpt is drawn at 1/HEAD_SCALE
HEAD_DIR = (1.0, -0.45, -0.12)

K_EGG, K_FOOT, K_HEAD = 1.5, 3.0, 3.0     # smooth-union radii; K_FOOT fills the crevice at the bed
K_KEY = 0.6
Z_G = 0.0
FAR = 6.0                   # stand-in field value well clear of every surface

DRAFT = "--draft" in sys.argv
SERPENT = "--serpent" in sys.argv   # leave the egg out: the serpent alone, coiled round where it was
H = 0.42 if DRAFT else 0.2  # marching-cubes voxel
HC = 1.0                    # coarse voxel
BAND = 2.4                  # exact evaluation within this distance of the coarse surface
TAUBIN = 8                  # smoothing passes: take marching-cubes stair-steps off cut edges
TARGET_FACES = 850_000


def smoothstep(a, b, x):
    t = np.clip((np.asarray(x, float) - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b + (a - b) * h - k * h * (1 - h)


def smax(a, b, k):
    return -smin(-a, -b, k)


def rounded(x):
    """Quarter-round edge profile: 0 at the edge, 1 once x reaches 1."""
    x = np.clip(x, 0.0, 1.0)
    return np.sqrt(1 - (1 - x) ** 2)


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def ellipsoid(x, y, z, a, b, c):
    k0 = np.sqrt((x / a) ** 2 + (y / b) ** 2 + (z / c) ** 2)
    k1 = np.sqrt((x / a / a) ** 2 + (y / b / b) ** 2 + (z / c / c) ** 2)
    return k0 * (k0 - 1) / np.maximum(k1, 1e-9)


class Egg:
    """Hugelschaffer egg of revolution about z, pointed end down, as an exact 2-D distance
    to its meridian."""

    def __init__(self):
        L, w = EGG_H, EGG_SHIFT
        t = np.linspace(0, np.pi, 20001)
        x = -0.5 * L * np.cos(t)                         # blunt end -L/2 .. pointed end +L/2
        q = np.sqrt(np.clip((L * L - 4 * x * x) / (L * L + 8 * w * x + 4 * w * w), 0, None))
        q *= 0.5 * EGG_W / q.max()
        z = Z_G - EGG_SINK + 0.5 * L - x
        m = np.r_[0, np.cumsum(np.hypot(np.diff(q), np.diff(z)))]
        self.mu_len = m[-1]
        self.mu = np.linspace(0, m[-1], int(m[-1] / 0.04) + 1)   # meridian length from the crown
        self.q = np.interp(self.mu, m, q)
        self.z = np.interp(self.mu, m, z)
        tq, tz = np.gradient(self.q, self.mu), np.gradient(self.z, self.mu)
        n = np.hypot(tq, tz)
        self.nq, self.nz = -tz / n, tq / n                  # outward normal in the meridian plane
        self.tree = cKDTree(np.c_[self.q, self.z])
        self.top, self.bottom = self.z[0], self.z[-1]

    def sdf2(self, q, z):
        d, i = self.tree.query(np.c_[q, z], distance_upper_bound=12.0, workers=-1)
        d = np.minimum(d, 12.0)
        inside = (z > self.bottom) & (z < self.top) & (q < np.interp(z, self.z[::-1], self.q[::-1]))
        return np.where(inside, -d, d)

    def sdf(self, P):
        return self.sdf2(np.hypot(P[:, 0], P[:, 1]), P[:, 2])

    def outward(self, P):
        """Unit outward normal of the egg at its point nearest each P."""
        q = np.hypot(P[:, 0], P[:, 1])
        _, i = self.tree.query(np.c_[q, P[:, 2]], workers=-1)
        th = np.arctan2(P[:, 1], P[:, 0])
        return np.c_[self.nq[i] * np.cos(th), self.nq[i] * np.sin(th), self.nz[i]]

    def at(self, mu, off):
        """(q, z) of the point `off` out along the normal from meridian position mu."""
        mu = np.clip(mu, 0.0, self.mu_len)
        return (np.interp(mu, self.mu, self.q) + off * np.interp(mu, self.mu, self.nq),
                np.interp(mu, self.mu, self.z) + off * np.interp(mu, self.mu, self.nz))

    def hug(self, z, off):
        """Radius q at height z where the egg's distance is `off` (bisection)."""
        lo, hi = np.zeros_like(z), np.full_like(z, 120.0)
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            out = self.sdf2(mid, z) > off
            hi = np.where(out, mid, hi)
            lo = np.where(out, lo, mid)
        return 0.5 * (lo + hi)


class Head:
    """Sculpted head in its own frame: f forward, l left, u up, origin at the back of the
    skull. Mirror-symmetric, so eyes and nostrils are placed once and used at |l|."""

    # The eye is built in its own skin frame: x along the head, y up, z out of the skin.
    EYE = (13.0, 1.8)                   # f, u of the eye's centre on the skull surface
    BALL = (2.6, 2.3, 2.1, 1.25)        # eyeball semi-axes x y z, and its centre's depth
    SOCKET = (3.8, 3.4, 1.6, 0.6)       # orbit carved round it: semi-axes, centre above the skin
    VESICA_H = 1.8                      # pupil: an exact vesica piscis this half-tall
    PUPIL_DEPTH = 0.5                   # floor of the vesica's horns, below the eyeball
    DISC_DEPTH, FUNNEL = 0.15, 0.35     # its inscribed circle: rim depth, extra fall to the centre
    SPIRAL = (1.5, 0.375, 0.3)          # vortex groove in that circle: turns, width, depth
    BROW = ((3.6, 1.6, 0.6, -0.25),     # supraocular shelf: x, y, radius, lift off the skin (sunk,
            (2.0, 2.9, 1.3, -0.5),      # so it stands 0.35-1.0 proud), swept from in front of the
            (0.0, 3.4, 1.6, -0.6),      # eye, over it, and back; the orbit cuts its lower edge
            (-2.4, 3.3, 1.6, -0.6),     # into the lid line
            (-4.8, 2.7, 1.2, -0.5),
            (-7.2, 1.6, 0.5, -0.15))
    NOSTRIL = (26.8, 0.5)
    EDGE = 0.2                          # round-over on cut edges, so marching cubes can't saw them
    BOX = ((-5.0, 33.0), (-12.0, 12.0), (-9.0, 9.0))

    def __init__(self, origin, forward):
        F = unit(forward)
        U = unit(np.array([0.0, 0.0, 1.0]) - F[2] * F)
        self.O, self.F, self.U = np.asarray(origin, float), F, U
        self.M = np.stack([F, np.cross(U, F), U])          # rows: f, l, u axes
        s = self._landmark(*self.EYE)
        n = self._normal(s)
        x = unit(np.array([1.0, 0.0, 0.0]) - n[0] * n)
        y = np.cross(x, n)
        self.eye_at, self.eye_axes = s, np.stack([x, y, n])
        brow = []
        for bx, by, r, lift in self.BROW:
            p = self._onto_skin(s + bx * x + by * y)
            brow.append(np.r_[p + lift * self._normal(p), r])
        self.brow = np.array(brow)
        sn = self._landmark(*self.NOSTRIL)
        self.nostril = sn + 0.15 * self._normal(sn)

    @staticmethod
    def base(f, l, u):
        skull = ellipsoid(f - 9.5, l, u - 0.6, 10.5, 8.8, 6.4)
        snout = ellipsoid(f - 19.2, l, u + 0.2, 9.2, 6.3, 4.7)
        jaw = ellipsoid(f - 8.5, l, u + 2.0, 9.5, 7.6, 4.4)
        h = smin(smin(skull, snout, 4.0), jaw, 2.0)
        return smax(h, u - (5.7 - 0.1 * f), 2.5)           # flat crown, falling to the snout

    def _landmark(self, f, u):
        lo, hi = 0.0, 14.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if self.base(np.array(f), np.array(mid), np.array(u)) > 0:
                hi = mid
            else:
                lo = mid
        assert lo > 2.0, f"landmark ({f}, {u}) is off the side of the head"
        return np.array([f, 0.5 * (lo + hi), u])

    def _normal(self, p, e=1e-3):
        g = [(self.base(*(p + d)) - self.base(*(p - d))) / (2 * e) for d in np.eye(3) * e]
        return unit(np.array(g, float))

    def _onto_skin(self, p):
        for _ in range(8):
            p = p - self.base(*p) * self._normal(p)
        return p

    def local(self, f, l, u):
        h = self.base(f, l, u)
        la = np.abs(l)

        p = np.stack([f, la, u], axis=-1)
        x, y, z = ((p - self.eye_at) @ self.eye_axes.T).T
        h = smin(h, capsules(p, self.brow), 1.0)
        a, b, c, lift = self.SOCKET
        h = smax(h, -ellipsoid(x, y, z - lift, a, b, c), 0.25)    # orbit
        a, b, c, sink = self.BALL
        ball = ellipsoid(x, y, z + sink, a, b, c)
        h = smin(h, ball, 0.12)                                   # crisp lid line round the ball
        h = smax(h, -(np.linalg.norm(p - self.nostril, axis=-1) - 0.75), 0.25)
        # mouth: a slit round the head, 1.1 mm deep
        half = 0.42 * smoothstep(3.5, 7.5, f)
        return smax(h, np.minimum(half - np.abs(u - self.lip(f)), h + 1.1), self.EDGE)

    @staticmethod
    def lip(f):
        """Height u of the mouth line, rising toward the gape."""
        return -1.5 + 0.9 * (1 - np.clip((f - 4.0) / 25.0, 0, 1)) ** 2

    def pupil(self, P):
        """Field of the solid cut from each eyeball, world units, negative inside. The pupil is
        an exact vesica piscis -- two circles of radius rho, each through the other's centre, so
        rho wide and rho*sqrt(3) tall -- and its inscribed circle is a vortex: a disc falling
        toward its centre, with a spiral groove winding inward from its rim to a pit."""
        loc = (P - self.O) @ self.M.T / HEAD_SCALE
        loc[:, 1] = np.abs(loc[:, 1])
        x, y, z = ((loc - self.eye_at) @ self.eye_axes.T).T
        rho = 2 * self.VESICA_H / math.sqrt(3)
        ves = np.maximum(np.hypot(x + rho / 2, y), np.hypot(x - rho / 2, y)) - rho
        a, r = rho / 2, np.hypot(x, y)
        turns, width, groove = self.SPIRAL
        pitch = a / turns
        t = (np.arctan2(y, x) / (2 * math.pi)) % 1.0
        k0 = np.round((a - r) / pitch - t)
        gap = np.full_like(r, np.inf)
        for k in (k0 - 1, k0, k0 + 1):                        # arms r = a - pitch * (t + k)
            run = t + k
            gap = np.where((run >= 0) & (run <= turns), np.minimum(gap, np.abs(r - a + pitch * run)), gap)
        disc = self.DISC_DEPTH + self.FUNNEL * (1 - r / a) + groove * np.clip((width / 2 - gap) / 0.05, 0, 1)
        rim = smoothstep(a - 0.04, a + 0.04, r)
        depth = (1 - rim) * disc + rim * self.PUPIL_DEPTH
        bx, by, bz, sink = self.BALL
        ball = ellipsoid(x, y, z + sink, bx, by, bz)
        return HEAD_SCALE * np.maximum(np.maximum(ves, -(ball + depth)), np.maximum(z - 1.6, -1.0 - z))

    def sdf(self, P):
        loc = (P - self.O) @ self.M.T / HEAD_SCALE
        (f0, f1), (l0, l1), (u0, u1) = self.BOX
        f, l, u = loc[:, 0], loc[:, 1], loc[:, 2]
        inn = (f > f0) & (f < f1) & (l > l0) & (l < l1) & (u > u0) & (u < u1)
        out = np.full(len(P), FAR)
        out[inn] = HEAD_SCALE * self.local(f[inn], l[inn], u[inn])
        return out


def extrude(d2, t, half, r):
    """A 2-D field d2 given thickness 2*half along t, all edges rounded by r."""
    x, y = d2 + r, np.abs(t) - half + r
    return np.minimum(np.maximum(x, y), 0) + np.hypot(np.maximum(x, 0), np.maximum(y, 0)) - r


def box2(x, y, cx, cy, hx, hy):
    qx, qy = np.abs(x - cx) - hx, np.abs(y - cy) - hy
    return np.minimum(np.maximum(qx, qy), 0) + np.hypot(np.maximum(qx, 0), np.maximum(qy, 0))


class Key:
    """A decorative skeleton key held crosswise in the jaws, like a stick. The bow stands on
    the serpent's right (toward the front view), the bit hangs on its left. Key frame, world mm:
    a along the shaft toward the bit, w up, t through the plate (the head's forward); origin
    where the shaft crosses the head's midplane on the mouth line.

    The bow is the Tripod of Life -- three circles, each through the other two's centres --
    inside a ring tangent to their three outer lobes."""

    BITE = 20.0                 # head-local f where the jaws close on the shaft
    SHAFT_R = 1.8
    PLATE = 1.3                 # half-thickness of bow and bit
    BOW_R, BAND = 4.6, 0.75     # bow circles: centreline radius, half band width
    COLLAR = (2.8, 1.1)         # throat collar between bow and shaft: radius, half-length
    GAP = 2.5                   # collar clear of the cheek
    BIT = (6.5, 7.0)            # bit length along the shaft, drop below it
    WARDS = ((1.9, 0.5, 2.8), (4.5, 0.5, 1.8))   # bottom slots: from bit start, half-width, depth
    TIP = 1.5                   # shaft beyond the bit

    def __init__(self, head):
        f, u = self.BITE, Head.lip(self.BITE)
        lo, hi = 0.0, 14.0
        for _ in range(50):                             # jaw half-width at the bite
            mid = 0.5 * (lo + hi)
            lo, hi = (lo, mid) if head.base(np.array(f), np.array(mid), np.array(u)) > 0 else (mid, hi)
        self.cheek = HEAD_SCALE * lo
        self.O = head.O + HEAD_SCALE * (f * head.F + u * head.U)
        self.A = np.stack([head.M[1], head.U, head.F])     # rows: a, w, t
        self.ring = self.BOW_R / math.sqrt(3)               # bow circle centres from the bow centre
        self.rim = self.ring + self.BOW_R                   # enclosing ring, tangent to the lobes
        reach = self.rim + self.BAND                        # bow edge toward the shaft
        self.collar_at = -(self.cheek + self.GAP + self.COLLAR[1])
        self.bow_at = self.collar_at - self.COLLAR[1] - reach
        self.bit_at = self.cheek + 1.5
        self.tip_at = self.bit_at + self.BIT[0] + self.TIP
        h = max(self.BIT[1], reach) + 1
        self.box = ((self.bow_at - reach - 1, self.tip_at + self.SHAFT_R + 1), (-h, h), (-4.0, 4.0))

    def sdf(self, P):
        q = (P - self.O) @ self.A.T
        (a0, a1), (w0, w1), (t0, t1) = self.box
        inn = (q[:, 0] > a0) & (q[:, 0] < a1) & (q[:, 1] > w0) & (q[:, 1] < w1) & (q[:, 2] > t0) & (q[:, 2] < t1)
        out = np.full(len(P), FAR)
        a, w, t = q[inn].T
        bow = np.full_like(a, np.inf)
        for th in (math.pi, math.pi / 3, -math.pi / 3):
            ring = np.abs(np.hypot(a - self.bow_at - self.ring * math.cos(th), w - self.ring * math.sin(th)) - self.BOW_R)
            bow = np.minimum(bow, ring - self.BAND)
        bow = np.minimum(bow, np.abs(np.hypot(a - self.bow_at, w) - self.rim) - self.BAND)
        key = extrude(bow, t, self.PLATE, 0.35)
        rc, hl = self.COLLAR
        key = smin(key, extrude(np.hypot(w, t) - rc, a - self.collar_at, hl, 0.5), 0.4)
        ac = np.clip(a, self.bow_at + self.rim, self.tip_at)
        key = smin(key, np.sqrt((a - ac) ** 2 + w * w + t * t) - self.SHAFT_R, 0.4)
        length, drop = self.BIT
        bit = box2(a, w, self.bit_at + length / 2, -drop / 2, length / 2, drop / 2)
        for at, hw, depth in self.WARDS:
            bit = np.maximum(bit, -box2(a, w, self.bit_at + at, -drop, hw, depth))
        bit = np.maximum(bit, -box2(a, w, self.bit_at + length, -0.55 * drop, 0.9, 0.6))   # end ward
        out[inn] = smin(key, extrude(bit, t, self.PLATE, 0.3), 0.4)
        return out

    def corners(self):
        (a0, a1), (w0, w1), (t0, t1) = self.box
        c = np.array([[a, w, t] for a in (a0, a1) for w in (w0, w1) for t in (t0, t1)])
        return self.O + c @ self.A


def capsules(p, path):
    """Distance to a chain of cones through path rows (x, y, z, radius)."""
    d = np.full(len(p), np.inf)
    for (a, ra), (b, rb) in zip(((q[:3], q[3]) for q in path[:-1]), ((q[:3], q[3]) for q in path[1:])):
        ab = b - a
        t = np.clip((p - a) @ ab / (ab @ ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(p - a - t[:, None] * ab, axis=1) - (ra + t * (rb - ra)))
    return d


def hermite(a, ta, b, tb, n):
    t = np.linspace(0, 1, n)[:, None]
    h00, h10 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t
    h01, h11 = -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
    return h00 * a + h10 * ta + h01 * b + h11 * tb


def resample(P, r, ds):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    si = np.linspace(0, s[-1], int(s[-1] / ds) + 1)
    return np.c_[[np.interp(si, s, P[:, k]) for k in range(3)]].T, np.interp(si, s, r), si


class Serpent:
    """Centreline from inside the skull to the tail tip, with radius, a belly-facing frame
    and the scale coordinate u (scale lengths from the head)."""

    def __init__(self, egg, head):
        two_pi = 2 * math.pi
        D = np.arange(0.0, 6 * two_pi, 0.002)                 # turns down from the root, radians
        zh = egg.at(MU_TOP + PITCH * D / two_pi, BODY_R - EMBED)[1]
        land = D[np.argmax(zh < Z_G + GROUND * BODY_R)]
        self.land_turns = land / two_pi
        end = land + TAIL_TURNS * two_pi
        D = D[D <= end]
        x = np.clip((D - (land - TAPER_TURNS * two_pi)) / (end - land + TAPER_TURNS * two_pi), 0, 1)
        r = TIP_R + (BODY_R - TIP_R) * (1 - x) ** TAPER_POW
        off = r * (1 - EMBED / BODY_R)
        z = smax(egg.at(MU_TOP + PITCH * D / two_pi, off)[1], Z_G + GROUND * r, 4.0)
        z += CLIMB * np.clip(1 - D / CLIMB_RAD, 0, 1) ** 2    # rise toward the neck
        q = egg.hug(z, off)
        # once down, keep the thinning tail from creeping under the egg, then curl it out
        q_land = np.interp(land, D, q)
        floor = q_land + CURL * smoothstep(end - 0.3 * two_pi, end, D) ** 1.5
        on = smoothstep(land - 0.5, land + 0.5, D)
        q = (1 - on) * q + on * smax(q, floor, 2.0)
        th = math.radians(ROOT_DEG) - D
        coil = np.c_[q * np.cos(th), q * np.sin(th), z]

        a, b = head.O + 6.0 * HEAD_SCALE * head.F, head.O
        root, t_root = coil[0], unit(coil[1] - coil[0])
        chord = np.linalg.norm(root - b)
        neck = hermite(b, -NECK_K[0] * chord * head.F, root, NECK_K[1] * chord * t_root, 400)
        inner = a + np.linspace(0, 1, 40)[:, None] * (b - a)
        P = np.r_[inner, neck[1:], coil[1:]]
        rr = np.r_[np.full(len(inner) + len(neck) - 1, BODY_R), r[1:]]

        ds = 0.15
        P, rr, s = resample(P, rr, ds)
        P[1:-1] = gaussian_filter1d(P, 1.5 / ds, axis=0, mode="nearest")[1:-1]
        P, rr, s = resample(P, rr, ds)
        rr *= 0.62 + 0.18 * smoothstep(0.0, 7.0, s) + 0.2 * smoothstep(6.0, 34.0, s)   # neck narrows
                                                     # behind the head and ends well inside it
        self.C, self.r, self.s, self.ds = P, rr, s, ds

        T = unit(np.gradient(P, axis=0))
        toward_egg = -egg.outward(P)
        on_ground = smoothstep(rr + 9.0, rr + 2.0, P[:, 2] - Z_G)[:, None]
        at_head = smoothstep(20.0, 8.0, s)[:, None]
        belly = (1 - at_head) * ((1 - on_ground) * toward_egg + on_ground * [0, 0, -1]) \
            + at_head * -head.U
        N = unit(belly - np.sum(belly * T, axis=1, keepdims=True) * T)
        N = gaussian_filter1d(N, 3.0 / ds, axis=0, mode="nearest")
        N = unit(N - np.sum(N * T, axis=1, keepdims=True) * T)
        self.T, self.N, self.B = T, N, np.cross(T, N)
        self.u = np.r_[0, np.cumsum(ds / (SCALE_LEN * np.maximum(rr[1:], 2.0)))]
        self.tree = cKDTree(P)
        self.egg_gap = egg.sdf(P) - rr              # < 0 where the tube sinks into the egg

    def clearance(self):
        """Smallest wall-to-wall gap between parts of the body more than 4R apart along it
        (negative: the tube runs into itself)."""
        pairs = self.tree.query_pairs(4 * BODY_R + 6.0, output_type="ndarray")
        i, j = pairs[:, 0], pairs[:, 1]
        far = np.abs(self.s[i] - self.s[j]) > 4 * R
        i, j = i[far], j[far]
        gap = np.linalg.norm(self.C[i] - self.C[j], axis=1) - self.r[i] - self.r[j]
        k = np.argmin(gap) if len(gap) else None
        return (gap[k], self.s[i[k]], self.s[j[k]]) if k is not None else (np.inf, 0, 0)

    def sdf(self, P, scales=True):
        dist, i = self.tree.query(P, distance_upper_bound=BODY_R + FAR, workers=-1)
        out = np.full(len(P), FAR)
        hit = np.isfinite(dist)
        i, p = i[hit], P[hit]
        w = p - self.C[i]
        t = np.einsum("ij,ij->i", w, self.T[i])
        lo = np.where(i == 0, 0.0, -self.ds)                 # both ends capped
        hi = np.where(i == len(self.C) - 1, 0.0, self.ds)
        t = np.clip(t, lo, hi)
        w -= t[:, None] * self.T[i]
        rho = np.linalg.norm(w, axis=1)
        sv = self.s[i] + t
        rr = np.interp(sv, self.s, self.r)
        d = rho - rr
        if scales:
            phi = np.arctan2(np.einsum("ij,ij->i", w, self.B[i]), np.einsum("ij,ij->i", w, self.N[i]))
            d -= scale_relief(np.interp(sv, self.s, self.u), phi, rr, sv)
        out[hit] = d
        return out


def scale_relief(u, phi, r, s):
    """Height of the scale coat. Dorsal scales are diamonds on a staggered lattice, each
    rising toward its free (tailward) edge, which is rounded; the edge it tucks under is a
    sharp step, so the coat reads as overlapping. The belly band carries transverse plates."""
    v = phi / (2 * math.pi) * ROWS
    best = None
    for j in (np.floor(v), np.floor(v) + 1):
        cu = np.round(u - 0.5 * j) + 0.5 * j
        du, dv = u - cu, v - j
        qd = 2 * np.abs(du) + np.abs(dv)
        if best is None:
            best = [qd, du, dv]
        else:
            take = qd < best[0]
            for k, val in enumerate((qd, du, dv)):
                best[k] = np.where(take, val, best[k])
    qd, du, dv = best
    eps = 0.07 + 0.2 * smoothstep(-0.2, 0.25, du)
    dorsal = (0.3 + 0.7 * (du + 0.5)) * rounded((1 - qd) / eps) * (1 - 0.3 * dv * dv)
    tv = u - np.floor(u)
    ventral = 0.8 * (0.3 + 0.7 * tv) * rounded((1 - tv) / 0.15) * rounded(tv / 0.05)
    a = np.abs(phi)
    h = np.where(a < BELLY, ventral, dorsal) * rounded(np.abs(a - BELLY) / SEAM)
    return SCALE_H * r * smoothstep(8.0, 13.0, s) * h


class Piece:
    """The print. Without the egg the serpent keeps exactly its coiled shape, as if wound round
    an egg that is not there; the egg still sets the geometry, it just isn't part of the solid."""

    def __init__(self, with_egg=True):
        self.with_egg = with_egg
        self.egg = Egg()
        self.head = Head((HEAD_AT[0], HEAD_AT[1], self.egg.top + HEAD_LIFT), HEAD_DIR)
        self.body = Serpent(self.egg, self.head)
        self.key = Key(self.head)

    def field(self, P, fine=True):
        d = self.body.sdf(P, fine)
        if self.with_egg:
            d = smin(self.egg.sdf(P), d, K_EGG + K_FOOT * smoothstep(9.0, 2.0, P[:, 2] - Z_G))
        d = smin(d, self.head.sdf(P), K_HEAD)
        d = smin(d, self.key.sdf(P), K_KEY)                     # the lips close round the shaft
        return np.maximum(d, Z_G - P[:, 2])

    def bounds(self):
        b = self.body
        lo = np.minimum(b.C.min(axis=0) - BODY_R, [-EGG_W / 2, -EGG_W / 2, Z_G])
        hi = np.maximum(b.C.max(axis=0) + BODY_R, [EGG_W / 2, EGG_W / 2, self.egg.top])
        corners = np.array([[f, l, u] for f in self.head.BOX[0] for l in self.head.BOX[1]
                            for u in self.head.BOX[2]])
        hc = np.r_[self.head.O + HEAD_SCALE * corners @ self.head.M, self.key.corners()]
        lo, hi = np.minimum(lo, hc.min(axis=0)), np.maximum(hi, hc.max(axis=0))
        return lo - 3.0, hi + 3.0


def lerp_matrix(fine, coarse):
    h = coarse[1] - coarse[0]
    i = np.clip(((fine - coarse[0]) / h).astype(int), 0, len(coarse) - 2)
    t = (fine - coarse[i]) / h
    W = np.zeros((len(fine), len(coarse)))
    W[np.arange(len(fine)), i] = 1 - t
    W[np.arange(len(fine)), i + 1] = t
    return W


def mesh(piece, h=H, hc=HC, band=BAND, slab=48):
    lo, hi = piece.bounds()
    cg = [np.arange(lo[k] - hc, hi[k] + 2 * hc, hc) for k in range(3)]
    G = np.stack(np.meshgrid(*cg, indexing="ij"), axis=-1).reshape(-1, 3)
    coarse = np.concatenate([piece.field(c, fine=False) for c in np.array_split(G, 8)])
    coarse = coarse.reshape([len(c) for c in cg])
    fg = [lo[k] + h * (np.arange(int((hi[k] - lo[k]) / h) + 1) + 0.37) for k in range(3)]
    Wx, Wy, Wz = (lerp_matrix(fg[k], cg[k]) for k in range(3))
    nz = len(fg[2])
    verts, faces, nv, nexact = [], [], 0, 0
    for k0 in range(0, nz - 1, slab):
        ks = np.arange(k0, min(k0 + slab, nz - 1) + 1)
        vol = np.tensordot(coarse, Wz[ks], axes=([2], [1]))
        vol = np.tensordot(Wy, vol, axes=([1], [1]))
        vol = np.tensordot(Wx, vol, axes=([1], [1]))
        idx = np.nonzero(np.abs(vol) < band)
        if len(idx[0]):
            pts = np.c_[fg[0][idx[0]], fg[1][idx[1]], fg[2][ks[idx[2]]]]
            vol[idx] = np.concatenate([piece.field(c) for c in
                                       np.array_split(pts, max(1, len(pts) // 1_500_000))])
            nexact += len(pts)
        vol[vol == 0] = 1e-7
        if vol.min() > 0 or vol.max() < 0:
            continue
        v, f, _, _ = marching_cubes(vol, 0.0, gradient_direction="descent")
        v[:, 2] += k0
        verts.append(v)
        faces.append(f + nv)
        nv += len(v)
    V, F = np.concatenate(verts), np.concatenate(faces)
    key, inv = np.unique(np.round(V, 6), axis=0, return_inverse=True)   # weld slab seams
    F = inv.reshape(-1)[F]
    F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 2] != F[:, 0])]
    V = lo + h * (key + 0.37)
    m = trimesh.Trimesh(V, F, process=False)
    parts = sorted(m.split(only_watertight=False), key=lambda p: -len(p.faces))
    crumb = max((p.extents.max() for p in parts[1:]), default=0.0)
    print(f"    field: {nexact / 1e6:.1f} M exact samples at {h} mm, {len(F) / 1e6:.2f} M triangles, "
          f"{len(parts) - 1} sub-voxel crumbs dropped (largest {crumb:.2f} mm)")
    m = parts[0]
    bed = np.abs(m.vertices[:, 2] - Z_G) < 1e-4
    filter_taubin(m, lamb=0.5, nu=0.53, iterations=TAUBIN)      # volume-preserving (Taubin 1995)
    m.vertices[:, 2] = np.where(bed, Z_G, np.maximum(m.vertices[:, 2], Z_G))
    return m


def decimate(m, faces):
    import pymeshlab
    if len(m.faces) <= faces:
        return m
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(vertex_matrix=m.vertices, face_matrix=m.faces))
    ms.meshing_decimation_quadric_edge_collapse(targetfacenum=faces, qualitythr=0.5,
                                                preservetopology=True, preservenormal=True,
                                                optimalplacement=True, planarquadric=True)
    c = ms.current_mesh()
    return trimesh.Trimesh(c.vertex_matrix(), c.face_matrix(), process=False)


def carve_pupils(piece, m, h=0.04, half=3.6):
    """Cut both vortex pupils out of the finished mesh with manifold3d. Each cutter is its own
    marching-cubes solid at h, fine enough for the 0.45 mm spiral that the body's voxel would
    alias away."""
    import manifold3d as m3
    hd = piece.head
    man = m3.Manifold(m3.Mesh(vert_properties=np.asarray(m.vertices, np.float32),
                              tri_verts=np.asarray(m.faces, np.uint32)))
    for side in (1.0, -1.0):
        c = hd.O + HEAD_SCALE * (hd.eye_at * [1.0, side, 1.0]) @ hd.M
        ax = c[:, None] + np.arange(-half, half + h / 2, h)
        n = ax.shape[1]
        d = hd.pupil(np.stack(np.meshgrid(*ax, indexing="ij"), axis=-1).reshape(-1, 3)).reshape(n, n, n)
        d[[0, -1]], d[:, [0, -1]], d[:, :, [0, -1]] = 1.0, 1.0, 1.0
        v, f, _, _ = marching_cubes(d, 0.0, spacing=(h, h, h), gradient_direction="descent")
        cut = trimesh.Trimesh(v + ax[:, 0], f, process=False)
        cut = max(cut.split(only_watertight=False), key=lambda p: len(p.faces))   # no specks
        cut = m3.Manifold(m3.Mesh(vert_properties=np.asarray(cut.vertices, np.float32),
                                  tri_verts=np.asarray(cut.faces, np.uint32)))
        assert cut.status() == m3.Error.NoError, cut.status()
        man = man - cut
    m = sg.to_trimesh(man)
    m.merge_vertices(digits_vertex=5)
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.remove_unreferenced_vertices()
    return m


def stability(m):
    """Tip angle: how far the piece can lean before its centre of mass leaves the footprint."""
    from scipy.spatial import ConvexHull
    from shapely.geometry import Point, Polygon
    flat = (m.vertices[m.faces][:, :, 2] < Z_G + 1e-3).all(axis=1)
    foot = m.vertices[m.vertices[:, 2] < Z_G + 1e-3][:, :2]
    hull = Polygon(foot[ConvexHull(foot).vertices])
    c = m.center_mass
    margin = hull.exterior.distance(Point(c[:2])) if hull.contains(Point(c[:2])) else -1.0
    return math.degrees(math.atan2(margin, c[2] - Z_G)), hull.area, m.area_faces[flat].sum(), c


def main():
    t0 = time.time()
    piece = Piece(with_egg=not SERPENT)
    b = piece.body
    print(f"    egg {EGG_H:.0f} x {EGG_W:.0f} mm, serpent {b.s[-1]:.0f} mm long, "
          f"lands after {b.land_turns:.2f} turns, min tube-egg gap {b.egg_gap.min():+.2f} mm")
    gap, s1, s2 = b.clearance()
    print(f"    closest approach of the body to itself {gap:.1f} mm (s = {s1:.0f} / {s2:.0f} mm)")
    m = mesh(piece)
    print(f"    meshed in {time.time() - t0:.0f} s")
    m = decimate(m, TARGET_FACES)
    V = m.vertices.astype(np.float32).astype(float)
    V[:, 2] = np.maximum(V[:, 2], Z_G)                 # decimation nudges a few microns below the bed
    m = trimesh.Trimesh(V, m.faces)                    # trimesh's own merge + degenerate drop
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    m = carve_pupils(piece, m)
    name = "orphic_serpent" if SERPENT else "orphic_egg"
    path = f"{name}_draft.stl" if DRAFT else f"{name}.stl"
    m.export(path)
    back = trimesh.load(path)
    ok = sg.check(back, path)
    tip, area, contact, c = stability(back)
    print(f"    on the bed {contact:.0f} mm2 within a {area:.0f} mm2 footprint, "
          f"centre of mass {c.round(1)}, tips at {tip:.1f} deg")
    print(f"    done in {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
