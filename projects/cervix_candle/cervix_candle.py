#!/usr/bin/env python3
"""The uterus as a candle holder, modelled on the user's reference photo (reference_photo.png)
of a white sculpture after the coronal-section diagram (reference_diagram.png). The 2D layout
is traced off the photo in photo_trace.py.

The front is flat, like a section through the organ, with softly rounded edges. The back is
domed: the outline is inflated as the union of its inscribed spheres, flattened to
AP_SCALE, so every part is round in section. The flat-topped tube tabs, the fundus band, the
walls sweeping down to the candle, the two lips with the canal between them and the cup round
their ends all follow the photo. The cavity is hollowed into the face like the inside of a
vessel cut in half: the same inflation, so it is deepest where it is widest.

The candle (0.23 in, 5.84 mm) stands down in the cervical canal, as in the photo. The canal is
the photo's S-curved slit between the lips: 4.2 mm where they close below the hollow, swelling
to 6.4 mm, pinching to 3.4 mm at the lips' knees and opening round their rounded ends into the
cup. Behind it runs a round channel, 6 mm behind the face. It is an hourglass: the channel
narrows smoothly to a 5.5 mm ring behind the canal's upper pinch, which grips the candle gently
as it is pressed down through, then opens to 6.2 mm below the ring. It ends in a rounded floor
that the candle's foot rests on. Both pinches are narrower than the candle, so it is held
at the ring and at its foot and can't come out forward. The lumen of each tube is a groove
along the fundus band's lower edge, running from the corner of the cavity out through the
tab's end.

Every cut (hollow, canal, channel, cup, grooves) is joined to its neighbours with a smooth
minimum and rounded into the face with a smooth maximum, so there are no ledges or sharp
inside corners anywhere.

The model is written in two poses:
    cervix_candle.stl          flat face down for printing (support only inside the hollow)
    cervix_candle_upright.stl  standing as displayed; dish.py takes the socket from this one

Run from this directory:  python cervix_candle.py
"""

import numpy as np
import trimesh
from scipy import ndimage as nd
from scipy.interpolate import PchipInterpolator
from skimage.measure import marching_cubes

import photo_trace as pt

HEIGHT = 100.0          # mm, base to the top of the fundus
AP_SCALE = 0.85         # back dome depth relative to a round section
CUT_Y = -2.0            # the flat front face (y = 0 is the plane the dome and hollow inflate about)
EDGE_R = 2.5            # roundover where the face meets the domed back
LIP_R = 3.0             # roundover on every cut's rim: the photo's bullnose
INNER_K = 1.2           # fillet where cuts meet each other
CANDLE_DEPTH = 6.0      # candle axis behind the face: deep enough that the lips close in front
CANDLE_D = 5.84         # the user's candles: 0.23 in
CHANNEL_R = 3.1        # candle channel radius down the canal (6.2 mm, the candle slides)
RING_R = 2.75           # ... narrowing at the internal orifice (5.5 mm): a gentle press fit
FLARE_R = 3.6           # ... and opening out above it into the hollow
VENT_R = 1.5            # radius of the tube lumen grooves
VOX = 0.2
FLAME = 25.0            # wick base to flame tip, generous


def sdf2(mask, sigma):
    """Signed distance in pixels (negative inside), smoothed to remove the pixel staircase."""
    d = nd.distance_transform_edt(~mask) - nd.distance_transform_edt(mask)
    return nd.gaussian_filter(d, sigma)


def polygon_sdf(poly, X, Z):
    """Exact signed distance (negative inside) from grid points to a closed polygon (x, z)."""
    d = np.full(X.shape, np.inf)
    inside = np.zeros(X.shape, bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        for (x0, z0), (x1, z1) in zip(poly, np.roll(poly, -1, axis=0)):
            dx, dz = x1 - x0, z1 - z0
            t = np.clip(((X - x0) * dx + (Z - z0) * dz) / (dx * dx + dz * dz), 0, 1)
            d = np.minimum(d, np.hypot(X - x0 - t * dx, Z - z0 - t * dz))
            inside ^= ((z0 > Z) != (z1 > Z)) & (X < x0 + (Z - z0) * dx / dz)
    return np.where(inside, -d, d)


def _envelope(f):
    """1D lower envelope of parabolas: g(q) = min_p((q - p)^2 + f(p)) (Felzenszwalb)."""
    n = len(f)
    v = np.zeros(n, int)
    zb = np.empty(n + 1)
    k, zb[0], zb[1] = 0, -np.inf, np.inf
    for q in range(1, n):
        while True:
            p = v[k]
            sx = ((f[q] + q * q) - (f[p] + p * p)) / (2 * q - 2 * p)
            if sx > zb[k]:
                break
            k -= 1
        k += 1
        v[k], zb[k], zb[k + 1] = q, sx, np.inf
    g = np.empty(n)
    k = 0
    for q in range(n):
        while zb[k + 1] < q:
            k += 1
        g[q] = (q - v[k]) ** 2 + f[v[k]]
    return g


def inflate(mask):
    """Half-depth (px) of the union of all balls inscribed in the mask, centred on its plane:
    h(q)^2 = max_p(D(p)^2 - |q - p|^2) over interior pixels p, D = distance to the edge. That is
    a squared distance transform of -D^2, done separably, so no skeleton can drop a branch."""
    D = nd.distance_transform_edt(mask)
    f = np.where(mask, -D * D, 1e9)
    g = np.apply_along_axis(_envelope, 0, f)
    g = np.apply_along_axis(_envelope, 1, g)
    return nd.gaussian_filter(np.sqrt(np.maximum(-g, 0)), 2.0)   # soften the pixel lattice


def smax(a, b, k):
    h = np.clip(0.5 - 0.5 * (b - a) / k, 0, 1)
    return b + (a - b) * h + k * h * (1 - h)


def smin(a, b, k):
    return -smax(-a, -b, k)


class Uterus:
    def __init__(self):
        sil, cav, grooves = pt.masks()
        ys, _ = np.nonzero(sil)
        self.s = HEIGHT / (pt.BASE - ys.min())                     # mm per photo pixel
        self.sil_px = sdf2(sil, 2.0)
        self.cav_px = sdf2(cav, 1.5)
        o = pt.canal_outline()
        self.canal = np.c_[(o[:, 0] - pt.AXIS) * self.s, self.z(o[:, 1])]
        # The base is a cut through the cervix, not its end: extend the neck straight down
        # before inflating, so the dome keeps its full depth all the way to the base.
        ext = np.vstack([sil, np.repeat(sil[-1:], 80, axis=0)])
        self.sil_h = inflate(ext)[: sil.shape[0]]
        self.cav_h = inflate(cav)
        self.sil, self.cav = sil, cav
        self.grooves = [np.c_[(g[:, 0] - pt.AXIS) * self.s, self.z(g[:, 1])] for g in grooves]
        self.z_ring = self.z(pt.ORIFICE)                           # the squeeze
        self.z_wick = self.z(pt.CANDLE_TOP)
        self.z_floor = self.z(pt.CANAL_FLOOR)                      # lowest point of the floor
        self.z_top = self.z_wick + 4.0                             # channel's rounded top end
        self.y_axis = CUT_Y + CANDLE_DEPTH
        self.half_w = (np.nonzero(sil.any(axis=0))[0] - pt.AXIS).__abs__().max() * self.s
        zr = self.z_ring
        self.channel_r = PchipInterpolator(
            [self.z_floor - 5, zr - 6, zr, zr + 5, self.z_top + 5],
            [CHANNEL_R, CHANNEL_R, RING_R, FLARE_R, FLARE_R])

    def z(self, py):
        return (pt.BASE - np.asarray(py, float)) * self.s

    def sample(self, field, x, z):
        r, c = pt.BASE - z / self.s, pt.AXIS + x / self.s
        return nd.map_coordinates(field, [r, c], order=1, mode="nearest") * self.s

    def floor_contact(self, d=CANDLE_D):
        """Height of a flat-footed candle of diameter d resting in the channel's rounded floor."""
        r = float(self.channel_r(self.z_floor))
        return self.z_floor + r - np.sqrt(max(r * r - (d / 2) ** 2, 0.0))

    def field(self):
        y_back = AP_SCALE * self.s * self.sil_h.max()
        off = 0.37 * VOX                                          # keep flat faces off grid nodes,
        xs = np.arange(-self.half_w - 2, self.half_w + 2 + VOX, VOX) + off   # or MC leaves holes
        ys = np.arange(CUT_Y - 2, y_back + 2 + VOX, VOX) + off
        zs = np.arange(-2, HEIGHT + 2 + VOX, VOX) + off
        X, Z = np.meshgrid(xs, zs, indexing="ij")
        d_sil = np.maximum(self.sample(self.sil_px, X, Z), -Z)    # flat base
        d_cav = self.sample(self.cav_px, X, Z)
        d_canal = polygon_sdf(self.canal, X, Z)
        h_sil = self.sample(self.sil_h, X, Z) * AP_SCALE
        h_cav = self.sample(self.cav_h, X, Z)

        groove = np.full(X.shape, np.inf)                         # distance to lumen centreline (x, z)
        for path in self.grooves:
            for (x0, z0), (x1, z1) in zip(path[:-1], path[1:]):
                dx, dz = x1 - x0, z1 - z0
                t = np.clip(((X - x0) * dx + (Z - z0) * dz) / (dx * dx + dz * dz), 0, 1)
                groove = np.minimum(groove, np.hypot(X - x0 - t * dx, Z - z0 - t * dz))

        # Candle channel: a capsule along the axis with an hourglass radius, rounded at both ends.
        zcl = np.clip(Z, self.z_floor + CHANNEL_R, self.z_top)
        r_ch = self.channel_r(zcl)
        dz_end = Z - zcl
        # Canal gap between the lips, in the traced outline: from the face back to the candle
        # axis (its back corners rounded, so where it swells wider than the channel it still
        # meets it smoothly), from the base up into the hollow (so nothing bridges in front of
        # the candle), easing below the floor along a long smoothstep ramp to a shallow groove
        # that runs on into the cup.
        ya = self.y_axis
        t = np.clip((self.z_floor - Z) / 8, 0, 1)
        back = ya + (CUT_Y + 2.5 - ya) * t * t * (3 - 2 * t)
        canal_z = np.maximum(Z - self.z_top, -1 - Z)              # z-extent, from the base up

        F = np.empty((len(xs), len(ys), len(zs)), np.float32)
        for j, y in enumerate(ys):
            dome = np.maximum(abs(y) - h_sil, d_sil)
            outer = smax(dome, CUT_Y - y, EDGE_R)                 # flat face, rounded edge
            hollow = np.maximum(abs(y) - h_cav, d_cav)
            lumen = np.hypot(groove, y - CUT_Y) - VENT_R          # half-round groove in the face
            channel = np.sqrt(X * X + (y - ya) ** 2 + dz_end ** 2) - r_ch
            canal = np.maximum(smax(d_canal, y - back, INNER_K), canal_z)
            cut = smin(smin(hollow, lumen, INNER_K), smin(channel, canal, INNER_K), INNER_K)
            F[:, j, :] = smax(outer, -cut, LIP_R)
        return F, (xs[0], ys[0], zs[0])

    def front_opening(self, F, origin, z):
        """Narrowest gap between the lips in front of the candle axis at height z (mm): the
        candle can only come out forward if this is wider than it."""
        x0, y0, z0 = origin
        k = int(round((z - z0) / VOX))
        widths = []
        for y in np.arange(CUT_Y - 0.5, self.y_axis, VOX):
            col = F[:, int(round((y - y0) / VOX)), k]
            xs = x0 + VOX * np.arange(len(col))
            # contiguous empty run through x = 0
            i0 = int(np.argmin(np.abs(xs)))
            if col[i0] <= 0:
                return 0.0                                       # solid in front of the axis
            lo = hi = i0
            while lo > 0 and col[lo - 1] > 0:
                lo -= 1
            while hi < len(col) - 1 and col[hi + 1] > 0:
                hi += 1
            widths.append(xs[hi] - xs[lo] + VOX)
        return min(widths)

    def report(self):
        zs = np.arange(self.z_wick, HEIGHT, 0.1)
        inside = self.sample(self.cav_px, np.zeros_like(zs), zs) < 0
        ceiling = zs[np.nonzero(inside)[0].max()]
        back = AP_SCALE * self.s * self.sil_h - self.s * self.cav_h
        foot = self.floor_contact()
        return {"mm_per_px": self.s, "width": 2 * self.half_w,
                "depth": AP_SCALE * self.s * self.sil_h.max() - CUT_Y,
                "ring_z": self.z_ring, "ring_dia": 2 * RING_R,
                "candle_foot_z": foot, "wick_z": self.z_wick,
                "candle_trim_len": self.z_wick - foot, "cavity_ceiling_z": ceiling,
                "flame_to_ceiling": ceiling - self.z_wick - FLAME,
                "back_wall_min": back[self.cav].min()}


def build():
    u = Uterus()
    F, origin = u.field()
    gaps = {z: u.front_opening(F, origin, z) for z in np.arange(u.z_floor, u.z_ring + 5, 1.0)}
    v, f, _, _ = marching_cubes(F, 0.0, spacing=(VOX,) * 3, gradient_direction="descent")
    m = trimesh.Trimesh(v + origin, f, process=True)
    m = max(m.split(only_watertight=False), key=lambda p: len(p.faces))
    trimesh.smoothing.filter_taubin(m, iterations=8)
    m.vertices[:, 2] = np.maximum(m.vertices[:, 2], 0.0)          # flat base exactly on z = 0
    m.export("cervix_candle_upright.stl")
    lying = m.copy()                                             # flat face down
    lying.vertices = np.c_[m.vertices[:, 0], -m.vertices[:, 2], m.vertices[:, 1] - CUT_Y]
    lying.vertices[:, 2] = np.maximum(lying.vertices[:, 2], 0.0)
    lying.export("cervix_candle.stl")
    return u, m, lying, gaps


if __name__ == "__main__":
    u, m, lying, gaps = build()
    for k, val in u.report().items():
        print(f"{k:24s} {val:7.2f}")
    print("front opening by height:", "  ".join(f"z{z:.0f}:{w:.1f}" for z, w in gaps.items()))
    for name, mesh in (("upright", m), ("print", lying)):
        print(name, "extents", np.round(mesh.extents, 1), "watertight", mesh.is_watertight,
              "volume cm3", round(mesh.volume / 1000, 1), "faces", len(mesh.faces))
