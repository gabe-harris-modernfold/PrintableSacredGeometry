#!/usr/bin/env python3
"""Mesh kit for the three voice-horn POCs (docs/voice-horn-pocs.md).

The one idea that carries this file: a printable shell is two radius FIELDS over
the same (z, theta) grid -- inner and outer -- so bore profile, helical threads,
QRS wells and rib arrays are all additive terms on a field, and the resulting
quad shell comes out watertight by construction (no booleans, no repair pass).

Ornament (thorn/dome fields, braille) is emitted as disjoint solids and unioned
once through manifold3d; lattices (pop cage, wind cowls, Fleece coupon) are strut
soups unioned the same way.
"""
import math
import os

import numpy as np
import trimesh

C_AIR = 343000.0                      # mm/s, speech-band sound speed
GOLDEN = math.radians(137.50776405)   # phyllotaxis divergence angle
WALL = 5.0                            # F1 sandwich total wall (2 shells + 40% Grid)
BED = 320.0                           # printer envelope, mm


# ----------------------------------------------------------------- shells

def clean(m):
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    return m


def shell(z, ri, ro, cap_start=True, cap_end=True, theta=None):
    """Watertight tube-with-wall. z:(nz,), ri/ro:(nz,nt) radii in mm.

    Inner-surface normals point into the bore, outer normals point out, and the
    two end annuli close the solid. `theta` may be non-uniform -- that is how
    the QRS dividers and the Halo ribs get crisp edges without paying for a fine
    grid everywhere."""
    z = np.asarray(z, float)
    ri, ro = np.asarray(ri, float), np.asarray(ro, float)
    nz, nt = ri.shape
    th = grid_theta(nt) if theta is None else np.asarray(theta, float)
    ct, st = np.cos(th)[None, :], np.sin(th)[None, :]
    zc = z[:, None] * np.ones((1, nt))
    Vi = np.stack([ri * ct, ri * st, zc], -1)
    Vo = np.stack([ro * ct, ro * st, zc], -1)
    V = np.concatenate([Vi, Vo], 1).reshape(-1, 3)          # station-major
    idx = lambda k, j, out: k * 2 * nt + (nt if out else 0) + (j % nt)

    k = np.arange(nz - 1)[:, None]
    j = np.arange(nt)[None, :]
    j1 = (j + 1) % nt
    F = []
    o00, o01 = idx(k, j, 1), idx(k, j1, 1)
    o10, o11 = idx(k + 1, j, 1), idx(k + 1, j1, 1)
    F += [np.stack([o00, o01, o10], -1), np.stack([o01, o11, o10], -1)]
    i00, i01 = idx(k, j, 0), idx(k, j1, 0)
    i10, i11 = idx(k + 1, j, 0), idx(k + 1, j1, 0)
    F += [np.stack([i00, i10, i01], -1), np.stack([i01, i10, i11], -1)]
    if cap_start:
        a, b = idx(0, j, 0), idx(0, j1, 0)
        c, d = idx(0, j, 1), idx(0, j1, 1)
        F += [np.stack([a, b, d], -1), np.stack([a, d, c], -1)]
    if cap_end:
        e = nz - 1
        a, b = idx(e, j, 0), idx(e, j1, 0)
        c, d = idx(e, j, 1), idx(e, j1, 1)
        F += [np.stack([a, d, b], -1), np.stack([a, c, d], -1)]
    F = np.concatenate([f.reshape(-1, 3) for f in F])
    return clean(trimesh.Trimesh(V, F, process=False))


def loft(sections, cap_start=True, cap_end=True):
    """Loft equal-length closed 3-D point rings into a tube. sections:(ns,nt,3).

    The Volute's duct is a swept vaulted section on a spiral path -- no axis of
    revolution to hang a radius field on -- so it is built this way instead."""
    P = np.asarray(sections, float)
    ns, nt, _ = P.shape
    V = P.reshape(-1, 3)
    k = np.arange(ns - 1)[:, None]
    j = np.arange(nt)[None, :]
    j1 = (j + 1) % nt
    a, b = k * nt + j, k * nt + j1
    c, d = (k + 1) * nt + j1, (k + 1) * nt + j
    F = [np.stack([a, b, c], -1).reshape(-1, 3),
         np.stack([a, c, d], -1).reshape(-1, 3)]
    extra = []
    if cap_start:
        f, v = fan(np.arange(nt), P[0], len(V) + len(extra), flip=True)
        F.append(f); extra.append(v)
    if cap_end:
        f, v = fan(np.arange(nt) + (ns - 1) * nt, P[-1], len(V) + len(extra),
                   flip=False)
        F.append(f); extra.append(v)
    if extra:
        V = np.vstack([V, np.asarray(extra)])
    return clean(trimesh.Trimesh(V, np.concatenate(F), process=False))


def fan(ring, ring_pts, hub_index, flip=False):
    """Fan a closed ring around a new centroid vertex.

    Fanning from a ring vertex instead would emit zero-area triangles wherever
    the outline samples a straight edge -- and those get culled, leaving holes."""
    n = len(ring)
    j = np.arange(n)
    tri = np.stack([np.full(n, hub_index), ring[j], ring[(j + 1) % n]], -1)
    return (tri[:, ::-1].copy() if flip else tri), np.asarray(ring_pts).mean(0)


def grid_theta(nt):
    return np.linspace(0.0, 2 * math.pi, nt, endpoint=False)


# ----------------------------------------------------------------- threads

def thread(z, th, z0, z1, pitch, depth, lead=4.0):
    """Trapezoidal single-lead helical ridge as an additive radius field (mm).

    Male and female halves call this with identical arguments, so the pair mates
    by construction; the female adds its clearance to the base radius instead of
    altering the profile."""
    Z, T = np.meshgrid(np.asarray(z, float), np.asarray(th, float), indexing="ij")
    phase = np.mod(Z - pitch * T / (2 * math.pi), pitch) / pitch
    ridge = np.clip(np.minimum(phase / 0.35, (1.0 - phase) / 0.35), 0.0, 1.0)
    win = np.clip(np.minimum((Z - z0) / lead, (z1 - Z) / lead), 0.0, 1.0)
    return depth * ridge * win


# --------------------------------------------------------- ornament fields

def phyllotaxis(prof_z, prof_r, n, weight=None, z_lo=None, z_hi=None):
    """Golden-angle sites on a surface of revolution, area-uniform x weight.

    Returns (points, normals, frac, area_per_site): frac is the site's
    fractional distance along the sampled meridian -- the argument RULE S2
    grades spike height on."""
    z, r = np.asarray(prof_z, float), np.asarray(prof_r, float)
    step = np.hypot(np.diff(z), np.diff(r))       # drop repeated stations first:
    keep = np.concatenate([[True], step > 1e-9])  # a zero-length one makes the
    z, r = z[keep], r[keep]                       # arc-length gradient NaN
    ds = np.hypot(np.diff(z), np.diff(r))
    zm, rm = 0.5 * (z[:-1] + z[1:]), 0.5 * (r[:-1] + r[1:])
    w = np.ones_like(zm) if weight is None else np.asarray(weight(zm), float)
    if z_lo is not None:
        w = np.where(zm >= z_lo, w, 0.0)
    if z_hi is not None:
        w = np.where(zm <= z_hi, w, 0.0)
    dA = 2 * math.pi * rm * ds * w
    cum = np.concatenate([[0.0], np.cumsum(dA)])
    if cum[-1] <= 0:
        raise ValueError("phyllotaxis: empty weighted band")
    s_edge = np.concatenate([[0.0], np.cumsum(ds)])
    s_hit = np.interp((np.arange(n) + 0.5) / n * cum[-1], cum, s_edge)
    zk, rk = np.interp(s_hit, s_edge, z), np.interp(s_hit, s_edge, r)
    dz = np.interp(s_hit, s_edge, np.gradient(z, s_edge))
    dr = np.interp(s_hit, s_edge, np.gradient(r, s_edge))
    thk = np.arange(n) * GOLDEN
    P = np.stack([rk * np.cos(thk), rk * np.sin(thk), zk], -1)
    nrm = np.stack([dz * np.cos(thk), dz * np.sin(thk), -dr], -1)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    frac = (s_hit - s_hit.min()) / max(np.ptp(s_hit), 1e-9)
    return P, nrm, frac, cum[-1] / n


def frusta(P, axis, h, r_base, r_tip=0.6, nseg=8, sink=1.5):
    """Disjoint thorn/dome solids (truncated cones), ready for one union."""
    P = np.asarray(P, float)
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    h = np.broadcast_to(np.asarray(h, float), (len(P),))
    r_base = np.broadcast_to(np.asarray(r_base, float), (len(P),))
    r_tip = np.broadcast_to(np.asarray(r_tip, float), (len(P),))
    sink = np.broadcast_to(np.asarray(sink, float), (len(P),))
    up = np.array([0.0, 0.0, 1.0])
    ph = grid_theta(nseg)
    V, F, off = [], [], 0
    for i in range(len(P)):
        a = axis[i]
        t = up if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u = np.cross(a, t); u /= np.linalg.norm(u)
        w = np.cross(a, u)
        ring = np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * w
        b, tp = P[i] - a * sink[i], P[i] + a * h[i]
        V.append(np.concatenate([b + r_base[i] * ring, tp + r_tip[i] * ring,
                                 [b], [tp]]))
        j = np.arange(nseg); j1 = (j + 1) % nseg
        F.append(np.concatenate([
            np.stack([off + j, off + j1, off + nseg + j1], -1),
            np.stack([off + j, off + nseg + j1, off + nseg + j], -1),
            np.stack([np.full(nseg, off + 2 * nseg), off + j1, off + j], -1),
            np.stack([np.full(nseg, off + 2 * nseg + 1), off + nseg + j,
                      off + nseg + j1], -1)]))
        off += 2 * nseg + 2
    return trimesh.Trimesh(np.concatenate(V), np.concatenate(F), process=False)


# --------------------------------------------------------- texture fields
#
# High-relief exterior skin (RULE S3: the exterior is acoustically free).
# These are additive (z, theta) fields like thread(): outward-only relief keeps
# the shell watertight by construction and keeps every unioned ornament solid
# rooted in the wall. Inward features (finger pits) are carved out of a bulge
# that was added first, so the 5 mm wall law is never thinned.

def smooth_band(z, z0, z1, ramp=6.0):
    """C1 window: 0 outside [z0, z1], 1 inside, cosine ramps at both ends."""
    z = np.asarray(z, float)
    t = np.clip(np.minimum(z - z0, z1 - z) / max(ramp, 1e-9), 0.0, 1.0)
    return 0.5 - 0.5 * np.cos(math.pi * t)


def _fade(t):
    return t * t * t * (t * (6.0 * t - 15.0) + 10.0)


def _vnoise(zc, tc, nt_cells, rng, ridged=True):
    """One octave of value noise on a lattice, periodic in the t direction.

    zc: cell coords (any shape), tc: cell coords with period nt_cells."""
    nz_cells = int(np.floor(zc.max())) + 2
    lat = rng.random((nz_cells + 1, nt_cells))
    i0 = np.clip(np.floor(zc).astype(int), 0, nz_cells - 1)
    j0 = np.floor(tc).astype(int) % nt_cells
    fz, ft = _fade(zc - np.floor(zc)), _fade(tc - np.floor(tc))
    j1 = (j0 + 1) % nt_cells
    v = ((lat[i0, j0] * (1 - fz) + lat[i0 + 1, j0] * fz) * (1 - ft)
         + (lat[i0, j1] * (1 - fz) + lat[i0 + 1, j1] * fz) * ft)
    return 1.0 - np.abs(2.0 * v - 1.0) if ridged else v


def ridge_noise(z, th, r_ref, z0, z1, amp_lo, amp_hi, cell_z=25.0, cell_arc=9.0,
                seed=0, octaves=3, ramp=6.0, sharp=1.6, grow=1.0):
    """Ridged-multifractal relief on a (z, theta) grid, periodic in theta.

    Crests are C0 creases (the |...| fold), valleys stay smooth -- jagged to the
    hand along the ridge lines, pleasant in between. Amplitude ramps amp_lo ->
    amp_hi across [z0, z1] (RULE S2: relief grows with the flare), windowed to
    zero outside it. r_ref converts theta to arc mm so cells are isotropic-ish."""
    z = np.asarray(z, float)
    th = np.asarray(th, float)
    rng = np.random.default_rng(seed)
    nt_cells = max(4, int(round(2 * math.pi * r_ref / cell_arc)))
    Z, T = np.meshgrid((z - z0) / cell_z, th / (2 * math.pi), indexing="ij")
    acc, wsum, w = np.zeros(Z.shape), 0.0, 1.0
    for o in range(octaves):
        k = 2 ** o
        acc += w * _vnoise(Z * k - Z.min() * k + 1.0, T * nt_cells * k,
                           nt_cells * k, rng)
        wsum += w
        w *= 0.5
    field = (acc / wsum) ** sharp
    x = np.clip((z - z0) / max(z1 - z0, 1e-9), 0.0, 1.0)
    amp = amp_lo + (amp_hi - amp_lo) * x ** grow
    return field * (amp * smooth_band(z, z0, z1, ramp))[:, None]


def voronoi_shards(z, th, r_of_z, z0, z1, spacing=22.0, ridge_w=3.2, amp=2.5,
                   seed=0, ramp=6.0, mask=None):
    """Crystalline plate field: flat Voronoi cells, raised jagged inter-plate
    ridges (height ~ closeness to the cell boundary, F2 - F1). Periodic in
    theta; distances use the local radius so the pattern stays isotropic on a
    cone. `mask(z, th) -> (nz, nt)` multiplies the field (keep-out zones)."""
    z = np.asarray(z, float)
    th = np.asarray(th, float)
    r = np.asarray(r_of_z, float)
    rng = np.random.default_rng(seed)
    n = max(8, int(round((z1 - z0) * 2 * math.pi * float(r.mean()) / spacing ** 2)))
    sz = z0 + (z1 - z0) * rng.random(n)
    st = 2 * math.pi * rng.random(n)
    out = np.zeros((len(z), len(th)))
    for a in range(0, len(z), 24):                # chunk rows: nz*nt*n floats
        b = min(a + 24, len(z))
        dz = z[a:b, None, None] - sz[None, None, :]
        dth = (th[None, :, None] - st[None, None, :] + math.pi) \
            % (2 * math.pi) - math.pi
        da = r[a:b, None, None] * dth
        d = np.sqrt(dz * dz + da * da)
        d.partition(1, axis=2)
        out[a:b] = np.clip(1.0 - (d[:, :, 1] - d[:, :, 0]) / ridge_w, 0.0, 1.0)
    field = amp * out ** 0.8 * smooth_band(z, z0, z1, ramp)[:, None]
    return field * mask(z, th) if mask is not None else field


def pit_field(z, th, pits):
    """Smooth finger/thumb pits: sum of C1 flat-bottomed depressions, POSITIVE
    depths (subtract from an outer-radius field). Each pit is a dict:
    z, th (centre), sz, sa (half-widths, mm along z and arc), r (local surface
    radius for the arc metric), depth."""
    z = np.asarray(z, float)
    th = np.asarray(th, float)
    out = np.zeros((len(z), len(th)))
    for p in pits:
        dth = (th - p["th"] + math.pi) % (2 * math.pi) - math.pi
        q = (((z - p["z"]) / p["sz"]) ** 2)[:, None] \
            + ((p["r"] * dth / p["sa"]) ** 2)[None, :]
        s = np.clip(1.0 - q ** 1.5, 0.0, 1.0)     # ^1.5 flattens the floor
        out += p["depth"] * s * s * (3.0 - 2.0 * s)
    return out


def pit_depth_at(zs, ths, pits):
    """pit_field sampled at scattered sites (for sinking ornament that lands
    inside a pit deep enough to stay rooted)."""
    zs, ths = np.asarray(zs, float), np.asarray(ths, float)
    out = np.zeros(len(zs))
    for p in pits:
        dth = (ths - p["th"] + math.pi) % (2 * math.pi) - math.pi
        q = ((zs - p["z"]) / p["sz"]) ** 2 + (p["r"] * dth / p["sa"]) ** 2
        s = np.clip(1.0 - q ** 1.5, 0.0, 1.0)
        out += p["depth"] * s * s * (3.0 - 2.0 * s)
    return out


def claws(P, axis, h, r_base, r_tip=0.18, curl=None, curl_amt=0.6, nseg=7,
          sink=1.5, split=0.45):
    """Two-segment sharp claw per site: a straight root frustum, then a tip
    frustum bent toward `curl` -- genuinely sharp (r_tip ~ the finest point a
    tapering PETG print can hold) where frusta() was a blunt stud.

    Returns [roots, tips] as TWO meshes: within each, bodies are disjoint, but
    every tip overlaps its root -- and manifold3d treats a single mesh whose
    own components intersect as invalid input. Union both onto the shell:
    H.union(shell, *claws(...))."""
    P = np.asarray(P, float)
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    h = np.broadcast_to(np.asarray(h, float), (len(P),)).astype(float)
    r_base = np.broadcast_to(np.asarray(r_base, float), (len(P),)).astype(float)
    if curl is None:
        curl = np.array([0.0, 0.0, 1.0])
    curl = np.broadcast_to(np.asarray(curl, float), P.shape)
    tipd = axis + curl_amt * curl
    tipd = tipd / np.linalg.norm(tipd, axis=1, keepdims=True)
    r_mid = r_tip + (r_base - r_tip) * (1.0 - split)
    mid = P + axis * (h * split)[:, None]
    root = frusta(P, axis, h * split + 0.2, r_base, r_tip=r_mid, nseg=nseg,
                  sink=sink)
    tip = frusta(mid, tipd, h * (1.0 - split), r_mid, r_tip=r_tip, nseg=nseg,
                 sink=0.6)
    return [root, tip]


BRAILLE = {
    "a": "1", "b": "12", "c": "14", "d": "145", "e": "15", "f": "124",
    "g": "1245", "h": "125", "i": "24", "j": "245", "k": "13", "l": "123",
    "m": "134", "n": "1345", "o": "135", "p": "1234", "q": "12345",
    "r": "1235", "s": "234", "t": "2345", "u": "136", "v": "1236",
    "w": "2456", "x": "1346", "y": "13456", "z": "1356", " ": "",
    "#": "3456", ".": "256", "-": "36", "/": "34",
}
_DIGIT = {"1": "a", "2": "b", "3": "c", "4": "d", "5": "e",
          "6": "f", "7": "g", "8": "h", "9": "i", "0": "j"}
_DOTPOS = {"1": (0, 2), "2": (0, 1), "3": (0, 0),
           "4": (1, 2), "5": (1, 1), "6": (1, 0)}


def braille_dots(text, dot=2.5, cell=6.1):
    """Grade-1 braille dot centres (x along the line, y up), mm. RULE S4."""
    out, x, num = [], 0.0, False
    for ch in text.lower():
        if ch.isdigit():
            if not num:
                out += _cell(BRAILLE["#"], x, dot); x += cell; num = True
            key = _DIGIT[ch]
        else:
            num = False
            key = ch if ch in BRAILLE else " "
        out += _cell(BRAILLE[key], x, dot)
        x += cell
    return np.asarray(out, float).reshape(-1, 2), x


def _cell(code, x0, dot):
    return [(x0 + _DOTPOS[c][0] * dot, _DOTPOS[c][1] * dot) for c in code]


def wrap_dots(dots, z0, prof_z, prof_r, th0=0.0):
    """Map planar braille coords onto a surface of revolution."""
    z = z0 + dots[:, 1]
    r = np.interp(z, prof_z, prof_r)
    r0 = float(np.interp(z0, prof_z, prof_r))
    th = th0 + dots[:, 0] / r0
    P = np.stack([r * np.cos(th), r * np.sin(th), z], -1)
    dr = np.gradient(np.asarray(prof_r, float), np.asarray(prof_z, float))
    dri = np.interp(z, prof_z, dr)
    nrm = np.stack([np.cos(th), np.sin(th), -dri], -1)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    return P, nrm


# ---------------------------------------------------------------- lattices

def struts(seg, r, nseg=6):
    """One prism per segment; overlapping prisms are resolved by union()."""
    ph = grid_theta(nseg)
    V, F, off = [], [], 0
    for p0, p1 in seg:
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        a = p1 - p0
        L = float(np.linalg.norm(a))
        if L < 1e-6:
            continue
        a = a / L
        t = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u = np.cross(a, t); u /= np.linalg.norm(u)
        w = np.cross(a, u)
        ring = np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * w
        V.append(np.concatenate([p0 + r * ring, p1 + r * ring, [p0], [p1]]))
        j = np.arange(nseg); j1 = (j + 1) % nseg
        F.append(np.concatenate([
            np.stack([off + j, off + j1, off + nseg + j1], -1),
            np.stack([off + j, off + nseg + j1, off + nseg + j], -1),
            np.stack([np.full(nseg, off + 2 * nseg), off + j1, off + j], -1),
            np.stack([np.full(nseg, off + 2 * nseg + 1), off + nseg + j,
                      off + nseg + j1], -1)]))
        off += 2 * nseg + 2
    return trimesh.Trimesh(np.concatenate(V), np.concatenate(F), process=False)


def tube_path(pts, r, nseg=6, closed=False):
    """Swept tube along a polyline, parallel-transported frames.

    Lattices are built from these rather than from per-edge prisms: one tube per
    hoop or rail keeps the boolean input count (and therefore the union) small,
    and each tube is individually watertight so manifold3d accepts it."""
    P = np.asarray(pts, float)
    if closed:
        P = np.vstack([P, P[:1]])
    seg = P[1:] - P[:-1]
    T = seg / np.linalg.norm(seg, axis=1, keepdims=True)
    Tv = np.vstack([T[:1], 0.5 * (T[:-1] + T[1:]), T[-1:]])
    Tv /= np.linalg.norm(Tv, axis=1, keepdims=True)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(Tv[0] @ ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    u = np.cross(Tv[0], ref); u /= np.linalg.norm(u)
    U = [u]
    for i in range(1, len(Tv)):                 # parallel transport
        u = u - Tv[i] * float(u @ Tv[i])
        n = np.linalg.norm(u)
        u = U[-1] if n < 1e-9 else u / n
        U.append(u)
    U = np.asarray(U)
    W = np.cross(Tv, U)
    ph = grid_theta(nseg)
    ring = (np.cos(ph)[None, :, None] * U[:, None, :]
            + np.sin(ph)[None, :, None] * W[:, None, :])
    V = P[:, None, :] + r * ring
    if closed:
        return _loft_closed(V[:-1])
    return loft(V, cap_start=True, cap_end=True)


def _loft_closed(V):
    """Loft a ring-of-rings that closes on itself in the sweep direction."""
    ns, nt, _ = V.shape
    k = np.arange(ns)[:, None]
    k1 = (k + 1) % ns
    j = np.arange(nt)[None, :]
    j1 = (j + 1) % nt
    a, b = k * nt + j, k * nt + j1
    c, d = k1 * nt + j1, k1 * nt + j
    F = np.concatenate([np.stack([a, b, c], -1).reshape(-1, 3),
                        np.stack([a, c, d], -1).reshape(-1, 3)])
    return clean(trimesh.Trimesh(V.reshape(-1, 3), F, process=False))


def union(*meshes):
    ms = [m for m in meshes if m is not None and len(m.faces)]
    if not ms:
        return None
    if len(ms) == 1:
        return ms[0]
    return trimesh.boolean.union(ms, engine="manifold")


def intersection(*meshes):
    return trimesh.boolean.intersection(list(meshes), engine="manifold")


def difference(a, *cutters):
    cs = [c for c in cutters if c is not None and len(c.faces)]
    return trimesh.boolean.difference([a] + cs, engine="manifold") if cs else a


def cyl(p0, p1, r, nseg=32):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    return trimesh.creation.cylinder(radius=r, segment=np.stack([p0, p1]),
                                     sections=nseg)


# ---------------------------------------------------------------- reporting

def unpinch(mesh, eps=6e-3):
    """Pull apart vertices that two surface sheets share exactly.

    manifold3d is happy to return a solid whose sheets touch at a point -- the
    volume is unambiguous -- but STL is float32 and the two copies survive the
    round trip as one vertex, which reads as a non-manifold edge in every slicer
    that checks. Nudging each copy 6 um into its own sheet keeps the geometry and
    loses the pinch."""
    v = np.asarray(mesh.vertices, float).copy()
    key = np.round(v, 6)
    order = np.lexsort(key.T)
    same = np.all(key[order[1:]] == key[order[:-1]], axis=1)
    if not same.any():
        return mesh
    groups, run = [], [order[0]]
    for i, s in enumerate(same):
        if s:
            run.append(order[i + 1])
        else:
            if len(run) > 1:
                groups.append(run)
            run = [order[i + 1]]
    if len(run) > 1:
        groups.append(run)
    f = np.asarray(mesh.faces)
    cent = mesh.triangles.mean(1)
    for g in groups:
        for vi in g:
            owns = np.any(f == vi, axis=1)
            if not owns.any():
                continue
            d = cent[owns].mean(0) - v[vi]
            n = np.linalg.norm(d)
            if n > 1e-12:
                v[vi] += eps * d / n
    out = trimesh.Trimesh(v, f, process=False)
    return out


def footprint(mesh):
    """Bed footprint of the part as printed: the minimum-area rectangle round its
    XY hull. Parts are free to rotate about z on the plate, so the axis-aligned
    box overstates what a diagonal part like the Volute bell actually needs."""
    from shapely.geometry import MultiPoint
    r = MultiPoint(mesh.vertices[:, :2]).convex_hull.minimum_rotated_rectangle
    c = np.asarray(r.exterior.coords)[:-1]
    d = np.linalg.norm(np.diff(np.vstack([c, c[:1]]), axis=0), axis=1)
    return float(max(d[0], d[1])), float(min(d[0], d[1]))


def report(name, mesh, path=None, note=""):
    if path:                       # report on the FILE, not the float64 original:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        unpinch(mesh).export(path)  # STL is float32, and rounding can pinch a mesh
        mesh = trimesh.load(path)
    fx, fy = footprint(mesh)
    ext = np.array([fx, fy, mesh.extents[2]])
    row = dict(part=name, faces=len(mesh.faces), x=ext[0], y=ext[1], z=ext[2],
               vol_cm3=mesh.volume / 1000.0, watertight=bool(mesh.is_watertight),
               winding=bool(mesh.is_winding_consistent),
               bodies=int(mesh.body_count),
               fits_bed=bool((ext <= BED + 1e-6).all()), note=note)
    if path:
        row["file"] = path
    return row


HDR = (f"{'part':<28}{'faces':>8}{'X':>7}{'Y':>7}{'Z':>7}{'cm3':>8}"
       f"{'wt':>4}{'wnd':>4}{'bod':>4}{'bed':>4}")


def print_table(rows):
    print(HDR)
    print("-" * len(HDR))
    for r in rows:
        print(f"{r['part']:<28}{r['faces']:>8}{r['x']:>7.1f}{r['y']:>7.1f}"
              f"{r['z']:>7.1f}{r['vol_cm3']:>8.1f}"
              f"{'Y' if r['watertight'] else 'n':>4}"
              f"{'Y' if r['winding'] else 'n':>4}{r['bodies']:>4}"
              f"{'Y' if r['fits_bed'] else 'NO':>4}")
    print("-" * len(HDR))
    bad = [r["part"] for r in rows
           if not (r["watertight"] and r["fits_bed"] and r["bodies"] == 1)]
    print("all parts watertight, single-body, inside the 320 bed" if not bad
          else "CHECK: " + ", ".join(bad))
    return bad
