"""The cascade screen: a staircase cone that forces the same water to detach ~21
times per module.

Geometry rules, all from params:
  * riser STEP = 10 mm  -- as short as STEP_MIN allows, because drops-in-air ~ 1/sqrt(STEP)
  * tread TREAD >= CAPILLARY_LEN, and each lip sits directly over the next tread's
    inner edge so every drop lands and re-forms instead of falling clear
  * lips are scalloped into travertine lobes purely to buy arc length (= drip sites)
    without moving the cone out toward the shell
  * notch teeth at NOTCH_P along the lip, tip = NOZZLE wide, so Tate gives the
    smallest drop the printer can make
"""

import math
import numpy as np
import trimesh
from shapely.geometry import LineString

import params as P


def _staircase(levels=P.LEVELS, r0=P.SCREEN_R_TOP, tread=P.TREAD, step=P.STEP,
               z_top=P.MOD_H, tread_slope=8.0):
    """(r,z) polyline top->bottom, plus the (r,z) of every lip corner.

    Each lip is directly above the next tread's inner edge -- that is what turns one
    fall into twenty-one."""
    drop = tread * math.tan(math.radians(tread_slope))
    pts, lips = [], []
    for k in range(levels):
        rk = r0 + k * tread
        zk = z_top - k * step
        pts.append((rk, zk))                    # tread inner
        pts.append((rk + tread, zk - drop))     # lip: outer, and slightly lower
        lips.append((rk + tread, zk - drop))
        pts.append((rk + tread, zk - step))     # straight down to the next tread
    pts.append((r0 + levels * tread, z_top - levels * step))
    return pts, lips


def _densify(pts, lips, fine=0.35, coarse=1.8, near=2.5):
    """Dense sampling only where the teeth need it."""
    L = np.asarray(lips)
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.asarray(a), np.asarray(b)
        seg = np.linalg.norm(b - a)
        if seg < 1e-9:
            continue
        mid = (a + b) / 2
        d = np.min(np.linalg.norm(L - mid, axis=1))
        h = fine if d < near else coarse
        n = max(1, int(math.ceil(seg / h)))
        for i in range(n):
            out.append(tuple(a + (b - a) * i / n))
    out.append(tuple(pts[-1]))
    return out


def _section(levels=P.LEVELS, **kw):
    """Closed wall cross-section in (r,z), plus a per-vertex lip weight."""
    pts, lips = _staircase(levels=levels, **kw)
    dense = _densify(pts, lips)
    poly = LineString(dense).buffer(P.SCREEN_W / 2, cap_style=2, join_style=2,
                                    mitre_limit=6.0)
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    sec = np.asarray(poly.exterior.coords)[:-1]
    L = np.asarray(lips)
    d = np.min(np.linalg.norm(sec[:, None, :] - L[None, :, :], axis=2), axis=1)
    w = np.clip(1.0 - d / 1.6, 0.0, 1.0) ** 0.6          # teeth live only at the lips
    w /= w.max()   # the buffer offsets every vertex off the lip corner by SCREEN_W/2,
                   # so without this the teeth build at ~82% of the depth asked for
    return sec, w, lips


def _tooth(u, tip_frac):
    """One notch: ramp out, a NOZZLE-wide flat tip, ramp back. u in [0,1)."""
    a, b = 0.5 - tip_frac / 2, 0.5 + tip_frac / 2
    return np.where(u < a, u / a,
           np.where(u < b, 1.0, (1.0 - u) / (1.0 - b)))


def _azimuth(n_teeth, tip_frac, extra=2):
    """Sample theta at the notch BREAKPOINTS, not on a uniform grid.

    A 0.4 mm tip inside a 6 mm pitch is 6.7% of a period -- a uniform grid fine
    enough to resolve it would need ~30 samples per tooth and put the mesh over
    9 M faces. Placing samples exactly at {root, tip start, tip end} gives an
    exact 0.4 mm tip for 5 samples per tooth."""
    a, b = 0.5 - tip_frac / 2, 0.5 + tip_frac / 2
    u = np.sort(np.concatenate([[0.0, a, b],
                                np.linspace(0, 1, extra, endpoint=False)[1:]]))
    U = (np.arange(n_teeth)[:, None] + u[None, :]).ravel()
    return U * 2 * math.pi / n_teeth, U % 1.0


def build(levels=P.LEVELS, tooth_d=1.4, verbose=True, **kw):
    sec, lipw, lips = _section(levels=levels, **kw)
    r_sec, z_sec = sec[:, 0], sec[:, 1]

    # arc length of one scalloped lip at the mean radius -> how many teeth fit
    r_mean = float(np.mean([p[0] for p in lips]))
    th = np.linspace(0, 2 * math.pi, 20000, endpoint=False)
    hf = (math.cos(math.pi / 6) / np.cos((th % (math.pi / 3)) - math.pi / 6)
          if P.SCREEN_HEX else np.ones_like(th))
    rr = r_mean * hf + P.LOBE_A * np.cos(P.LOBES * th)
    dr = np.gradient(rr, th)
    lip_len = float(np.trapezoid(np.sqrt(rr ** 2 + dr ** 2), th))
    n_teeth = max(24, int(round(lip_len / P.NOTCH_P)))
    tip_frac = P.LIP_W / (lip_len / n_teeth)

    T, u = _azimuth(n_teeth, tip_frac)
    nt = len(T)
    hexf = (math.cos(math.pi / 6) / np.cos((T % (math.pi / 3)) - math.pi / 6)
            if P.SCREEN_HEX else np.ones_like(T))
    scallop = P.LOBE_A * np.cos(P.LOBES * T)
    teeth = _tooth(u, tip_frac) * tooth_d

    # R(theta, s) = staircase radius + travertine lobe + notch (lips only)
    R = (r_sec[None, :] * hexf[:, None] + scallop[:, None]
         + lipw[None, :] * teeth[:, None])
    X = R * np.cos(T)[:, None]
    Y = R * np.sin(T)[:, None]
    Z = np.broadcast_to(z_sec[None, :], R.shape)
    V = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)

    ns = len(sec)
    i = np.arange(nt)[:, None]
    j = np.arange(ns)[None, :]
    a = i * ns + j
    b = ((i + 1) % nt) * ns + j
    c = ((i + 1) % nt) * ns + (j + 1) % ns
    d = i * ns + (j + 1) % ns
    F = np.concatenate([np.stack([a, b, c], -1).reshape(-1, 3),
                        np.stack([a, c, d], -1).reshape(-1, 3)])

    m = trimesh.Trimesh(vertices=V, faces=F, process=False)
    m.merge_vertices()
    if m.volume < 0:
        m.invert()
    if verbose:
        rate_t = P.per_site_rate(P.Q_TRICKLE, n_teeth)
        rate_d = P.per_site_rate(P.Q_DUMP, n_teeth)
        print(f"  lip arc length   {lip_len:7.1f} mm  (plain circle would be "
              f"{2*math.pi*r_mean:.0f})")
        print(f"  notches / level  {n_teeth:7d}   pitch {lip_len/n_teeth:.2f} mm  "
              f"tip {tip_frac*lip_len/n_teeth:.2f} mm  depth {tooth_d} mm "
              f"(tread {P.TREAD:.2f}, so drops land mid-tread)")
        print(f"  azimuth samples  {nt:7d}   {nt/n_teeth:.0f} per tooth, "
              f"placed on the notch breakpoints")
        print(f"  detach events    {levels} levels x {n_teeth} = "
              f"{levels*n_teeth} drip sites")
        print(f"  per-site rate    {rate_t:6.2f}/s trickle   {rate_d:5.2f}/s dump   "
              f"(jet limit {P.JET_LIMIT:.1f}) "
              f"{'OK' if rate_d < P.JET_LIMIT else '*** STREAMS ***'}")
        print(f"  mesh             {len(m.faces):,} faces   watertight="
              f"{m.is_watertight}  winding={m.is_winding_consistent}")
        print(f"  extents          {np.round(m.extents,1)}   volume "
              f"{m.volume/1000:.0f} cm3")
    return m, dict(n_teeth=n_teeth, lip_len=lip_len, sites=levels * n_teeth)


if __name__ == "__main__":
    print("cascade screen, one module:")
    m, info = build()
    m.export("cascade_screen.stl")
    print("  wrote cascade_screen.stl")
