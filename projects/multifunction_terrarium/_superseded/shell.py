"""Outer shell: a sawtooth hexagonal wall carrying eye-aimed viewing ports.

Why a sawtooth. Aiming a port at the viewer means tilting it back, which makes the
wall narrow as it rises. Integrated over the whole object that is 171 mm of inward
excursion -- the thing would taper to nothing and lose its footprint. So each aimed
port band is followed by a riser leaning back OUT by the same amount. Net silhouette
vertical, footprint constant at 310, and every band's normal points between a standing
and a seated eye (worst-case incidence 14 deg, against a 45 deg Fresnel budget).

Each band's lower edge is a sharp convex crease: a Gibbs pinning line, i.e. a drip lip,
with a 70 mm free fall to the next one.
"""

import math
import numpy as np
import trimesh
from shapely.geometry import LineString

import params as P

BANDS_PER_MODULE = 3
TARGET_OVERHANG = 40.0         # deg from vertical, with margin on P.OVERHANG_OK
SHELL_W = 2.4                  # structural frame; the panes carry no load
FRAME_SIDE = 14.0              # rib left either side of a port, per face
FRAME_END  = 6.0               # rib above and below a port


def band_schedule(iters=4):
    """Solve each band's port/riser split.

    The port height is not free: aiming a band at tilt t costs port_h*tan(t) of
    radius, and the riser has to give exactly that back within TARGET_OVERHANG.
    So riser_h = port_h*tan(t)/tan(overhang), and the band heights within a module
    are whatever makes them sum to MOD_H. Tilt depends on height and height on tilt,
    so iterate -- it converges in two passes.

    Consequence worth knowing: the steeply-aimed bands at the bottom get short
    windows and tall risers, the near-vertical bands at the top get tall windows.
    The object's biggest panes end up where the eye rests anyway."""
    tan_o = math.tan(math.radians(TARGET_OVERHANG))
    out = []
    for mod in range(P.N_MODULE):
        h = [P.MOD_H / BANDS_PER_MODULE] * BANDS_PER_MODULE   # seed: equal periods
        for _ in range(iters):
            z, tilts = 0.0, []
            for b in range(BANDS_PER_MODULE):
                zc = P.TABLE_H + mod * P.MOD_H + z + h[b] / 2
                tilts.append(P.band_tilt(zc)[0])
                z += h[b]
            k = [1.0 / (1.0 + math.tan(math.radians(t)) / tan_o) for t in tilts]
            ks = sum(k)
            h = [P.MOD_H * ki / ks / (1.0 / (1.0 + math.tan(math.radians(t)) / tan_o))
                 * (1.0 / (1.0 + math.tan(math.radians(t)) / tan_o)) for ki, t in zip(k, tilts)]
            h = [P.MOD_H * ki / ks for ki in k]
        z = 0.0
        for b in range(BANDS_PER_MODULE):
            zc = P.TABLE_H + mod * P.MOD_H + z + h[b] / 2
            tilt, worst = P.band_tilt(zc)
            port_h = h[b] / (1.0 + math.tan(math.radians(tilt)) / tan_o)
            out.append(dict(module=mod, band=b, z_local=z, period=h[b],
                            port_h=port_h, riser_h=h[b] - port_h,
                            z_world=zc, tilt=tilt, worst=worst,
                            dr=port_h * math.tan(math.radians(tilt))))
            z += h[b]
    return out


def _profile(mod):
    """(r,z) polyline of one module's outer face, bottom to top, plus port spans."""
    sched = [s for s in band_schedule() if s["module"] == mod]
    r = P.HEX_R
    pts = [(r, 0.0)]
    ports = []
    z = 0.0
    for s in sched:
        r_top = r - s["dr"]
        ports.append(dict(z0=z, z1=z + s["port_h"], r0=r, r1=r_top, **s))
        z += s["port_h"]
        pts.append((r_top, z))              # aimed band: leans in
        z += s["riser_h"]
        pts.append((r, z))                  # riser: leans back out
    return pts, ports, sched


def _hexf(t):
    return math.cos(math.pi / 6) / np.cos((t % (math.pi / 3)) - math.pi / 6)


def _sweep(pts, n_face=14, wall=SHELL_W):
    """Sweep a closed (r,z) section around a hexagon."""
    dense = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        n = max(1, int(math.ceil(np.linalg.norm(b - a) / 2.0)))
        for i in range(n):
            dense.append(tuple(a + (b - a) * i / n))
    dense.append(tuple(pts[-1]))
    poly = LineString(dense).buffer(wall / 2, cap_style=2, join_style=2, mitre_limit=6)
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    sec = np.asarray(poly.exterior.coords)[:-1]

    # azimuth: land samples exactly on the six corners, n_face between
    T = np.concatenate([np.linspace(k * math.pi / 3, (k + 1) * math.pi / 3,
                                    n_face, endpoint=False) for k in range(6)])
    hf = _hexf(T)
    R = sec[:, 0][None, :] * hf[:, None]
    X, Y = R * np.cos(T)[:, None], R * np.sin(T)[:, None]
    Z = np.broadcast_to(sec[:, 1][None, :], R.shape)
    V = np.stack([X, Y, Z], -1).reshape(-1, 3)

    nt, ns = len(T), len(sec)
    i = np.arange(nt)[:, None]; j = np.arange(ns)[None, :]
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
    return m


def _aperture(port, face):
    """One trapezoidal cutter, lying in the tilted plane of a port on one hex face."""
    ang = face * math.pi / 3 + math.pi / 6          # face-centre azimuth
    nx, ny = math.cos(ang), math.sin(ang)
    tx, ty = -ny, nx                                # along the face
    a = math.cos(math.pi / 6)                       # face radius = a * circumradius
    z0, z1 = port["z0"] + FRAME_END, port["z1"] - FRAME_END
    f = ((z0 - port["z0"]) / port["port_h"], (z1 - port["z0"]) / port["port_h"])
    r0 = (port["r0"] + (port["r1"] - port["r0"]) * f[0]) * a
    r1 = (port["r0"] + (port["r1"] - port["r0"]) * f[1]) * a
    # half-widths: the face is 2*a*R*tan(30) wide at radius r
    w0 = r0 * math.tan(math.pi / 6) - FRAME_SIDE
    w1 = r1 * math.tan(math.pi / 6) - FRAME_SIDE
    inn, out = -6.0, 12.0                           # cut clean through the wall
    quad = []
    for dr in (inn, out):
        quad += [(nx * (r0 + dr) - tx * w0, ny * (r0 + dr) - ty * w0, z0),
                 (nx * (r0 + dr) + tx * w0, ny * (r0 + dr) + ty * w0, z0),
                 (nx * (r1 + dr) + tx * w1, ny * (r1 + dr) + ty * w1, z1),
                 (nx * (r1 + dr) - tx * w1, ny * (r1 + dr) - ty * w1, z1)]
    V = np.array(quad)
    F = np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],
                  [1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]])
    m = trimesh.Trimesh(vertices=V, faces=F, process=False)
    if m.volume < 0:
        m.invert()
    return m


def build(mod=0, verbose=True):
    pts, ports, sched = _profile(mod)
    body = _sweep(pts)
    cutters = [_aperture(p, f) for p in ports for f in range(6)]
    shell = trimesh.boolean.difference([body] + cutters, engine="manifold")
    if verbose:
        print(f"  module {mod}: {len(ports)} bands x 6 faces = {len(cutters)} ports")
        for p in ports:
            lean = math.degrees(math.atan2(p['dr'], p['riser_h']))
            glass = p['port_h'] - 2 * FRAME_END
            print(f"    band z {p['z0']:5.0f}-{p['z1']:5.0f}  tilt {p['tilt']:5.1f}  "
                  f"incid {p['worst']:4.1f}  port {p['port_h']:5.1f} (glass {glass:4.1f})  "
                  f"riser {p['riser_h']:5.1f}  overhang {lean:4.1f} "
                  f"{'OK' if lean < P.OVERHANG_OK else '*** SUPPORT ***'}")
        print(f"  mesh {len(shell.faces):,} faces  watertight={shell.is_watertight}  "
              f"extents {np.round(shell.extents,1)}  volume {shell.volume/1000:.0f} cm3")
    return shell, ports


if __name__ == "__main__":
    print("outer shell:")
    for mod in range(P.N_MODULE):
        m, ports = build(mod)
        m.export(f"shell_module_{mod}.stl")
    print("  wrote shell_module_{0,1,2}.stl")
