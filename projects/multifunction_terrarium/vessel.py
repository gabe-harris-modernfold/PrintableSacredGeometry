"""The vessel: three printable shell modules, the reservoir pan, and the living-bed
tray with its bell siphon. Ovoid crown over a 12-sided tapered body."""

import math
import numpy as np
import trimesh
from shapely.geometry import LineString
import params as P
import spiral

FRAME_SIDE, FRAME_END = 13.0, 9.0
PORT_BANDS = {0: 2, 1: 2, 2: 1}          # per module


def _sweep(pts, wall, n_face=P.N_FACE, per_face=9, close=True):
    dense = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        n = max(1, int(math.ceil(np.linalg.norm(b - a) / 2.5)))
        dense += [tuple(a + (b - a) * i / n) for i in range(n)]
    dense.append(tuple(pts[-1]))
    poly = LineString(dense).buffer(wall / 2, cap_style=2, join_style=2, mitre_limit=6)
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    sec = np.asarray(poly.exterior.coords)[:-1]
    T = np.concatenate([np.linspace(k * 2 * math.pi / n_face,
                                    (k + 1) * 2 * math.pi / n_face,
                                    per_face, endpoint=False) for k in range(n_face)])
    R = sec[:, 0][None, :] * P.hexf(T, n_face)[:, None]
    V = np.stack([R * np.cos(T)[:, None], R * np.sin(T)[:, None],
                  np.broadcast_to(sec[:, 1][None, :], R.shape)], -1).reshape(-1, 3)
    nt, ns = len(T), len(sec)
    i = np.arange(nt)[:, None]; j = np.arange(ns)[None, :]
    a = i * ns + j; b = ((i + 1) % nt) * ns + j
    c = ((i + 1) % nt) * ns + (j + 1) % ns; d = i * ns + (j + 1) % ns
    F = np.concatenate([np.stack([a, b, c], -1).reshape(-1, 3),
                        np.stack([a, c, d], -1).reshape(-1, 3)])
    m = trimesh.Trimesh(vertices=V, faces=F, process=False)
    m.merge_vertices()
    if m.volume < 0: m.invert()
    return m


def _port(face, z0, z1, tilt_r0, tilt_r1):
    """Trapezoidal aperture cutter on one face."""
    ang = (face + 0.5) * 2 * math.pi / P.N_FACE
    nx, ny = math.cos(ang), math.sin(ang)
    tx, ty = -ny, nx
    a = math.cos(math.pi / P.N_FACE)
    r0, r1 = tilt_r0 * a, tilt_r1 * a
    w0 = r0 * math.tan(math.pi / P.N_FACE) - FRAME_SIDE
    w1 = r1 * math.tan(math.pi / P.N_FACE) - FRAME_SIDE
    if w0 <= 2 or w1 <= 2: return None
    V = []
    for dr in (-9.0, 14.0):
        V += [(nx*(r0+dr)-tx*w0, ny*(r0+dr)-ty*w0, z0),
              (nx*(r0+dr)+tx*w0, ny*(r0+dr)+ty*w0, z0),
              (nx*(r1+dr)+tx*w1, ny*(r1+dr)+ty*w1, z1),
              (nx*(r1+dr)-tx*w1, ny*(r1+dr)-ty*w1, z1)]
    F = np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],
                  [1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]])
    m = trimesh.Trimesh(vertices=np.array(V), faces=F, process=False)
    if m.volume < 0: m.invert()
    return m


def shell(mod, verbose=True):
    z0, z1 = P.SPLIT_Z[mod], P.SPLIT_Z[mod + 1]
    zz = np.linspace(z0, z1, 60)
    pts = list(zip(P.vessel_r(zz), zz))
    body = _sweep(pts, P.SHELL_W)
    cutters, ports = [], []
    nb = PORT_BANDS[mod]
    span = (z1 - z0) / nb
    for b in range(nb):
        a0 = z0 + b * span + FRAME_END
        a1 = z0 + (b + 1) * span - FRAME_END
        if mod == 2:                                    # crown: shorter, higher up
            a0, a1 = z0 + 14, z0 + 96
        for f in range(P.N_FACE):
            c = _port(f, a0, a1, float(P.vessel_r(a0)), float(P.vessel_r(a1)))
            if c is not None:
                cutters.append(c); ports.append((f, a0, a1))
    if mod == 2:
        # Apex penetration for the condenser tap's intake stem. ONE small sealed
        # hole: the collector cone seats against the inner face below it, so only
        # the 11.8 mm stem crosses the shell. (With the tap's root outside the
        # vessel this had to be ~40 mm wide to clear two diverged branches, which
        # would have vented the crown and killed the condensate the tap collects.)
        r_stem = spiral.BORE * 4 / 2 + spiral.GWALL             # 5.9
        # The boss must start ABOVE the cone's seat (z=539.07) or it hangs into the
        # crown and lands inside the cone: that collision is 2.1 cm3, and it is the
        # boss that moves, not the seat -- the seat height IS the contact.
        z_boss = spiral.crown_in_z(spiral.SEAT_R) + 0.4
        body = trimesh.boolean.union(
            [body, trimesh.creation.cylinder(radius=r_stem + 5.1, height=15,
                                             sections=32)
                   .apply_translation([0, 0, z_boss + 7.5])], engine="manifold")
        cutters.append(trimesh.creation.cylinder(radius=r_stem + 0.4, height=30,
                                                 sections=32)
                       .apply_translation([0, 0, z_boss + 8.0]))
    out = trimesh.boolean.difference([body] + cutters, engine="manifold")
    if verbose:
        print(f"  shell module {mod} (z {z0:.0f}-{z1:.0f}): {len(ports)} ports, "
              f"{len(out.faces):,} faces, watertight={out.is_watertight}, "
              f"{np.round(out.extents,1)}, {out.volume/1000:.0f} cm3")
    return out, ports


def reservoir(verbose=True):
    """Pan holding the standing water, with a foot and a pump bay."""
    zz = np.linspace(0, P.Z_RES + 16, 40)
    outer = P.vessel_r(zz) - P.SHELL_W - 1.2
    pts = list(zip(outer, zz))
    pan = _sweep(pts, 2.6)
    floor = _sweep([(outer[0], 0.0), (2.0, 0.0)], 3.2)
    bay = trimesh.creation.cylinder(radius=26, height=P.Z_RES,
                                    sections=32).apply_translation([0, 0, P.Z_RES/2])
    m = trimesh.boolean.union([pan, floor], engine="manifold")
    m = trimesh.boolean.difference([m, bay], engine="manifold")
    vol = math.pi * (outer.mean() * math.cos(math.pi/P.N_FACE)) ** 2 * P.Z_RES / 1e6
    if verbose:
        print(f"  reservoir: holds ~{vol:.2f} L to the waterline at z={P.Z_RES:.0f}, "
              f"{len(m.faces):,} faces, watertight={m.is_watertight}")
    return m


def bell_siphon(bell_r=26.0, stand_r=9.0, h=50.0):
    """Standpipe + bell + snorkel: a genuine relaxation oscillator, no timer.

    h is the standpipe height and therefore the flood depth. It must leave the bell
    cap (h+30) BELOW the bed-tray rim, or the tray overflows before the siphon ever
    trips. At h=54 the cap lands at 203 against a 205 rim."""
    parts = []
    parts.append(trimesh.creation.annulus(r_min=stand_r, r_max=stand_r+2.2,
                                          height=h, sections=40)
                 .apply_translation([0, 0, h/2]))
    parts.append(trimesh.creation.annulus(r_min=stand_r+7, r_max=stand_r+9.5,
                                          height=7, sections=40)
                 .apply_translation([0, 0, h - 3.5]))        # flare
    parts.append(trimesh.creation.annulus(r_min=bell_r, r_max=bell_r+2.4,
                                          height=h+22, sections=48)
                 .apply_translation([0, 0, (h+22)/2 + 7]))
    parts.append(trimesh.creation.cylinder(radius=bell_r+2.4, height=2.4, sections=48)
                 .apply_translation([0, 0, h + 29 + 1.2]))    # bell cap
    parts.append(trimesh.creation.annulus(r_min=3.0, r_max=4.6, height=52, sections=20)
                 .apply_translation([bell_r - 8, 0, 26]))     # snorkel
    for k in range(6):                                        # inlet feet
        a = k * math.pi / 3
        parts.append(trimesh.creation.box(extents=[7, 3.2, 7])
                     .apply_translation([(bell_r+1)*math.cos(a),
                                         (bell_r+1)*math.sin(a), 3.5]))
    return trimesh.boolean.union(parts, engine="manifold")


def bed_tray(verbose=True):
    """Living bed: >=30 deg drain slope over an air gap, so it trickles and never ponds."""
    r_out = float(P.vessel_r(P.Z_BED)) - P.SHELL_W - 1.2
    zz_top = P.Z_BED
    slope = math.radians(32)
    r_in = 30.0
    drop = (r_out - r_in) * math.tan(slope)
    pts = [(r_out, zz_top), (r_in, zz_top - drop), (r_in, zz_top - drop - 14),
           (r_out, zz_top - 8)]
    cone = _sweep(pts, 2.4)
    m = trimesh.boolean.union([cone, bell_siphon().apply_translation(
        [0, 0, zz_top - drop - 14])], engine="manifold")
    if verbose:
        print(f"  bed tray: r {r_in:.0f}-{r_out:.0f}, {math.degrees(slope):.0f} deg slope, "
              f"drop {drop:.0f} mm, {len(m.faces):,} faces")
    return m


if __name__ == "__main__":
    print("vessel:")
    for k in range(3):
        m, _ = shell(k); m.export(f"vessel_shell_{k}.stl")
    reservoir().export("vessel_reservoir.stl")
    bed_tray().export("vessel_bedtray.stl")
    print("  wrote vessel_shell_{0,1,2}.stl, vessel_reservoir.stl, vessel_bedtray.stl")
