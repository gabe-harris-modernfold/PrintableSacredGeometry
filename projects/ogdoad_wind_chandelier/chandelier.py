#!/usr/bin/env python3
"""The Hymn of the Ogdoad -- a wind-blown Helmholtz chandelier, CAD.

Nine pigeon-whistle voices tuned to the nine Solfeggio tones, hung as a geocentric
cosmos. Earth (174 Hz), the largest gourd, sits at the centre on the hub. Around it a
wheel carries the seven planetary spheres in the Chaldean order Poimandres climbs
(Moon 285, Mercury 396, Venus 417, Sun 528, Mars 639, Jupiter 741, Saturn 852) and the
eighth sphere, the Ogdoad (963). They sit as radial rays in four tiers that climb the
wheel. Each tier's pair is split left/right so the wheel balances, which makes the
ascent zig-zag across the axis like the serpents of a caduceus. A dove-tail vane on a
boom keeps every intake facing the wind, a crown with a comb of holes takes the
hanging swivel, and a plumb-bob filled with sand or steel shot hangs below as ballast.

Coordinates (assembly): the wheel lies in the XZ plane, z up. The wind blows toward +y,
so every intake faces -y. Each gourd has its own local frame: z runs along its axis,
from the foot to the neck; y runs downwind; x = y cross z.

Acoustic sizes come from voices.py. This script builds the solids, re-measures each
cavity volume against the design (and so each pitch), checks the assembly clearances
and the bed fit, exports print-oriented STLs and renders previews.

Run from this directory:  python chandelier.py [--no-preview]
"""

import math
import os
import sys

import numpy as np
import trimesh
from build123d import (Align, Axis, Box, Cylinder, Location, Plane, Polygon, Pos, Rot,
                       export_stl, extrude, revolve)

import voices as vc

HERE = os.path.dirname(os.path.abspath(__file__))
C, MIN, MAX = Align.CENTER, Align.MIN, Align.MAX
BED = 320.0
PETG = 1.27                     # g/cm^3

# ---------------------------------------------------------------- frame parameters (mm)

T_F = 10.0                      # wheel thickness along y (y in [-5, 5]); prints rear-face down
HUB_R = 24.0
HUB_Y0 = -15.0                  # hub boss reaches forward to here; Earth's peg starts here
SPOKE_W = 9.0
RING_W = 7.0
PEG_X, PEG_Y, PEG_L = 7.0, 8.0, 7.5     # satellite foot peg (tangential x radial-normal y)
EPEG = 10.0                     # Earth's square peg
EPEG_L = 12.0
L_BOSS = 16.0                   # Earth's rear boss, from its bulb surface to the hub
BOSS_HW = 11.0                  # boss half width / half height
BOOM = 12.0                     # vane boom section
BOOM_IN = 14.0                  # boom depth in the hub
BOOM_FREE = 220.0               # boom span between the hub and the fin sleeve
SLEEVE = 40.0
FIN_T = 2.0
LUG_W = 16.0                    # crown / pendant lug width
PIN_D = 3.4                     # M3 through-bolts
GAP = 6.0                       # minimum clearance between neighbouring voices
RF_MIN = 105.0                  # open the wheel past the packing minimum: an airier rose,
                                # less wake between rays, ~0.55 m across the ray tips

# the climb: four tiers of two, each at +-22.5 deg off the vertical axis
TIERS = [(247.5, 292.5), (202.5, 337.5), (157.5, 22.5), (112.5, 67.5)]

# ---------------------------------------------------------------- helpers


def yz_prism(pts, x0, x1):
    """Polygon in the (y, z) plane extruded over x in [x0, x1]."""
    face = Plane.YZ.offset(x0) * Polygon(*[tuple(p) for p in pts], align=None)
    return extrude(face, amount=x1 - x0)


def box(x0, x1, y0, y1, z0, z1):
    return Pos((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2) * Box(x1 - x0, y1 - y0, z1 - z0)


def cyl_z(r, z0, z1, x=0.0, y=0.0):
    return Pos(x, y, z0) * Cylinder(r, z1 - z0, align=(C, C, MIN))


def cyl_y(r, y0, y1, x=0.0, z=0.0):
    return Pos(x, (y0 + y1) / 2, z) * Rot(90, 0, 0) * Cylinder(r, y1 - y0)


def rev_rz(pts):
    return revolve(Plane.XZ * Polygon(*[tuple(p) for p in pts], align=None), Axis.Z)


def arc_pts(cy, cz, r, a0, a1, n=10):
    a = np.radians(np.linspace(a0, a1, n))
    return [(cy + r * math.cos(t), cz + r * math.sin(t)) for t in a]


def frame_of(theta):
    """Location placing a satellite's local frame on the wheel ray at angle theta (deg)."""
    t = math.radians(theta)
    z = (math.cos(t), 0.0, math.sin(t))
    x = (math.sin(t), 0.0, -math.cos(t))
    return z, x


# ---------------------------------------------------------------- the voice: gourd + head


def flat_dist(d):
    """D-flat on the spigot (and in the head socket), at local +x: keys the head's
    orientation so its intake can only face upwind. Leaves >= 1 mm of wall at the flat."""
    return max(vc.D_FLAT * d["rno"], d["rb"] + 1.0)


def gourd(v, socket=True):
    d = v["dims"]
    outer, inner = v["outer"], v["inner"]
    pts = [(0.0, 0.0)] + [tuple(p) for p in outer] + [tuple(p) for p in inner[::-1]] + \
          [(0.0, vc.FLOOR_T)]
    body = rev_rz(pts)
    # D-flat on the spigot
    zf0 = d["z_top"] - vc.SPIG_H
    fd = flat_dist(d)
    body -= box(fd, d["rno"] + 2, -d["rno"] - 2, d["rno"] + 2, zf0, d["z_top"] + 1)
    if socket:
        # blind foot socket for the wheel's peg, offset toward +y: one way in only
        body -= box(-vc.SOCK_X / 2, vc.SOCK_X / 2, vc.SOCK_OFF - vc.SOCK_Y / 2,
                    vc.SOCK_OFF + vc.SOCK_Y / 2, -1, vc.SOCK_D)
    return body


def earth_boss(v):
    """Rear corbel on Earth's lower bulb (local +y), with the hub peg's socket.
    Its underside is chamfered at 45 deg so the gourd still prints upright unsupported."""
    d = v["dims"]
    R, t, c1 = d["R"], d["t"], d["c1"]
    y0, y1 = R - t / 2, R + L_BOSS
    zt, zb = c1 + BOSS_HW, c1 - BOSS_HW
    prof = [(y0, zt), (y0, zb - (y1 - y0)), (y1, zb), (y1, zt)]
    boss = yz_prism(prof, -BOSS_HW, BOSS_HW)
    sock = EPEG + 0.4
    boss -= box(-sock / 2, sock / 2, y1 - EPEG_L - 1, y1 + 1, c1 - sock / 2, c1 + sock / 2)
    return boss


def head_geometry(v):
    """Head stations, with the windway / ramp stretched so the head covers its skirt."""
    d = v["dims"]
    hs = dict(v["head"])
    ys = hs["ys"]
    # big voices: the skirt is longer than the jet path; lengthen the windway and ramp so
    # the intake and ramp end flush with the skirt (no ledge in front of the intake)
    if hs["y_front"] > -ys:
        hs["y_front"] = -ys
        hs["y_ww0"] = -ys + vc.R_LIP + vc.T_BLK
    if hs["y_rear"] < ys:
        hs["y_rear"] = ys
        hs["z_ramp"] = hs["z_tip"] + (ys - hs["yL"]) * math.tan(math.radians(vc.LAB_RAMP))
    hs["flat"] = flat_dist(d)
    return hs


def head(v):
    """Whistle head, local frame of its gourd shifted so z=0 is the spigot top.
    Prints on its -x side, so the jet's whole profile (windway, window, labium) lies in
    the bed plane and comes out crisp."""
    b = v["b"]
    hs = head_geometry(v)
    xh, ys = hs["xh"], hs["ys"]
    z_f, z_r, z_bt = hs["z_f"], hs["z_r"], hs["z_bt"]
    y0, yL, yf, yw0, yr = hs["y0"], hs["yL"], hs["y_front"], hs["y_ww0"], hs["y_rear"]
    R = vc.R_LIP
    # deck + skirt: covers the whole neck bore and holds the socket
    solid = box(-xh, xh, -ys, ys, -vc.SPIG_H, z_f)
    # the jet body: intake nose, block over the windway, open window, labium ramp
    Rn = R + vc.T_BLK
    prof = ([(yf, 0.0), (yr, 0.0), (yr, hs["z_ramp"]), (yL, hs["z_tip"]), (yL, z_f),
             (y0, z_f), (y0, z_bt)]
            + arc_pts(yw0, z_r + R, Rn, 90, 180, 10))
    solid += yz_prism(prof, -xh, xh)
    # windway with bellmouth lips
    ww = ([(yf - 1.0, z_f - R)]
          + arc_pts(yw0, z_f - R, R, 180, 90, 10)
          + [(y0 + 0.3, z_f), (y0 + 0.3, z_r)]
          + arc_pts(yw0, z_r + R, R, 270, 180, 10)
          + [(yf - 1.0, z_r + R)])
    solid -= yz_prism(ww, -b / 2, b / 2)
    # plenum + window, down to the bore
    solid -= yz_prism(vc.plenum_poly(hs), -b / 2, b / 2)
    # D socket over the spigot
    rs = hs["rso"] + vc.FIT
    sock = cyl_z(rs, -vc.SPIG_H - 1, 0.0) - box(hs["flat"] + vc.FIT, rs + 1, -rs - 1, rs + 1,
                                                 -vc.SPIG_H - 2, 1)
    solid -= sock
    return solid, hs


def cavity(v, hs):
    """The air: gourd interior + plenum. Measured to confirm the pitch."""
    d = v["dims"]
    air = rev_rz([(0.0, vc.FLOOR_T)] + [tuple(p) for p in v["inner"]] + [(0.0, d["z_top"])])
    pl = Pos(0, 0, d["z_top"]) * yz_prism(vc.plenum_poly(hs), -v["b"] / 2, v["b"] / 2)
    return air + pl


# ---------------------------------------------------------------- layout


def solve_layout(vs):
    """Assign each tier's pair to left/right for balance, then find the smallest wheel
    radius at which no two voices come within GAP of each other or of Earth."""
    sats = vs[1:]
    pairs = [(sats[2 * k], sats[2 * k + 1]) for k in range(4)]
    # mass moment about the vertical axis, with each voice's centroid ~ R_f + c1 + ...
    best = None
    for mask in range(16):
        mom = 0.0
        assign = []
        for k, (a, b) in enumerate(pairs):
            th_l, th_r = TIERS[k]
            if mask >> k & 1:
                a_th, b_th = th_r, th_l
            else:
                a_th, b_th = th_l, th_r
            assign += [(a, a_th), (b, b_th)]
            for vv, th in ((a, a_th), (b, b_th)):
                rc = 100 + vv["dims"]["c1"] * 1.4
                mom += vv["mass"] * rc * math.cos(math.radians(th))
        if best is None or abs(mom) < abs(best[0]):
            best = (mom, assign)
    assign = best[1]

    def samples(vv, inner_r=0.0):
        """Axis stations and radii of a voice incl. its head (as a disc of radius ~xh)."""
        o = vv["outer"]
        zz = np.linspace(0, o[-1, 1], 60)
        rr = np.interp(zz, o[:, 1], o[:, 0])
        hs = head_geometry(vv)
        zh = o[-1, 1] + np.linspace(0, hs["z_bt"], 4)
        rh = np.full(4, math.hypot(hs["xh"], max(-hs["y_front"], hs["y_rear"])))
        return np.concatenate([zz, zh]), np.concatenate([rr, rh])

    def clear(Rf):
        balls = []                      # (centre, radius) spheres approximating each voice
        for vv, th in assign:
            z_l, _ = frame_of(th)
            zz, rr = samples(vv)
            for z, r in zip(zz, rr):
                balls.append((np.array(z_l) * (Rf + z), r, vv["name"]))
        # Earth: its axis is vertical at y = yE, foot at z = -c1
        e = vs[0]
        yE = HUB_Y0 - (e["dims"]["R"] + L_BOSS)
        zz, rr = samples(e)
        for z, r in zip(zz, rr):
            balls.append((np.array([0.0, yE, z - e["dims"]["c1"]]), r, "Earth"))
        P = np.array([b[0] for b in balls])
        Rr = np.array([b[1] for b in balls])
        names = np.array([b[2] for b in balls])
        dist = np.linalg.norm(P[:, None] - P[None], axis=2) - Rr[:, None] - Rr[None]
        diff = names[:, None] != names[None]
        return dist[diff].min()

    Rf = 60.0
    while clear(Rf) < GAP:
        Rf += 1.0
    Rf = max(math.ceil(Rf / 5.0) * 5.0, RF_MIN)
    return assign, Rf


# ---------------------------------------------------------------- the wheel and its fittings


def wheel(assign, Rf, vs):
    """Hub, eight spokes, the Ogdoad ring, foot pads + pegs, crown and pendant lugs.
    Prints flat, rear face (+y) down."""
    y0, y1 = -T_F / 2, T_F / 2
    rring = Rf - 7.0
    w = cyl_y(HUB_R, HUB_Y0, y1)
    ring = cyl_y(rring + RING_W / 2, y0, y1) - cyl_y(rring - RING_W / 2, y0 - 1, y1 + 1)
    w += ring
    for vv, th in assign:
        z_l, x_l = frame_of(th)
        loc = Location(Plane(origin=(0, 0, 0), x_dir=x_l, z_dir=z_l))
        rf = vv["dims"]["rf"]
        spoke = box(-SPOKE_W / 2, SPOKE_W / 2, y0, y1, HUB_R - 2, Rf - 4)
        pad = box(-rf, rf, y0, y1, Rf - 8, Rf)
        # peg flush with the rear face, so it prints on the bed: centred at y = +1
        peg = box(-PEG_X / 2, PEG_X / 2, y1 - PEG_Y, y1, Rf - 1, Rf + PEG_L)
        w += loc * (spoke + pad + peg)
    # crown lug (top) and pendant lug (bottom)
    # (the pendant lug stays short: the Moon's bottom cone flares toward it)
    for th, reach, hole in ((90.0, Rf + 26, Rf + 18), (270.0, Rf + 6, Rf + 0.5)):
        z_l, x_l = frame_of(th)
        loc = Location(Plane(origin=(0, 0, 0), x_dir=x_l, z_dir=z_l))
        lug = box(-LUG_W / 2, LUG_W / 2, y0, y1, rring, reach)
        lug -= cyl_y(PIN_D / 2 if th == 90 else 2.2, y0 - 1, y1 + 1, 0, hole)
        w += loc * lug
    # Earth's peg (forward) and the boom socket (rear)
    w += box(-EPEG / 2, EPEG / 2, HUB_Y0 - EPEG_L, HUB_Y0, -EPEG / 2, EPEG / 2)
    s = BOOM + 0.4
    w -= box(-s / 2, s / 2, y1 - BOOM_IN, y1 + 1, -s / 2, s / 2)
    return w


def crown(Rf):
    """A comb that straddles the top lug: one M3 bolt through both, and a row of hanging
    holes along y so the swivel can be put exactly above the centre of mass."""
    zl = Rf + 18.0                  # bolt height
    sl = T_F + 0.4
    xt = 8.0
    top = zl + 30.0
    arch = arc_pts(0, top - 70, 72, 54, 126, 24)          # right -> left over the top
    prof = [(arch[-1][0], arch[-1][1]), (-30, zl - 6), (-sl / 2 - 6, zl - 14),
            (sl / 2 + 6, zl - 14), (30, zl - 6)] + arch[:-1]
    plate = yz_prism(prof, -xt / 2, xt / 2)
    plate -= box(-xt, xt, -sl / 2, sl / 2, zl - 15, zl + 8.6)           # slot over the lug
    plate -= cyl_y(PIN_D / 2, -50, 50, 0, zl)                           # bolt
    holes = []
    for k, yy in enumerate(np.arange(-30, 30.1, 6.0)):
        # the arch top: z of the arc at this y, minus 6 mm
        zz = (top - 70) + math.sqrt(72 ** 2 - yy ** 2) - 6.5
        plate -= Pos(0, yy, zz) * Rot(0, 90, 0) * Cylinder(2.1, 2 * xt)
        holes.append((yy, zz))
    return plate, holes


def plumb_bob():
    """Hollow plumb-bob ballast, printed point-up on its flat top. Fill with dry sand
    (~0.3 kg) or steel shot (~0.9 kg) and close with the stopper."""
    R, H, t, tip = 40.0, 46.0, 2.4, 3.0
    s2 = math.sqrt(2)
    outer = [(0, 0), (R, 0), (R, H), (tip, H + R - tip), (0, H + R - tip)]
    zi = H - (s2 - 1) * t                      # inner 45-deg cone, offset t from the outer
    inner = [(0, t), (R - t, t), (R - t, zi), (tip, H + R - tip - s2 * t + (s2 - 1) * 0),
             (0, H + R - tip - s2 * t)]
    body = rev_rz(outer) - rev_rz(inner)
    body -= cyl_z(8.2, -1, t + 1)                                      # fill port
    return body


def stopper():
    """Press-fit plug for the bob's fill port, with the hanging eye."""
    plug = cyl_z(8.05, 0, 9) + cyl_z(11, 9, 12)
    eye = Pos(0, 0, 12) * (box(-3, 3, -9, 9, 0, 18) + Pos(0, 0, 18) * Rot(0, 90, 0)
                           * Cylinder(9, 6))
    eye -= Pos(0, 0, 30) * Rot(0, 90, 0) * Cylinder(3.6, 10)
    return plug + eye


def boom():
    yA = T_F / 2 - BOOM_IN
    yB = T_F / 2 + BOOM_FREE + SLEEVE
    return box(-BOOM / 2, BOOM / 2, yA, yB, -BOOM / 2, BOOM / 2)


def fin():
    """Dove-tail vane: a fan of seven feathers on a sleeve that slides over the boom end.
    The pigeon carries its whistle on its tail; here the tail steers the whistles."""
    yS = T_F / 2 + BOOM_FREE                 # sleeve front
    so, si = BOOM / 2 + 2.0, (BOOM + 0.4) / 2
    sleeve = box(-so, so, yS, yS + SLEEVE, -so, so)
    bore = box(-si, si, yS - 1, yS + SLEEVE + 1, -si, si)
    ax, ay = yS - 20.0, 0.0                  # fan apex (y, z)
    half = 32.0
    n = 7
    edge = []
    angs = np.linspace(-half, half, n + 1)
    for j in range(n):
        a0, a1 = angs[j], angs[j + 1]
        am = 0.5 * (a0 + a1)
        edge.append((ax + 186 * math.cos(math.radians(a0)), ay + 186 * math.sin(math.radians(a0))))
        for a, r in ((a0 + 0.3 * (a1 - a0), 200), (am, 206), (a0 + 0.7 * (a1 - a0), 200)):
            edge.append((ax + r * math.cos(math.radians(a)), ay + r * math.sin(math.radians(a))))
    edge.append((ax + 186 * math.cos(math.radians(half)), ay + 186 * math.sin(math.radians(half))))
    prof = [(yS, -so)] + edge + [(yS, so)]
    x0 = -so
    plate = yz_prism(prof, x0, x0 + FIN_T)
    # feather shafts, raised on the sleeve side
    ribs = None
    for j in range(n):
        am = math.radians(0.5 * (angs[j] + angs[j + 1]))
        r0, r1 = 46.0, 196.0
        p0 = (ax + r0 * math.cos(am), ay + r0 * math.sin(am))
        p1 = (ax + r1 * math.cos(am), ay + r1 * math.sin(am))
        nrm = (-math.sin(am) * 0.7, math.cos(am) * 0.7)
        q = [(p0[0] + nrm[0], p0[1] + nrm[1]), (p1[0], p1[1]), (p0[0] - nrm[0], p0[1] - nrm[1])]
        r = yz_prism(q, x0 + FIN_T - 0.01, x0 + FIN_T + 1.2)
        ribs = r if ribs is None else ribs + r
    return plate + ribs + sleeve - bore


# ---------------------------------------------------------------- assembly


def place_voice(part, v, th, Rf):
    """Local gourd frame -> assembly."""
    if v["name"] == "Earth":
        yE = HUB_Y0 - (v["dims"]["R"] + L_BOSS)
        return Pos(0, yE, -v["dims"]["c1"]) * part
    z_l, x_l = frame_of(th)
    org = tuple(Rf * c for c in z_l)
    return Location(Plane(origin=org, x_dir=x_l, z_dir=z_l)) * part


def build(vs):
    assign, Rf = solve_layout(vs)
    print(f"wheel radius (foot pads) R_f = {Rf:.0f} mm")
    parts = {}
    meas = []
    for v in vs:
        g = gourd(v, socket=v["name"] != "Earth")
        if v["name"] == "Earth":
            g += earth_boss(v)
        h, hs = head(v)
        cav = cavity(v, hs)
        Vm = cav.volume
        # the pitch from the measured volume, scaling the model (lower bulb share kept)
        meas.append((v["name"], v["freq"], v["V"], Vm, None))
        parts[v["name"]] = dict(gourd=g, head=h, hs=hs, v=v)
    th_of = {vv["name"]: th for vv, th in assign}
    th_of["Earth"] = None
    parts["_wheel"] = wheel(assign, Rf, vs)
    parts["_crown"], holes = crown(Rf)
    parts["_bob"] = plumb_bob()
    parts["_stopper"] = stopper()
    parts["_boom"] = boom()
    parts["_fin"] = fin()
    return parts, assign, Rf, th_of, meas, holes


def assembled(parts, th_of, Rf):
    """Every part in assembly position (bob/stopper hung 170 mm under the pendant lug)."""
    out = {}
    for name, th in th_of.items():
        p = parts[name]
        v = p["v"]
        out[name + " gourd"] = place_voice(p["gourd"], v, th, Rf)
        out[name + " head"] = place_voice(Pos(0, 0, v["dims"]["z_top"]) * p["head"], v, th, Rf)
    out["wheel"] = parts["_wheel"]
    out["crown"] = parts["_crown"]
    out["boom"] = parts["_boom"]
    out["fin"] = parts["_fin"]
    zb = -(Rf + 0.5) - 190.0
    out["bob"] = Pos(0, 0, zb) * Rot(180, 0, 0) * parts["_bob"]
    out["stopper"] = Pos(0, 0, zb) * Rot(180, 0, 0) * parts["_stopper"]
    return out


def mass_props(asm):
    """Mass (g) and centre of mass of the printed assembly (ballast excluded)."""
    m_tot, mom = 0.0, np.zeros(3)
    for k, p in asm.items():
        if k in ("bob", "stopper"):
            continue
        m = p.volume / 1000 * PETG
        c = p.center()
        m_tot += m
        mom += m * np.array([c.X, c.Y, c.Z])
    return m_tot, mom / m_tot


def clearance_check(asm):
    """Distance between parts that should not touch. Joints that are meant to mate are
    skipped (gourd-head, gourd-wheel, Earth-wheel, boom-wheel, boom-fin, crown-wheel)."""
    keys = list(asm)
    mates = set()
    for k in keys:
        if k.endswith(" gourd"):
            n = k[:-6]
            mates |= {(k, n + " head"), (k, "wheel")}
    mates |= {("boom", "wheel"), ("boom", "fin"), ("crown", "wheel"), ("bob", "stopper")}
    worst = (1e9, None)
    ok = True
    bbs = {k: asm[k].bounding_box() for k in keys}
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if (a, b) in mates or (b, a) in mates:
                continue
            A, B = bbs[a], bbs[b]
            sep = max(A.min.X - B.max.X, B.min.X - A.max.X, A.min.Y - B.max.Y,
                      B.min.Y - A.max.Y, A.min.Z - B.max.Z, B.min.Z - A.max.Z)
            if sep > 4.0:
                continue
            dist = asm[a].distance_to(asm[b])
            if dist < worst[0]:
                worst = (dist, (a, b))
            if dist < 3.0:
                print(f"  CLASH {a} / {b}: {dist:.2f} mm")
                ok = False
    print(f"  closest non-mating parts: {worst[1]} at {worst[0]:.1f} mm")
    return ok


def mating_check(parts, th_of, Rf):
    """Mating parts must not interpenetrate (beyond float noise)."""
    ok = True
    asm = assembled(parts, th_of, Rf)
    for name in th_of:
        g, h = asm[name + " gourd"], asm[name + " head"]
        iv = (g & h).volume
        iw = (g & asm["wheel"]).volume
        if iv > 0.5 or iw > 0.5:
            print(f"  INTERFERENCE {name}: head {iv:.2f} mm3, wheel {iw:.2f} mm3")
            ok = False
    for a, b in (("boom", "wheel"), ("boom", "fin"), ("crown", "wheel")):
        iv = (asm[a] & asm[b]).volume
        if iv > 0.5:
            print(f"  INTERFERENCE {a}/{b}: {iv:.2f} mm3")
            ok = False
    return ok


# ---------------------------------------------------------------- export


def to_bed(part, rot=None):
    if rot is not None:
        part = rot * part
    part = part.clean()
    bb = part.bounding_box()
    return Pos(-(bb.min.X + bb.max.X) / 2, -(bb.min.Y + bb.max.Y) / 2, -bb.min.Z) * part


def export(part, name, rot=None):
    path = os.path.join(HERE, "stl", name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    export_stl(to_bed(part, rot), path, tolerance=0.01, angular_tolerance=0.12)
    m = trimesh.load(path)
    m.vertices = m.vertices + 0.0
    m.export(path)
    m = trimesh.load(path)
    bodies = len(m.split(only_watertight=False))
    ext = m.extents
    ok = m.is_watertight and m.volume > 0 and bodies == 1 and max(ext[:2]) <= BED and ext[2] <= BED
    print(f"  {'OK ' if ok else 'BAD'} {name:26s} {ext[0]:6.1f} x {ext[1]:6.1f} x {ext[2]:6.1f} mm  "
          f"watertight={m.is_watertight} bodies={bodies} {m.volume / 1000 * PETG:6.1f} g")
    return ok


def export_all(parts):
    ok = True
    side = Rot(0, -90, 0)                   # head: local +x (the D-flat) up, -x on the bed
    for i, (f, name, _) in enumerate(vc.VOICES):
        p = parts[name]
        tag = f"{i}_{name.lower()}_{f}hz"
        ok &= export(p["gourd"], f"{tag}_gourd.stl")
        ok &= export(p["head"], f"{tag}_head.stl", side)
    ok &= export(parts["_wheel"], "wheel.stl", Rot(-90, 0, 0))
    ok &= export(parts["_crown"], "crown.stl", Rot(0, -90, 0))
    ok &= export(parts["_bob"], "plumb_bob.stl")
    ok &= export(parts["_stopper"], "bob_stopper.stl")
    ok &= export(parts["_boom"], "vane_boom.stl")
    ok &= export(parts["_fin"], "vane_fin.stl", Rot(0, -90, 0))
    return ok


# ---------------------------------------------------------------- previews


def tess(part, tol=0.15):
    import pyvista as pv
    v, t = part.tessellate(tol, 0.25)
    v = np.array([[p.X, p.Y, p.Z] for p in v])
    f = np.hstack([np.full((len(t), 1), 3), np.asarray(t)]).ravel()
    return pv.PolyData(v, f)


GOURD_COL = "#d9c7a3"
HEAD_COL = "#3d6fb6"
FRAME_COL = "#b08a3c"


def previews(asm, parts):
    import pyvista as pv
    pv.OFF_SCREEN = True
    meshes = {k: tess(p) for k, p in asm.items()}

    def col(k):
        if k.endswith("gourd"):
            return GOURD_COL
        if k.endswith("head"):
            return HEAD_COL
        if k in ("bob", "stopper"):
            return "#7d7d78"
        return FRAME_COL

    def scene(out, cam, size=(1400, 1400), zoom=1.0, bg="#f2f2f0"):
        pl = pv.Plotter(off_screen=True, window_size=size)
        pl.set_background(bg)
        for k, m in meshes.items():
            pl.add_mesh(m, color=col(k), smooth_shading=True, specular=0.25)
        pl.camera_position = cam
        pl.camera.zoom(zoom)
        pl.screenshot(os.path.join(HERE, out))
        pl.close()
        print("  wrote", out)

    scene("preview_front.png", [(0, -1700, -60), (0, 0, -60), (0, 0, 1)], zoom=1.1)
    scene("preview_three_quarter.png", [(-900, -1150, 420), (0, 60, -10), (0, 0, 1)], zoom=1.2)
    scene("preview_side.png", [(1500, 80, 30), (0, 120, -20), (0, 0, 1)], zoom=1.15)

    head_section(parts["Earth"])


def head_section(p):
    """2D section of a whistle head on its mid-plane (x = 0), drawn from the same polygons
    the solids are cut from: the jet's whole path in one picture."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MPoly
    v = p["v"]
    hs = p["hs"]
    d = v["dims"]
    R, Rn = vc.R_LIP, vc.R_LIP + vc.T_BLK
    yf, yw0, y0, yL, yr = hs["y_front"], hs["y_ww0"], hs["y0"], hs["yL"], hs["y_rear"]
    z_f, z_r, z_bt = hs["z_f"], hs["z_r"], hs["z_bt"]
    ext = ([(yf, 0.0), (yr, 0.0), (yr, hs["z_ramp"]), (yL, hs["z_tip"]), (yL, z_f), (y0, z_f),
            (y0, z_bt)] + arc_pts(yw0, z_r + R, Rn, 90, 180, 16))
    skirt = [(-hs["ys"], -vc.SPIG_H), (hs["ys"], -vc.SPIG_H), (hs["ys"], z_f), (-hs["ys"], z_f)]
    ww = ([(yf - 1.0, z_f - R)] + arc_pts(yw0, z_f - R, R, 180, 90, 16)
          + [(y0 + 0.3, z_f), (y0 + 0.3, z_r)] + arc_pts(yw0, z_r + R, R, 270, 180, 16)
          + [(yf - 1.0, z_r + R)])
    pl = vc.plenum_poly(hs)
    rs = hs["rso"] + vc.FIT
    sock = [(-rs, -vc.SPIG_H - 0.01), (rs, -vc.SPIG_H - 0.01), (rs, 0), (-rs, 0)]
    neck = [(-d["rno"], -vc.SPIG_H - 6), (d["rno"], -vc.SPIG_H - 6), (d["rno"], 0), (-d["rno"], 0)]
    bore = [(-d["rb"], -vc.SPIG_H - 7), (d["rb"], -vc.SPIG_H - 7), (d["rb"], 0.01), (-d["rb"], 0.01)]
    fig, ax = plt.subplots(figsize=(12, 6.2), dpi=130)
    blue, sand, air = "#3d6fb6", "#d9c7a3", "#f2f2f0"
    for poly, col in ((skirt, blue), (ext, blue), (sock, air), (neck, sand), (bore, air),
                      (ww, air), (pl, air)):
        ax.add_patch(MPoly(poly, closed=True, fc=col, ec="none"))
    ax.annotate("", xy=(y0 - 0.2, (z_f + z_r) / 2), xytext=(yf - 6, (z_f + z_r) / 2),
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.6))
    ax.text(yf - 6, z_r + 1.0, "wind", color="#c0392b", fontsize=11)
    notes = [((y0 + yL) / 2, z_f + 3.6, f"window {vc.W_CUT} mm"),
             (yL + 4.5, hs["z_tip"] + 0.2, "labium (30 deg)"),
             ((yw0 + y0) / 2, z_r + 0.45, f"windway {vc.H_WW} mm"),
             (0.0, 0.9, "plenum"), (0.0, -6.5, "neck bore to the cavity")]
    for x, z, t in notes:
        ax.text(x, z, t, ha="center", fontsize=10, color="#222")
    ax.set_xlim(yf - 8, yr + 2)
    ax.set_ylim(-vc.SPIG_H - 7, z_bt + 5)
    ax.set_aspect("equal")
    ax.set_xlabel("y, along the wind (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title(f"Whistle head, section on the window's mid-plane ({v['name']}, {v['freq']} Hz)")
    ax.set_facecolor(air)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "preview_head_section.png"))
    plt.close(fig)
    print("  wrote preview_head_section.png")


# ---------------------------------------------------------------- main


def main():
    vs = vc.design()
    parts, assign, Rf, th_of, meas, holes = build(vs)
    ok = True
    # This confirms the CAD builds the cavity the solver sized; it is not a pitch check.
    # The pitch itself is checked independently by fem_check.py (3-D acoustic FEM).
    print("\nCavity volumes, design vs built solid (volume agreement only; pitch: fem_check.py):")
    for name, f, Vd, Vm, _ in meas:
        err = Vm / Vd - 1
        flag = abs(err) < 1e-3
        ok &= flag
        print(f"  {'OK ' if flag else 'BAD'} {name:8s} {f:4d} Hz  V {Vd / 1000:8.2f} / {Vm / 1000:8.2f} cm3"
              f"  ({100 * err:+.3f} %)")
    print("\nLayout (tier: left / right, the climb in order):")
    for vv, th in assign:
        side = "left " if math.cos(math.radians(th)) < 0 else "right"
        print(f"  {vv['name']:8s} {vv['freq']:4d} Hz at {th:6.1f} deg ({side})")

    asm = assembled(parts, th_of, Rf)
    print("\nMating parts:")
    ok &= mating_check(parts, th_of, Rf)
    print("Clearances:")
    ok &= clearance_check(asm)
    m, cg = mass_props(asm)
    print(f"\nPrinted mass {m:.0f} g, centre of mass ({cg[0]:+.1f}, {cg[1]:+.1f}, {cg[2]:+.1f}) mm")
    best = min(holes, key=lambda h: abs(h[0] - cg[1]))
    print(f"  hang from the crown hole at y = {best[0]:+.0f} mm "
          f"(nearest to the centre of mass; +y is toward the vane)")
    bbm = None
    for k, p in asm.items():
        bb = p.bounding_box()
        bbm = (np.array([bb.min.X, bb.min.Y, bb.min.Z]), np.array([bb.max.X, bb.max.Y, bb.max.Z])) \
            if bbm is None else (np.minimum(bbm[0], [bb.min.X, bb.min.Y, bb.min.Z]),
                                 np.maximum(bbm[1], [bb.max.X, bb.max.Y, bb.max.Z]))
    span = bbm[1] - bbm[0]
    print(f"  assembled envelope {span[0]:.0f} (w) x {span[1]:.0f} (d) x {span[2]:.0f} (h) mm "
          f"incl. the bob on a 190 mm cord")

    layout = dict(
        R_f=Rf, hub_y0=HUB_Y0, frame_t=T_F, ring_r=Rf - 7.0, boom_free=BOOM_FREE,
        earth_y=HUB_Y0 - (vs[0]["dims"]["R"] + L_BOSS), earth_z0=-vs[0]["dims"]["c1"],
        angles={vv["name"]: th for vv, th in assign}, climb=[vv["name"] for vv, _ in assign],
        mass_g=round(m, 0), cg=[round(x, 1) for x in cg], hang_hole_y=best[0],
        envelope=[round(x) for x in span],
    )
    vc.export_json(vs, layout)

    print("\nExporting STLs (print orientation):")
    ok &= export_all(parts)
    if "--no-preview" not in sys.argv:
        print("\nRendering previews:")
        previews(asm, parts)
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
