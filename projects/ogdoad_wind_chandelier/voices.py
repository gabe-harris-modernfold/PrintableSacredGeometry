#!/usr/bin/env python3
"""The Hymn of the Ogdoad -- acoustic design of the nine wind-blown voices.

Each voice is a pigeon-whistle (gelin) style vessel flute: a printed gourd is the
Helmholtz cavity, and a whistle head on its neck holds the window. Wind enters a
bellmouth intake, leaves a thin windway as a jet, crosses the window and splits on a
sharp labium. The jet's flapping pumps air in and out of the window at the cavity's
resonance.

This module is the single source of truth for the CAD (chandelier.py) and the web page
(voices.json). It has no CAD dependency: the gourd profiles and every volume are
computed here exactly (frustum sums, polygon areas), and chandelier.py re-measures
the built solids against them.

Model
-----
* Helmholtz with a coupled waist. A bottle gourd is two bulbs joined by a waist, which
  is itself a short neck. The voice is solved as the 2-DOF chain
      outside --M_mouth-- V_upper --M_waist-- V_lower
  and the lower root is the voice. The waist is kept wide (area >> window), so the
  second mode sits far above and the lumped picture holds.
* Window inertance: M = rho * L' / A with A = W*b and L' = t_p + 2*delta, where delta is
  the low-frequency radiation-mass end correction of a uniform-velocity W x b piston in
  a baffle, delta = J / (2 pi A), J being the rectangle's self-potential integral. For a
  circle this gives 8a/(3 pi) = 0.849a; a slit's correction is much smaller than the
  equal-area circle's, which is why it is computed rather than assumed.
* Jet drive (wind bands). Jet speed U = ETA * wind. A jet across a window of length W
  drives the cavity while the jet Strouhal number St = f W / U lies in a band: too slow
  (St > ST_ON) and the jet's travel delay is too long for gain; optimum near ST_OPT
  (jet-delay criterion, Fletcher & Rossing; recorder and ocarina playing data put
  St ~ 0.15-0.25); too fast (St < ST_OFF) and it overblows into hiss. Because every
  voice shares one cut-up W, each voice's band scales with its pitch: light air wakes
  Earth, and the Ogdoad only sings in a fresh breeze. These band edges are a model and
  are meant to be checked with a fan and an anemometer.

Run:  python voices.py          (prints the table, writes voices.json)
"""

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- air and drive

T_AIR_C = 20.0
C_AIR = 331.3 * math.sqrt(1 + T_AIR_C / 273.15)        # 343.2 m/s at 20 C
RHO = 1.204                                             # kg/m^3 at 20 C
ETA = 0.85              # jet speed / wind speed (stagnation-fed intake, minus losses)
ST_ON, ST_OPT, ST_OFF = 0.40, 0.22, 0.10                # Strouhal band of the jet drive
# Window end-correction factor, K_END(b) = K_END0 + K_END1 * ln(b / 10 mm), fitted to the 3-D
# FEM of the built heads (fem_check.py). It is < 1 because the uniform-velocity piston formula
# overstates a slit's radiation mass and the inner side opens into a confined plenum.
K_END0, K_END1 = 0.960, 0.0375

# ---------------------------------------------------------------- the nine voices

# (Hz, name, sphere) -- geocentric order: Earth at the centre, then the seven planetary
# spheres in the Chaldean order Poimandres climbs, then the eighth sphere, the Ogdoad.
VOICES = [
    (174, "Earth", "centre"),
    (285, "Moon", "1st sphere"),
    (396, "Mercury", "2nd sphere"),
    (417, "Venus", "3rd sphere"),
    (528, "Sun", "4th sphere"),
    (639, "Mars", "5th sphere"),
    (741, "Jupiter", "6th sphere"),
    (852, "Saturn", "7th sphere"),
    (963, "Ogdoad", "8th sphere (fixed stars)"),
]
SOLFEGE = {396: "UT", 417: "RE", 528: "MI", 639: "FA", 741: "SOL", 852: "LA", 963: "SI"}

# ---------------------------------------------------------------- whistle head (mm)

W_CUT = 2.4             # window length along the jet (the "cut-up"), shared by all voices
H_WW = 0.70             # windway height (jet thickness); W/h = 3.4
T_P = 1.2               # windway-floor slab = the window's neck length
L_WW = 8.0              # windway length (lets the jet settle into a sheet)
R_LIP = 2.5             # bellmouth lip radii (upper and lower)
T_BLK = 1.4             # block material above the upper lip
PLEN_H = 2.0            # plenum height between the spigot top and the slab
PLEN_FRONT = 2.0        # plenum reaches this far ahead of the window, under the slab
PLEN_REAR = 7.0         # ... and this far behind the labium tip
LAB_DOWN = 15.0         # labium underside angle below the jet axis (deg)
LAB_RAMP = 15.0         # labium ramp angle above the jet axis (deg); wedge = 30 deg
LAB_TIP = 0.35          # labium tip height above the window plane, in windway heights
L_RAMP = 8.0            # ramp length behind the tip
SPIG_H = 8.0            # neck spigot that the head's socket slides over
NECK_FREE = 2.0         # exposed neck between the gourd shoulder and the head skirt
FIT = 0.15              # radial glue clearance, socket over spigot
T_SKIRT = 1.6           # head socket wall
T_SIDE = 1.8            # head side walls beside the windway
D_FLAT = 0.72           # D-flat at local +x, at this fraction of the spigot radius

# ---------------------------------------------------------------- gourd body

BETA = 0.26            # window width b as a fraction of the gourd's bulb diameter
RHO2_K = 0.72           # upper bulb radius / lower bulb radius
WAIST_V = 0.42          # waist V-vertex radius / lower bulb radius
WAIST_F = 0.12          # waist fillet radius / lower bulb radius
FOOT_K, FOOT_MIN = 0.42, 9.5     # foot radius = max(FOOT_K * R, FOOT_MIN)
FLOOR_T = 10.0          # solid foot: holds the 8 mm blind socket for the frame peg
SOCK_X, SOCK_Y, SOCK_D = 7.4, 8.4, 8.0     # foot socket (frame peg 7 x 8 + clearance)
SOCK_OFF = 1.0          # socket centre offset toward local +y: keys the orientation
T_NECK = 1.6            # neck wall


def wall_t(R):
    """Shell thickness grows with size: stiff walls keep the cavity's Q up."""
    return float(np.clip(0.032 * R + 0.8, 1.6, 2.8))


def bore_r(b):
    return max(b / 2.0, 4.0) + 0.8


def arc(cx, cz, r, a0, a1, n):
    """Points on a circle; angles in degrees measured from +r toward +z."""
    a = np.radians(np.linspace(a0, a1, n))
    return np.stack([cx + r * np.cos(a), cz + r * np.sin(a)], axis=1)


def gourd_dims(R, b):
    """Every station of the bottle-gourd profile, from the lower bulb radius R (mm)."""
    rho2 = RHO2_K * R
    rf = max(FOOT_K * R, FOOT_MIN)
    c1 = R * math.sqrt(2) - rf                  # bottom 45-deg cone meets the foot at z=0
    rv = WAIST_V * R
    zv = c1 + R * math.sqrt(2) - rv             # waist V vertex
    c2 = zv + rho2 * math.sqrt(2) - rv
    rb = bore_r(b)
    rno = rb + T_NECK
    z_sh = c2 + rho2 * math.sqrt(2) - rno       # top 45-deg cone reaches the neck
    z_top = z_sh + NECK_FREE + SPIG_H
    return dict(R=R, rho2=rho2, rf=rf, c1=c1, rv=rv, zv=zv, c2=c2, rwf=WAIST_F * R,
                rb=rb, rno=rno, z_sh=z_sh, z_top=z_top, t=wall_t(R))


def gourd_profile(d, inner=False, n=28):
    """(r, z) polyline from the foot (or floor) up the side to the neck top.

    Lower bulb: 45-deg cone, sphere arc between polar angles 135 and 45 deg, 45-deg cone.
    Waist: concave fillet in the 90-deg V of the two cones. Upper bulb likewise, then a
    45-deg cone into a vertical neck. Every surface stays within 45 deg of vertical, so
    the gourd prints upright with no supports. The inner surface is the same construction
    with each sphere shrunk by the wall t (offsetting a 45-deg cone tangent to a sphere
    gives the cone tangent to the offset sphere), the same fillet centre with radius
    r + t, the floor raised to FLOOR_T and the neck bore at rb."""
    t = d["t"] if inner else 0.0
    R, rho2 = d["R"] - t, d["rho2"] - t
    c1, c2, zv = d["c1"], d["c2"], d["zv"]
    s2 = math.sqrt(2)
    rv = d["rv"] - t * s2
    rwf = d["rwf"] + t
    pts = []
    if inner:
        pts.append((0.0, c1 - R * s2))           # bottom cone apex; clipped at the floor below
    else:
        pts.append((d["rf"], 0.0))
    # lower sphere arc: polar 135 -> 45 deg, i.e. angle from +r of -45 -> +45 deg
    pts += [tuple(p) for p in arc(0.0, c1, R, -45, 45, n)]
    # waist fillet: centre (rv + rwf*sqrt2, zv); tangent to both 45-deg lines
    cx = rv + rwf * s2
    pts += [tuple(p) for p in arc(cx, zv, rwf, 225, 135, max(8, n // 3))]
    pts += [tuple(p) for p in arc(0.0, c2, rho2, -45, 45, n)]
    rn = d["rb"] if inner else d["rno"]
    zs = c2 + rho2 * s2 - rn                     # top cone meets the neck radius
    pts.append((rn, zs))
    pts.append((rn, d["z_top"]))
    pts = np.asarray(pts)
    if inner:
        pts = split_at(pts, FLOOR_T)[1]
    return pts


def check_profile(d):
    """The construction is valid only if each cone segment is long enough for its fillet
    and the inner surface stays inside the outer one."""
    R, rho2, rv, rwf = d["R"], d["rho2"], d["rv"], d["rwf"]
    assert rho2 - rv * math.sqrt(2) >= rwf - 1e-9, "upper cone too short for the waist fillet"
    assert R - rv * math.sqrt(2) >= rwf - 1e-9, "lower cone too short for the waist fillet"
    assert d["rf"] < R * math.sin(math.pi / 4), "foot wider than the bottom cone"
    ri = d["rv"] - d["t"] * math.sqrt(2) + (math.sqrt(2) - 1) * (d["rwf"] + d["t"])
    # the line model counts the waist's mass, so it may pinch below the bore; keep it open
    # enough to print cleanly and to stay a minor constriction
    assert ri >= max(0.75 * d["rb"], 3.5), "waist too narrow"


def revolved_volume(pts):
    """Volume enclosed by revolving a monotone-z (r, z) polyline about z (frustum sums)."""
    r, z = pts[:, 0], pts[:, 1]
    dz = np.diff(z)
    return float(np.sum(math.pi * dz * (r[:-1] ** 2 + r[:-1] * r[1:] + r[1:] ** 2) / 3.0))


def split_at(pts, zc):
    """Split a monotone-z polyline at height zc, interpolating the cut point."""
    i = int(np.searchsorted(pts[:, 1], zc))
    p0, p1 = pts[i - 1], pts[i]
    f = (zc - p0[1]) / (p1[1] - p0[1])
    pc = p0 + f * (p1 - p0)
    return np.vstack([pts[:i], pc]), np.vstack([pc, pts[i:]])


# ---------------------------------------------------------------- whistle head geometry

def head_stations(b, rso):
    """Key heights/positions of the head, local frame: z up the gourd axis (z=0 at the
    spigot top), +y along the jet (downwind), x across the window."""
    z_pc = PLEN_H                       # plenum ceiling (slab underside)
    z_f = z_pc + T_P                    # window plane = windway floor
    z_r = z_f + H_WW                    # windway roof
    z_tip = z_f + LAB_TIP * H_WW
    y0 = -W_CUT / 2.0                   # windway exit / window front edge
    yL = W_CUT / 2.0                    # labium tip
    y_ww0 = y0 - L_WW                   # windway entry
    y_front = y_ww0 - R_LIP - T_BLK     # nose: lips are quarter-annuli of thickness T_BLK
    z_bt = z_r + 2 * R_LIP + T_BLK      # block top (its rear face shelters the window)
    y_rear = yL + L_RAMP
    z_ramp = z_tip + L_RAMP * math.tan(math.radians(LAB_RAMP))
    # skirt half-length in y: around the socket, and past the plenum so it is closed below
    ys = max(rso + T_SKIRT, W_CUT / 2 + PLEN_REAR + 1.2, W_CUT / 2 + PLEN_FRONT + 1.2)
    xh = max(b / 2.0 + T_SIDE, rso + T_SKIRT)
    return dict(z_pc=z_pc, z_f=z_f, z_r=z_r, z_tip=z_tip, y0=y0, yL=yL, y_ww0=y_ww0,
                y_front=y_front, z_bt=z_bt, y_rear=y_rear, z_ramp=z_ramp, ys=ys, xh=xh,
                rso=rso)


def plenum_poly(hs):
    """yz section of the plenum + window void (extruded over x in [-b/2, b/2]).
    Its top edge runs from the windway exit, over the window, to the labium tip and down
    the labium underside; it is clipped at z=0, where it meets the neck bore."""
    tdown = math.tan(math.radians(LAB_DOWN))
    y_pf = hs["y0"] - PLEN_FRONT
    y_pr = hs["yL"] + PLEN_REAR
    z_under = hs["z_tip"] - (y_pr - hs["yL"]) * tdown     # underside height at the rear wall
    z_under = max(z_under, 0.6)
    return np.array([
        (y_pf, 0.0), (y_pr, 0.0), (y_pr, z_under), (hs["yL"], hs["z_tip"]),
        (hs["yL"], hs["z_f"]), (hs["y0"], hs["z_f"]), (hs["y0"], hs["z_pc"]),
        (y_pf, hs["z_pc"]),
    ])


def poly_area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


# ---------------------------------------------------------------- acoustics

def rect_self_potential(a, b, n=400):
    """J = integral over a W x b rectangle twice of 1/|r - r'| (mm^3).

    Reduced to 4 * int_0^a int_0^b (a-u)(b-v) / sqrt(u^2+v^2) dv du, which in polar
    coordinates has no singularity; integrated with Gauss-Legendre on two triangles."""
    xg, wg = np.polynomial.legendre.leggauss(n)

    def tri(A, B):                    # triangle with legs A (along u) and B (along v)
        phi_max = math.atan2(B, A)
        ph = 0.5 * phi_max * (xg + 1)
        wp = 0.5 * phi_max * wg
        tot = 0.0
        for p, w in zip(ph, wp):
            rmax = A / math.cos(p)
            rr = 0.5 * rmax * (xg + 1)
            wr = 0.5 * rmax * wg
            u, v = rr * math.cos(p), rr * math.sin(p)
            tot += w * np.sum(wr * (A - u) * (B - v))     # r dr dphi / r = dr dphi
        return tot

    # the integrand is symmetric in (u, a) <-> (v, b); split the rectangle along its diagonal
    return 4.0 * (tri(a, b) + tri(b, a))


def end_correction(W, b):
    """Uniform-velocity baffled-piston end correction of a W x b window (mm)."""
    return rect_self_potential(W, b) / (2 * math.pi * W * b)


def line_admittance(inner, f, n=600):
    """Acoustic admittance looking down into the gourd from the neck top (m^3/(Pa s)).

    The interior is a 1-D acoustic line (plane waves in slices of the revolved profile),
    closed by the rigid floor. Each slice carries both compliance and mass, so the neck
    bore's air counts as mass in series with the window (a lumped model that treats the
    bore as volume alone puts these voices 2.5-4 semitones flat), and the gourd's length
    (k L up to ~0.6 for Earth) is handled without a lumped approximation."""
    z = np.linspace(inner[0, 1], inner[-1, 1], n + 1)
    r = np.interp(z, inner[:, 1], inner[:, 0])
    A = math.pi * (0.5 * (r[:-1] + r[1:])) ** 2 * 1e-6
    dz = np.diff(z) * 1e-3
    k = 2 * math.pi * f / C_AIR
    p, U = 1.0 + 0j, 0j                         # floor: rigid, U = 0
    for Ai, L in zip(A, dz):
        Z0 = RHO * C_AIR / Ai
        cs, sn = math.cos(k * L), math.sin(k * L)
        p, U = cs * p - 1j * Z0 * sn * U, -1j * sn / Z0 * p + cs * U
    return -U / p


def mouth_length(b):
    """Effective neck length of the window (mm): slab + calibrated end corrections.
    The baffled-piston correction (uniform velocity) is scaled by K_END, fitted to the
    3-D FEM of the built heads (fem_check.py), which also absorbs the bore-to-slot
    junction under the plenum."""
    return T_P + k_end(b) * 2 * end_correction(W_CUT, b)


def k_end(b):
    return K_END0 + K_END1 * math.log(b / 10.0)


def line_pitch(inner, b, V_plen, f_guess):
    """Lowest root of  j w M_mouth + 1 / (Y_gourd + j w C_plenum) = 0."""
    A_m = W_CUT * b * 1e-6
    L_m = mouth_length(b) * 1e-3

    def X(f):
        w = 2 * math.pi * f
        Y = line_admittance(inner, f) + 1j * w * V_plen * 1e-9 / (RHO * C_AIR ** 2)
        return (1j * w * RHO * L_m / A_m + 1 / Y).imag

    # X rises through zero at each resonance and jumps from +inf to -inf at each pole (a
    # gourd anti-resonance). Scan upward from well below the note and take the first
    # rising zero: that is the fundamental. (A plain bisection can lock onto a pole.)
    grid = np.geomspace(0.2 * f_guess, 4.0 * f_guess, 160)
    xs = [X(f) for f in grid]
    for i in range(len(grid) - 1):
        if xs[i] < 0 <= xs[i + 1]:
            lo, hi = grid[i], grid[i + 1]
            break
    else:
        raise ValueError("no resonance in the scan window")
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if X(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def voice_acoustics(R, f_target=None):
    """Everything about one voice from its bulb radius R (mm)."""
    b = BETA * 2 * R
    d = gourd_dims(R, b)
    outer = gourd_profile(d)
    inner = gourd_profile(d, inner=True)
    rso = d["rno"]
    hs = head_stations(b, rso)
    low, up = split_at(inner, d["zv"])
    V_low = revolved_volume(low)                 # starts at the flat floor
    V_up = revolved_volume(up)                   # ends at the neck top (spigot bore incl.)
    pl = plenum_poly(hs)
    V_plen = poly_area(pl) * b
    V_tot = V_low + V_up + V_plen
    A_m = W_CUT * b
    delta = end_correction(W_CUT, b)
    L_m = mouth_length(b)
    # single-volume Helmholtz, window only: what the naive formula would say
    f_single = C_AIR / (2 * math.pi) * math.sqrt(A_m * 1e-6 / (V_tot * 1e-9 * L_m * 1e-3))
    f = line_pitch(inner, b, V_plen, f_single)
    rw = d["rv"] - d["t"] * math.sqrt(2) + (math.sqrt(2) - 1) * (d["rwf"] + d["t"])
    # bore + taper inertance per unit density (1/m), for reporting
    zz = np.linspace(d["c2"] + (d["rho2"] - d["t"]) / math.sqrt(2), inner[-1, 1], 300)
    rr = np.interp(zz, inner[:, 1], inner[:, 0])
    L_bore = float(np.trapezoid(1 / (math.pi * rr ** 2), zz)) * 1e3     # (1/mm) -> (1/m)
    return dict(R=R, b=b, dims=d, outer=outer, inner=inner, head=hs, plenum=pl,
                V_low=V_low, V_up=V_up, V_plen=V_plen, V=V_tot, A=A_m, delta=delta,
                L_eff=L_m, r_waist=rw, L_bore=L_bore, f=f, f_single=f_single)


def solve_voice(f_target):
    """Bisect the bulb radius so the voice's fundamental hits f_target."""
    lo, hi = 8.0, 120.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if voice_acoustics(mid, f_target)["f"] > f_target:
            lo = mid
        else:
            hi = mid
    v = voice_acoustics(0.5 * (lo + hi), f_target)
    check_profile(v["dims"])
    return v


def wind_band(f):
    """Wind speeds (m/s): onset, optimum, overblow, for one voice."""
    k = f * W_CUT * 1e-3 / ETA
    return k / ST_ON, k / ST_OPT, k / ST_OFF


def roughness(f1, f2):
    """Plomp-Levelt roughness of two equal pure tones (Sethares' parameterisation)."""
    lo, hi = min(f1, f2), max(f1, f2)
    s = 0.24 / (0.021 * lo + 19.0)
    x = s * (hi - lo)
    return math.exp(-3.5 * x) - math.exp(-5.75 * x)


def shell_mass_g(v, density=1.27):
    """PETG mass of gourd shell + head (g), from exact profile volumes."""
    d = v["dims"]
    V_out = revolved_volume(v["outer"]) + 0.0
    V_in = revolved_volume(v["inner"])
    hs = v["head"]
    head_vol = (hs["y_rear"] - hs["y_front"]) * 2 * hs["xh"] * (hs["z_bt"] + SPIG_H) * 0.55
    return density * (V_out - V_in + head_vol) / 1000.0


def design():
    out = []
    for f, name, sphere in VOICES:
        v = solve_voice(f)
        v.update(freq=f, name=name, sphere=sphere, solfege=SOLFEGE.get(f, ""))
        v["band"] = wind_band(f)
        v["mass"] = shell_mass_g(v)
        out.append(v)
    return out


def table(vs):
    print(f"c = {C_AIR:.1f} m/s at {T_AIR_C:.0f} C; window cut-up W = {W_CUT} mm, windway h = {H_WW} mm")
    print(f"{'voice':8s} {'Hz':>4s} {'R':>6s} {'b':>6s} {'V cm3':>8s} {'d_end':>6s} {'L_eff':>6s} "
          f"{'waist':>6s} {'f line':>7s} {'1-vol':>7s} {'H mm':>6s} {'g':>5s}  wind on/opt/off m/s")
    for v in vs:
        H = v["dims"]["z_top"] + v["head"]["z_bt"]
        on, opt, off = v["band"]
        print(f"{v['name']:8s} {v['freq']:4d} {v['R']:6.2f} {v['b']:6.2f} {v['V'] / 1000:8.2f} "
              f"{v['delta']:6.2f} {v['L_eff']:6.2f} {v['r_waist']:6.2f} {v['f']:7.2f} "
              f"{v['f_single']:7.1f} {H:6.1f} {v['mass']:5.0f}  {on:4.1f} / {opt:4.1f} / {off:4.1f}")


def export_json(vs, layout=None, path=os.path.join(HERE, "voices.json")):
    data = dict(
        c_air=C_AIR, rho=RHO, eta=ETA, st=[ST_ON, ST_OPT, ST_OFF], W=W_CUT, h=H_WW, k_end=[K_END0, K_END1],
        t_p=T_P, temp_c=T_AIR_C,
        voices=[dict(
            freq=v["freq"], name=v["name"], sphere=v["sphere"], solfege=v["solfege"],
            R=round(v["R"], 3), b=round(v["b"], 5), V_cm3=round(v["V"] / 1000, 5),
            V_lower_cm3=round(v["V_low"] / 1000, 5), V_plen_cm3=round(v["V_plen"] / 1000, 6),
            A_mm2=round(v["A"], 5),
            delta_mm=round(v["delta"], 5), L_eff_mm=round(v["L_eff"], 5),
            waist_r_mm=round(v["r_waist"], 4), L_bore_per_m=round(v["L_bore"], 2),
            f_model=round(v["f"], 3),
            f_single_volume=round(v["f_single"], 2),
            wind=[round(x, 2) for x in v["band"]], mass_g=round(v["mass"], 1),
            outer=[[round(r, 2), round(z, 2)] for r, z in v["outer"]],
            inner=[[round(r, 4), round(z, 4)] for r, z in v["inner"]],
            head=dict((k, round(val, 3)) for k, val in v["head"].items()),
        ) for v in vs],
    )
    if layout:
        data["layout"] = layout
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1)
    print("wrote", os.path.basename(path))


if __name__ == "__main__":
    # sanity: the quadrature reproduces the closed form for a rectangle a x b,
    # J = 2/3 (a^3 + b^3 - (a^2+b^2)^1.5) + 2 a^2 b asinh(b/a) + 2 a b^2 asinh(a/b)
    a, bb = 2.4, 20.0
    J_sq = rect_self_potential(a, bb)
    J_exact = (2 / 3 * (a ** 3 + bb ** 3 - (a * a + bb * bb) ** 1.5)
               + 2 * a * a * bb * math.asinh(bb / a) + 2 * a * bb * bb * math.asinh(a / bb))
    print(f"rectangle self-potential: {J_sq:.4f} vs exact {J_exact:.4f}  "
          f"(rel err {abs(J_sq / J_exact - 1):.1e})")
    vs = design()
    table(vs)
    print("\nPure-tone roughness (Plomp-Levelt, 0..1 of the worst case), pairs above 0.05:")
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            r = roughness(vs[i]["freq"], vs[j]["freq"]) / 0.18
            if r > 0.05:
                print(f"  {vs[i]['freq']} + {vs[j]['freq']}: {r:.2f}")
    export_json(vs)
