#!/usr/bin/env python3
"""POC-B VOLUTE -- printable STLs for the cycle-10 KNEE.

ACOUSTICS, taken verbatim from pareto_c10.csv (the KNEE of G vs Gamma):
    r_t 23.9624   r_m 118.8258   path L 871.1085 mm   spiral pitch k 0.10574
    S(x) = pi r_t^2 e^{mx},  m = ln((r_m/r_t)^2)/L
    -> f_c 100 Hz, G 17.6 dB, Gamma 0.212.  None of it is touched below.

LAYOUT, re-derived here. The c10 packing gate compares each winding against the
spacing to the turn OUTSIDE it, `gap = r_c(e^{2 pi k} - 1)`, while the binding
neighbour is the turn INSIDE, whose spacing is smaller by e^{2 pi k} -- and the
half-width that must fit beside it is the outer turn's, not the inner one's. Run
the honest check (`python volute.py --pack`) and the published layout
(r1 45.52, H 95.4->245.5, fold 0.6745) overlaps itself by 14.8 mm.

The fix costs nothing acoustically, because S(x), H(x) and L set the 1-D chain
and the fold only sets where the duct is bent: raise the channel height grading
to 130 -> 290 mm and the spiral start radius to 52 mm. Same f_c, same G, same
fold fraction, same one-segment bell -- and the corrected gate passes with the
turns' walls clearing by 1.4 mm, which is what the gate was always trying to say.

Cross-section: flat floor, vertical sides, 45-degree pointed vault. The vault is
why the body prints in one piece with no support and no lid -- and it holds S(x)
EXACTLY: width solves w^2/4 - wH + S = 0, so the area lost to the vault is given
back by the width. (The doc's micro-perforated lid absorber is moot here: the
flat lids it was designed to damp no longer exist.)

Run from this directory:  python volute.py [--no-skin] [--pack]
"""
import argparse
import math

import numpy as np
import trimesh

import horn_lib as H

OUT = "stl/"

RT, RM, LPATH, KSP = 23.9624, 118.8258, 871.1085, 0.10574      # published c10
H0, H1, R1, FFOLD = 130.0, 290.0, 52.0, 0.6745                 # re-laid out
WALL, PLATE, R_HUB = 5.0, 5.0, 15.0
PLUG_LEN, PLUG_WALL, PLUG_CLR = 18.0, 3.0, 0.50
BITE = 0.6            # how far a plug root bites into its own wall, so the two
                      # solids overlap instead of merely touching along a face
CUP_R, CUP_LEN = 22.5, 40.0                                    # fixed lip cup
TURN_DEG = 62.0                                                # cup swing, see mouthpiece()
MFL = math.log((RM / RT) ** 2) / LPATH
Q = 4                                                          # outline density
SEC = 32 * Q                                                   # points per section

# the Skin: sharp claws along the vault, and the cradle grip carved into the
# first turn's roof (x is path length, mm); the shell surfaces stay clean
GRIP_X = (85.0, 178.0)             # exposed on every side: the second turn only
GRIP_RISE, GRIP_DEPTH = 1.6, 1.6   # shadows x < 63, the bell run x < ~150+
# (x centre, fraction up the right roof face, half-width mm along x, across face)
# -- the roof face is only ~20 mm wide here, so the pits are finger-PAD sized,
# not the Clarion's full finger grooves
FINGERS = ((95.0, 0.52, 9.0, 6.5), (118.0, 0.56, 9.0, 6.5),
           (141.0, 0.60, 9.0, 6.5), (164.0, 0.64, 9.0, 6.5))
THUMB = (126.0, 0.18, 11.0, 8.0)


def S(x):
    return math.pi * RT * RT * np.exp(MFL * np.asarray(x, float))


def Hx(x):
    return H0 + (H1 - H0) * np.asarray(x, float) / LPATH


def wv(x):
    """Channel width that holds S(x) exact under a 45-degree vaulted ceiling."""
    h = Hx(x)
    return 2 * h - 2 * np.sqrt(np.maximum(h * h - S(x), 0.0))


def rc(x):
    return R1 + KSP * np.asarray(x, float) / math.sqrt(1 + KSP ** 2)


def th_of(x):
    return np.log(rc(x) / R1) / KSP


X_FOLD = FFOLD * LPATH
TH_FOLD = float(th_of(X_FOLD))


# ------------------------------------------------------------------ sections

def outline(w, h, wall=0.0):
    """Duct section in local (u = lateral, v = up); wall > 0 offsets outward.

    Fixed point budget per edge so consecutive sections loft without twisting;
    Q scales the budget so the skin fields resolve (~1-3 mm along the roofs)."""
    w2 = 0.5 * w + wall
    hs = (h - 0.5 * w) + (wall * math.sqrt(2.0) if wall else 0.0)
    v0 = -wall
    apex = hs + w2
    e = lambda a, b, n: np.linspace(a, b, n, endpoint=False)
    u = np.concatenate([e(-w2, w2, 8 * Q), np.full(4 * Q, w2), e(w2, 0.0, 8 * Q),
                        e(0.0, -w2, 8 * Q), np.full(4 * Q, -w2)])
    v = np.concatenate([np.full(8 * Q, v0), e(v0, hs, 4 * Q), e(hs, apex, 8 * Q),
                        e(apex, hs, 8 * Q), e(hs, v0, 4 * Q)])
    return np.stack([u, v], -1)


def outline_circle(r):
    """SEC points round a circle, matched to outline()'s perimeter fractions."""
    ref = outline(2 * r, 2 * r)
    seg = np.linalg.norm(np.diff(np.vstack([ref, ref[:1]]), axis=0), axis=1)
    frac = np.concatenate([[0.0], np.cumsum(seg)[:-1]]) / seg.sum()
    a = -0.5 * math.pi + 2 * math.pi * frac
    return np.stack([r * np.cos(a), r * np.sin(a)], -1)


# ---------------------------------------------------------------- the skin
#
# The Volute has no axis of revolution, so its texture is applied to the loft
# SECTIONS: each outline point moves along its own 2-D outward normal by a
# field of (path position x, perimeter position p) -- ridge crests that run
# ALONG the duct like growth lines, a serrated keel on the vault apex, and the
# grip carved into the first turn's roof. Only exposed faces move: the roofs
# everywhere, the walls only on the bell (between the body's turns there is
# 1.4 mm of air and nothing more -- the packing gate stays untouched).

def _outline_normals(o):
    """Per-point outward 2-D normals of a closed outline (periodic)."""
    d = np.roll(o, -1, 0) - np.roll(o, 1, 0)
    n = np.stack([d[:, 1], -d[:, 0]], -1)
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)


def _sstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _keel(x, x0, x1, h0, h1, pitch=14.0):
    """Serrated sawtooth along the apex, teeth growing toward the mouth."""
    x = np.asarray(x, float)
    tooth = 1.0 - np.abs(2.0 * ((x / pitch) % 1.0) - 1.0)
    amp = (h0 + (h1 - h0) * np.clip((x - x0) / max(x1 - x0, 1e-9), 0, 1))
    return amp * tooth * H.smooth_band(x, x0, x1, 10.0)


def texture_secs(x, secs, ridges=None, keel=None, grip=False, walls=False,
                 seed=31):
    """Displace loft sections: ridges=(x0, x1, amp_lo, amp_hi); keel=(x0, x1,
    h_lo, h_hi); grip carves the cradle handprint; walls extends the ridge
    field down the side walls. Any of them may be None."""
    x = np.asarray(x, float)
    n, sec_n = len(x), len(secs[0])
    p = (np.arange(sec_n) + 0.5) / sec_n
    if ridges is not None:
        rx0, rx1, a_lo, a_hi = ridges
        F = H.ridge_noise(x, 2 * math.pi * p, r_ref=80.0, z0=rx0, z1=rx1,
                          amp_lo=a_lo, amp_hi=a_hi, cell_z=30.0, cell_arc=11.0,
                          seed=seed)
    else:
        F = np.zeros((n, sec_n))
    kf = _keel(x, *keel) if keel is not None else np.zeros(n)
    out = []
    for i, xi in enumerate(x):
        o = np.array(secs[i], float)
        xc = float(np.clip(xi, 0.0, LPATH))
        w2 = 0.5 * float(wv(xc)) + WALL
        hs = (float(Hx(xc)) - 0.5 * float(wv(xc))) + WALL * math.sqrt(2.0)
        apex = hs + w2
        u, v = o[:, 0], o[:, 1]
        g = np.clip((v - hs) / max(apex - hs, 1e-9), 0.0, 1.0)
        wgt = _sstep(g / 0.18)                        # roofs, fading in off the
        if walls:                                     # shoulder; bell walls too,
            wgt = np.maximum(wgt, _sstep((v - 1.0) / 5.0)   # kept off the plate
                             * (v < hs - 1.0))              # notch at the floor
        disp = F[i] * wgt
        if grip:
            gw = float(H.smooth_band(np.array([xi]), GRIP_X[0], GRIP_X[1], 10.0)[0])
            right = (u > 1e-6) & (v > hs)
            face = (apex - hs) * math.sqrt(2.0)       # roof face length, mm
            gs = _sstep(np.minimum(g / 0.22, (0.95 - g) / 0.15))
            pit = np.zeros(sec_n)
            for xc_, gc, sz_, sa_ in FINGERS + (THUMB,):
                q = ((xi - xc_) / sz_) ** 2 + ((g - gc) * face / sa_) ** 2
                s = np.clip(1.0 - q ** 1.5, 0.0, 1.0)
                pit += GRIP_DEPTH * s * s * (3.0 - 2.0 * s)
            bulge = GRIP_RISE * gw * gs
            disp = np.where(right, disp * (1.0 - gw) + bulge - np.minimum(pit, bulge),
                            disp)
        o = o + _outline_normals(o) * disp[:, None]
        if kf[i] > 0:                                 # the keel rides the apex
            j = int(np.argmax(v))
            for dj, f in ((0, 1.0), (1, 0.55), (-1, 0.55), (2, 0.2), (-2, 0.2)):
                o[(j + dj) % sec_n, 1] += f * kf[i]
        out.append(o)
    return out


def spiral_frames(x):
    """Centre, lateral (radial) and up vectors for each station on the spiral."""
    t = th_of(x)
    r = rc(x)
    ct, st = np.cos(t), np.sin(t)
    ctr = np.stack([r * ct, r * st, np.zeros_like(r)], -1)
    lat = np.stack([ct, st, np.zeros_like(r)], -1)
    return ctr, lat


def straight_frames(x, x0):
    """Tangential run past the fold: the bell leaves on the spiral's tangent.

    Its sections keep the RADIAL lateral axis of the spiral sweep, not one
    perpendicular to the tangent -- otherwise every section would be rotated by
    the spiral's 5.6 degree pitch angle and the bell would not mate with the
    body's port. The cost is that sections are sheared 5.6 deg off normal, which
    scales the acoustic area by cos 5.6 = 0.995."""
    t0 = float(th_of(x0))
    r0 = float(rc(x0))
    tan = np.array([KSP * math.cos(t0) - math.sin(t0),
                    KSP * math.sin(t0) + math.cos(t0)])
    tan /= np.linalg.norm(tan)
    lat2 = np.array([math.cos(t0), math.sin(t0)])          # radial, as swept
    p0 = np.array([r0 * math.cos(t0), r0 * math.sin(t0)])
    s = (np.asarray(x, float) - x0)[:, None]
    ctr = np.concatenate([p0 + s * tan, np.zeros((len(s), 1))], 1)
    lat = np.repeat(np.concatenate([lat2, [0.0]])[None, :], len(s), 0)
    return ctr, lat


def body_frames(x):
    """Frames of the body's own duct: spiral inside the fold, tangent outside.

    Both plugs are built on THIS, not on their own part's axis -- over 18 mm of
    insertion the spiral's lateral axis turns 9 degrees, and a plug swept on a
    straight axis drives its corners straight into the duct wall."""
    x = np.atleast_1d(np.asarray(x, float))
    ctr, lat = spiral_frames(np.clip(x, 0.0, X_FOLD))
    for i, xi in enumerate(x):
        if xi < 0.0 or xi > X_FOLD:
            c, l = straight_frames(np.array([xi]), 0.0 if xi < 0 else X_FOLD)
            ctr[i], lat[i] = c[0], l[0]
    return ctr, lat


def sweep(x, frames, wall=0.0, sections=None):
    """Place duct sections along a path -> (ns, SEC, 3) for loft()."""
    ctr, lat = frames
    out = []
    for i, xi in enumerate(np.atleast_1d(x)):
        o = outline(float(wv(xi)), float(Hx(xi)), wall) if sections is None \
            else sections[i]
        p = (ctr[i][None, :] + o[:, :1] * lat[i][None, :]
             + o[:, 1:2] * np.array([0.0, 0.0, 1.0])[None, :])
        out.append(p)
    return np.asarray(out)


def tube(sections):
    m = H.loft(sections)
    if m.volume < 0:
        m.invert()
    return m


# ------------------------------------------------------------------- checks

def pack_report():
    """The corrected gate, published layout vs built layout."""
    def worst(r1, h0, h1, xf, wall=WALL, n=400):
        S_ = lambda x: math.pi * RT * RT * math.exp(MFL * x)
        Hf = lambda x: h0 + (h1 - h0) * x / LPATH
        wf = lambda x: 2 * Hf(x) - 2 * math.sqrt(max(Hf(x) ** 2 - S_(x), 0.0))
        tf = math.log((r1 + KSP * xf / math.sqrt(1 + KSP ** 2)) / r1) / KSP
        bad = -1e9
        for t in np.linspace(0, max(tf - 2 * math.pi, 0.0), n):
            ra, rb = r1 * math.exp(KSP * t), r1 * math.exp(KSP * (t + 2 * math.pi))
            xa = (ra - r1) * math.sqrt(1 + KSP ** 2) / KSP
            xb = (rb - r1) * math.sqrt(1 + KSP ** 2) / KSP
            if xb > xf:
                break
            bad = max(bad, 0.5 * wf(xa) + 0.5 * wf(xb) + 2 * wall - (rb - ra))
        return bad, tf / 2 / math.pi
    for name, r1, h0, h1, ff in (("published c10", 45.5170, 95.4316, 245.5498, 0.6745),
                                 ("built here   ", R1, H0, H1, FFOLD)):
        b, turns = worst(r1, h0, h1, ff * LPATH)
        print(f"  {name}: r1 {r1:5.1f}  H {h0:5.1f}->{h1:5.1f}  fold {ff:.4f}  "
              f"turns {turns:4.2f}  worst wall-to-wall {-b:+6.1f} mm "
              f"{'OK' if b <= 0 else 'OVERLAP'}")
    print(f"  f_c {MFL * 1000 * 343 / (4 * math.pi):.0f} Hz   throat A {S(0):.0f}   "
          f"mouth A {S(LPATH):.0f} mm2   path {LPATH:.0f} mm "
          f"(fold {X_FOLD:.0f} + bell {LPATH - X_FOLD:.0f})")


# -------------------------------------------------------------------- parts

def spiral_body(skin=True):
    """Folded duct + base plate, built as sub-turn lofts so nothing self-
    intersects, then unioned; the duct void is subtracted last."""
    per = 1 + int(math.ceil(X_FOLD / 2.0))
    x = np.linspace(0.0, X_FOLD, per)           # ONE sweep: the corrected gate
    secs = [outline(float(wv(xi)), float(Hx(xi)), WALL) for xi in x]
    if skin:
        secs = texture_secs(x, secs, grip=True)   # just the cradle hand print
    outer = tube(sweep(x, spiral_frames(x), sections=secs))  # keeps the turns
    # apart, so the outer surface never self-intersects, needs no per-turn union
    r_max = float(rc(X_FOLD) + 0.5 * wv(X_FOLD)) + WALL + 1.0
    plate = H.cyl((0, 0, -PLATE), (0, 0, -0.5), r_max, 256)
    hub = H.cyl((0, 0, -PLATE - 1), (0, 0, 1), R_HUB, 64)
    body = H.union(outer, plate)
    body = H.difference(body, hub)

    xv = np.concatenate([[-10.0, -5.0], np.linspace(0, X_FOLD, 2 * per),
                         [X_FOLD + 5.0, X_FOLD + 10.0]])   # just enough to punch
    # through the end caps: a longer backward run would graze the outer turn
    sec = [outline(float(wv(np.clip(xi, 0.0, X_FOLD))),
                   float(Hx(np.clip(xi, 0.0, X_FOLD)))) for xi in xv]
    void = tube(sweep(xv, body_frames(xv), sections=sec))
    body = H.difference(body, void, *neighbour_notches())
    if skin:
        body = H.union(body, *volute_skin())
    return body


def neighbour_notches(clear=0.3):
    """Cut the base plate back where the bell and the mouthpiece land.

    The plate is a full disc, so without this it runs under both neighbours and
    their floors would have to be shaved to 0.5 mm to clear it. Cutting the plate
    instead leaves every floor at full thickness and turns the joint into a notch
    that locates the part."""
    out = []
    xb = np.linspace(X_FOLD, X_FOLD + 150.0, 20)   # out past the plate rim
    cb, lb = straight_frames(xb, X_FOLD)
    out.append(tube(sweep(xb, (cb, lb),
                          sections=[outline(float(wv(xi)), float(Hx(xi)),
                                            WALL + clear) for xi in xb])))
    ctr, lat, sec_o = mouthpiece_path(grow=clear)
    out.append(tube(sweep(np.arange(len(ctr)), (ctr, lat), sections=sec_o)))
    return out


def volute_skin():
    """Phyllotaxis along the spiral path: small studs on the inner whorl, sharp
    claws curling downstream through the outer turn (RULE S2/S3). The cradle
    grip window stays bare -- that is where the hand goes."""
    n = 340
    x = np.linspace(8.0, X_FOLD - 8.0, n)
    ctr, lat = spiral_frames(x)
    w2o = 0.5 * wv(x) + WALL                       # outer half-width
    hso = (Hx(x) - 0.5 * wv(x)) + WALL * math.sqrt(2.0)     # outer shoulder
    u = w2o * 0.85 * (2 * ((np.arange(n) * 0.6180339887) % 1.0) - 1.0)
    keep = ~((x > GRIP_X[0] - 8.0) & (x < GRIP_X[1] + 8.0) & (u > 0))
    x, ctr, lat, w2o, hso, u = (a[keep] for a in (x, ctr, lat, w2o, hso, u))
    frac = (x - 8.0) / (X_FOLD - 16.0)
    P = ctr + lat * u[:, None]
    P[:, 2] = hso + (w2o - np.abs(u))              # the 45-degree roof itself
    s = np.sign(u)[:, None] / math.sqrt(2.0)
    axis = lat * s + np.array([[0.0, 0.0, 1 / math.sqrt(2.0)]])
    t = th_of(x)
    tan = np.stack([KSP * np.cos(t) - np.sin(t),
                    KSP * np.sin(t) + np.cos(t), np.zeros_like(t)], -1)
    tan /= np.linalg.norm(tan, axis=1, keepdims=True)
    h = 0.8 + 7.2 * (np.exp(2.2 * frac) - 1) / (math.exp(2.2) - 1)
    rb = np.clip(h / 7.0, 0.7, 2.2)
    return H.claws(P, axis, h, rb, r_tip=0.18, curl=tan, curl_amt=0.5, sink=2.5)


def bell_skin():
    """Claws along the bell's vault: armored where it speaks, growing toward
    the mouth, curling with the flow."""
    n = 150
    x = np.linspace(X_FOLD + 55.0, LPATH - 15.0, n)
    ctr, lat = straight_frames(x, X_FOLD)
    w2o = 0.5 * wv(x) + WALL
    hso = (Hx(x) - 0.5 * wv(x)) + WALL * math.sqrt(2.0)
    u = w2o * 0.85 * (2 * ((np.arange(n) * 0.6180339887) % 1.0) - 1.0)
    P = ctr + lat * u[:, None]
    P[:, 2] = hso + (w2o - np.abs(u))
    s = np.sign(u)[:, None] / math.sqrt(2.0)
    axis = lat * s + np.array([[0.0, 0.0, 1 / math.sqrt(2.0)]])
    t0 = float(TH_FOLD)
    tan = np.array([KSP * math.cos(t0) - math.sin(t0),
                    KSP * math.sin(t0) + math.cos(t0), 0.0])
    tan /= np.linalg.norm(tan)
    frac = (x - x[0]) / (x[-1] - x[0])
    h = 2.0 + 6.0 * frac
    rb = np.clip(h / 7.0, 0.7, 2.2)
    return H.claws(P, axis, h, rb, r_tip=0.18, curl=tan, curl_amt=0.5, sink=2.5)


def hollow(ctr, lat, sec_out, sec_in, ext=(1.5, 1.5)):
    """Shell between two section lists, with the void run past both ends so the
    duct opens instead of sealing itself into a closed cavity."""
    idx = np.arange(len(ctr))
    outer = tube(sweep(idx, (ctr, lat), sections=sec_out))
    cv = np.vstack([ctr[:1] + (ctr[0] - ctr[1]) * ext[0], ctr,
                    ctr[-1:] + (ctr[-1] - ctr[-2]) * ext[1]])
    lv = np.vstack([lat[:1], lat, lat[-1:]])
    void = tube(sweep(np.arange(len(cv)), (cv, lv),
                      sections=[sec_in[0]] + list(sec_in) + [sec_in[-1]]))
    return H.difference(outer, void)


def bell(skin=True):
    """The unrolled 283.6 mm of path: straight, vaulted, floor-down, one piece.

    Its inlet is a plug landing inside the body's port -- at that joint the two
    turns' walls touch, so there is no room outside for a sleeve or a flange."""
    x = np.linspace(X_FOLD, LPATH, 90)
    ctr, lat = straight_frames(x, X_FOLD)
    sec_o = [outline(float(wv(xi)), float(Hx(xi)), WALL) for xi in x]
    sec_i = [outline(float(wv(xi)), float(Hx(xi))) for xi in x]
    m = hollow(ctr, lat, sec_o, sec_i, ext=(6.0, 1.5))

    xp = np.linspace(X_FOLD - PLUG_LEN, X_FOLD + 3.0, 9)
    cp, lp = body_frames(xp)
    o, i_ = plug_sections(X_FOLD, xp)
    m = H.union(m, hollow(cp, lp, o, i_, ext=(2.0, 0.5)))
    if skin:
        m = H.union(m, *bell_skin())
    return m


def plug_sections(x0, xp, own=+1):
    """Lap-joint spigot: inside its own part it bites into the local bore wall so
    the two solids overlap; where it protrudes it follows the mating duct's own
    taper, stepped in by the fit clearance."""
    o, i_ = [], []
    for xi in xp:
        w, h = float(wv(xi)), float(Hx(xi))       # follow the taper: a constant
        s_ = (xi - x0) * own                      # section would jam going in,
        if s_ >= 0:                               # since the duct narrows upstream
            d = BITE                              # own part: bite into its wall
        else:
            d = -PLUG_CLR * min(1.0, -s_ / 2.0)   # protruding: fit clearance
        o.append(outline(w, h, d))
        i_.append(outline(w, h, d - PLUG_WALL))
    return o, i_


def mouthpiece():
    """Cup adapter: the modelled 40 mm lip taper, curved in the horizontal plane
    so the dia-45 cup lands in the clear centre of the shell instead of running
    into the second turn. Plugs into the body's throat port; pull it off to
    print a different cup."""
    ctr, lat, sec_o, sec_i, w0, h0, step = mouthpiece_path(full=True)
    m = hollow(ctr, lat, sec_o, sec_i, ext=(1.5, 1.5))

    xp = np.linspace(-3.0, PLUG_LEN, 9)                       # 3 mm inside -> tip
    cp, lp = body_frames(xp)
    o, i_ = plug_sections(0.0, xp[::-1], own=-1)
    plug = hollow(cp[::-1], lp[::-1], o, i_, ext=(0.5, 2.0))
    return H.union(m, plug)


def mouthpiece_path(full=False, grow=0.0):
    n = 26
    t = np.linspace(0.0, 1.0, n)
    back = np.array([0.0, -1.0])            # port-plane normal at theta = 0, out
    turn = math.radians(TURN_DEG) * t
    dirs = -np.stack([back[0] * np.cos(turn) - back[1] * np.sin(turn),
                      back[0] * np.sin(turn) + back[1] * np.cos(turn)], -1)
    step = CUP_LEN / (n - 1)
    pts = [np.array([R1, 0.0])]
    for i in range(1, n):
        pts.append(pts[-1] - step * dirs[i])                  # back out of the port
    ctr = np.concatenate([np.asarray(pts), np.zeros((n, 1))], 1)
    lat = np.stack([dirs[:, 1], -dirs[:, 0], np.zeros(n)], -1)   # radial at root

    w0, h0 = float(wv(0.0)), float(Hx(0.0))
    circ = outline_circle(CUP_R) + np.array([0.0, 0.5 * h0])
    circ_o = outline_circle(CUP_R + WALL + grow) + np.array([0.0, 0.5 * h0])
    rect, rect_o = outline(w0, h0), outline(w0, h0, WALL + grow)
    blend = 0.5 - 0.5 * np.cos(math.pi * t)
    sec_i = [rect * (1 - s) + circ * s for s in blend]
    sec_o = [rect_o * (1 - s) + circ_o * s for s in blend]
    if full:
        return ctr, lat, sec_o, sec_i, w0, h0, step
    return ctr, lat, sec_o


def build(skin=True):
    print("VOLUTE (cycle-10 KNEE acoustics, re-laid-out fold)")
    pack_report()
    rows = []
    b = spiral_body(skin=skin)
    rows.append(H.report("volute_body", b, OUT + "volute_body.stl"))
    rows.append(H.report("volute_bell", bell(skin=skin), OUT + "volute_bell.stl"))
    rows.append(H.report("volute_mouthpiece", mouthpiece(),
                         OUT + "volute_mouthpiece.stl"))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skin", action="store_true",
                    help="add the ornament pass; default is plain smooth bodies")
    ap.add_argument("--pack", action="store_true", help="packing check only")
    a = ap.parse_args()
    if a.pack:
        pack_report()
    else:
        rows = build(skin=a.skin)
        print()
        H.print_table(rows)
