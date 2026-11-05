#!/usr/bin/env python3
"""Evolutionary simulation of the Plumb Octahedron (Concept A): tuning and theta, as built.

The rev-1 intent is checked against a physical model of the built object:

* Static equilibrium of the hanging frame, solved per state by Newton on 8 unknowns: the cross (hub, tubes,
  tip nodes) hangs from the hook on 4 heaven strings and may translate and tilt; the bottom node hangs on
  4 earth strings and carries the bob. Waist strings are internal to the rigid cross.
* Steel strings with real axial stiffness (EA/L) in series with optional springs; tension is stretch x
  stiffness, so a 4-string star is statically indeterminate exactly as it will be in the workshop.
* Install tuning simulated as a builder would do it: shot in the bob, trim in the tips, anchor adjusters,
  iterated until every string reads its target at the install temperature.
* Monte Carlo manufacturing tolerances with common random numbers (gauge, density, modulus, saddle and
  tube positions, frame masses, spring rates, PETG expansion, Strouhal law, creep).
* Environment: three seasons x a daily cycle, sun on one side, PETG creep after months, gust drag.
* Wind: Weibull speeds, Williamson Strouhal shedding, lock-in on the nearest stiff-string harmonic
  (inharmonicity included).

scipy's differential evolution searches design + tuning plan for the closest match to the intent: fixed
fundamentals, just intervals at the harmonics the wind sounds, theta beats (4-8 Hz) for as much of the
sounding time as possible, a hub drone on the root; inside mass, string-stress, spring and bed limits.

    python evolve_tuning.py --check     # rev 1, no tolerances: compare with plumb_octahedron.py
    python evolve_tuning.py             # score rev 1, evolve, write evolution_*.json / .csv
"""
import argparse
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import gamma

G = 9.81
NU_AIR = 1.5e-5
ROOT = 146.832                                   # D3
JUST = np.array([1.0, 6 / 5, 3 / 2, 9 / 5])      # minor seventh chord
HEAVEN = 9 / 8                                   # heaven chord a whole tone above
CATALOG = np.array([.008, .009, .010, .011, .012, .013, .014, .015, .016, .017, .018, .019, .020,
                    .022, .024, .026])           # plain steel gauges, inch
HUB_R, TIP_R = 125.0, 55.0
S_TIP, S_END = 30.0, 25.0                        # saddle offsets from node centres, mm
BREAK_MPA = 2600.0
M_HUB, M_TIPNODE, M_BOTTOM = 0.45, 0.12, 0.08
RHAT = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
EZ = np.array([0.0, 0.0, 1.0])
LOCK_SIGMA = 0.06                                # lock-in width, natural-log units
N_MC = 12
SEASONS = (6.0, 14.0, 22.0)
HOURS = (0, 4, 8, 12, 16, 20)
U_GRID = np.arange(0.25, 12.01, 0.25)
STEPS = np.array([1e-3, 1e-3, 1e-3, 1e-7, 1e-7, 1e-3, 1e-3, 1e-3])
LIMITS = np.array([20, 20, 20, 0.02, 0.02, 20, 20, 20])
BOUNDS = [(500, 690), (-0.3, 2.0), (-0.3, 1.3), (0.2, 2.5), (0.0, 30.0), (1.0, 14.0)] + [(0, 15)] * 8 \
    + [(0, 1), (0, 1)]                            # frame: 0 PETG, 1 carbon cores; tensioner: 0 springs, 1 steelyard
INTEGRALITY = [False] * 6 + [True] * 10
REV1 = np.array([609.5 + S_TIP + S_END, 6.0, math.log10(2.0), 0.8, 20.0, 1.0, 10, 7, 4, 2, 10, 7, 4, 2, 0, 0], float)
LEVER = 8.0                                       # steelyard advantage: poise weight = tension / 8
POISE_ARM = 120.0                                 # mm, poise working radius


# ---------------------------------------------------------------- tolerances and environment
def tolerances(n, seed=7, zero=False):
    rng = np.random.default_rng(seed)
    z = 0.0 if zero else 1.0
    return dict(
        dg=z * rng.normal(0, 0.00008, (n, 12)),              # gauge error, inch
        rho=7850 * (1 + z * rng.normal(0, 0.004, (n, 12))),
        E=2.0e5 * (1 + z * rng.normal(0, 0.03, (n, 12))),    # N/mm^2
        sad=z * rng.normal(0, 0.7, (n, 12)),                 # saddle-sum error, mm
        dr=z * rng.normal(0, 0.5, (n, 4)),                   # tip radius error, mm
        m_hub=M_HUB * (1 + z * rng.normal(0, 0.1, n)),
        m_tipnode=M_TIPNODE * (1 + z * rng.normal(0, 0.1, (n, 4))),
        m_tube=1 + z * rng.normal(0, 0.05, (n, 4)),
        m_bottom=M_BOTTOM * (1 + z * rng.normal(0, 0.1, n)),
        kfac=1 + z * rng.normal(0, 0.08, (n, 12)),
        cte=65e-6 if zero else rng.uniform(55e-6, 75e-6, n),
        cte_cf=9e-6 if zero else rng.uniform(7e-6, 11e-6, n),     # carbon core + PETG tip node
        st=1 + z * rng.normal(0, 0.04, n),
        creep=rng.uniform(0, 0.0012, (n, 4)) * (0 if zero else 1),
        sun=rng.uniform(2, 8, n) * (0 if zero else 1),
    )


def states():
    """(temperature, daylight 0..1, aged 0/1) for three seasons x daily cycle x fresh/aged."""
    out = []
    for mean in SEASONS:
        for h in HOURS:
            t = mean + 7 * math.cos(2 * math.pi * (h - 15) / 24)
            day = max(0.0, math.cos(2 * math.pi * (h - 12) / 24))
            for aged in (0, 1):
                out.append((t, day, aged))
    return np.array(out)


def weibull_weights(mean):
    c = mean / gamma(1 + 1 / 2)
    pdf = (2 / c) * (U_GRID / c) * np.exp(-(U_GRID / c) ** 2)
    return pdf / pdf.sum()


def strouhal(re):
    return np.where(re < 178, -3.3265 / np.maximum(re, 1) + 0.1816 + 1.6e-4 * re, 0.19)


# ---------------------------------------------------------------- statics
def tips_world(x, tipR):
    c = x[..., 0:3]
    w = np.stack([x[..., 3], x[..., 4], np.zeros_like(x[..., 3])], -1)[..., None, :]
    arm = tipR[..., :, None] * RHAT
    arm = arm + np.cross(w, arm)
    return c[..., None, :] + arm, arm


def residual(x, P):
    tips, arm = tips_world(x, P["tipR"])
    p = x[..., 5:8]
    vu = -tips
    du = np.linalg.norm(vu, axis=-1)
    tu = P["k"][..., 0:4] * np.maximum(du - P["ell"][..., 0:4], 0)
    vl = p[..., None, :] - tips
    dl = np.linalg.norm(vl, axis=-1)
    tl = P["k"][..., 4:8] * np.maximum(dl - P["ell"][..., 4:8], 0)
    fu = tu[..., None] * vu / du[..., None]
    fl = tl[..., None] * vl / dl[..., None]
    ftip = fu + fl - (P["m_tip"][..., :, None] * G) * EZ
    wtube = -(P["m_tube"][..., :, None] * G) * EZ
    fc = ftip.sum(-2) + wtube.sum(-2) - (P["m_hub"][..., None] * G) * EZ + P["drag_c"]
    mc = np.cross(arm, ftip).sum(-2) + np.cross(arm * P["mid"][..., None, None], wtube).sum(-2)
    fn = -fl.sum(-2) - (P["m_bob"][..., None] * G) * EZ + P["drag_b"]
    r = np.concatenate([fc, mc[..., 0:2] / 400.0, fn], -1)
    return r, tu, tl, du, dl, tips, p


def solve(P, x, iters=14):
    Pb = {k: v[:, None] for k, v in P.items()}
    ok = np.zeros(len(x), bool)
    for _ in range(iters):
        xx = np.repeat(x[:, None, :], 9, axis=1)
        xx[:, 1:, :] += np.diag(STEPS)[None]
        r = residual(xx, Pb)[0]
        jac = np.swapaxes((r[:, 1:, :] - r[:, :1, :]) / STEPS[None, :, None], 1, 2)
        jac = jac + np.eye(8)[None] * 1e-9
        try:
            dx = np.linalg.solve(jac, -r[:, 0, :][..., None])[..., 0]
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(jac[0], -r[0, 0], rcond=None)[0][None].repeat(len(x), 0)
        scale = np.max(np.abs(dx) / LIMITS, axis=1, keepdims=True)
        x = x + dx / np.maximum(1, scale)
        ok = np.max(np.abs(r[:, 0, :]), axis=1) < 1e-6
        if ok.all():
            break
    out = residual(x, P)
    ok = np.max(np.abs(out[0]), axis=1) < 1e-4
    return x, out, ok


# ---------------------------------------------------------------- one design
def decode(v):
    edge, lkg, lkw, delta, t_set, n_c = v[:6]
    ge = CATALOG[np.clip(np.round(v[6:10]).astype(int), 0, len(CATALOG) - 1)]
    gh = CATALOG[np.clip(np.round(v[10:14]).astype(int), 0, len(CATALOG) - 1)]
    frame = int(round(v[14])) if len(v) > 14 else 0
    tens = int(round(v[15])) if len(v) > 15 else 0
    kg = 0.02 if tens else 10 ** lkg               # a steelyard is a near-zero-rate tensioner
    return dict(edge=edge, kg=kg, kw=10 ** lkw, delta=delta, t_set=t_set, n_c=n_c, ge=ge, gh=gh,
                frame=frame, tens=tens)


def build(D, tol, n, shift=0.0):
    """Per-sample constants for one design: strings 0-3 heaven, 4-7 earth, 8-11 waist.
    shift transposes every target by that many cents (tunability checks)."""
    g_in = np.concatenate([D["gh"], D["ge"], D["ge"]])[None, :] + tol["dg"][:n]
    d_m = g_in * 0.0254
    area = math.pi * (d_m * 1000) ** 2 / 4                       # mm^2
    mu = tol["rho"][:n] * math.pi * d_m ** 2 / 4                 # kg/m
    E = tol["E"][:n]
    sad = np.concatenate([np.full(8, S_TIP + S_END), np.full(4, 2 * S_TIP)])[None, :] + tol["sad"][:n]
    ks = np.concatenate([np.full(8, D["kg"]), np.full(4, D["kw"])])[None, :] * tol["kfac"][:n]
    ell0 = np.full((n, 12), D["edge"])
    k = 1 / (1 / ks + ell0 / (E * area))
    r0 = D["edge"] / math.sqrt(2) + tol["dr"][:n]
    tube_len = D["edge"] / math.sqrt(2) - HUB_R - TIP_R
    m_tube = 1.27e-3 * math.pi / 4 * (30 ** 2 - 26 ** 2) * tube_len / 1000 * tol["m_tube"][:n] \
        + (0.03 if D.get("frame") else 0.0)
    mid = (HUB_R + tube_len / 2) / (D["edge"] / math.sqrt(2))

    # targets with inharmonicity pre-compensation at harmonic n_c
    L = (D["edge"] - S_TIP - S_END) / 1000
    def bcoef(gin, f):
        dm = gin * 0.0254
        T = (2 * L * f) ** 2 * 7850 * math.pi * dm ** 2 / 4
        return math.pi ** 3 * 2e11 * dm ** 4 / (64 * T * L ** 2)
    fe = ROOT * JUST
    fh = ROOT * HEAVEN * JUST
    be = np.array([bcoef(g, f) for g, f in zip(D["ge"], fe)])
    bh = np.array([bcoef(g, f) for g, f in zip(D["gh"], fh)])
    nc = D["n_c"]
    fe = fe * math.sqrt(1 + be[0] * nc ** 2) / np.sqrt(1 + be * nc ** 2)
    fh = fh * math.sqrt(1 + bh[0] * nc ** 2) / np.sqrt(1 + bh * nc ** 2)
    fstar = np.concatenate([fh, fe, fe + D["delta"]]) * 2 ** (shift / 1200)
    return dict(d_m=d_m, area=area, mu=mu, E=E, sad=sad, ks=ks, k=k, r0=r0, tube_len=tube_len,
                m_tube=m_tube, mid=mid, fstar=fstar)


def tension_targets(S, dist):
    lv = (dist - S["sad"]) / 1000
    return (2 * lv * S["fstar"]) ** 2 * S["mu"]


def install(D, tol, n, shift=0.0):
    S = build(D, tol, n, shift)
    h = D["edge"] / math.sqrt(2)
    x = np.zeros((n, 8))
    x[:, 2], x[:, 7] = -h, -2 * h
    tipR = S["r0"]
    dist = np.full((n, 12), D["edge"])
    span = np.sqrt(tipR ** 2 + np.roll(tipR, -1, axis=1) ** 2)
    dist[:, 8:] = span
    tstar = tension_targets(S, dist)
    ell = dist - tstar / S["k"]
    # steelyards sit at the poles: heaven poises hang from the hook (outside the frame), earth poises
    # are part of the bob, so neither loads the tip nodes
    tipnode = tol["m_tipnode"][:n].copy()
    m_bob = (tstar[:, 4:8] * math.cos(math.pi / 4)).sum(1) / G - tol["m_bottom"][:n]
    frame = m_bob + tol["m_bottom"][:n] + tol["m_hub"][:n] + S["m_tube"].sum(1) + tipnode.sum(1)
    trim = ((tstar[:, 0:4] * math.cos(math.pi / 4)).sum(1) / G - frame) / 4
    ell[:, 8:] = span - tstar[:, 8:] / S["k"][:, 8:]

    def forward(y, x0, reps):
        """y: (n*reps, 10) = 8 gravity anchor lengths, bob mass, trim per tip."""
        rr = lambda a: np.repeat(a, reps, axis=0)
        nn = len(y)
        P = dict(tipR=rr(tipR), ell=np.concatenate([y[:, :8], rr(ell[:, 8:])], 1), k=rr(S["k"]),
                 m_tip=rr(tipnode) + y[:, 9:10], m_tube=rr(S["m_tube"]),
                 m_hub=rr(tol["m_hub"][:n]), m_bob=y[:, 8] + rr(tol["m_bottom"][:n]),
                 mid=np.full(nn, S["mid"]), drag_c=np.zeros((nn, 3)), drag_b=np.zeros((nn, 3)))
        return solve(P, x0)

    # the builder's tuning passes as a coupled Newton step: anchors, shot and trim together
    steps = np.concatenate([np.clip(0.3 / S["k"][:, :8], 0.005, 50), np.full((n, 1), 0.01),
                            np.full((n, 1), 0.005)], 1)                       # (n, 10), ~0.3 N per step
    y = np.concatenate([ell[:, :8], m_bob[:, None], trim[:, None]], 1)
    step_lim = np.concatenate([6.0 / S["k"][:, :8], np.full((n, 1), 2.0), np.full((n, 1), 0.5)], 1)
    for _ in range(14):
        yy = np.repeat(y[:, None, :], 11, axis=1)
        yy[:, 1:, :] += steps[:, :, None] * np.eye(10)[None]
        xs, (r, tu, tl, du, dl, tips, p), ok = forward(yy.reshape(n * 11, 10), np.repeat(x, 11, 0), 11)
        T8 = np.concatenate([tu, tl], 1).reshape(n, 11, 8)
        d8 = np.concatenate([du, dl], 1).reshape(n, 11, 8)
        lv = (d8 - S["sad"][:, None, :8]) / 1000
        F = T8 - (2 * lv * S["fstar"][:8]) ** 2 * S["mu"][:, None, :8]
        jac = np.swapaxes((F[:, 1:, :] - F[:, :1, :]) / steps[:, :, None], 1, 2)       # (n, 8, 10)
        dy = -(np.linalg.pinv(jac) @ F[:, 0, :, None])[..., 0]
        scale = np.max(np.abs(dy) / step_lim, axis=1, keepdims=True)
        y = y + dy / np.maximum(1, scale)
        x = xs.reshape(n, 11, 8)[:, 0]
        if np.abs(F[:, 0, :]).max() < 1e-4:
            break
    x, (r, tu, tl, du, dl, tips, p), ok = forward(y, x, 1)
    ell[:, :8] = y[:, :8]
    m_bob, trim = y[:, 8], y[:, 9]
    m_tip = tipnode + trim[:, None]
    tw = S["k"][:, 8:] * np.maximum(span - ell[:, 8:], 0)
    T = np.concatenate([tu, tl, tw], 1)
    return S, dict(x=x, ell=ell, m_bob=m_bob, trim=trim, m_tip=m_tip, T=T, ok=ok, jac=jac,
                   dist=np.concatenate([du, dl, span], 1))


def tunability(D, tol, n, S, I):
    """Can a person retune it?
    * whole harp +/-200 cents (shot and trim, or every poise);
    * any one gravity string +/-100 cents on its own adjuster or poise, the others held;
    * resolution: cents per turn of an M5 x 0.8 adjuster, or per 4 mm of poise.
    Limits: springs - anchor travel 12 mm, spring stretch 40 mm, 40 cents per turn;
            steelyard - poise travel 20 mm for one string and 40 mm for the whole harp."""
    steel = bool(D["tens"])
    lim_one, lim_all = (20.0, 40.0) if steel else (12.0, 12.0)
    v, out = 0.0, {}
    for shift in (-200.0, 200.0):
        St, It = install(D, tol, n, shift)
        mass = It["m_bob"] + M_BOTTOM + M_HUB + St["m_tube"].sum(1) + It["m_tip"].sum(1)             + (It["T"][:, 0:4].sum(1) / (LEVER * G) if steel else 0.0)
        stress = (It["T"] / (BREAK_MPA * St["area"])).max()
        sprung = slice(8, 12) if steel else slice(0, 12)
        ext = (It["T"][:, sprung] / St["ks"][:, sprung]).max()
        if steel:
            travel = float((np.abs(It["T"][:, :8] / I["T"][:, :8] - 1) * POISE_ARM).max())
        else:
            travel = float(np.percentile(np.abs(It["ell"] - I["ell"]), 95))
        ok = bool(It["ok"].all())
        v += (0 if ok else 3) + max(0, mass.max() - 22) / 2 + max(0, stress - 0.45) / 0.05             + max(0, ext - 40) / 5 + max(0, -It["trim"].min()) / 0.1 + max(0, travel - lim_all) / 3
        out["transpose_" + ("up" if shift > 0 else "down")] = dict(
            ok=ok, mass_kg=float(mass.max()), stress=float(stress), spring_ext_mm=float(ext),
            travel_mm=travel, bob_kg=float(It["m_bob"].mean()), trim_kg=float(It["trim"].mean()))
    T8 = I["T"][:, :8]
    ratio = np.array([2 ** (2 / 12) - 1, 2 ** (-2 / 12) - 1])
    tw = I["T"][:, 8:]
    ext_w = float((tw * 2 ** (2 / 12) / S["ks"][:, 8:]).max())
    if steel:
        one = float(np.abs(ratio).max() * POISE_ARM)                    # one poise, nothing else moves
        ext_g = 0.0
        per_turn = np.full_like(T8, 600 * math.log2(1 + 4 / POISE_ARM))
    else:
        dT = np.zeros((n, 8, 16))
        for i in range(8):
            dT[:, i, 2 * i] = T8[:, i] * ratio[0]
            dT[:, i, 2 * i + 1] = T8[:, i] * ratio[1]
        dy = np.linalg.pinv(I["jac"]) @ dT                              # coupled adjusters, linear
        one = float(np.percentile(np.abs(dy[:, :8, :]).max(1), 95))
        ext_g = float(((T8[:, :, None] + dT) / S["ks"][:, :8, None]).max())
        diag = np.abs(np.einsum("nii->ni", I["jac"][:, :, :8]))
        per_turn = 600 * np.log2(1 + 0.8 * diag / T8)
    worst_turn = float(per_turn.max())
    v += max(0, one - lim_one) / 3 + max(0, max(ext_g, ext_w) - 40) / 5 + max(0, worst_turn - 40) / 10
    out.update(tensioner="steelyard" if steel else "springs", single_string_travel_mm=one,
               cents_per_turn_max=worst_turn, cents_per_turn_median=float(np.median(per_turn)),
               spring_ext_after_mm=max(ext_g, ext_w), violation=float(v))
    return v, out


def pitch(S, T, dist, idx=slice(None)):
    lv = (dist - S["sad"][..., idx]) / 1000
    f1 = np.sqrt(np.maximum(T, 1e-9) / S["mu"][..., idx]) / (2 * lv)
    B = math.pi ** 3 * S["E"][..., idx] * 1e6 * S["d_m"][..., idx] ** 4 / (64 * np.maximum(T, 1e-9) * lv ** 2)
    return f1, B


def evaluate(v, tol, wind_mean=3.0, n=N_MC, detail=False, tune_req=False):
    D = decode(v)
    pen = 0.0
    if D["edge"] / math.sqrt(2) - HUB_R - TIP_R > 310:
        pen += 5 * (D["edge"] / math.sqrt(2) - HUB_R - TIP_R - 310)
    S, I = install(D, tol, n)
    if not I["ok"].all() or (I["T"] <= 0).any():
        return 1e4 if not detail else (1e4, {})
    # install accuracy (should be ~0) and feasibility
    f1, _ = pitch(S, I["T"], I["dist"])
    install_err = np.abs(1200 * np.log2(f1 / S["fstar"])).max()
    pen += 20 * max(0, install_err - 0.5)
    trim = I["trim"]
    pen += 300 * np.maximum(0, -trim).max()
    total_mass = I["m_bob"] + M_BOTTOM + M_HUB + S["m_tube"].sum(1) + I["m_tip"].sum(1)         + (I["T"][:, 0:4].sum(1) / (LEVER * G) if D["tens"] else 0.0)
    pen += 20 * max(0, total_mass.max() - 20)
    stress = (I["T"] / (BREAK_MPA * S["area"])).max()
    pen += 300 * max(0, stress - 0.45)
    sprung = slice(8, 12) if D["tens"] else slice(0, 12)       # steelyards have no spring to overstretch
    ext = (I["T"][:, sprung] / S["ks"][:, sprung]).max()
    pen += 1.0 * max(0, ext - 40)
    ref_S, ref_I = install(D, tolerances(1, zero=True), 1)
    travel = 0.0 if D["tens"] else np.percentile(np.abs(I["ell"] - ref_I["ell"]), 95)
    pen += 2 * max(0, travel - 12)
    tune_info = {}
    if tune_req:
        tv, tune_info = tunability(D, tol, n, S, I)
        pen += 20 * tv

    # operating states
    st = states()
    ns = len(st)
    t_env = st[:, 0][None, :] - D["t_set"]
    asym = np.array([1.0, 0.0, -1.0, 0.0])
    cte = np.broadcast_to(np.asarray(tol["cte"]) if np.ndim(tol["cte"]) else np.full(len(tol["st"]), tol["cte"]),
                          (len(tol["st"]),))[:n]
    if D["frame"]:
        cte = np.broadcast_to(np.asarray(tol["cte_cf"]), (len(tol["st"]),))[:n]
    creep_scale = 0.1 if D["frame"] else 1.0
    dT_tube = t_env[..., None] + (tol["sun"][:n, None, None] * st[:, 1][None, :, None] * asym)
    tipR = S["r0"][:, None, :] * (1 + cte[:, None, None] * dT_tube) \
        - (creep_scale * tol["creep"][:n, None, :] * st[:, 2][None, :, None]) * S["tube_len"]
    ell = I["ell"][:, None, :] * (1 + 12e-6 * t_env[..., None])
    nb = n * ns
    rep = lambda a: np.repeat(a, ns, axis=0)
    tipR_f, ell_f = tipR.reshape(nb, 4), ell.reshape(nb, 12)
    base = dict(tipR=tipR_f, ell=ell_f, k=rep(S["k"]), m_tip=rep(I["m_tip"]), m_tube=rep(S["m_tube"]),
                m_hub=rep(tol["m_hub"][:n]), m_bob=rep(I["m_bob"] + tol["m_bottom"][:n]),
                mid=np.full(nb, S["mid"]))
    drag_b = np.array([0.5 * 1.2 * 100 * 0.9 * 0.031, 0, 0])
    drag_c = np.array([0.5 * 1.2 * 100 * 0.9 * 0.060, 0, 0])
    res = []
    for gust in (0.0, 1.0):
        P = dict(base, drag_c=np.tile(drag_c * gust, (nb, 1)), drag_b=np.tile(drag_b * gust, (nb, 1)))
        x, out, ok = solve(P, rep(I["x"]))
        res.append((out, ok))
    (r0_, tu0, tl0, du0, dl0, _, _), ok0 = res[0]
    (r1_, tu1, tl1, du1, dl1, _, _), ok1 = res[1]
    fail = 1 - (ok0 & ok1).mean()
    pen += 1000 * fail
    span = np.sqrt(tipR_f ** 2 + np.roll(tipR_f, -1, axis=1) ** 2)
    tw = base["k"][:, 8:] * np.maximum(span - ell_f[:, 8:], 0)
    Sb = {k_: rep(v_) for k_, v_ in S.items() if isinstance(v_, np.ndarray) and v_.ndim == 2}
    Sb["fstar"] = np.tile(S["fstar"], (nb, 1))

    u = U_GRID[None, :, None]
    gust_frac = (U_GRID / 10.0) ** 2
    Tg0 = np.concatenate([tu0, tl0], 1)
    Tg1 = np.concatenate([tu1, tl1], 1)
    Tg = Tg0[:, None, :] + (Tg1 - Tg0)[:, None, :] * gust_frac[None, :, None]       # (nb, U, 8)
    T = np.concatenate([Tg, np.repeat(tw[:, None, :], len(U_GRID), 1)], 2)            # (nb, U, 12)
    dist = np.concatenate([du0, dl0, span], 1)[:, None, :]
    lv = (dist - Sb["sad"][:, None, :]) / 1000
    f1 = np.sqrt(T / Sb["mu"][:, None, :]) / (2 * lv)
    Bc = math.pi ** 3 * Sb["E"][:, None, :] * 1e6 * Sb["d_m"][:, None, :] ** 4 / (64 * T * lv ** 2)
    d = Sb["d_m"][:, None, :]
    re = u * d / NU_AIR
    stf = np.repeat(tol["st"][:n], ns)[:, None, None]
    fs = strouhal(re) * u / d * stf
    n0 = fs / f1
    n1 = np.maximum(1, np.round(n0 / np.sqrt(1 + Bc * n0 ** 2)))
    best_n, best_e = n1, np.full(n1.shape, np.inf)
    for dn in (-1, 0, 1):
        cand = np.maximum(1, n1 + dn)
        e = np.abs(np.log(fs / (cand * f1 * np.sqrt(1 + Bc * cand ** 2))))
        better = e < best_e
        best_n, best_e = np.where(better, cand, best_n), np.where(better, e, best_e)
    w = np.exp(-best_e ** 2 / (2 * LOCK_SIGMA ** 2)) * (re >= 47)
    fn = best_n * f1 * np.sqrt(1 + Bc * best_n ** 2)
    pw = weibull_weights(wind_mean)[None, :]

    # E1: fundamentals of the 8 gravity strings against their targets, gust-weighted
    e1_all = np.abs(1200 * np.log2(f1[..., :8] / Sb["fstar"][:, None, :8]))
    e1 = (e1_all.mean(2) * pw).sum(1).mean()
    # E2: intervals between strings of one chord sounding on the same harmonic
    num = den = 0.0
    for base_i in (0, 4):
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = base_i + i, base_i + j
                same = best_n[..., a] == best_n[..., b]
                ww = w[..., a] * w[..., b] * same * pw
                err = np.abs(1200 * np.log2((fn[..., a] / fn[..., b]) / (JUST[i] / JUST[j])))
                num += (ww * err).sum()
                den += ww.sum()
    e2 = num / max(den, 1e-9)
    # coherence: the chord keeps its shape only if a star's strings sound on the same harmonic
    co_num = co_den = 0.0
    for base_i in (0, 4):
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = base_i + i, base_i + j
                ww = w[..., a] * w[..., b] * pw
                co_num += (ww * (best_n[..., a] == best_n[..., b])).sum()
                co_den += ww.sum()
    coherence = co_num / max(co_den, 1e-9)
    # theta: waist i against earth i on the same harmonic
    hit = np.zeros(w.shape[:2])
    for i in range(4):
        a, b = 8 + i, 4 + i
        same = best_n[..., a] == best_n[..., b]
        beat = np.abs(fn[..., a] - fn[..., b])
        both = (w[..., a] * w[..., b] > 0.25) & same
        hit = np.maximum(hit, (both & (beat >= 4) & (beat <= 8)).astype(float))
    sounding = (w[..., :8] > 0.25).any(2).astype(float)
    cover = (hit * pw).sum() / max((sounding * pw).sum(), 1e-9)
    cover_all = (hit * pw).sum() / pw.sum() / nb
    # E3: hub drone against the root earth string
    t_abs = np.repeat(st[:, 0][None, :], n, 0).reshape(-1)
    f_hub = ROOT * np.sqrt((t_abs + 273.15) / (D["t_set"] + 273.15))
    e3 = np.abs(1200 * np.log2(f_hub / f1[:, 0, 4])).mean()

    J = e1 + W_E2 * e2 + 0.15 * e3 + 50 * (1 - cover) + 30 * (1 - coherence) + pen
    if not detail:
        return float(J)
    by_u = (hit.reshape(n, ns, -1).mean((0, 1)), sounding.reshape(n, ns, -1).mean((0, 1)))
    beat_day = []
    i3 = int(np.argmin(np.abs(U_GRID - 3.0)))
    fn_r = fn.reshape(n, ns, len(U_GRID), 12)
    for s_i, (t, day, aged) in enumerate(st):
        if aged == 0 and abs(t - st[st[:, 2] == 0][:, 0].mean()) < 99:
            beat_day.append((float(t), float(np.abs(fn_r[0, s_i, i3, 8] - fn_r[0, s_i, i3, 4]))))
    info = dict(J=float(J), e1_cents=float(e1), e2_cents=float(e2), e3_cents=float(e3), coherence=float(coherence),
                theta_cover_sounding=float(cover), theta_cover_all=float(cover_all), penalty=float(pen),
                install_err_cents=float(install_err), total_mass_kg=float(total_mass.mean()),
                bob_kg=float(I["m_bob"].mean()), trim_kg=float(trim.mean()),
                max_stress=float(stress), max_spring_ext_mm=float(ext), adjuster_travel95_mm=float(travel),
                cents_per_0p1mm=float(1200 * math.log2(1 + S["k"][0, 4] * 0.1 / (2 * I["T"][0, 4]))),
                gust_cents_10ms=float(np.abs(1200 * np.log2(np.sqrt(Tg1 / Tg0)))[:, 4:8].max(1).mean()),
                solve_fail=float(fail), design={k_: (v_.tolist() if isinstance(v_, np.ndarray) else float(v_))
                                                for k_, v_ in D.items()},
                tension_N=I["T"][0].round(2).tolist(), f1_Hz=f1[0, 0].round(2).tolist(),
                theta_by_u=dict(u=U_GRID.tolist(), hit=by_u[0].tolist(), sounding=by_u[1].tolist()),
                beat_vs_temp_3ms=beat_day, tunability=tune_info)
    return float(J), info


# ---------------------------------------------------------------- run
TOL = None
TUNE = True
W_E2 = 0.5                                        # weight on interval error (v4 raises it)


def fitness(v):
    global TOL
    if TOL is None:
        TOL = tolerances(N_MC)
    try:
        return evaluate(v, TOL, tune_req=TUNE)
    except (FloatingPointError, ValueError, np.linalg.LinAlgError):
        return 1e4


def check():
    tol = tolerances(1, zero=True)
    v = REV1.copy()
    D = decode(v)
    S, I = install(D, tol, 1)
    f1, _ = pitch(S, I["T"], I["dist"])
    print("rev 1, no tolerances, vibrating length 609.5 mm")
    print("  tensions N:", I["T"][0].round(2))
    print("  f1 Hz     :", f1[0].round(2))
    print("  targets   :", S["fstar"].round(2))
    print(f"  bob {I['m_bob'][0]:.2f} kg, trim {I['trim'][0]:.3f} kg per tip, converged {I['ok'].all()}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--maxiter", type=int, default=45)
    ap.add_argument("--popsize", type=int, default=10)
    ap.add_argument("--wind-mean", type=float, default=3.0)
    ap.add_argument("--tag", default="v2")
    ap.add_argument("--w-interval", type=float, default=0.5)
    ap.add_argument("--seed-from", default="")
    a = ap.parse_args()
    global W_E2
    W_E2 = a.w_interval
    np.seterr(all="ignore")
    if a.check:
        check()
        return
    tol = tolerances(N_MC)
    fresh = tolerances(48, seed=99)            # out-of-sample check, never seen by the evolution
    t0 = time.time()
    skip = ("theta_by_u", "beat_vs_temp_3ms")
    j1, rev1 = evaluate(REV1.copy(), tol, a.wind_mean, detail=True, tune_req=True)
    rev1b_x = REV1.copy()
    rev1b_x[1] = math.log10(2.0)
    j1b, rev1b = evaluate(rev1b_x, tol, a.wind_mean, detail=True, tune_req=True)
    print(f"rev 1 as drawn: J {j1:.2f}; rev 1 + 2 N/mm springs: J {j1b:.2f}")
    print(json.dumps({k: v for k, v in rev1b.items() if k not in skip}, indent=1))
    sys.stdout.flush()
    history = []

    def cb(xk, convergence=None):
        history.append(np.array(xk, float))
        print(f"gen {len(history)}: best so far x={np.round(xk, 3).tolist()}", flush=True)

    res = differential_evolution(fitness, BOUNDS, integrality=INTEGRALITY, popsize=a.popsize,
                                 maxiter=a.maxiter, tol=1e-4, seed=11, polish=False, updating="deferred",
                                 workers=1, callback=cb, init="sobol", x0=(np.array(json.load(open(a.seed_from))["x"]) if a.seed_from else np.r_[rev1b_x[:14], 1, 1]).clip(
                                     [b[0] for b in BOUNDS], [b[1] for b in BOUNDS]))
    jb, best = evaluate(res.x, tol, a.wind_mean, detail=True, tune_req=True)
    print(f"evolved: J {jb:.2f}")
    print(json.dumps({k: v for k, v in best.items() if k not in skip}, indent=1))
    _, val_b = evaluate(res.x, fresh, a.wind_mean, n=48, detail=True, tune_req=True)
    _, val_1b = evaluate(rev1b_x, fresh, a.wind_mean, n=48, detail=True, tune_req=True)
    sens = {}
    for wm in (2.0, 4.0):
        sens[str(wm)] = dict(evolved=evaluate(res.x, tol, wm, detail=True, tune_req=True)[1]["theta_cover_sounding"],
                             rev1b=evaluate(rev1b_x, tol, wm, detail=True, tune_req=True)[1]["theta_cover_sounding"])
    print("validation (48 fresh samples):",
          {k: round(val_b[k], 3) for k in ("J", "e1_cents", "e2_cents", "theta_cover_sounding")},
          {k: round(val_1b[k], 3) for k in ("J", "e1_cents", "e2_cents", "theta_cover_sounding")})
    rows = []
    for g_i, xk in enumerate(history, 1):
        jk, ik = evaluate(xk, tol, a.wind_mean, detail=True, tune_req=True)
        rows.append((g_i, jk, ik["e1_cents"], ik["e2_cents"], ik["theta_cover_sounding"], ik["coherence"]))
    with open(f"evolution_history_{a.tag}.csv", "w") as fh:
        fh.write("generation,fitness,fund_err_cents,interval_err_cents,theta_cover,coherence\n")
        for r in rows:
            fh.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    with open(f"evolution_result_{a.tag}.json", "w") as fh:
        json.dump(dict(w_interval=W_E2, rev1=rev1, rev1b=rev1b, evolved=best, validation=dict(evolved=val_b, rev1b=val_1b),
                       wind_sensitivity=sens, x=res.x.tolist(), nfev=int(res.nfev), wind_mean=a.wind_mean,
                       n_mc=N_MC, seconds=time.time() - t0), fh, indent=1)
    print(f"done in {time.time() - t0:.0f} s, {res.nfev} evaluations")


if __name__ == "__main__":
    main()
