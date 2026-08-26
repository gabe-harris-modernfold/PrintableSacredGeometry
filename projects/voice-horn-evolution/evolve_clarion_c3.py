"""Forge cycle 3 for the Clarion family. Change vs cycle 2 (evolve_clarion_c2.py):

The epilarynx chain is modeled in. The transfer-matrix chain is now
  glottal Norton source (Z_s = 4 MPa*s/m^3)
    -> lip cup (fixed: dia 45 -> tube, 40 mm linear taper)
    -> epilaryngeal tube (GENES: radius r_e, length L_e)  <- the squillo element, S15.11
    -> area step into the horn throat (ratio A_t/A_e now VISIBLE to the optimizer)
    -> horn bore (genes r_t, r_m, L, n, bulge A, x0)  -> piston radiation load.

8 genes; new constraint r_t >= r_e (no expansion-then-choke). Reflection objective is
now measured at the lip-cup entry - where the mouth actually is.
Objectives (3): -G_speech, Gamma_cup, P_cost = vol_L + (segments-1).

Run from this directory: python evolve_clarion_c3.py
Outputs: pareto_c3.csv, evolution_c3.png
"""
import numpy as np
from scipy.special import j1, struve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6
FREQS = np.array([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000], float)
WF = np.array([0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2])
WF = WF / WF.sum()
K = 2 * np.pi * FREQS / C
R_CUP, L_CUP, NC = 22.5, 40.0, 16
NT, NH = 12, 140
#        r_e    L_e    r_t    r_m    L      n     A      x0
XL = np.array([8.0, 10.0, 10.0, 60.0, 120.0, 0.35, -0.30, 0.20])
XU = np.array([20.0, 60.0, 40.0, 155.0, 580.0, 4.00, 0.30, 0.80])
RNG = np.random.default_rng(7)
UH = np.linspace(0, 1, NH + 1)


def chain_from_genes(X):
    """Returns segment areas Sseg (N,S) m^2, segment lengths dseg (N,S) m,
    mouth radius a_m (N,1) m, horn radii r_h (N,NH+1) mm, total length (N,) mm."""
    re_, Le, rt, rm, L, n, A, x0 = [X[:, i:i + 1] for i in range(8)]
    uc = np.linspace(0, 1, NC + 1)[None, :]
    r_cup = R_CUP + (re_ - R_CUP) * uc
    r_h = rt + (rm - rt) * UH[None, :] ** n \
        + A * (rm - rt) * np.exp(-((UH[None, :] - x0) / 0.12) ** 2) * np.sin(np.pi * UH[None, :])
    r_h = np.maximum(r_h, 1.0)
    mid = lambda r: 0.5 * (r[:, :-1] + r[:, 1:])
    S_cup = np.pi * (mid(r_cup) / 1000.0) ** 2
    S_tub = np.pi * np.repeat(re_ / 1000.0, NT, axis=1) ** 2
    S_hrn = np.pi * (mid(r_h) / 1000.0) ** 2
    Sseg = np.concatenate([S_cup, S_tub, S_hrn], axis=1)
    d_cup = np.full_like(S_cup, L_CUP / NC / 1000.0)
    d_tub = np.repeat(Le / NT / 1000.0, NT, axis=1)
    d_hrn = np.repeat(L / NH / 1000.0, NH, axis=1)
    dseg = np.concatenate([d_cup, d_tub, d_hrn], axis=1)
    Ltot = (L_CUP + Le + L)[:, 0]
    return Sseg, dseg, r_h[:, -1:] / 1000.0, r_h, Ltot


def evaluate_chain(Sseg, dseg, a_m, r_h, Ltot, re_mm, rt_mm):
    N, S = Sseg.shape
    x = 2 * K[None, :] * a_m
    R1 = 1 - 2 * j1(x) / x
    X1 = 2 * struve(1, x) / x
    Sm = np.pi * a_m ** 2
    Zrad = RHO * C / Sm * (R1 + 1j * X1)
    p = Zrad.astype(complex)
    Uv = np.ones_like(p)
    for s in range(S - 1, -1, -1):
        Zc = RHO * C / Sseg[:, s:s + 1]
        kd = K[None, :] * dseg[:, s:s + 1]
        cs, sn = np.cos(kd), np.sin(kd)
        p, Uv = cs * p + 1j * Zc * sn * Uv, 1j * sn / Zc * p + cs * Uv
    Zin = p / Uv
    src = ZS / (ZS + Zin)
    W = 0.5 * np.real(Zrad) * np.abs(src / Uv) ** 2
    Wref = 0.5 * RHO * C * K ** 2 / (4 * np.pi)
    Q = (K[None, :] * a_m) ** 2 / np.maximum(R1, 1e-9)
    G = 10 * np.log10(np.maximum(W / Wref[None, :], 1e-12)) \
        + 10 * np.log10(np.maximum(Q / 2.0, 1e-9))
    Gs = (G * WF[None, :]).sum(1)
    Scup = np.pi * (R_CUP / 1000.0) ** 2
    Gam = np.abs((Zin - RHO * C / Scup) / (Zin + RHO * C / Scup)).mean(1)
    r_eq = np.sqrt(Sseg / np.pi) * 1000.0                     # mm, per segment
    dr = np.diff(r_eq, axis=1)
    dxm = 1000.0 * 0.5 * (dseg[:, :-1] + dseg[:, 1:])
    Alat = (2 * np.pi * 0.5 * (r_eq[:, :-1] + r_eq[:, 1:])
            * np.sqrt(dxm ** 2 + dr ** 2)).sum(1)
    vol_L = Alat * 5.0 / 1e6
    segs = np.ceil((Ltot + 20.0) / 290.0)
    Pcost = vol_L + (segs - 1.0)
    ripple = G[:, 2:].max(1) - G[:, 2:].min(1)
    choke = np.maximum(0.6 * rt_mm, 8.0) - r_h.min(1)
    stepv = re_mm - rt_mm                                     # r_t >= r_e
    cv = np.maximum(0, ripple - 12.0) + np.maximum(0, choke) + np.maximum(0, stepv)
    F = np.stack([-Gs, Gam, Pcost], 1)
    return F, cv, G


def evaluate(X):
    Sseg, dseg, a_m, r_h, Ltot = chain_from_genes(X)
    return evaluate_chain(Sseg, dseg, a_m, r_h, Ltot, X[:, 0], X[:, 2])


def nds(F, cv):
    N = len(F)
    dominates = np.zeros((N, N), bool)
    for i in range(N):
        ci = cv[i]
        both_feas = (ci == 0) & (cv == 0)
        dom_f = (F[i][None, :] <= F).all(1) & (F[i][None, :] < F).any(1)
        dominates[i] = ((ci < cv) & (cv > 0)) | (both_feas & dom_f) | ((ci == 0) & (cv > 0))
    fronts, rank = [], np.full(N, -1)
    n_dom = dominates.sum(0)
    cur = np.where(n_dom == 0)[0].tolist()
    fr = 0
    while cur:
        fronts.append(cur)
        rank[cur] = fr
        nxt = []
        for i in cur:
            for j in np.where(dominates[i])[0]:
                n_dom[j] -= 1
                if n_dom[j] == 0:
                    nxt.append(j)
        cur, fr = nxt, fr + 1
    return fronts, rank


def crowding(F, idx):
    n, m = len(idx), F.shape[1]
    dist = np.zeros(n)
    for k in range(m):
        o = np.argsort(F[idx, k])
        f = F[idx, k][o]
        span = max(f[-1] - f[0], 1e-12)
        dist[o[0]] = dist[o[-1]] = np.inf
        dist[o[1:-1]] += (f[2:] - f[:-2]) / span
    return dist


def sbx(a, b, eta=15):
    u = RNG.random(a.shape)
    beta = np.where(u <= 0.5, (2 * u) ** (1 / (eta + 1)), (1 / (2 * (1 - u))) ** (1 / (eta + 1)))
    return (np.clip(0.5 * ((1 + beta) * a + (1 - beta) * b), XL, XU),
            np.clip(0.5 * ((1 - beta) * a + (1 + beta) * b), XL, XU))


def polymut(x, eta=20, pm=1 / 8):
    y = x.copy()
    mask = RNG.random(x.shape) < pm
    u = RNG.random(x.shape)
    delta = np.where(u < 0.5, (2 * u) ** (1 / (eta + 1)) - 1, 1 - (2 * (1 - u)) ** (1 / (eta + 1)))
    return np.clip(np.where(mask, y + delta * (XU - XL), y), XL, XU)


POP, NGEN = 64, 75
X = XL + RNG.random((POP, 8)) * (XU - XL)
F, cv, _ = evaluate(X)
trace = []
for gen in range(NGEN):
    fronts, rank = nds(F, cv)
    crowd = np.zeros(len(F))
    for fr in fronts:
        crowd[fr] = crowding(F, fr)
    def pick():
        i, j = RNG.integers(0, len(F), 2)
        if rank[i] != rank[j]:
            return i if rank[i] < rank[j] else j
        return i if crowd[i] > crowd[j] else j
    kids = []
    while len(kids) < POP:
        c1, c2 = sbx(X[pick()], X[pick()])
        kids.append(polymut(c1)); kids.append(polymut(c2))
    Xk = np.array(kids[:POP])
    Fk, cvk, _ = evaluate(Xk)
    Xall, Fall, cvall = np.vstack([X, Xk]), np.vstack([F, Fk]), np.concatenate([cv, cvk])
    fronts, rank_all = nds(Fall, cvall)
    crowd_all = np.zeros(len(Fall))
    for fr in fronts:
        crowd_all[fr] = crowding(Fall, fr)
    order = sorted(range(len(Fall)), key=lambda i: (rank_all[i], -crowd_all[i]))
    keep = np.array(order[:POP])
    X, F, cv = Xall[keep], Fall[keep], cvall[keep]
    feas = cv == 0
    trace.append((-F[feas, 0]).max() if feas.any() else np.nan)

fronts, rank = nds(F, cv)
pf = np.array(fronts[0]); pf = pf[cv[pf] == 0]
Fp, Xp = F[pf], X[pf]
step_ratio = (Xp[:, 2] / Xp[:, 0]) ** 2
print(f"cycle 3 front: {len(pf)} feasible of {POP}")
print(f"G range: {(-Fp[:,0]).min():.1f} .. {(-Fp[:,0]).max():.1f} dB | "
      f"Gam: {Fp[:,1].min():.3f} .. {Fp[:,1].max():.3f}")
print(f"epilarynx tube on the front: r_e {Xp[:,0].min():.1f}..{Xp[:,0].max():.1f} mm, "
      f"L_e {Xp[:,1].min():.0f}..{Xp[:,1].max():.0f} mm")
print(f"area step A_t/A_e on the front: {step_ratio.min():.1f} .. {step_ratio.max():.1f}, "
      f"median {np.median(step_ratio):.1f}")
print(f"members with step >= 6 (squillo decoupling): {(step_ratio >= 6).sum()}/{len(pf)}")

# hand Clarion, FULL chain (Rev-1 geometry): cup 22.5->11, tube r11 x 30, exp horn 26.5->125 x 280
Xhand = np.array([[11.0, 30.0, 26.5, 125.0, 280.0, 1.0, 0.0, 0.5]])
Sseg_h, dseg_h, a_h, rh_h, Lt_h = chain_from_genes(Xhand)
r_exp = 26.5 * (125.0 / 26.5) ** UH
Sh = np.pi * ((0.5 * (r_exp[:-1] + r_exp[1:])) / 1000.0) ** 2
Sseg_h[:, NC + NT:] = Sh[None, :]
Fh, cvh, Gh = evaluate_chain(Sseg_h, dseg_h, a_h, r_exp[None, :], Lt_h, np.array([11.0]),
                             np.array([26.5]))
print(f"hand Clarion FULL CHAIN: G={-Fh[0,0]:.1f} dB (bare-bore c2 rescore was 16.1) "
      f"Gam={Fh[0,1]:.3f} cv={cvh[0]:.2f} step={(26.5/11)**2:.1f}")
domd = ((Fp <= Fh[0]).all(1) & (Fp < Fh[0]).any(1)).sum()
print(f"front members dominating hand Clarion: {domd}")

gmax = pf[np.argmax(-Fp[:, 0])]
cmin = pf[np.argmin(Fp[:, 2])]
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
knee = pf[np.argmax(np.minimum(gn, rn))]
picks = [("MAX-GAIN", gmax, "#7a1f1f"), ("KNEE", knee, "#B45309"), ("MIN-COST", cmin, "#1D4ED8")]
for nm, i, _ in picks:
    sr = (X[i, 2] / X[i, 0]) ** 2
    print(f"{nm} genes: {np.round(X[i],2)}  G={-F[i,0]:.1f} Gam={F[i,1]:.3f} "
          f"Pcost={F[i,2]:.2f} step={sr:.1f} f_qw~{C/(4*X[i,1]/1000)/1000:.1f}kHz")

np.savetxt("pareto_c3.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:], step_ratio[:, None]]),
           delimiter=",", header="re_mm,Le_mm,rt_mm,rm_mm,L_mm,n,bulgeA,bulgeX,"
           "G_speech_dB,refl,print_cost,step_ratio", comments="")

# ------------------------------- figure -------------------------------
BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
INK = "#1F2937"
fig = plt.figure(figsize=(11.8, 8.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28, left=0.07, right=0.97,
                      top=0.87, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(Fp[:, 1], -Fp[:, 0], c=Fp[:, 2], cmap="viridis", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("print cost (L + extra segments)")
ax.plot(Fh[0, 1], -Fh[0, 0], "*", ms=18, color="#B45309", mec=INK, mew=0.8)
ax.annotate("hand CLARION, full chain\n(cup + tube + 1:5.8 step)", (Fh[0, 1], -Fh[0, 0]),
            (Fh[0, 1] + 0.05, -Fh[0, 0] - 3.0), fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.8))
for nm, i, col in picks:
    ax.plot(F[i, 1], -F[i, 0], "o", ms=9, mfc="none", mec=col, mew=2)
    ax.annotate(nm, (F[i, 1], -F[i, 0]), (F[i, 1] + 0.015, -F[i, 0] + 0.7),
                fontsize=7.5, color=col)
ax.set_xlabel("mean reflection at the LIP CUP (phonation-load proxy)")
ax.set_ylabel("speech-weighted forward gain (dB)")
ax.set_title("cycle-3 front: the epilarynx chain is in the model", fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[0, 1])
for nm, i, col in picks:
    re_, Le, rt, rm, L, n, A, x0 = X[i]
    xs_c = np.linspace(-L_CUP - Le, -Le, NC + 1)
    ax.plot(xs_c, np.linspace(R_CUP, re_, NC + 1), color=col, lw=1.6)
    ax.plot(xs_c, -np.linspace(R_CUP, re_, NC + 1), color=col, lw=1.6)
    ax.plot([-Le, 0], [re_, re_], color=col, lw=1.6)
    ax.plot([-Le, 0], [-re_, -re_], color=col, lw=1.6)
    rr = rt + (rm - rt) * UH ** n + A * (rm - rt) * np.exp(-((UH - x0) / 0.12) ** 2) * np.sin(np.pi * UH)
    ax.plot([0, 0], [re_, rr[0]], color=col, lw=1.6)
    ax.plot([0, 0], [-re_, -rr[0]], color=col, lw=1.6)
    ax.plot(UH * L, rr, color=col, lw=2, label=f"{nm} (re={re_:.0f}, Le={Le:.0f}, step {((rt/re_)**2):.1f})")
    ax.plot(UH * L, -rr, color=col, lw=2)
ax.plot(UH * 280, r_exp, color="#B45309", lw=1.5, ls="--", label="hand CLARION")
ax.plot(UH * 280, -r_exp, color="#B45309", lw=1.5, ls="--")
ax.set_aspect("equal")
ax.set_xlabel("axial position (mm; device throat at 0)")
ax.set_ylabel("radius (mm)")
ax.legend(fontsize=6.6, loc="upper left", frameon=False)
ax.set_title("full chains: cup, tube, step, bore", fontsize=9.5, loc="left")

ax = fig.add_subplot(gs[1, 0])
for nm, i, col in picks:
    _, _, Gi = evaluate(X[i][None, :])
    ax.semilogx(FREQS, Gi[0], color=col, lw=1.9, marker="o", ms=3.5, label=nm)
ax.semilogx(FREQS, Gh[0], color="#B45309", lw=1.6, ls="--", marker="s", ms=3.5,
            label="hand CLARION")
ax.axvspan(2400, 3400, color="#B45309", alpha=0.10)
ax.text(2500, ax.get_ylim()[0] + 1.0, "squillo band", fontsize=7.5, color="#B45309")
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("forward gain vs bare voice (dB)")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("third-octave gain with the tube in the chain", fontsize=9.5, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
sc = ax.scatter(Xp[:, 1], step_ratio, c=-Fp[:, 0], cmap="magma_r", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("G_speech (dB)")
ax.axhline(6, color="#2c5c3b", lw=1.2, ls="--")
ax.text(12, 6.4, "1:6 decoupling threshold (S0.4)", fontsize=7.6, color="#2c5c3b")
ax.plot(30, (26.5 / 11) ** 2, "*", ms=16, color="#B45309", mec=INK, mew=0.8)
ax.text(31.5, (26.5 / 11) ** 2 - 0.9, "hand CLARION", fontsize=7.5, color="#B45309")
ax.set_xlabel("epilarynx tube length L_e (mm)")
ax.set_ylabel("area step  A_t / A_e")
ax.set_title("what evolution decided about the squillo:\ntube length vs decoupling step, colored by gain",
             fontsize=9, loc="left")
ax.grid(alpha=0.15)

fig.suptitle("THE FORGE, cycle 3  -  the epilarynx chain modeled in\n"
             f"8 genes, {POP*NGEN+POP} tier-1 evaluations, seed 7; source at the lip cup",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("evolution_c3.png", dpi=150)
print("wrote pareto_c3.csv, evolution_c3.png")
