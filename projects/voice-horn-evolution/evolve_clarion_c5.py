"""Forge cycle 5 for the Clarion family - PROJECTOR LANE, tier-1.5 directivity correction.

Change vs cycle 3: the on-axis directivity uses the tier-2-derived cap

    Q_eff = min( Q_piston, Q_cone ),   Q_cone = 2 / (1 - cos(theta_exit))

with theta_exit the bore wall angle at the mouth. Tier-2 FEM (tier2_knee_c4.py) showed the
uniform-piston Q overestimates DI by up to 9 dB above 2 kHz for wide mouths; the flare-cone
coverage cap reproduces the observed saturation within ~2 dB. Everything else is cycle 3:
same chain (cup, tube, step, bore), glottal source, genes, constraints, objectives, budget.

Calibration printed: corrected tier-1 G_speech of the cycle-4 KNEE vs tier-2's measured value.
Prediction under test: Q_cone rewards gentle exit angles, so evolution should abandon
steep-exit profiles (c3 chose n ~ 2.5) and move toward tractrix-character bores.

Run from this directory: python evolve_clarion_c5.py
Outputs: pareto_c5.csv, evolution_c5.png
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


def evaluate_chain(Sseg, dseg, a_m, r_h, Ltot, re_mm, rt_mm, L_mm):
    N, S = Sseg.shape
    x = 2 * K[None, :] * a_m
    R1 = 1 - 2 * j1(x) / x
    X1 = 2 * struve(1, x) / x
    Sm = np.pi * a_m ** 2
    Zrad = RHO * C / Sm * (R1 + 1j * X1)
    # tier-1.5: exit wall angle from the last bore segment; flare-cone directivity cap
    slope = (r_h[:, -1:] - r_h[:, -2:-1]) / (L_mm[:, None] / NH)  # mm/mm, (N,1)
    theta_exit = np.arctan(slope)                                  # (N,1)
    Qcone = 2.0 / np.maximum(1.0 - np.cos(theta_exit), 1e-6)
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
    Qp = (K[None, :] * a_m) ** 2 / np.maximum(R1, 1e-9)      # piston directivity factor
    Q = np.minimum(Qp, Qcone)                                 # tier-1.5 cap [from tier 2]
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
    return F, cv, G, theta_exit[:, 0]


def evaluate(X):
    Sseg, dseg, a_m, r_h, Ltot = chain_from_genes(X)
    return evaluate_chain(Sseg, dseg, a_m, r_h, Ltot, X[:, 0], X[:, 2], X[:, 4])


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
F, cv, _, _ = evaluate(X)
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
    Fk, cvk, _, _ = evaluate(Xk)
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
_, _, _, TH = evaluate(Xp)
TH = np.degrees(TH)
print(f"cycle 5 front: {len(pf)} feasible of {POP}")
print(f"G range: {(-Fp[:,0]).min():.1f} .. {(-Fp[:,0]).max():.1f} dB | "
      f"Gam: {Fp[:,1].min():.3f} .. {Fp[:,1].max():.3f}")
print(f"flare exponent n on the front: {Xp[:,5].min():.2f} .. {Xp[:,5].max():.2f}, "
      f"median {np.median(Xp[:,5]):.2f}   (cycle 3 chose ~2.5-2.9)")
print(f"exit wall angle: {TH.min():.0f} .. {TH.max():.0f} deg, median {np.median(TH):.0f} deg")
print(f"members with n < 1 (tractrix-character): {(Xp[:,5] < 1).sum()}/{len(pf)}")

# calibration: cycle-4 KNEE under the corrected evaluator vs tier-2 FEM truth
XC4 = np.array([[20.0, 24.59, 24.12, 154.1, 451.25, 2.53, -0.05, 0.67]])
Fc4, cvc4, Gc4, thc4 = evaluate(XC4)
G_T2 = np.array([24.8, 19.0, 18.5, 16.1, 15.6, 16.6, 15.7, 15.6, 11.8, 10.3, 7.8, 5.6])
Gs_t2 = (G_T2 * WF).sum()
print(f"CALIBRATION - c4 KNEE: corrected tier-1 G_speech = {-Fc4[0,0]:.1f} dB | "
      f"tier-2 FEM = {Gs_t2:.1f} dB | uncorrected tier-1 was 20.1 dB")

# c3 KNEE rescored under correction
XC3 = np.array([[19.98, 31.34, 21.38, 155.0, 484.17, 2.87, -0.01, 0.75]])
Fc3, _, _, _ = evaluate(XC3)
print(f"c3 KNEE rescored under correction: G = {-Fc3[0,0]:.1f} dB (was 19.7 uncorrected)")

# hand Clarion full chain under correction
Xhand = np.array([[11.0, 30.0, 26.5, 125.0, 280.0, 1.0, 0.0, 0.5]])
Sseg_h, dseg_h, a_h, rh_h, Lt_h = chain_from_genes(Xhand)
r_exp = 26.5 * (125.0 / 26.5) ** UH
Sh = np.pi * ((0.5 * (r_exp[:-1] + r_exp[1:])) / 1000.0) ** 2
Sseg_h[:, NC + NT:] = Sh[None, :]
Fh, cvh, Gh, thh = evaluate_chain(Sseg_h, dseg_h, a_h, r_exp[None, :], Lt_h,
                                  np.array([11.0]), np.array([26.5]), np.array([280.0]))
print(f"hand Clarion (corrected, full chain): G = {-Fh[0,0]:.1f} dB, cv = {cvh[0]:.2f}")

gmax = pf[np.argmax(-Fp[:, 0])]
cmin = pf[np.argmin(Fp[:, 2])]
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
knee = pf[np.argmax(np.minimum(gn, rn))]
picks = [("MAX-GAIN", gmax, "#7a1f1f"), ("KNEE", knee, "#B45309"), ("MIN-COST", cmin, "#1D4ED8")]
for nm, i, _ in picks:
    _, _, _, th_i = evaluate(X[i][None, :])
    print(f"{nm} genes: {np.round(X[i],2)}  G={-F[i,0]:.1f} Gam={F[i,1]:.3f} "
          f"Pcost={F[i,2]:.2f} n={X[i,5]:.2f} exit={np.degrees(th_i[0]):.0f}deg")

np.savetxt("pareto_c5.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:], TH[:, None]]),
           delimiter=",", header="re_mm,Le_mm,rt_mm,rm_mm,L_mm,n,bulgeA,bulgeX,"
           "G_speech_dB,refl,print_cost,exit_deg", comments="")

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
ax.plot(Fc4[0, 1], -Fc4[0, 0], "X", ms=10, color="#6D28D9", alpha=0.85)
ax.annotate("c4 KNEE, corrected", (Fc4[0, 1], -Fc4[0, 0]),
            (Fc4[0, 1] + 0.03, -Fc4[0, 0] - 1.4), fontsize=7.5, color="#6D28D9",
            arrowprops=dict(arrowstyle="->", color="#6D28D9", lw=0.7))
if cvh[0] == 0:
    ax.plot(Fh[0, 1], -Fh[0, 0], "*", ms=17, color="#B45309", mec=INK, mew=0.8)
    ax.annotate("hand CLARION (corrected)", (Fh[0, 1], -Fh[0, 0]),
                (Fh[0, 1] + 0.04, -Fh[0, 0] - 1.6), fontsize=7.5,
                arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.7))
for nm, i, col in picks:
    ax.plot(F[i, 1], -F[i, 0], "o", ms=9, mfc="none", mec=col, mew=2)
    ax.annotate(nm, (F[i, 1], -F[i, 0]), (F[i, 1] + 0.012, -F[i, 0] + 0.4),
                fontsize=7.5, color=col)
ax.set_xlabel("mean reflection at the lip cup")
ax.set_ylabel("speech-weighted forward gain (dB, corrected)")
ax.set_title("cycle-5 front: projector lane under the tier-1.5 lens", fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[0, 1])
for nm, i, col in picks:
    re_, Le, rt_, rm_, L_, n_, A_, x0_ = X[i]
    rr = rt_ + (rm_ - rt_) * UH ** n_ + A_ * (rm_ - rt_) * np.exp(-((UH - x0_) / 0.12) ** 2) * np.sin(np.pi * UH)
    ax.plot(UH * L_, rr, color=col, lw=2, label=f"{nm} (n={n_:.2f}, exit {np.degrees(np.arctan((rr[-1]-rr[-2])/(L_/NH))):.0f} deg)")
    ax.plot(UH * L_, -rr, color=col, lw=2)
rr4 = 24.12 + (154.1 - 24.12) * UH ** 2.53 + (-0.05) * (154.1 - 24.12) * np.exp(-((UH - 0.67) / 0.12) ** 2) * np.sin(np.pi * UH)
ax.plot(UH * 451.25, rr4, color="#6D28D9", lw=1.4, ls="--", label="c4 KNEE (n=2.53, dashed)")
ax.plot(UH * 451.25, -rr4, color="#6D28D9", lw=1.4, ls="--")
ax.set_aspect("equal")
ax.set_xlabel("axial position (mm)"); ax.set_ylabel("radius (mm)")
ax.legend(fontsize=6.6, loc="upper left", frameon=False)
ax.set_title("did the correction bend the bores toward tractrix character?", fontsize=9.5,
             loc="left")

ax = fig.add_subplot(gs[1, 0])
for nm, i, col in picks:
    _, _, Gi, _ = evaluate(X[i][None, :])
    ax.semilogx(FREQS, Gi[0], color=col, lw=1.9, marker="o", ms=3.5, label=nm)
ax.semilogx(FREQS, Gc4[0], color="#6D28D9", lw=1.6, ls="--", marker="s", ms=3.5,
            label="c4 KNEE, corrected tier-1")
ax.semilogx(FREQS, G_T2, color=INK, lw=0, marker="D", ms=5, label="c4 KNEE, tier-2 FEM (truth)")
ax.axvspan(1000, 4000, color="#B45309", alpha=0.07)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("forward gain (dB)")
ax.legend(fontsize=7, frameon=False)
ax.set_title("calibration: corrected tier-1 (dashed) against the FEM measurements (diamonds)",
             fontsize=9, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
sc = ax.scatter(Xp[:, 3], TH, c=-Fp[:, 0], cmap="magma_r", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("G_speech (dB)")
try:
    c3f = np.loadtxt("pareto_c3.csv", delimiter=",", skiprows=1)
    n3, rm3, L3 = c3f[:, 5], c3f[:, 3], c3f[:, 4]
    th3 = np.degrees(np.arctan((rm3 - c3f[:, 2]) * n3 / L3))
    ax.scatter(rm3, th3, s=20, facecolors="none", edgecolors="#8a8578", lw=0.8,
               label="cycle-3 front (uncorrected)")
    ax.legend(fontsize=7.2, frameon=False, loc="upper left")
except OSError:
    pass
ax.set_xlabel("mouth radius r_m (mm)")
ax.set_ylabel("exit wall angle (deg)")
ax.set_title("the collimation move: exit angle vs mouth,\ncycle 5 (filled) vs cycle 3 (hollow)",
             fontsize=9, loc="left")
ax.grid(alpha=0.15)

fig.suptitle("THE FORGE, cycle 5  -  projector lane under the tier-1.5 directivity cap\n"
             f"Q_eff = min(Q_piston, 2/(1-cos theta_exit)); {POP*NGEN+POP} evaluations, seed 7",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("evolution_c5.png", dpi=150)
print("wrote pareto_c5.csv, evolution_c5.png")
