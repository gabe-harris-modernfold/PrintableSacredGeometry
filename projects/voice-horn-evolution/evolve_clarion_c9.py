"""Forge cycle 9 - CLARION deepening pass. Richer bore + triple budget.

Bore becomes a 3-segment piecewise-exponential (each segment its own flare constant):
genes [r_e, L_e, r_t, r_a, r_b, r_m, L] (7). Chain, glottal source, tier-1.5 directivity
cap, objectives and constraints identical to cycle 5. Pop 96 x 120 gens, seed 7.
Success test: dominate or push past the cycle-5 KNEE (G 18.1 dB, Gam 0.231, Pcost 2.10).
Run from this directory -> pareto_c9.csv, evolution_c9.png
"""
import numpy as np
from scipy.special import j1, struve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6
FREQS = np.array([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000], float)
WF = np.array([0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2]); WF /= WF.sum()
K = 2 * np.pi * FREQS / C
R_CUP, L_CUP, NC = 22.5, 40.0, 16
NT, NH = 12, 141   # NH divisible by 3
UH = np.linspace(0, 1, NH + 1)
#              r_e   L_e   r_t   r_a    r_b    r_m    L
XL = np.array([8.0, 10.0, 13.0, 13.0, 13.0, 60.0, 200.0])
XU = np.array([20.0, 60.0, 40.0, 160.0, 160.0, 155.0, 600.0])
RNG = np.random.default_rng(7)


def bore(X):
    rt, ra, rb, rm = X[:, 2:3], X[:, 3:4], X[:, 4:5], X[:, 5:6]
    u = UH[None, :]
    r = np.where(u <= 1 / 3, rt * (ra / rt) ** (3 * u),
        np.where(u <= 2 / 3, ra * (rb / ra) ** (3 * u - 1),
                 rb * (rm / rb) ** (3 * u - 2)))
    return np.maximum(r, 1.0)


def evaluate(X):
    re_, Le, L = X[:, 0:1], X[:, 1:2], X[:, 6:7]
    r_h = bore(X)
    mid = lambda r: 0.5 * (r[:, :-1] + r[:, 1:])
    r_cup = R_CUP + (re_ - R_CUP) * np.linspace(0, 1, NC + 1)[None, :]
    Sseg = np.concatenate([np.pi * (mid(r_cup) / 1000) ** 2,
                           np.pi * np.repeat(re_ / 1000, NT, 1) ** 2,
                           np.pi * (mid(r_h) / 1000) ** 2], 1)
    dseg = np.concatenate([np.full((len(X), NC), L_CUP / NC / 1000),
                           np.repeat(Le / NT / 1000, NT, 1),
                           np.repeat(L / NH / 1000, NH, 1)], 1)
    a_m = r_h[:, -1:] / 1000
    x = 2 * K[None, :] * a_m
    R1 = 1 - 2 * j1(x) / x
    X1 = 2 * struve(1, x) / x
    Zrad = RHO * C / (np.pi * a_m ** 2) * (R1 + 1j * X1)
    slope = (r_h[:, -1:] - r_h[:, -2:-1]) / (L / NH)
    Qcone = 2.0 / np.maximum(1 - np.cos(np.arctan(slope)), 1e-6)
    Q = np.minimum((K[None, :] * a_m) ** 2 / np.maximum(R1, 1e-9), Qcone)
    p, Uv = Zrad.astype(complex), np.ones_like(Zrad, dtype=complex)
    for s in range(Sseg.shape[1] - 1, -1, -1):
        Zc = RHO * C / Sseg[:, s:s + 1]
        kd = K[None, :] * dseg[:, s:s + 1]
        cs, sn = np.cos(kd), np.sin(kd)
        p, Uv = cs * p + 1j * Zc * sn * Uv, 1j * sn / Zc * p + cs * Uv
    Zin = p / Uv
    src = ZS / (ZS + Zin)
    W = 0.5 * np.real(Zrad) * np.abs(src / Uv) ** 2
    Wref = 0.5 * RHO * C * K ** 2 / (4 * np.pi)
    G = 10 * np.log10(np.maximum(W / Wref[None, :], 1e-12)) + 10 * np.log10(Q / 2)
    Gs = (G * WF[None, :]).sum(1)
    Acup = np.pi * (R_CUP / 1000) ** 2
    Gam = np.abs((Zin - RHO * C / Acup) / (Zin + RHO * C / Acup)).mean(1)
    dr = np.diff(np.sqrt(Sseg / np.pi) * 1000, axis=1)
    dxm = 1000 * 0.5 * (dseg[:, :-1] + dseg[:, 1:])
    Alat = (2 * np.pi * 0.5 * (np.sqrt(Sseg[:, :-1] / np.pi) + np.sqrt(Sseg[:, 1:] / np.pi))
            * 1000 * np.sqrt(dxm ** 2 + dr ** 2)).sum(1)
    Pc = Alat * 5 / 1e6 + (np.ceil((L_CUP + Le[:, 0] + L[:, 0] + 20) / 290) - 1)
    ripple = G[:, 2:].max(1) - G[:, 2:].min(1)
    cv = np.maximum(0, ripple - 12) + np.maximum(0, np.maximum(0.6 * X[:, 2], 8) - r_h.min(1)) \
        + np.maximum(0, X[:, 0] - X[:, 2])
    return np.stack([-Gs, Gam, Pc], 1), cv, G


def nds(F, cv):
    N = len(F)
    dom = np.zeros((N, N), bool)
    for i in range(N):
        ci = cv[i]
        bf = (ci == 0) & (cv == 0)
        df = (F[i][None] <= F).all(1) & (F[i][None] < F).any(1)
        dom[i] = ((ci < cv) & (cv > 0)) | (bf & df) | ((ci == 0) & (cv > 0))
    fronts, rank = [], np.full(N, -1)
    nd = dom.sum(0)
    cur = np.where(nd == 0)[0].tolist(); fr = 0
    while cur:
        fronts.append(cur); rank[cur] = fr
        nxt = []
        for i in cur:
            for j in np.where(dom[i])[0]:
                nd[j] -= 1
                if nd[j] == 0:
                    nxt.append(j)
        cur, fr = nxt, fr + 1
    return fronts, rank


def crowding(F, idx):
    dist = np.zeros(len(idx))
    for k in range(F.shape[1]):
        o = np.argsort(F[idx, k]); f = F[idx, k][o]
        span = max(f[-1] - f[0], 1e-12)
        dist[o[0]] = dist[o[-1]] = np.inf
        dist[o[1:-1]] += (f[2:] - f[:-2]) / span
    return dist


POP, NGEN, NV = 96, 120, 7
X = XL + RNG.random((POP, NV)) * (XU - XL)
F, cv, _ = evaluate(X)
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
        a, b = X[pick()], X[pick()]
        uu = RNG.random(NV)
        beta = np.where(uu <= .5, (2 * uu) ** (1 / 16), (1 / (2 * (1 - uu))) ** (1 / 16))
        for ch in (.5 * ((1 + beta) * a + (1 - beta) * b), .5 * ((1 - beta) * a + (1 + beta) * b)):
            m = RNG.random(NV) < 1 / NV * 1.6
            um = RNG.random(NV)
            dl = np.where(um < .5, (2 * um) ** (1 / 21) - 1, 1 - (2 * (1 - um)) ** (1 / 21))
            kids.append(np.clip(np.where(m, ch + dl * (XU - XL), ch), XL, XU))
    Xk = np.array(kids[:POP])
    Fk, ck, _ = evaluate(Xk)
    Xa, Fa, ca = np.vstack([X, Xk]), np.vstack([F, Fk]), np.concatenate([cv, ck])
    fronts, ra = nds(Fa, ca)
    cra = np.zeros(len(Fa))
    for fr in fronts:
        cra[fr] = crowding(Fa, fr)
    keep = np.array(sorted(range(len(Fa)), key=lambda i: (ra[i], -cra[i]))[:POP])
    X, F, cv = Xa[keep], Fa[keep], ca[keep]

fronts, _ = nds(F, cv)
pf = np.array(fronts[0]); pf = pf[cv[pf] == 0]
Fp, Xp = F[pf], X[pf]
print(f"cycle 9 front: {len(pf)} feasible | G {(-Fp[:,0]).min():.1f}..{(-Fp[:,0]).max():.1f} dB "
      f"| Gam {Fp[:,1].min():.3f}..{Fp[:,1].max():.3f}")
# cycle-5 KNEE reference in this evaluator: rebuild its profile as 3-seg approx of its power law
# exact rescoring: represent c5 KNEE (n=1.74, A=-0.14, x0=0.80) radii at u = 1/3, 2/3
u3 = np.array([1 / 3, 2 / 3])
rc5 = 21.69 + (143.58 - 21.69) * u3 ** 1.74 - 0.14 * (143.58 - 21.69) * \
    np.exp(-((u3 - 0.80) / 0.12) ** 2) * np.sin(np.pi * u3)
Xc5 = np.array([[19.60, 16.47, 21.69, rc5[0], rc5[1], 143.58, 501.26]])
Fc5, cvc5, _ = evaluate(Xc5)
print(f"c5 KNEE (3-seg approx): G={-Fc5[0,0]:.1f} Gam={Fc5[0,1]:.3f} Pc={Fc5[0,2]:.2f} cv={cvc5[0]:.1f}")
domd = ((Fp <= Fc5[0]).all(1) & (Fp < Fc5[0]).any(1)).sum()
print(f"front members dominating c5 KNEE: {domd}")
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
best = pf[np.argmax(np.minimum(gn, rn))]
print(f"c9 KNEE genes {np.round(X[best],1)}  G={-F[best,0]:.1f} Gam={F[best,1]:.3f} Pc={F[best,2]:.2f}")
gmax = pf[np.argmax(-Fp[:, 0])]
print(f"c9 MAX-GAIN {np.round(X[gmax],1)}  G={-F[gmax,0]:.1f} Gam={F[gmax,1]:.3f}")
np.savetxt("pareto_c9.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:]]), delimiter=",",
           header="re,Le,rt,ra,rb,rm,L,G_dB,refl,Pcost", comments="")

BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.7),
                        gridspec_kw=dict(left=.06, right=.98, top=.80, bottom=.13, wspace=.24))
ax = axs[0]
sc = ax.scatter(Fp[:, 1], -Fp[:, 0], c=Fp[:, 2], cmap="viridis", s=40, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=.02).set_label("print cost")
if cvc5[0] == 0:
    ax.plot(Fc5[0, 1], -Fc5[0, 0], "*", ms=16, color="#B45309", mec="#1F2937", mew=.8)
    ax.annotate("cycle-5 KNEE", (Fc5[0, 1], -Fc5[0, 0]), (Fc5[0, 1] + .02, -Fc5[0, 0] - 1),
                fontsize=8, arrowprops=dict(arrowstyle="->", color="#5c584d", lw=.8))
ax.plot(F[best, 1], -F[best, 0], "o", ms=10, mfc="none", mec="#B45309", mew=2.2)
ax.set_xlabel("reflection at lip cup"); ax.set_ylabel("G_speech (dB, corrected)")
ax.set_title("cycle-9 front (3-segment exponential bore, pop 96 x 120)", fontsize=9.3, loc="left")
ax.grid(alpha=.15)
ax = axs[1]
for i, nm, col in [(best, "c9 KNEE", "#B45309"), (gmax, "c9 MAX-GAIN", "#7a1f1f")]:
    rr = bore(X[i][None, :])[0]
    ax.plot(UH * X[i, 6], rr, color=col, lw=2, label=f"{nm} (L={X[i,6]:.0f})")
    ax.plot(UH * X[i, 6], -rr, color=col, lw=2)
rr5 = 21.69 + (143.58 - 21.69) * UH ** 1.74 - 0.14 * (143.58 - 21.69) * \
    np.exp(-((UH - .80) / .12) ** 2) * np.sin(np.pi * UH)
ax.plot(UH * 501, rr5, color="#8a8578", lw=1.5, ls="--", label="c5 KNEE")
ax.plot(UH * 501, -rr5, color="#8a8578", lw=1.5, ls="--")
ax.set_aspect("equal"); ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.set_xlabel("axial (mm)"); ax.set_ylabel("radius (mm)")
ax.set_title("what three flare constants buy", fontsize=9.3, loc="left")
fig.suptitle("THE FORGE, cycle 9  -  Clarion deepening: piecewise-exponential bore",
             fontsize=12, fontweight="bold", color="#1F2937", x=.06, ha="left")
fig.savefig("evolution_c9.png", dpi=150)
print("wrote pareto_c9.csv, evolution_c9.png")
