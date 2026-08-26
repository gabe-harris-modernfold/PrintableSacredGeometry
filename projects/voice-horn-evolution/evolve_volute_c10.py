"""Forge cycle 7 - the VOLUTE family enters the Forge.

1-D chain: lip cup (40 mm, 22.5 -> r_t equiv) -> exponential rectangular channel
S(x) = S_t e^{mx} folded on a log spiral -> oval bell (equivalent piston radiation,
tier-1.5 directivity cap). Glottal Norton source, corrected evaluator throughout.

Genes (5): r_t (equiv, mm), r_m (equiv, mm; oval bell), L (path, mm),
           k (spiral pitch), H (channel height, mm).

The Volute's signature constraint is GEOMETRIC - the packing gate:
  channel width w(x) = S(x)/H must fit the gap between successive spiral turns,
  gap(x) = r_c(x) (e^{2 pi k} - 1),  r_c(x) = r1 + k x / sqrt(1+k^2),  r1 = 35 mm,
  over the folded 80% of the path (the last 20% unrolls as the bell), with 4 mm walls,
  and the body must stay inside dia 310 x height 300.

Objectives (3): -G_speech (corrected), Gamma_cup, P_cost.
Run from this directory: python evolve_volute_c7.py -> pareto_c7.csv, evolution_c7.png
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
R_CUP, L_CUP, NC, NB = 22.5, 40.0, 12, 120
# cycle 10: freed genes - spiral start r1, fold fraction, height-graded channel H0->H1
#              r_t   r_m    L       k     H0     H1     r1    ffold
XL = np.array([13.0, 50.0, 300.0, 0.06, 40.0, 80.0, 25.0, 0.60])
XU = np.array([30.0, 123.0, 1000.0, 0.30, 120.0, 280.0, 60.0, 0.90])
RNG = np.random.default_rng(7)
UB = np.linspace(0, 1, NB + 1)


def evaluate(X):
    Fo, cvs, Gs_all = [], [], []
    for g in X:
        rt, rm, L, kk, H0, H1, R1g, ff = g
        Hx = H0 + (H1 - H0) * UB          # channel height along path
        St, Sm = np.pi * rt ** 2, np.pi * rm ** 2          # mm^2
        m = np.log(Sm / St) / L                            # 1/mm
        S_b = St * np.exp(m * UB * L)
        r_eq = np.sqrt(S_b / np.pi)
        r_cup = np.linspace(R_CUP, rt, NC + 1)
        r_all = np.concatenate([0.5 * (r_cup[:-1] + r_cup[1:]),
                                np.sqrt(0.5 * (S_b[:-1] + S_b[1:]) / np.pi)])
        S_all = np.pi * (r_all / 1000.0) ** 2
        d_all = np.concatenate([np.full(NC, L_CUP / NC), np.full(NB, L / NB)]) / 1000.0
        a_m = r_eq[-1] / 1000.0
        # radiation + corrected directivity
        x2 = 2 * K * a_m
        Rr = 1 - 2 * j1(x2) / x2
        Xr = 2 * struve(1, x2) / x2
        Zrad = RHO * C / (np.pi * a_m ** 2) * (Rr + 1j * Xr)
        slope = (r_eq[-1] - r_eq[-2]) / (L / NB)
        Qcone = 2.0 / max(1.0 - np.cos(np.arctan(slope)), 1e-6)
        Q = np.minimum((K * a_m) ** 2 / np.maximum(Rr, 1e-9), Qcone)
        p, Uv = Zrad.astype(complex), np.ones_like(Zrad, dtype=complex)
        for s in range(len(S_all) - 1, -1, -1):
            Zc = RHO * C / S_all[s]
            kd = K * d_all[s]
            p, Uv = np.cos(kd) * p + 1j * Zc * np.sin(kd) * Uv, \
                    1j * np.sin(kd) / Zc * p + np.cos(kd) * Uv
        Zin = p / Uv
        src = ZS / (ZS + Zin)
        W = 0.5 * np.real(Zrad) * np.abs(src / Uv) ** 2
        Wref = 0.5 * RHO * C * K ** 2 / (4 * np.pi)
        G = 10 * np.log10(np.maximum(W / Wref, 1e-12)) + 10 * np.log10(Q / 2.0)
        Gsp = (G * WF).sum()
        A_cup = np.pi * (R_CUP / 1000.0) ** 2
        Gam = np.abs((Zin - RHO * C / A_cup) / (Zin + RHO * C / A_cup)).mean()
        # packing gate over the folded fraction, with local height
        xs = UB * L
        fold = xs <= ff * L
        w = S_b / Hx                                       # channel width, mm
        r_c = R1g + kk * xs / np.sqrt(1 + kk ** 2)
        gap = r_c * (np.exp(2 * np.pi * kk) - 1)
        cv = np.maximum(0, (w + 4.0 - gap)[fold]).sum() / 10.0
        cv += max(0, r_c[fold][-1] + w[fold][-1] / 2 + 2 - 155.0) / 5.0
        ripple = G[2:].max() - G[2:].min()
        cv += max(0, ripple - 12.0)
        # print cost: channel shell + lids + bell segment
        per = 2 * (w + Hx)                                 # mm perimeter
        vol_L = (per.mean() * L * 5.0) / 1e6 + 0.7
        Fo.append([-Gsp, Gam, vol_L + 1.0])
        cvs.append(cv)
        Gs_all.append(G)
    return np.array(Fo), np.array(cvs), np.array(Gs_all)


def nds(F, cv):
    N = len(F)
    dom = np.zeros((N, N), bool)
    for i in range(N):
        ci = cv[i]
        bf = (ci == 0) & (cv == 0)
        df = (F[i][None, :] <= F).all(1) & (F[i][None, :] < F).any(1)
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
    n, mm = len(idx), F.shape[1]
    dist = np.zeros(n)
    for k2 in range(mm):
        o = np.argsort(F[idx, k2]); f = F[idx, k2][o]
        span = max(f[-1] - f[0], 1e-12)
        dist[o[0]] = dist[o[-1]] = np.inf
        dist[o[1:-1]] += (f[2:] - f[:-2]) / span
    return dist


POP, NGEN, NV = 64, 100, 8
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
        beta = np.where(uu <= 0.5, (2 * uu) ** (1 / 16), (1 / (2 * (1 - uu))) ** (1 / 16))
        for ch in (0.5 * ((1 + beta) * a + (1 - beta) * b),
                   0.5 * ((1 - beta) * a + (1 + beta) * b)):
            msk = RNG.random(NV) < 0.2
            um = RNG.random(NV)
            dl = np.where(um < 0.5, (2 * um) ** (1 / 21) - 1, 1 - (2 * (1 - um)) ** (1 / 21))
            kids.append(np.clip(np.where(msk, ch + dl * (XU - XL), ch), XL, XU))
    Xk = np.array(kids[:POP])
    Fk, cvk, _ = evaluate(Xk)
    Xa, Fa, ca = np.vstack([X, Xk]), np.vstack([F, Fk]), np.concatenate([cv, cvk])
    fronts, ra = nds(Fa, ca)
    cra = np.zeros(len(Fa))
    for fr in fronts:
        cra[fr] = crowding(Fa, fr)
    keep = np.array(sorted(range(len(Fa)), key=lambda i: (ra[i], -cra[i]))[:POP])
    X, F, cv = Xa[keep], Fa[keep], ca[keep]

fronts, rank = nds(F, cv)
pf = np.array(fronts[0]); pf = pf[cv[pf] == 0]
Fp, Xp = F[pf], X[pf]
print(f"cycle 10 (VOLUTE v2) front: {len(pf)} feasible of {POP}")
print(f"G: {(-Fp[:,0]).min():.1f}..{(-Fp[:,0]).max():.1f} dB | Gam: "
      f"{Fp[:,1].min():.3f}..{Fp[:,1].max():.3f}")
fc = C * np.log((Xp[:, 1] / Xp[:, 0]) ** 2) / (Xp[:, 2] / 1000) / (4 * np.pi)
print(f"cutoffs f_c on the front: {fc.min():.0f}..{fc.max():.0f} Hz | path L: "
      f"{Xp[:,2].min():.0f}..{Xp[:,2].max():.0f} mm | pitch k: {Xp[:,3].min():.2f}..{Xp[:,3].max():.2f}")
# hand Volute with its true height-graded intent: H 55 -> 295, r1 35, fold 0.8
Fh, cvh, Gh = evaluate(np.array([[26.5, 122.5, 492.0, 0.154, 55.0, 295.0, 35.0, 0.80]]))
# cycle-7 KNEE rescored in this evaluator (fixed H expressed as flat grade)
Fk7, cvk7, _ = evaluate(np.array([[22.14, 123.0, 506.05, 0.14, 200.0, 200.0, 35.0, 0.80]]))
print(f"c7 KNEE rescored: G={-Fk7[0,0]:.1f} Gam={Fk7[0,1]:.3f} cv={cvk7[0]:.2f}")
print(f"hand VOLUTE: G={-Fh[0,0]:.1f} dB Gam={Fh[0,1]:.3f} cv={cvh[0]:.2f} "
      f"f_c={C*np.log((122.5/26.5)**2)/0.492/(4*np.pi):.0f} Hz")
domd = ((Fp <= Fh[0]).all(1) & (Fp < Fh[0]).any(1)).sum() if cvh[0] == 0 else -1
print(f"front members dominating hand Volute: {domd}")
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
best = pf[np.argmax(np.minimum(gn, rn))]
print(f"KNEE: genes {np.round(X[best],2)}  G={-F[best,0]:.1f} Gam={F[best,1]:.3f} "
      f"f_c={C*np.log((X[best,1]/X[best,0])**2)/(X[best,2]/1000)/(4*np.pi):.0f} Hz")
np.savetxt("pareto_c10.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:], fc[:, None]]), delimiter=",",
           header="rt,rm,L,k,H0,H1,r1,ffold,G_dB,refl,Pcost,fc_Hz", comments="")

BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.8),
                        gridspec_kw=dict(left=0.06, right=0.98, top=0.80, bottom=0.12,
                                         wspace=0.26))
ax = axs[0]
sc = ax.scatter(Fp[:, 1], -Fp[:, 0], c=fc, cmap="viridis", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("cutoff f_c (Hz)")
if cvh[0] == 0:
    ax.plot(Fh[0, 1], -Fh[0, 0], "*", ms=17, color="#1D4ED8", mec="#1F2937", mew=0.8)
    ax.annotate("hand VOLUTE (Rev 1)", (Fh[0, 1], -Fh[0, 0]),
                (Fh[0, 1] + 0.03, -Fh[0, 0] - 1.2), fontsize=8,
                arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.8))
ax.plot(F[best, 1], -F[best, 0], "o", ms=10, mfc="none", mec="#1D4ED8", mew=2.2)
ax.set_xlabel("mean reflection at the lip cup")
ax.set_ylabel("G_speech (dB, corrected)")
ax.set_title("cycle-7 front, colored by cutoff - the fold trades\nlength for packing, "
             "and packing is the binding gate", fontsize=9.3, loc="left")
ax.grid(alpha=0.15)
ax = axs[1]
ax.scatter(Xp[:, 2], Xp[:, 3], c=-Fp[:, 0], cmap="magma_r", s=42, edgecolors="none")
kline = np.linspace(0.08, 0.30, 50)
ax.plot((155 - 35) * np.sqrt(1 + kline ** 2) / kline / 0.8, kline, color="#7a1f1f", lw=1.2,
        ls="--")
ax.text(660, 0.27, "spiral capacity limit\n(fold must fit dia 310)", fontsize=7.4,
        color="#7a1f1f")
ax.plot(492, 0.154, "*", ms=15, color="#1D4ED8", mec="#1F2937", mew=0.8)
ax.set_xlabel("path length L (mm)"); ax.set_ylabel("spiral pitch k")
cb = fig.colorbar(ax.collections[0], ax=ax, pad=0.02); cb.set_label("G_speech (dB)")
ax.set_title("where the front lives in (L, k): the packing gate\nshapes the family",
             fontsize=9.3, loc="left")
ax.grid(alpha=0.15)
fig.suptitle("THE FORGE, cycle 10  -  VOLUTE v2: freed fold, start radius, height-graded channel",
             fontsize=12, fontweight="bold", color="#1F2937", x=0.06, ha="left")
fig.savefig("evolution_c10.png", dpi=150)
print("wrote pareto_c10.csv, evolution_c10.png")
