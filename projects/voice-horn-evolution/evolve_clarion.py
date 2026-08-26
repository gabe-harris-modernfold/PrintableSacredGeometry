"""Tier-1 evolutionary search for the Clarion family (axisymmetric voice horn).

Implements the simulation-screened multi-objective protocol from
docs/voice-horn-pocs.md ("The Forge"): a 1-D Webster transfer-matrix evaluator
(staircase of cylindrical duct segments, baffled-piston radiation load) driven
by a self-contained NSGA-II (no external optimizer dependency, per
docs/PYTHON_TOOLING.md policy).

Genes (6): throat radius rt (mm), mouth radius rm (mm), length L (mm),
flare exponent n, mid-bore bulge amplitude A, bulge position x0.
Bore: r(u) = rt + (rm-rt)*u^n + A*(rm-rt)*exp(-((u-x0)/0.12)^2)*sin(pi*u)

Objectives (5, all minimized):
  -G_speech   speech-weighted forward gain vs bare (baffled monopole) source
  -DI_14      mean directivity index, 1-4 kHz third-octave bands
   Gamma_bar  mean input reflection |(Zin-Zc_throat)/(Zin+Zc_throat)|  (phonation load proxy)
   vol_L      material volume (lateral area x 5 mm F1 sandwich wall), litres
   P_print    print segments + length penalty

Hard constraints (violations dominate, per Law 15.13 / S15.3 discipline):
   spectral ripple (500 Hz - 4 kHz) <= 12 dB;  no bore choke below max(0.6*rt, 8 mm).

Run from this directory:  python evolve_clarion.py
Outputs: pareto.csv, evolution.png (bare relative filenames).
"""
import numpy as np
from scipy.special import j1, struve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
FREQS = np.array([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000], float)
WF = np.array([0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2])
WF = WF / WF.sum()
K = 2 * np.pi * FREQS / C
NSEG = 160
UGRID = np.linspace(0, 1, NSEG + 1)
XL = np.array([9.0, 60.0, 120.0, 0.35, -0.30, 0.20])
XU = np.array([30.0, 155.0, 580.0, 4.00, 0.30, 0.80])
RNG = np.random.default_rng(7)


def bore_from_genes(X):
    rt, rm, L, n, A, x0 = [X[:, i:i + 1] for i in range(6)]
    u = UGRID[None, :]
    r = rt + (rm - rt) * u ** n
    r = r + A * (rm - rt) * np.exp(-((u - x0) / 0.12) ** 2) * np.sin(np.pi * u)
    return r, X[:, 2]


def evaluate_profiles(r, L):
    """r: (N, NSEG+1) radii mm; L: (N,) lengths mm. Returns objectives, cv, G(f)."""
    r = np.maximum(r, 1.0)
    rt = r[:, 0:1]
    a_m = r[:, -1:] / 1000.0
    Sseg = np.pi * ((0.5 * (r[:, :-1] + r[:, 1:])) / 1000.0) ** 2
    St = np.pi * (rt / 1000.0) ** 2
    Sm = np.pi * a_m ** 2
    d = (L[:, None] / 1000.0) / NSEG
    x = 2 * K[None, :] * a_m
    R1 = 1 - 2 * j1(x) / x
    X1 = 2 * struve(1, x) / x
    Zrad = RHO * C / Sm * (R1 + 1j * X1)
    p = Zrad.astype(complex)
    Uv = np.ones_like(p)
    kd = K[None, :] * d
    cs, sn = np.cos(kd), np.sin(kd)
    for s in range(NSEG - 1, -1, -1):
        Zc = RHO * C / Sseg[:, s:s + 1]
        p, Uv = cs * p + 1j * Zc * sn * Uv, 1j * sn / Zc * p + cs * Uv
    Zin = p / Uv
    W = 0.5 * np.real(Zrad) * np.abs(1.0 / Uv) ** 2          # unit throat volume velocity
    Wref = 0.5 * RHO * C * K ** 2 / (4 * np.pi)              # baffled monopole, same source
    Q = (K[None, :] * a_m) ** 2 / np.maximum(R1, 1e-9)       # piston directivity factor
    G = 10 * np.log10(np.maximum(W / Wref[None, :], 1e-12)) \
        + 10 * np.log10(np.maximum(Q / 2.0, 1e-9))
    Gs = (G * WF[None, :]).sum(1)
    DI = (10 * np.log10(np.maximum(Q[:, 5:] / 2.0, 1e-9))).mean(1)
    Gam = np.abs((Zin - RHO * C / St) / (Zin + RHO * C / St)).mean(1)
    dr = np.diff(r, axis=1)
    dx = (L[:, None] / NSEG) * np.ones_like(dr)
    Alat = (2 * np.pi * 0.5 * (r[:, :-1] + r[:, 1:]) * np.sqrt(dx ** 2 + dr ** 2)).sum(1)
    vol_L = Alat * 5.0 / 1e6
    segs = np.ceil((L + 40.0) / 290.0)
    Pp = segs + L / 280.0 * 0.2
    ripple = G[:, 2:].max(1) - G[:, 2:].min(1)
    choke = np.maximum(0.6 * r[:, 0], 8.0) - r.min(1)
    cv = np.maximum(0, ripple - 12.0) + np.maximum(0, choke)
    F = np.stack([-Gs, -DI, Gam, vol_L, Pp], 1)
    return F, cv, G


def evaluate(X):
    r, L = bore_from_genes(X)
    return evaluate_profiles(r, L)


# ------------------------- NSGA-II (self-contained) -------------------------
def nds(F, cv):
    N = len(F)
    dominates = np.zeros((N, N), bool)
    for i in range(N):
        fi, ci = F[i], cv[i]
        better_c = (ci < cv) & (cv > 0)
        both_feas = (ci == 0) & (cv == 0)
        dom_f = (F[i][None, :] <= F).all(1) & (F[i][None, :] < F).any(1)
        dominates[i] = better_c | (both_feas & dom_f) | ((ci == 0) & (cv > 0))
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
    c1 = 0.5 * ((1 + beta) * a + (1 - beta) * b)
    c2 = 0.5 * ((1 - beta) * a + (1 + beta) * b)
    return np.clip(c1, XL, XU), np.clip(c2, XL, XU)


def polymut(x, eta=20, pm=1 / 6):
    y = x.copy()
    mask = RNG.random(x.shape) < pm
    u = RNG.random(x.shape)
    delta = np.where(u < 0.5, (2 * u) ** (1 / (eta + 1)) - 1, 1 - (2 * (1 - u)) ** (1 / (eta + 1)))
    y = np.where(mask, y + delta * (XU - XL), y)
    return np.clip(y, XL, XU)


POP, NGEN = 64, 75
X = XL + RNG.random((POP, 6)) * (XU - XL)
F, cv, _ = evaluate(X)
trace = []
for gen in range(NGEN):
    fronts, rank = nds(F, cv)
    crowd = np.zeros(len(F))
    for fr in fronts:
        crowd[fr] = crowding(F, fr)
    # tournament selection
    def pick():
        i, j = RNG.integers(0, len(F), 2)
        if rank[i] != rank[j]:
            return i if rank[i] < rank[j] else j
        return i if crowd[i] > crowd[j] else j
    kids = []
    while len(kids) < POP:
        a, b = X[pick()], X[pick()]
        c1, c2 = sbx(a, b)
        kids.append(polymut(c1))
        kids.append(polymut(c2))
    Xk = np.array(kids[:POP])
    Fk, cvk, _ = evaluate(Xk)
    Xall = np.vstack([X, Xk]); Fall = np.vstack([F, Fk]); cvall = np.concatenate([cv, cvk])
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
pf = np.array(fronts[0])
pf = pf[cv[pf] == 0]
Fp, Xp = F[pf], X[pf]
print(f"final front: {len(pf)} feasible non-dominated of {POP}")
print(f"G_speech range on front: {(-Fp[:,0]).min():.1f} .. {(-Fp[:,0]).max():.1f} dB")
print(f"reflection range: {Fp[:,2].min():.3f} .. {Fp[:,2].max():.3f}")
print(f"volume range: {Fp[:,3].min():.2f} .. {Fp[:,3].max():.2f} L")

# hand-derived Clarion baseline: exponential bore, rt=26.5, rm=125, L=280
r_cl = 26.5 * (125.0 / 26.5) ** UGRID
Fc, cvc, Gc = evaluate_profiles(r_cl[None, :], np.array([280.0]))
print(f"hand Clarion: G={-Fc[0,0]:.1f} dB  DI={-Fc[0,1]:.1f} dB  Gam={Fc[0,2]:.3f} "
      f"vol={Fc[0,3]:.2f} L  cv={cvc[0]:.2f}")
# is it dominated by any front member?
domd = ((Fp <= Fc[0]).all(1) & (Fp < Fc[0]).any(1)).sum()
print(f"front members dominating the hand Clarion: {domd}")

# pick three named candidates: max gain, knee (min dist to utopia in gain/refl), min volume
gmax = pf[np.argmax(-Fp[:, 0])]
vmin = pf[np.argmin(Fp[:, 3])]
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 2].max() - Fp[:, 2]) / max(np.ptp(Fp[:, 2]), 1e-9)
knee = pf[np.argmax(np.minimum(gn, rn))]
picks = [("MAX-GAIN", gmax, "#7a1f1f"), ("KNEE", knee, "#B45309"), ("MIN-VOL", vmin, "#1D4ED8")]
for nm, i, _ in picks:
    print(nm, "genes:", np.round(X[i], 2), " G=%.1f Gam=%.3f vol=%.2fL" %
          (-F[i, 0], F[i, 2], F[i, 3]))

np.savetxt("pareto.csv",
           np.hstack([Xp, -Fp[:, :2], Fp[:, 2:]]), delimiter=",",
           header="rt_mm,rm_mm,L_mm,n,bulgeA,bulgeX,G_speech_dB,DI14_dB,refl,vol_L,print_pen",
           comments="")

# ------------------------------- figure -------------------------------
BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
INK = "#1F2937"
fig = plt.figure(figsize=(11.8, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.28, left=0.07, right=0.97,
                      top=0.88, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(Fp[:, 2], -Fp[:, 0], c=Fp[:, 3], cmap="viridis", s=38, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("material volume (L)")
ax.plot(Fc[0, 2], -Fc[0, 0], "*", ms=18, color="#B45309", mec=INK, mew=0.8)
ax.annotate("hand-derived CLARION\n(Rev 1, no optimizer)", (Fc[0, 2], -Fc[0, 0]),
            (0.45, -Fc[0, 0] - 4.5), fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.8))
for nm, i, col in picks:
    ax.plot(F[i, 2], -F[i, 0], "o", ms=9, mfc="none", mec=col, mew=2)
    ax.annotate(nm, (F[i, 2], -F[i, 0]), (F[i, 2] + 0.02, -F[i, 0] + 0.8),
                fontsize=7.5, color=col)
ax.set_xlabel("mean input reflection  (phonation-load proxy, S15.2)")
ax.set_ylabel("speech-weighted forward gain (dB)")
ax.set_title("the Pareto front: gain vs back-pressure, colored by material",
             fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[0, 1])
for nm, i, col in picks:
    rr, LL = bore_from_genes(X[i][None, :])
    xx = np.linspace(0, LL[0], NSEG + 1)
    ax.plot(xx, rr[0], color=col, lw=2, label=f"{nm}  (L={LL[0]:.0f}, rm={rr[0,-1]:.0f})")
    ax.plot(xx, -rr[0], color=col, lw=2)
ax.plot(np.linspace(0, 280, NSEG + 1), r_cl, color="#B45309", lw=1.6, ls="--",
        label="hand CLARION (exp, dashed)")
ax.plot(np.linspace(0, 280, NSEG + 1), -r_cl, color="#B45309", lw=1.6, ls="--")
ax.set_aspect("equal")
ax.set_xlabel("axial position (mm)"); ax.set_ylabel("radius (mm)")
ax.legend(fontsize=7, loc="upper left", frameon=False)
ax.set_title("evolved bores vs the hand-derived exponential", fontsize=9.5, loc="left")

ax = fig.add_subplot(gs[1, 0])
for nm, i, col in picks:
    rr, LL = bore_from_genes(X[i][None, :])
    _, _, Gi = evaluate_profiles(rr, LL)
    ax.semilogx(FREQS, Gi[0], color=col, lw=1.9, marker="o", ms=3.5, label=nm)
ax.semilogx(FREQS, Gc[0], color="#B45309", lw=1.6, ls="--", marker="s", ms=3.5,
            label="hand CLARION")
ax.axvspan(1000, 4000, color="#B45309", alpha=0.07)
ax.text(1900, ax.get_ylim()[0] + 1, "intelligibility band", fontsize=7.5, color="#B45309")
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("forward gain vs bare voice (dB)")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("third-octave gain curves (tier-1 model)", fontsize=9.5, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
ax.plot(np.arange(1, NGEN + 1), trace, color=INK, lw=2)
ax.set_xlabel("generation"); ax.set_ylabel("best feasible G_speech (dB)")
ax.set_title("convergence: NSGA-II, pop 64, seed 7\n"
             f"{POP*NGEN + POP} tier-1 evaluations (Webster TM, 160 segments, 12 bands)",
             fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

fig.suptitle("THE FORGE, first cycle  -  evolved Clarion family vs the hand-derived design\n"
             "tier-1 screen only: coiling, Fleece and skin are frozen; the bore alone evolves",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("evolution.png", dpi=150)
print("wrote pareto.csv, evolution.png")
