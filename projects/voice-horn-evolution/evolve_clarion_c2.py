"""Forge cycle 2 for the Clarion family. Changes vs cycle 1 (evolve_clarion.py):

1. Throat floor raised: r_t >= 13 mm (anatomical epilaryngeal coupling, S15.11's 1:6 step).
2. Finite glottal source: Norton flow source U_g with internal impedance Z_s = 4 MPa*s/m^3
   [PARAM given, coarse - mid-range of 30-100 cgs-ohm glottal resistance in phonation].
   Delivered throat flow U_t = U_g * Z_s / (Z_s + Z_in); the tiny-throat free lunch closes.
3. Objectives reduced to the three load-bearing ones (many-objective fix recorded in Rev 6):
   -G_speech, Gamma_bar, P_cost = vol_L + (print segments - 1).

Also re-scores the cycle-1 picks and the hand-derived Clarion under the glottal source.
Run from this directory: python evolve_clarion_c2.py
Outputs: pareto_c2.csv, evolution_c2.png
"""
import numpy as np
from scipy.special import j1, struve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6                       # glottal source impedance, Pa*s/m^3
FREQS = np.array([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000], float)
WF = np.array([0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2])
WF = WF / WF.sum()
K = 2 * np.pi * FREQS / C
NSEG = 160
UGRID = np.linspace(0, 1, NSEG + 1)
XL = np.array([13.0, 60.0, 120.0, 0.35, -0.30, 0.20])   # throat floor raised 9 -> 13
XU = np.array([30.0, 155.0, 580.0, 4.00, 0.30, 0.80])
RNG = np.random.default_rng(7)


def bore_from_genes(X):
    rt, rm, L, n, A, x0 = [X[:, i:i + 1] for i in range(6)]
    u = UGRID[None, :]
    r = rt + (rm - rt) * u ** n
    r = r + A * (rm - rt) * np.exp(-((u - x0) / 0.12) ** 2) * np.sin(np.pi * u)
    return r, X[:, 2]


def evaluate_profiles(r, L):
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
    src = ZS / (ZS + Zin)                        # Norton delivery factor (U_g = 1)
    W = 0.5 * np.real(Zrad) * np.abs(src / Uv) ** 2
    Wref = 0.5 * RHO * C * K ** 2 / (4 * np.pi)  # bare voice: monopole, Zrad_ref << ZS
    Q = (K[None, :] * a_m) ** 2 / np.maximum(R1, 1e-9)
    G = 10 * np.log10(np.maximum(W / Wref[None, :], 1e-12)) \
        + 10 * np.log10(np.maximum(Q / 2.0, 1e-9))
    Gs = (G * WF[None, :]).sum(1)
    Gam = np.abs((Zin - RHO * C / St) / (Zin + RHO * C / St)).mean(1)
    dr = np.diff(r, axis=1)
    dx = (L[:, None] / NSEG) * np.ones_like(dr)
    Alat = (2 * np.pi * 0.5 * (r[:, :-1] + r[:, 1:]) * np.sqrt(dx ** 2 + dr ** 2)).sum(1)
    vol_L = Alat * 5.0 / 1e6
    segs = np.ceil((L + 40.0) / 290.0)
    Pcost = vol_L + (segs - 1.0)
    ripple = G[:, 2:].max(1) - G[:, 2:].min(1)
    choke = np.maximum(0.6 * r[:, 0], 8.0) - r.min(1)
    cv = np.maximum(0, ripple - 12.0) + np.maximum(0, choke)
    F = np.stack([-Gs, Gam, Pcost], 1)
    return F, cv, G


def evaluate(X):
    r, L = bore_from_genes(X)
    return evaluate_profiles(r, L)


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


def polymut(x, eta=20, pm=1 / 6):
    y = x.copy()
    mask = RNG.random(x.shape) < pm
    u = RNG.random(x.shape)
    delta = np.where(u < 0.5, (2 * u) ** (1 / (eta + 1)) - 1, 1 - (2 * (1 - u)) ** (1 / (eta + 1)))
    return np.clip(np.where(mask, y + delta * (XU - XL), y), XL, XU)


POP, NGEN = 64, 75
X = XL + RNG.random((POP, 6)) * (XU - XL)
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
print(f"cycle 2 front: {len(pf)} feasible non-dominated of {POP} "
      f"(3 objectives - selection pressure restored)" if len(pf) < POP else
      f"cycle 2 front: {len(pf)} of {POP} non-dominated")
print(f"G range: {(-Fp[:,0]).min():.1f} .. {(-Fp[:,0]).max():.1f} dB | "
      f"Gam: {Fp[:,1].min():.3f} .. {Fp[:,1].max():.3f} | Pcost: {Fp[:,2].min():.2f} .. {Fp[:,2].max():.2f}")
print(f"throat radius across front: {Xp[:,0].min():.1f} .. {Xp[:,0].max():.1f} mm "
      f"(floor = 13.0); at floor (<13.5): {(Xp[:,0] < 13.5).sum()}/{len(pf)}")

# re-score hand Clarion and cycle-1 picks under the glottal source
r_cl = 26.5 * (125.0 / 26.5) ** UGRID
Fc, cvc, Gc = evaluate_profiles(r_cl[None, :], np.array([280.0]))
C1 = np.array([[9.00, 155.00, 534.64, 2.53, -0.01, 0.57],    # c1 MAX-GAIN
               [9.21, 145.48, 580.00, 2.40, -0.01, 0.63],    # c1 KNEE
               [14.93, 60.00, 120.00, 1.68, -0.12, 0.20]])   # c1 MIN-VOL
Fc1, _, _ = evaluate(C1)
print(f"hand Clarion rescored: G={-Fc[0,0]:.1f} dB (was 16.5 under ideal source) "
      f"Gam={Fc[0,1]:.3f}")
print(f"c1 MAX-GAIN rescored: G={-Fc1[0,0]:.1f} dB (was 28.2 - the tiny-throat mirage)")
print(f"c1 KNEE rescored:     G={-Fc1[1,0]:.1f} dB (was 27.5)")
domd = ((Fp <= Fc[0]).all(1) & (Fp < Fc[0]).any(1)).sum()
print(f"front members dominating hand Clarion: {domd}")

gmax = pf[np.argmax(-Fp[:, 0])]
cmin = pf[np.argmin(Fp[:, 2])]
gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
knee = pf[np.argmax(np.minimum(gn, rn))]
picks = [("MAX-GAIN", gmax, "#7a1f1f"), ("KNEE", knee, "#B45309"), ("MIN-COST", cmin, "#1D4ED8")]
for nm, i, _ in picks:
    print(nm, "genes:", np.round(X[i], 2),
          " G=%.1f Gam=%.3f Pcost=%.2f" % (-F[i, 0], F[i, 1], F[i, 2]))

np.savetxt("pareto_c2.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:]]), delimiter=",",
           header="rt_mm,rm_mm,L_mm,n,bulgeA,bulgeX,G_speech_dB,refl,print_cost", comments="")

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
                      top=0.87, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(Fp[:, 1], -Fp[:, 0], c=Fp[:, 2], cmap="viridis", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("print cost (L + extra segments)")
ax.plot(Fc[0, 1], -Fc[0, 0], "*", ms=18, color="#B45309", mec=INK, mew=0.8)
ax.annotate("hand CLARION\n(rescored, glottal source)", (Fc[0, 1], -Fc[0, 0]),
            (Fc[0, 1] + 0.06, -Fc[0, 0] - 3.2), fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.8))
for j, nm in [(0, "c1 MAX-GAIN"), (1, "c1 KNEE")]:
    ax.plot(Fc1[j, 1], -Fc1[j, 0], "X", ms=10, color="#7a1f1f", alpha=0.6)
    ax.annotate(f"{nm}\nrescored", (Fc1[j, 1], -Fc1[j, 0]),
                (Fc1[j, 1] + 0.05, -Fc1[j, 0] + (1.2 if j else -2.6)), fontsize=7,
                color="#7a1f1f", arrowprops=dict(arrowstyle="->", color="#7a1f1f",
                                                 lw=0.7, alpha=0.6))
for nm, i, col in picks:
    ax.plot(F[i, 1], -F[i, 0], "o", ms=9, mfc="none", mec=col, mew=2)
    ax.annotate(nm, (F[i, 1], -F[i, 0]), (F[i, 1] + 0.015, -F[i, 0] + 0.7),
                fontsize=7.5, color=col)
ax.set_xlabel("mean input reflection (phonation-load proxy)")
ax.set_ylabel("speech-weighted forward gain (dB)")
ax.set_title("cycle-2 front (glottal source, throat floor 13 mm)\n"
             "X = cycle-1 winners re-scored under the honest source", fontsize=9, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[0, 1])
for nm, i, col in picks:
    rr, LL = bore_from_genes(X[i][None, :])
    xx = np.linspace(0, LL[0], NSEG + 1)
    ax.plot(xx, rr[0], color=col, lw=2,
            label=f"{nm} (rt={rr[0,0]:.0f}, rm={rr[0,-1]:.0f}, L={LL[0]:.0f})")
    ax.plot(xx, -rr[0], color=col, lw=2)
ax.plot(np.linspace(0, 280, NSEG + 1), r_cl, color="#B45309", lw=1.6, ls="--",
        label="hand CLARION (dashed)")
ax.plot(np.linspace(0, 280, NSEG + 1), -r_cl, color="#B45309", lw=1.6, ls="--")
ax.set_aspect("equal")
ax.set_xlabel("axial position (mm)"); ax.set_ylabel("radius (mm)")
ax.legend(fontsize=6.8, loc="upper left", frameon=False)
ax.set_title("cycle-2 bores: phonatable throats by constraint", fontsize=9.5, loc="left")

ax = fig.add_subplot(gs[1, 0])
for nm, i, col in picks:
    rr, LL = bore_from_genes(X[i][None, :])
    _, _, Gi = evaluate_profiles(rr, LL)
    ax.semilogx(FREQS, Gi[0], color=col, lw=1.9, marker="o", ms=3.5, label=nm)
ax.semilogx(FREQS, Gc[0], color="#B45309", lw=1.6, ls="--", marker="s", ms=3.5,
            label="hand CLARION")
ax.axvspan(1000, 4000, color="#B45309", alpha=0.07)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("forward gain vs bare voice (dB)")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("third-octave gain, glottal source (tier-1)", fontsize=9.5, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
try:
    c1front = np.loadtxt("pareto.csv", delimiter=",", skiprows=1)
    ax.hist(c1front[:, 0], bins=np.arange(8, 32, 1), color="#7a1f1f", alpha=0.45,
            label="cycle-1 front (floor 9 mm)")
except OSError:
    pass
ax.hist(Xp[:, 0], bins=np.arange(8, 32, 1), color="#1D4ED8", alpha=0.55,
        label="cycle-2 front (floor 13 mm)")
ax.axvline(9, color="#7a1f1f", ls=":", lw=1); ax.axvline(13, color="#1D4ED8", ls=":", lw=1)
ax.axvline(26.5, color="#B45309", ls="--", lw=1.4)
ax.text(26.8, ax.get_ylim()[1] * 0.85, "hand\nCLARION", fontsize=7.5, color="#B45309")
ax.set_xlabel("throat radius r_t (mm)"); ax.set_ylabel("front members")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("where evolution puts the throat, by cycle\n"
             "(does the glottal source lift the population off the floor?)",
             fontsize=9, loc="left")

fig.suptitle("THE FORGE, cycle 2  -  glottal source impedance + anatomical throat floor\n"
             f"Z_s = 4 MPa*s/m^3; 3 objectives; {POP*NGEN+POP} tier-1 evaluations, seed 7",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("evolution_c2.png", dpi=150)
print("wrote pareto_c2.csv, evolution_c2.png")
