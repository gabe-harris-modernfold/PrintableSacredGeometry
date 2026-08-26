"""Forge cycle 8 - the HALO enters the Forge: in-situ chord tuning.

The Halo is not a horn; it is a lumped acoustic network. Chamber (fixed dia 160 x 200,
V_ch = 4.02 L) with a radiating annular slot, four Helmholtz spheres shunted off the chamber
node, driven by the glottal Norton source through the inlet. The four resonators share the
chamber and PULL EACH OTHER - closed-form isolated tuning (Law 14.6) lands off-pitch in situ.
Evolution tunes the coupled system.

Genes (7): V1..V4 (cm^3), neck radius a_n (mm, shared), neck length L_n (mm, shared),
           slot width g_s (mm).
Targets: just chord 220 / 330 / 440 / 550 Hz, shared t60 = 0.30 s (Law 15.13a).

Network (e^{j omega t}): chamber node pressure P; branches to ground:
  chamber compliance   C_ch = V_ch / (rho c^2)
  slot                 M_s = rho L_s,eff / A_s   + radiation R_s = rho c k^2 / (2 pi) (baffled)
  resonator i          M_i = rho L_eff / A_n ,  C_i = V_i / (rho c^2) ,
                       R_i = 2 L_eff sqrt(2 mu rho omega) / (A_n a_n)   (neck viscous)
Metrics from the radiated-power spectrum: peak nearest each target -> f_i (cents error),
Q_i -> t60_i = 2.2 Q_i / f_i.

Objectives (3): tuning error (cents, sum), ring error sum(|t60_i - 0.3|/0.3),
                -level (mean radiated dB at the four targets).
Run from this directory: python evolve_halo_c8.py -> pareto_c8.csv, evolution_c8.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C, MU = 1.2, 343.0, 1.81e-5
ZS = 4.0e6
V_CH = np.pi * 0.08 ** 2 * 0.200                    # m^3 (dia 160 x 200)
C_CH = V_CH / (RHO * C ** 2)
FT = np.array([220.0, 330.0, 440.0, 550.0])
FRQ = np.linspace(150, 700, 1101)
OM = 2 * np.pi * FRQ
# cycle 11: per-sphere neck radii (constant-t60 wants bigger necks on higher spheres)
#              V1    V2    V3    V4   a1   a2   a3   a4   L_n   g_s
XL = np.array([150., 60., 40., 25., 4.0, 4.0, 4.0, 4.0, 15.0, 4.0])
XU = np.array([700., 400., 250., 160., 16.0, 16.0, 16.0, 16.0, 50.0, 16.0])
RNG = np.random.default_rng(7)


def network(g):
    V = g[:4] * 1e-6
    a_i = g[4:8] * 1e-3
    L_n, g_s = g[8] * 1e-3, g[9] * 1e-3
    A_s = np.pi * 0.160 * g_s                       # annular slot area
    L_s = 0.005 + 1.7 * np.sqrt(A_s / np.pi) * 0.5  # 5 mm wall + end correction
    Y = 1j * OM * C_CH                               # chamber compliance admittance
    Zs_branch = 1j * OM * RHO * L_s / A_s + RHO * C * (OM / C) ** 2 / (2 * np.pi)
    Y = Y + 1 / Zs_branch
    Zres = []
    for Vi, an in zip(V, a_i):
        A_n = np.pi * an ** 2
        L_eff = L_n + 1.7 * an
        Ci = Vi / (RHO * C ** 2)
        Mi = RHO * L_eff / A_n
        Ri = 2 * L_eff * np.sqrt(2 * MU * RHO * OM) / (A_n * an)
        Zi = Ri + 1j * OM * Mi + 1 / (1j * OM * Ci)
        Zres.append(Zi)
        Y = Y + 1 / Zi
    P = 1.0 / (Y + 1 / ZS)                           # Norton source U_g = 1
    U_slot = P / Zs_branch
    W = 0.5 * np.abs(U_slot) ** 2 * RHO * C * (OM / C) ** 2 / (2 * np.pi)
    return P, W, Zres


def metrics(g):
    P, W, Zres = network(g)
    WdB = 10 * np.log10(np.maximum(W, 1e-30))
    cents, rings, lvls = 0.0, 0.0, []
    for i, ft in enumerate(FT):
        # resonance of branch i in situ: peak of branch velocity |P/Zi|
        Ui = np.abs(P / Zres[i])
        win = (FRQ > ft * 0.7) & (FRQ < ft * 1.35)
        if not win.any():
            return 1e6, 1e6, -1e6, None
        j = np.argmax(np.where(win, Ui, 0))
        f_i = FRQ[j]
        half = Ui[j] / np.sqrt(2)
        lo = j
        while lo > 0 and Ui[lo] > half:
            lo -= 1
        hi = j
        while hi < len(FRQ) - 1 and Ui[hi] > half:
            hi += 1
        bw = max(FRQ[hi] - FRQ[lo], 0.5)
        Qi = f_i / bw
        t60 = 2.2 * Qi / f_i
        cents += abs(1200 * np.log2(f_i / ft))
        rings += abs(t60 - 0.30) / 0.30
        lvls.append(np.interp(ft, FRQ, WdB))
    return cents, rings, float(np.mean(lvls)), WdB


def evaluate(X):
    Fo = []
    for g in X:
        c_, r_, l_, _ = metrics(g)
        Fo.append([c_, r_, -l_])
    return np.array(Fo), np.zeros(len(X))


def nds(F):
    N = len(F)
    dom = np.zeros((N, N), bool)
    for i in range(N):
        dom[i] = (F[i][None, :] <= F).all(1) & (F[i][None, :] < F).any(1)
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
    n, m = len(idx), F.shape[1]
    dist = np.zeros(n)
    for k in range(m):
        o = np.argsort(F[idx, k]); f = F[idx, k][o]
        span = max(f[-1] - f[0], 1e-12)
        dist[o[0]] = dist[o[-1]] = np.inf
        dist[o[1:-1]] += (f[2:] - f[:-2]) / span
    return dist


POP, NGEN, NV = 64, 100, 10
X = XL + RNG.random((POP, NV)) * (XU - XL)
F, _ = evaluate(X)
for gen in range(NGEN):
    fronts, rank = nds(F)
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
            msk = RNG.random(NV) < 0.15
            um = RNG.random(NV)
            dl = np.where(um < 0.5, (2 * um) ** (1 / 21) - 1, 1 - (2 * (1 - um)) ** (1 / 21))
            kids.append(np.clip(np.where(msk, ch + dl * (XU - XL), ch), XL, XU))
    Xk = np.array(kids[:POP])
    Fk, _ = evaluate(Xk)
    Xa, Fa = np.vstack([X, Xk]), np.vstack([F, Fk])
    fronts, ra = nds(Fa)
    cra = np.zeros(len(Fa))
    for fr in fronts:
        cra[fr] = crowding(Fa, fr)
    keep = np.array(sorted(range(len(Fa)), key=lambda i: (ra[i], -cra[i]))[:POP])
    X, F = Xa[keep], Fa[keep]

fronts, rank = nds(F)
pf = np.array(fronts[0])
Fp, Xp = F[pf], X[pf]
print(f"cycle 11 (HALO v2) front: {len(pf)} of {POP}")
print(f"tuning error: {Fp[:,0].min():.0f}..{Fp[:,0].max():.0f} cents total | "
      f"ring error: {Fp[:,1].min():.2f}..{Fp[:,1].max():.2f} | level: "
      f"{(-Fp[:,2]).min():.1f}..{(-Fp[:,2]).max():.1f} dB")
# references: closed-form hand Halo, and cycle-8 evolved (shared necks)
Xhand = np.array([412., 183., 103., 66., 10., 10., 10., 10., 30.0, 8.0])
ch_, rh_, lh_, Wh = metrics(Xhand)
print(f"hand HALO (closed-form): tuning {ch_:.0f} cents, ring err {rh_:.2f}, level {lh_:.1f} dB")
Xc8 = np.array([402., 349., 162., 102., 11.4, 11.4, 11.4, 11.4, 22.0, 6.6])
c8c, c8r, c8l, _ = metrics(Xc8)
print(f"c8 evolved (shared necks): tuning {c8c:.0f} cents, ring err {c8r:.2f}, level {c8l:.1f} dB")
best = pf[np.argmin(Fp[:, 0] + 30 * Fp[:, 1])]
cb_, rb_, lb_, Wb = metrics(X[best])
print(f"EVOLVED HALO v2: V = {np.round(X[best,:4],0)} cm3, necks a = "
      f"{np.round(X[best,4:8],1)} mm, L_n {X[best,8]:.0f}, slot {X[best,9]:.1f} mm")
print(f"  tuning {cb_:.0f} cents, ring err {rb_:.2f}, level {lb_:.1f} dB")
print(f"  neck radii ascend with pitch? {np.all(np.diff(X[best,4:8]) > -0.5)}")
np.savetxt("pareto_c11.csv", np.hstack([Xp, Fp]), delimiter=",",
           header="V1,V2,V3,V4,a1,a2,a3,a4,L_n,slot,cents,ring_err,neg_level", comments="")

BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.8),
                        gridspec_kw=dict(left=0.06, right=0.98, top=0.80, bottom=0.12,
                                         wspace=0.24))
ax = axs[0]
sc = ax.scatter(Fp[:, 0], Fp[:, 1], c=-Fp[:, 2], cmap="viridis", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("level at chord (dB)")
ax.plot(ch_, rh_, "*", ms=17, color="#6D28D9", mec="#1F2937", mew=0.8)
ax.annotate("hand HALO\n(closed-form Law 14.6)", (ch_, rh_), (ch_ + 8, rh_ + 0.12),
            fontsize=8, arrowprops=dict(arrowstyle="->", color="#5c584d", lw=0.8))
ax.plot(F[best, 0], F[best, 1], "o", ms=10, mfc="none", mec="#6D28D9", mew=2.2)
ax.text(F[best, 0] + 4, F[best, 1], "EVOLVED", fontsize=8, color="#6D28D9")
ax.set_xlabel("total tuning error (cents)")
ax.set_ylabel("ring-time error  sum|t60 - 0.3|/0.3")
ax.set_title("cycle-8 front: in-situ chord tuning vs the shared halo", fontsize=9.3, loc="left")
ax.grid(alpha=0.15)
ax = axs[1]
ax.plot(FRQ, Wh - Wh.max(), color="#8a8578", lw=1.6, ls="--",
        label=f"hand (closed-form): {ch_:.0f} cents off")
ax.plot(FRQ, Wb - Wb.max(), color="#6D28D9", lw=2.0,
        label=f"evolved in situ: {cb_:.0f} cents off")
for ft in FT:
    ax.axvline(ft, color="#B45309", lw=0.8, ls=":")
ax.text(FT[0] - 42, 2, "targets", fontsize=7.4, color="#B45309")
ax.set_ylim(-45, 5)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("radiated power (dB re max)")
ax.legend(fontsize=7.6, frameon=False)
ax.set_title("the chord, in situ: coupled network spectrum", fontsize=9.3, loc="left")
ax.grid(alpha=0.15)
fig.suptitle("THE FORGE, cycle 11  -  HALO v2: per-sphere necks chase the shared ring time",
             fontsize=12, fontweight="bold", color="#1F2937", x=0.06, ha="left")
fig.savefig("evolution_c11.png", dpi=150)
print("wrote pareto_c11.csv, evolution_c11.png")
