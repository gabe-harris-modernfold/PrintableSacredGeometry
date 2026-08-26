"""Forge cycle 6 - EVOLVE THE INTERIOR TEXTURE (the Fleece) on the fixed cycle-5 KNEE bore.

The bore is frozen; the liner is the genome:
    s0    in [0.00, 0.40]  liner start, fraction of bore length (wind mandate: cover >= 60%)
    d0    in [3, 10] mm    liner depth at start   (wind mandate: >= 3 mm everywhere)
    d1    in [3, 12] mm    liner depth at mouth end of Fleece zone
    R_f   in [5, 60] Rayl  face-sheet flow resistance
    phi   in [0.85, 0.97]  porosity (cavity compliance ~ phi)

Wall model (locally reacting, rigid-backed - the PESSIMISTIC model that failed RULE W3
for the naive liner): Z_w = R_f + i (rho c / phi) cot(k d(x)).
Tier-1.5 lined-duct propagation in covered segments:
    k_z^2 = k^2 - 2 i k beta / a ,  beta = rho c / Z_w ,  Zc' = (rho c / S) (k / k_z)
with the cycle-5 directivity cap Q_eff = min(Q_piston, Q_cone).

Objectives (3): -G_speech, ripple(>=500 Hz), Gamma_cup.
After the search: tier-2 FEM verification of the knee pick vs the naive W4 liner vs rigid.

Run from this directory: python evolve_fleece_c6.py
Outputs: pareto_c6.csv, evolution_c6.png
"""
import numpy as np
from scipy.special import j1, struve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6
GEN = dict(re=19.60, Le=16.47, rt=21.69, rm=143.58, L=501.26, n=1.74, A=-0.14, x0=0.80)
R_CUP, L_CUP, NC = 22.5, 40.0, 16
NT, NH = 12, 140
FRQ = np.geomspace(300, 4200, 24)
WB = np.interp(np.log(FRQ), np.log([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500,
                                    3150, 4000]),
               [0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2])
WB = WB / WB.sum()
K = 2 * np.pi * FRQ / C
UH = np.linspace(0, 1, NH + 1)
XL = np.array([0.00, 3.0, 3.0, 5.0, 0.85])
XU = np.array([0.40, 10.0, 12.0, 60.0, 0.97])
RNG = np.random.default_rng(7)

# fixed chain geometry (mm)
re_, Le, rt, rm, L = GEN["re"], GEN["Le"], GEN["rt"], GEN["rm"], GEN["L"]
nfl, Abu, x0 = GEN["n"], GEN["A"], GEN["x0"]
r_cup = np.linspace(R_CUP, re_, NC + 1)
r_h = rt + (rm - rt) * UH ** nfl + Abu * (rm - rt) * np.exp(-((UH - x0) / 0.12) ** 2) * np.sin(np.pi * UH)
mid = lambda r: 0.5 * (r[:-1] + r[1:])
r_all = np.concatenate([mid(r_cup), np.full(NT, re_), mid(r_h)])          # mm, per segment
S_all = np.pi * (r_all / 1000.0) ** 2
d_all = np.concatenate([np.full(NC, L_CUP / NC), np.full(NT, Le / NT),
                        np.full(NH, L / NH)]) / 1000.0
u_seg = np.concatenate([np.full(NC, -1.0), np.full(NT, -1.0),
                        (np.arange(NH) + 0.5) / NH])                       # bore fraction
a_m = r_h[-1] / 1000.0
A_cup = np.pi * (R_CUP / 1000.0) ** 2
slope = (r_h[-1] - r_h[-2]) / (L / NH)
QCONE = 2.0 / (1.0 - np.cos(np.arctan(slope)))
x2 = 2 * K * a_m
R1 = 1 - 2 * j1(x2) / x2
X1 = 2 * struve(1, x2) / x2
ZRAD = RHO * C / (np.pi * a_m ** 2) * (R1 + 1j * X1)
QP = np.minimum((K * a_m) ** 2 / np.maximum(R1, 1e-9), QCONE)
WREF = 0.5 * RHO * C * K ** 2 / (4 * np.pi)


def evaluate(X):
    """X: (N,5) liner genomes. Tier-1.5 lined-duct chain, vectorized over freq."""
    out_F, out_G = [], []
    for g in X:
        s0, d0, d1, Rf, phi = g
        p, Uv = ZRAD.astype(complex).copy(), np.ones_like(ZRAD, dtype=complex)
        for s in range(len(S_all) - 1, -1, -1):
            a_loc = r_all[s] / 1000.0
            if u_seg[s] >= s0:
                frac = (u_seg[s] - s0) / max(1 - s0, 1e-6)
                d = (d0 + (d1 - d0) * frac) / 1000.0
                Zw = Rf + 1j * (RHO * C / phi) / np.tan(K * d)
                beta = RHO * C / Zw
                kz = np.sqrt(K ** 2 - 2j * K * beta / a_loc)
                kz = np.where(kz.imag > 0, kz, np.conj(kz))   # decay, not growth
                Zc = (RHO * C / S_all[s]) * (K / kz)
                kd = kz * d_all[s]
            else:
                Zc = RHO * C / S_all[s]
                kd = K * d_all[s]
            cs, sn = np.cos(kd), np.sin(kd)
            p, Uv = cs * p + 1j * Zc * sn * Uv, 1j * sn / Zc * p + cs * Uv
        Zin = p / Uv
        src = ZS / (ZS + Zin)
        W = 0.5 * np.real(ZRAD) * np.abs(src / Uv) ** 2
        G = 10 * np.log10(np.maximum(W / WREF, 1e-12)) + 10 * np.log10(QP / 2.0)
        Gs = (G * WB).sum()
        rip = G[FRQ >= 500].max() - G[FRQ >= 500].min()
        Gam = np.abs((Zin - RHO * C / A_cup) / (Zin + RHO * C / A_cup)).mean()
        out_F.append([-Gs, rip, Gam])
        out_G.append(G)
    return np.array(out_F), np.zeros(len(X)), np.array(out_G)


def nds(F, cv):
    N = len(F)
    dom = np.zeros((N, N), bool)
    for i in range(N):
        dom[i] = (F[i][None, :] <= F).all(1) & (F[i][None, :] < F).any(1)
    fronts, rank = [], np.full(N, -1)
    nd = dom.sum(0)
    cur = np.where(nd == 0)[0].tolist()
    fr = 0
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


POP, NGEN = 48, 50
X = XL + RNG.random((POP, 5)) * (XU - XL)
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
        uu = RNG.random(a.shape)
        beta = np.where(uu <= 0.5, (2 * uu) ** (1 / 16), (1 / (2 * (1 - uu))) ** (1 / 16))
        for ch in (0.5 * ((1 + beta) * a + (1 - beta) * b),
                   0.5 * ((1 - beta) * a + (1 + beta) * b)):
            mask = RNG.random(5) < 0.25
            um = RNG.random(5)
            delta = np.where(um < 0.5, (2 * um) ** (1 / 21) - 1, 1 - (2 * (1 - um)) ** (1 / 21))
            kids.append(np.clip(np.where(mask, ch + delta * (XU - XL), ch), XL, XU))
    Xk = np.array(kids[:POP])
    Fk, cvk, _ = evaluate(Xk)
    Xall, Fall = np.vstack([X, Xk]), np.vstack([F, Fk])
    fronts, rank_all = nds(Fall, np.zeros(len(Fall)))
    crowd_all = np.zeros(len(Fall))
    for fr in fronts:
        crowd_all[fr] = crowding(Fall, fr)
    keep = np.array(sorted(range(len(Fall)), key=lambda i: (rank_all[i], -crowd_all[i]))[:POP])
    X, F = Xall[keep], Fall[keep]

fronts, rank = nds(F, np.zeros(POP))
pf = np.array(fronts[0])
Fp, Xp = F[pf], X[pf]
print(f"cycle 6 front: {len(pf)} of {POP}")
print(f"G: {(-Fp[:,0]).min():.1f}..{(-Fp[:,0]).max():.1f} dB | ripple: "
      f"{Fp[:,1].min():.1f}..{Fp[:,1].max():.1f} dB | Gam: {Fp[:,2].min():.3f}..{Fp[:,2].max():.3f}")

# references in the same tier-1.5 model: rigid (liner off ~ d->3mm? no: use s0 -> 1 unreachable;
# emulate rigid by Rf=60, phi=0.85, d=3 minimal? Proper rigid: bypass liner entirely)
F_rigid, _, G_rigid = evaluate(np.array([[1.1, 3, 3, 60, 0.9]]))   # s0>1: liner never engages
F_naive, _, G_naive = evaluate(np.array([[0.0, 4, 8, 20, 0.90]]))  # the hand-specified W4 liner
print(f"rigid   : G={-F_rigid[0,0]:.1f}  rip={F_rigid[0,1]:.1f}")
print(f"naive W4: G={-F_naive[0,0]:.1f}  rip={F_naive[0,1]:.1f}  "
      f"(IL={-F_rigid[0,0]+F_naive[0,0]:+.2f} dB)")

gn = (-Fp[:, 0] - (-Fp[:, 0]).min()) / max(np.ptp(-Fp[:, 0]), 1e-9)
rn = (Fp[:, 1].max() - Fp[:, 1]) / max(np.ptp(Fp[:, 1]), 1e-9)
best = pf[np.argmax(np.minimum(gn, rn))]
s0b, d0b, d1b, Rfb, phib = X[best]
print(f"EVOLVED TEXTURE: start {s0b:.2f}L, depth {d0b:.1f}->{d1b:.1f} mm, "
      f"R_f {Rfb:.0f} Rayl, phi {phib:.2f}")
print(f"  G={-F[best,0]:.1f} dB (IL={-F_rigid[0,0]+F[best,0]:+.2f}), rip={F[best,1]:.1f} dB "
      f"(rigid {F_rigid[0,1]:.1f}), Gam={F[best,2]:.3f}")
np.savetxt("pareto_c6.csv", np.hstack([Xp, -Fp[:, :1], Fp[:, 1:]]), delimiter=",",
           header="s0,d0_mm,d1_mm,Rf_rayl,phi,G_speech_dB,ripple_dB,refl", comments="")

# ---------------- tier-2 FEM verification of the evolved texture ----------------
print("\ntier-2 FEM verification (rigid vs naive vs evolved)...")
import subprocess, json, io, os
# reuse the tier2_fleece machinery inline (mesh once)
exec(io.open("tier2_fleece_c5.py", encoding="utf-8").read().split("A_cup = np.pi * R_CUP ** 2")[0]
     .replace('FRQ = np.geomspace(300, 4200, 30)', 'FRQ_FEM = np.geomspace(300, 4200, 24)')
     .replace("import matplotlib\nmatplotlib.use(\"Agg\")\n", "")
     .replace('GEN = dict(re=19.60', 'GEN = dict(re=19.60'), globals())
from scipy.sparse.linalg import splu
A_cupF = np.pi * R_CUP ** 2          # R_CUP now metres from the exec'd script
U0 = A_cupF
configs = {"RIGID": None,
           "NAIVE": (0.0, 4e-3, 8e-3, 20.0, 0.90),
           "EVOLVED": (s0b, d0b * 1e-3, d1b * 1e-3, Rfb, phib)}
GF = {}
for cfg, prm in configs.items():
    Gc = []
    for f in FRQ_FEM:
        k = 2 * np.pi * f / C
        om = 2 * np.pi * f
        A = (S - k ** 2 * M - (1j * k - 1.0 / R_FAR) * Bf).tocsc().astype(complex)
        if prm is not None:
            s0_, dA, dB_, Rf_, ph_ = prm
            for Bm, (fs_, za, zb) in zip(B_FL + B_FR, FL_SEGS + FR_SEGS):
                zmid = 0.5 * (za + zb)
                ufrac = (zmid - z_step) / (z_mouth - z_step)
                if ufrac < s0_:
                    continue
                dloc = dA + (dB_ - dA) * (ufrac - s0_) / max(1 - s0_, 1e-6)
                Zw = Rf_ + 1j * (RHO * C / ph_) / np.tan(k * dloc)
                A = A - (1j * om * RHO / Zw) * Bm
        b = (-1j * om * RHO) * Ld.astype(complex)
        p = splu(A).solve(b)
        p_avg = (Bd @ p).sum() / Bd.sum()
        src = ZS / (ZS + p_avg / U0)
        p_ax = complex((probe @ p)[0]) * src
        Gc.append(20 * np.log10(abs(p_ax) / (om * RHO * U0 / (2 * np.pi * 0.5))))
    GF[cfg] = np.array(Gc)
    print(f"  {cfg}: G_speech(FEM) = {(GF[cfg]*WB).sum():.1f} dB")
ilN = ((GF["RIGID"] - GF["NAIVE"]) * WB).sum()
ilE = ((GF["RIGID"] - GF["EVOLVED"]) * WB).sum()
hi = FRQ_FEM >= 500
ripR = GF["RIGID"][hi].max() - GF["RIGID"][hi].min()
ripN = GF["NAIVE"][hi].max() - GF["NAIVE"][hi].min()
ripE = GF["EVOLVED"][hi].max() - GF["EVOLVED"][hi].min()
print(f"FEM insertion loss: naive {ilN:+.2f} dB | evolved {ilE:+.2f} dB (RULE W3 budget 0.4)")
print(f"FEM ripple >=500Hz: rigid {ripR:.1f} | naive {ripN:.1f} | evolved {ripE:.1f} dB")

# ---------------- figure ----------------
BGc = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BGc, "axes.facecolor": BGc, "savefig.facecolor": BGc,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
INK = "#1F2937"
fig = plt.figure(figsize=(11.8, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.28, left=0.07, right=0.97,
                      top=0.87, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(Fp[:, 1], -Fp[:, 0], c=Fp[:, 3] if Fp.shape[1] > 3 else Xp[:, 3],
                cmap="viridis", s=42, edgecolors="none")
fig.colorbar(sc, ax=ax, pad=0.02).set_label("face resistance R_f (Rayl)")
ax.plot(F_rigid[0, 1], -F_rigid[0, 0], "s", ms=10, color="#8a8578")
ax.text(F_rigid[0, 1] + 0.1, -F_rigid[0, 0], "rigid bore", fontsize=7.6)
ax.plot(F_naive[0, 1], -F_naive[0, 0], "X", ms=11, color="#7a1f1f")
ax.text(F_naive[0, 1] + 0.1, -F_naive[0, 0], "naive W4 liner", fontsize=7.6, color="#7a1f1f")
ax.plot(F[best, 1], -F[best, 0], "o", ms=11, mfc="none", mec="#0F766E", mew=2.4)
ax.text(F[best, 1] + 0.1, -F[best, 0] - 0.35, "EVOLVED", fontsize=8, color="#0F766E")
ax.set_xlabel("spectral ripple >= 500 Hz (dB)")
ax.set_ylabel("G_speech (dB, tier-1.5 lined)")
ax.set_title("cycle-6 front: the liner's gain/smoothness trade", fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[0, 1])
ub = np.linspace(0, 1, 200)
dn = np.where(ub >= 0.0, 4 + 4 * ub, np.nan)
ax.plot(ub * GEN["L"], dn, color="#7a1f1f", lw=1.8, ls="--",
        label="naive W4: 4->8 mm, R 20, phi 0.90, full bore")
de = np.where(ub >= s0b, d0b + (d1b - d0b) * (ub - s0b) / max(1 - s0b, 1e-6), np.nan)
ax.plot(ub * GEN["L"], de, color="#0F766E", lw=2.4,
        label=f"evolved: start {s0b:.2f}L, {d0b:.1f}->{d1b:.1f} mm, R {Rfb:.0f}, phi {phib:.2f}")
ax.set_xlabel("position along bore (mm)"); ax.set_ylabel("liner depth d(x) (mm)")
ax.set_ylim(0, 13)
ax.legend(fontsize=7.4, frameon=False)
ax.set_title("what evolution did to the texture", fontsize=9.5, loc="left")
ax.grid(alpha=0.15)

ax = fig.add_subplot(gs[1, 0])
ax.semilogx(FRQ_FEM, GF["RIGID"], color="#8a8578", lw=2.2, label="rigid (FEM)")
ax.semilogx(FRQ_FEM, GF["NAIVE"], color="#7a1f1f", lw=1.6, ls="--",
            label=f"naive liner (FEM, IL {ilN:+.2f} dB)")
ax.semilogx(FRQ_FEM, GF["EVOLVED"], color="#0F766E", lw=2.0,
            label=f"evolved liner (FEM, IL {ilE:+.2f} dB)")
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("on-axis gain (dB)")
ax.legend(fontsize=7.4, frameon=False)
ax.set_title("tier-2 verification of the evolved texture", fontsize=9.5, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
cats = ["rigid", "naive W4", "EVOLVED"]
ils = [0.0, ilN, ilE]
rips = [ripR, ripN, ripE]
xpos = np.arange(3)
b1 = ax.bar(xpos - 0.18, ils, width=0.34, color="#7a1f1f", alpha=0.75,
            label="insertion loss (dB)")
b2 = ax.bar(xpos + 0.18, rips, width=0.34, color="#1D4ED8", alpha=0.65,
            label="ripple >=500 Hz (dB)")
ax.axhline(0.4, color="#7a1f1f", lw=1.0, ls=":")
ax.text(2.35, 0.5, "RULE W3\nbudget", fontsize=7, color="#7a1f1f")
ax.set_xticks(xpos); ax.set_xticklabels(cats)
ax.legend(fontsize=7.6, frameon=False)
ax.set_title("the verdict, at FEM fidelity (pessimistic wall model)", fontsize=9.5, loc="left")
ax.grid(alpha=0.15, axis="y")

fig.suptitle("THE FORGE, cycle 6  -  the interior metamaterial texture evolves\n"
             f"liner genome on the frozen cycle-5 KNEE; {POP*NGEN+POP} tier-1.5 lined-duct "
             "evaluations + FEM verification",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("evolution_c6.png", dpi=150)
print("wrote pareto_c6.csv, evolution_c6.png")
