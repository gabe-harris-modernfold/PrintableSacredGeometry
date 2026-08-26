"""Tier-2 verification of the Forge cycle-4 KNEE candidate.

2-D axisymmetric Helmholtz FEM (scikit-fem + gmsh) of the full chain - lip cup, epilaryngeal
tube, area step, bore - radiating through the mouth plane into a baffled half-space bounded
by a first-order (BGT-1) radiation arc at R_far = 1.0 m from the mouth center.

Physics (e^{-i omega t}): u = grad(p)/(i omega rho); drive plane z=0 with piston velocity
v0 = 1 m/s (Neumann g = -i omega rho v0 with outward normal); rigid walls natural;
far arc dp/dn = (ik - 1/R) p. Axisymmetric weight r in every form.

The Norton glottal source (Z_s = 4 MPa*s/m^3) is applied identically in both tiers through
the computed input impedance, so the comparison isolates the STRUCTURE model, not the source.

Outputs: printed band-by-band comparison, tier2_knee.png.
Run from this directory: python tier2_knee_c4.py
"""
import numpy as np
from scipy.special import j1, struve
from scipy.sparse.linalg import splu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6
FREQS = np.array([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000], float)
GENES = dict(re=20.0, Le=24.59, rt=24.12, rm=154.1, L=451.25, n=2.53, A=-0.05, x0=0.67)
R_CUP, L_CUP = 22.5e-3, 40.0e-3
R_FAR = 1.0

# ---------------- geometry (meridian, metres; x = z, y = r) ----------------
re_, Le = GENES["re"] * 1e-3, GENES["Le"] * 1e-3
rt, rm, L = GENES["rt"] * 1e-3, GENES["rm"] * 1e-3, GENES["L"] * 1e-3
nfl, Abu, x0 = GENES["n"], GENES["A"], GENES["x0"]
z_tube0, z_step = L_CUP, L_CUP + Le
z_mouth = z_step + L
u = np.linspace(0, 1, 91)
r_bore = rt + (rm - rt) * u ** nfl + Abu * (rm - rt) * np.exp(-((u - x0) / 0.12) ** 2) * np.sin(np.pi * u)
z_bore = z_step + u * L

# ---------------- mesh with gmsh ----------------
import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("knee")
LC_T, LC_B, LC_F = 2.5e-3, 5.0e-3, 16.0e-3
P = gmsh.model.geo.addPoint
pts = [P(0, 0, 0, LC_T), P(0, R_CUP, 0, LC_T), P(z_tube0, re_, 0, LC_T),
       P(z_step, re_, 0, LC_T), P(z_step, rt, 0, LC_T)]
bore_ids = [P(z_bore[i], r_bore[i], 0, LC_B) for i in range(1, len(u))]
p_baffle = P(z_mouth, R_FAR, 0, LC_F)
p_axis_far = P(z_mouth + R_FAR, 0, 0, LC_F)
p_center = P(z_mouth, 0, 0, LC_F)
Ln = gmsh.model.geo.addLine
lines = [Ln(pts[0], pts[1]), Ln(pts[1], pts[2]), Ln(pts[2], pts[3]), Ln(pts[3], pts[4])]
prev = pts[4]
for pid in bore_ids:
    lines.append(Ln(prev, pid)); prev = pid
lines.append(Ln(prev, p_baffle))
arc = gmsh.model.geo.addCircleArc(p_baffle, p_center, p_axis_far)
lines.append(arc)
lines.append(Ln(p_axis_far, pts[0]))
loop = gmsh.model.geo.addCurveLoop(lines)
surf = gmsh.model.geo.addPlaneSurface([loop])
gmsh.model.geo.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 16.0e-3)
gmsh.model.mesh.generate(2)
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
coords = np.array(ncoords).reshape(-1, 3)[:, :2]
idmap = {t: i for i, t in enumerate(ntags)}
etypes, etags, enodes = gmsh.model.mesh.getElements(2)
tris = np.array(enodes[0]).reshape(-1, 3)
tris = np.vectorize(idmap.get)(tris)
gmsh.finalize()
used = np.unique(tris)                      # drop nodes unused by any triangle
remap = -np.ones(len(coords), dtype=np.int64)
remap[used] = np.arange(len(used))
coords = coords[used]
tris = remap[tris]
print(f"mesh: {len(coords)} nodes, {len(tris)} triangles")

# ---------------- scikit-fem assembly ----------------
from skfem import MeshTri, Basis, FacetBasis, ElementTriP2, asm, BilinearForm, LinearForm
from skfem.helpers import dot, grad

mesh = MeshTri(coords.T.astype(np.float64).copy(), tris.T.astype(np.int64).copy())
elem = ElementTriP2()
basis = Basis(mesh, elem)
fac_drive = mesh.facets_satisfying(lambda x: x[0] < 1e-9)
fac_far = mesh.facets_satisfying(
    lambda x: np.abs(np.sqrt((x[0] - z_mouth) ** 2 + x[1] ** 2) - R_FAR) < 3e-3)
fb_drive = FacetBasis(mesh, elem, facets=fac_drive)
fb_far = FacetBasis(mesh, elem, facets=fac_far)

@BilinearForm
def stiff(uu, vv, w):
    return w.x[1] * dot(grad(uu), grad(vv))

@BilinearForm
def mass(uu, vv, w):
    return w.x[1] * uu * vv

@BilinearForm
def bmass(uu, vv, w):
    return w.x[1] * uu * vv

@LinearForm
def bload(vv, w):
    return w.x[1] * vv

S = asm(stiff, basis)
M = asm(mass, basis)
Bf = asm(bmass, fb_far)
Ld = asm(bload, fb_drive)
Bd = asm(bmass, fb_drive)
print(f"dofs: {S.shape[0]}")

# probes
def make_probe(points):
    try:
        return basis.probes(points)
    except Exception:
        locs = basis.doflocs
        idx = [np.argmin((locs[0] - p[0]) ** 2 + (locs[1] - p[1]) ** 2)
               for p in points.T]
        import scipy.sparse as sp
        Pm = sp.lil_matrix((points.shape[1], S.shape[0]))
        for row, i in enumerate(idx):
            Pm[row, i] = 1.0
        return Pm.tocsr()

pt_axis = np.array([[z_mouth + 0.5], [1e-4]])
probe_axis = make_probe(pt_axis)
theta = np.radians(np.arange(0, 91, 5))
pol_pts = np.vstack([z_mouth + 0.8 * np.cos(theta), np.maximum(0.8 * np.sin(theta), 1e-4)])
probe_pol = make_probe(pol_pts)

A_cup = np.pi * R_CUP ** 2
U0 = A_cup * 1.0  # v0 = 1 m/s
G2, DI2, GAM2, POL2 = [], [], [], []
for f in FREQS:
    k = 2 * np.pi * f / C
    om = 2 * np.pi * f
    A = (S - k ** 2 * M - (1j * k - 1.0 / R_FAR) * Bf).tocsc().astype(complex)
    b = (-1j * om * RHO * 1.0) * Ld.astype(complex)
    p = splu(A).solve(b)
    # input impedance: area-averaged drive pressure / volume velocity
    p_avg = (Bd @ p).sum() / Bd.sum()            # r-weighted mean over drive plane
    Zin = p_avg / U0                             # acoustic impedance, Pa/(m^3/s)
    src = ZS / (ZS + Zin)                         # Norton glottal delivery factor
    Z0 = RHO * C / A_cup
    GAM2.append(abs((Zin - Z0) / (Zin + Z0)))
    p_ax = complex((probe_axis @ p)[0]) * src
    p_ref = om * RHO * U0 / (2 * np.pi * 0.5)     # baffled monopole at 0.5 m, same U0
    G2.append(20 * np.log10(abs(p_ax) / p_ref))
    pol = np.abs(np.asarray(probe_pol @ p).ravel())
    POL2.append(pol / pol[0])
    Iavg = np.trapezoid(pol ** 2 * np.sin(theta), theta)
    DI2.append(10 * np.log10(2 * pol[0] ** 2 / max(Iavg, 1e-30)))
G2, DI2, GAM2 = np.array(G2), np.array(DI2), np.array(GAM2)
p_field = p  # 4 kHz field (last solve) for the plot

# ---------------- tier-1 (same chain, same source) for comparison ----------------
NSEG_C, NSEG_T, NSEG_H = 16, 12, 140
r_cup_t1 = np.linspace(R_CUP, re_, NSEG_C + 1)
uh = np.linspace(0, 1, NSEG_H + 1)
r_h_t1 = rt + (rm - rt) * uh ** nfl + Abu * (rm - rt) * np.exp(-((uh - x0) / 0.12) ** 2) * np.sin(np.pi * uh)
mid = lambda r: 0.5 * (r[:-1] + r[1:])
Ssegs = np.concatenate([np.pi * mid(r_cup_t1) ** 2, np.pi * re_ ** 2 * np.ones(NSEG_T),
                        np.pi * mid(r_h_t1) ** 2])
dsegs = np.concatenate([np.full(NSEG_C, L_CUP / NSEG_C), np.full(NSEG_T, Le / NSEG_T),
                        np.full(NSEG_H, L / NSEG_H)])
a_m = r_h_t1[-1]
G1, DI1, GAM1 = [], [], []
for f in FREQS:
    k = 2 * np.pi * f / C
    x2 = 2 * k * a_m
    R1 = 1 - 2 * j1(x2) / x2
    X1 = 2 * struve(1, x2) / x2
    Zrad = RHO * C / (np.pi * a_m ** 2) * (R1 + 1j * X1)
    p_, U_ = complex(Zrad), 1.0 + 0j
    for s in range(len(Ssegs) - 1, -1, -1):
        Zc = RHO * C / Ssegs[s]
        kd = k * dsegs[s]
        p_, U_ = np.cos(kd) * p_ + 1j * Zc * np.sin(kd) * U_, \
                 1j * np.sin(kd) / Zc * p_ + np.cos(kd) * U_
    Zin = p_ / U_
    src = ZS / (ZS + Zin)
    W = 0.5 * np.real(Zrad) * abs(src / U_) ** 2
    Wref = 0.5 * RHO * C * k ** 2 / (4 * np.pi)
    Q = (k * a_m) ** 2 / max(R1, 1e-9)
    G1.append(10 * np.log10(W / Wref) + 10 * np.log10(Q / 2))
    DI1.append(10 * np.log10(Q))                 # full-sphere DI; hemisphere limit = 3 dB
    Z0 = RHO * C / A_cup
    GAM1.append(abs((Zin - Z0) / (Zin + Z0)))
G1, DI1, GAM1 = np.array(G1), np.array(DI1), np.array(GAM1)

print("\nband   G_t1    G_t2    dG     DI_t1  DI_t2   Gam1   Gam2")
for i, f in enumerate(FREQS):
    print(f"{f:5.0f}  {G1[i]:6.1f}  {G2[i]:6.1f}  {G2[i]-G1[i]:+5.1f}  "
          f"{DI1[i]:5.1f}  {DI2[i]:5.1f}  {GAM1[i]:5.2f}  {GAM2[i]:5.2f}")
WF = np.array([0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2]); WF /= WF.sum()
print(f"\nG_speech tier-1 = {(G1*WF).sum():.1f} dB   tier-2 = {(G2*WF).sum():.1f} dB")
print(f"Gam_bar  tier-1 = {GAM1.mean():.3f}       tier-2 = {GAM2.mean():.3f}")

# ---------------- figure ----------------
BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
INK = "#1F2937"
fig = plt.figure(figsize=(11.8, 8.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28, left=0.07, right=0.97,
                      top=0.88, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
ax.semilogx(FREQS, G1, color="#8a8578", lw=1.8, ls="--", marker="s", ms=4,
            label="tier-1 (Webster TM + piston load)")
ax.semilogx(FREQS, G2, color="#7a1f1f", lw=2.2, marker="o", ms=4.5,
            label="tier-2 (axisym. FEM, true radiation)")
ax.axvspan(1000, 4000, color="#B45309", alpha=0.07)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("on-axis gain vs bare voice at 0.5 m (dB)")
ax.legend(fontsize=7.6, frameon=False)
ax.set_title("cycle-4 KNEE: tier-1 prediction vs tier-2 verification", fontsize=9.5, loc="left")
ax.grid(True, which="both", alpha=0.15)

ax = fig.add_subplot(gs[0, 1], projection="polar")
cols = {1000: "#c9b18a", 2000: "#c78a4b", 3150: "#b45309", 4000: "#7a1f1f"}
for i, f in enumerate(FREQS):
    if f in cols:
        dB = 20 * np.log10(np.maximum(POL2[i], 1e-3))
        ax.plot(theta, np.clip(dB, -30, 0) + 30, color=cols[f], lw=1.9, label=f"{f:.0f} Hz")
        ka = 2 * np.pi * f / C * a_m
        arg = ka * np.sin(theta)
        Dpist = np.where(np.abs(arg) < 1e-9, 1.0, 2 * j1(np.where(arg == 0, 1, arg)) / np.where(arg == 0, 1, arg))
        dBp = 20 * np.log10(np.maximum(np.abs(Dpist), 1e-3))
        ax.plot(theta, np.clip(dBp, -30, 0) + 30, color=cols[f], lw=0.9, ls=":")
ax.set_thetamin(0); ax.set_thetamax(90)
ax.set_rticks([10, 20, 30]); ax.set_yticklabels(["-20", "-10", "0 dB"], fontsize=7)
ax.legend(loc="lower left", fontsize=7, frameon=False, bbox_to_anchor=(-0.12, -0.12))
ax.set_title("tier-2 polars (solid) vs tier-1 piston model (dotted)", fontsize=9)

ax = fig.add_subplot(gs[1, :])
spl = 20 * np.log10(np.maximum(np.abs(p_field[:basis.mesh.p.shape[1]]), 1e-6))
tp = ax.tripcolor(mesh.p[0], mesh.p[1], mesh.t.T, spl, cmap="magma", shading="gouraud")
ax.plot(np.concatenate([[0, 0, z_tube0, z_step, z_step], z_bore]),
        np.concatenate([[0, R_CUP, re_, re_, rt], r_bore]), color="w", lw=1.2)
ax.plot([z_mouth, z_mouth], [rm, R_FAR], color="w", lw=1.2)
cb = fig.colorbar(tp, ax=ax, pad=0.01)
cb.set_label("|p| (dB re 1 Pa, v0 = 1 m/s)")
ax.set_xlim(-0.02, z_mouth + 0.75); ax.set_ylim(0, 0.55)
ax.set_aspect("equal")
ax.set_xlabel("z (m)"); ax.set_ylabel("r (m)")
ax.set_title("the 4 kHz field: cup, tube, step, bore, and the beam leaving the mouth "
             "(meridian half-plane)", fontsize=9.5, loc="left")

fig.suptitle("TIER 2  -  axisymmetric FEM verification of the cycle-4 KNEE "
             f"({S.shape[0]} dofs, gmsh + scikit-fem, BGT-1 arc at 1.0 m)",
             fontsize=12.5, fontweight="bold", color=INK, x=0.07, ha="left")
fig.savefig("tier2_knee.png", dpi=150)
print("wrote tier2_knee.png")
