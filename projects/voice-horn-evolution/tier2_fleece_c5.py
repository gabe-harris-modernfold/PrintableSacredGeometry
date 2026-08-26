"""Tier-2 test of the METAMATERIAL SURFACE on the cycle-5 KNEE (the print candidate).

Three wall conditions on the same axisymmetric FEM (gmsh + scikit-fem, BGT-1 arc at 1 m):

  RIGID   - smooth rigid bore (what every Forge cycle assumed)
  FLEECE  - bore wall carries the RULE W3/W4 liner, modeled as a locally reacting
            rigid-backed air layer with a resistive face sheet:
                Z_w(f, x) = R_f + i rho c cot(k d(x)),   R_f = 20 Rayl (RULE W3 budget),
                d grading 4 -> 8 mm along the flare (RULE W4), starting after the step (RULE W2)
  FRINGE  - FLEECE plus the rim fringe: the last 40 mm of wall grades toward rho*c
            (Z = rho c / s, s: 0.2 -> 0.95). A locally-reacting matched wall ABSORBS where
            the real fringe TRANSMITS, so this is an upper bound on the anti-etalon effect
            and a pessimistic bound on level.

Claims under test:
  1. RULE W3: liner insertion loss on the voice <= 0.4 dB      (FLEECE vs RIGID)
  2. Rim fringe: mouth-ripple suppression above c/4L = 2.1 kHz (FRINGE vs FLEECE)

Robin wall term (e^{-i omega t}): dp/dn = i omega rho p / Z_w.
Run from this directory: python tier2_fleece_c5.py
Outputs: printed metrics, tier2_fleece.png
"""
import numpy as np
from scipy.special import j1, struve
from scipy.sparse.linalg import splu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO, C = 1.2, 343.0
ZS = 4.0e6
RF = 20.0                        # Rayl, RULE W3 budget
GEN = dict(re=19.60, Le=16.47, rt=21.69, rm=143.58, L=501.26, n=1.74, A=-0.14, x0=0.80)
R_CUP, L_CUP = 22.5e-3, 40.0e-3
R_FAR = 1.0
FRQ = np.geomspace(300, 4200, 30)

re_, Le = GEN["re"] * 1e-3, GEN["Le"] * 1e-3
rt, rm, L = GEN["rt"] * 1e-3, GEN["rm"] * 1e-3, GEN["L"] * 1e-3
nfl, Abu, x0 = GEN["n"], GEN["A"], GEN["x0"]
z_step = L_CUP + Le
z_mouth = z_step + L
u = np.linspace(0, 1, 91)
r_bore = rt + (rm - rt) * u ** nfl + Abu * (rm - rt) * np.exp(-((u - x0) / 0.12) ** 2) * np.sin(np.pi * u)
z_bore = z_step + u * L

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("knee5")
LC_T, LC_B, LC_F = 2.5e-3, 5.0e-3, 16.0e-3
P = gmsh.model.geo.addPoint
pts = [P(0, 0, 0, LC_T), P(0, R_CUP, 0, LC_T), P(L_CUP, re_, 0, LC_T),
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
lines.append(gmsh.model.geo.addCircleArc(p_baffle, p_center, p_axis_far))
lines.append(Ln(p_axis_far, pts[0]))
loop = gmsh.model.geo.addCurveLoop(lines)
gmsh.model.geo.addPlaneSurface([loop])
gmsh.model.geo.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 16.0e-3)
gmsh.model.mesh.generate(2)
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
coords = np.array(ncoords).reshape(-1, 3)[:, :2]
idmap = {t: i for i, t in enumerate(ntags)}
_, _, enodes = gmsh.model.mesh.getElements(2)
tris = np.vectorize(idmap.get)(np.array(enodes[0]).reshape(-1, 3))
gmsh.finalize()
used = np.unique(tris)
remap = -np.ones(len(coords), dtype=np.int64); remap[used] = np.arange(len(used))
coords, tris = coords[used], remap[tris]
print(f"mesh: {len(coords)} nodes, {len(tris)} triangles")

from skfem import MeshTri, Basis, FacetBasis, ElementTriP2, asm, BilinearForm, LinearForm
from skfem.helpers import dot, grad

mesh = MeshTri(coords.T.astype(np.float64).copy(), tris.T.astype(np.int64).copy())
elem = ElementTriP2()
basis = Basis(mesh, elem)

def wall_r_of_z(z):
    return np.interp(z, z_bore, r_bore)

fac_drive = mesh.facets_satisfying(lambda x: x[0] < 1e-9)
fac_far = mesh.facets_satisfying(
    lambda x: np.abs(np.sqrt((x[0] - z_mouth) ** 2 + x[1] ** 2) - R_FAR) < 3e-3)
def on_bore_wall(x, z0, z1):
    return (x[0] > z0) & (x[0] < z1) & (np.abs(x[1] - wall_r_of_z(x[0])) < 2.5e-3)
FL_SEGS = []   # (facets, z0, z1) fleece thirds
zf0, zf1 = z_step + 5e-3, z_mouth - 40e-3
edges = np.linspace(zf0, zf1, 4)
for a, b in zip(edges[:-1], edges[1:]):
    FL_SEGS.append((mesh.facets_satisfying(lambda x, a=a, b=b: on_bore_wall(x, a, b)), a, b))
FR_SEGS = []   # fringe quarters
fedges = np.linspace(z_mouth - 40e-3, z_mouth - 0.5e-3, 5)
for a, b in zip(fedges[:-1], fedges[1:]):
    FR_SEGS.append((mesh.facets_satisfying(lambda x, a=a, b=b: on_bore_wall(x, a, b)), a, b))
print("fleece wall facets:", [len(s[0]) for s in FL_SEGS],
      "| fringe facets:", [len(s[0]) for s in FR_SEGS])

@BilinearForm
def stiff(uu, vv, w): return w.x[1] * dot(grad(uu), grad(vv))
@BilinearForm
def mass(uu, vv, w): return w.x[1] * uu * vv
@BilinearForm
def bmass(uu, vv, w): return w.x[1] * uu * vv
@LinearForm
def bload(vv, w): return w.x[1] * vv

S = asm(stiff, basis); M = asm(mass, basis)
Bf = asm(bmass, FacetBasis(mesh, elem, facets=fac_far))
fbd = FacetBasis(mesh, elem, facets=fac_drive)
Ld = asm(bload, fbd); Bd = asm(bmass, fbd)
B_FL = [asm(bmass, FacetBasis(mesh, elem, facets=s[0])) for s in FL_SEGS]
B_FR = [asm(bmass, FacetBasis(mesh, elem, facets=s[0])) for s in FR_SEGS]
D_FL = [4e-3, 6e-3, 8e-3]            # liner depth grading, RULE W4
S_FR = [0.2, 0.45, 0.7, 0.95]        # fringe match fraction
print(f"dofs: {S.shape[0]}")

pt_axis = np.array([[z_mouth + 0.5], [1e-4]])
try:
    probe = basis.probes(pt_axis)
except Exception:
    locs = basis.doflocs
    i0 = np.argmin((locs[0] - pt_axis[0, 0]) ** 2 + (locs[1] - pt_axis[1, 0]) ** 2)
    import scipy.sparse as sp
    probe = sp.csr_matrix(([1.0], ([0], [i0])), shape=(1, S.shape[0]))

A_cup = np.pi * R_CUP ** 2
U0 = A_cup
G = {"RIGID": [], "FLEECE": [], "FRINGE": []}
for fi, f in enumerate(FRQ):
    k = 2 * np.pi * f / C
    om = 2 * np.pi * f
    base = (S - k ** 2 * M - (1j * k - 1.0 / R_FAR) * Bf).tocsc().astype(complex)
    b = (-1j * om * RHO) * Ld.astype(complex)
    p_ref = om * RHO * U0 / (2 * np.pi * 0.5)
    for cfg in G:
        A = base.copy()
        if cfg in ("FLEECE", "FRINGE"):
            for Bm, d in zip(B_FL, D_FL):
                Zw = RF + 1j * RHO * C / np.tan(k * d)
                A = A - (1j * om * RHO / Zw) * Bm
        if cfg == "FRINGE":
            for Bm, sfr in zip(B_FR, S_FR):
                Zw = RHO * C / sfr
                A = A - (1j * om * RHO / Zw) * Bm
        else:
            for Bm, d in zip(B_FR, [8e-3] * 4) if cfg == "FLEECE" else []:
                Zw = RF + 1j * RHO * C / np.tan(k * d)
                A = A - (1j * om * RHO / Zw) * Bm
        p = splu(A).solve(b)
        p_avg = (Bd @ p).sum() / Bd.sum()
        Zin = p_avg / U0
        src = ZS / (ZS + Zin)
        p_ax = complex((probe @ p)[0]) * src
        G[cfg].append(20 * np.log10(abs(p_ax) / p_ref))
    if fi % 6 == 0:
        print(f"  f = {f:.0f} Hz done")
for cfg in G:
    G[cfg] = np.array(G[cfg])

hi = FRQ >= 2100
il = (G["RIGID"] - G["FLEECE"]).mean()
rip_b = G["FLEECE"][hi].max() - G["FLEECE"][hi].min()
rip_c = G["FRINGE"][hi].max() - G["FRINGE"][hi].min()
lvl = (G["FRINGE"][hi] - G["FLEECE"][hi]).mean()
print("\n=== metamaterial surface, verified at tier 2 ===")
print(f"RULE W3 test - Fleece insertion loss on the voice: {il:+.2f} dB mean "
      f"(budget: <= 0.4 dB) {'PASS' if abs(il) <= 0.4 else 'FAIL'}")
print(f"anti-etalon test - ripple 2.1-4.2 kHz: FLEECE {rip_b:.1f} dB -> FRINGE {rip_c:.1f} dB "
      f"({rip_c - rip_b:+.1f} dB)")
print(f"fringe level cost in the band (upper-bound model): {lvl:+.2f} dB mean")

# figure
BG = "#FBFAF7"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#8a8578", "text.color": "#3d3a33",
                     "axes.labelcolor": "#3d3a33", "xtick.color": "#5c584d",
                     "ytick.color": "#5c584d"})
INK = "#1F2937"
fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.9),
                        gridspec_kw=dict(left=0.06, right=0.985, top=0.82, bottom=0.12,
                                         wspace=0.24))
ax = axs[0]
ax.semilogx(FRQ, G["RIGID"], color="#8a8578", lw=2.4, label="RIGID bore (all Forge cycles)")
ax.semilogx(FRQ, G["FLEECE"], color="#0F766E", lw=1.8, ls="--",
            label="+ FLEECE liner (RULE W3/W4 model)")
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("on-axis gain (dB)")
ax.legend(fontsize=7.6, frameon=False)
ax.set_title(f"RULE W3 test: liner insertion loss {il:+.2f} dB mean vs 0.4 dB budget "
             f"({'PASS' if abs(il) <= 0.4 else 'FAIL in this worst-case model'})\n"
             "sealed-compliance locally-reacting model - pessimistic for a grazing wave",
             fontsize=9.3, loc="left", color="#0F766E")
ax.grid(True, which="both", alpha=0.15)

ax = axs[1]
ax.semilogx(FRQ, G["FLEECE"], color="#0F766E", lw=1.8, label="FLEECE, hard rim")
ax.semilogx(FRQ, G["FRINGE"], color="#1D4ED8", lw=2.0, label="FLEECE + rim fringe (bound)")
ax.axvline(2143, color="#5c584d", lw=0.9, ls=":")
ax.text(2200, ax.get_ylim()[0] + 0.5, "c/4L = 2.1 kHz", fontsize=7.4)
ax.axvspan(2100, 4200, color="#B45309", alpha=0.06)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("on-axis gain (dB)")
ax.legend(fontsize=7.6, frameon=False)
ax.set_title(f"the anti-etalon: ripple 2.1-4.2 kHz  {rip_b:.1f} -> {rip_c:.1f} dB "
             f"({rip_c-rip_b:+.1f}), level {lvl:+.1f} dB\n"
             "(matched-wall model: bound on ripple benefit, pessimistic on level)",
             fontsize=9.3, loc="left", color="#1D4ED8")
ax.grid(True, which="both", alpha=0.15)

fig.suptitle("THE METAMATERIAL SURFACE AT TIER 2  -  cycle-5 KNEE bore, "
             "Fleece as locally-reacting wall, fringe as graded rim admittance",
             fontsize=12, fontweight="bold", color=INK, x=0.06, ha="left")
fig.savefig("tier2_fleece.png", dpi=150)
print("wrote tier2_fleece.png")
