"""FEM check of the FLOOR texture - the thinnest lawful Fleece (d = 3 mm uniform,
phi = 0.85, R_f = 60 Rayl): the minimum-compliance corner the wind mandate permits.
Run from this directory: python tier2_floor.py
"""
import io
import numpy as np

src = io.open("tier2_fleece_c5.py", encoding="utf-8").read()
head = src.split("A_cup = np.pi * R_CUP ** 2")[0].replace(
    "FRQ = np.geomspace(300, 4200, 30)", "FRQ_FEM = np.geomspace(300, 4200, 24)")
exec(head, globals())
from scipy.sparse.linalg import splu

WB = np.interp(np.log(FRQ_FEM), np.log([315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
                                        2500, 3150, 4000]),
               [0.5, 0.5, 0.7, 0.7, 0.8, 1.0, 1.2, 1.5, 1.5, 1.5, 1.3, 1.2])
WB = WB / WB.sum()
U0 = np.pi * R_CUP ** 2
configs = {"RIGID": None, "FLOOR": (0.0, 3e-3, 3e-3, 60.0, 0.85)}
G = {}
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
                dloc = dA + (dB_ - dA) * ufrac
                Zw = Rf_ + 1j * (RHO * C / ph_) / np.tan(k * dloc)
                A = A - (1j * om * RHO / Zw) * Bm
        b = (-1j * om * RHO) * Ld.astype(complex)
        p = splu(A).solve(b)
        p_avg = (Bd @ p).sum() / Bd.sum()
        srcf = ZS / (ZS + p_avg / U0)
        p_ax = complex((probe @ p)[0]) * srcf
        Gc.append(20 * np.log10(abs(p_ax) / (om * RHO * U0 / (2 * np.pi * 0.5))))
    G[cfg] = np.array(Gc)
il = ((G["RIGID"] - G["FLOOR"]) * WB).sum()
hi = FRQ_FEM >= 500
ripR = G["RIGID"][hi].max() - G["RIGID"][hi].min()
ripF = G["FLOOR"][hi].max() - G["FLOOR"][hi].min()
print(f"FLOOR texture (3 mm / phi 0.85 / R 60): insertion loss {il:+.2f} dB "
      f"(RULE W3 budget 0.4 -> {'PASS' if abs(il) <= 0.4 else 'FAIL'})")
print(f"ripple: rigid {ripR:.1f} dB, floor {ripF:.1f} dB")
np.savetxt("floor_fem.csv", np.column_stack([FRQ_FEM, G["RIGID"], G["FLOOR"]]),
           delimiter=",", header="f_hz,G_rigid_dB,G_floor_dB", comments="")
