#!/usr/bin/env python3
"""Candle plume vs. PETG: temperature at the rotor and the frame, and a rotor speed
estimate. Run from this directory after pyramid.py:  python physics.py

Plume model: Heskestad's turbulent-plume correlations for each tea light, virtual
origin at the flame base, Gaussian radial profiles, four plumes superposed
linearly. A candle plume is laminar over the first ~15 cm, so the correlation is
a *lower* bound on centreline temperature there and about right once the plumes
have transitioned and begun to merge -- the regime the rotor and collar live in.
Measure before trusting it (see README).
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyramid as P

# ---- air and candle
T_INF = 293.0          # K
RHO, CP, G = 1.20, 1.005, 9.81       # cp in kJ/(kg K): Heskestad's constants take kW and kJ
Q_TEALIGHT = 40.0      # W total (14 g paraffin over ~4 h)
CHI_CONV = 0.70        # convective fraction (soot-poor flame)
QC = CHI_CONV * Q_TEALIGHT / 1000.0   # kW, per candle
Z_FLAME_BASE = P.PLATE_T + 16.0       # tea-light wax surface, world z (mm)

K_T = 9.1 * (T_INF / (G * CP**2 * RHO**2)) ** (1 / 3)      # Heskestad temperature constant
K_U = 3.4 * (G / (CP * RHO * T_INF)) ** (1 / 3)            # velocity constant

def plume(z_m, r_m):
    """(dT [K], w [m/s]) of one tea-light plume at height z (m above flame base), radius r."""
    z = np.maximum(z_m, 0.02)
    dT0 = K_T * QC ** (2 / 3) * z ** (-5 / 3)
    u0 = K_U * QC ** (1 / 3) * z ** (-1 / 3)
    b = 0.12 * np.sqrt((T_INF + dT0) / T_INF) * z            # Gaussian radius of dT
    dT = dT0 * np.exp(-(r_m / b) ** 2)
    w = u0 * np.exp(-(r_m / (1.1 * b)) ** 2)
    return dT, w

def field(x_mm, y_mm, z_world_mm):
    """Superpose the four candles. x,y,z in world mm."""
    z = (z_world_mm - Z_FLAME_BASE) / 1000.0
    dT = np.zeros(np.broadcast(x_mm, y_mm).shape); w = np.zeros_like(dT)
    for i in range(4):
        a = math.radians(90 * i)
        cx, cy = P.R_CANDLE * math.cos(a), P.R_CANDLE * math.sin(a)
        r = np.hypot(x_mm - cx, y_mm - cy) / 1000.0
        d, u = plume(z, r)
        dT += d; w += u
    return dT, w

if __name__ == "__main__":
    z_blade = P.rotor_hub_z() + P.SLOT_Z
    print(f"blade plane z = {z_blade:.0f} mm  ({z_blade - Z_FLAME_BASE:.0f} mm above wax surface)")

    # --- temperature map on the blade plane
    g = np.linspace(-130, 130, 261)
    X, Y = np.meshgrid(g, g)
    dT, W = field(X, Y, z_blade)
    R = np.hypot(X, Y)
    print(f"blade plane: peak {T_INF - 273 + dT.max():.0f} C over a candle, "
          f"hub (r<18) {T_INF - 273 + dT[R < 18].max():.0f} C, tip (r=120) {T_INF - 273 + dT[np.abs(R-120)<1].max():.0f} C")

    # --- temperature along the frame: collar and posts
    probes = {"mast top": (0, 0, P.COLLAR_Z1), "mast @ z=150": (0, 0, 150.0), "mast @ z=250": (0, 0, 250.0)}
    if P.FRAME == "posts":
        probes.update({"post @ z=150": tuple(P.post_point(0, 150)), "post @ z=250": tuple(P.post_point(0, 250))})
    for name, (x, y, z) in probes.items():
        d, _ = field(np.array(x, float), np.array(y, float), z)
        print(f"  {name:22s} {T_INF - 273 + float(d):5.0f} C  (Heskestad; laminar core may run hotter)")

    # --- rotor speed: blade-element balance with a flat plate, C_N = 1.6 sin(alpha)
    th = math.radians(P.PITCH_DEG)
    r = np.linspace(P.SLOT_INNER_R + 10, P.ROTOR_R, 200) / 1000.0
    chord = np.interp(r * 1000, [P.SLOT_INNER_R, 0.82 * P.ROTOR_R, P.ROTOR_R],
                      [P.ROOT_CHORD, P.TIP_CHORD, 0.3 * P.TIP_CHORD]) / 1000.0
    phi = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    ROTOR_MASS = 0.079                        # kg, solid hub + 8 blades (pyramid.py report)
    FRICTION = 0.3 * ROTOR_MASS * G * 0.05e-3  # mu * m g * contact-patch radius (~50 um dome-in-vee)
    I = P.N_BLADE * 5.6e-3 * (0.5 * (P.ROTOR_R / 1000) ** 2) + 0.5 * 34.5e-3 * (P.HUB_D / 2000) ** 2

    def rotor_at(zb):
        """(mean updraft, starting torque, steady omega, peak blade-plane temp) for blade plane at world z=zb."""
        w_r = np.array([field(rr * 1000 * np.cos(phi), rr * 1000 * np.sin(phi), zb)[1].mean() for rr in r])
        def net_torque(omega):
            v = omega * r
            V = np.hypot(w_r, v)
            alpha = np.arctan2(w_r, v) - th
            drive = 0.5 * RHO * V ** 2 * chord * 1.6 * np.sin(alpha) * math.sin(th) * r
            Re = np.maximum(V * chord / 1.5e-5, 50)
            drag = 0.5 * RHO * V ** 2 * chord * 2 * (1.33 / np.sqrt(Re)) * math.cos(th) * r
            return P.N_BLADE * np.trapezoid(drive - drag, r)
        om = np.linspace(0.01, 20, 1000)
        tq = np.array([net_torque(o) for o in om]) - FRICTION
        idx = np.where(tq < 0)[0]
        om_ss = om[idx[0]] if len(idx) else om[-1]
        peak = T_INF - 273 + float(field(np.array(P.R_CANDLE, float), np.array(0.0), zb)[0])
        return w_r.mean(), net_torque(0.0), om_ss, peak

    w_mean, T0, om_ss, _ = rotor_at(z_blade)
    print(f"updraft at blades: mean {w_mean:.2f} m/s, over a candle {W.max():.2f} m/s")
    print(f"starting torque {T0*1e3:.3f} mN.m vs pivot friction {FRICTION*1e3:.3f} mN.m")
    print(f"steady speed ~ {om_ss * 60 / (2 * math.pi):.0f} rpm; inertia {I*1e6:.0f} g.cm2, "
          f"spin-up time constant ~ {I * om_ss / T0:.1f} s")
    print()
    print("blade-height sweep (move Z_APEX in pyramid.py to change it):")
    print("  blades above wax   peak air T over candle   start torque   steady rpm")
    for zb_above in (160, 200, 240, 287, 320):
        w_m, t0, o, peak = rotor_at(Z_FLAME_BASE + zb_above)
        print(f"  {zb_above:8.0f} mm        {peak:6.0f} C                {t0*1e3:6.3f} mN.m   {o*60/(2*math.pi):5.0f}")

    # --- figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    im = ax[0].pcolormesh(X, Y, T_INF - 273 + dT, cmap="inferno", shading="auto", vmin=20, vmax=90)
    ax[0].add_patch(plt.Circle((0, 0), P.ROTOR_R, fill=False, ec="w", lw=1))
    ax[0].add_patch(plt.Circle((0, 0), P.HUB_D / 2, fill=False, ec="w", lw=1))
    for i in range(4):
        a = math.radians(90 * i)
        ax[0].add_patch(plt.Circle((P.R_CANDLE * math.cos(a), P.R_CANDLE * math.sin(a)), P.TEALIGHT_D / 2,
                                   fill=False, ec="c", lw=1, ls="--"))
    ax[0].set_aspect("equal"); ax[0].set_title(f"air temperature at the blade plane, z = {z_blade:.0f} mm")
    ax[0].set_xlabel("mm"); fig.colorbar(im, ax=ax[0], label="degC")
    zz = np.linspace(60, 340, 200)
    on_axis = [T_INF - 273 + float(field(np.array(0.0), np.array(0.0), z)[0]) for z in zz]
    over_candle = [T_INF - 273 + float(field(np.array(P.R_CANDLE, float), np.array(0.0), z)[0]) for z in zz]
    ax[1].plot(on_axis, zz, label="pyramid axis (mast, hub)")
    ax[1].plot(over_candle, zz, label="directly over a candle (blade path)")
    ax[1].axhline(z_blade, color="k", lw=0.8, ls=":"); ax[1].text(95, z_blade + 3, "blades")
    ax[1].axhline(P.COLLAR_Z1, color="0.6", lw=0.8); ax[1].text(95, P.COLLAR_Z1 - 12, "mast top / collar")
    ax[1].axvline(70, color="r", lw=1, ls="--"); ax[1].text(71, 80, "PETG softens", color="r", rotation=90)
    ax[1].set_xlim(20, 160); ax[1].set_xlabel("degC"); ax[1].set_ylabel("height above base, mm")
    ax[1].set_title("air temperature vs height (4 tea lights)"); ax[1].legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig("plume_temperature.png", dpi=130)
    print("wrote plume_temperature.png")
