#!/usr/bin/env python3
"""Singing flame (Higgins 1802): how low can a candle-class diffusion flame drive a tube?

Same linear energy balance as rijke_limit.py, but the heater is the flame itself.
Heat-release response to acoustic velocity at the flame is a flame transfer function
    Q'/Q0 = F(w) u'/u0 ,   F(w) = n0 * exp(-i w tau_d) / (1 + i w tau_c)
-- a gain n0, a pure convective delay tau_d (fuel-to-flame-front transport) and a
first-order low-pass with time constant tau_c (flame-area / mixing response). For a
buoyant laminar diffusion flame tau_c is of the order of the flame's own residence
time, which is also what sets the ~10 Hz flicker; tau_d is a few ms.

Rayleigh driving with the flame at x_h:
    W_gain = (gamma-1) P^2 Q0 |F| sin(phi) m(x_h) / (2 rho^2 c^3 u0)
    phi = w tau_d + atan(w tau_c)            (drive is positive only for 0 < phi < pi)
    m(x_h) = |sin(k x_h) cos(k x_h)| = |sin(2 k x_h)| / 2 for open-open (max 1/2 at L/4);
    for a quarter-wave tube (acoustically closed bottom, open top) m = |sin(2 k x_h)|/2 too,
    with L = c/4f and the flame at L/2 of that tube; radiation from one end only.

Losses: Stokes-layer wall dissipation, open-end radiation, and (switchable, on by default)
the acoustic energy the mean draft carries out of the open ends. Near a pressure node the
energy is all kinetic, <e> = P^2/(4 rho c^2), so each open end loses W = u0 A <e>; against
the stored energy E = P^2 A L/(4 rho c^2) that is a decay rate of u0/L per open end.
It is an order-of-magnitude term (the proper treatment is the Mach-number correction to
the open-end reflection coefficient plus the entropy wave the flame sheds), and it turns
out to be the same size as the wall loss for the long, wide, low tubes this study cares
about. Still neglected: flame-anchoring damping. The result remains an optimistic floor.

Mean draft: chimney balance with all of Q0 released into the gas,
    u0 = (2 g L Q0 / (K rho A cp T_inf))^(1/3).
Because W_gain ~ Q0 / u0 ~ Q0^(2/3) and the convective loss ~ u0 ~ Q0^(1/3), the onset
condition is solved in closed form (a quadratic in Q0^(1/3)); no root finder.

Tube radius has a floor set by the flame: a tea light needs R >= 25 mm, four of them
R >= 70 mm; a single 100 W wick flame fits in R >= 35 mm.

    python singing_flame.py            tables
    python singing_flame.py --figure   also writes singing_flame.png
"""
import sys
import numpy as np

GAMMA, RHO0, C0, T_INF, PR, NU0 = 1.4, 1.20, 343.0, 293.0, 0.71, 1.5e-5
CP, G = 1005.0, 9.81
K_LOSS = 2.0
CONV_LOSS = True      # include mean-flow convection of acoustic energy out of the open ends
R_GRID = np.linspace(0.02, 0.20, 181)


def threshold(f, R, tau_c, tau_d, n0=1.0, mode="open", T_gas=450.0, conv=None):
    """Flame power [W] at onset. R may be an array. Returns (Q0_min, L, u0, phi)."""
    conv = CONV_LOSS if conv is None else conv
    R = np.asarray(R, dtype=float)
    c = C0 * np.sqrt(T_gas / T_INF)
    rho = RHO0 * T_INF / T_gas
    nu = NU0 * (T_gas / T_INF) ** 1.7
    w = 2 * np.pi * f
    k = w / c
    L = c / (2 * f) if mode == "open" else c / (4 * f)
    n_ends = 2 if mode == "open" else 1
    A = np.pi * R ** 2
    alpha = (1 / (R * c)) * np.sqrt(nu * w / 2) * (1 + (GAMMA - 1) / np.sqrt(PR))
    E_per_P2 = A * L / (4 * rho * c ** 2)
    loss0 = 2 * alpha * c * E_per_P2 + n_ends * A * (k * R) ** 2 / (8 * rho * c)
    phi = w * tau_d + np.arctan(w * tau_c)
    if not (0 < phi < np.pi):
        return np.full_like(R, np.inf), L, np.full_like(R, np.nan), phi
    Fmag = n0 / np.sqrt(1 + (w * tau_c) ** 2)
    m = 0.5                                   # flame at the optimum station

    # u0 = a Q0^(1/3);  gain = b Q0^(2/3);  conv loss = n_ends A u0 /(4 rho c^2) = d Q0^(1/3)
    a = (2 * G * L / (K_LOSS * rho * A * CP * T_INF)) ** (1 / 3)
    b = (GAMMA - 1) * Fmag * np.sin(phi) * m / (2 * rho ** 2 * c ** 3 * a)
    d = n_ends * A * a / (4 * rho * c ** 2) if conv else 0.0
    # b x^2 - d x - loss0 = 0,  x = Q0^(1/3)
    x = (d + np.sqrt(d ** 2 + 4 * b * loss0)) / (2 * b)
    Q0 = x ** 3
    return Q0, L, a * x, phi


def best_R(f, R_min, **kw):
    Rs = R_GRID[R_GRID >= R_min]
    Q, L, u0, phi = threshold(f, Rs, **kw)
    i = int(np.argmin(Q))
    return Rs[i], (Q[i], L, u0[i], phi)


def onset_curve(fs, R_min, **kw):
    return np.array([best_R(f, R_min, **kw)[1][0] for f in fs])


def f_min(budget, R_min, fs=None, **kw):
    fs = np.linspace(5, 400, 396) if fs is None else fs
    Q = onset_curve(fs, R_min, **kw)
    ok = fs[Q <= budget]
    return (ok.min(), Q) if len(ok) else (None, Q)


def tables():
    c = C0 * np.sqrt(450 / 293)

    def fmt(f, mode):
        if not f:
            return f"{'no onset':>18}"
        L = c / (2 * f) if mode == "open" else c / (4 * f)
        return f"{f:6.0f} Hz {L:5.2f} m"

    print(f"mean-flow convective loss {'ON' if CONV_LOSS else 'OFF'}")
    print("1. Low-pass flame (buoyant laminar diffusion flame): tau_c response time, tau_d delay, n0 gain.")
    print("   Minimum onset power over all frequencies, open-open tube, R >= 20 mm.")
    print(f"{'tau_c [ms]':>10} {'tau_d [ms]':>10} {'n0':>4} {'min Q0 [W]':>11} {'at f [Hz]':>9}")
    fs = np.linspace(5, 300, 60)
    for tau_c in (5e-3, 10e-3, 20e-3):
        for tau_d in (1e-3, 3e-3, 6e-3):
            Q = onset_curve(fs, 0.02, tau_c=tau_c, tau_d=tau_d, n0=0.7)
            i = int(np.argmin(Q))
            print(f"{tau_c*1e3:10.0f} {tau_d*1e3:10.0f} {0.7:4.1f} {Q[i]:11.0f} {fs[i]:9.0f}")

    budgets = (40, 100, 160, 400, 1000)
    hdr = f"{'tau_d [ms]':>10} {'n0':>4} " + " ".join(f"{b:>15} W" for b in budgets)
    print()
    print("2. n-tau flame (pure delay, flat gain -- a gas-jet flame whose feed is acoustically modulated, Higgins' case).")
    print("   Lowest self-excited tone vs flame power, R >= 20 mm.")
    print(hdr)
    for tau_d in (3e-3, 5e-3, 8e-3, 12e-3):
        for n0 in (0.3, 1.0):
            _, Q = f_min(1e9, 0.02, tau_c=0.0, tau_d=tau_d, n0=n0, mode="open")
            fs_ = np.linspace(5, 400, 396)
            row = [fmt(fs_[Q <= b].min() if (Q <= b).any() else None, "open") for b in budgets]
            print(f"{tau_d*1e3:10.0f} {n0:4.1f} " + " ".join(row))

    print()
    print("3. Same n-tau flame, quarter-wave tube (acoustically closed bottom with a small air inlet, open top).")
    print(hdr)
    for tau_d in (5e-3, 8e-3, 12e-3):
        for n0 in (0.3, 1.0):
            _, Q = f_min(1e9, 0.02, tau_c=0.0, tau_d=tau_d, n0=n0, mode="quarter")
            fs_ = np.linspace(5, 400, 396)
            row = [fmt(fs_[Q <= b].min() if (Q <= b).any() else None, "quarter") for b in budgets]
            print(f"{tau_d*1e3:10.0f} {n0:4.1f} " + " ".join(row))


def figure(path="singing_flame.png"):
    """Onset power vs fundamental for the three heater models, radius optimised at every point."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import rijke_limit

    fs = np.geomspace(8, 400, 120)
    curves = []
    # gauze: rijke_limit with the 0.71 ms lag calibrated against the 1 m / 150 W textbook tube
    Qg = np.array([rijke_limit.best_radius(f, tau_fixed=0.71e-3)[1][0] for f in fs])
    curves.append(("gauze heater, 0.71 ms lag (Rijke)", Qg, "C0", "-"))
    curves.append(("candle-class flame, low-pass FTF (tau_c 10 ms, tau_d 3 ms, n 0.7)",
                   onset_curve(fs, 0.02, tau_c=10e-3, tau_d=3e-3, n0=0.7), "C3", "-"))
    curves.append(("gas flame, n-tau FTF (tau 5 ms, n 1), open-open",
                   onset_curve(fs, 0.02, tau_c=0.0, tau_d=5e-3, n0=1.0, mode="open"), "C2", "-"))
    curves.append(("gas flame, n-tau FTF, quarter-wave (closed bottom)",
                   onset_curve(fs, 0.02, tau_c=0.0, tau_d=5e-3, n0=1.0, mode="quarter"), "C2", "--"))
    curves.append(("gas flame, quarter-wave, convective loss off",
                   onset_curve(fs, 0.02, tau_c=0.0, tau_d=5e-3, n0=1.0, mode="quarter", conv=False), "C2", ":"))

    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    for label, Q, col, ls in curves:
        ax.loglog(fs, Q, color=col, ls=ls, lw=1.8 if ls == "-" else 1.4, label=label)
    for W, txt in ((40, "one tea light, 40 W"), (160, "four tea lights, 160 W"), (1000, "1 kW gas burner")):
        ax.axhline(W, color="0.5", lw=0.8, ls="-.")
        ax.text(fs[-1], W * 1.08, txt, ha="right", va="bottom", fontsize=8, color="0.35")
    ax.axvspan(20, 400, color="0.93", zorder=0)
    ax.text(21, 3.2e5, "audible", fontsize=8, color="0.4")
    ax.set_xlim(fs[0], fs[-1]); ax.set_ylim(10, 5e5)
    ax.set_xlabel("fundamental [Hz]  (open-open tube length = c/2f, gas at 450 K)")
    ax.set_ylabel("heater power at onset [W]")
    ax.set_title("Lowest tone the heat can drive: linear onset, radius optimised at each frequency")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    tables()
    if "--figure" in sys.argv:
        figure()
