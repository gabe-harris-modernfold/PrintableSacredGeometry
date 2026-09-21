#!/usr/bin/env python3
"""Lowest self-excited tone from candle heat: linear onset of a Rijke tube.

Model (Heckl 1990; Matveev & Culick 2003 energy form). Open-open tube, length L,
radius R, gauze heater at x_h = L/4. Fundamental p'(x,t) = P sin(kx) cos(wt),
u'(x,t) = -(P/rho c) cos(kx) sin(wt).

Driving (Rayleigh):  heat-release fluctuation follows the acoustic velocity at the
heater with King's-law sensitivity and a time lag tau (Lighthill: tau ~ 0.2 d/u0):
    Q'(t) = (chi Q0 / 2 u0) u'(x_h, t - tau)
    W_gain = (gamma-1)/(rho c^2) * <p' Q'>  = (gamma-1) P^2 chi Q0 sin(w tau) / (4 rho^2 c^3 u0)
    (positive for x_h in the lower half of the tube -- the Rijke condition)

Stored energy:  E = P^2 V / (4 rho c^2),  V = pi R^2 L

Losses:
  wall (viscous + thermal Stokes layers):  alpha = (1/(R c)) sqrt(nu w / 2) (1 + (gamma-1)/sqrt(Pr))
      W_wall = 2 alpha c E
  radiation, two unflanged open ends, kR << 1:  W_rad = P^2 pi R^2 (kR)^2 / (4 rho c)
  heater grid pressure-drop damping is neglected -> the result is an optimistic floor.

Mean draft u0 through the gauze from the chimney balance
    rho u0 A cp dT = Q0 ,  (1/2) K rho u0^2 = rho g L dT / T_inf   (K = grid + end losses)
so u0 = (2 g L Q0 / (K rho A cp T_inf))^(1/3).

Onset when W_gain >= W_wall + W_rad; P^2 cancels, leaving a threshold heater power
Q0_min(f, R). For each frequency the radius is optimised, and the lowest frequency
whose threshold fits under the candle budget is reported.
"""
import numpy as np

GAMMA, RHO0, C0, T_INF, PR, NU0 = 1.4, 1.20, 343.0, 293.0, 0.71, 1.5e-5
CP, G = 1005.0, 9.81
CHI = 0.6          # velocity-sensitive fraction of the gauze heat transfer (King's law, Re ~ 10)
K_LOSS = 2.5       # chimney flow-loss coefficient (grid + inlet + outlet)
D_WIRE = 0.4e-3    # gauze wire, m


def threshold_power(f, R, T_gas=420.0, d_wire=D_WIRE, tau_fixed=None, chi=CHI):
    """Heater power [W] at onset for fundamental f in a tube of radius R.
    Solves W_gain(Q0) = W_loss exactly (u0 and tau both depend on Q0 through the draft).
    tau_fixed: use a heater with a fixed response lag (e.g. a flame) instead of Lighthill's 0.2 d/u0.
    Returns (Q0_min, L, u0, tau)."""
    from scipy.optimize import brentq
    c = C0 * np.sqrt(T_gas / T_INF)
    rho = RHO0 * T_INF / T_gas
    nu = NU0 * (T_gas / T_INF) ** 1.7
    L = c / (2 * f)
    w = 2 * np.pi * f
    k = w / c
    A = np.pi * R ** 2
    alpha = (1 / (R * c)) * np.sqrt(nu * w / 2) * (1 + (GAMMA - 1) / np.sqrt(PR))
    E_per_P2 = A * L / (4 * rho * c ** 2)
    loss = 2 * alpha * c * E_per_P2 + A * (k * R) ** 2 / (4 * rho * c)

    def u0_of(Q0):
        return (2 * G * L * Q0 / (K_LOSS * rho * A * CP * T_INF)) ** (1 / 3)

    def tau_of(Q0):
        return tau_fixed if tau_fixed is not None else 0.2 * d_wire / u0_of(Q0)

    def gain(Q0):
        return (GAMMA - 1) * chi * Q0 * np.sin(np.minimum(w * tau_of(Q0), np.pi / 2)) / (4 * rho ** 2 * c ** 3 * u0_of(Q0))

    g = lambda lq: np.log(gain(np.exp(lq))) - np.log(loss)
    try:
        lq = brentq(g, np.log(1e-3), np.log(1e12))
    except ValueError:
        return np.inf, L, np.nan, np.nan
    Q0 = np.exp(lq)
    return Q0, L, u0_of(Q0), tau_of(Q0)


def best_radius(f, **kw):
    Rs = np.linspace(0.01, 0.15, 141)
    res = [threshold_power(f, R, **kw) for R in Rs]
    i = int(np.argmin([r[0] for r in res]))
    return Rs[i], res[i]


def table(title, freqs, **kw):
    print(title)
    print(f"{'f [Hz]':>7} {'L [m]':>6} {'R_opt [mm]':>10} {'u0 [m/s]':>9} {'tau [ms]':>9} {'Q0_min [W]':>12}")
    for f in freqs:
        R, (Q, L, u0, tau) = best_radius(f, **kw)
        print(f"{f:7.0f} {L:6.2f} {R*1e3:10.0f} {u0:9.2f} {tau*1e3:9.2f} {Q:12.4g}")
    print()


def f_min(budget, **kw):
    fs = np.linspace(10, 600, 296)
    Q = np.array([best_radius(f, **kw)[1][0] for f in fs])
    ok = fs[Q <= budget]
    return ok.min() if len(ok) else None


if __name__ == "__main__":
    F = (20, 30, 50, 75, 100, 150, 200, 300, 500)
    table(f"A. classic gauze heater, {D_WIRE*1e3:.1f} mm wire, lag = 0.2 d/u0 (Lighthill), chi={CHI}, K={K_LOSS}, gas 420 K", F)
    table("B. 3 mm rod heater (longer convective lag)", F, d_wire=3e-3)
    table("C. heater with a fixed 8 ms response lag (a diffusion flame's own time constant)", F, tau_fixed=8e-3)
    print("lowest self-excited frequency under a heater budget:")
    print(f"{'budget [W]':>11} {'A gauze':>9} {'B 3mm rod':>10} {'C flame 8ms':>12}")
    for budget in (30, 60, 100, 200, 500, 2000):
        row = [f_min(budget), f_min(budget, d_wire=3e-3), f_min(budget, tau_fixed=8e-3)]
        print(f"{budget:11.0f} " + " ".join(f"{(x if x else float('nan')):9.0f}" if i == 0 else f"{(x if x else float('nan')):10.0f}" if i == 1 else f"{(x if x else float('nan')):12.0f}" for i, x in enumerate(row)))
