#!/usr/bin/env python3
"""Quartz Flash Windscrew -- the numbers behind the build.

Run from this directory:  python physics.py
Prints the design table and writes physics.png.

Chain of models (assumptions in README.md):
  1. Rotor     Liam-type conical Archimedes screw, linear torque-coefficient curve
               Ct(lam) = Ct0 (1 - lam/lam0), conservative small-rotor values.
  2. Tumbler   loose quartz stones in an icosahedron on the shaft.  Froude number
               Fr = omega^2 r / g says when they cascade (Fr < 1) and when they are
               pinned to the walls (Fr > 1) -- which is every wind above cut-in.
  3. Paddle    a stationary paddle on a fixed sleeve inside the icosahedron meets each
               pinned stone once per revolution.  Strike + friction re-acceleration
               costs ~ eta m (omega r)^2 / 2 per stone per rev, so the load power goes
               as omega^3, like the wind's: the operating tip-speed ratio is the same
               at every wind speed.
  4. Impact    Hertz sphere-on-sphere impact: peak force, contact radius, peak
               pressure, and the rim tension that starts a cone crack.
  5. Sliding   flash temperature of a sliding contact, Tian & Kennedy (1994)
               circular-source formula with Blok heat partition; smooth Hertz contact
               and plastic asperities (p = H) both evaluated.
  6. Photons   Planck emission from the hot real-contact area, scotopic-weighted, into
               a dark-adapted pupil, against the Hecht-Shlaer-Pirenne threshold -- swept
               over emissivity, because the visible emissivity of a micron-thick hot
               skin of clear quartz is the most uncertain number here.
"""
import math
import numpy as np

RHO_AIR = 1.20             # kg/m^3 (ridge, ~200 m ASL, 15 C)
G = 9.81
H_PL, C_L, K_B = 6.626e-34, 2.998e8, 1.381e-23

# rotor (matches windscrew.py)
R_ROT = 0.150
A_ROT = math.pi * R_ROT**2
CT0, LAM0 = 0.20, 2.0      # -> Cp_max = Ct0 lam0 / 4 = 0.10 (independent Liam tests: ~0.1-0.2)
T_BEARING = 4e-3           # N m, four shielded 608s

# tumbler (matches windscrew.py)
N_STONES = 6               # one per outer vertex pocket of the face-axis icosahedron
M_STONE = 0.012            # kg (a ~20 mm tumbled quartz stone)
R_POCKET = 0.055           # m, centroid radius of a pinned stone (icosahedron Rc 70)
ETA = 1.5                  # strike (~1) + friction re-acceleration (~0.5) per 1/2 m v^2
REPOSE = math.radians(35)  # charge angle of repose against the paddle at rest
SLIDE_LEN = 0.010          # m, scrape length per paddle pass

# quartz (alpha, room temperature, directional averages)
E_Q, NU_Q = 95e9, 0.08
K_Q = 7.0
KAPPA_Q = K_Q / (2650 * 740)
H_Q = 10e9
MU_Q = 0.5
SIGMA_F = 70e6
R_STONE = 0.009
T_AMB = 288.0
T_ALPHA_BETA = 846.0
T_CEIL = 1950.0

# eye
PUPIL_D = 7e-3
VIEW_L = 0.30              # m, eye to the icosahedron at night
T_PETG = 0.85              # transmission of one clear 2.4 mm PETG wall
HSP = (54, 148)

_SC = np.array([[400, 9.29e-3], [420, 9.66e-2], [440, 0.328], [460, 0.567], [480, 0.793],
                [500, 0.982], [507, 1.0], [520, 0.935], [540, 0.650], [560, 0.329],
                [580, 0.121], [600, 3.32e-2], [620, 7.37e-3], [640, 1.50e-3],
                [660, 3.13e-4], [680, 7.15e-5], [700, 1.78e-5]])

def ct(lam):
    return CT0 * np.clip(1 - lam / LAM0, 0, None)

def rotor_torque(v, lam):
    return 0.5 * RHO_AIR * A_ROT * R_ROT * v**2 * ct(lam)

def strike_energy(omega):
    return ETA * 0.5 * M_STONE * (omega * R_POCKET)**2

def load_coeff():
    return N_STONES * ETA * M_STONE * R_POCKET**2 / (4 * math.pi * R_ROT**3)

def operating_lambda():
    """n eta m r^2 lam^2 / (4 pi R^3) = 1/2 rho A Ct0 (1 - lam/lam0); v has cancelled."""
    a, b = load_coeff(), 0.5 * RHO_AIR * A_ROT * CT0
    return (-(b / LAM0) + math.sqrt((b / LAM0)**2 + 4 * a * b)) / (2 * a)

def operating_point(v):
    lam = np.linspace(1e-3, LAM0, 4000)
    om = lam * v / R_ROT
    load = N_STONES * strike_energy(om) / (2 * math.pi) + T_BEARING
    net = rotor_torque(v, lam) - load
    i = np.where(np.diff(np.sign(net)) < 0)[0]
    if len(i) == 0:
        return 0.0, 0.0
    return lam[i[0]], lam[i[0]] * v / R_ROT

def charge_torque():
    return N_STONES * M_STONE * G * R_POCKET * math.sin(REPOSE) * 0.6

def cut_in():
    return math.sqrt((T_BEARING + charge_torque()) / (0.5 * RHO_AIR * A_ROT * R_ROT * CT0))

def hertz_impact(v_n, m=M_STONE, r1=R_STONE, r2=R_STONE):
    e_star = E_Q / (2 * (1 - NU_Q**2))
    r_star = r1 * r2 / (r1 + r2)
    d = (15 * m * v_n**2 / (16 * e_star * math.sqrt(r_star)))**0.4
    f = 4 / 3 * e_star * math.sqrt(r_star) * d**1.5
    a = math.sqrt(r_star * d)
    p0 = 3 * f / (2 * math.pi * a**2)
    return f, a, p0, (1 - 2 * NU_Q) / 3 * p0

def hertz_static(n, r1=R_STONE, r2=R_STONE):
    e_star = E_Q / (2 * (1 - NU_Q**2))
    return (3 * n * (r1 * r2 / (r1 + r2)) / (4 * e_star))**(1 / 3)

def flash_rise(q, a, v):
    """Tian & Kennedy (1994), uniform circular source, flux q, Pe = v a / 2 kappa on the
    moving body; Blok partition against a body that sees the source stationary."""
    pe = v * a / (2 * KAPPA_Q)
    t_moving = 2 * q * a / (K_Q * math.sqrt(math.pi * (1.273 + pe)))
    t_still = 2 * q * a / (K_Q * math.sqrt(math.pi * 1.273))
    return 1 / (1 / t_moving + 1 / t_still)

def flash_T(q, a, v):
    return min(T_AMB + flash_rise(q, a, v), T_CEIL)

def photons_to_eye(T, area, tau, eps, dist=VIEW_L):
    lam_nm = np.linspace(400, 700, 301)
    lam = lam_nm * 1e-9
    mp = 2 * math.pi * C_L / lam**4 / np.expm1(H_PL * C_L / (lam * K_B * T))
    w = np.interp(lam_nm, _SC[:, 0], _SC[:, 1])
    n_scot = np.trapezoid(mp * w, lam) * area * tau * eps
    frac = (PUPIL_D**2 / 4) / dist**2 * T_PETG                  # Lambertian, on-axis
    return n_scot * frac

def slide_state(v):
    lam, om = operating_point(v)
    vs = om * R_POCKET
    fc = M_STONE * om**2 * R_POCKET
    return om, vs, fc

if __name__ == "__main__":
    lam_op = operating_lambda()
    w_c = math.sqrt(G / R_POCKET)
    print("== rotor + tumbler load ==")
    print(f"rotor swept area {A_ROT*1e4:.0f} cm2, Cp_max (model) {CT0*LAM0/4:.3f}")
    print(f"operating tip-speed ratio (drag-free, wind-independent): {lam_op:.2f}  (runaway {LAM0:.1f})")
    print(f"cut-in wind (bearings + stones resting on the paddle): {cut_in():.1f} m/s")
    print(f"Froude 1 (stones pin to the walls) at omega = {w_c:.1f} rad/s = {w_c*30/math.pi:.0f} rpm"
          f" = {w_c*R_ROT/lam_op:.1f} m/s wind  -> without the paddle nothing collides above that")
    print()
    print(" wind   lam   rpm   Fr   strikes/s  stone m/s  E_strike mJ  P_load W  "
          "F_pin N  Hertz F N  p0 GPa  rim tension/strength")
    for v in (4, 6, 8, 10, 12, 15, 20, 25):
        lam, om = operating_point(v)
        vs = om * R_POCKET
        e = strike_energy(om)
        f, a, p0, st = hertz_impact(vs)
        print(f"{v:5.0f} {lam:5.2f} {om*30/math.pi:5.0f} {om**2*R_POCKET/G:5.1f} "
              f"{N_STONES*om/(2*math.pi):9.1f} {vs:10.2f} {e*1e3:11.1f} "
              f"{N_STONES*om/(2*math.pi)*e:9.2f} {M_STONE*om**2*R_POCKET:8.1f} {f:10.0f} "
              f"{p0/1e9:7.2f} {st/SIGMA_F:10.0f}x")

    print("\n== sliding flash temperature (quartz on quartz; quartz on PETG cannot run hot) ==")
    print("  v m/s   Hertz a um  T_smooth K   T_asperity(a=3um) K   T_asperity(a=10um) K")
    for v in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
        n = 1.5
        a_h = hertz_static(n)
        q_h = MU_Q * n / (math.pi * a_h**2) * v
        q_asp = MU_Q * H_Q * v
        print(f"{v:7.1f} {a_h*1e6:11.0f} {flash_T(q_h, a_h, v):11.0f} "
              f"{flash_T(q_asp, 3e-6, v):20.0f} {flash_T(q_asp, 10e-6, v):21.0f}")
    print(f"  (hand rubbing ~0.2 m/s; alpha-beta inversion {T_ALPHA_BETA:.0f} K; ceiling {T_CEIL:.0f} K)")

    print(f"\n== thermal photons at the cornea per scrape, 8 m/s wind, eye {VIEW_L*100:.0f} cm away ==")
    om, vs, fc = slide_state(8)
    area, tau = fc / H_Q, SLIDE_LEN / vs
    T = flash_T(MU_Q * H_Q * vs, 3e-6, vs)
    print(f"  real contact {area*1e12:.0f} um2 at {T:.0f} K for {tau*1e3:.1f} ms")
    for eps in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        ns = photons_to_eye(T, area, tau, eps)
        flag = "seen" if ns > HSP[1] else ("threshold" if ns > HSP[0] else "below")
        print(f"  emissivity {eps:7.0e}: {ns:10.3g} scotopic-equivalent photons  ({flag})")
    eps_need = next(e for e in np.logspace(-8, 0, 800) if photons_to_eye(T, area, tau, e) > HSP[0])
    print(f"  emissivity needed for the HSP floor: {eps_need:.1e}  -> a ~3 um hot skin must absorb "
          f">= {eps_need/3e-6/100:.1f} /cm in the visible")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    winds = np.linspace(2, 25, 200)
    ops = [operating_point(v) for v in winds]
    ax[0, 0].plot(winds, [o[1] * 30 / math.pi for o in ops], color="#3a6ea5")
    ax[0, 0].set_xlabel("wind m/s"); ax[0, 0].set_ylabel("rpm", color="#3a6ea5")
    axb = ax[0, 0].twinx()
    axb.plot(winds, [o[1]**2 * R_POCKET / G for o in ops], color="#c0552b")
    axb.axhline(1, color="#c0552b", ls=":"); axb.set_ylabel("Froude number (pinned above 1)", color="#c0552b")
    ax[0, 0].axvline(cut_in(), color="0.6", ls="--"); ax[0, 0].text(cut_in() + 0.3, 50, "cut-in", fontsize=8)
    ax[0, 0].set_title("Operating point: stones are pinned at every running speed")

    lam = np.linspace(0, LAM0, 200)
    ax[0, 1].plot(lam, ct(lam) * lam, color="#3a6ea5", label="rotor Cp(lam)")
    ax[0, 1].plot(lam, load_coeff() * lam**3 / (0.5 * RHO_AIR * A_ROT), color="#c0552b",
                  label="paddle demand / (1/2 rho A v^3)")
    ax[0, 1].axvline(lam_op, color="0.5", ls=":")
    ax[0, 1].set_ylim(0, 0.15); ax[0, 1].set_xlabel("tip-speed ratio")
    ax[0, 1].legend(); ax[0, 1].set_title("Matching: wind speed cancels")

    vv = np.logspace(-1, 1, 100)
    for a_asp, ls in ((1e-6, "--"), (3e-6, "-"), (10e-6, "-.")):
        ax[1, 0].plot(vv, [flash_T(MU_Q * H_Q * v, a_asp, v) for v in vv], ls,
                      color="#c0552b", label=f"asperity a={a_asp*1e6:.0f} um")
    a_h = hertz_static(1.5)
    ax[1, 0].plot(vv, [flash_T(MU_Q * 1.5 / (math.pi * a_h**2) * v, a_h, v) for v in vv],
                  color="#3a6ea5", label="smooth Hertz, 1.5 N")
    ax[1, 0].axhline(T_ALPHA_BETA, color="0.5", ls=":"); ax[1, 0].text(0.11, 880, "alpha-beta 573 C", fontsize=8)
    ax[1, 0].axvspan(0.1, 0.3, color="0.9"); ax[1, 0].text(0.105, 1700, "hand\nrubbing", fontsize=8)
    lo, hi = slide_state(4)[1], slide_state(15)[1]
    ax[1, 0].axvspan(lo, hi, color="#f3d9c8", alpha=0.6); ax[1, 0].text(lo * 1.05, 400, "windscrew\n4-15 m/s", fontsize=8)
    ax[1, 0].set_xscale("log"); ax[1, 0].set_xlabel("sliding speed m/s"); ax[1, 0].set_ylabel("flash T, K")
    ax[1, 0].legend(fontsize=8); ax[1, 0].set_title("Flash temperature (Tian-Kennedy, capped at softening)")

    eps = np.logspace(-7, 0, 100)
    for v, col in ((5, "#7a9cc6"), (8, "#3a6ea5"), (15, "#1d3f66")):
        om, vs, fc = slide_state(v)
        T = flash_T(MU_Q * H_Q * vs, 3e-6, vs)
        ax[1, 1].plot(eps, [photons_to_eye(T, fc / H_Q, SLIDE_LEN / vs, e) for e in eps],
                      color=col, label=f"{v} m/s wind")
    ax[1, 1].axhspan(*HSP, color="#f3d9c8"); ax[1, 1].text(2e-7, 170, "HSP threshold", fontsize=8)
    ax[1, 1].set_xscale("log"); ax[1, 1].set_yscale("log")
    ax[1, 1].set_xlabel("visible emissivity of the hot skin")
    ax[1, 1].set_ylabel(f"photons at cornea / scrape ({VIEW_L*100:.0f} cm)")
    ax[1, 1].legend(fontsize=8); ax[1, 1].set_title("Thermal channel: seen only if the skin absorbs")
    fig.tight_layout()
    fig.savefig("physics.png", dpi=110)
    print("\nwrote physics.png")
