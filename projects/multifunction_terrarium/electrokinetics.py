"""Charge separation in the drip system: what is real, what is not, and what it
costs the living volume.

Short version. Streaming POTENTIAL is not the kilovolt mechanism -- conductivity,
not geometry, sets its ceiling, and in terrarium water it is 2 mV. The kilovolts in
the literature come from Kelvin-dropper INDUCTION, which is self-amplifying and
does work here. But it collides head-on with the drop-maximising geometry: the
0.4 mm notch tips that make small drops are also 2000 corona points at 0.6 kV.
"""

import math
import params as P

E0, ER = 8.854e-12, 80.0
EPS = E0 * ER
ETA = 1.0e-3
FARADAY = 96485.0
ZETA_PETG = -0.030          # PET/PETG in neutral water, -20..-40 mV
E_CRIT = 3.0e6              # V/m, dry air

SIGMA = {"distilled": 5.5e-6, "condensate": 1.0e-3, "terrarium water": 5.0e-2}


def debye(ionic_M):
    return 0.304e-9 / math.sqrt(ionic_M)

def streaming_current(bore_mm, dPdL, zeta=ZETA_PETG):
    A = math.pi * (bore_mm * 1e-3 / 2) ** 2
    return EPS * abs(zeta) / ETA * A * dPdL

def streaming_potential(dP, sigma, zeta=ZETA_PETG):
    """Note the sigma in the denominator: this is why geometry cannot save it."""
    return EPS * abs(zeta) * dP / (ETA * sigma)

def rayleigh_limit(drop_d_mm):
    R = drop_d_mm * 1e-3 / 2
    return 8 * math.pi * math.sqrt(E0 * P.SIGMA * R ** 3)

def drop_capacitance(drop_d_mm):
    return 4 * math.pi * E0 * (drop_d_mm * 1e-3 / 2)

def kelvin_ceiling(drop_d_mm):
    """Induced charge q = C_drop.V. The drop disintegrates when q hits Rayleigh,
    so DROP SIZE sets the maximum voltage the machine can reach."""
    return rayleigh_limit(drop_d_mm) / drop_capacitance(drop_d_mm)

def corona_onset(tip_r_mm):
    return E_CRIT * tip_r_mm * 1e-3

def kelvin_growth(f_drops, drop_d_mm, C_collector):
    """Kelvin droppers grow exponentially, not linearly."""
    return f_drops * drop_capacitance(drop_d_mm) / C_collector

def ozone_per_day(V, I, volume_L, g_per_kWh=75.0):
    """Corona makes ozone. In a sealed living volume it accumulates."""
    kWh = V * I * 24 / 1000.0          # V*I = W; x24 h = Wh; /1000 = kWh
    mg = kWh * g_per_kWh * 1000.0
    mg_m3 = mg / (volume_L / 1000.0)
    return mg, mg_m3 / 2.14          # ~2.14 mg/m3 per ppm for O3


def report():
    D = P.DROP_D
    print("=" * 74)
    print("1. STREAMING POTENTIAL -- ruled out as the kV mechanism")
    dP = P.RHO * P.GRAV * 0.5
    for k, s in SIGMA.items():
        print(f"   {k:18s} sigma {s:8.1e} S/m  ->  V = {streaming_potential(dP,s):9.3f} V")
    print("   V = eps.zeta.dP/(eta.sigma). Conductivity is in the denominator, so no")
    print("   channel geometry rescues it. A living terrarium is the worst case.")

    print()
    print("2. STREAMING CURRENT -- small but usable as a SEED")
    Is = streaming_current(2.5, P.RHO * P.GRAV)
    print(f"   {Is*1e9:.2f} nA per 2.5 mm bore  x16 = {16*Is*1e9:.1f} nA")
    print("   BUT only the closed splitter tubes generate it. The spiral gutters are")
    print("   open C-sections with a free surface -- no confined dP, so ~no streaming.")

    print()
    print("3. KELVIN INDUCTION -- this is the real mechanism, and it self-amplifies")
    Cd = drop_capacitance(D)
    f = P.Q_TRICKLE / 3600 * 1e-3 / (P.DROP_V * 1e-9)
    for C in (5e-12, 20e-12):
        g = kelvin_growth(f, D, C)
        t = math.log(6000 / 1e-3) / g
        print(f"   C={C*1e12:4.0f} pF: growth {g:6.2f} /s (e-fold {1/g*1e3:5.0f} ms), "
              f"seed to 6 kV in {t:.1f} s")
    print(f"   drop capacitance {Cd*1e15:.1f} fF, {f:.0f} drops/s total")

    print()
    print("4. THE VOLTAGE CEILING IS SET BY DROP SIZE")
    for d in (1.5, P.DROP_D, 3.0, 4.8):
        print(f"   {d:4.2f} mm drop: Rayleigh {rayleigh_limit(d)*1e12:6.0f} pC "
              f"-> ceiling {kelvin_ceiling(d)/1000:5.2f} kV")
    print("   Above that the drop disintegrates and carries no more charge.")

    print()
    print("5. THE COLLISION: sharp tips make small drops AND corona at 0.6 kV")
    for r, lab in ((0.2, "0.4 mm notch tip (2000 of them)"), (2.0, "4 mm rounded bead"),
                   (10.0, "20 mm sphere")):
        print(f"   {lab:32s} corona onset {corona_onset(r)/1000:5.2f} kV")
    print(f"   The drop-maximising geometry is 2000 corona points. It cannot hold kV.")

    print()
    print("6. pH -- the acidity worry is unfounded, and self-limiting")
    q = Is / (P.Q_TRICKLE / 3600 * 1e-3 / 16 / (P.DROP_V * 1e-9))
    n = q / FARADAY
    c = n / (P.DROP_V * 1e-9 * 1e3)
    print(f"   {q*1e12:.0f} pC/drop = {n:.2e} mol H+ in {P.DROP_V:.1f} uL "
          f"= {c:.2e} mol/L excess")
    print(f"   pH 7.00 -> {-math.log10(1e-7+c):.2f}. You cannot strip more H+ than")
    print("   self-ionisation supplies, and substrate buffers what little there is.")

    print()
    print("7. OZONE -- this IS the hazard, not acid")
    for vol in (20, 60):
        mg, ppm = ozone_per_day(6000, 1e-8, vol)
        flag = "SAFE" if ppm < 0.1 else f"{ppm/0.1:.0f}x OVER LIMIT"
        print(f"   6 kV x 10 nA, 24 h -> {mg:.3f} mg O3 in {vol:2d} L = {ppm:5.2f} ppm  {flag}")
    print("   0.1 ppm is the human 8 h limit; invertebrates are no less sensitive.")
    print("   Springtails, isopods, worms and moss share that air. Corona must not")
    print("   happen inside the sealed living volume.")

    print()
    print("8. ENERGY -- harmless to a person")
    E = 0.5 * 5e-12 * 6000 ** 2
    print(f"   0.5.C.V^2 = {E*1e6:.0f} uJ per discharge. Perception threshold ~1 mJ.")
    print("   Startling at worst. The livestock, not the operator, is the concern.")

    print()
    print("9. THE PRIZE: a fourth oscillator, on the same law")
    C = 5e-12
    g = kelvin_growth(f, D, C)
    T = math.log(6000 / 1e-3) / g
    spans = [("spark discharge", 1e-6), ("drop pinch-off", 0.2),
             ("KELVIN CYCLE", T), ("terrace cell", 3.0),
             ("bell siphon", 600.0), ("solar day", 86400.0)]
    for lab, t in spans:
        print(f"   {lab:18s} {t:10.2e} s")
    dec = math.log10(spans[-1][1] / spans[0][1])
    print(f"   accumulate -> threshold -> discharge -> reset, over {dec:.1f} DECADES")
    print("   (was 5.5 with water alone). Charge is a fourth substrate, and it is a")
    print("   working part, not an engraving -- which is what BRIEF.md demands.")
    print("=" * 74)


if __name__ == "__main__":
    report()
