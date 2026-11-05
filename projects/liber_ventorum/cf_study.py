#!/usr/bin/env python3
"""Carbon-fibre cores in the Plumb Octahedron's tubes: does the simulation want them?

A pultruded carbon rod inside each clear PETG tube carries the compression, so the tip positions follow
carbon (about 0-1 um/m/K) instead of PETG (55-75 um/m/K). The PETG tip nodes still sit between rod end and
string saddle, so the effective expansion of the tip radius is about 7-11 um/m/K, and creep falls ~10x.

    python cf_study.py            # PETG vs carbon for rev 1 + springs and the v2 design, then re-tune
                                  # the six continuous settings for carbon (gauges held at v2's)
"""
import json
import math
import time

import numpy as np
from scipy.optimize import differential_evolution

import evolve_tuning as e

KEYS = ("J", "e1_cents", "e2_cents", "e3_cents", "coherence", "theta_cover_sounding", "gust_cents_10ms",
        "total_mass_kg", "max_stress", "max_spring_ext_mm", "adjuster_travel95_mm")


def carbon(tol, seed):
    rng = np.random.default_rng(seed)
    out = dict(tol)
    n = len(tol["st"])
    out["cte"] = rng.uniform(7e-6, 11e-6, n)
    out["creep"] = tol["creep"] * 0.1
    return out


def summary(info):
    return {k: round(info[k], 3) for k in KEYS}


def main():
    np.seterr(all="ignore")
    tol, fresh = e.tolerances(e.N_MC), e.tolerances(48, seed=99)
    tol_cf, fresh_cf = carbon(tol, 5), carbon(fresh, 6)
    v2 = np.array(json.load(open("evolution_result_v2.json"))["x"])
    rev1b = e.REV1.copy()
    rev1b[1] = math.log10(2.0)
    rows = {}
    for name, x in (("rev1b", rev1b), ("v2", v2)):
        for frame, t in (("petg", fresh), ("carbon", fresh_cf)):
            _, info = e.evaluate(x, t, n=48, detail=True)
            rows[f"{name}_{frame}"] = summary(info)
            print(name, frame, rows[f"{name}_{frame}"], flush=True)

    gauges = v2[6:]
    history = []

    def f(c):
        try:
            return e.evaluate(np.r_[c, gauges], tol_cf)
        except (FloatingPointError, ValueError, np.linalg.LinAlgError):
            return 1e4

    t0 = time.time()
    res = differential_evolution(f, e.BOUNDS[:6], popsize=10, maxiter=30, tol=1e-4, seed=3, polish=False,
                                 updating="deferred", workers=1, init="sobol", x0=v2[:6],
                                 callback=lambda xk, convergence=None: history.append(float(f(xk))))
    xcf = np.r_[res.x, gauges]
    _, info_cf = e.evaluate(xcf, fresh_cf, n=48, detail=True)
    rows["v3_carbon"] = summary(info_cf)
    print("v3 carbon", rows["v3_carbon"], f"{time.time() - t0:.0f} s", flush=True)
    json.dump(dict(rows=rows, x_v3=xcf.tolist(), design_v3=info_cf["design"], history=history,
                   theta_by_u=info_cf["theta_by_u"], beat_vs_temp_3ms=info_cf["beat_vs_temp_3ms"]),
              open("cf_study.json", "w"), indent=1)


if __name__ == "__main__":
    main()
