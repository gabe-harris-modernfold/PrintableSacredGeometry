"""Re-render only the shedding figure and re-report its cadence.

    python refresh_shedding.py --run run_ra1e7

The full figures.py pass spends most of its time fitting vortex cores, which does
not depend on the spectral analysis at all. This regenerates just that one figure
after a change to the frequency estimator, and prints the numbers the writeup quotes.
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import body_model as bm
import figures as F
import vortex_analysis as va


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run_ra1e7")
    args = ap.parse_args()

    d, meta = F.load(args.run)
    res = F.fig_shedding(d, meta, os.path.join(args.run, "fig_shedding.png"))
    if not res:
        print("no usable probe record")
        return

    st = va.strouhal(res["freq_hz"], bm.BODY_HEIGHT, meta["u_ff"])
    st_meas = va.strouhal(1 / 5.0, bm.BODY_HEIGHT, meta["u_ff"])
    out = dict(
        window_s=[round(res["t_start"], 1), round(res["t_end"], 1)],
        shedding_period_s=round(res["period_s"], 2),
        shedding_period_peakcount_s=round(float(res["period_peakcount_s"]), 2),
        slow_meander_period_s=round(float(res["period_slow_s"]), 1),
        n_peaks=res["n_peaks"],
        strouhal_sim=round(st, 3),
        strouhal_measured=round(st_meas, 3),
        ratio_sim_over_measured=round(st / st_meas, 2),
        buoyant_time_s=round(meta["t_buoy"], 2),
        shed_intervals_per_buoyant_time=round(res["period_s"] / meta["t_buoy"], 2),
    )
    print(json.dumps(out, indent=2))

    path = os.path.join(args.run, "shedding.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote", path)


if __name__ == "__main__":
    main()
