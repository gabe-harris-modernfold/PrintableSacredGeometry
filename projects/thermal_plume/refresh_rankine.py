"""Re-render the Rankine/Oseen figure from cached fits.

    python refresh_rankine.py --run run_ra1e7

Fitting several hundred cores is the expensive part of figures.py, and it does not
change when the figure's layout or wording does. figures.py writes core_fits.json;
this reads it back and re-renders in a second.
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figures as F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run_ra1e7")
    args = ap.parse_args()

    path = os.path.join(args.run, "core_fits.json")
    if not os.path.exists(path):
        print(f"{path} not found — run figures.py once to produce it")
        return
    with open(path) as fh:
        fits = json.load(fh)
    fits = [f for f in fits if f.get("rankine") and f.get("oseen")]
    print(f"{len(fits)} cached core fits")

    stats = F._rankine_stats(fits)
    print(json.dumps(stats, indent=2))
    F.fig_rankine(fits, os.path.join(args.run, "fig_rankine.png"))
    print("wrote", os.path.join(args.run, "fig_rankine.png"))


if __name__ == "__main__":
    main()
