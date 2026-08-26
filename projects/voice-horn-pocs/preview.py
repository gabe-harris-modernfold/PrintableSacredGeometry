#!/usr/bin/env python3
"""Preview renders of the POC parts (lib/render_mesh z-buffer rasteriser)."""
import os
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "lib"))
import render_mesh as R


def shot(paths, out, az=32, el=18, w=1300, h=900, section=False):
    layers = []
    for p, mat in paths:
        m = trimesh.load(p)
        v, f = np.asarray(m.vertices), np.asarray(m.faces)
        if section:
            c = v[f].mean(1)
            f = f[c[:, 1] <= 0.5]
        layers.append((v, f, mat))
    R.render(layers, az, el, w=w, h=h, out=out)
    print("wrote", out)


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        sys.exit("usage: preview.py out.png [--section] file.stl[:mat] ...")
    out, args = args[0], args[1:]
    sec = "--section" in args
    args = [a for a in args if a != "--section"]
    shot([(a.split(":")[0], a.split(":")[1] if ":" in a else "frame")
          for a in args], out, section=sec)
