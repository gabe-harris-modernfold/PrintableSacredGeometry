#!/usr/bin/env python3
"""Preview renders of the assembled pyramid and of the rotor alone.
Run after pyramid.py, from this directory:  python preview.py"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import render_mesh
from render_mesh import render
import pyramid as P

render_mesh.MAT["candle"] = {'col': (236, 226, 200), 'amb': 0.5, 'spec': 0.3, 'shin': 20.0}
render_mesh.MAT["steel"] = {'col': (200, 204, 210), 'amb': 0.3, 'spec': 0.9, 'shin': 60.0}

def layers_of(world):
    out = []
    for name, m in world:
        v, f = np.asarray(m.vertices, float), np.asarray(m.faces, np.int64)
        if name.startswith("blade"):
            mat = "pad"
        elif name.startswith("tealight"):
            mat = "candle"
        elif name == "steel_pin":
            mat = "steel"
        elif name in ("rotor_hub", "apex_collar", "pin_carrier"):
            mat = "glass"
        else:
            mat = "frame"
        out.append((v, f, mat))
    return out

if __name__ == "__main__":
    parts = {"base_plate": P.base_plate(), "pin_carrier": P.pin_carrier(),
             "rotor_hub": P.rotor_hub(), "blade": P.blade()}
    if P.FRAME == "posts":
        parts["apex_collar"] = P.apex_collar()
    else:
        parts["mast"] = P.mast()
    _, post_len = P.post()
    world, z_hub = P.assembled(parts, post_len)
    render(layers_of(world), az=35, el=18, w=1000, h=1200, out="preview_assembly.png")
    render(layers_of(world), az=0, el=89, w=1000, h=1000, out="preview_top.png")
    rotor = [(n, m) for n, m in world if n.startswith(("rotor", "blade", "steel", "pin_carrier", "apex", "mast"))]
    render(layers_of(rotor), az=30, el=25, w=1000, h=800, out="preview_rotor.png")
    render(layers_of(rotor), az=0, el=2, w=1200, h=500, out="preview_rotor_side.png")
    bl = P.blade_print_orientation(parts["blade"])
    render(layers_of([("blade", bl)]), az=-90, el=55, w=1000, h=600, out="preview_blade.png")
