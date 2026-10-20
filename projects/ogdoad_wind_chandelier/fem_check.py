#!/usr/bin/env python3
"""Independent check of a voice's pitch: 3D acoustic FEM of the real air domain.

The air is the gourd's cavity + head plenum, joined through the window to a sphere of
outside air centred on the window. Walls are rigid (natural Neumann); the sphere's surface
is pressure-release (p = 0). The lowest eigenvalue of  K p = k^2 M p  is the Helmholtz
mode with the outside radiation mass truncated at radius Rs. That truncation removes an
inertance ~ rho / (2 pi Rs), so 1/f^2 is linear in 1/Rs; two radii extrapolate to free air.

The windway is filled (a long 0.7 mm duct in parallel carries next to no flow at the
resonance) and the glue gap of the socket is closed, as glue does. Nothing else is
simplified: plenum, labium underside, neck bore, waist and the 45-deg cones are all as built.

Run:  python fem_check.py [voice names...]     (default: Earth Sun Ogdoad)
"""

import math
import os
import sys
import tempfile
import time

import gmsh
import meshio
import numpy as np
from scipy.sparse.linalg import eigsh
from build123d import Pos, Sphere, export_step

import voices as vc
import chandelier as ch

C_MM = vc.C_AIR * 1000.0
R_OVERRIDE = {}
ALG3D = 10               # gmsh 3-D algorithm: 10 = HXT; 1 = Delaunay (fallback)
H_BORE_K = 0.25
H_PLEN = 0.4             # mesh size in the plenum + window column (mm)          # mesh size in the neck bore, as a fraction of its radius


def air_solid(v, Rs):
    d = v["dims"]
    hs = ch.head_geometry(v)
    zt = d["z_top"]
    xh, ys = hs["xh"], hs["ys"]
    R, Rn = vc.R_LIP, vc.R_LIP + vc.T_BLK
    # the head's outer envelope (no windway, socket filled), minus its plenum
    env = ch.box(-xh, xh, -ys, ys, -vc.SPIG_H, hs["z_f"])
    prof = ([(hs["y_front"], 0.0), (hs["y_rear"], 0.0), (hs["y_rear"], hs["z_ramp"]),
             (hs["yL"], hs["z_tip"]), (hs["yL"], hs["z_f"]), (hs["y0"], hs["z_f"]),
             (hs["y0"], hs["z_bt"])]
            + ch.arc_pts(hs["y_ww0"], hs["z_r"] + R, Rn, 90, 180, 10))
    env += ch.yz_prism(prof, -xh, xh)
    env -= ch.yz_prism(vc.plenum_poly(hs), -v["b"] / 2, v["b"] / 2)
    env = Pos(0, 0, zt) * env
    outer = ch.rev_rz([(0.0, 0.0)] + [tuple(p) for p in v["outer"]] + [(0.0, zt)])
    centre = (0.0, 0.0, zt + hs["z_f"])
    sph = Pos(*centre) * Sphere(Rs)
    air = (sph - env - outer) + ch.cavity(v, hs)
    return air, centre, hs


def mesh_air(v, Rs, h_win, h_mid, h_far):
    air, centre, hs = air_solid(v, Rs)
    tmp = tempfile.mkdtemp()
    step = os.path.join(tmp, "air.step")
    export_step(air, step)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    zt = v["dims"]["z_top"]
    zc = centre[2]
    rb = v["dims"]["rb"]
    # size grows linearly with distance from the window: fine where the flow squeezes
    # through the window, coarse where the pressure is uniform (deep cavity, far air)
    f1 = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        f1, "F", f"{h_win} + {h_mid}*Sqrt(x*x + y*y + (z-({zc}))*(z-({zc})))")
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_far)
    # the neck bore and the taper below it carry the plug's extra mass: resolve them too
    f2 = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(f2, "VIn", H_BORE_K * rb)
    gmsh.model.mesh.field.setNumber(f2, "VOut", 1e3)
    for k, val in (("XMin", -rb - 2), ("XMax", rb + 2), ("YMin", -rb - 2), ("YMax", rb + 2),
                   ("ZMin", zt - 14 - 0.6 * rb), ("ZMax", zt + 0.5)):
        gmsh.model.mesh.field.setNumber(f2, k, val)
    # the 2 mm plenum under the window and the window column itself
    f3 = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(f3, "VIn", H_PLEN)
    gmsh.model.mesh.field.setNumber(f3, "VOut", 1e3)
    for k, val in (("XMin", -v["b"] / 2 - 0.5), ("XMax", v["b"] / 2 + 0.5), ("YMin", hs["y0"] - 3),
                   ("YMax", hs["yL"] + 8), ("ZMin", zt - 0.5), ("ZMax", zt + hs["z_f"] + 1.5)):
        gmsh.model.mesh.field.setNumber(f3, k, val)
    fm = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fm, "FieldsList", [f1, f2, f3])
    gmsh.model.mesh.field.setAsBackgroundMesh(fm)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", ALG3D)
    gmsh.model.mesh.generate(3)
    msh = os.path.join(tmp, "air.msh")
    gmsh.write(msh)
    gmsh.finalize()
    m = meshio.read(msh)
    tets = np.vstack([c.data for c in m.cells if c.type == "tetra"])
    return m.points, tets, np.array(centre)


def gourd_volume_ratio(points, tets, v):
    """Meshed / exact volume of the gourd's interior. Straight-sided tets cut the corners of
    the curved walls, so a coarse mesh loses a few % of the cavity and reads sharp; the
    pitch is corrected by sqrt of this ratio (checked: two mesh densities then agree to
    within 4 cents, where uncorrected they differ by ~9)."""
    zt, inn, t = v["dims"]["z_top"], v["inner"], v["dims"]["t"]
    P = points[tets]
    vol = np.abs(np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                           P[:, 3] - P[:, 0])) / 6
    Cc = P.mean(1)
    rc = np.hypot(Cc[:, 0], Cc[:, 1])
    inside = ((Cc[:, 2] < zt) & (Cc[:, 2] > inn[0, 1])
              & (rc < np.interp(Cc[:, 2], inn[:, 1], inn[:, 0]) + t / 2))
    return vol[inside].sum() / vc.revolved_volume(inn)


def solve(points, tets, centre, Rs, order=2):
    from skfem import Basis, BilinearForm, ElementTetP1, ElementTetP2, MeshTet, asm, condense
    from skfem.helpers import dot, grad
    used = np.unique(tets)
    remap = -np.ones(len(points), int)
    remap[used] = np.arange(len(used))
    mesh = MeshTet(points[used].T, remap[tets].T)
    el = ElementTetP2() if order == 2 else ElementTetP1()
    basis = Basis(mesh, el)
    K = asm(BilinearForm(lambda u, w, _: dot(grad(u), grad(w))), basis)
    M = asm(BilinearForm(lambda u, w, _: u * w), basis)
    # pressure-release on the outer sphere: boundary facets whose midpoints lie on it
    # (every vertex exactly on the sphere AND the facet facing radially: a test on the
    # midpoint's distance alone also catches cavity-wall facets that happen to lie Rs away)
    bf = mesh.boundary_facets()
    V = mesh.p[:, mesh.facets[:, bf]] - centre[:, None, None]       # 3 x 3 x nf
    dist = np.linalg.norm(V, axis=0)
    nrm = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0], axis=0)
    rad = V.mean(axis=1)
    cosang = np.abs((nrm * rad).sum(0)) / (np.linalg.norm(nrm, axis=0) * np.linalg.norm(rad, axis=0))
    on = (np.abs(dist - Rs).max(axis=0) < 1e-3 * Rs) & (cosang > 0.95)
    D = basis.get_dofs(bf[on]).all()
    Kc, Mc, _, I = condense(K, M, D=D)
    k_guess = 2 * math.pi * 50 / C_MM
    w, _ = eigsh(Kc, k=2, M=Mc, sigma=k_guess ** 2 * 0.1, which="LM")
    k = math.sqrt(sorted(w)[0])
    return k * C_MM / (2 * math.pi), mesh.t.shape[1], len(I)


def check(name, radii=(36.0, 60.0), h=(0.18, 0.16, 8.0), order=2):
    f = dict((n, ff) for ff, n, _ in vc.VOICES)[name]
    v = vc.solve_voice(f)
    v.update(name=name)
    if R_OVERRIDE.get(name):                       # check a given gourd size instead
        v = vc.voice_acoustics(R_OVERRIDE[name], f)
        v.update(name=name)
    h_far = h[2] or float(np.clip(0.12 * v["R"], 3.0, 9.0))
    res = []
    for Rs in radii:
        t0 = time.time()
        pts, tets, centre = mesh_air(v, Rs, h[0], h[1], h_far)
        fr, ne, nd = solve(pts, tets, centre, Rs, order)
        fr *= math.sqrt(gourd_volume_ratio(pts, tets, v))
        res.append((Rs, fr))
        print(f"  {name:8s} Rs={Rs:5.1f}  f={fr:8.2f} Hz  ({ne} tets, {nd} dofs, {time.time() - t0:.0f}s)")
    # extrapolate 1/f^2 linearly in 1/Rs to free air
    x = np.array([1 / r for r, _ in res])
    y = np.array([1 / fr ** 2 for _, fr in res])
    a, b0 = np.polyfit(x, y, 1)
    f_inf = 1 / math.sqrt(b0)
    cents = 1200 * math.log2(f_inf / f)
    cm = 1200 * math.log2(f_inf / v["f"])
    print(f"  {name:8s} R={v['R']:.2f}  target {f} Hz | model {v['f']:.1f} | FEM free-air {f_inf:.1f} Hz "
          f"({cents:+.0f} c vs target, {cm:+.0f} c vs model)", flush=True)
    return f, v["f"], f_inf


if __name__ == "__main__":
    names = sys.argv[1:] or ["Ogdoad", "Sun", "Earth"]
    for n in names:
        check(n)
