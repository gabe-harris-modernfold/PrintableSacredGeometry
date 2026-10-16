#!/usr/bin/env python3
"""Self-weight sag of the spoke-less desk sculpture (negative_space_desk.py).

With no spokes the ribbon is a coil spring held only where it enters the two pole
caps, and the centre rod is all that holds the caps apart and square. Model: 3-D
Euler-Bernoulli frame -- the ribbon along its mid-thickness centreline, the rod up
the axis, each cap a stiff link from the rod to the ribbon's end -- clamped at the
eye where it hangs, gravity only. Ribbon elements carry the local band section
(width w across, thickness t radial) scaled by the share of plastic any holes
leave there. Short-term PETG modulus; creep under a constant load adds to this
over months. Linear: treat sags past a few percent of the size as indicative.

    python sag_check.py                          # the sizes in negative_space_desk.py
    python sag_check.py --thick 6 --rod 8 5 4 3  # any combination, mm"""

import argparse

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import negative_space_desk as nd
import preview_render as pr

E, NU = 2000.0, 0.38            # PETG, MPa (N/mm^2), short-term
G = E / (2.0 * (1.0 + NU))
WEIGHT_N_MM3 = pr.PETG_G_MM3 * 1e-3 * 9.81
N_EL, N_ROD = 900, 80
CAP_STIFF = 1000.0              # caps: short, thick plates -> effectively rigid links


def frame_k(L, A, Iy, Iz, J, stiff=1.0):
    e, g = E * stiff, G * stiff
    k = np.zeros((12, 12))
    k[0, 0] = k[6, 6] = e * A / L
    k[0, 6] = -e * A / L
    k[3, 3] = k[9, 9] = g * J / L
    k[3, 9] = -g * J / L
    a, b, c, d = 12 * e * Iz / L ** 3, 6 * e * Iz / L ** 2, 4 * e * Iz / L, 2 * e * Iz / L
    k[1, 1] = k[7, 7] = a
    k[1, 7] = -a
    k[1, 5] = k[1, 11] = b
    k[5, 7] = k[7, 11] = -b
    k[5, 5] = k[11, 11] = c
    k[5, 11] = d
    a, b, c, d = 12 * e * Iy / L ** 3, 6 * e * Iy / L ** 2, 4 * e * Iy / L, 2 * e * Iy / L
    k[2, 2] = k[8, 8] = a
    k[2, 8] = -a
    k[2, 4] = k[2, 10] = -b
    k[4, 8] = k[8, 10] = b
    k[4, 4] = k[10, 10] = c
    k[4, 10] = d
    return np.triu(k) + np.triu(k, 1).T


def axes(x, y=None):
    """Element axes as rows: x along it, y the section's width direction."""
    if y is None:
        y = np.cross([0.0, 0.0, 1.0] if abs(x[2]) < 0.9 else [1.0, 0.0, 0.0], x)
    y = y - y.dot(x) * x
    y /= np.linalg.norm(y)
    return np.vstack([x, y, np.cross(x, y)])


def sag(rib, thick, rod_d, place=None, normal=None, top=None, stand=None):
    """place(lat, lon, offset) -> ribbon points `offset` mm out from the surface it
    wraps and normal(points) -> that surface's unit normals; default the sphere.
    top: z of the rod's upper end (default nd.TUBE_TOP). Hangs clamped at the eye,
    or with stand = a z level, stands instead: ribbon nodes below it are pinned."""
    top = nd.TUBE_TOP if top is None else top
    lat_all, _ = pr.lat_lon(rib.tv)
    on = np.nonzero(pr.R * np.cos(lat_all) >= nd.CAP_R)[0]
    tv = np.linspace(rib.tv[on[0]], rib.tv[on[-1]], N_EL + 1)
    lat, lon = pr.lat_lon(tv)
    if place is None:
        rib_xyz = pr.sphere_pt(lat, lon, pr.R + 0.5 * thick)
    else:
        rib_xyz = place(lat, lon, 0.5 * thick)

    tc = 0.5 * (tv[:-1] + tv[1:])
    w = rib.width(tc)
    v = (np.arange(40) + 0.5) / 40 - 0.5
    frac = np.zeros(N_EL)                        # plastic share of each element
    for s in np.linspace(-0.4, 0.4, 5):
        ts = tc + s * (tv[1] - tv[0])
        frac += 1.0 - rib.holes(ts[:, None], v[None, :] * rib.width(ts)[:, None]).mean(axis=1)
    frac /= 5

    # Rod nodes from the bottom cap up to the eye, with a node at the top cap.
    z_rod = np.unique(np.r_[np.linspace(rib_xyz[-1, 2], top, N_ROD), rib_xyz[0, 2]])
    nodes = np.vstack([rib_xyz, np.column_stack([np.zeros((len(z_rod), 2)), z_rod])])
    rod0 = N_EL + 1
    top_cap = rod0 + int(np.searchsorted(z_rod, rib_xyz[0, 2]))

    # (i, j, A, Iy, Iz, J, stiffness factor, carries weight, ribbon element index)
    els = []
    for e in range(N_EL):
        A = frac[e] * w[e] * thick
        els.append((e, e + 1, A, frac[e] * w[e] * thick ** 3 / 12.0,
                    frac[e] * thick * w[e] ** 3 / 12.0,
                    frac[e] * w[e] * thick ** 3 / 3.0 * (1.0 - 0.63 * thick / max(w[e], thick)),
                    1.0, True, e))
    r = 0.5 * rod_d
    rod_I = np.pi * r ** 4 / 4.0
    for k in range(len(z_rod) - 1):
        els.append((rod0 + k, rod0 + k + 1, np.pi * r * r, rod_I, rod_I, 2 * rod_I,
                    1.0, True, None))
    for rod_node, rib_node in [(top_cap, 0), (rod0, N_EL)]:
        els.append((rod_node, rib_node, 100.0, 1000.0, 1000.0, 2000.0, CAP_STIFF, False, None))

    n_dof = 6 * len(nodes)
    rows, cols, vals = [], [], []
    F = np.zeros(n_dof)
    kept = []
    for i, j, A, Iy, Iz, J, stiff, heavy, e in els:
        x = nodes[j] - nodes[i]
        L = np.linalg.norm(x)
        x /= L
        y = None
        if e is not None:                        # ribbon: width runs across, n x tangent
            mid = 0.5 * (nodes[i] + nodes[j])
            n = mid / np.linalg.norm(mid) if normal is None else normal(mid[None])[0]
            y = np.cross(n, x)
        T = np.kron(np.eye(4), axes(x, y))
        k_loc = frame_k(L, A, Iy, Iz, J, stiff)
        dof = np.r_[6 * i:6 * i + 6, 6 * j:6 * j + 6]
        rows.append(np.repeat(dof, 12))
        cols.append(np.tile(dof, 12))
        vals.append((T.T @ k_loc @ T).ravel())
        if heavy:
            F[[6 * i + 2, 6 * j + 2]] -= 0.5 * WEIGHT_N_MM3 * A * L
        kept.append((T, k_loc, dof, A, Iy, Iz, e, L))
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n_dof, n_dof))
    eye = len(nodes) - 1                         # hang point: clamped
    fixed = np.arange(6 * eye, 6 * eye + 6)
    if stand is not None:                        # on a table: pin what rests on it
        down = np.nonzero(rib_xyz[:, 2] < stand)[0]
        fixed = (6 * down[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(n_dof), fixed)
    u = np.zeros(n_dof)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])

    # Peak fibre stress per ribbon / rod element: axial + both bending moments.
    sig_rib, sig_rod = np.zeros(N_EL), 0.0
    for T, k_loc, dof, A, Iy, Iz, e, L in kept:
        f = k_loc @ (T @ u[dof])
        my, mz = max(abs(f[4]), abs(f[10])), max(abs(f[5]), abs(f[11]))
        if e is not None:
            sig_rib[e] = abs(f[0]) / A + my * 0.5 * thick / Iy + mz * 0.5 * w[e] / Iz
        elif A < 50.0 * rod_d:                   # rod (not a cap link)
            sig_rod = max(sig_rod, abs(f[0]) / A + np.hypot(my, mz) * r / rod_I)
    disp = u.reshape(-1, 6)[:, :3]
    i_dn = np.argmin(disp[:N_EL + 1, 2])
    rod_bow = np.hypot(disp[rod0:, 0], disp[rod0:, 1]).max()
    mass = sum(WEIGHT_N_MM3 / 9.81 * 1e3 * A * L for _, _, _, A, _, _, e, L in kept
               if A < 50.0 * max(rod_d, 1.0) or e is not None)
    i_s = int(np.argmax(sig_rib))
    return dict(sag=-disp[i_dn, 2], t_sag=tv[i_dn], rod_bow=rod_bow, sigma=sig_rib[i_s],
                t_sigma=tc[i_s], sig_rod=sig_rod, mass=mass)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thick", type=float, nargs="+", default=[nd.THICK])
    ap.add_argument("--rod", type=float, nargs="+", default=[2 * nd.ROD_R])
    args = ap.parse_args()
    rib = nd.make_ribbon()
    date = lambda t: pr.datetime(pr.YEAR, 1, 1) + pr.timedelta(days=int(t * pr.DAYS))
    for t in args.thick:
        for d in args.rod:
            r = sag(rib, t, d)
            print(f"ribbon {t:.1f} mm, rod {d:.1f} mm: {r['mass']:.0f} g; ribbon drops "
                  f"{r['sag']:.1f} mm at {date(r['t_sag']):%b %d}, rod bows {r['rod_bow']:.1f} "
                  f"mm; peak stress ribbon {r['sigma']:.2f} MPa, rod {r['sig_rod']:.2f} MPa")


if __name__ == "__main__":
    main()
