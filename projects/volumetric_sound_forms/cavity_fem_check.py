"""The page's slowly-varying-duct ("adiabatic") mode solver against a finite-element solution.

The cavity form's modes come from −(a² Z′)′ + κ² Z = k² a² Z on each cross-section mode J_m(κρ)cos(mθ):
exact for a straight cylinder, an approximation for a cone frustum. Here the same frustum is solved
without that approximation: for each azimuthal order m, the axisymmetric Helmholtz eigenproblem in (r, z)

    ∫ (∇u·∇v + m²/r² u v) r dA = k² ∫ u v r dA

with quadratic triangles (scikit-fem), rigid wall and floor (natural), p = 0 on the axis for m ≥ 1, and at
the top either rigid (closed) or p = 0 (the open mouth, as the page models it). Eigenvalues are compared
mode by mode, k·a0 up to --kmax.

    python cavity_fem_check.py                    # the page's default frustum, open and closed
    python cavity_fem_check.py --cr1 30 --kmax 8  # a stronger taper
Exit status 1 if any mode differs by more than --tol percent (default 2)."""
import argparse
import sys

import numpy as np
from scipy.sparse.linalg import eigsh
from skfem import Basis, BilinearForm, ElementTriP2, MeshTri
from skfem.helpers import dot, grad

import sound_form_reference as R


def adiabatic(P, m, kmax):
    M = R.cavity_modes(P)
    ks = [md["k"] for md in M["modes"] if md["m"] == m and md["k"] <= kmax]
    return np.sort(np.array(ks))


def fem(P, m, kmax, nr=40):
    G = R.cavity_geom(P)
    a0, Hn = G["a0"], G["Hn"]
    b0, b1 = G["cr0"] / a0, G["cr1"] / a0
    nz = int(np.ceil(nr * Hn / max(b0, b1)))
    mesh = MeshTri.init_tensor(np.linspace(0, 1, nr + 1), np.linspace(0, Hn, nz + 1))
    p = mesh.p.copy()
    p[0] = p[0] * (b0 + (b1 - b0) * p[1] / Hn)
    mesh = MeshTri(p, mesh.t)
    basis = Basis(mesh, ElementTriP2())

    @BilinearForm
    def stiff(u, v, w):
        r = w.x[0]
        return dot(grad(u), grad(v)) * r + (m * m / np.maximum(r, 1e-12)) * u * v

    @BilinearForm
    def mass(u, v, w):
        return u * v * w.x[0]

    K = stiff.assemble(basis).tocsc()
    Mm = mass.assemble(basis).tocsc()
    fixed = []
    if m > 0:
        fixed.append(basis.get_dofs(lambda x: x[0] < 1e-9).all())
    if G["open"]:
        fixed.append(basis.get_dofs(lambda x: x[1] > Hn - 1e-9).all())
    if fixed:
        D = np.unique(np.concatenate(fixed))
        keep = np.setdiff1d(np.arange(K.shape[0]), D)
        K, Mm = K[keep][:, keep], Mm[keep][:, keep]
    vals = eigsh(K, k=60, M=Mm, sigma=0.05, which="LM", return_eigenvectors=False)
    k = np.sqrt(np.clip(np.sort(vals), 0, None))
    return k[(k > 1e-3) & (k <= kmax)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cr0", type=float, default=50)
    ap.add_argument("--cr1", type=float, default=40)
    ap.add_argument("--ch", type=float, default=120)
    ap.add_argument("--kmax", type=float, default=9.0)
    ap.add_argument("--tol", type=float, default=2.0, help="percent")
    a = ap.parse_args()
    worst = 0.0
    for top in ("open", "closed"):
        P = dict(R.defaults(), cr0=a.cr0, cr1=a.cr1, ch=a.ch, ctop=top)
        G = R.cavity_geom(P)
        print(f"\nfrustum r {a.cr0:g} -> {a.cr1:g} mm, height {a.ch:g} mm, {top} top; k*a0 up to {a.kmax:g} (a0 = {G['a0']:g} mm)")
        for m in (0, 1, 2, 3):
            ka, kf = adiabatic(P, m, a.kmax), fem(P, m, a.kmax + 0.5)
            kf = kf[kf <= a.kmax + 0.5]
            # pair each adiabatic mode with the nearest finite-element mode
            err = np.array([(k - kf[np.argmin(np.abs(kf - k))]) / kf[np.argmin(np.abs(kf - k))] * 100 for k in ka])
            w = float(np.abs(err).max()) if len(err) else 0.0
            worst = max(worst, w)
            print(f"  m={m}: {len(ka):2d} modes (FEM {int((kf <= a.kmax).sum()):2d})  max |err| {w:5.2f} %  mean {np.abs(err).mean() if len(err) else 0:5.2f} %   "
                  + " ".join(f"{e:+.2f}" for e in err[:8]))
    ok = worst <= a.tol
    print(f"\nworst {worst:.2f} % against a {a.tol:g} % limit:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
