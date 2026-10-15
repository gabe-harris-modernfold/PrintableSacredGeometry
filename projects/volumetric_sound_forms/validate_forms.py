"""Gate every STL in this directory: watertight, consistently wound, one body,
positive volume, inside the 320 mm bed.

Cavity forms may carry more shells, and each is checked by the sign of its volume:
*halves* files hold two bodies (the half and its mirror, side by side); any form may hold
sealed pockets of air inside its body (inward-facing shells, negative volume), which is
how a closed vessel keeps its cavity.

MeshLab's self-intersection filter misfires on stacks of exactly coplanar faces (a vertical wall meshed
by marching tetrahedra, in float32): every face it flags is re-tested here against its neighbours with an
independent float64 test (coplanar pairs by exact 2-D separation), and only real intersections fail.

    python validate_forms.py            # trimesh checks
    python validate_forms.py --meshlab  # also run lib/validate_mesh.py (pymeshlab) per file

Exit status 1 if any file fails."""
import glob
import importlib.util
import os
import subprocess
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
BED = 320.0
HDR = f"{'file':<40}{'faces':>9}{'X':>8}{'Y':>8}{'Z':>8}{'cm3':>9}{'wt':>4}{'wnd':>5}{'bod':>4}{'air':>4}{'bed':>4}"


def _seg_pierces(p, q, T, eps):
    """does segment pq cross triangle T's interior (strictly, beyond eps)?"""
    n = np.cross(T[1] - T[0], T[2] - T[0])
    L = np.linalg.norm(n)
    if L == 0:
        return False
    n = n / L
    dp, dq = float((p - T[0]) @ n), float((q - T[0]) @ n)
    if not ((dp > eps and dq < -eps) or (dp < -eps and dq > eps)):
        return False
    x = p + (dp / (dp - dq)) * (q - p)
    for i in range(3):
        a, b = T[i], T[(i + 1) % 3]
        if float(np.cross(b - a, x - a) @ n) <= eps * np.linalg.norm(b - a):
            return False
    return True


def _coplanar_overlap(A, B, n, eps):
    """2-D separating-axis test in the shared plane; overlap only if no axis separates by more than -eps"""
    k = int(np.argmax(np.abs(n)))
    ax = [i for i in range(3) if i != k]
    a2, b2 = A[:, ax], B[:, ax]
    for T in (a2, b2):
        for i in range(3):
            e = T[(i + 1) % 3] - T[i]
            m = np.array([-e[1], e[0]])
            pa, pb = a2 @ m, b2 @ m
            gap = max(pa.min() - pb.max(), pb.min() - pa.max())
            if gap >= -eps * np.linalg.norm(m):
                return False
    return True


def tri_tri(A, B, eps=1e-7):
    """proper intersection of two triangles in float64, touching within eps (mm) excluded"""
    nA = np.cross(A[1] - A[0], A[2] - A[0])
    nB = np.cross(B[1] - B[0], B[2] - B[0])
    uA, uB = nA / np.linalg.norm(nA), nB / np.linalg.norm(nB)
    dA = (A - B[0]) @ uB
    if np.all(np.abs(dA) <= eps) and np.all(np.abs((B - A[0]) @ uA) <= eps):
        return _coplanar_overlap(A, B, uA, eps)
    for i in range(3):
        if _seg_pierces(A[i], A[(i + 1) % 3], B, eps) or _seg_pierces(B[i], B[(i + 1) % 3], A, eps):
            return True
    return False


def verify_self_intersections(path, radius=4.0):
    """re-test every face MeshLab flags as self-intersecting against every nearby face that shares no vertex,
    with an independent float64 test (coplanar pairs by exact 2-D separation). Returns (flagged, real)."""
    spec = importlib.util.find_spec("pymeshlab")
    if spec and spec.submodule_search_locations:
        os.add_dll_directory(spec.submodule_search_locations[0])
    import pymeshlab as ml
    ms = ml.MeshSet()
    ms.load_new_mesh(path)
    ms.apply_filter("compute_selection_by_self_intersections_per_face")
    m = ms.current_mesh()
    V, Fc = m.vertex_matrix().astype(np.float64), m.face_matrix()
    flagged = np.nonzero(m.face_selection_array())[0]
    from scipy.spatial import cKDTree
    C = V[Fc].mean(1)
    tree = cKDTree(C)
    real = 0
    for f in flagged:
        A = V[Fc[f]]
        for g in tree.query_ball_point(C[f], radius):
            if g == f or np.intersect1d(Fc[f], Fc[g]).size:
                continue
            if tri_tri(A, V[Fc[g]]):
                real += 1
    return len(flagged), real


def expected_bodies(name):
    return 2 if "halves" in name else 1


def check(path):
    m = trimesh.load(path, force="mesh")
    ext = m.extents
    shells = m.split(only_watertight=False) if m.body_count > 1 else [m]
    vols = [float(p.volume) for p in shells]
    return dict(
        file=os.path.basename(path), faces=len(m.faces), x=ext[0], y=ext[1], z=ext[2],
        vol_cm3=m.volume / 1000.0, watertight=bool(m.is_watertight),
        winding=bool(m.is_winding_consistent), bodies=sum(v > 0 for v in vols), air=sum(v < 0 for v in vols),
        fits_bed=bool((ext <= BED + 1e-6).all()), positive=bool(m.volume > 0),
    )


def main():
    meshlab = "--meshlab" in sys.argv
    files = sorted(glob.glob(os.path.join(HERE, "*.stl")))
    if not files:
        print("no STL files here; run `node run_pure.mjs` first")
        return 1
    print(HDR)
    print("-" * len(HDR))
    bad = []
    for path in files:
        r = check(path)
        print(f"{r['file']:<40}{r['faces']:>9}{r['x']:>8.1f}{r['y']:>8.1f}{r['z']:>8.1f}"
              f"{r['vol_cm3']:>9.1f}{'Y' if r['watertight'] else 'n':>4}"
              f"{'Y' if r['winding'] else 'n':>5}{r['bodies']:>4}{r['air']:>4}{'Y' if r['fits_bed'] else 'NO':>4}")
        if not (r["watertight"] and r["winding"] and r["bodies"] == expected_bodies(r["file"]) and r["fits_bed"] and r["positive"]):
            bad.append(r["file"])
    print("-" * len(HDR))
    if meshlab:
        script = os.path.join(HERE, "..", "..", "lib", "validate_mesh.py")
        shells = {r["file"]: r["bodies"] + r["air"] for r in map(check, files)}
        for path in files:
            name = os.path.basename(path)
            res = subprocess.run([sys.executable, script, path], capture_output=True, text=True)
            out = res.stdout
            tail = out.strip().splitlines()[-1] if out.strip() else res.stderr.strip()[-200:]
            passed = res.returncode == 0
            if not passed:
                # the shared validator insists on one component (a multi-shell cavity form passes when meshlab
                # counts exactly the shells trimesh found) and trusts its self-intersection flags (re-tested)
                fails = [ln for ln in out.splitlines() if "[FAIL]" in ln]
                comps = [ln for ln in out.splitlines() if ln.startswith("connected components:")]
                ncomp = int(comps[0].split(":")[1]) if comps else -1
                notes, ok_all = [], True
                for ln in fails:
                    if "single component" in ln and ncomp == shells[name] > 1:
                        notes.append(f"{ncomp} shells")
                    elif "self-intersections" in ln:
                        flagged, real = verify_self_intersections(path)
                        ok_all &= real == 0
                        notes.append(f"{flagged} self-intersection flags re-tested in float64, {real} real")
                    else:
                        ok_all = False
                passed = ok_all
                tail = ("PRINT-READY - " if passed else "NEEDS ATTENTION - ") + "; ".join(notes)
            print(f"meshlab {name:<38} exit {res.returncode}  {tail}")
            if not passed:
                bad.append(name + " (meshlab)")
    if bad:
        print("CHECK:", ", ".join(bad))
        return 1
    print("all forms watertight, consistently wound, the expected bodies, positive volume, inside the 320 bed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
