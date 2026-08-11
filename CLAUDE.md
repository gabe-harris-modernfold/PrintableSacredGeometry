# PrintableSacredGeometry

Parametric sacred-geometry solids and optical/physics studies, output as printable STLs and
self-contained HTML visualizers.

## Python tooling

A full numeric / CAD / simulation stack is **already installed and verified** — NumPy, SciPy,
Matplotlib, SymPy, QuTiP, Kwant, build123d, CadQuery, trimesh, numpy-stl, PyVista, pymeshfix, gmsh,
meshio, SDyPy, SDynPy, FLife, Vibration Toolbox, scikit-fem, MuJoCo.

**Read [docs/PYTHON_TOOLING.md](docs/PYTHON_TOOLING.md) before adding a dependency or writing an
analysis/CAD script.** It lists exact versions, which library to reach for per task, the API gotchas
that will otherwise bite (QuTiP 5 keyword-only `c_ops`, pymeshfix int32 requirement, Kwant's band
window and source-build recipe, headless Matplotlib/PyVista), and two PyPI packages that must never be
installed here (`pychrono`, `PyOMA` — both are unrelated to the libraries their names suggest, and
`PyOMA` breaks the NumPy install).

Pinned set: [requirements-analysis.txt](requirements-analysis.txt). Interpreter is CPython 3.13
(Windows Store build), packages in per-user site-packages — no venv, no conda.

## Hardware constraint

Printer bed is **320 × 320 × 320 mm, PETG**. Size printable parts to fit and prefer stackable or
segmented modules over one oversized print. Verify with `trimesh.load(path).extents` before calling a
part done, and check `is_watertight`.
