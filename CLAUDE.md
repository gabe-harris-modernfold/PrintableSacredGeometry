# PrintableSacredGeometry

Parametric sacred-geometry solids and optical/physics studies, output as printable STLs and
self-contained HTML visualizers.

## Repository layout

- `projects/<name>/` — one directory per buildable thing: its generator scripts, STLs and preview
  renders live together. **New work gets a new project directory here**, not a file at the root.
- `lib/` — shared modules: `mesh_kit`, `render_mesh`, `ray_optics`, `validate_mesh`.
- `web/` — source for the HTML visualizers. `site/` is the Netlify publish dir; `deploy.ps1` copies
  `web/*.html` into it, so **edit `web/`, never `site/`** (except `site/index.html`, which is the
  menu page and lives there only).
- `docs/` — tooling notes and research writeups. `assets/` — brand art not tied to one build.

Build scripts read and write **bare relative filenames**, so run them from inside their own project
directory (`cd projects/boerdijk_coxeter_tetrahelix && python tetrahelix.py`); that is where their output lands. A
script needing `lib/` adds it to `sys.path` itself, relative to its own file:

```python
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
```

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
