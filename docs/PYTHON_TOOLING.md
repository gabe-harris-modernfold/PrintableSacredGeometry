# Python tooling available in this project

Verified working on **2026-08-11**. Interpreter: **CPython 3.13.14** (Windows Store build), packages in
the per-user site-packages — no venv, no conda. Run everything as plain `python foo.py` from the repo root.

```bash
python -c "import sys; print(sys.version)"
```

Everything marked **installed** here was **functionally smoke-tested**, not just imported: a real solve, a
real mesh export, a real spectrum, a real conductance plateau. If something marked installed doesn't work,
the environment changed — don't assume the doc is aspirational.

Libraries that are *not* installed are listed too, with their status and the reason — see the
[quantum ecosystem map](#quantum-ecosystem--what-to-reach-for-and-whats-here) and
[Not available — and why](#not-available--and-why). Check those before pip-installing anything: three
PyPI names in this space resolve to packages that are **not** what they appear to be.

Reproduce the install with:

```bash
python -m pip install -r requirements-analysis.txt
```

---

## Core numeric layer

| Library | Version | Use it for |
|---|---|---|
| **NumPy** | 2.4.6 | arrays, linear algebra, FFT. Base layer for everything else. |
| **SciPy** | 1.18.0 | `signal` (filters, Welch PSD, spectrograms, peak finding), `optimize` (fitting, minimization, root finding), `integrate` (`solve_ivp`, quadrature), `special` (Bessel, Legendre, elliptic), `spatial`, `interpolate`. |
| **Matplotlib** | 3.11.1 | beam intensity plots, phase maps, spectra, Bode plots, field slices, sensor logs. |
| **SymPy** | 1.14.0 | symbolic derivation before numerics — ABCD/transfer matrices, dispersion relations, Snell/critical-angle algebra, analytic sanity checks. |
| **pandas** | 3.0.5 | tabular sensor logs, CSV/Excel I/O. |
| **h5py** | 3.16.0 | HDF5 for large datasets. |

**Headless plotting is mandatory** — there is no display in agent runs:

```python
import matplotlib
matplotlib.use("Agg")   # before pyplot import
import matplotlib.pyplot as plt
```

## Quantum / quantum optics

| Library | Version | Use it for |
|---|---|---|
| **QuTiP** | 5.3.1 | Hamiltonians, density matrices, open-system dynamics (`mesolve`), driven cavities, time-dependent drives, Bloch spheres. |

**QuTiP 5 API change:** `c_ops` and `e_ops` are keyword-only now.

```python
res = qt.mesolve(H, psi0, tlist, c_ops=[...], e_ops=[...])   # v5 — correct
res = qt.mesolve(H, psi0, tlist, [...], [...])               # v4 style — TypeError
```

## Mesoscopic quantum transport

| Library | Version | Use it for |
|---|---|---|
| **Kwant** | 1.5.0 | Tight-binding transport on arbitrary geometries: conductance, S-matrix, scattering states, LDOS, band structure, edge states, nanostructures. Bundles **tinyarray** 1.2.5 (its small-array backend, used for spinful/multi-orbital onsite terms). |

Verified end to end: a clean W=10 wire gives conductance quantized in exact integer plateaus 0 → 4.
`kwant.continuum.discretize` (SymPy-backed k·p → tight-binding), `kwant.kpm`, `kwant.operator.Density`
/`Current`, `physics.Bands`, `physics.magnetic_gauge`, and all of `kwant.plotter` (`map`, `density`,
`current`, `bands`, `spectrum`, `plot`) work headless.

```python
import matplotlib; matplotlib.use("Agg")
import kwant, numpy as np
lat = kwant.lattice.square(1.0, norbs=1)
syst = kwant.Builder()
syst[(lat(x, y) for x in range(30) for y in range(10))] = 4.0   # onsite = 4t
syst[lat.neighbors()] = -1.0
lead = kwant.Builder(kwant.TranslationalSymmetry((-1.0, 0)))
lead[(lat(0, y) for y in range(10))] = 4.0
lead[lat.neighbors()] = -1.0
syst.attach_lead(lead); syst.attach_lead(lead.reversed())
f = syst.finalized()
g = kwant.smatrix(f, 0.8).transmission(1, 0)
```

Three things to know:

- **Energy window.** With onsite `4t` and hopping `-t` the band occupies **[0, 8t]**, so probe energies
  like `0.05 … 1.2`. Asking for a negative energy returns transmission `0.0` with no error — that is
  correct physics below the band edge, not a broken solver. Easy to misread as a failure.
- **No MUMPS — SciPy sparse fallback.** Kwant warns `MUMPS is not available … performance can be very
  poor` on every import. Correctness is unaffected; large systems will just be slow. MUMPS needs
  conda/WSL. Suppress with `warnings.filterwarnings("ignore", category=RuntimeWarning)` if it's noise.
- **`magnetic_gauge` needs one field per lead**, not just the bulk: `gauge(bulk_B, lead0_B, lead1_B)`.

**Kwant was built from source here — it is not a wheel off PyPI.** `pip install kwant` **fails** on this
machine: the sdist ships pre-generated Cython C from the NumPy 1.x era, and NumPy 2 moved `subarray` out
of `PyArray_Descr` (`error C2039`). The working recipe, if it ever needs rebuilding:

```bash
python -m pip install "cython>=3.0" wheel
python -m pip download --no-binary :all: --no-deps --no-build-isolation -d . kwant
tar -xzf kwant-1.5.0.tar.gz && cd kwant-1.5.0
rm kwant/_system.c kwant/operator.c
INCLUDE="$(python -c 'import numpy;print(numpy.get_include())');$INCLUDE" python setup.py --cython bdist_wheel
python -m pip install dist/kwant-1.5.0-cp313-cp313-win_amd64.whl
```

Deleting the stale `.c` files and passing `--cython` is the crux — it regenerates them with Cython 3.2,
which emits NumPy-2-correct accessors. Requires the MSVC build tools (VS 2019 BuildTools are present).

### Quantum ecosystem — what to reach for, and what's here

Two of the nine common quantum libraries are installed. The rest are **evaluated but deliberately not
installed**; all were checked against this environment and **none would downgrade NumPy, SciPy, or
QuTiP**, so any of them can be added on request without risking the stack.

| Library | Here? | Best for | Install note |
|---|---|---|---|
| **QuTiP** | ✅ 5.3.1 | Open quantum systems, quantum optics, spin systems, driven/dissipative dynamics. Schrödinger / master-equation / stochastic solvers. | — |
| **Kwant** | ✅ 1.5.0 | Mesoscopic transport, tight-binding devices, conductance, edge states, nanostructures. | Source-built; see recipe above. |
| **Qiskit** | ❌ | Gate-model circuits, transpilation, error mitigation, IBM Quantum hardware execution. | `pip install qiskit qiskit-aer` — clean. |
| **PennyLane** | ❌ | Differentiable/parameterized circuits, hybrid quantum–classical, QML, variational chemistry. | `pip install pennylane` — clean. |
| **Cirq** | ❌ | NISQ circuits where device constraints and circuit structure matter; noisy simulation. | `pip install cirq` — clean. |
| **OpenFermion** | ❌ | Fermionic operators, electronic structure → qubit mappings, lattice models, materials. | `pip install openfermion` — clean. |
| **NetKet** | ❌ | Neural quantum states, variational ground states and dynamics for large many-body systems. | Pulls JAX 0.11. **jaxlib on Windows is CPU-only** — no CUDA/TPU, which is most of NetKet's appeal. Use WSL2 if you need the acceleration. |
| **TeNPy** | ❌ | Tensor networks in Python — approachable DMRG/TEBD for strongly correlated systems. | `pip install physics-tenpy` — clean. The Python answer for MPS/MPO work. |
| **ITensor** | ❌ **n/a** | High-performance tensor networks, MPS/MPO/DMRG/TEBD with index bookkeeping handled for you. | **Not a Python library — it is Julia (and C++).** Nothing to pip install. Requires a Julia toolchain. In Python, reach for **TeNPy** or **quimb** instead. |

`quimb` is the other strong Python tensor-network option (also not installed, also clean) — worth
considering alongside TeNPy if a task needs tensor networks plus general quantum-information tooling.

Rough routing for this repo's interests: **QuTiP** for anything with dissipation, driving, or cavity
coupling; **Kwant** for spatial/geometric tight-binding structure and transport; **SymPy** first for any
closed-form derivation before either. Circuit frameworks (Qiskit/Cirq/PennyLane) are only relevant if the
work turns toward gate-model algorithms, which nothing here currently does.

**Kwant's own test suite: 368 passed, 8 failed, 102 skipped.** The 8 failures are drift in Kwant's *test
helpers* against Matplotlib 3.11 / SciPy 1.18 / SymPy 1.14 (e.g. tests calling `scipy.stats.linregress`
with a packed 2D array), **not** defects in Kwant. I re-checked all 13 affected features directly at
runtime and every one works. Don't chase these.

## Parametric CAD → STL/STEP

| Library | Version | Use it for |
|---|---|---|
| **build123d** | 0.11.1 | Python-first parametric CAD. Preferred for new geometry work in this repo: enclosures, sensor mounts, lattices, jigs, brackets, Platonic/sacred-geometry solids. Direct `export_stl` / `export_step`. |
| **CadQuery** | 2.8.0 | Mature dimension-driven CAD on the same OCCT kernel (`cadquery-ocp` 7.9.3). Use when you want its fluent selector API or high-quality STEP/AMF/3MF export with mesh-tolerance control. |
| **SolidPython2** | 2.1.3 | Generates OpenSCAD `.scad` source from Python (declarative CSG). **Note:** the `openscad` binary is *not* installed, so this can author SCAD but cannot render to STL here. Import name is `solid2`, not `solid`. |
| **ezdxf** | 1.4.4 | DXF read/write for 2D profiles and laser/CNC output. |
| **lib3mf** | 2.5.0 | 3MF container read/write (pulled in by build123d). |

```python
from build123d import BuildPart, Box, Cylinder, Mode, export_stl
with BuildPart() as p:
    Box(20, 20, 10)
    Cylinder(4, 20, mode=Mode.SUBTRACT)
export_stl(p.part, "part.stl")
```

```python
import cadquery as cq
r = cq.Workplane("XY").box(20, 20, 8).edges("|Z").fillet(2)
cq.exporters.export(r, "part.step")
cq.exporters.export(r, "part.stl", tolerance=0.01, angularTolerance=0.1)
```

## Mesh inspection, repair, and visualization

| Library | Version | Use it for |
|---|---|---|
| **trimesh** | 5.0.0 | The workhorse for existing STLs: load, `is_watertight`, `volume`, `bounds`, `moment_inertia`, transforms, booleans, section/slice, convex hull. Handles STL/OBJ/PLY/GLTF/3MF. |
| **numpy-stl** | 4.0.0 | Fast direct vertex/triangle edits on ASCII or binary STL via vectorized NumPy. Import name is `stl`. |
| **PyVista** | 0.48.4 | Interrogate and render meshes, slice planes, scalar fields. Use off-screen mode in agent runs. |
| **pymeshfix** | 0.18.1 | Repair holes, self-intersections, non-manifold edges on closed triangular surfaces. **Coarsens CAD tessellation — use as a last resort, not routinely.** |
| **meshio** | 5.3.5 | Convert between mesh/FEM formats (`.msh`, `.vtu`, `.xdmf`, …). |
| **gmsh** | 4.15.2 | Generate 2D/3D (tet) meshes from OCC geometry — the meshing front end for FEM work. |

**Printability check** — always verify watertightness before calling an STL done:

```python
import trimesh
m = trimesh.load("part.stl")
assert m.is_watertight, "not printable"
print(m.volume, m.bounds, m.extents)   # extents must fit the 320x320x320 bed
```

**pymeshfix wants strict dtypes** — trimesh's tracked arrays are rejected:

```python
v = np.ascontiguousarray(m.vertices, dtype=np.float64)
f = np.ascontiguousarray(m.faces, dtype=np.int32)     # int32, not int64
vc, fc = pymeshfix.clean_from_arrays(v, f)
```

**PyVista off-screen:**

```python
import pyvista as pv
pv.OFF_SCREEN = True
pv.read("part.stl").plot(screenshot="view.png")
```

## Vibration, modal analysis, fatigue

| Library | Version | Use it for |
|---|---|---|
| **SDyPy** | 0.5.1 | Structural-dynamics ecosystem. `sdypy.EMA` for experimental modal analysis — FRF estimation, pole fitting, stabilization diagrams. Bundles `pyFRF` 1.4.0 and `pyuff` 2.5.6 (UFF file I/O). |
| **SDynPy** | 0.23.0 | Sandia's structural-dynamics package: geometry/data objects, test-geometry handling, file I/O, signal-processing utilities. Emits a harmless `numpy.ndarray size changed` ABI warning on import — it works. |
| **FLife** | 2.2.2 | Vibration-fatigue life from measured/simulated PSDs — Dirlik, Tovo-Benasciutti, narrowband, rainflow. |
| **Vibration Toolbox** | 0.6.10 | SDOF/MDOF mass–spring–damper models for validating analytic results. |

```python
import FLife
sd = FLife.SpectralData(input={"data": x, "dt": 1e-4})   # dict input; tuple form is deprecated
life = FLife.Dirlik(sd).get_life(C=1.8e22, k=7.3)
```

**Local patch applied:** `vibration_toolbox/vibesystem.py:9` called `plt.style.use('seaborn-white')`,
which Matplotlib ≥3.6 removed — it crashed on import. Patched to `'seaborn-v0_8-white'`. **If you
reinstall or upgrade `vibration_toolbox`, the patch is lost and the import breaks again**; reapply it.

## FEM / PDE

| Library | Version | Use it for |
|---|---|---|
| **scikit-fem** | 12.0.2 | Pure-Python FEM: Poisson/electrostatics, heat conduction, linear elasticity, custom weak forms in 1D/2D/3D. Pairs with `gmsh` + `meshio` for real geometry. |

> **FEniCSx is NOT installed and cannot be pip-installed on Windows** — `fenics-dolfinx` publishes no
> Windows wheels. Reaching it requires WSL2, Docker, or conda-forge, none of which are set up here.
> `scikit-fem` is the substitute in place: it covers the same weak-form-based electrostatics / heat /
> deformation / coupled-continuum work at smaller scale, in-process, with no MPI. Escalate to
> FEniCSx-in-WSL only if a model genuinely needs distributed-memory scale or DOLFINx-specific features.

```python
from skfem import *
from skfem.helpers import dot, grad
m = MeshTri().refined(4)
b = Basis(m, ElementTriP1())
A = asm(BilinearForm(lambda u, v, w: dot(grad(u), grad(v))), b)
f = asm(LinearForm(lambda v, w: 1.0 * v), b)
x = solve(*condense(A, f, D=m.boundary_nodes()))
```

## Multibody / rigid-body dynamics

| Library | Version | Use it for |
|---|---|---|
| **MuJoCo** | 3.11.0 | Rigid bodies, joints, constraints, springs, dampers, actuators, contacts and impacts. Fast, deterministic, pip-clean. |

> **PyChrono is NOT installed.** The real Project Chrono Python bindings are distributed via conda-forge
> only. The `pychrono` package on PyPI is a 10 kB stub that is **not** Project Chrono — do not install it.
> MuJoCo is the substitute in place and covers the same mechanism-simulation ground (joints, springs,
> bushings-as-soft-constraints, contact). Use conda or WSL if a task specifically requires Chrono's FEA
> or granular-material modules.

## Not available — and why

| Requested | Status | Path if truly needed |
|---|---|---|
| **FEniCSx** | No Windows wheels exist | WSL2 / Docker / conda-forge. Use `scikit-fem` first. |
| **PyChrono** | conda-forge only; PyPI `pychrono` is an unrelated 10 kB stub — **do not install it** | conda or WSL. Use `mujoco` first. |
| **PyOMA / pyoma2** | `pyoma2` has no Python 3.13 release. **Warning:** `pip install PyOMA` fetches an unrelated *bioinformatics* package (Orthologous Matrix API) that also **downgrades NumPy to 1.26 and breaks this environment** | Output-only/operational modal analysis is reachable today via `sdypy.EMA` + `scipy.signal` (SSI-COV / FDD by hand), or install `pyoma2` under a separate Python 3.12 interpreter. |
| **ITensor** | **Not a Python library** — it is Julia (and C++) | See the [quantum ecosystem map](#quantum-ecosystem--what-to-reach-for-and-whats-here) above. |
| **Qiskit, PennyLane, Cirq, OpenFermion, NetKet, TeNPy, quimb** | Not installed — evaluated, none would break the stack | See the [quantum ecosystem map](#quantum-ecosystem--what-to-reach-for-and-whats-here) above for per-library install notes. |

---

## Notes for agents working in this repo

- **Bed limit is 320 × 320 × 320 mm, PETG.** Check `mesh.extents` against it before declaring a part done; prefer stackable/segmented modules over one oversized print.
- **Prefer build123d** for new parametric parts, **trimesh** for anything that reads or measures an existing STL. Existing repo scripts (`flower_of_life_spheres.py`, `enhanced_render.py`, `gem_trace_reference.py`, …) predate this install — check what they already import before adding a dependency.
- **Derive symbolically, then go numeric.** For optics work in `gem-cut-light-paths.html` / `gem_trace_reference.py`, use SymPy to establish the relation (critical angle, ABCD product, dispersion) and assert the closed form, then implement in NumPy and check the two agree.
- Set `matplotlib.use("Agg")` and `pv.OFF_SCREEN = True` in any script an agent runs.
- Console scripts (`isympy.exe`, `jupyter.exe`, …) are installed to a directory **not on PATH**. Invoke via `python -m <module>` instead.
- Jupyter, notebook, PyQt5/PyQt6, VTK and trame came in as transitive dependencies. They work but nothing here needs them; don't build on them casually.
