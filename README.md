# PrintableSacredGeometry

Parametric sacred-geometry solids and optical/physics studies, output as printable STLs and
self-contained HTML visualizers.

## Layout

| Path | Contents |
| --- | --- |
| `projects/<name>/` | One directory per buildable thing: generator scripts, STLs, preview renders. |
| `lib/` | Shared modules imported by the project scripts. |
| `web/` | Source for the self-contained HTML visualizers. Deploy source. |
| `site/` | Netlify publish directory. `index.html` menu, `robots.txt`, previews, and copies of `web/`. |
| `docs/` | Tooling notes and long-form research writeups. |
| `assets/` | Brand art and reference images not tied to a single build. |

## Running a build script

Scripts read and write **bare relative filenames**, so run them from inside their own project
directory — that is where their STLs and preview PNGs land:

```bash
cd projects/tetrahelix && python tetrahelix.py
```

Scripts that need `lib/` put it on `sys.path` themselves, relative to their own location, so this
works without setting `PYTHONPATH`.

## Shared modules (`lib/`)

| Module | Purpose |
| --- | --- |
| `mesh_kit.py` | Mesh assembly and STL writing (`Mesh`, `tube`, `write_stl`). |
| `render_mesh.py` | Offscreen preview rendering used by the `preview_*.py` scripts. |
| `ray_optics.py` | Ray/vector helpers for the optical studies. |
| `validate_mesh.py` | Watertightness and geometry checks. |

## Projects

| Project | What it is |
| --- | --- |
| [crystal_lattices](projects/crystal_lattices) | The 7 crystal systems, 14 Bravais lattices and the three cubic packings, as printable models. |
| [cymatics_dish](projects/cymatics_dish) | Cymatics petri dish for a 12" driver, plus skirt, support ring and pattern studies. |
| [dodecahedron_mirror_frame](projects/dodecahedron_mirror_frame) | 12-piece frame holding twelve 10 cm hexagonal mirror tiles. |
| [e8_torus_lattice](projects/e8_torus_lattice) | The 240 roots of E8 as a printable 150 mm horn-torus strut lattice. |
| [equilateral_frame](projects/equilateral_frame) | Space frame for the (3,7) phyllotactic equilateral vortex tube. |
| [faraday_disc](projects/faraday_disc) | 200 mm tactile Faraday-instability teaching disc. |
| [fcc_negative_space](projects/fcc_negative_space) | Interstitial voids of an FCC close-packing of spheres, and its connected "sponge" variant. |
| [five_intersecting_tetrahedra](projects/five_intersecting_tetrahedra) | 150 mm compound of five woven tetrahedral frames. |
| [flower_of_life](projects/flower_of_life) | Flower of Life built from spheres. |
| [gem_cut_light_paths](projects/gem_cut_light_paths) | Pure-numpy auditor for the gem-cut light-path visualizer in `web/`. |
| [horn_torus](projects/horn_torus) | The torus where tube radius equals revolution radius, plus spindle and vortex variants. |
| [incense-vortex](projects/incense-vortex) | Fully passive helical smoke vortex — no fans, no electronics. |
| [magnetic_vortex](projects/magnetic_vortex) | Parametric magnetic-vortex relief, clockwise and counter-clockwise. |
| [misc_prints](projects/misc_prints) | Standalone STL/3MF prints with no generator script in the repo. |
| [rainbow_tray](projects/rainbow_tray) | Hex-tile prism tray and tabletop scene (OpenSCAD). |
| [sacred_geometry_plates](projects/sacred_geometry_plates) | The numbered 2D plate series (circle → icosahedron) and its OpenSCAD generator. |
| [scrying_pool](projects/scrying_pool) | Golden-ratio scrying pool: flared half-cosine wall, sound-collector rim, stand, print segments. |
| [stellated_platonics](projects/stellated_platonics) | First stellation of each Platonic solid, at 150 mm. |
| [terrarium](projects/terrarium) | Multifunction terrarium — brief and idea notes. |
| [tetrahelix](projects/tetrahelix) | Boerdijk–Coxeter tetrahelix: plain, braided three-ply, ribbon and stellated variants. |

## Web visualizers

Sources live in `web/`. To publish:

```bash
pwsh ./make-previews.ps1   # refresh site/previews/*.webp from web/*.html
pwsh ./deploy.ps1          # copy web/*.html into site/ and push to Netlify
```

New pages must be added to `site/index.html` by hand.

## Constraints

Printer bed is **320 × 320 × 320 mm, PETG**. See [CLAUDE.md](CLAUDE.md) for the full build
constraints and [docs/PYTHON_TOOLING.md](docs/PYTHON_TOOLING.md) for the installed analysis stack.
