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
| `assets/` | Brand art and reference images not tied to a single build, indexed in [assets/README.md](assets/README.md). |

## Running a build script

Scripts read and write **bare relative filenames**, so run them from inside their own project
directory — that is where their STLs and preview PNGs land:

```bash
cd projects/boerdijk_coxeter_tetrahelix && python tetrahelix.py
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
| [boerdijk_coxeter_tetrahelix](projects/boerdijk_coxeter_tetrahelix) | Boerdijk–Coxeter tetrahelix: plain, braided three-ply, ribbon and stellated variants. |
| [bravais_crystal_lattices](projects/bravais_crystal_lattices) | The 7 crystal systems, 14 Bravais lattices and the three cubic packings, as printable models. |
| [cymatics_petri_dish](projects/cymatics_petri_dish) | Cymatics petri dish for a 12" driver, plus skirt, support ring and pattern studies. |
| [dodecahedron_mirror_frame](projects/dodecahedron_mirror_frame) | 12-piece frame holding twelve 10 cm hexagonal mirror tiles. |
| [e8_torus_lattice](projects/e8_torus_lattice) | The 240 roots of E8 as a printable 150 mm horn-torus strut lattice. |
| [equilateral_vortex_tube_frame](projects/equilateral_vortex_tube_frame) | Space frame for the (3,7) phyllotactic equilateral vortex tube. |
| [faraday_instability_disc](projects/faraday_instability_disc) | 200 mm tactile Faraday-instability teaching disc. |
| [fcc_negative_space](projects/fcc_negative_space) | Interstitial voids of an FCC close-packing of spheres, and its connected "sponge" variant. |
| [five_intersecting_tetrahedra](projects/five_intersecting_tetrahedra) | 150 mm compound of five woven tetrahedral frames. |
| [flesh-sheath-standing-wave](projects/flesh-sheath-standing-wave) | Research note on universality classes and renormalization-group flow. No build yet. |
| [flower_of_life_spheres](projects/flower_of_life_spheres) | Flower of Life built from spheres. |
| [gem_cut_light_paths](projects/gem_cut_light_paths) | Pure-numpy auditor for the gem-cut light-path visualizer in `web/`. |
| [golden_ratio_scrying_pool](projects/golden_ratio_scrying_pool) | Golden-ratio scrying pool: flared half-cosine wall, sound-collector rim, stand, print segments. |
| [horn_torus_variants](projects/horn_torus_variants) | The torus where tube radius equals revolution radius, plus spindle and vortex variants. |
| [incense_vortex](projects/incense_vortex) | Fully passive helical smoke vortex — no fans, no electronics. |
| [incense_vortex_ring](projects/incense_vortex_ring) | Hand-pumped smoke-ring generator for the same 100 mm tube, with seven interchangeable orifice shapes. |
| [indirect-perception](projects/indirect-perception) | The eye's spectral window: photopigment fits to CVRL data for `web/indirect-perception.html`. |
| [magnetic_vortex_relief](projects/magnetic_vortex_relief) | Parametric magnetic-vortex relief, clockwise and counter-clockwise. |
| [misc_standalone_prints](projects/misc_standalone_prints) | Standalone STL/3MF prints with no generator script in the repo. |
| [multifunction_terrarium](projects/multifunction_terrarium) | Multifunction terrarium — brief and idea notes. |
| [negative-space](projects/negative-space) | Data ribbons wound around invisible solids (sphere, cube, icosahedron, pyramid), the data encoded as holes. Also builds the assets for `web/indirect-perception.html`. |
| [ogdoad_wind_chandelier](projects/ogdoad_wind_chandelier) | The Hymn of the Ogdoad: nine wind-blown pigeon-whistle voices hung as a geocentric cosmos. |
| [rainbow_prism_tray](projects/rainbow_prism_tray) | Hex-tile prism tray and tabletop scene (OpenSCAD). |
| [sacred_correspondant_forms](projects/sacred_correspondant_forms) | Research report on printed acoustic geometries, cymatics and water resonance. No build yet. |
| [sacred_geometry_plates](projects/sacred_geometry_plates) | The numbered 2D plate series (circle → icosahedron) and its OpenSCAD generator. |
| [stellated_platonic_solids](projects/stellated_platonic_solids) | First stellation of each Platonic solid, at 150 mm. |
| [SubConsiousness](projects/SubConsiousness) | Research note on evidence for influencing nonconscious processing. No build yet. |
| [thermal_plume](projects/thermal_plume) | Thermal lattice-Boltzmann study of the convective plume off a standing person, for `web/human-thermal-plume.html`. |
| [voice-horn-evolution](projects/voice-horn-evolution) | NSGA-II and FEM search that evolved the voice-horn genomes ("The Forge" in `docs/voice-horn-pocs.md`). |
| [voice-horn-pocs](projects/voice-horn-pocs) | Printable parts for the three acoustic-metamaterial voice horns, built from the evolved genomes. |
| [volumetric_sound_forms](projects/volumetric_sound_forms) | Browser tool: mp3/wav → printable STL, as a spectral cylinder or as the cavity form (the clip's standing-wave pressure in a rigid vessel, lobes fused to the shell), plus the node/numpy/FEM harness that verifies the page. |
| [weihnachtspyramide](projects/weihnachtspyramide) | Candle-driven Christmas pyramid: four tea lights turn an eight-blade wheel. |

## Project timeline

First commit of each project's files (renames followed), oldest first:

| Date | Projects added |
| --- | --- |
| 2026-06-17 | [horn_torus_variants](projects/horn_torus_variants), [magnetic_vortex_relief](projects/magnetic_vortex_relief) |
| 2026-07-27 | [incense_vortex](projects/incense_vortex) |
| 2026-07-29 | [fcc_negative_space](projects/fcc_negative_space), [flower_of_life_spheres](projects/flower_of_life_spheres), [misc_standalone_prints](projects/misc_standalone_prints), [rainbow_prism_tray](projects/rainbow_prism_tray), [sacred_geometry_plates](projects/sacred_geometry_plates), [golden_ratio_scrying_pool](projects/golden_ratio_scrying_pool) |
| 2026-08-04 | [dodecahedron_mirror_frame](projects/dodecahedron_mirror_frame) |
| 2026-08-06 | [gem_cut_light_paths](projects/gem_cut_light_paths) |
| 2026-08-11 | [equilateral_vortex_tube_frame](projects/equilateral_vortex_tube_frame), [five_intersecting_tetrahedra](projects/five_intersecting_tetrahedra), [boerdijk_coxeter_tetrahelix](projects/boerdijk_coxeter_tetrahelix) |
| 2026-08-12 | [cymatics_petri_dish](projects/cymatics_petri_dish) |
| 2026-08-13 | [bravais_crystal_lattices](projects/bravais_crystal_lattices), [faraday_instability_disc](projects/faraday_instability_disc) |
| 2026-08-18 | [e8_torus_lattice](projects/e8_torus_lattice), [multifunction_terrarium](projects/multifunction_terrarium) |
| 2026-08-20 | [stellated_platonic_solids](projects/stellated_platonic_solids) |
| 2026-08-26 | [SubConsiousness](projects/SubConsiousness), [voice-horn-evolution](projects/voice-horn-evolution), [voice-horn-pocs](projects/voice-horn-pocs) |
| 2026-09-21 | [thermal_plume](projects/thermal_plume), [weihnachtspyramide](projects/weihnachtspyramide) |
| 2026-09-25 | [flesh-sheath-standing-wave](projects/flesh-sheath-standing-wave), [sacred_correspondant_forms](projects/sacred_correspondant_forms), [volumetric_sound_forms](projects/volumetric_sound_forms) |
| 2026-10-14 | [indirect-perception](projects/indirect-perception), [negative-space](projects/negative-space) |
| 2026-10-16 | [incense_vortex_ring](projects/incense_vortex_ring) |
| 2026-10-20 | [ogdoad_wind_chandelier](projects/ogdoad_wind_chandelier) |

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
