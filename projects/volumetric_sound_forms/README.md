# Volumetric Sound Forms

Browser tool: upload an mp3 or wav, get a printable STL of the sound in one of three forms.

- **Waveform solid** (the default). The waveform an audio editor draws, spun about its time axis:
  every window of the clip is a thin disc as wide as that moment is loud and as thick as its pitch is
  low, threaded on a spindle with pointed tips. Lies whole along its axis, or prints as two halves on their cut faces with no
  supports. `preview_waveform_solid.png` shows the spoken-vowel demo both ways.
- **Spectral cylinder.** The spectrogram stood up: time climbs the axis, frequency wraps the
  circumference, loudness per band raises the wall. Solid or hollow (wall + floor, open top).
- **Cavity form.** The clip played into a rigid cone-frustum vessel by a small speaker in its wall.
  The time-averaged standing-wave pressure is computed from the vessel's acoustic modes, the loudest
  share of the air is made solid, and those lobes are fused to a thin shell of the vessel. Whole or as
  two halves laid on their cut faces; open mouth or closed lid; the whole clip or one moment of it.

Everything runs client-side; nothing is uploaded.

The page is [`web/volumetric-sound-forms.html`](../../web/volumetric-sound-forms.html)
(one self-contained file, zero dependencies, per house rule). This directory holds the
verification harness for it.

## How the waveform solid is made

1. **Windows.** The trimmed span is cut into `wn` equal windows (fins mode) or one station per
   0.3 mm of `wlen` (waveform mode). A window's height is its peak `|s|` (`wmeas=peak`, what an editor
   draws) or its RMS, divided by the loudest window and raised to `wgam`.
2. **Pitch → thickness** (`wpitch=1`, fins mode). YIN per window: a frame of `W + τmax` samples
   centred on the window (`τmax = ⌈sr/50⌉`, `W = 2τmax`), difference function from one FFT
   cross-correlation plus prefix sums of x², cumulative-mean normalisation, first dip under 0.15
   refined by a parabola, search 50–1200 Hz. Windows under 2 % of the loudest, or whose deepest dip
   stays above 0.3, are unpitched and take the log-interpolated pitch of their voiced neighbours.
   Thickness is proportional to the period: the median voiced pitch gets `√(wfin · 0.8·spacing)`,
   clamped to `[wfin, 0.8·spacing]`. A clip with no pitched window keeps the fixed fin. The track is
   cached per samples and windows, so shape-only sliders don't recompute it.
3. **Profile.** Radius against axial position, pole to pole. Fins mode: per window a trapezoid centred
   in it, rim `wfin` (or the pitch thickness) thick at `rc + (wrad − rc)·a`, flanks leaning `wtap`° and capped at the window edge,
   so the interior profile is strictly increasing along the axis (no zero-area faces on a cut). Windows
   under 0.1 mm leave bare core (`wcore`). Ends close in cones of length `wtip` (0 = flat). Waveform mode
   joins the stations directly.
4. **Lathe.** `latheMesh`: one ring of `wseg` vertices per interior profile point, quad bands between
   rings, fans to the poles. Halves: θ ∈ [0, π] in `wseg/2 + 1` columns plus the cut face on z = 0; a flat
   end becomes a half-disc fanned from its θ = 0 corner. Closed and oriented by construction; the volume
   equals the stack of regular-polygon frusta (`latheVolume`) to float32 precision.
5. **Placement.** Whole: lying along x, resting on its largest disc (shown like the reference
   photograph). Halves: both cut faces down, side by side along y, 4 mm apart. All surface normals of a
   half point up or sideways, so the Overhang cell reads 0°.
6. The demo is `synthVoice`: a 2.4 s synthetic "ah" (ten harmonics through two formants, f0
   `165·2^(−0.9 t/2.4) + 3 sin(2π·5.5 t)` Hz, falling 0.9 of an octave,
   a long main envelope and a short second breath, a 16 ms random roughness from a 32-bit LCG).
   The 100 Hz tone stays one button away. The waveform solid skips the STFT; the other forms request it
   when chosen.

## How the spectral cylinder is made

1. `OfflineAudioContext.decodeAudioData` decodes any browser-supported codec at 44.1 kHz
   (a hand-written RIFF/WAVE reader is the fallback); stereo is folded to mono.
2. An inline Blob-URL worker (built from the page's own `/* @pure-begin */ … /* @pure-end */`
   block) runs the STFT: `nt` rows over the trimmed clip, Hann window of `fft` samples, up to 8
   sub-frames averaged per row when rows are far apart.
3. Bins → `nf` bands (log or linear between `fmin`/`fmax`), power summed per band, empty bands
   interpolated, dB against the clip maximum, floor / gamma / contrast, then the layout
   (`mirror` = no seam, `plain` = one seam, `transpose` = time around) and a wrapped box blur.
4. `meshCylinder` places every vertex on a ring grid and closes it with fans (solid) or an inner
   sheet + floor fan + rim annulus (hollow, inner wall = min-filtered radial offset, window radii
   rounded up so the wall is never thinner than set). The optional slope limit (`slope` 45/55/65°)
   fills a cone under every ledge so no face leans past that angle. Every edge is shared by exactly
   two triangles by construction; `solidCheck` reports V/E/F, Euler, signed volume, odd and reversed edges.
5. `writeSTL` streams binary STL in mm; fit-to-bed scales uniformly to the chosen limit (≤ 320).

## How the cavity form is made

1. **Modes.** Radius `a(z) = cr0 + (cr1 − cr0) z/ch`, lengths in units of `a0 = max(cr0, cr1)`.
   Cross-section modes `J_m(κρ) cos mθ` with κ a zero of `J_m′` (rigid wall; Bessel functions by
   Miller's recurrence, zeros by bisection). Up the height, the slowly-varying-duct equation
   `−(a² Z′)′ + κ² Z = k² a² Z` (Webster's horn equation at κ = 0) in finite volumes, rigid floor,
   rigid lid or `Z = 0` at an open mouth; eigenpairs by Sturm bisection and inverse iteration, for
   every mode with `k·a0 ≤ 22` (≈ 670 for the default vessel).
2. **Tuning.** `ctune` is the audio frequency that lands on the first sideways mode (κ = 1.8412) of
   the widest cross-section; the vessel then holds its lowest mode up to ≈ 12 × `ctune`.
3. **Drive.** A speaker on the wall at θ = 0, height `csrc·ch`, volume acceleration following the
   clip. Mode j's share of the time-averaged p² is `φ_j(src)² / (N_j² k_j⁴) · Σ_bins P(f) L(f/f_j)`,
   with the resonance `L(r) = 1/((1 − r²)² + (r/Q)²)` integrated exactly over each FFT bin (closed-form
   antiderivative). Cross terms between modes are dropped. Moment mode uses one STFT row; its fill
   scales with that row's ringing energy against the clip average.
4. **Field.** p² on a voxel grid (`cvox` mm; odd x/y counts through the axis; z layers on the floor's
   and lid's mid-planes), one radial table per group, Gaussian smoothing σ = `csmooth` mm,
   normalised to the loudest point in the cavity. The level is the quantile that makes the loudest
   `cfill` % of the cavity's air solid.
5. **Solid.** `F = max(lobes, shell)` in mm, lobes clipped to the vessel's envelope, wall ≥ 1.5 voxels.
   6-connected islands (lobes touching nothing) are removed; sealed pockets of quiet air are kept as
   voids. Halves: `min(F, y)` on the y ≥ 0 rows, so the cut is one plane.
6. **Mesh.** Marching tetrahedra (6 Kuhn tetrahedra per cube, per-layer shared edge vertices), closed
   and oriented by construction; halves rotated onto their cut faces and mirrored side by side.
   Genus and shell count are reported; the existing `solidCheck`, `fitToBed`, `writeSTL` are reused.

Parameters live in the URL hash (never the audio). `window.__selfTest()`, `__cavitySelfTest()`,
`__derived`, `__loadArrayBuffer(ab, name)`, `__exportSTL()`, `__runDemo('voice' | 'tone')` and `__whenIdle()` are
the test hooks. The load-time self-test covers the chart and the cheap cavity maths; the heavier
cavity checks run in the worker after first paint (the footer reports both).

## Verification

```
node run_pure.mjs                      # slices the pure block out of the HTML, runs both self-tests,
                                       # writes js_cylinder_*.stl (7), js_cavity_{open,closed,halves,moment}.stl,
                                       # js_wave_{fins,fins_halves,flat_rms_fixed_halves,waveform}.stl,
                                       # demo_field.*, demo_cavity.json, demo_cavity_field.f32, demo_wave.json,
                                       # demo_voice.f32, and checks the bare vessel against its closed-form volume
python validate_forms.py --meshlab     # trimesh + lib/validate_mesh.py gate over every *.stl here
python sound_form_reference.py --write-test-wav test_tone.wav --compare
                                       # numpy re-implementation. Chart: field, band map, surface and every
                                       # vertex (≈1e-7 / 4e-6 mm). Cavity: J_m' zeros vs scipy (1e-14), axial
                                       # eigenvalues vs LAPACK (1e-11), modal weights (1e-12), p² with exact
                                       # Bessel functions (1e-6), fill thresholds, and each cavity STL's volume
                                       # vs VTK marching cubes on the numpy implicit function (< 1 %). Waveform:
                                       # windows, YIN by direct sums vs the page's FFT (same voiced windows,
                                       # pitch to 1e-15), pitch vs the voice's own contour (< 1.2 %), thickness,
                                       # profile and every triangle identical (0.0 mm), volume vs the exact
                                       # faceted lathe (< 1e-6), every disc's rim vs its window's peak
python cavity_fem_check.py             # the duct approximation vs an axisymmetric finite-element solution
                                       # (scikit-fem, P2): ≤ 0.8 % for the default 4.8° wall, open and closed;
                                       # --cr1 30 (9.5°) ≤ 1.2 %; --cr1 20 (14°) ≈ 2.1 %
python sound_form_reference.py --wav test_tone.wav --layout mirror --hollow   # py_*.stl from a WAV
python render_waveform_preview.py      # preview_waveform_solid.png from the js_wave_*.stl
```

`validate_forms.py` expects two bodies in `*halves*` files (cavity and waveform alike) and accepts sealed pockets (inward shells,
negative volume). MeshLab's self-intersection filter misfires on stacks of exactly coplanar faces (any
vertical stretch of wall meshed by marching tetrahedra in float32): every face it flags is re-tested
against its neighbours with an independent float64 triangle–triangle test, and only real intersections
fail. The test is itself checked on two overlapping spheres, where it finds 372 real intersecting pairs.

Browser: serve the repo root (`python -m http.server 8317`), open
`http://localhost:8317/web/volumetric-sound-forms.html`, check the footer self-test line (the waveform solid opens by default; `#form=chart`, `#form=cavity` for the others), open
`#form=cavity`, then
`fetch('/projects/volumetric_sound_forms/test_tone.mp3').then(r=>r.arrayBuffer()).then(ab=>__loadArrayBuffer(ab,'t.mp3'))`.
`test_tone.mp3` was encoded from the test WAV with ffmpeg (`-codec:a libmp3lame -b:a 128k`).
A cavity rebuild at the defaults takes ≈ 0.5–0.8 s (91 × 91 × 107 voxels, ≈ 530 k triangles).

Generated files (`js_*.stl`, `py_*.stl`, `demo_*`, `test_tone.wav`) are git-ignored; the commands
above regenerate them in a couple of minutes.

## Publishing

`pwsh ./make-previews.ps1 -Only volumetric-sound-forms`, check the card in `site/index.html`, then
`pwsh ./deploy.ps1` (copies `web/*.html` into `site/` and runs `netlify deploy --prod`).

## Not built yet

Spherical cavity, spectral sphere, waveform vase, playback with a sweeping cursor, sound ring; a
second, longer FFT for the lowest bands of the chart, and a cap on sub-frames for very long clips.
