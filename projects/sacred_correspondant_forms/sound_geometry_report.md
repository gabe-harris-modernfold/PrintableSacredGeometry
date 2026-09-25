# Shape, Sound, and Energy: 3D-Printed Acoustic Geometries, Sacred Geometry, Cymatics, and Water Resonance

**Research report — compiled 11 August 2026**

## Source Labeling

Every claim below carries an inline link to a page that was actually fetched during this research. Each source is tagged as one of:

- **[MAINSTREAM SCIENCE]** — peer-reviewed papers, university/lab publications, engineering references, encyclopedic technical entries.
- **[CYMATICS DEMO]** — demonstration, instrument-maker, or practitioner documentation of visible sound patterns (sits between the two categories; often empirical but not peer-reviewed).
- **[ESOTERIC / OCCULT]** — Hermetic, metaphysical, radionic, vortex-mathematics, orgonomic, or free-energy literature and practitioner sites.

Esoteric frameworks are presented on their own terms, as the beliefs, practices, and claimed results of their proponents. Scientific material is presented on its own terms. No attempt is made here to adjudicate between them.

---

## SECTION 1 — 3D-Printed Acoustic Geometries (Mainstream Acoustic Science)

### 1.1 Horns and Waveguides: Impedance Matching and Directivity

A horn is an acoustic impedance transformer: it couples the high mechanical impedance of a small radiating diaphragm to the low impedance of open air. The governing relation is Webster's horn equation, and the horn's behavior follows from how its cross-sectional area S(x) expands along the axis.

#### Profile families and their formulas:

| Profile | Area function | Key parameter |
|---------|----------------|----------------|
| Exponential | S(x) = S₀·e^{mx} | flare constant m (m⁻¹); cutoff f_c = mc/4π |
| Conical | S(x) = π(r₀ + x·tan θ)² | half-angle θ |
| Tractrix | tractrix curve to mouth radius | mouth radius sets loading |

The exponential horn has a hard cutoff frequency below which it stops loading the driver; the same reference gives a worked example of an exponential horn with f_c = 80 Hz. A separate practical limit is set by mouth size: f_min ≈ c / (2 × mouth diameter).

**Directivity** is governed by the dimensionless ka* (wavenumber × mouth radius):

- ka < 1 → essentially omnidirectional radiation
- ka ≈ 2–3 → useful pattern control
- ka > 5 → strong beaming / narrowing

#### Printable design ranges

The same design tool exposes the parameter space a desktop printer can actually realize:
- Throat diameter: 10–100 mm (a 1-inch compression-driver throat = 25.4 mm)
- Flare constant: 1–50 m⁻¹
- Conical half-angle: 5–85°
- Horn length: 10–500 mm

Two concrete worked designs:

- **PA conical horn:** 30° half-angle → 60° nominal coverage, usable 2–16 kHz.
- **Studio oblate-spheroid waveguide (OSWG):** 25.4 mm throat, 150–250 mm length, 40–50° half-angle, mouth ≈ 180–220 mm, giving pattern control down to ~1.2 kHz.

Those numbers are directly buildable on a 250 mm-bed FDM machine; the 180–220 mm mouth is the binding constraint.

### 1.2 Helmholtz Resonators: Cavity and Neck Geometry

The Helmholtz resonator — a cavity of volume V₀ with a neck of cross-sectional area A and length L — is the simplest geometry-defined acoustic resonance and the workhorse of 3D-printed absorbers.

**Canonical formula:**

```
f_H = (v / 2π) √(A / (V₀ · L_eq))
```

with the end-corrected neck length L_eq = L_n + 0.3 D (D = neck diameter) and adiabatic index γ = 1.4 for air.

An equivalent thermodynamic form used in acoustic-materials engineering is:

```
f = (1 / 2π) √(γP₀ / ρ) · √(A_n / (V_c · L_n))
```

with ρ = 1.2 kg/m³ for air at standard conditions.

**Design consequence for printing:** Because f scales as √(A/V₀L), halving the cavity volume raises the resonance by only √2 (≈ +41%), while the neck is the sensitive lever. This is why printed absorber arrays vary neck geometry, not cavity size, to tune across a band.

### 1.3 Acoustic Metamaterials: Labyrinthine Channels, Coiling, and Phononic Crystals

Space-coiling / labyrinthine unit cells fold a long acoustic path into a subwavelength footprint, producing a large effective refractive index. Duke University's Xie, Cummer et al. (Applied Physics Letters 103, 201906, 2013) fabricated tapered labyrinthine cells in ABS by FDM and reported:

- Horn-port cells with inner paths rotated ~30°, port aperture 6.9 mm (4.3× the channel width)
- Anisotropic spiral coiling coefficient g ≈ 3.1
- Isotropic spiral g ≈ 3.95, wall thickness 0.5 mm
- Twisted spiral: height 5.08 cm, twist 18°, cross-section 25 × 25 mm, g ≈ 3.7
- Effective index relation: n_eff ≈ g·2πc/(ωD)

#### 3D-printed labyrinthine absorber panel (Wunderlich space-filling curve)

An arXiv study designed and tested a printed panel of quarter-wave labyrinthine unit cells targeting 800–1300/1400 Hz, centered on 1 kHz:

| Parameter | Value |
|-----------|-------|
| Unit-cell channel width w | 6 mm |
| Unit-cell thickness t | 9.6 mm |
| Unit-cell footprint a | 19.6–19.7 mm |
| Aperture height h | 0.5 mm |
| Aperture area S_a | 36 mm² |
| Tuning sensitivity | df/dt = −42.7 Hz/mm, df/da = −145.3 Hz/mm |
| Analytical perfect absorption | 1370 Hz, bandwidth Δf = 110 Hz |
| Measured absorption peak | ~1480 Hz, α = 0.95, Δf ≈ 390 Hz |
| Panel cells | a₁ = 19.5 mm, a₂ = 26.0 mm; macro-cell 156 × 156 mm; 20 macro-cells; panel 628 × 784 mm |
| Subwavelength ratio | λ ≈ 250 mm ≈ 26× the unit-cell thickness |
| Impedance-tube test range | 100–5000 Hz in 2 Hz steps |

Those two tuning derivatives (−42.7 Hz/mm and −145.3 Hz/mm) are the most directly actionable numbers in this section for a maker: they tell you exactly how much your printed cell drifts per millimetre of dimensional error.

#### Review-level metamaterial numbers:

- At 20 kHz, λ = 17.25 mm, setting the largest practical "meta-atom"; subwavelength cells run around λ/10.
- A 3D acoustic cloak was demonstrated at 3 kHz; a tunable cloak spanned 3000–6000 Hz.
- A phononic crystal showed a pass band near ~2 kHz and a band gap near ~4 kHz.
- A Helmholtz-based metasurface made by vat photopolymerization achieved full 2π phase coverage in π/4 steps at a thickness of λ/30.
- A gradient coiled metamaterial achieved pressure magnification ≈ 80×, versus 50× for a straight channel — i.e., geometry alone bought a 60% gain in concentrated pressure.
- A Hilbert fractal metasurface with 400 μm walls absorbed across 225–1175 Hz.

#### Acoustic metalens, printed in resin

A 2024 study from Universitas Gadjah Mada printed a metalens on an Anycubic Photon Mono X SLA printer using Anycubic 405 nm resin at 50 μm layer thickness. The lens comprised 22 unit cells (11 distinct types, mirrored, i.e. 11 symmetric pairs) built from Helmholtz resonators and labyrinthine chambers, with fixed dimensions S₀ = 22 mm, w₀ = 1 mm, L₀ = 108.6 mm, achieving full 2π phase shift. Testing used a waveguide of two glass plates (1000 × 500 × 3 mm), 8 × 5 W speakers (40 × 70 mm), and an LT Lutron SL-4013 meter (31.5–8000 Hz, 30–130 dB), driven at 1000 Hz, with a further design discussed at 4 kHz.

### 1.4 Real 3D-Printed Acoustic Concentrators and Manipulators

#### The "acoustic spanner" (Nature Scientific Reports, 2019)

This is the most directly replicable printed-geometry-concentrator project found. Researchers 3D-printed a PLA funnel-shaped concentrator on an Ultimaker 2, designed in Fusion 360, with 8 sound-guides channeling 8 × Visaton K28WP 28 mm 1 W speakers arranged on a 10 cm-diameter circle down to a 3 cm funnel outlet. Driven at 650 Hz with orbital angular momentum states ℓ = +1, 0, −1, it rotated a ~3 g packing peanut (40 × 30 mm) and ~2 mm graphite flakes floating on water. Drive electronics: Arduino Due DDS, PWM at 328 kHz, filtered at 12.5 kHz; design and firmware are public at the AcousticSpanner GitHub repository — **[MAINSTREAM SCIENCE]**.

#### TinyLev acoustic levitator

The standard open acoustic-levitation build operates at 40 kHz (λ = 8.65 mm at 25 °C), levitating samples up to ~4 mm (λ/2), using 72 × 10 mm transducers in rings of 6/12/18 on a 6 cm radius of curvature, driven at ~20 V (up to 70 Vpp) for about 10 W; the larger BigLev variant reached densities of 6.5 g/cm³.

#### Acoustic holograms

Melde, Fischer and colleagues (Nature 537:518, 2016; DOI 10.1038/nature19755) 3D-printed plastic relief plates whose local thickness encodes phase delay, projecting arbitrary sound fields that assembled particles into a Picasso dove and a holographic stack rendering the characters "1", "2", "3". This is, in effect, a printed shape acting as a programmable field-shaping device — the mainstream-physics analogue of "form directing energy."

---

## SECTION 2 — Sacred Geometry, Hermetic and Occult Concepts of Shape Amplifying Energy

Everything in this section is presented as the framework, claim, or reported measurement of its proponents.

### 2.1 Rodin Coils, Vortex-Based Mathematics, and Toroidal Windings

#### The 1-2-4-8-7-5 doubling circuit

Marko Rodin's Vortex-Based Mathematics (VBM) reduces numbers to digital roots. Repeated doubling produces the closed loop 1 → 2 → 4 → 8 → 7 → 5 → (1), while 3 and 6 oscillate against each other and 9 stands as the invisible axis / source — the "emanation" from which the doubling circuit unfolds. Fibonacci digital roots form a 24-step cycle summing to 108; Randy Powell's extension is called the ABHA torus.

#### Claimed performance figures from Rodin Aerodynamics

The primary proponent document states that:
- "crude versions of the Rodin Coil, created and tested by leading scientists and engineers, show 60% more efficiency than anything presently used in antennas, computer research, or life-saving medical devices"
- "the Rodin Torus Coil tested at Hewlett Packard and determined to be creating 62.5% greater magnetic output than the present day standard wound electrical coils"
- the coil's magnetism permits "a minimum 20% reduction in copper"
- a prototype ceiling fan tests "600% more efficient than a box fan"
- the coil is "60% more sensitive than any antenna existing"

Geometrically, it describes the central figure in the torus hole as "a polygon of 24 sides for each completed wrap," with "six different self-referencing axis configurations" and "perfect mathematical coherence on all six axes".

#### Independent maker bench data

A widely circulated test report gives each Rodin winding as 18 turns of #26 wire, 4.5 Ω, with a load inductor of 1800 turns of #28 wire, 128 Ω, driven with a 1200 kHz square wave at 5 V and 10 V, producing a 125 V spike; a motor test at 12 VDC / 10,000 RPM yielded 220 VAC, with the author's claim that the setup "shows OU by as much as 1700%".

#### Bearden / JLN Labs framing

Jean-Louis Naudin's 1997 replication cites Tom Bearden's 1996 "Introduction to Rodin Coil design," reporting that the coil behaves inductively on its center axis, that inserting a ferrite rod changes its inductance, and framing the effect in terms of a curl-free magnetic vector potential (A-potential) and the Aharonov–Bohm effect.

#### Peer-reviewed magnetometry of a Rodin coil

A study in Revista Mexicana de Física 70, 030901 (May–June 2024) mapped a Rodin coil's field with FEM plus a scanning magnetometer. 

**Geometry:** Dodecagonal form, 36 poles, winding angle 150°, 12 poles at 30° separation, toroid external diameter 9.95 cm, aperture 2.9 cm. The field was scanned over a 16 × 16 × 8 cm volume in 0.5 cm steps using an ATmega2560 controller.

### 2.2 Metatron's Cube and the Platonic Solids

Metatron's Cube is constructed from the 13 circles of the Fruit of Life, connected by 78 lines — exactly n(n−1)/2 for n = 13 — and is held to contain all five Platonic solids:
- Tetrahedron (4 faces / 4 vertices)
- Cube (6/8)
- Octahedron (8/6)
- Icosahedron (20/12)
- Dodecahedron (12/20)

In the popular esoteric literature the solids carry elemental correspondences and the figure is explicitly linked to 432 Hz tuning.

### 2.3 Vogel-Cut Quartz: Facet Angles and Cut Theory

Marcel Vogel (IBM research chemist, later crystal researcher) originated a specific quartz cut whose angles are its defining claim.

#### Core angles:

- **Universal / "female" end:** Cut to 52°, justified as "the interlattice angle of quartz is 90 degrees minus the basal angle which is 38 degrees. Therefore the interlattice angle is 52 degrees," and "this is also the angle of the Cheops pyramid."
- **"Male" end:** "Formed with a more acute degree angle of around 60 degrees that is resonant with the individual."
- **Golden ratio:** Φ = 1.618 is invoked as "same as the Great Pyramid."

The primary Vogel monograph gives the pyramid-angle figure to the second: the large end cut to 51°51'51", the Great Pyramid of Cheops angle.

#### Water-geometry rationale invoked in the cut theory:

- Water's intermolecular bond angle: 104.5°
- Tetrahedral angle: 109.5°
- Flat hexagon interior angle: 120°
- Phase angle in a six-molecule ring: 60° (360/6)
- For a seven-molecule ring: 51.43° (360/7)
- Hexagonal structure ratio: √3 = 1.732
- Water "bending angle": typically 20°–60°, most commonly 26°, second most common 54°

#### Facet counts

Vogel crystals are cut with 4, 6, 8, 12, 13, or 24 facets, the first form having been a 4-sided bipyramid, with the doctrine that "the more facets the more amplification," and the design attributed to a vision of the Tree of Life. A practitioner authentication guide lists 4, 6, 8, 10, 12, 13, and 16 facets as authentic counts and specifies the male end must be more acute than the female end.

#### Vogel's water experiments and reported readings

Vogel circulated 1000 cc of water through a Pyrex tube coiled seven clockwise (right-hand) turns, tubing ¾ inch diameter, coil diameter approximately six inches, chamber sealed at 70 °F, drawing a 50 cc sample per pass for pH, conductivity, and "Omega 5" readings. He reported that:

- Treated water developed a new set of UV absorption bands around 3000–4000 angstroms
- "The magnetic moment of this structured water is increased by 0.07 gauss"
- In the radionic framing, "the crystal is cut to resonate radionically to the rate of water which is 454" — described as "the radionic signature, or the fundamental numerical identity of the information band of the water itself" (no Hz unit is given for 454)
- A larger coil experiment of 250 turns on a 6-inch coil is also documented

### 2.4 Pyramid, Cone, and Torus: Orgone, Pyramid Power, Radionics

#### Wilhelm Reich's orgone accumulator

Invented in 1940 and described in *The Orgone Energy Accumulator* (Orgone Institute Press, Rangeley, 1951, pp. 12, 16–19), the device is built of alternating organic layers on the outside and metallic layers on the inside, where "one layer" in Reich's counting means two physical layers. Reich rated:
- Plastic as a good absorber and wood as a poor one
- Specified only iron for medical use
- Gave 40–50% relative humidity as optimal
- Held the accumulator to be more efficient nearer the equator and at higher altitude

#### Reported orgone measurements

A 2011 orgonomy-conference review by Maglione and Mazzocchi compiles measurement geometry and results: accumulators described as 3-fold, 4-fold, 5-fold, 8-fold, 15-fold, and 24-fold; a 3-fold cubic accumulator of 10 cm size; a Geiger–Müller counter at 50 cm from the accumulator walls and 50 cm above the floor; control runs at ~1 m outside; thermometers with 0.03 °C accuracy. The paper reports:

- An "electroscope discharge many times longer (up to 9) within an orgone accumulator compared to that measured outside"
- The ratio falling toward 1 when humidity exceeded 70%
- Longer discharge times below 30% humidity
- A transition zone at 57–65% relative humidity
- A claimed "halo of orgone energy" extending 7.5–8.5 m, 30–35 m, or between 7.5 and 42.5 m depending on the case and season, peaking around 3.5 m from the walls

#### Pyramid power

The modern tradition traces to Antoine Bovis in the 1930s, and to Karel Drbal's Czechoslovak patent #91,304, applied for in 1949 and published in 1959, titled the "Pharaoh's shaving device" (razor-blade sharpener). It reached the West through:
- Ostrander & Schroeder's 1970 *Psychic Discoveries Behind the Iron Curtain*
- Max Toth and Patrick Flanagan's 1970s *Pyramid Power*
- In the 1975–76 season the Toronto Maple Leafs placed a pyramid by their bench (players sitting under it for 4 minutes)
- Martin Gardner treated the topic in *Scientific American* in June 1974

The Drbal patent text itself is available in full (Patent 91,304 PDF — **[ESOTERIC / OCCULT]** primary document).

### 2.5 Scalar Wave Theory and Geometric Coil / Antenna Designs

Konstantin Meyl's scalar-wave program is built around specific hardware:
- **Demo-Kit generator:** operates at 4–8 MHz
- **Experimental-Kit:** spans 135 kHz–10 MHz, using pancake coil pairs dimensioned at half or twice the wire length for resonance

Meyl claims:
- The Tesla magnifying transmitter propagated scalar waves at approximately 1.5c
- Dates his discovery of the "potential vortex" to 1990
- His own paper appeared in the *Journal of Scientific Exploration* 15(2), 2001 (Meyl, JSE 15(2) PDF — **[ESOTERIC / OCCULT]**)

---

## SECTION 3 — Cymatics

### 3.1 Chladni Plates: History and Physics of Nodal Patterns

Ernst Chladni (1756–1827) and Sophie Germain (1776–1831) established the modern study of plate nodal figures.

Square plate mode shapes follow cos(mπx/L)·cos(nπy/L), with frequency scaling as f ∝ m² + n². Chladni's empirical law for circular plates is f ~ (m + 2n)².

#### Fully worked circular-plate example with measured mode frequencies

For a circular stainless baking pan, 0.5 mm thick, 0.229 m diameter (a = 0.115 m), E = 205 GPa, ρ = 7700 kg/m³, μ = 0.3, giving flexural rigidity D = Eh³/12(1−μ²) = 2.490 N·m and mode frequencies from f_n = (λ²/2πa²)·√(D/ρh):

| Frequency | Mode shape |
|-----------|-----------|
| 110 Hz | "oil-can" (fundamental umbrella mode) |
| 230 Hz | 1 diametral line |
| 400 Hz | 2 diametral lines |
| 495 Hz | 1 circular node |
| 600 Hz | 3 diametral lines |
| 720 Hz | 1 circular + 1 diametral |
| 980 Hz | 1 circular + 2 diametral |

This is the single most useful table in this report for diagramming: it maps a specific Hz to a specific named nodal geometry on a specifiable, printable-scale plate.

#### Square plate resonance sweep

A classroom demonstration drove a square Chladni plate from 932 → 1459 Hz, passing through resonances at 955 Hz, 1055 Hz, and 1377 Hz.

#### PASCO plate resonance catalogue

The Oklahoma junior-lab study lists measured plate resonances at 174, 197, 254, 345, 396, 400, 490, 795, 950, 1060, 1400, 1725, 1900, 2105, 2260, and 2700 Hz. For its circular plate it identifies:
- n = 4 at 112 Hz and 270 Hz (predicted 255 Hz)
- n = 5 at 2700 Hz
- n = 6 at 4000 Hz
- m = 7 at 5750 Hz
- (2,2) mode at 450 Hz

#### Chladni patterns as a particle-manipulation tool

A Nature Communications paper (article 12764) used chromatic-note-addressed Chladni modes for micro-manipulation, spanning C6 = 1.047 kHz to A#10 = 29.83 kHz — 59 chromatic notes, of which 34 are C-major. Specific working notes cited include:
- F8 = 5588 Hz
- C#6 = 1109 Hz
- A6 = 1760 Hz
- Range F6 = 1397 Hz to E9 = 10548 Hz

Plate resonance peaks occurred at 11.9, 19.6, 23.9, 27.1, 34.5, and 43.5 kHz, characterized with a 200 Hz–50 kHz chirp. At 29.83 kHz, λ = 16 mm giving a 4 mm node-to-antinode spacing, with positioning resolution of 4 μm.

### 3.2 Hans Jenny and the Coining of "Cymatics"

Hans Jenny (1904–1972) coined the term cymatics and published *Kymatik* volume 1 in 1967 and volume 2 in 1972. His method included:
- Quartz sand on a 60 cm black drum membrane excited by singing through a cardboard pipe
- The tonoscope for rendering the voice visible

His summary observation: low tones produce simple figures, high tones produce complex ones (Wikipedia, Hans Jenny (cymatics) — **[MAINSTREAM SCIENCE]** / biographical). Jenny's full book is available as a scan (*Cymatics: A Study of Wave Phenomena and Vibration*, Monoskop PDF — **[CYMATICS DEMO]** / primary source).

### 3.3 The CymaScope and Modern Water-Surface Cymatics

The instrument lineage runs from:
- A brass plate excited with a violin bow
- Through Margaret Watts-Hughes' Eidophone and Mary D. Waller's dry-ice technique
- To John Stuart Reid's CymaScope Pro, which imprints sound onto a water surface for imaging

#### Peer-reviewed-adjacent CymaScope study

Rupert Sheldrake and Merlin Sheldrake published "Determinants of Faraday Wave-Patterns in Water Samples Oscillated Vertically at a Range of Frequencies from 50–200 Hz" in the October 2017 edition of *Water Journal*, using the CymaScope with precise control of acoustic energy into the visualising cell. They tested patterns at six temperatures (5, 10, 15, 20, 25, 30 °C) for three sample frequencies and reported that water temperature:
- "did not alter the overall morphology of the pattern"
- "had a negligible effect on the degree to which patterns were expressed"

#### Systematic cuvette-geometry study (2025)

Geffen, Val Baker, Gentile and Oomen (WATER, DOI 10.14294/WATER.2025.3, published 9 September 2025) classified Faraday-wave resonances across cuvette geometries using a CymaScope Pro MF+ with a voice-coil-motor-driven piston, a 32-diode light ring, and a Canon EOS 5D Mk IV.

| Item | Value |
|------|-------|
| Instrument compatibility range | ~250 Hz |
| Investigated frequency range | 8–256 Hz (16–256 Hz for rectangular, square, circular cuvettes) |
| Triangular cuvette | 50 mm length, 2.1 ml water |
| Rectangular cuvette | 24.25 × 15 mm, 1.1 ml |
| Circular cuvette | 24.25 mm diameter, 2.3 ml |
| Square cuvette | 25.2 mm length, 3.3 ml |
| Circular cuvette at 40 Hz | 14-fold and 4-fold patterns compared |
| Reference example | 69 Hz (Figure 22) |
| "Unidentified anomaly" | 128.4 Hz at −19.4 dB — pattern forms in one trial, not another |
| Observed rotational motion | 10°, 45°, and 360°+ |

They report that pattern formation requires:
- A compatible reference frequency
- Sufficient reference amplitude (which increases with frequency for every cuvette shape)
- Sufficient Time to Full Expression (TFE), which decreases with both amplitude and frequency

Shifting ±a few Hz from a reference frequency either reproduces the same pattern or yields a completely new one.

### 3.4 Water-Drop Cymatics (Lauterwasser) and Popular Demonstrations

Alexander Lauterwasser's work vibrating a single water drop between 40 and 120 Hz produced sequential polygonal figures — triangles, squares, pentagons, hexagons, heptagons — as frequency rose. An ARTE documentary clip on his work references:
- A sand "dance" around 11,000 Hz
- A figure at 380 Hz (ARTE / Lauterwasser video — **[CYMATICS DEMO]**)
- A French-language book excerpt captions a water sound-image at 1,330 Hertz (Calaméo book excerpt — **[CYMATICS DEMO]**)

Common maker-level demonstrations include water on a speaker cone swept 20–220 Hz (YouTube water-on-speaker demo — **[CYMATICS DEMO]**).

### 3.5 Esoteric Interpretations of Cymatics

Within the esoteric reading, cymatic figures are treated as evidence that sound is the formative principle of creation and that specific tones generate specific sacred-geometric forms. Practitioner demonstrations circulate showing:
- The Schumann resonance at 7.83 Hz applied to a thin layer of distilled water (Sacred Geometrics, 7.83 Hz water cymatics post — **[ESOTERIC / OCCULT]**)
- 528 Hz water cymatics (528 Hz cymatics video — **[ESOTERIC / OCCULT]**)

The link between Metatron's Cube and 432 Hz noted in Section 2.2 belongs to the same interpretive current.

---

## SECTION 4 — Water, Frequency, and "Fluid Mirrors"

### 4.1 Faraday Waves: The Physics Underneath Water Cymatics

When a fluid layer is oscillated vertically, standing surface waves appear at half the driving frequency (subharmonic response), forming stripes, hexagons, squares, and quasiperiodic patterns. Michael Faraday described the phenomenon in an appendix to an 1831 *Philosophical Transactions* paper.

#### The canonical experimental survey

Kudrolli and Gollub (Physica D 97, 1996, pp. 133–154) systematically mapped primary patterns and secondary instabilities using:
- A 32 cm Delrin container with silicone oil 0.3 cm deep
- Drive range: 20 Hz to over 200 Hz at accelerations up to 15 g via a 2200 N shaker
- Temperature-controlled to 25.0 ± 0.1 °C
- Control parameter ε = (a − a_c)/a_c

The companion theoretical work notes that minimization of the associated Lyapunov functional "predicts a primary bifurcation to a standing wave pattern of square symmetry for capillary-dominated surface waves, in agreement with experiments" (arXiv patt-sol/9607006 — **[MAINSTREAM SCIENCE]**).

#### Practical rule for builders

A printed cuvette or plate driven at f will show surface structure at f/2, and the pattern's symmetry (square vs. hexagonal vs. stripe) is set jointly by:
- Depth
- Viscosity
- Drive amplitude
- Boundary shape

This is exactly why the cuvette-geometry study in §3.3 found different fold-symmetries at the same 40 Hz.

### 4.2 "Fluid Mirrors": Named Installations and Instruments

Lachlan Turczan (b. 1992) is the principal contemporary artist working the fluid-mirror concept. His piece *Speaking Water* (2020) is a 1 × 1 × 1 ft work of water, aluminum, and a transducer, using a parabolic mirror and subsonic frequencies to form a liquid lens. His exhibition *Lensing Water* ran 25 July – 14 September 2025, comprising:
- Sympathetic Resonance
- Contact Lens
- Watergram
- Annular Optic

This followed the *Lucida* (I–VI) series made with Google for Milan Design Week 2025.

#### Liquid-mirror telescopes

Liquid-mirror telescopes are the engineering realization of the same idea: a rotating dish of mercury forms a true paraboloid:

- **Large Zenith Telescope:** 6 m across, rotating at ~8.5 rpm, cost about $1M, decommissioned in 2016
- **Proposed successors:** ALPACA (8 m) and LAMA (66 × 6.15 m)
- **Historical:** Henry Skey built the first laboratory liquid-mirror telescope in 1872, following Capocci in 1850
- **Current:** The International Liquid Mirror Telescope is a 4 m instrument at Devasthal, Nainital (Wikipedia, International Liquid Mirror Telescope — **[MAINSTREAM SCIENCE]**)

### 4.3 Masaru Emoto's Claims (Presented as His Claims)

Masaru Emoto (22 July 1943 – 17 October 2014) claimed that human consciousness, words, and music alter the crystalline structure of freezing water. He published:
- *Messages from Water* (1999)
- *The Hidden Messages in Water* (2004)
- A 2008 paper in the *Journal of Scientific Exploration*
- He declined James Randi's $1 million challenge in 2003

#### Emoto's stated methodology (from his own institute)

- 50 or 100 Petri dishes per sample
- Typically no more than 15% of dishes yielding "beautiful" crystals
- A 17-pattern scoring table scored 0–100 (100 = "Beautiful", 80 = "Hexagonal", 0 = "None")
- A 90-second exposure to an MRA-generated "Hado of immunity"
- Music exposures including Bach's *Air on the G String* and Chopin's *Étude Op. 10 No. 3*

### 4.4 Structured / "Fourth-Phase" and Hexagonal Water Claims

Gerald Pollack's EZ (exclusion zone) or "fourth phase" water is described as:
- Extending millions of molecular layers out from a hydrophilic surface
- Being negatively charged
- Being built by UV, visible, and infrared light — with near-infrared the most capable

Marcel Vogel's structured-water claims (§2.3) belong to the same family:
- New UV bands at 3000–4000 Å
- A +0.07 gauss shift in magnetic moment after seven right-hand-turn circulation around a charged, 52°-cut crystal

### 4.5 Specific Frequencies Practitioners Associate with Water

#### The Solfeggio set

Joseph Puleo derived the core six frequencies in the 1990s from the Book of Numbers 7:12–83:
- 396, 417, 528, 639, 741, 852 Hz

Later extended with:
- 174, 285, and 963 Hz

Leonard Horowitz popularized 528 Hz as the "love frequency," co-authoring the 1999 book with Puleo. The same source notes A440 was adopted in 1939 and that historical tuning pitches ranged 415–466 Hz.

#### 432 vs 440 Hz

- **Scientific pitch:** Fixes C4 = 256 Hz and correspondingly A4 = 430.54 Hz, traced to Sauveur (1713)
- **Verdi (1874):** Used the A = 435 Hz diapason normal
- **Schiller Institute:** Has campaigned since 1988 for A = 432

The Verdi/432 controversy has received academic treatment in:
- *Journal of Popular Music Studies* (Perfect Pitch: 432 Hz Music and the Promise of Frequency — **[MAINSTREAM SCIENCE]**)
- *Nineteenth-Century Music Review* (Tuning Sounds in Italy, 1750–1885 — **[MAINSTREAM SCIENCE]**)

#### A laboratory test of 432 vs 440 Hz on water

A 2022 aquaphotomics pilot study (Molecules 27(19):6332, DOI 10.3390/molecules27196332) exposed:
- **Samples:** Purified water and Yunosato Gold mineral water
- **Frequencies:** 432 Hz and 440 Hz tones
- **Equipment:** YAMAHA MOTIF XF8 synthesizer, BOSE L1 speakers, samples a metre from the source
- **Duration:** 10 min exposures
- **Data:** 2 replicates × 4 vials per frequency × 5 consecutive scans = 80 spectra, MicroNIR 1700-ES

The authors report that "even though there was only 8 Hz difference in applied sound frequencies, the change of absorbance at water absorbance bands was specific for each frequency and also water-type-dependent," and note prior literature claiming 432 Hz music "can significantly decrease the heart rate and slightly decrease mean respiratory rate values".

A practitioner page presents side-by-side water images attributed to the 432 vs 440 Hz contrast, describing "the beautiful symmetry of the 432 Hertz musical note versus the complete lack of harmony" of 440 (Crystal Blue, Healthy, Natural Frequencies — **[ESOTERIC / OCCULT]**). A cymatics-in-hertz practitioner reference notes that "in 432 Hz tuning, 38.18 Hz = note D#1/Eb1, with a wavelength of 9.03 meters" (Satori Rei, Cymatics in Hertz — **[ESOTERIC / OCCULT]**).

---

## KEY FREQUENCIES REFERENCE TABLE

Every specific Hz value located across all four sections, with what it is claimed or shown to do.

| Frequency | What it does / is associated with | Category | Source |
|-----------|-----------------------------------|----------|--------|
| 7.83 Hz | Schumann resonance applied to thin-layer distilled water in cymatic demo | Esoteric-occult | Sacred Geometrics |
| 8–256 Hz | Full investigated range for Faraday-wave cuvette classification study | Cymatics demo | Geffen et al., WATER |
| 16–256 Hz | Range used for rectangular, square, and circular cuvettes | Cymatics demo | Geffen et al., WATER |
| 20–200+ Hz | Kudrolli & Gollub Faraday-wave survey drive range (silicone oil, 0.3 cm deep) | Mainstream science | Physica D 97:133 |
| 20–220 Hz | Popular water-on-speaker cymatics sweep | Cymatics demo | YouTube demo |
| 38.18 Hz | D#1/Eb1 in 432 Hz tuning; λ = 9.03 m | Esoteric-occult | Satori Rei |
| 38–689 Hz | Audible music range in a cited water-exposure study | Mainstream science | Molecules 27:6332 |
| 40 Hz | Circular cuvette: 14-fold and 4-fold patterns compared at same frequency | Cymatics demo | Geffen et al., WATER |
| 40–120 Hz | Lauterwasser single water drop → triangles, squares, pentagons, hexagons, heptagons | Cymatics demo | Generationen Tandem |
| 50–200 Hz | Sheldrake & Sheldrake CymaScope Faraday-wave determinants study | Cymatics demo | CymaScope |
| 69 Hz | Reference frequency example, cuvette pattern classification | Cymatics demo | Geffen et al., WATER |
| 80 Hz | Worked exponential-horn cutoff f_c example | Mainstream science | Maker Portal |
| 100–5000 Hz | Impedance-tube test range (2 Hz steps) for printed labyrinthine absorber | Mainstream science | arXiv 2110.05026 |
| 110 Hz | Circular plate (0.229 m, 0.5 mm steel): "oil-can" fundamental mode | Mainstream science | Vibrationdata |
| 112 Hz | Circular Chladni plate n = 4 mode | Mainstream science | OU junior lab |
| 128.4 Hz (at −19.4 dB) | "Unidentified anomaly": pattern forms in one trial, not another | Cymatics demo | Geffen et al., WATER |
| 174 Hz | PASCO plate resonance; also extended-Solfeggio tone | Mainstream science / esoteric-occult | OU junior lab / Muzen Cosmos |
| 197 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 225–1175 Hz | Absorption band of Hilbert-fractal metasurface (400 μm walls) | Mainstream science | PMC8226526 |
| 230 Hz | Circular plate: 1 diametral nodal line | Mainstream science | Vibrationdata |
| ~250 Hz | Compatibility ceiling of the CymaScope instrument | Cymatics demo | Geffen et al., WATER |
| 254 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 256 Hz | Scientific pitch C4 (→ A4 = 430.54 Hz) | Mainstream science | Scientific pitch |
| 270 Hz | Circular plate n = 4 mode (predicted 255 Hz) | Mainstream science | OU junior lab |
| 285 Hz | Extended Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 345 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 380 Hz | Lauterwasser cymatic figure (ARTE documentary) | Cymatics demo | ARTE clip |
| 396 Hz | Solfeggio tone; also a measured PASCO plate resonance | Esoteric-occult / mainstream science | Muzen Cosmos / OU junior lab |
| 400 Hz | Circular plate: 2 diametral lines; also PASCO plate resonance | Mainstream science | Vibrationdata / OU junior lab |
| 415–466 Hz | Historical range of concert tuning pitch before standardization | Esoteric-occult source, historical claim | Muzen Cosmos |
| 417 Hz | Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 430.54 Hz | A4 under scientific pitch (C4 = 256 Hz) | Mainstream science | Scientific pitch |
| 432 Hz | Schiller Institute tuning campaign (since 1988); linked to Metatron's Cube; lab-tested NIR effect on water distinct from 440 Hz | Esoteric-occult + mainstream science | Scientific pitch / Gaia / Molecules 27:6332 |
| 435 Hz | Verdi's 1874 diapason normal | Mainstream science | Scientific pitch |
| 440 Hz | ISO 16 standard A4 (standardized 1955; A440 adopted 1939); comparison tone in water study | Mainstream science | Molecules 27:6332 / Muzen Cosmos |
| 450 Hz | Circular plate (2,2) mode | Mainstream science | OU junior lab |
| 490 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 495 Hz | Circular plate: 1 circular nodal line | Mainstream science | Vibrationdata |
| 528 Hz | Solfeggio "love frequency" (Horowitz); water cymatics demonstrations | Esoteric-occult | Muzen Cosmos / 528 Hz demo |
| 600 Hz | Circular plate: 3 diametral lines | Mainstream science | Vibrationdata |
| 639 Hz | Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 650 Hz | Drive frequency of the 3D-printed PLA "acoustic spanner" (OAM ℓ = ±1) | Mainstream science | Sci. Rep. 9:43046 |
| 720 Hz | Circular plate: 1 circular + 1 diametral | Mainstream science | Vibrationdata |
| 741 Hz | Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 795 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 800–1300/1400 Hz | Target absorption band of 3D-printed Wunderlich labyrinthine panel | Mainstream science | arXiv 2110.05026 |
| 852 Hz | Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 932–1459 Hz | Square Chladni plate demonstration sweep | Mainstream science | ChemEd X |
| 950 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 955 Hz | Square plate resonance in sweep | Mainstream science | ChemEd X |
| 963 Hz | Extended Solfeggio tone | Esoteric-occult | Muzen Cosmos |
| 980 Hz | Circular plate: 1 circular + 2 diametral | Mainstream science | Vibrationdata |
| 1000 Hz | Drive frequency for SLA-printed 22-cell acoustic metalens focusing test | Mainstream science | J. Eng. Appl. Technol. 5(1) |
| 1 kHz | Center design frequency of printed labyrinthine absorber panel | Mainstream science | arXiv 2110.05026 |
| 1047 Hz (C6) | Lowest chromatic note used in Chladni particle-manipulation platform | Mainstream science | Nat. Commun. 12764 |
| 1055 Hz | Square plate resonance in sweep | Mainstream science | ChemEd X |
| 1060 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 1109 Hz (C#6) | Working note in Chladni particle manipulation | Mainstream science | Nat. Commun. 12764 |
| ~1.2 kHz | Lowest pattern-control frequency of a 25.4 mm-throat studio OSWG | Mainstream science | Maker Portal |
| 1330 Hz | Captioned water sound-image ("Image sonore d'eau, 1 330 Hertz") | Cymatics demo | Calaméo excerpt |
| 1370 Hz | Analytical perfect-absorption frequency of printed labyrinthine cell (Δf = 110 Hz) | Mainstream science | arXiv 2110.05026 |
| 1377 Hz | Square plate resonance in sweep | Mainstream science | ChemEd X |
| 1397 Hz (F6) | Lower bound of a working note range in Chladni manipulation | Mainstream science | Nat. Commun. 12764 |
| 1400 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| ~1480 Hz | Measured absorption peak of printed labyrinthine cell, α = 0.95, Δf ≈ 390 Hz | Mainstream science | arXiv 2110.05026 |
| 1725 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 1760 Hz (A6) | Working note in Chladni particle manipulation | Mainstream science | Nat. Commun. 12764 |
| 1900 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| ~2 kHz | Phononic-crystal pass band | Mainstream science | PMC8226526 |
| 2–16 kHz | Usable band of a 30°-half-angle PA conical horn (60° coverage) | Mainstream science | Maker Portal |
| 2105 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 2260 Hz | PASCO plate resonance | Mainstream science | OU junior lab |
| 2700 Hz | PASCO plate resonance; circular plate n = 5 mode | Mainstream science | OU junior lab |
| 3 kHz | Demonstrated 3D acoustic cloak operating frequency | Mainstream science | PMC8226526 |
| 3000–6000 Hz | Tunable acoustic cloak operating range | Mainstream science | PMC8226526 |
| ~4 kHz | Phononic-crystal band gap; alternative metalens design frequency; circular plate n = 6 | Mainstream science | PMC8226526 / J. Eng. Appl. Technol. / OU junior lab |
| 5588 Hz (F8) | Working note in Chladni particle manipulation | Mainstream science | Nat. Commun. 12764 |
| 5750 Hz | Circular plate m = 7 mode | Mainstream science | OU junior lab |
| 10548 Hz (E9) | Upper bound of a working note range in Chladni manipulation | Mainstream science | Nat. Commun. 12764 |
| ~11,000 Hz | Sand "dance" regime in Lauterwasser documentary | Cymatics demo | ARTE clip |
| 11.9 / 19.6 / 23.9 / 27.1 / 34.5 / 43.5 kHz | Measured plate resonance peaks in Chladni manipulation platform | Mainstream science | Nat. Commun. 12764 |
| 12.5 kHz | Anti-alias filter corner in acoustic-spanner drive electronics | Mainstream science | Sci. Rep. 9:43046 |
| 20 kHz | λ = 17.25 mm; sets the largest practical acoustic meta-atom | Mainstream science | PMC8226526 |
| 29.83 kHz (A#10) | Highest note used in Chladni manipulation; λ = 16 mm, 4 mm node spacing | Mainstream science | Nat. Commun. 12764 |
| 40 kHz | TinyLev acoustic levitator operating frequency; λ = 8.65 mm at 25 °C | Mainstream science | Rev. Sci. Instrum. 88:085105 |
| 200 Hz – 50 kHz | Chirp used to characterize manipulation plate resonances | Mainstream science | Nat. Commun. 12764 |
| 135 kHz – 10 MHz | Meyl "Experimental-Kit" scalar-wave generator range | Esoteric-occult | Rex Research / Meyl |
| 328 kHz | PWM switching frequency in acoustic-spanner Arduino Due drive | Mainstream science | Sci. Rep. 9:43046 |
| 1200 kHz (1.2 MHz) | Square-wave drive used in maker Rodin-coil overunity test | Esoteric-occult | Rex Research Rodin tests |
| 4–8 MHz | Meyl "Demo-Kit" scalar-wave generator range | Esoteric-occult | Rex Research / Meyl |

---

## Key Angles and Geometric Constants

(Companion reference for diagrams)

| Value | Meaning | Category | Source |
|-------|---------|----------|--------|
| 51°51'51" | Vogel crystal large/female end cut = Great Pyramid of Cheops angle | Esoteric-occult | MarcelVogel.pdf |
| 52° | Quartz interlattice angle = Vogel female-end cut (90° − 38° basal) | Esoteric-occult | marcelvogel.org |
| ~60° | Vogel male-end cut angle ("resonant with the individual") | Esoteric-occult | marcelvogel.org |
| 104.5° | Water intermolecular bond angle | Esoteric-occult framing of a standard value | marcelvogel.org |
| 109.5° | Tetrahedral angle | Esoteric-occult framing | marcelvogel.org |
| 120° / 60° / 51.43° | Hexagon interior / 6-ring phase angle (360/6) / 7-ring phase angle (360/7) | Esoteric-occult | marcelvogel.org |
| 1.618 (Φ), 1.732 (√3) | Golden ratio and hexagonal ratio invoked in Vogel cut theory | Esoteric-occult | marcelvogel.org |
| 150° winding angle, 36 poles, 30° pole separation | Measured Rodin coil geometry (dodecagon, OD 9.95 cm, aperture 2.9 cm) | Mainstream science | Rev. Mex. Fís. 70:030901 |
| 24-sided polygon; 6 self-referencing axes | Rodin torus geometry per proponent literature | Esoteric-occult | Rodin Aerodynamics |
| 13 circles, 78 lines, 5 Platonic solids | Metatron's Cube construction | Esoteric-occult | Metatron Collective |
| 30° inner-path rotation, 18° twist | Duke tapered labyrinthine metamaterial cell geometry | Mainstream science | Xie & Cummer, APL 103:201906 |
| 5–85° half-angle, 10–100 mm throat, 10–500 mm length | Printable horn/waveguide parameter envelope | Mainstream science | Maker Portal |

---

## Notes on Coverage and Gaps

- **Vogel's radionic rate "454":** No source located gives a Hz value; the source explicitly states no unit. Recorded as unitless.

- **Rodin Aerodynamics:** The proponent document states no frequencies in Hz and no mm/cm dimensions — only that the coil is "much more sensitive to a particular frequency of applied current". Recorded as n.a. for those fields.

- **Maglione & Mazzocchi orgonomy review:** States no frequencies in Hz; its quantitative content is temperature, distance, humidity, and electroscope-discharge ratios.

- **Emoto's methodology:** Specifies exposure durations and scoring scales but no Hz values for water crystal formation; the frequency associations in Section 4.5 come from the Solfeggio and 432 Hz literatures, not from Emoto directly.
