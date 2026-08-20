# Terrarium — the idea

**One sentence:** the terrarium is an **instrument for making ambient gradients visible**,
in which geometry is the only actuator — it never pushes, it sets boundary conditions and
waits, and the sacred geometry is what the physics *produces*, not what we engrave.

**Material is fixed: clear PETG.** That is not a finish note; it moves the physics three
ways: (1) the shell is only *clear* where the wall is one extrusion wide, so "thin wall,
mostly viewing ports" is **enforced by the material**, not chosen; (2) clear PETG passes
sunlight but is opaque across the 8–13 µm thermal window — the enclosure is a
**greenhouse**, which is where the driving ΔT comes from; (3) at θ ≈ 70° it fogs with
dropwise condensation, which is either the death of visibility or a harvest, depending on
the dome's shape. All three are handled below.

---

## 1. The reframe

A subtle force is one whose energy is already in the room — a fraction of a kelvin, a
surface-tension gradient, sunlight through a wall. It cannot be harnessed by making it
bigger. It is expressed by putting it somewhere it *has no choice but to become a pattern*.

Geometry has exactly three ways to do that. Everything below is one of these three.

| Amplifier | What geometry does | Gain |
|---|---|---|
| **Threshold** | Poise the system at a bifurcation. An imperceptible input crosses it and produces a **discrete visible event**. | ∞ at the point of instability |
| **Accumulation** | Collect a nanolitre-per-second over hours into one bead that falls. | ~10⁴–10⁶ in time |
| **Selection** | A dimension picks one mode out of broadband noise. Depth sets cell size; slot width sets rise height. | Turns noise into a lattice |

## 2. The subtle forces actually available, with honest numbers

### The greenhouse ΔT — the prime mover

Clear PETG transmits solar/visible in and blocks thermal IR out. In any sunlit room the
glazed headspace runs warm while the wet moss bed sits near its depressed wet-bulb
temperature. Net: a **3–8 K internal gradient, for free, every day**, dying at night and
returning with the sun — which is the correct behaviour for the thesis, not a fault.

(Radiative sky cooling — a fin seeing the night sky at 3–8 K below ambient — is real
physics but **does not work through PETG or glass**: both are opaque exactly where the sky
is cold. It is an outdoor option only, demoted to §8.)

### Bénard–Marangoni — hexagons that build themselves

A shallow water layer with a vertical ΔT tiles itself into hexagonal cells above
Marangoni number ≈ 80. For a **2 mm** layer:

```
Ma ≈ 2100 · ΔT[K]   →   supercritical at ΔT ≈ 0.04 K
```

Evaporative cooling of the free surface alone supplies far more than 0.04 K, so the honest
statement is stronger than the number: **the dish cannot not convect.** The only failure
mode is chemical — organic surfactants flatten Marangoni flow dead. So the dish must run
on **distilled condensate**, which the shell produces for free (see fog, below).

- Cell pitch ≈ 2–3× depth → **4–6 mm hexagons** in a 2 mm layer. Flower-of-Life scale,
  convected into existence because the dish is 2 mm deep.
- The depth is a mode selector, so the depth must be **held by geometry**: an overflow
  weir locks the layer at 2 mm regardless of drip rate. The constant is kept by a shape,
  not by attention.
- Made visible by settled spore dust (lycopodium); cells creep at ~0.1 mm/s.

### Rayleigh–Bénard in air — the chamber picks its own convection cell

```
Ra = gβΔT d³/(να);   air, d = 50 mm, ΔT = 1 K  →  Ra ≈ 12,400   (critical 1708)
```

A 50 mm headspace at one degree is 7× supercritical — it convects on its own, in cells
about as wide as the gap is tall. Change the chamber height and you change the pattern.

### Capillarity — water climbing on shape alone, derated for PETG

Jurin's law with the **real contact angle**. PETG is not glass: θ ≈ 70°, cos θ ≈ 0.34.

```
h = 2σcosθ/(ρgw)      w = 0.4 mm (one nozzle width)  →  h ≈ 13 mm
                      w = 0.15 mm (tuned gap)        →  h ≈ 33 mm
```

(The often-quoted 37 mm at 0.4 mm assumes a perfectly wetting wall.) Two recoveries, both
legitimate: **layer grooves** aligned with the rise act as sub-slot capillaries and lift
higher than the nominal gap predicts — measure, don't assume; and a **flame/plasma pass**
takes PETG to θ ≈ 20°, restoring ~34 mm at 0.4 mm. Either way the slot width is a design
variable that sets lift height, and "watering pipes" become passive capillary slots.

### Conical transport — the cactus-spine effect

A droplet on a cone of half-angle 2–5° is pushed by its own Laplace-pressure asymmetry
toward the **wider** end, unaided, at ~mm/s. Consequence for us: spines that deliver
condensate *downward* must **widen downward** — a print-orientation constraint, and a
happy one (fat end on the bed, needle tip up).

### Fog — the enemy converted

Dropwise condensation on θ ≈ 70° walls scatters light and blinds every viewing port. The
geometric answer is to make the dome **crowned with meridian gutters**: droplets grow,
merge, hit run-off size, and are steered off the view panels into drip edges that feed the
spines that feed the dish. **Anti-fog and distillate harvest are the same part.** This is
the function-stacking requirement earned, not asserted.

### Rayleigh–Plateau — where drops decide their own spacing

A film on a fibre beads up at λ ≈ 8.89 × fibre radius (the free-jet constant is 9.02 R —
same instability, different boundary). A wettable filament across a chamber turns dew into
a rhythm nobody chose.

### Plateau's laws — 120°, always

Three films meet at exactly 120°; four edges at cos⁻¹(−⅓) ≈ 109.47°. Any wet junction in
this object will insist on it.

## 3. The stack

```
        sun ──► clear PETG shell (greenhouse: solar in, thermal IR trapped)
                    │
        headspace warm · wet bed cool  ──►  the working ΔT (3–8 K, daily)
                    │
        fog beads on the dome ──► meridian gutters ──► drip edges
                    │  (distilled, surfactant-free)
                    ▼
        spines (widening downward) walk each droplet to
                    ▼
        the Marangoni dish — 2 mm deep, weir-locked ──► HEXAGONS
                    │  overflow drips at a threshold
                    ▼
        bead-fibre span (λ ≈ 8.9 r) ──► dew spaces itself
                    ▼
        capillary slots (w = 0.4 mm → ~13 mm lift) water the bed
                    ▼
        living bed · 50 mm air gap above ──► Ra ≈ 12,400: it convects
```

One gradient, made by the shell itself, expressing as five different lattices — each with
a geometric constant attached: 0.04 K, 4–6 mm, 1708, 8.89, 120°, 13 mm.

## 4. The water circuit — the siphons stay

Explicit decision, because [BRIEF.md](BRIEF.md) lists the bell siphon as worth keeping
whatever the form: **the siphons remain the water circuit.** The solar pump (settled,
~0.5 W) lifts to the top chamber; three cells flood and drain through bell siphons at
volumes **1:2:3**, beating on a common period of 6T. That is the *loud clock* — the
multichamber water levels, water flowing through, of the original ask.

The subtle layer of §3 **rides on top of it** and runs on nothing: the siphons move the
bulk water; the instabilities pattern the margins. Two registers, one object — the audible
pulse and the almost-invisible lattice, and the demonstration is that both obey the same
accumulate → threshold → discharge → reset law: drip (~s), siphon (~10 min, 1:2:3),
convective overturn (~min), solar day (24 h).

What stays demoted is the siphon-as-*bellows* (Tesla-valve air rectification, Re ≈ 240 in
a 6 mm duct). Sound physics, wrong register — it belongs in an object about pumping.

## 5. Requirements map (against BRIEF.md, line by line)

| Requirement | Where it lands |
|---|---|
| Multichamber water levels, water flowing | Bell siphons, 1:2:3, pump-fed (§4) |
| Watering pipes | Capillary slots — passive, width sets lift (§2) |
| Thin wall, mostly viewing ports | Enforced by clear PETG: clarity *requires* single-extrusion walls |
| Organic tubes running vertically | The spine-and-fibre columns and capillary risers — they do the transport, they are not ridges |
| Trays with snap-fit domes | Each chamber = tray + clear crowned dome; the dome's crown *is* the fog-harvester, so the lift-off part earns its shape |
| Occult / hermetic through working parts | The constants: hexagons from a depth, 120° from films, spacing from a radius |
| Permaculture function stacking | Condensate = anti-fog + distillate for the dish + watering, one flow, three jobs |
| Vibration & correspondence | Same threshold law at four timescales, all stopwatch-able (§4) |
| ≤ 320 × 320 mm | Constraint stands; nothing above forces it larger — chamber count and dome crowns are the variables to fit it |

## 6. Why the form language stops being an open question

BRIEF.md says the form language "needs a reference, not another guess." The reference is:
**every dimension must be traceable to a physical constant.** Depth 2 mm because Ma > 80.
Gap 50 mm because Ra > 1708. Slot 0.4 mm because Jurin. Weir height because the depth is a
mode selector. Cone angle 3° because that walks a droplet. Wall one extrusion wide because
that is where PETG is clear. Nothing is styled.

## 7. What must be proven first

1. **Surfactant lifetime.** Bench a 2 mm distilled dish inside a humid, spore-laden box.
   Hours before contamination flattens the hexagons? If short, the dish needs sacrificial
   skimming geometry (weir overflow already helps — surface film exits over the weir).
2. **Real capillary rise on printed clear PETG.** θ plus layer-groove effects, measured,
   before any slot dimension is trusted.
3. **Greenhouse ΔT survey.** Log headspace vs bed temperature in the intended room over a
   sunny day and a cloudy one. If it never clears ~1 K, the Ra chamber shrinks or the
   claim softens.
4. **Does the dome actually shed?** Fog run-off size vs crown slope on printed clear PETG.
   If drops pin on layer lines and the ports stay blind, the geometry has failed at its
   first job: being seen through.
5. **Visibility.** Tracer strategy (lycopodium, mist, one dyed filament) decided *before*
   the geometry, not after.

## 8. Demoted, not deleted

- **Siphon-bellows + Tesla valves** — a real pump, wrong register here.
- **Radiative sky fin** — real cold, but blind behind PETG or glass; outdoor variant only,
  or under a polyethylene window if the terrarium ever lives outside.

## 9. Still excluded

Knudsen channels (need sub-700 nm; PETG floor is 0.4 mm — three orders short) and acoustic
streaming (60 dB gurgle → 5 × 10⁻⁵ m/s particle velocity, streaming scales as its square).
Correct physics, wrong regime. Claiming them would be decoration pretending to be
mechanism.
