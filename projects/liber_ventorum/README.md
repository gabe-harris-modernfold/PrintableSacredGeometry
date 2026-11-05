# Liber Ventorum: a hanging octahedral wind harp

Research and design of a centre-radiating aeolian harp that hangs from a hook, steadied by a counterweight.
It is built from clear PETG tubes, steel music strings and Helmholtz chambers, with the octahedron as its
form. Started 2026-11-05. The work is still open-ended and pursuing design goals; this is not yet a build
spec.

**Resume here.** The main artifact is the research document. This folder holds the scripts, renders and
simulation results behind it. `notes/` holds the full conversation and everything the user pasted in.

## Where everything lives

| What | Where |
| --- | --- |
| Research document (Claude Doc, the main artifact) | [Liber Ventorum — Research Toward a Wind Harp](https://claude.ai/artifact/VuZTgwy8xGD7Ekpxursbcy), doc id `ea169013-8f61-42bc-b1c5-28ef216a318e` |
| Readable conversation, 2026-11-05 (through the request to save everything; run v4's result is recorded here and in the doc) | [notes/conversation_2026-11-05.md](notes/conversation_2026-11-05.md) |
| The user's pasted briefs and the Pd patch screenshot | [notes/sources/](notes/sources/) |
| Full session export (transcript, metadata; 21 MB, not in the repo) | `C:\Users\Harrisg1\Downloads\session-export-1793901871619.zip` |
| Four concept massings (A–D) | `concepts.py` → `concept_a.png` … `concept_d.png`, `concepts_lineup.png` |
| Concept A developed, rev 1 | `plumb_octahedron.py` → `plumb_octahedron_rev1.png` |
| Evolutionary simulation of the built harp | `evolve_tuning.py`, `cf_study.py`, `evolution_*_v1…v4` |

The research document's sections, in order:
- purpose and legend;
- Aeolus, what is stated and what is inferred;
- the 13-string Pd model;
- Physics I and Physics II (beating);
- the magician's reading, and the magic of beating;
- rhymes between the two voices;
- theta and trance;
- emerging design goals (Status dropdowns);
- design proposals A–D;
- Concept A developed (rev 1);
- the evolutionary simulation, including carbon fibre and tunability;
- open questions, kin in this repository, and sources.

## How the user works

The user is a self-described poetic mystic. Every finding is read in two voices.
- **Physicist:** marks each claim established (●), inferred (◐) or speculative (○).
- **Magician:** reads the same facts through correspondence, scale invariance, geometry and archetypal pairs
  (vessel and current, threshold and passage, seed and tree, centre and circumference, descent and return).

Resonance is kept distinct from proof. The full frame is in
[notes/sources/brief_poetic_mystic_frame.md](notes/sources/brief_poetic_mystic_frame.md).

## Decisions so far (2026-11-05)

1. **Two-voice reading.** Jerram's *Aeolus* was read for every correspondence to the wind; no design yet.
2. **Pd patch decoded.** The user's screenshot is a 13-string aeolian harp model: 0.020 in strings, 24 in long,
   with Mass_Den as a volume density. Its "harmonic 3 for all strings" only reproduces with a
   Reynolds-dependent Strouhal law (Williamson's fit).
3. **Research document created.** Acoustic beating and the two humming briefs were added: the 1–5 Hz felt
   throb, bone conduction, hum gently, 85 dBA caution.
4. **Theta and trance.** The harp can carry theta-rate (4–8 Hz) beats; evidence that this entrains the brain
   is weak. Hard cautions: no flicker, gentle levels.
5. **Physical brief.** The user set it: centre-radiating, hook, counterweight, clear PETG tubes, steel strings,
   Helmholtz chambers, octahedron. Four concepts followed: A Plumb Octahedron, B Axis Monochord,
   C Tonos Tensegrity, D Star of Valihas.
6. **User chose A.** Rev 1: 12 strings on the edges of an octahedron.
   - Earth strings D3 F3 A3 C4 are tensioned by an 11.8 kg bob.
   - Heaven strings E3 G3 B3 D4 are tensioned by the bob plus trim; together they spell D Dorian.
   - Waist strings run on springs.
   - The hub is a Helmholtz chamber on D3.
7. **Evolutionary simulation requested.** The user asked to make the tuning and theta "as close as possible
   to the design implementation".
   - **Run v1** gamed its score (strings on different harmonics, chord coherence 4 %) and was rejected.
   - **Run v2** added a coherence term: about 6 % better than rev 1 with springs on fresh builds.
8. **Carbon fibre allowed.** Carbon rods inside the clear tubes take the compression. Rev 1's theta rose from
   41 % to 62 %.
9. **Tunability required.** The tests:
   - transpose the whole harp ±2 semitones;
   - move one string ±1 semitone on its own adjuster;
   - at most 40 cents per turn.

   Springs fail the one-string test because the strings are coupled. Steelyard levers at the poles pass it.
   **Run v3** chose carbon cores and steelyards.
10. **Run v4** weighted chord intervals three times as much. Intervals improved from 20 to 9.7 cents and theta
    rose to 69 %, but chord coherence fell from 44 % to 26 %. This is the physical trade-off; choosing between
    v3 and v4 is the user's call.

## The two candidate designs (checked on 48 fresh simulated builds)

Both use carbon cores and steelyards, and both pass the tunability tests.

| | Run v3: the chord holds | Run v4: intervals first |
| --- | --- | --- |
| Interval error at sounding harmonics | 20.2 cents | 9.7 cents |
| Chord on one harmonic | 44 % | 26 % |
| Sounding time with a theta beat | 62 % | 69 % |
| Drift / gust bend | 0.2 / 5 cents | 0.2 / 3 cents |
| Octahedron, node to node | 506 mm | 548 mm |
| Earth strings D F A C | .018 .024 .009 .015 | .016 .018 .014 .011 |
| Heaven strings E G B D | .026 .020 .015 .011 | .019 .018 .014 .012 |
| Waist springs, offset | 7.9 N/mm, +1.2 Hz | 1.9 N/mm, +1.2 Hz |
| Set at, bob, trim | 12 °C, 10.7 kg, 0.59 kg | 2 °C, 9.4 kg, 0.73 kg |
| Hanging mass | 16.9 kg | 16.0 kg |

Full numbers: `evolution_result_v3.json`, `evolution_result_v4.json`, and the research doc's "third run" and
"fourth run" sections.

## Run v3 in detail

| Choice | Value |
| --- | --- |
| Frame | Clear PETG tubes with carbon rods inside; octahedron 506 mm node to node (about 451 mm sounding) |
| Gravity tensioners | Steelyards, 8:1. Heaven levers at the top node (poises hang from the hook); earth levers at the bottom node (poises are part of the bob) |
| Earth strings D F A C | .018 .024 .009 .015 in, at 21.5, 62.1, 14.6 and 51.4 N |
| Heaven strings E G B D | .026 .020 .015 .011 in, at 56.9, 45.0, 51.0 and 44.3 N |
| Waist strings | Earth gauges, 7.9 N/mm springs, 1.2 Hz above their partners |
| Install temperature, bob, trim | 12 °C, 10.7 kg, 0.59 kg per tip |
| Results | Drift 0.2 cents; theta in 62 % of the sounding time; chord on one harmonic 44 %; gust bend 5 cents; hub drone 21 cents from the root; **interval error 20 cents** (rev 1: 13); 16.9 kg hanging, including 2.5 kg of poises at the hook |

## Key findings to carry forward

- **A four-string star meeting at one node is statically indeterminate.** Without compliance, 0.1 mm of anchor
  travel moves a string about 100 cents, so the harp is untunable. Springs (about 2 N/mm) fix tuning; steelyards
  also decouple the strings.
- **Gravity tuning holds.** Weight-set tension ignores heat and creep. Drift is about 1.3 cents with springs
  and 0.1–0.2 cents with steelyards.
- **The chord holds less often than the equal-tension theory says** (26 % for rev 1). Thin strings shed at a
  lower Strouhal number at low Reynolds numbers.
- **Coherence and interval accuracy pull against each other.** Coherence wants thinner strings for higher
  notes; equal string stiffness (just intervals at high harmonics) wants thicker ones.
- **Theta tops out at roughly 50–65 % of the sounding time.** Carbon helps by stopping the waist strings from
  drifting with temperature, but it removes their sun voice.
- **The hub drone cannot track the strings.** Air in the chamber shifts about 3 cents per °C; the drone sits
  20–25 cents from the root on average unless the neck compensates.
- **Model limits.** Lock-in is a 6 % band; there is no string coupling through the nodes, no turbulence and no
  hub FEM; frame masses are estimates; the wind is Weibull with a 3 m/s mean.

## Open items

- [ ] Choose between v3 (the chord holds) and v4 (intervals first), or map the trade-off between them with a few
      structured gauge sets.
- [ ] Tune the hub separately at the site's mean temperature; the model currently ties it to the strings' setting
      temperature.
- [ ] Design the steelyard levers: ratio, graduated arm, poise, pivots, how they sit at the poles.
- [ ] CAD the printed parts within the 320 mm bed and check watertightness: hub halves, tip and pole nodes,
      bob shell, tube end fittings for the carbon rods.
- [ ] Check the hub pitch with FEM (the Ogdoad lesson); decide whistle mouth or humming port; consider a
      temperature-compensating neck.
- [ ] Optionally keep one PETG-only "sun pair" of waist strings, so the sun stays audible.
- [ ] Measure the site wind at the hanging height; get the Pd patch's `.pd` file.
- [ ] Bench tests: string on a skin (centre vs off-centre warble), sand figures, sunlit vs shaded tubes.

## How to rerun

Run from this folder. Windows multiprocessing is blocked in the agent sandbox, so the searches run serially:
about 0.2–0.6 s per design, 20–50 minutes per run.

```bash
python concepts.py                       # four concept massings + lineup
python plumb_octahedron.py               # rev-1 numbers + labelled render
python evolve_tuning.py --check          # statics model reproduces rev 1
python evolve_tuning.py --tag v3 --popsize 6 --maxiter 35
python evolve_tuning.py --tag v4 --popsize 6 --maxiter 35 --w-interval 1.5 --seed-from evolution_result_v3.json
python cf_study.py                       # PETG vs carbon comparison
```

Constraints from the repo: PETG, 320 × 320 × 320 mm bed, segmented parts; the matched sacred-geometry set
is sized from `lib/sacred_grid` (R = 15 mm).
