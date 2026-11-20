# The Shadow Solids

A modern table of magical correspondences for psychological occultism and shadow work. The five
Platonic solids, simplest to most complex, each hold one Jungian archetype with two faces:

| Solid | Integrated | Shadow | Quality |
|---|---|---|---|
| Tetrahedron | The Child | The Eternal Child (*puer aeternus*) | Ignition |
| Hexahedron (Cube) | The Persona | The Mask | Stability |
| Octahedron | The Inner Other (anima / animus) | The Projected Other | Reciprocity |
| Icosahedron | The Sage (mana-personality) | The Inflated Guru | Flow |
| Dodecahedron | The Self | The Divided Self | Wholeness |

For each, the table gives integrated and shadow traits (with the defence mechanisms named), core drive
and fear, a somatic cue for both states, a balanced and a shadow hex colour with the shift between them,
a five-sense anchor plus a "shadow tell" for each sense, one or two cultural objects from 2006–2026, and
a transmutation practice with a sigil.

**Excluded by design:** Hebrew letters, Kabbalah, Golden Dawn frameworks, traditional planetary and
elemental attributions, classical god names and Renaissance correspondences (planetary metals, classical
incense). Each solid's quality is read from its geometry alone: faces, corners, symmetry and dual.

## Files

- The page is [`web/the-shadow-solids.html`](../../web/the-shadow-solids.html). It is one self-contained
  HTML file with no build step. All of its content sits in a single JSON block (`<script id="ss-data">`)
  that drives the hero, the master table, the five rooms, the closing diagram and the Markdown export.
- [`check_solids.py`](check_solids.py) validates that block and writes the Markdown reference. Run it
  from this directory after editing the page:

  ```
  python check_solids.py
  ```

  It recomputes vertices, edges, faces, duals, circumsphere fill and sphericity with SciPy's convex hull
  and compares them with the page. It checks the quantitative colour-shift sentences against the hex values
  in HSL, scans every string on the page for the excluded vocabulary (with a short allow-list for literal
  modern uses such as "sea air"), confirms each cultural object falls inside 2006–2026, and enforces the
  brief's structure: 2–3 traits, five senses, 1–2 objects, a transmutation of 1–2 sentences.
- [`the-shadow-solids.md`](the-shadow-solids.md) is the generated reference: the master Markdown table,
  the sensory and ritual breakdown for each solid, and the closing note. Do not edit it by hand.

## What is on the page

1. **The five.** Five wireframe solids drawn live from their exact coordinates, with one switch that turns
   all five to their shadow state. Each shadow has its own motion: the tetrahedron races and wobbles, the
   cube freezes and thickens, the octahedron throws a ghost of itself in the complementary colour (the
   projection as an afterimage), the icosahedron swells, lifts and glows cold, and the dodecahedron's
   edges fly apart and drop out.
2. **The master table.** All ten columns of the brief, with the solid pinned while the table scrolls
   sideways. *Copy table as Markdown* puts the master table on the clipboard; *Download the full
   reference* saves the whole document as `.md`.
3. **Five rooms**, one per solid. A slider runs from integrated to shadow and moves the solid, the live
   colour readout and the emphasis of the text. *Hear the anchor* synthesises the sound anchor with Web
   Audio and crossfades it into the shadow tell as the slider moves:
   - Tetrahedron: a square-wave major-third arpeggio at 128 BPM, then a loop of notification pings.
   - Cube: a sine perfect fifth with a click at 60 BPM, then office mains hum and buzz.
   - Octahedron: two sawtooth voices in unison under a 66 BPM pulse, pulled 3 Hz apart so they beat.
   - Icosahedron: brown-noise rain under a triangle-wave perfect fourth, then a pressured hush and a thin
     high whine.
   - Dodecahedron: a sine drone on the harmonic series of 55 Hz, then the same drone stuttering with
     dropouts.
4. **How the geometry mirrors the growth.** Edges climb 6, 12, 12, 30, 30. One self-dual seed is followed
   by two dual pairs, and this matches Jung's 1928 sequence (persona, anima/animus, mana-personality, self)
   with the Child (1940) as the seed. Faces change kind only twice, and the pentagons carry φ. The
   icosahedron is the most spherical but the dodecahedron fills the most of its sphere.

Non-English terms (persona, anima, animus, puer aeternus, mana) take a tap-to-translate underline, and
there is a *Show all translations* toggle. The care note says plainly that this is a reflective practice,
not therapy.

## Printing

Nothing here is a new print. Each room links to the printable solids in
[From One Circle](../../web/from-one-circle.html#solids) as objects to hold during the practice.
