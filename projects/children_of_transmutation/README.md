# The Children of Transmutation

An infographic catalogue and timeline of what humans have turned into what. It covers energy changed
from one form to another, elements changed into other elements, the laboratory alchemy that sat between
them, and the inner transmutations of the spiritual traditions. The three kinds of change make different
claims, so every entry carries a label saying which kind it makes:

- **Demonstrated**: shown by instrument, record or excavation, and repeatable.
- **Measured in practitioners**: an effect on the body has been measured; what it means is still open.
- **Practised and taught**: a real discipline, described in its own terms.
- **Reported**: testimony, a vision, or a later account.
- **Claimed**: asserted but never shown.
- **Refuted**: tested and failed, or exposed as fraud.

The page is [`web/the-children-of-transmutation.html`](../../web/the-children-of-transmutation.html).
It is one self-contained HTML file with no build step. Nothing here is printable; this folder holds the
tooling that checks the page.

## What is on the page

0. **The emblem**, *The Transformer*, beside the opening. The source is
   [`the-transformer.png`](the-transformer.png) in this folder; the page uses a 700 px WebP at
   `site/assets/transmutation/the-transformer.webp`. Its crow, dove and phoenix are hotspots that open
   the nigredo, albedo and rubedo rungs of the colours ladder, and the ouroboros opens the
   "Nothing lost, all transformed" card. The image loads only from `site/`, where `assets/` lives.

   Five more engraved plates explain the hardest concepts. They were generated from
   [image-prompts.md](image-prompts.md); the sources are the `illus-*.png` files here and the WebPs
   are in `site/assets/transmutation/`.
   - *The two roads* (IV): a single step versus a chain through heat.
   - *The tithe of heat* (V): the losses at each station of the coal chain.
   - *The near miss* (VI): lead into gold at the LHC.
   - *The valley of iron* (VI): fusion and fission both falling toward iron.
   - *One heat, two readings* (VII): the measured effect versus the tradition's reading of it.

   Every caption links its parts to the matching chart, and clicking any plate enlarges it.
1. **Born of transmutation**: a ring chart of the body's mass by where its atoms were made. About 76%
   came from exploding massive stars, 14% from dying low-mass stars and 10% from Big Bang hydrogen,
   with traces from white dwarfs and neutron-star mergers. Composition follows ICRP Reference Man and
   Emsley; stellar origins follow Johnson (2019).
2. **Three kinds of change, and the Art between them**, with the six-label legend.
3. **The chronicle**: 163 dated entries in four lanes (energy, matter, the Art, the inner work).
   - The default view is a log scale of years before present, running from the first fire people kept
     (about 1 million years ago) to lead turned into gold at the LHC (2025).
   - Three zoomed linear views cover 3500 BC–1500, 1500–1900 and 1900 to now.
   - Lanes and labels can be filtered.
   - A detail card and a scrolling list show every entry.
4. **The octagon of forms**: 8 forms of energy and 37 conversions over 30 directed pairs.
   - Solid arrows are single physical steps.
   - Dashed arrows are chains (fuel to motion by an engine, current to light by a filament). Each
     dashed arrow is drawn so that its midpoint falls exactly on the form it passes through, Heat:
     the control point is C = 2V − (A+B)/2.
   - The table under the octagon lists every conversion.
5. **The inheritance of a joule**: five chains from source to lamp (coal, sunlight, uranium, falling
   water, bread), drawn as one tapering ribbon. The heat lost at each step peels away below it, and
   you can switch the lamp between an LED and a filament bulb.
6. **The ladder of elements**: proton number plotted against year.
   - Dots: the 29 elements people made.
   - Arrows: 15 landmark transmutations. Gold arrows end at gold.
   - Reference lines: the iron peak, lead and uranium.
7. **Ladders of refinement**: twelve traditions' maps of inner transmutation, each stretched to the
   same height.
   - Non-English rung names (32 rungs, plus two column headers) are underlined with dots. Tapping one
     shows its English; *Show all translations* reveals them all.
   - Each tradition gets a practitioner's voice and a physicist's voice.
   - Under the ladders is a strip of the practices whose bodily effects have been measured.
8. **Where the lanes touch**: twelve correspondence cards, each labelled as one of four kinds of link:
   inheritance (history), measured bridge, rhyme (resonance, not evidence) or contrast. Each card says
   what recurs and what the change of scale alters.

## Data and checks

All content lives in JSON blocks at the bottom of the page: `d-sources`, `d-body`, `d-events`,
`d-forms`, `d-conv`, `d-chains`, `d-elements`, `d-ladders` and `d-touch`. Edit those blocks, then run:

```bash
cd projects/children_of_transmutation && python check_data.py
```

The check confirms that ids are unique and that every cross-reference resolves (sources, chronicle
events, octagon pairs, element arrows, ladders, cards, chains, page links). It also checks that lanes
and labels are known values, that the body origins sum to 100, and that each chain has one start, sane
efficiencies and a lamp at the end. It exits non-zero on any failure.

One gotcha: in `d-elements.synth`, `s` is the element *symbol*. That element's sources go in `src`.

## Fact-check notes (2026-11-05)

An independent web fact-check corrected six points:

- Kozhevnikov's 38.3 °C is an armpit temperature, reached with breathing and visualisation together.
- The body composition now uses ICRP and Emsley's figures.
- ALICE's 86 billion gold nuclei is an estimate for the whole LHC, worked out from ALICE's measurement.
- Six of the seven Sufi nafs names have Qur'anic roots.
- Volta's letter was dated March 1800 and read to the Royal Society in June.
- Holborn Viaduct and Pearl Street were the first *coal-fired* central stations.

The Uncompahgre Ute quartz-rattle story could not be traced to an ethnographic source and is worded as
"said to". The Oklo reactors are dated as 1.7–2 billion years ago.
