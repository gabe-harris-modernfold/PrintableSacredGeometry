# Language of the Plume

A reading table for incense smoke, and a search for a sign language the smoke can carry. You sit
at a round table on a ridge top with a burning cone at its centre. A communicator sends signs
into the plume, and the ridge's wind tries to wreck them. The page learns how *you* read the smoke
and which signals you favour. Then it evolves a set of signs whose reading lands on what was meant.

The page is [`web/language-of-the-plume.html`](../../web/language-of-the-plume.html). It is a
single self-contained HTML file with no build step. This folder holds the tools that test and
picture it.

## What is on the page

- **Read.** One unbroken sitting with no start, clip or submit. Now and then a sign drifts up in
  the smoke. Touch what you saw, or touch the table where a Stone leaned. A touch while the sign
  shows, or up to 12 s after it, counts as your reading of that sign. A touch in a quiet stretch
  counts as a reading of the plain wind. A pick is kept after 9 s left alone, and nothing on screen
  changes when a sign starts or ends.
- **Blind.** You are never told what was sent. Your answers record how you see; they are not
  corrected toward a key. They are saved in the artifact's database collection `trials` (or in
  this browser when the page is opened elsewhere).
- **Evolve.** A genetic search over burn time, dial size, forms, hold length, beat and opening
  count. Each candidate language is scored against the reader model fitted from your answers.
- **Free.** Send any sign by hand and play with the time of day, gusts and burn.

## The model in one paragraph

The smoke is a set of particles riding an Ornstein–Uhlenbeck ridge wind with divergence-free eddies.
It is drawn as threads: each strand is one smooth curve whose darkness is mass over area. The smoke
stays smooth for only a few centimetres above the cone, then breaks up. A tall straight column is
rare: it needs a calm spell at dusk (about 2–5% of the time) or the communicator's grip (Still,
Fork, Tower). The **dominant air arc** on the table's rim is a slowly forgetting tally of the
natural wind's bearing alone. It never includes the communicator's push, so it cannot give a sign
away. It shows where a Stone would land most often if nobody steered it.

The reader model fits from your answers:
- a form-confusion matrix, by EM, which separates "the air wrecked it" from "it was misread";
- the sharpness of your direction reading, as a likelihood over the physical lean samples;
- your puff-count noise;
- how many forms you can hold at once. The communicator's vocabulary changes from chapter to
  chapter, so the option count varies, and that variation is what makes this limit fittable.

## Files

| File | What it does |
|---|---|
| `core.py` | Copies sections 1–6 of the page (the model, with no DOM in it) into `plume_core.js` so Node can run exactly the page's code. It generates that file, so never edit `plume_core.js` by hand. |
| `check_plume.js` | Checks the behaviour the page must keep: a short smooth column, a rare tall thread, Tower still reaching about 55 cm, no false two-streams, a Fork read as two streams, a Stone landing on its hour, recovery of the stand-in reader, evolution beating the first guess, and step cost. Exit code = number of failures. Takes about 2 min. |
| `shots.py` | Headless-Edge preview renders of named scenes into `previews/`. |
| `previews/` | `sitting`, `picked`, `air_day`, `air_night`, `column`, `tower`, `fork`, `ring`, `light`. |

```bash
cd projects/language_of_the_plume
python core.py && node check_plume.js
python shots.py            # or: python shots.py fork air_day
```

The page also has a self-test in the browser: open it with `#selftest` to append the results.
`#evolve` and `#free` open those modes directly.

## Things that bite

- **Headless Edge hangs** if it shares a profile with an open Edge window. `shots.py` gives it a
  throwaway `--user-data-dir`.
- **Don't add per-particle random walk to the smoke.** It zigzags the threads. Spread is carried as
  a width (`w2`), and the stream detector counts only smooth smoke (`tb < 0.8`).
- **Saved readings use the localStorage key `smoke-speech-lab.trials.v1`**, from the page's first
  name. Keep it, or earlier readings stop loading. Trial records from v1 are renumbered on load.
  Records from v3 on carry `flow: true`, with `rt` measured from the start of the sign.
- **Not yet on the site or in the menu.** `deploy.ps1` copies every `web/*.html` page to `site/`,
  so the next deploy will publish it. Add a card to `site/index.html` when it should be listed.
