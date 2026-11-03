# Language of the Plume

A reading table for incense smoke, and a trainer for reading it. You sit at a round table on a
ridge top with a burning cone at its centre and watch. Now and then something is sent into the
plume; the rest is the ridge's own air. The page learns how *you* read the smoke, above all when
you feel a presence in it and when you feel only weather, and tells you, sitting by sitting, how
that reading grows.

The page is [`web/language-of-the-plume.html`](../../web/language-of-the-plume.html). It is a
single self-contained HTML file with no build step. This folder holds the tools that test and
picture it.

## What is on the page

- **Read.** One unbroken sitting with no start, clip or submit. Now and then a sign drifts up in
  the smoke. Touch what you saw, or touch the table where a Stone leaned, and say how it felt:
  sent, or like the wind. A touch while the sign shows, or up to 12 s after it, counts as your
  reading of that sign. A touch in a quiet stretch counts as a reading of the plain wind. A pick is
  kept after 9 s left alone, and nothing on screen changes when a sign starts or ends.
- **Blind in the moment, told between sittings.** You are never told what a given moment was. When
  a sitting ends (every 4–7 signs, as the hour and the cone move on) the dock tells the sitting as
  a whole: how many sent moments felt sent, how many quiet stretches did, how many forms were named
  right, and the one confusion that recurred. The presence score and the fitted reader model move
  only then. Readings are saved in the artifact's database collection `trials` (or in this browser
  when the page is opened elsewhere).
- **Presence.** The headline: d′ and criterion from the felt ratings (sent = 4–5, wind = 1–2, and
  an unrated "Just wind" pick counts as felt wind) against the truth of each moment, with the
  log-linear correction. It needs six rated sent moments and six rated quiet stretches to show.
- **Codebook.** The sender's side of the table: a genetic search over burn time, dial size, forms,
  hold length, beat and opening count, scored against the reader model fitted from your answers.
  It shows what a sender would have to do to be read by you.
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
| `check_plume.js` | Checks the behaviour the page must keep: a short smooth column, a rare tall thread, Tower still reaching about 55 cm, no false two-streams, a Fork read as two streams, a Stone landing on its hour, recovery of the stand-in reader, the presence score on the stand-in and on a toy reader, evolution beating the first guess, and step cost. Exit code = number of failures. Takes about 2 min. |
| `shots.py` | Headless-Edge preview renders of named scenes into `previews/`. |
| `previews/` | `sitting`, `told` (a sitting's end: the card in the dock), `picked`, `air_day`, `air_night`, `column`, `tower`, `fork`, `ring`, `light`. |

```bash
cd projects/language_of_the_plume
python core.py && node check_plume.js
python shots.py            # or: python shots.py fork air_day
```

The page also has a self-test in the browser: open it with `#selftest` to append the results.
`#evolve` and `#free` open those modes directly.

## Where this is going

The page began as a channel-coding lab (`smoke-speech-lab.html`): a sender with a grip on the air,
the wind as the adversary, a code evolved for bits per cone. The intent is the inverse, a trainer
for reading a plume as the speech of something that lives in the wind. Phase 1 (the page as it is
now) made the felt sense the headline and moved every word the table says back to the end of a
sitting. Still to do:

- **Phase 2, the model.** The sender becomes latent structure in the wind process itself and
  `WORLD.k` (the grip) is retired; direction is read relative to the seat (toward, away, left,
  right) instead of twelve hours; puff-count coding and the prime opening go; difficulty adapts by
  staircase; the saved seeds let the same weather be replayed weeks apart for test-retest.
- **Phase 3, lexicon and lenses.** The notes cluster into a personal lexicon the reader can name
  and promote to words; the smoke-reading traditions become lenses (a sourced writeup in `docs/`
  comes first); the genetic search is retargeted from the sender's code to reading protocols,
  scored on this reader's consistency and discrimination.
- **Phase 4, Field mode.** The same table and dock with no simulation: readings of a real cone,
  then a camera tracker that yields the same instrument readouts.

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
