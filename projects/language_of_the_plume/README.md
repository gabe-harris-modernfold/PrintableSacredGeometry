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
  the smoke. Touch what you saw, or touch the side of the table the smoke leaned to (toward you,
  away, left or right), and say how it felt: sent, or like the wind. A touch while the sign shows, or up to 12 s after it, counts as your
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
- **Air.** Seven levels, from a soft wind fully committed to its signs to a gusty one that only
  half means them (`LEVELS`: gust multiplier and hold). After a sitting in which the felt sense told
  sent from wind by a margin of 0.5 or more, the air moves up one; at 0.15 or less, down one. The
  level lives in this browser (`smoke-speech-lab.practice`).
- **Sittings that come round again.** Every sitting follows from its seed, so after a week an old
  sitting can return under the same weather and the same signs, unannounced, every fourth sitting.
  Agreement between the two readings of each sign is the reliability figure in "The reader".
- **Lexicon and lenses.** Every note you write gathers under the cue that drew it (the form, with
  its side for a Stone); a word that comes back in another note of the same cue is offered as your
  word for it, and you can keep one. A lens is a dictionary from what you saw to what a tradition
  made of it, chosen in the dock; after each kept reading the hint says what the lens made of it.
  It speaks of what you saw, never of what was sent, so it can speak at once. "Your own" is the
  lexicon itself. The traditions' dictionaries are lifted from
  [docs/smoke-reading-traditions.md](../../docs/smoke-reading-traditions.md), every rule with its
  source. The lens board shows, per lens, your presence score and how steady your reading was when
  a sitting came round again, once six readings have been kept under it.
- **The readers of smoke.** Every known way of reading smoke, in five families (omen, rite,
  catalogue, signal, science), as cards and on one timeline of the people who did it, from the Old
  Babylonian tablets to Morton, Taylor and Turner. Each item links its source; a card whose
  tradition left rules for what was seen offers itself as a lens. All of it is lifted from the
  writeup by `lift_traditions.py`, which refuses to run if any item cites a source missing from the
  table.
- **Protocol.** Your side of the table: a genetic search over which forms to attend to, whether to
  read the Stone's side, when to sit and which stretch of the cone to listen through, scored on
  the reader model fitted from your answers. The wind sends every form; a form let pass is read as
  weather. Scored as bits learned (mutual information) by default, or as clean marks.
- **Free.** Send any sign by hand and play with the time of day, gusts and burn.

## The model in one paragraph

The smoke is a set of particles riding an Ornstein–Uhlenbeck ridge wind with divergence-free eddies.
It is drawn as threads: each strand is one smooth curve whose darkness is mass over area. The smoke
stays smooth for only a few centimetres above the cone, then breaks up. A tall straight column is
rare: it needs a calm spell at dusk (about 2–5% of the time) or the wind holding its breath (Still,
Tower). A sign is a mode the wind enters, not a push laid over it (`Wind.drive`): for a few
seconds the mean turns toward a bearing on a short memory while its own habit fades by the form's
`QUIET`, the gusts die down, and when the sign ends the wind wanders back on its slow clock. The
**air arc** on the rim is a slowly forgetting tally of where the one wind, signs and all, has been
carrying the smoke: its habit. A lean that turns against the habit and holds is what a reader
learns to notice.

The reader model fits from your answers:
- a form-confusion matrix, by EM, which separates "the air wrecked it" from "it was misread";
- the sharpness of your direction reading, as a likelihood over the physical lean samples,
  integrated over the sector you chose (four seat bearings now; twelve hours in older readings);
- how many forms you can hold at once. The communicator's vocabulary changes from chapter to
  chapter, so the option count varies, and that variation is what makes this limit fittable.

## Files

| File | What it does |
|---|---|
| `core.py` | Copies sections 1–6 of the page (the model, with no DOM in it) into `plume_core.js` so Node can run exactly the page's code. It generates that file, so never edit `plume_core.js` by hand. |
| `lift_traditions.py` | Copies the JSON block (lenses, forms, timeline) and the Sources table of `docs/smoke-reading-traditions.md` into the page between the `TRADITIONS` markers. Run it whenever the writeup changes, before `core.py`. Never edit between the markers by hand. |
| `check_plume.js` | Checks the behaviour the page must keep: a short smooth column, a rare tall thread, Tower still reaching about 55 cm, no false two-streams, a Fork read as two streams, a Stone leaning where the wind turned, a Still that drops the air and lets it back slowly, recovery of the stand-in reader, the presence score on the stand-in and on a toy reader, reliability on a toy pair, the lexicon's clustering and a lens's matching on toy notes, every lifted tradition rule in the page's cue vocabulary with a resolvable source, the protocol search beating the first protocol, and step cost. Exit code = number of failures. Takes about 2 min. |
| `shots.py` | Headless-Edge preview renders of named scenes into `previews/`. |
| `previews/` | `sitting`, `told` (a sitting's end: the card in the dock), `readers` (the timeline and the forms), `card` (cropped by `shots.py card` into `site/previews/language-of-the-plume.webp`), `picked`, `air_day`, `air_night`, `column`, `tower`, `fork`, `ring`, `light`. |

```bash
cd projects/language_of_the_plume
python lift_traditions.py && python core.py && node check_plume.js
python shots.py            # or: python shots.py fork air_day
```

The page also has a self-test in the browser: open it with `#selftest` to append the results.
`#evolve` and `#free` open those modes directly.

## Where this is going

The page began as a channel-coding lab (`smoke-speech-lab.html`): a sender with a grip on the air,
the wind as the adversary, a code evolved for bits per cone. The intent is the inverse, a trainer
for reading a plume as the speech of something that lives in the wind. Phase 1 made the felt sense
the headline and moved every word the table says back to the end of a sitting. Phase 2 made the
sender a mode of the wind itself and retired the grip, read direction from the seat, dropped the
puff codes and the prime opening, and added the air levels and the sittings that come round
again. Phase 3 (the page as it is now) built the lexicon from the notes, made the smoke-reading
traditions lenses on the strength of a sourced writeup in `docs/`, and retargeted the genetic
search from the sender's code to the reader's protocol. Still to do:

- **Phase 4, Field mode.** The same table and dock with no simulation: readings of a real cone,
  then a camera tracker that yields the same instrument readouts.

## Things that bite

- **Headless Edge hangs** if it shares a profile with an open Edge window. `shots.py` gives it a
  throwaway `--user-data-dir`.
- **Don't add per-particle random walk to the smoke.** It zigzags the threads. Spread is carried as
  a width (`w2`), and the stream detector counts only smooth smoke (`tb < 0.8`).
- **Saved readings use the localStorage key `smoke-speech-lab.trials.v1`**, from the page's first
  name. Keep it, or earlier readings stop loading. Trial records from v1 are renumbered on load.
  Records from v3 on carry `flow: true`, with `rt` measured from the start of the sign. Records
  from v4 (phase 2) carry `resp.dir` (0 away, 1 right, 2 toward you, 3 left) with `nDir: 4`, plus
  `hold`, `sit`, `epi` and `replay`. Older records carry `resp.hour` and `k`; the fit reads them
  with twelve-hour sectors and full hold, and they can never come round again. Records from v5
  (phase 3) add `lens`, the lens in use (null for none). Your words live in the artifact's
  collection `lexicon` (the latest row per cue wins, so only `add` is ever used) or under
  `smoke-speech-lab.trials.v1.lexicon` in this browser.
- **On the site menu** (`site/index.html`, with `site/previews/language-of-the-plume.webp`). The page
  itself reaches `site/` only when `deploy.ps1` copies `web/*.html` across.
