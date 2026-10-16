// Runs the page's DOM-free core under node and writes the demo forms as STL.
//
//   node run_pure.mjs                 -> js_cylinder_<layout>[_hollow].stl, demo_field.f32, demo_field.json,
//                                        js_cavity_<open|closed|halves|moment>.stl, demo_cavity.json, demo_cavity_field.f32,
//                                        js_wave_*.stl, demo_wave.json, demo_voice.f32
//   node run_pure.mjs some.wav        -> also js_<basename>.stl from a real WAV via VSF.parseWav
//
// The block between /* @pure-begin */ and /* @pure-end */ in
// ../../web/volumetric-sound-forms.html is sliced out verbatim and evaluated,
// so what runs here is byte-for-byte what runs in the browser and its worker.
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const htmlPath = path.join(here, '..', '..', 'web', 'volumetric-sound-forms.html');
const html = fs.readFileSync(htmlPath, 'utf8');
const a = html.indexOf('/* @pure-begin */'), b = html.indexOf('/* @pure-end */');
if (a < 0 || b <= a){ console.error('pure markers not found in', htmlPath); process.exit(2); }
vm.runInThisContext(html.slice(a, b), { filename: 'vsf-pure.js' });
const VSF = globalThis.VSF;
if (!VSF){ console.error('VSF did not attach to globalThis'); process.exit(2); }

let failures = 0;
const check = (c, msg) => { if (!c){ failures++; console.error('FAIL:', msg); } };

const st = VSF.selfTestPure();
console.log(`selfTestPure: ${st.asserts} asserts, ${st.failures} failures`);
failures += st.failures;
const t0c = Date.now(), stc = VSF.selfTestCavity();
console.log(`selfTestCavity: ${stc.asserts} asserts, ${stc.failures} failures (${Date.now() - t0c} ms)`);
failures += stc.failures;

const sr = 44100, samples = VSF.synthDemo(sr, 6);
const base = VSF.defaults(); base.nt = 96; base.nf = 48;
const runs = [];
const summary = { sr, seconds: 6, P: base, runs: [] };
const chartRuns = [];
for (const layout of ['mirror', 'plain', 'transpose']) for (const hollow of [0, 1]) chartRuns.push({ layout, hollow, slope: 0 });
chartRuns.push({ layout: 'mirror', hollow: 1, slope: 45 });          // the support-free pass
for (const { layout, hollow, slope } of chartRuns){
  const P = { ...base, layout, hollow, slope };
  const t = Date.now();
  const r = VSF.pipeline(samples, sr, P);
  const s = r.stats, name = `js_cylinder_${layout}${hollow ? '_hollow' : ''}${slope ? '_slope' + slope : ''}.stl`;
  const buf = VSF.writeSTL(r.mesh.pos, r.mesh.idx, 'Promethean Theoria - volumetric sound form - demo - ' + layout + (hollow ? ' hollow' : ''));
  fs.writeFileSync(path.join(here, name), Buffer.from(buf));
  check(s.badEdges === 0, `${name}: ${s.badEdges} edges not shared by exactly two faces`);
  check(s.badDir === 0, `${name}: ${s.badDir} reversed edges`);
  check(s.euler === 2, `${name}: Euler ${s.euler}`);
  check(s.vol > 0, `${name}: volume ${s.vol}`);
  check(s.fits, `${name}: exceeds the bed`);
  check(buf.byteLength === 84 + 50 * s.F, `${name}: STL byte length`);
  if (slope) check(s.over <= slope + 1e-3, `${name}: overhang ${s.over.toFixed(2)} deg past the ${slope} deg limit`);
  console.log(`${name.padEnd(34)} V ${String(s.V).padStart(7)}  F ${String(s.F).padStart(7)}  Euler ${s.euler}  odd ${s.badEdges}  rev ${s.badDir}  vol ${(s.vol / 1000).toFixed(1).padStart(7)} cm3  over ${s.over.toFixed(0).padStart(2)} deg  size ${s.size.map(v => v.toFixed(1)).join(' x ')}  ${Date.now() - t} ms`);
  runs.push({ layout, hollow, r });
  summary.runs.push({ name, layout, hollow, slope, over: s.over, nT: r.nT, nF: r.nF, nR: r.nR, nC: r.nC, nRi: r.mesh.nRi, nFoot: s.nFoot, dz: s.dz,
    V: s.V, E: s.E, F: s.F, euler: s.euler, vol: s.vol, badEdges: s.badEdges, badDir: s.badDir, size: s.size, scale: s.scale });
}
// the mirror/solid run's field and band map, for the numpy comparison
const ref = runs[0].r;
fs.writeFileSync(path.join(here, 'demo_field.f32'), Buffer.from(ref.field.buffer, ref.field.byteOffset, ref.field.byteLength));
fs.writeFileSync(path.join(here, 'demo_surface.f32'), Buffer.from(ref.S.buffer, ref.S.byteOffset, ref.S.byteLength));
summary.field = { nT: ref.nT, nF: ref.nF, N: ref.N, hop: ref.hop, K: ref.K, binToBand: Array.from(ref.map.binToBand), count: Array.from(ref.map.count) };
summary.surface = { nR: ref.nR, nC: ref.nC };
fs.writeFileSync(path.join(here, 'demo_field.json'), JSON.stringify(summary, null, 1));

// ── the cavity form: the 6 s demo in the default frustum, four ways ──
{
  const P0 = { ...VSF.defaults(), form: 'cavity', cvox: 1.5 };
  const an = VSF.analyse(samples, sr, P0);
  const cruns = [
    { name: 'js_cavity_open.stl', P: { ...P0 } },
    { name: 'js_cavity_closed.stl', P: { ...P0, ctop: 'closed' } },
    { name: 'js_cavity_halves.stl', P: { ...P0, chalf: 1 } },
    { name: 'js_cavity_moment.stl', P: { ...P0, cmode: 'moment', ct: 3.2 } },
  ];
  const cav = { sr, seconds: 6, P: P0, runs: [] };
  for (const run of cruns){
    const t = Date.now(), r = VSF.cavityPipeline(an, sr, run.P), s = r.stats;
    const buf = VSF.writeSTL(r.mesh.pos, r.mesh.idx, 'Promethean Theoria - volumetric sound form - demo - cavity');
    fs.writeFileSync(path.join(here, run.name), Buffer.from(buf));
    const minShells = s.half ? 2 : 1, maxShells = minShells + s.sealed;
    check(s.badEdges === 0 && s.badDir === 0, `${run.name}: ${s.badEdges} odd, ${s.badDir} reversed edges`);
    check(s.vol > 0 && s.fits && s.shells >= minShells && s.shells <= maxShells, `${run.name}: vol ${s.vol} fits ${s.fits} shells ${s.shells} of ${minShells}..${maxShells}`);
    check((s.euler & 1) === 0 && Number.isInteger(s.genus) && s.genus >= 0, `${run.name}: Euler ${s.euler}`);
    // coarse sanity only: a voxel count is a biased estimator on a wall 1.5 voxels thick; the precise checks are the
    // bare shell against its closed form (below) and the numpy/VTK isosurface in sound_form_reference.py --compare
    check(Math.abs(s.voxVol / s.vol - 1) < 0.15, `${run.name}: mesh volume ${s.vol.toFixed(0)} vs voxel count ${s.voxVol.toFixed(0)}`);
    check(buf.byteLength === 84 + 50 * s.F, `${run.name}: STL byte length`);
    console.log(`${run.name.padEnd(34)} V ${String(s.V).padStart(7)}  F ${String(s.F).padStart(7)}  Euler ${s.euler}  genus ${s.genus}  shells ${s.shells} (sealed ${s.sealed})  odd ${s.badEdges}  rev ${s.badDir}  vol ${(s.vol / 1000).toFixed(1).padStart(6)} cm3 (voxels ${(s.voxVol / 1000).toFixed(1)})  islands ${s.islands}  size ${s.size.map(v => v.toFixed(1)).join(' x ')}  ${Date.now() - t} ms`);
    const gq = r.field.grid;
    cav.runs.push({ name: run.name, P: run.P, stats: { ...s, min: s.min, max: s.max, size: s.size }, drive: { row: r.drive.row, gain: r.drive.gain, bmax: r.drive.bmax, s: r.drive.s, window: r.drive.window },
      grid: { nx: gq.nx, ny: gq.ny, nz: gq.nz, c: gq.c, below: gq.below, vx: gq.vx, vz: gq.vz, x0: gq.x0, y0: gq.y0, z0: gq.z0, tw: gq.tw, flr: gq.flr, zRim: gq.zRim, zTop: gq.zTop } });
    if (run.name === 'js_cavity_open.stl'){
      const m = r.modes, d = r.drive, g = r.field.grid;
      cav.modes = { a0: m.a0, Hn: m.Hn, nz: m.nz, h: m.h, mMax: m.mMax, count: m.count,
        groups: m.groups.map(q => ({ m: q.m, n: q.n, kappa: q.kappa, C: q.C, Jw: q.Jw, modes: q.modes })),
        m: Array.from(m.m), g: Array.from(m.g), k: Array.from(m.k), N: Array.from(m.N) };
      cav.drive = { s: d.s, Q: d.Q, bmax: d.bmax, df: d.df, gain: d.gain, row: d.row, window: d.window,
        f: Array.from(d.f), D: Array.from(d.D), phis2: Array.from(d.phis2), w: Array.from(d.w) };
      cav.grid = { nx: g.nx, ny: g.ny, nz: g.nz, c: g.c, below: g.below, vx: g.vx, vz: g.vz, x0: g.x0, y0: g.y0, z0: g.z0, tw: g.tw, flr: g.flr, zRim: g.zRim, zTop: g.zTop };
      cav.field = { max: r.field.max, gain: r.field.gain, smooth: r.field.smooth };
      const p2 = r.field.p2;
      fs.writeFileSync(path.join(here, 'demo_cavity_field.f32'), Buffer.from(p2.buffer, p2.byteOffset, p2.byteLength));
    }
  }
  // the bare vessel (no lobe reaches 99 %) against the closed-form volume of a conical tube on a floor disc
  {
    const P = { ...P0, cfill: 0, cvox: 1.2 }, r = VSF.cavityPipeline(an, sr, P), g = r.field.grid;
    let tube = 0; const n = 4000;
    for (let i = 0; i < n; i++){ const a = P.cr0 + (P.cr1 - P.cr0) * (i + 0.5) / n; tube += Math.PI * ((a + g.tw) ** 2 - a * a) * P.ch / n; }
    const exact = tube + Math.PI * (P.cr0 + g.tw) ** 2 * P.cflr, rel = r.stats.vol / exact - 1;
    check(Math.abs(rel) < 0.02, `bare vessel volume ${(r.stats.vol / 1000).toFixed(2)} vs closed form ${(exact / 1000).toFixed(2)} cm3`);
    console.log(`bare vessel (no lobes)                mesh ${(r.stats.vol / 1000).toFixed(2)} cm3 vs closed form ${(exact / 1000).toFixed(2)} cm3 (${(rel * 100).toFixed(2)} %)`);
    cav.bare = { vol: r.stats.vol, exact };
  }
  fs.writeFileSync(path.join(here, 'demo_cavity.json'), JSON.stringify(cav));
}

// ── the waveform solid: the spoken-vowel demo four ways (three with thickness by pitch); its samples, windows,
//    pitch track and profiles for numpy ──
{
  const x = VSF.synthVoice(sr, VSF.VOICE_SECONDS);
  fs.writeFileSync(path.join(here, 'demo_voice.f32'), Buffer.from(x.buffer, x.byteOffset, x.byteLength));
  const W0 = VSF.defaults();
  const wruns = [
    { name: 'js_wave_fins.stl', P: { ...W0 } },
    { name: 'js_wave_fins_halves.stl', P: { ...W0, whalf: 1 } },
    { name: 'js_wave_flat_rms_fixed_halves.stl', P: { ...W0, whalf: 1, wtip: 0, wn: 48, wmeas: 'rms', wgam: 0.6, wpitch: 0 } },
    { name: 'js_wave_waveform.stl', P: { ...W0, wmode: 'wave', t0: 0.5, t1: 0.8 } },
  ];
  const wav = { sr, seconds: VSF.VOICE_SECONDS, runs: [] };
  for (const run of wruns){
    const t = Date.now(), r = VSF.waveShape(x, sr, run.P), s = r.stats;
    const buf = VSF.writeSTL(r.mesh.pos, r.mesh.idx, 'Promethean Theoria - volumetric sound form - demo - waveform');
    fs.writeFileSync(path.join(here, run.name), Buffer.from(buf));
    const bodies = s.half ? 2 : 1;
    check(s.badEdges === 0 && s.badDir === 0, `${run.name}: ${s.badEdges} odd, ${s.badDir} reversed edges`);
    check(s.shells === bodies && s.euler === 2 * bodies, `${run.name}: ${s.shells} shells, Euler ${s.euler}`);
    check(s.vol > 0 && s.fits && Math.abs(s.vol / s.exact - 1) < 1e-5, `${run.name}: volume ${s.vol} vs exact lathe ${s.exact}, fits ${s.fits}`);
    if (s.half) check(s.over < 0.5 && Math.abs(s.min[2]) < 1e-6, `${run.name}: halves overhang ${s.over.toFixed(3)} deg, z from ${s.min[2]}`);
    check(buf.byteLength === 84 + 50 * s.F, `${run.name}: STL byte length`);
    if (run.P.wmode !== 'wave' && run.P.wpitch) check(s.pitched && s.voiced > 0.8 * s.n, `${run.name}: pitched ${s.pitched}, ${s.voiced} of ${s.n} windows voiced`);
    const pitchNote = s.pitched ? `  pitch ${s.fLo.toFixed(0)}-${s.fHi.toFixed(0)} Hz (${s.voiced} voiced)` : '';
    console.log(`${run.name.padEnd(34)} V ${String(s.V).padStart(7)}  F ${String(s.F).padStart(7)}  Euler ${s.euler}  shells ${s.shells}  odd ${s.badEdges}  rev ${s.badDir}  vol ${(s.vol / 1000).toFixed(1).padStart(6)} cm3 (lathe ${((s.vol / s.exact - 1)).toExponential(1)})  over ${s.over.toFixed(0).padStart(2)} deg  ${s.mode === 'wave' ? s.n + ' stations' : s.fins + ' fins'}${pitchNote}  size ${s.size.map(v => v.toFixed(1)).join(' x ')}  ${Date.now() - t} ms`);
    wav.runs.push({ name: run.name, P: run.P, s0: r.amp.s0, s1: r.amp.s1, raw: Array.from(r.amp.raw), amp: Array.from(r.amp.a),
      pitch: r.pitch ? { f: Array.from(r.pitch.f), voiced: Array.from(r.pitch.voiced), dip: Array.from(r.pitch.dip), fMed: r.pitch.fMed } : null,
      thick: r.thick ? Array.from(r.thick.t) : null, X: Array.from(r.prof.X), R: Array.from(r.prof.R), stats: s });
  }
  fs.writeFileSync(path.join(here, 'demo_wave.json'), JSON.stringify(wav));
}

if (process.argv[2]){
  const wavPath = process.argv[2];
  const ab = fs.readFileSync(wavPath);
  const w = VSF.parseWav(ab.buffer.slice(ab.byteOffset, ab.byteOffset + ab.byteLength));
  const P = { ...VSF.defaults(), nt: 192, nf: 96 };
  const r = VSF.pipeline(w.samples, w.sr, P);
  const s = r.stats, name = 'js_' + path.basename(wavPath).replace(/\.[^.]+$/, '') + '.stl';
  fs.writeFileSync(path.join(here, name), Buffer.from(VSF.writeSTL(r.mesh.pos, r.mesh.idx)));
  check(s.badEdges === 0 && s.badDir === 0 && s.euler === 2 && s.vol > 0, `${name}: not closed`);
  console.log(`${name.padEnd(34)} ${w.sr} Hz ${w.channels} ch ${(w.samples.length / w.sr).toFixed(1)} s -> F ${s.F}  Euler ${s.euler}  vol ${(s.vol / 1000).toFixed(1)} cm3  size ${s.size.map(v => v.toFixed(1)).join(' x ')}`);
}
console.log(failures ? `${failures} FAILURE(S)` : 'all forms closed, oriented, inside the bed (halves: two bodies; sealed cavity: one body, its cavity a shell of its own)');
process.exit(failures ? 1 : 0);
