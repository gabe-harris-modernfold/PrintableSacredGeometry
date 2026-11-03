/* Invariants the smoke and the reader model must keep. Run after `python core.py`:
     node check_plume.js          (about two minutes)
   Each line prints the measured value, the allowed range, and ok or FAIL; the exit code is the
   number of failures. The ranges are the behaviour the user asked for, not tuning targets. */
"use strict";
const m=require('./plume_core.js');
let fails=0;
function check(name,v,lo,hi){const ok=v>=lo&&v<=hi;if(!ok)fails++;
  console.log(`${ok?'ok  ':'FAIL'} ${name}: ${typeof v==='number'?+v.toFixed(3):v}  [${lo}..${hi}]`);}
const sim=(env,seed,signs)=>{const s=new m.SmokeSim();s.reset({env,gm:1,p:0.35,seed});if(signs)s.schedule(signs);return s;};

// 1. The smooth column is short most of the time; the tall thread is a rare dusk event.
for(const env of ['day','dusk','night']){const s=sim(env,21),Ls=[];let tall=0;
  for(let i=0;i<300*60;i++){s.step(1/60);if(i%30===0){Ls.push(s.L);if(s.L>30)tall++;}}
  Ls.sort((a,b)=>a-b);const med=Ls[Ls.length>>1],frac=tall/Ls.length;
  check(`${env}: median smooth column (cm)`,med,1,9);
  check(`${env}: share of time above 30 cm`,frac,0,env==='dusk'?0.15:0.02);}

// 2. Tower: the communicator can still stand a tall column on purpose.
{const s=sim('dusk',4,[{form:m.F.TOWER,t0:1,T:6}]);let maxL=0;for(let i=0;i<9*60;i++){s.step(1/60);maxL=Math.max(maxL,s.L);}
  check('Tower: tallest column (cm)',maxL,40,70);}

// 3. Streams: the plain wind never reads as two streams; a Fork sign mostly does.
for(const env of ['day','dusk','night']){let two=0,n=0;for(const seed of [1,2,3,4,5]){const s=sim(env,seed);
  for(let i=0;i<40*60;i++){s.step(1/60);if(i>300&&i%6===0){n++;if(s.read.streams===2)two++;}}}
  check(`${env}: false two-streams`,two/n,0,0.02);}
{let a=0,b=0;for(const seed of [5,6,7,8,9]){const s=sim('dusk',seed,[{form:m.F.FORK,t0:4,T:6}]);
  for(let i=0;i<14*60;i++){s.step(1/60);if(i%6===0&&s.t>5.5&&s.t<10){a++;if(s.read.streams===2)b++;}}}
  check('Fork sign read as two streams',b/a,0.6,1);}

// 4. A Stone leans where the wind turned; a Still drops the air, and the wind lets go slowly afterwards.
{const s=sim('dusk',5,[{form:m.F.STONE,angle:Math.PI/2,hold:1,t0:1,T:6}]);for(let i=0;i<5*60;i++)s.step(1/60);
  check('Stone to the right leans right (0 away, 1 right, 2 toward, 3 left)',m.dirOf(s.read.bearing),1,1);}
{const s=sim('day',8,[{form:m.F.STILL,t0:1,T:6,hold:1}]);for(let i=0;i<5*60;i++)s.step(1/60);const during=Math.hypot(s.wind.m[0],s.wind.m[1]);
  for(let i=0;i<14*60;i++)s.step(1/60);const after=Math.hypot(s.wind.m[0],s.wind.m[1]);
  check("Still at midday: the wind's mean 4 s in (cm/s)",during,0,10);check("Still at midday: the wind's mean 12 s after it ends (cm/s)",after,18,60);}

// 5. The reader model recovers a known stand-in, and evolution beats the first guess.
{const st=[];for(let i=0;i<300;i++){const tr=m.makeTrial(100+i);st.push({...tr,resp:m.standInRespond(tr,m.rng32(tr.seed^0x9e37))});}
  const f=m.fitObserver(st);check('stand-in direction kappa (true 6)',f.kappa,3,14);
  // presence: the stand-in feels the signs it can see and never the wind; an unrated "Just wind" pick is a felt wind
  const pr=m.presence(st);check("stand-in presence d'",pr.d,1.5,6);check('stand-in false alarms',pr.fa,0,0.05);
  const toy=[];for(let i=0;i<12;i++)toy.push({truth:i<6?m.F.RING:m.F.WIND,resp:{form:i<6?m.F.RING:m.F.WIND,conf:i<6?5:null}});
  check("toy reader d' from 6 rated hits and 6 unrated wind picks",m.presence(toy).d,2.5,5);
  // reliability: a sign read twice under the same seed, weeks apart
  const rel=[{seed:1,sit:10,epi:0,at:1,truth:m.F.RING,resp:{form:m.F.RING,conf:5}},{seed:1,sit:20,epi:0,at:2,truth:m.F.RING,resp:{form:m.F.RING,conf:1}},
    {seed:1,sit:10,epi:1,at:1,truth:m.F.FORK,resp:{form:m.F.FORK,conf:5}},{seed:1,sit:20,epi:1,at:2,truth:m.F.FORK,resp:{form:m.F.CURL,conf:5}}];
  const R=m.reliability(rel);check('reliability: pairs found',R.n,2,2);check('reliability: same form',R.form,1,1);check('reliability: same felt',R.felt,1,1);
  // lexicon: a note's word must come back in another note of the same cue to count
  const lx=m.lexicon([{resp:{form:m.F.FORK,dir:null,note:'a choice, two roads'}},{resp:{form:m.F.FORK,dir:null,note:'the choice again'}},{resp:{form:m.F.FORK,dir:null,note:'roads roads roads'}},{resp:{form:m.F.STONE,dir:2,note:'welcome'}}]);
  check('lexicon: notes on the fork',lx[0].n,3,3);check('lexicon: words that come back on the fork',lx[0].recur.length,2,2);check('lexicon: top word count',lx[0].recur[0].n,2,2);
  const lens={entries:[{form:'stone',dir:'toward',means:'welcome'},{form:'stone',means:'a lean'},{form:'fork',means:'two roads'}]};
  check('lens: stone toward matches its own rule',m.lensMatch(lens,{form:m.F.STONE,dir:2,env:'dusk',p:0.3}).means==='welcome'?1:0,1,1);
  check('lens: stone left falls to the general rule',m.lensMatch(lens,{form:m.F.STONE,dir:3,env:'dusk',p:0.3}).means==='a lean'?1:0,1,1);
  check('lens: ring has no rule',m.lensMatch(lens,{form:m.F.RING,dir:null,env:'dusk',p:0.3})===null?1:0,1,1);
  const ev=new m.Evolver({M:f.M,kappa:f.kappa,C:f.C,beta:f.beta},7);for(let i=0;i<40;i++)ev.step();
  check('evolved bits per cone (at least 20, and 15% over the first guess)',ev.champ.e.fit,Math.max(20,ev.base.fit*1.15),1e9);}

// 6. The traditions as lifted from docs/smoke-reading-traditions.md: every rule speaks the page's cue
//    vocabulary, every item cites a source in the table, the timeline runs in order in known families.
{const T=m.TRADITIONS,S=T.sources,bad=[];
  for(const l of T.lenses)for(const e of l.entries){
    if(e.form&&!m.FORM_KEY.includes(e.form))bad.push(l.id+': form '+e.form);
    if(e.dir&&(!m.DIR_KEY.includes(e.dir)||e.form!=='stone'))bad.push(l.id+': dir '+e.dir);
    if(e.env&&!['day','dusk','night'].includes(e.env))bad.push(l.id+': env '+e.env);
    if(e.burn&&!['early','mid','late'].includes(e.burn))bad.push(l.id+': burn '+e.burn);
    if(!e.means||!S[e.src])bad.push(l.id+': means/src');}
  for(const x of [...T.lenses,...T.forms,...T.timeline])if(x.src&&!S[x.src])bad.push('unresolved '+x.src);
  for(const f of T.forms){if(!m.FAMILY[f.family])bad.push(f.id+': family '+f.family);if(f.lens&&!T.lenses.some(l=>l.id===f.lens))bad.push(f.id+': lens '+f.lens);}
  for(let i=0;i<T.timeline.length;i++){const t=T.timeline[i];if(!m.FAMILY[t.family])bad.push(t.who+': family');if(i&&t.year<T.timeline[i-1].year)bad.push(t.who+': order');}
  if(bad.length)console.log('  '+bad.slice(0,8).join('\n  '));
  check('traditions: rules, sources, families and order that do not fit',bad.length,0,0);
  check('traditions: lenses lifted',T.lenses.length,1,99);}

// 7. Cost: one stage step stays well inside a frame.
{const s=sim('dusk',3),t0=Date.now();for(let i=0;i<60*60;i++)s.step(1/60);check('sim ms per step',(Date.now()-t0)/3600,0,4);}

console.log(fails?`${fails} failed`:'all ok');process.exitCode=fails;
