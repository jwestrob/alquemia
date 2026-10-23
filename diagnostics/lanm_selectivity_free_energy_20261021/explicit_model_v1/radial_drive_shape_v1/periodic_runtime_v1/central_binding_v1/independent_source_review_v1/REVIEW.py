"""Saved source/input binding review only; no preparation or scientific replay."""
import datetime,hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent;PINS={};CHECKS=[]
def pin(p):
 p=Path(p).resolve()
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 return {'path':str(p),'sha256':h}
def verify(p):
 if p['path'] not in PINS:assert pin(p['path'])==p,p;PINS[p['path']]=p
 else:assert PINS[p['path']]==p
 return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def check(n,v):
 assert v,n
 CHECKS.append({'label':n,'passed':True})
def split(text):
 base=[];response=[]
 for line in text.splitlines(keepends=True):
  words=line.split('#',1)[0].split()
  if words and words[0].lower() in ('polar-local-response','polar-radial-drive'):response.append(words)
  else:base.append(line)
 return ''.join(base),response
def main():
 wall=time.monotonic();cpu=time.process_time();rp=pin(P/'READY.json');r=read(rp)
 check('Exact requested frozen READY',rp['sha256']=='dde404f4c9e09e72f16763b9a2107b8291db56b8663f1cab16fb7d7dfb967185')
 for p in r.values():
  if isinstance(p,dict) and set(p)=={'path','sha256'}:verify(p)
 t=read(r['required_source_review_template']);check('Exact receipt template',t==r['required_source_review'])
 s=read(r['inputs']);d=read(r['root_declaration']);cfg=read(r['configuration']);prep=read(r['preparation']);m=read(r['manifest']);c=read(r['composition']);meta=read(r['physics_metadata'])
 for p in s.values():verify(p)
 check('Root exact parent pins',all(s[k]==p for k,p in d['parents'].items()))
 check('Source/configuration/result binding',cfg['implementation']==prep['implementation']==r['implementation'] and cfg['inputs']==r['inputs'] and prep['configuration']==r['configuration'])
 check('35 saved preparation checks pass without execution',prep['status']=='prepared_periodic_central_binding_only' and prep['error'] is None and len(prep['checks'])==r['checks_passed']==35 and all(v==0 for v in prep['counts'].values()) and prep['unexpected_input_mismatches']==[])
 check('Candidate remains declared20605 terminal5d992',r['candidate_key']==d['candidate_key']==c['candidate_key']==prep['candidate_key']=='5d992dbffb9fbe646102782cb90b53182deff7e384838ffd809683e1a187fa78' and r['candidate']==s['candidate']==c['candidate']==m['model']['candidate'])
 admitted=read(s['fit_root_admission']);candidate=read(s['candidate'])
 check('Literal vector and existing native admission retained',admitted['selected_artifact']==s['candidate'] and c['parameters']==candidate['parameters'] and candidate['status']=='complete')
 check('Model/result native coefficients and both radial slots exact',m['model']['coefficients']==prep['coefficients']==c['coefficients']==meta['coefficients'] and m['model']['radial']==prep['radial']==c['radial']==meta['radial'] and set(c['coefficients'])==set(c['radial'])=={'900','901'})
 check('Four exact roles',set(m['rows'])==set(c['roles'])=={'central','BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO'} and m['rows']==c['roles'])
 old=read(s['baseline_manifest'])['rows']['central']['inputs'];mapping=read(s['BIO_mapping']);mapped=next(x for x in mapping['aqueous_endpoints'] if x['id']=='La_aqueous_BIO')
 check('Original mapped geometry and BIO library authority',old['xyz']['sha256']==mapped['files']['xyz']['sha256'] and old['parameters']['sha256']==mapped['files']['parameters']['sha256']==mapping['BIO_parameters_preserved_byte_for_byte']['sha256'])
 texts={}
 for sid,row in m['rows'].items():
  for p in row['inputs'].values():verify(p)
  check('Unchanged xyz/prm and truthful provenance '+sid,row['inputs']['xyz']['sha256']==old['xyz']['sha256'] and row['inputs']['parameters']['sha256']==old['parameters']['sha256'] and row['source_inputs']==row['inputs'] and row['baseline_source_inputs']==old and row['cell']==[49.324]*3)
  texts[sid]=verify(row['inputs']['key']).read_text()
 base,active=split(texts['central'])
 check('Identical nonresponse bytes across all4 roles',all(split(text)[0]==base for text in texts.values()))
 check('Both absent rows contain zero response records',not split(texts['BYPASS_OLD_ABSENT'])[1] and not split(texts['BYPASS_CLONE_ABSENT'])[1] and texts['BYPASS_OLD_ABSENT']==texts['BYPASS_CLONE_ABSENT'])
 records={(v[0].lower(),v[1]):v for v in active}
 check('Both old and radial ACTIVE records explicit',len(active)==len(records)==4 and set(records)=={(kind,metal) for kind in ('polar-local-response','polar-radial-drive') for metal in ('900','901')})
 check('ACTIVE native literal records match mapped metadata',all(records[('polar-local-response',metal)][2]=='402' and list(map(float,records[('polar-local-response',metal)][3:]))==c['coefficients'][metal] and records[('polar-radial-drive',metal)][2]=='402' and float(records[('polar-radial-drive',metal)][3])==c['radial'][metal] for metal in ('900','901')))
 zero=split(texts['BYPASS_CLONE_ZERO'])[1]
 check('ZERO bypass has all4 explicit zero records',len(zero)==4 and {(v[0].lower(),v[1]) for v in zero}==set(records) and all(v[2]=='402' and all(float(x)==0. for x in v[3:]) for v in zero))
 prior=verify(old['key']).read_text()
 check('Original POLPAIR text unchanged', [l for l in prior.splitlines() if l.lower().startswith('polpair ')]==[l for l in base.splitlines() if l.lower().startswith('polpair ')])
 pairs=[v for line in base.splitlines() if (v:=line.split()) and v[0].lower()=='vdwpr']
 check('Both candidate pair values common to all4 roles',len(pairs)==2 and all(v[1:3]==['73',metal] and list(map(float,v[3:]))==c['effective_a_b_c_d_R_epsilon_POLPAIR'][name][4:6] for v,(metal,name) in zip(pairs,[('900','La'),('901','Dy')])))
 check('Fixed cutoff12 settings and type counts',m['model']['settings']==prep['settings']==d['settings'] and prep['type_counts']=={'900':1,'402':3980,'403':7960,'405':7,'414':10})
 oldmeta=read(s['old_physics_metadata']);keep=('type_parameters','mutual_scales','pair_overrides','metal_site_index','sites','water_count','all_alpha_positive')
 check('Only intrinsic metadata inherited; no old H/PME',all(meta[k]==oldmeta[k] for k in keep) and meta['initialized_PME'] is None and meta['H_evaluated'] is False)
 check('All new state/build/runtime gates unbound',m['status']=='prepared_periodic_v2_central_inputs_not_executable' and all(m[k] is None for k in ('build','build_admission','runtime_release','initializer','certificate','central_parity_admission','central_scout_admission')) and r['probe_source'] is r['probe_data'] is r['initialized_state'] is None)
 check('Engineering only and explicit kinetic mass convention',c['status']=='prepared_only_not_executable' and c['common_metal_mass_Da']==138.91 and c['mass_convention_requires_dynamics_decision'] is True and c['physical_accuracy_qualified'] is False and prep['scientific_admission'] is False and prep['fresh_V2_chain_required'] is True)
 report='Independent radial periodic central binding source review: PASS\n\nThe binding uses exactly admitted20605 terminal5d992, never the new1/16 witness. The source maps normalized u/v to native c/d using the pinned original formulas, appends literal e900/e901, and preserves both unused Dy slots. Both candidate VDWPR values are identical across ACTIVE and all3 bypasses. The absent keys omit all response records; ZERO retains all4 explicit zero records. POLPAIR and other nonresponse settings remain common.\n\nFour xyz/prm hashes match the existing mapped11958-site La box and BIO library:3980 waters,7Na,10Cl,49.324A orthogonal cell. No geometry/topology or intrinsic alpha/damping/salt record is changed. Old physics metadata is limited to intrinsic properties; no old H, initialized PME, sparse matrix or certificate is admitted. Cutoff12 V2 build/init/certificate/parity/scout and probes remain unbound. The common138.91Da metal mass is retained explicitly as a later dynamics choice.\n\nRead35 saved preparation checks and literal key/provenance records without rerunning preparation, h/H, native, FFT, assembly, factor or solver. This admits source/input correspondence only, not runtime stability, physical accuracy, production or sampling.\n'
 (HERE/'REPORT.txt').write_text(report)
 result=dict(t);result.update(ready=rp,physics_metadata=r['physics_metadata'],report=pin(HERE/'REPORT.txt'),reviewer_implementation=pin(__file__),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,verified_pins=list(PINS.values()),findings=[],source_only=True,preparation_repeated=False,native_admission=False,physical_accuracy_qualified=False,counts={'saved_preparation_checks_read':35,'native_calls':0,'builds':0,'FFT_replays':0,'assemblies':0,'factorizations':0,'solver_calls':0,'new_preparations':0},cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1})
 with (HERE/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(CHECKS),'pins':len(PINS),'cost':result['cost']}))
if __name__=='__main__':main()
