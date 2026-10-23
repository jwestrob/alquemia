#!/usr/bin/env python3
"""Narrow immutable-source and actual-preparation joins; no native/control-suite replay."""
import ast,hashlib,json,re,time
from pathlib import Path
ROOT=Path.cwd();REL=Path('reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_controls_v1');DB=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL;WB=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL;OUT=WB/'independent_review_v1';OUT.mkdir(exist_ok=False);checks=[];t=time.monotonic();cpu=time.process_time()
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def ck(q,n):
 checks.append({'check':n,'passed':bool(q)});assert q,n
def verify(p):ck(pin(p['path'])==p,'pin:'+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,q):Path(p).write_text(json.dumps(q,indent=2)+'\n')
ready=json.loads((DB/'READY.json').read_text())
for v in ready.values():
 if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
p=read(ready['panel']);prep=read(ready['preparation']);cfg=read(p['configuration']);parent=read(p['parent_panel']);scout=read(p['scout_collection']);rev=read(p['scout_actual_review']);controls=read(ready['saved_controls']);failed=read(p['failed_tight_controls'])
ck(ready['implementation']==p['adapter']==prep['implementation']==cfg['implementation'],'source joins')
ck(prep['panel']==ready['panel'] and prep['status']=='passed_tighter_controls_preparation' and prep['new_native_calls']==0,'actual zero-call preparation')
ck(controls['status']=='passed_saved_controls' and controls['implementation']==p['adapter'] and controls['panel']==ready['panel'] and controls['all_passed'] and len(controls['checks'])==22 and all(controls['checks'].values()) and controls['new_native_calls']==0,'existing22 controls inspected only')
ck(rev['status']=='passed_tighter_dy_acetate_recovery_actual_review' and rev['collection']==p['scout_collection'] and rev['implementation']==p['recovery_adapter'] and rev['diagnostic_limit_passed'],'actual tighter3 review join')
ck(failed['status']=='local_controls_failed' and failed['complete']==18 and not failed['all_controls_passed'],'failed original preserved')
ck(p['assets']==parent['assets'] and p['orca']==parent['orca'] and p['controls']==parent['controls'],'original basis/binary/canonical evidence unchanged')
for k in ['tighter_adapter','recovery_adapter','parser','native_input_adapter']:verify(p[k])
for aset in p['assets'].values():
 for x in aset.values():
  if isinstance(x,dict) and set(x)=={'path','sha256'}:verify(x)
ck(p['tier']=={'TCutPairs':1e-6,'TCutPNO':1e-8,'TCutMKN':.0001,'TCutDO':.005},'exact tighter recipe')
old={c['cell_id']:c for c in parent['cells'] if c['stage']=='controls'};new={c['cell_id']:c for c in p['cells']};reuse=set(p['reused_cell_ids']);newids=set(p['new_cell_ids'])
ck(len(new)==18 and set(new)==set(old) and len(reuse)==3 and len(newids)==15 and newids.isdisjoint(reuse) and newids|reuse==set(new),'full18 finite split')
ck(reuse=={'control__cc__acetate_bi__Dy__'+r for r in ['complex','ion_ghost_ligand','ligand_ghost_ion']},'exact reused roles')
ck(p['maximum_new_native_calls']==15 and p['declared_controls']==18 and p['qualification_check_count']==35 and p['ranks']==8 and p['nodes']==1,'declared limits')
for cid,c in new.items():
 before=old[cid];e={k:v for k,v in before.items() if k!='template_input'};e['tier']='tighter';ck({k:v for k,v in c.items() if k!='template_input'}==e,cid+'/literal original metadata')
 original=verify(before['template_input']).read_text();actual=verify(c['template_input']).read_text();expected=original.replace('%pal nprocs 4 end','%pal nprocs 8 end').replace(' TCutPairs 1e-05\n',' TCutPairs 1e-06\n').replace(' TCutPNO 1e-07\n',' TCutPNO 1e-08\n').replace(' TCutMKN 0.001\n',' TCutMKN 0.0001\n');ck(actual==expected,cid+'/only declared deck delta')
 ck('! RHF DLPNO-CCSD(T1) TightPNO def2-TZVP def2-TZVP/C NoFrozenCore VeryTightSCF NoAutostart PrintBasis\n' in actual and ' STol 1e-8\n' in actual and ' TriTolE 1e-8\n' in actual and ' TCutDO 0.005\n' in actual,cid+'/method and stops')
 pos=actual.split('* xyz ',1)[1].splitlines();ck(pos[0]==str(c['charge'])+' 1',cid+'/charge effective singlet');coords=[x.split() for x in pos[1:] if x.strip()!='*'];ck(len(coords)==len(c['symbols']),cid+'/all centers')
 ck([x[0] for x in coords]==[s+(':' if g else '') for s,g in zip(c['symbols'],c['ghost_mask'])],cid+'/ghost mask');ck([[float(v) for v in x[1:]] for x in coords]==c['positions_A'],cid+'/unchanged coordinates')
 ck(('DelECP '+c['metal']) in actual if c['ghost_mask'][0] else ('DelECP '+c['metal']) not in actual,cid+'/metal ghost ECP');ck('NewAuxCGTO '+c['metal'] in actual and 'NewGTO '+c['metal'] in actual,cid+'/metal orbital and AuxC retained')
 ck(sum(c['expected_effective_nuclear_charges'])-c['charge']==c['expected_electrons']==c['expected_correlated_electrons'],cid+'/electron and all-electron correlation accounting')
 ck(sum(c['expected_basis_functions_by_center'])==c['expected_basis_functions'],cid+'/basis inventory')
 ck(c['physical_multiplicity_of_parent']==(1 if c['metal']=='La' else 6),cid+'/physical versus effective state')
 if cid in reuse:
  q=p['reused'][cid];stored=scout['rows'][q['original_scout_cell_id']];ck({k:v for k,v in q.items() if k not in ['reused_from_tighter_scout','original_scout_cell_id']}==stored,cid+'/exact admitted reused object')
  for key in ['result','input','output','start']:verify(q[key])
  ck(verify(q['input']).read_text()==actual,cid+'/byte-identical reused deck')
 else:ck(prep['input_deltas'][cid]['new_native_required'],cid+'/new required')
 # Canonical values are frozen from the preparation's18 raw replays; verify producer identities without rerunning parsers.
 canon=p['controls'][c['canonical_id']]
 for key in ['result','output']:
  if key in canon:verify(canon[key])
source=verify(p['adapter']).read_text();tree=ast.parse(source);fn={n.name:ast.get_source_segment(source,n) for n in tree.body if isinstance(n,ast.FunctionDef)}
ck('deck=tight.deck;parse=tight.parse' in source,'frozen deck/parser aliases')
ck('old.control_score(p,rows)' in fn['score'] and "len(q['checks'])==35" in fn['score'] and "len(rows)==18" in fn['score'],'direct original full35 scoring')
ck("cell_id in p['new_cell_ids']" in fn['run'] and 'not out.exists()' in fn['run'],'no reused starts/retry')
ck('nodes==1' in fn['resources'] and 'tight.resources(env)' in fn['resources'],'one-node actual eight-rank resource contract')
ck('recovery.admit_replacement' in fn['reused_rows'] and 'tight.admit_result' in fn['reused_rows'],'producer-specific reused admissions')
ck('len(good)==18' in fn['collect'] and 'cluster_release_authorized=False' in fn['collect'] and 'physical_accuracy_qualified=False' in fn['collect'],'complete denominator no cluster/accuracy release')
wrapper=verify(ready['run_wrapper']).read_text();ck('#SBATCH --nodes=1' in wrapper and '#SBATCH --ntasks=8' in wrapper and '#SBATCH --cpus-per-task=1' in wrapper and not re.search(r'^#SBATCH\s+--(?:mem|time)',wrapper,re.M),'one node eight threads/default memory-time')
inv=read(ready['result_inventory']);ck(set(inv)==newids,'only15 new result inventory')
for cid,path in inv.items():ck(Path(path)==Path(p['execution_root'])/'tasks'/cid/'execution/RESULT.json',cid+'/one immutable path')
orig=verify(p['native_input_adapter']).read_text();score=ast.get_source_segment(orig,next(n for n in ast.parse(orig).body if isinstance(n,ast.FunctionDef) and n.name=='control_score'))
ck('<=1e-6' in score and score.count('<=.25')==3 and "r['aux_C']" in score and "paired[d][i]-paired['water'][i]" in score,'unchanged35 categories/signs/tolerances inspected')
save(OUT/'CHECKS.json',checks);receipt={'status':'passed_tighter_controls_source_review','implementation':p['adapter'],'panel':ready['panel'],'preparation':ready['preparation'],'ready':pin(DB/'READY.json'),'declaration':p['declaration'],'configuration':p['configuration'],'scout_actual_review':p['scout_actual_review'],'scout_collection':p['scout_collection'],'original_failed_controls':p['failed_tight_controls'],'saved_controls':ready['saved_controls'],'review_implementation':pin(__file__),'checks':pin(OUT/'CHECKS.json'),'check_count':len(checks),'new_native_calls':0,'new_test_suite_runs':0,'existing_control_checks_read':22,'new_declared':15,'reused':3,'full_declared':18,'original_qualification_checks':35,'blockers':[],'limitations':['The35 numerical/canonical-overlap gates are unchanged, not physical accuracy guarantees.','Effective restricted frozen-f/ECP states remain distinct from physical Dy state.','Collector native cost totals include admitted good rows only; root must add scheduler costs for all failed/unavailable jobs.','No cluster release, new thresholds, core/basis change, or automatic retry.'],'review_cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu}}
save(OUT/'REVIEW.json',receipt);save(DB/'independent_review_v1/REVIEW.json',receipt);print(json.dumps({'review':pin(OUT/'REVIEW.json'),'checks':len(checks)},indent=2))
