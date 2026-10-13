#!/usr/bin/env python3
"""Two-source MxaF solvent recovery using frozen preparation, components and bands."""
import argparse,copy,fcntl,hashlib,importlib.util,json,os,subprocess,time
from collections import Counter
from pathlib import Path

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
STAGE=ROOT/'diagnostics/orca_parameter_audit_20261012/slurm_resumption/GFN2_STAGE.py'
STAGE_SHA='b7ac3b133876f4cc9fa521beaee9a539e6213cb0f60a7206145b2e28c119dc1f'
SCOPE_SHA='6e31585165092f10a8294d31c7def13331fe0e6196949b5a3eb94dc75cf59ccb'
RECOVERY_SHA='14a14cbd8acb561fd255e582043928bed378a2676d64f831e7b7b8deedd26037'
IDS=[f'p12293-pqq-la_model__conditioned_Ca__seed-1_sample-{i}__envelope_e88ef960dad9c35e' for i in (1,4)]
PROTOCOL='MxaF_same_policy_q0_restart_corrected_GFN2_recovery_v1'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stage():
 assert sha(STAGE)==STAGE_SHA,'frozen Slurm executor changed'
 spec=importlib.util.spec_from_file_location('mxaf_frozen_solvent_stage',STAGE);s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
 placement=s.allocation()
 actual=subprocess.run(['scontrol','show','job','-o',placement['job_id']],text=True,capture_output=True,check=True).stdout
 fields=dict(t.split('=',1) for t in actual.split() if '=' in t)
 assert fields['JobState']=='RUNNING' and int(fields['NumCPUs'])==12 and int(fields['NumNodes'])==1
 os.environ.update(s.THREAD_ENV)
 return s,s.driver(),placement

def prepare(g,scope_path,collection_path,primary_path,base):
 f,r=g.f,g.r;scope=f.read(scope_path);native=f.read(collection_path);recovery=f.read(native['manifest'])
 assert f.pin(scope_path)['sha256']==SCOPE_SHA and native['terminal'] is True
 assert native['manifest']['sha256']==RECOVERY_SHA and recovery['protocol_id']=='MxaF_Ca_conditioned_same_policy_La_q0_restart_v1'
 assert len(scope['rows'])==scope['counts']['declared_Ca_sources']==125
 primary=f.read(primary_path);prior=f.read(primary['manifest']);old_native=f.read(prior['MACE_collection'])
 assert len(primary['rows'])==125 and primary['counts']['operational']['unavailable']==18
 assert prior['scope']==f.pin(scope_path) and prior['model']==scope['model']==recovery['model']
 assert recovery['primary_collection']==prior['MACE_collection'] and recovery['primary_manifest']==old_native['manifest']
 assert recovery['pool_settings']==prior['pool_settings'] and recovery['settings']==f.read(old_native['manifest'])['settings']
 assert primary['reference']==prior['reference']==scope['reference']
 ref=f.read(scope['reference']['source']);assert ref['calibration']==scope['reference']['frozen_calibration']
 assert [c['case_id'] for c in native['cases']]==IDS
 cache={}
 for task in prior['tasks']:
  path=task['reuse'] or Path(primary['manifest']['path']).parent/'points'/task['id']/'RESULT.json'
  if not isinstance(path,dict) and not path.exists():continue
  p=path if isinstance(path,dict) else f.pin(path);raw=f.read(p)
  assert r.key(task)==r.key(raw['task'])
  cache.setdefault(r.key(task),p)
 m={key:copy.deepcopy(prior[key]) for key in ('runtime','model','settings','pool_settings','reference','scope')}
 m.update(protocol_id=PROTOCOL,driver=f.pin(g.__file__),helpers=f.pin(g.RP),implementation=f.pin(f.__file__),adapter=f.pin(__file__),
          executor=f.pin(STAGE),MACE_collection=f.pin(collection_path),primary_collection=f.pin(primary_path),primary_manifest=primary['manifest'],
          cases=[],tasks=[],evaluation_kind='technical_NODE_FAIL_recovery_same_policy_q0',execution_disabled=False,
          declared_full_source_denominator=125,recovery_source_denominator=2,new_calibration=False,new_Ca_triples=False)
 identities={};reuses={}
 for cid in IDS:
  c=next(x for x in native['cases'] if x['case_id']==cid);old=next(x for x in old_native['cases'] if x['case_id']==cid)
  assert set(c['candidates'])==set(g.ORDER) and c['native_pool_status'] in ('available','unavailable')
  source=next(x for x in scope['rows'] if x.get('case_id')==cid);oldrow=next(x for x in primary['rows'] if x['case_id']==cid)
  assert oldrow['decisions']['operational']['outcome']=='unavailable' and source['execution_status']=='new_native_pool_required'
  for key in ('source_case_id','selection_id'):assert c[key]==old[key]==source[key]
  admission=c['source']['origin_row']['source'];assert admission['state_signature']==source['state_signature'] and admission['maps']==source['maps']
  assert c['source']['origin_row']==old['source']['origin_row'] and c['selection']==old['selection']
  assert c['source']['origin_row']['source_preparation']==source['source_preparation']
  assert c['technical_recovery']['manifest']==native['manifest'] and c['proposals']['Ca']==old['proposals']['Ca']
  for z in ('Ca','La'):
   for q in ('origin','adaptive_Ca'):assert c['matrix'][z][q]==old['matrix'][z][q]
  row=copy.deepcopy(source);row.update(execution_status='new_component_tasks',native_pool_status=c['native_pool_status'],native_source_record=c,
                                     matrix={z:{} for z in ('Ca','La')},paired_coordinate_max_error_A={})
  geometry={};reuses[cid]=0
  for z in ('Ca','La'):
   charge=f.read(source['state_signature'])['endpoint_charges'][z]
   for q in g.ORDER:
    cell=g.verified_cell(c['matrix'][z].get(q),z,scope);row['matrix'][z][q]=cell
    if cell is None:continue
    assert c['candidates'][q]['status']=='available' and cell['charge']==charge
    geometry[z,q]=f.coordinates(cell['xyz']);cell['inputs']={}
    for medium in ('vacuum','water'):
     key=(cell['xyz']['sha256'],cell['charge'],0,medium)
     if key not in identities:
      digest=hashlib.sha256(json.dumps([key,m['runtime'],{k:m['settings'][k] for k in r.SCIENTIFIC_SETTINGS}],sort_keys=True).encode()).hexdigest()
      task={'id':'input_'+digest,'protein':source['protein'],'context':cid,'geometry':q,'metal':z,'medium':medium,
            'xyz':cell['xyz'],'charge':charge,'multiplicity':1,'uhf':0,'valence_electrons':cell['valence_electrons'],'atom_count':cell['atom_count'],
            'reference':None,'scout':False,'reuse':cache.get(key)}
      if task['reuse']:r.verify_raw(f.read(task['reuse']),task,m|{'_base':str(base)})
      else:assert q=='adaptive_La','recovery would unexpectedly recompute an old origin/Ca candidate'
      identities[key]=task['id'];m['tasks'].append(task)
     cell['inputs'][medium]=identities[key]
     if q in ('origin','adaptive_Ca'):
      assert cache.get(key)==oldrow['components'][z+'__'+q]['inputs'][medium]
      assert oldrow['components'][z+'__'+q]['accepted'];reuses[cid]+=1
  assert reuses[cid]==8
  for q in g.ORDER:
   if all((z,q) in geometry for z in ('Ca','La')):
    ca,la=geometry['Ca',q],geometry['La',q];assert ca[0][1:]==la[0][1:] and len(ca[0])==len(la[0])
    residual=max(abs(a-b) for x,y in zip(ca[1],la[1]) for a,b in zip(x,y));assert residual<=m['pool_settings']['coordinate_copy_tolerance_A']
    row['paired_coordinate_max_error_A'][q]=residual
  for z in ('Ca','La'):
   for q in g.ORDER:
    if (z,q) in geometry:assert geometry[z,q][0]==geometry[z,'origin'][0]
  m['cases'].append(row)
 assert len(m['tasks'])<=24 and sum(t['reuse'] is None for t in m['tasks'])<=8
 base.mkdir(parents=True,exist_ok=False);f.write(base/'manifest.json',m)
 summary=g.dry_run(base);summary.update(reused_old_inputs_by_source=reuses,old_cache_inputs=len(cache),primary_unavailable_preserved=18)
 f.write(base/'PREPARATION.json',summary)
 replay=g.collect(base,'BEFORE_NEW_CALLS.json')
 for row in replay['rows']:
  old=next(x for x in primary['rows'] if x['case_id']==row['case_id'])
  assert row['pool']['origin']==old['pool']['origin'] and row['decisions']['origin']==old['decisions']['origin']
  for z in ('Ca','La'):
   for q in ('origin','adaptive_Ca'):assert row['components'][z+'__'+q]==old['components'][z+'__'+q]
 f.write(base/'REUSE_CHECK.json',{'status':'pass','old_components_exact':8,'old_raw_inputs_reused':16,'origin_scores_exact':2,'new_molecular_calls':0})
 return m

def overlay(g,primary_path,recovery_path,dest):
 f=g.f;old=f.read(primary_path);new=f.read(recovery_path);newrows={r['case_id']:r for r in new['rows']}
 assert len(old['rows'])==125 and set(newrows)==set(IDS) and old['reference']==new['reference']
 rows=[{'primary':r,'technical_recovery':newrows.get(r['case_id'])} for r in old['rows']]
 counts={p:dict(Counter((r['technical_recovery'] or r['primary'])['decisions'][p]['outcome'] for r in rows)) for p in ('operational','mathematical','origin')}
 result={'protocol_id':PROTOCOL+'_primary_preserving_overlay','primary':f.pin(primary_path),'technical_recovery':f.pin(recovery_path),'rows':rows,
         'primary_counts':old['counts'],'annotated_recovery_counts':counts,'declared_sources':125,'original_unavailable':18,
         'recovery_sources':2,'reference':old['reference'],'new_calibration':False,'production_changed':False,
         'interpretation':'Primary fields remain immutable; technical_recovery is a separate same-policy repeat after NODE_FAIL. No omitted unavailable sources.'}
 f.write(dest,result);return result

def run(args):
 s,g,placement=stage();base=args.directory.resolve();base.parent.mkdir(parents=True,exist_ok=True)
 with (base.parent/('.'+base.name+'.RECOVERY.lock')).open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);started=time.time()
  if not (base/'manifest.json').exists():prepare(g,args.scope,args.collection,args.primary,base)
  m=g.f.read(base/'manifest.json')
  assert m['adapter']==g.f.pin(__file__) and m['executor']==g.f.pin(STAGE) and m['primary_collection']==g.f.pin(args.primary)
  assert m['scope']==g.f.pin(args.scope) and m['MACE_collection']==g.f.pin(args.collection)
  assert g.f.read(base/'REUSE_CHECK.json')['status']=='pass','exact archived component replay required before new calls'
  result=s.run(args.scope,args.collection,base)
  recovery=base/'FINAL_COLLECTION.json'
  # Every invocation preserves prior collections. Resolve the collection actually returned.
  if g.f.read(recovery)!=result:
   matches=[p for p in base.glob('SLURM_ATTEMPT_*_COLLECTION.json') if g.f.read(p)==result]
   assert len(matches)==1;recovery=matches[0]
  name='PRIMARY_RECOVERY_OVERLAY.json' if not (base/'PRIMARY_RECOVERY_OVERLAY.json').exists() else recovery.stem+'_OVERLAY.json'
  combined=overlay(g,args.primary,recovery,base/name);finished=time.time()
  g.f.write(base/(recovery.stem+'_RECOVERY_RECEIPT.json'),{'adapter':g.f.pin(__file__),'collection':g.f.pin(recovery),'overlay':g.f.pin(base/name),
           'placement':placement,'elapsed_seconds_including_preparation':finished-started,'allocated_CPU_seconds_including_preparation':12*(finished-started),
           'primary_unavailable_unchanged':18,'new_backend_calls':result['cost'],'scientific_acceptance':'per-source frozen gates and common-pool availability'})
  return combined

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('scope','collection','primary','directory'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();answer=run(a);print(json.dumps({'primary':answer['primary_counts'],'recovery_annotated':answer['annotated_recovery_counts']}))
