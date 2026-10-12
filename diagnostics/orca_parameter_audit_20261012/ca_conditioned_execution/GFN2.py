"""Qualified solvent components for the frozen Ca-source transfer; no new fit."""
from pathlib import Path
import argparse,copy,fcntl,hashlib,importlib.util,json,math,sys,time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import context_factorial_gfn2 as f
RP=ROOT/'diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py'
spec=importlib.util.spec_from_file_location('expanded_pool_component_helpers',RP);r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ORDER=['origin','adaptive_Ca','adaptive_La']

def verified_cell(cell,z,scope):
 if cell is None or cell.get('status') not in ('complete','available'):return None
 for field in ('xyz','MACE'):assert f.pin(cell[field]['path'])==cell[field]
 receipt=f.read(cell['MACE']);assert receipt['status'] in ('complete','computed')
 state=receipt.get('state',receipt.get('input_state_check'))
 charge=cell['charge'];multiplicity=cell['multiplicity']
 assert charge==int(charge) and multiplicity==1
 assert state['charge']==charge and state['spin_multiplicity']==multiplicity and state['head']=='omol'
 if receipt.get('request'):
  assert receipt['model']==scope['model'];request=f.read(receipt['request'])
 else:
  m=f.read(receipt['manifest']);assert m['model']==scope['model']
  request=next(t for t in m['tasks'] if t['task_id']==receipt['task_id'])
 assert request['xyz']['sha256']==cell['xyz']['sha256'] and request['charge']==charge
 assert request.get('multiplicity',request.get('spin_multiplicity'))==multiplicity
 symbols,positions=f.coordinates(cell['xyz'])
 assert symbols[0]==z and state['atoms']==len(symbols)
 assert receipt['energy_eV']==cell['energy_eV'] and math.isfinite(cell['energy_eV'])
 electrons=sum(f.VALENCE[a] for a in symbols)-charge;assert electrons%2==0
 if receipt.get('forces'):assert f.pin(receipt['forces']['path'])==receipt['forces']
 return {'xyz':cell['xyz'],'MACE':cell['MACE'],'MACE_eV':cell['energy_eV'],'charge':charge,'multiplicity':multiplicity,
         'atom_count':len(symbols),'valence_electrons':electrons}

def verify_prior(source,manifest):
 collection=f.read(source['prior_collection'])
 old=next(row for row in collection['rows'] if row['case_id']==source['case_id'])
 assert old==source['prior_row']
 for component in old['components'].values():
  assert component['accepted']
  for medium,pin in component['inputs'].items():
   raw=f.read(pin);r.verify_raw(raw,raw['task']|{'reuse':pin},manifest)
   assert raw['backends']['xtb']['energy_hartree']==component['primary_'+medium+'_energy_hartree']
 return old

def prepare(scope_path,collection_path,output,fixture=False):
 scope=f.read(scope_path);data=f.read(collection_path);assert data['terminal'] is True,'MACE collection must be terminal before solvent stage'
 assert 'manifest' in data;native_manifest=f.read(data['manifest'])
 pool_settings=native_manifest['settings'] if fixture else native_manifest['pool_settings']
 assert pool_settings['candidate_order']==ORDER and pool_settings['origin_selection_tolerance_kcal_mol']==.1
 assert pool_settings['coordinate_copy_tolerance_A']==1e-12 and pool_settings['equivalent_cell_tolerance_kcal_mol']==.1
 if not fixture:
  assert native_manifest['scientific_scope']==f.pin(scope_path) and native_manifest['model']==scope['model']
 ref=f.read(scope['reference']['source']);assert ref['calibration']==scope['reference']['frozen_calibration']
 source_manifest=f.read(ref['manifest']);assert source_manifest['model']==scope['model']
 assert scope['counts']['declared_Ca_sources']==125 and len(scope['rows'])==125
 incoming={c['case_id']:c for c in data['cases']};assert len(incoming)==len(data['cases'])
 cache={}
 for old in scope['rows']:
  if old['execution_status']!='reuse_actual_corrected_pool':continue
  for component in old['prior_row']['components'].values():
   assert component['accepted']
   for p in component['inputs'].values():
    raw=f.read(p);cache[r.key(raw['task'])]=p
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
 m={'protocol_id':'Nikasha_corrected_expanded_Ca_GFN2_transfer_v1','driver':f.pin(__file__),'helpers':f.pin(RP),'implementation':f.pin(f.__file__),
    'scope':f.pin(scope_path),'MACE_collection':f.pin(collection_path),'runtime':source_manifest['runtime'],'model':scope['model'],
    'settings':source_manifest['settings'],'pool_settings':pool_settings,'reference':scope['reference'],'tasks':[],'cases':[],
    'evaluation_kind':'archived_parser_replay_not_new_science' if fixture else 'frozen_Ca_source_transfer',
    'execution_disabled':bool(fixture),'declared_full_source_denominator':125,'new_calibration':False,'new_Ca_triples':False}
 identities={}
 for source in scope['rows']:
  if fixture and source.get('case_id') not in incoming:continue
  row=copy.deepcopy(source)
  if source['execution_status']=='excluded':m['cases'].append(row);continue
  if source['execution_status']=='reuse_actual_corrected_pool' and not fixture:m['cases'].append(row);continue
  cid=source['case_id'];c=incoming[cid]
  for field in ('source_case_id','selection_id'):assert c[field]==source[field]
  assert set(c['candidates'])==set(ORDER) and c['native_pool_status'] in ('available','unavailable')
  admission=c['source']['origin_row']['source']
  assert admission['state_signature']==source['state_signature'] and admission['maps']==source['maps']
  assert c['source']['origin_row']['source_preparation']==source['source_preparation']
  row['execution_status']='new_component_tasks';row['native_pool_status']=c['native_pool_status'];row['native_source_record']=c
  row['matrix']={'Ca':{},'La':{}};symbols_by_cell={};row['paired_coordinate_max_error_A']={}
  for z in ('Ca','La'):
   source_origin=f.read(source['state_signature'])['endpoint_charges'][z]
   for q in ORDER:
    cell=verified_cell(c['matrix'][z].get(q),z,scope);row['matrix'][z][q]=cell
    if cell is None:continue
    assert c['candidates'][q]['status']=='available'
    assert cell['charge']==source_origin
    symbols,coords=f.coordinates(cell['xyz']);symbols_by_cell[z,q]=(symbols,coords)
    cell['inputs']={}
    for medium in ('vacuum','water'):
     key=(cell['xyz']['sha256'],cell['charge'],0,medium)
     if key not in identities:
      digest=hashlib.sha256(json.dumps([key,m['runtime'],{k:m['settings'][k] for k in r.SCIENTIFIC_SETTINGS}],sort_keys=True).encode()).hexdigest()
      t={'id':'input_'+digest,'protein':source['protein'],'context':cid,'geometry':q,'metal':z,'medium':medium,
        'xyz':cell['xyz'],'charge':cell['charge'],'multiplicity':1,'uhf':0,'valence_electrons':cell['valence_electrons'],
        'atom_count':cell['atom_count'],'reference':None,'scout':False,'reuse':cache.get(key)}
      identities[key]=t['id'];m['tasks'].append(t)
     cell['inputs'][medium]=identities[key]
  for q in ORDER:
   if all((z,q) in symbols_by_cell for z in ('Ca','La')):
    ca,la=[symbols_by_cell[z,q] for z in ('Ca','La')];assert ca[0][1:]==la[0][1:]
    residual=max(abs(a-b) for ca_atom,la_atom in zip(ca[1],la[1]) for a,b in zip(ca_atom,la_atom))
    assert residual<=pool_settings['coordinate_copy_tolerance_A']
    row['paired_coordinate_max_error_A'][q]=residual
  for z in ('Ca','La'):
   if (z,'origin') not in symbols_by_cell:continue
   for q in ORDER:
    if (z,q) in symbols_by_cell:assert symbols_by_cell[z,q][0]==symbols_by_cell[z,'origin'][0]
  m['cases'].append(row)
 if not fixture:
  assert len(m['cases'])==125 and len(incoming)==106
  assert set(incoming)=={c['case_id'] for c in scope['rows'] if c['execution_status']=='new_native_pool_required'}
 assert len({t['id'] for t in m['tasks']})==len(m['tasks']) and len(m['tasks'])<=1272
 f.write(out/'manifest.json',m);summary=dry_run(out);f.write(out/'PREPARATION.json',summary)
 return summary

def dry_run(base):
 m=f.read(base/'manifest.json');assert f.pin(__file__)==m['driver'] and f.pin(RP)==m['helpers'] and f.pin(f.__file__)==m['implementation']
 scope=f.read(m['scope']);data=f.read(m['MACE_collection']);assert data['terminal'] is True
 runtime=f.read(m['runtime'])
 for key in ('xtb','xtb_parameter_file','tblite_python'):assert f.pin(runtime[key]['path'])==runtime[key]
 for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:assert f.pin(p['path'])==p
 for item in runtime['tblite_dynamic'].values():
  for p in item['libraries']:assert f.pin(p['path'])==p
 for t in m['tasks']:
  f.coordinates(t['xyz'])
  if t['reuse']:r.verify_raw(f.read(t['reuse']),t,m|{'_base':str(base)})
 for case in m['cases']:
  if case['execution_status']=='reuse_actual_corrected_pool':verify_prior(case,m)
 assert m['reference']==scope['reference']
 summary={'status':'prepared_no_new_scientific_calls','manifest':f.pin(base/'manifest.json'),'evaluation_kind':m['evaluation_kind'],
  'source_rows':len(m['cases']),'declared_full_source_denominator':125,'unique_inputs_per_backend':len(m['tasks']),
  'actual_reused_inputs_per_backend':sum(t['reuse'] is not None for t in m['tasks']),
  'fresh_inputs_per_backend':sum(t['reuse'] is None for t in m['tasks']),'frozen_calibration':m['reference']['frozen_calibration'],
  'execution_disabled':m['execution_disabled'],'production_changed':False}
 print(json.dumps(summary,indent=2));return summary

def decision(value,bands,expected):
 if value is None:return {'decision':None,'outcome':'unavailable'}
 predicted='Ca' if value<=bands['Ca_max'] else 'La' if value>=bands['La_min'] else 'inconclusive'
 return {'decision':predicted,'outcome':'inconclusive' if predicted=='inconclusive' else 'correct' if predicted==expected else 'wrong'}

def collect(base,name):
 m=f.read(base/'manifest.json');tasks={t['id']:t for t in m['tasks']};rows=[]
 for c in m['cases']:
  row={k:c.get(k) for k in ('source_case_id','case_id','selection_id','protein','biological_group','expected_class','source_conditioning_metal','preparation_status','preparation_reason','execution_status')}
  if c['execution_status']=='reuse_actual_corrected_pool':
   old=verify_prior(c,m);row.update(pool=old['pool'],components=old['components'],prior_collection=c['prior_collection'])
  elif c['execution_status']=='excluded':row.update(pool=r.choose({z:{q:None for q in ORDER} for z in ('Ca','La')}),components={})
  else:
   matrix={z:{} for z in ('Ca','La')};components={}
   for z in ('Ca','La'):
    for q in ORDER:
     cell=c['matrix'][z][q]
     if cell is None:matrix[z][q]=None;components[z+'__'+q]={'accepted':False,'reason':'unavailable_MACE_cell'};continue
     component=r.component(base,*[tasks[cell['inputs'][medium]] for medium in ('vacuum','water')]);components[z+'__'+q]=component
     matrix[z][q]={'MACE_eV':cell['MACE_eV'],'GFN2_vacuum_hartree':component['primary_vacuum_energy_hartree'],
                  'GFN2_ALPB_hartree':component['primary_water_energy_hartree']} if component['accepted'] else None
   pool=r.choose(matrix)
   if c['native_pool_status']!='available':pool.update(status='unavailable',operational=None,mathematical=None)
   row.update(pool=pool,components=components,native_pool_status=c['native_pool_status'],paired_coordinate_max_error_A=c['paired_coordinate_max_error_A'])
  row['decisions']={policy:decision(row['pool'][policy]['R_model_kcal_mol'] if row['pool'][policy] is not None else None,
       ref['bands'],c['expected_class']) for policy,ref in m['reference']['frozen_calibration'].items()}
  rows.append(row)
 usage={'new_primary_calls':0,'new_independent_calls':0,'CPU_seconds':0.,'summed_call_wall_seconds':0.,'maximum_process_RSS_KiB':0,
        'started_calls_without_execution_receipt':[],'execution_receipts':[]}
 for t in m['tasks']:
  if t['reuse']:continue
  for backend in ('xtb','tblite'):
   directory=base/'points'/t['id']/backend;receipt=directory/'EXECUTION.json'
   if not receipt.exists():
    if (directory/'START.json').exists():usage['started_calls_without_execution_receipt'].append(f.pin(directory/'START.json'))
    continue
   rec=f.read(receipt);usage['execution_receipts'].append(f.pin(receipt))
   usage['new_primary_calls' if backend=='xtb' else 'new_independent_calls']+=1
   usage['CPU_seconds']+=rec['cpu_user_seconds']+rec['cpu_system_seconds'];usage['summed_call_wall_seconds']+=rec['wall_seconds']
   usage['maximum_process_RSS_KiB']=max(usage['maximum_process_RSS_KiB'],rec['peak_RSS_KiB'])
 common=[row for row in rows if all(row['pool'][p] is not None for p in ('operational','origin'))]
 usage['cost_coverage']='incomplete_interrupted_calls' if usage['started_calls_without_execution_receipt'] else 'all_started_backend_calls_receipted'
 result={'protocol_id':m['protocol_id'],'manifest':f.pin(base/'manifest.json'),'evaluation_kind':m['evaluation_kind'],
  'rows':rows,'counts':{p:dict(Counter(row['decisions'][p]['outcome'] for row in rows)) for p in ('operational','mathematical','origin')},
  'declared_sources':len(rows),'declared_biological_groups':len({x['biological_group'] for x in rows}),
  'supported_preparation_groups':len({x['biological_group'] for x in rows if x['preparation_status']=='prepared'}),
  'common_origin_operational_coverage':{'sources':len(common),'biological_groups':len({x['biological_group'] for x in common}),
     'counts':{p:dict(Counter(row['decisions'][p]['outcome'] for row in common)) for p in ('operational','origin')}},
  'reference':m['reference'],'cost':usage,'new_calibration':False,'new_Ca_triples':False,'production_changed':False}
 f.write(base/name,result);print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2));return result

def execute(base,workers):
 assert workers==14,'this finite stage is authorized for14 one-thread workers'
 with (base/'EXECUTION.lock').open('a') as lock:
  fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB);dry_run(base);m=f.read(base/'manifest.json')
  assert not m['execution_disabled'],'archived parser replay cannot launch science'
  runtime=f.read(m['runtime']);pending,complete,interrupted=r.pending_tasks(base,[t for t in m['tasks'] if t['reuse'] is None])
  attempt=1+len(list(base.glob('ATTEMPT_*_START.json')));stem=f'ATTEMPT_{attempt:03d}'
  f.write(base/(stem+'_START.json'),{'started_epoch':time.time(),'workers':14,'threads_per_worker':1,'pending_inputs':len(pending),'completed_reused':complete,'interrupted_unavailable':interrupted})
  def run(t):
   try:f.execute_point(t,m,base,runtime)
   except Exception as exc:return {'id':t['id'],'reason':str(exc)}
   return None
  with ThreadPoolExecutor(max_workers=workers) as pool:errors=[e for e in pool.map(run,pending) if e is not None]
  f.write(base/(stem+'_FINISH.json'),{'finished_epoch':time.time(),'exceptions':errors})
  name='FINAL_COLLECTION.json' if not (base/'FINAL_COLLECTION.json').exists() else stem+'_COLLECTION.json'
  result=collect(base,name)
  if errors:raise RuntimeError(f'{len(errors)} task exceptions; retained in {stem}_FINISH.json')
  return result

if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True)
 q=s.add_parser('prepare');q.add_argument('--scope',required=True);q.add_argument('--collection',required=True);q.add_argument('--output',required=True);q.add_argument('--archive-parser-fixture',action='store_true')
 for op in ('dry-run','execute','collect'):
  q=s.add_parser(op);q.add_argument('--directory',required=True)
  if op=='execute':q.add_argument('--workers',type=int,required=True)
  if op=='collect':q.add_argument('--output-name',required=True)
 a=p.parse_args()
 if a.op=='prepare':prepare(a.scope,a.collection,a.output,a.archive_parser_fixture)
 elif a.op=='dry-run':dry_run(Path(a.directory).resolve())
 elif a.op=='execute':execute(Path(a.directory).resolve(),a.workers)
 else:collect(Path(a.directory).resolve(),a.output_name)
