"""Consumed three-source La/Dy electronic exchange; separate frozen-core model."""
import argparse,json,os,shutil
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,cache_key,xyz,HA_TO_KCAL
from affordable_workflow import dry_run,execute
import metal_environment_frozen_f as model
PROTOCOL='nikasha_LaDy_compact_exchange_frozencore_v1'
REPAIRED_PROTOCOL='nikasha_LaDy_compact_CHrepair_v1'
def prepared_sources(path):
 m=read_json(path)
 if m['protocol_id']!='nikasha_compact_CH_repair_preparation_v1':raise InvalidArtifact('repair protocol differs')
 sources={s['source_id']:s for s in m['sources']}
 if set(sources)!=set(SOURCES):raise InvalidArtifact('repair source matrix differs')
 for s in sources.values():
  admission=read_json(verify(s['reviewed_admission']));verify(s['source_map'])
  if admission.get('admitted') is not True or admission.get('source_id')!=s['source_id']:raise InvalidArtifact('reviewed repair admission missing')
  if not all(s['checks'].values()):raise InvalidArtifact('repair invariants failed')
  for metal in ('La','Dy'):
   new=xyz(verify(s['xyz'][metal]));old=xyz(verify(s['original_xyz'][metal]))
   indices=set(s['repaired_H_indices'])
   if len(new)!=len(old) or [x[0] for x in new]!=[x[0] for x in old]:raise InvalidArtifact('repair inventory differs')
   for i,(a,b) in enumerate(zip(old,new)):
    if i in indices:
     if a[0]!='H':raise InvalidArtifact('repair changed nonhydrogen')
    elif a!=b:raise InvalidArtifact('repair changed frozen coordinates')
 return sources
ROOT=Path(__file__).resolve().parents[1]
SOURCES={'Hans8DQ2':'prepared_v1/endpoints/Hans_EF3__','Hans8FNR':'dy_transfer_v1/endpoints/Hans_EF3__','Mex8FNS':'prepared_v1/endpoints/Mex_EF3__'}
def validate(path):
 m=read_json(path)
 if m['protocol_id'] not in (PROTOCOL,REPAIRED_PROTOCOL):raise InvalidArtifact('protocol differs')
 repaired=m['protocol_id']==REPAIRED_PROTOCOL
 prepared=prepared_sources(verify(m['prepared_inputs'])) if repaired else None
 gate=read_json(verify(m['force_gate']))
 if not (gate.get('force_consistency') or {}).get('pass_declared_tolerance'):raise InvalidArtifact('Dy force qualification missing')
 for p in m['implementation'].values():verify(p)
 expected={s+'_'+metal for s in SOURCES for metal in ('La','Dy')}
 if not repaired:expected.remove('Hans8DQ2_Dy')
 if {t['task_id'] for t in m['tasks']}!=expected:raise InvalidArtifact('finite cell matrix differs')
 for t in m['tasks']:
  a=m['assets'][t['metal']]
  if repaired and t['xyz']['sha256']!=prepared[t['source']]['xyz'][t['metal']]['sha256']:raise InvalidArtifact('repaired coordinates differ')
  if model.state(verify(t['xyz']))!=t['electronic_state']:raise InvalidArtifact('state differs')
  if verify(t['input']).read_text()!=model.input_text(verify(a['basis']).read_text(),verify(a['aux']).read_text()):raise InvalidArtifact('input differs')
  bare={k:v for k,v in t.items() if k!='cache_key'}
  if t['cache_key']!=cache_key(dict(task=bare,protocol=m['protocol_id'],implementation=m['implementation'],assets=a,orca=m['orca'])):raise InvalidArtifact('cache differs')
 return dry_run(path)
def prepare(a):
 gate=read_json(a.force_gate)
 if not (gate.get('force_consistency') or {}).get('pass_declared_tolerance'):raise InvalidArtifact('force check must pass before preparation')
 origin=model.collect(a.origin_manifest);om=read_json(a.origin_manifest)
 if origin['rows']['Dy_origin']['status']!='complete':raise InvalidArtifact('origin reuse unavailable')
 prepared=prepared_sources(a.prepared_inputs) if a.prepared_inputs else None
 protocol=REPAIRED_PROTOCOL if prepared else PROTOCOL
 w=Path(a.output).resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir();pins={}
 for p in Path(__file__).resolve().parent.glob('*.py'):shutil.copyfile(p,impl/p.name);pins[p.name]=record(impl/p.name)
 shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
 for n in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[n]=record(impl/n)
 assets={'La':dict(basis=record(a.la_basis),aux=record(a.la_aux)),'Dy':dict(basis=om['basis'],aux=om['aux'])}
 m=dict(protocol_id=protocol,assets=assets,agreement=record(a.plan),force_gate=record(a.force_gate),orca=om['orca'],implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':a.ranks,'concurrent_tasks':a.workers},tasks=[],reused_origin=origin)
 if prepared:
  m.pop('reused_origin');m['prepared_inputs']=record(a.prepared_inputs)
 for source,prefix in SOURCES.items():
  paths={metal:ROOT/'workspaces/lanm_series_followup_20260923'/(prefix+metal)/'core.xyz' for metal in ('La','Dy')}
  if prepared:paths={metal:verify(prepared[source]['xyz'][metal]) for metal in ('La','Dy')}
  at={metal:xyz(p) for metal,p in paths.items()}
  if len(at['La'])!=len(at['Dy']) or not np.array_equal([x[1:] for x in at['La']],[x[1:] for x in at['Dy']]):raise InvalidArtifact('paired coordinates differ')
  if [x[0].replace('La','Ln') for x in at['La']]!=[x[0].replace('Dy','Ln') for x in at['Dy']]:raise InvalidArtifact('paired atom identities differ')
  for metal in ('La','Dy'):
   if source=='Hans8DQ2' and metal=='Dy' and not prepared:
    if record(paths[metal])['sha256']!=om['source_xyz']['sha256']:raise InvalidArtifact('origin reuse geometry differs')
    continue
   name=source+'_'+metal;d=w/name;d.mkdir();shutil.copyfile(paths[metal],d/'core.xyz');b=assets[metal]
   (d/'endpoint.inp').write_text(model.input_text(verify(b['basis']).read_text(),verify(b['aux']).read_text()))
   t=dict(task_id=name,source=source,metal=metal,charge=-1,multiplicity=1,electronic_state=model.state(d/'core.xyz'),xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient')
   t['cache_key']=cache_key(dict(task=t,protocol=protocol,implementation=pins,assets=b,orca=m['orca']));m['tasks'].append(t)
 write_new(w/'manifest.json',m);return validate(w/'manifest.json')
def collect(path):
 from ggr_sensitivity import executed
 validate(path);m,rows=executed(path)
 for t in m['tasks']:
  r=rows[t['task_id']]
  if r['status']=='complete':
   try:r.update(model.parse(t,{'electronic_state':t['electronic_state']}))
   except (ValueError,OSError) as e:r.update(status='invalid',reason=str(e),energy_hartree=None)
 if 'reused_origin' in m:rows['Hans8DQ2_Dy']=m['reused_origin']['rows']['Dy_origin']
 contrasts={}
 for s in ('Hans8DQ2','Hans8FNR'):
  names=[s+'_Dy',s+'_La','Mex8FNS_Dy','Mex8FNS_La']
  contrasts[s]=sum(c*rows[n]['energy_hartree'] for c,n in zip((1,-1,-1,1),names))*HA_TO_KCAL if all(rows[n]['status']=='complete' for n in names) else None
 return dict(protocol_id=m['protocol_id'],rows=rows,conditional_exchange_kcal_mol=contrasts,sign='positive supports greater relative La preference in Hans',biological_evidence='protein-level qualitative; no EF3 labels',affinity=None,full_hybrid_qualified=False,collector=record(__file__))
def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('prepare')
 for k in ('origin-manifest','force-gate','la-basis','la-aux','plan','output'):a.add_argument('--'+k,required=True)
 a.add_argument('--prepared-inputs')
 a.add_argument('--ranks',type=int,required=True);a.add_argument('--workers',type=int,required=True)
 for op in ('dry-run','execute','collect'):
  a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
 a=p.parse_args()
 if a.op=='prepare':r=prepare(a)
 elif a.op=='collect':r=collect(a.manifest)
 elif a.op=='dry-run':r=validate(a.manifest)
 else:
  validate(a.manifest);m=read_json(a.manifest);os.environ['METAL_ENV_WORKERS']=str(m['execution_resources']['concurrent_tasks']);r=execute(a.manifest)
 if a.op!='prepare' and a.output:write_new(a.output,r)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
