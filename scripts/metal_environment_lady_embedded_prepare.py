"""Preparation-only four-cell frozen-f La/Dy embedded electronic response.

No execute/collect operation: reference qualification and integration belong to
its owning runner. Existing physical-spin preparation and fields are immutable.
"""
import argparse,json,re,shutil
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,cache_key,xyz
from affordable_workflow import dry_run
from metal_environment_lanm_reference import check_config
from metal_environment_reference import ORCA
from metal_environment_electronic_state import NUMBERS
PROTOCOL='nikasha_lady_frozen4f_embedded_response_preparation_v1'
METHOD='PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad'

def state(path,metal,charge):
 atoms=xyz(path);syms=[a[0] for a in atoms]
 if type(charge) is not int or metal not in ('La','Dy') or syms.count(metal)!=1 or sum(s in ('La','Dy') for s in syms)!=1 or set(syms)-{'La','Dy','C','H','N','O','S'}:raise InvalidArtifact('unsupported target/ligands/charge')
 core={'La':46,'Dy':55}[metal];total=sum(NUMBERS[s] for s in syms)-charge;explicit=total-core
 if explicit<0 or explicit%2:raise InvalidArtifact('restricted valence parity differs')
 return dict(metal=metal,oxidation_state_hypothesis=3,charge=charge,physical_multiplicity=1 if metal=='La' else 6,physical_f_occupation=0 if metal=='La' else 9,effective_multiplicity=1,ecp_core_electrons=core,all_electron_count=total,explicit_electrons=explicit,alpha_electrons=explicit//2,beta_electrons=explicit//2,spin_orbit_included=False,representation='spin_free_frozen_f_valence_model',physical_state_qualified=False)

def asset_text(metal,basis,aux):
 b=verify(basis).read_text();a=verify(aux).read_text();core={'La':46,'Dy':55}[metal]
 if re.findall(r'NewGTO\s+(\w+)',b,re.I)!=[metal] or re.findall(r'NewECP\s+(\w+)',b,re.I)!=[metal] or re.findall(r'N_core\s+(\d+)',b,re.I)!=[str(core)] or re.findall(r'NewAuxJGTO\s+(\w+)',a,re.I)!=[metal]:raise InvalidArtifact('literal orbital/ECP/AuxJ asset identity differs')
 return b,a

def input_text(charge,basis,aux,trah=False,seeded=False):
 guess='Guess MORead\n MOInp "initial.gbw"' if seeded else 'Guess PModel'
 return f'! {METHOD}{" TRAH" if trah else ""}\n%scf\n {guess}\nend\n%method\n DoEQ false\nend\n%basis\n{basis}\n{aux}\nend\n%pointcharges "environment.pc"\n* xyzfile {charge} 1 core.xyz\n'

def check_seed(seed_pin,target,method,orca,assets):
 from metal_environment_orbital_seed import validate_seed,task
 s=read_json(verify(seed_pin));original=verify(s['target_manifest'])
 checked=validate_seed(verify(s['source_manifest']),verify(s['source_collection']),s['source_task_id'],original,s['target_task_id'])
 if any(s.get(k)!=v for k,v in checked.items()):raise InvalidArtifact('seed validation record differs')
 verify(s['staged_gbw'])
 if s['staged_gbw']['sha256']!=s['source_gbw']['sha256']:raise InvalidArtifact('staged orbital copy differs')
 m=read_json(original);t=task(m,target['task_id'])
 if m['method']!=method or m['orca']!=orca:raise InvalidArtifact('seed target method differs')
 for k in ('xyz','pointcharges','core_mapping','boundary_mapping'):
  if t[k]['sha256']!=target[k]['sha256']:raise InvalidArtifact('seed target geometry or mapping differs')
 if t['electronic_state']!=target['electronic_state']:raise InvalidArtifact('seed target state differs')
 for kind in ('basis','aux'):
  if m['assets'][t['metal']][kind]['sha256']!=assets[t['metal']][kind]['sha256']:raise InvalidArtifact('seed target basis differs')
 return s

def validate(path):
 m=read_json(path)
 trah=m.get('solver')=='TRAH'
 origins=m.get('origin_metals')
 if origins is not None and (trah or origins not in [['La'],['Dy'],['La','Dy']]):raise InvalidArtifact('invalid declared origin subset')
 protocol=PROTOCOL+'_origin_subset_v1' if origins is not None else PROTOCOL+'_trah_scout_v1' if trah else PROTOCOL
 seeded=m.get('initial_guess')=='MORead'
 if m.get('initial_guess','PModel') not in ('PModel','MORead'):raise InvalidArtifact('unknown initial guess')
 if seeded:protocol+='_moread_v1'
 if m.get('solver','default') not in ('default','TRAH') or m['protocol_id']!=protocol or m['method']!=METHOD:raise InvalidArtifact('method identity differs')
 c=check_config(read_json(verify(m['inputs'])))
 if m['source_id']!=c['source_id']:raise InvalidArtifact('source differs')
 for p in m['implementation'].values():verify(p)
 verify(m['agreement']);verify(m['orca'])
 expected={(z,'A') for z in origins} if origins is not None else {('Dy','A')} if trah else {(z,q) for z in ('La','Dy') for q in ('A','B')}
 if len(m['tasks'])!=len(expected) or {(t['metal'],t['configuration']) for t in m['tasks']}!=expected:raise InvalidArtifact('finite matrix differs')
 if c['configurations']['A']['pointcharges']['sha256']!=c['configurations']['B']['pointcharges']['sha256']:raise InvalidArtifact('fixed field changed between A/B')
 for t in m['tasks']:
  conf=c['configurations'][t['configuration']];e=conf['endpoints'][t['metal']]
  if t['xyz']['sha256']!=e['xyz']['sha256'] or t['pointcharges']['sha256']!=conf['pointcharges']['sha256'] or t['charge']!=e['charge'] or t['multiplicity']!=1 or t['physical_multiplicity']!=e['multiplicity']:raise InvalidArtifact('source geometry/field/state differs')
  verify(t['pointcharges']);st=state(verify(t['xyz']),t['metal'],t['charge'])
  if st!=t['electronic_state']:raise InvalidArtifact('state metadata differs')
  b,a=asset_text(t['metal'],m['assets'][t['metal']]['basis'],m['assets'][t['metal']]['aux'])
  if seeded:
   s=check_seed(t['orbital_seed'],t,m['method'],m['orca'],m['assets'])
   if verify(t['initial_gbw']).name!='initial.gbw' or t['initial_gbw']['sha256']!=s['source_gbw']['sha256']:raise InvalidArtifact('task initial orbital differs')
  elif 'orbital_seed' in t or 'initial_gbw' in t:raise InvalidArtifact('undeclared orbital seed')
  if verify(t['input']).read_text()!=input_text(t['charge'],b,a,trah,seeded):raise InvalidArtifact('literal input differs')
  bare={k:v for k,v in t.items() if k!='cache_key'}
  if t['cache_key']!=cache_key(dict(task=bare,protocol=m['protocol_id'],method=METHOD,inputs=m['inputs'],assets=m['assets'],orca=m['orca'],implementation=m['implementation'])):raise InvalidArtifact('cache identity differs')
 return dry_run(path)

def prepare(inputs,plan,output,assets,mpi_ranks,workers,trah_scout=False,origin_metals=None,seed_records=None):
 if origin_metals is not None and (trah_scout or origin_metals not in [['La'],['Dy'],['La','Dy']] or workers>len(origin_metals)):raise InvalidArtifact('invalid origin subset or excess workers')
 if trah_scout and workers!=1:raise InvalidArtifact('TRAH scout has one endpoint')
 protocol=PROTOCOL+'_origin_subset_v1' if origin_metals is not None else PROTOCOL+'_trah_scout_v1' if trah_scout else PROTOCOL
 if seed_records is not None:protocol+='_moread_v1'
 if type(mpi_ranks) is not int or mpi_ranks<1 or type(workers) is not int or not 1<=workers<=4:raise InvalidArtifact('invalid four-cell layout')
 c=check_config(read_json(inputs))
 if c['configurations']['A']['pointcharges']['sha256']!=c['configurations']['B']['pointcharges']['sha256']:raise InvalidArtifact('fixed field must be identical')
 for z in ('La','Dy'):asset_text(z,assets[z]['basis'],assets[z]['aux'])
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
 for p in sorted(Path(__file__).resolve().parent.glob('*.py')):shutil.copyfile(p,impl/p.name);pins[p.name]=record(impl/p.name)
 shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
 for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
 pinned={}
 for z in ('La','Dy'):
  pinned[z]={}
  for kind,pin in assets[z].items():dest=out/f'{z}_{kind}.inc';shutil.copyfile(verify(pin),dest);pinned[z][kind]=record(dest)
 m=dict(protocol_id=protocol,method=METHOD,source_id=c['source_id'],inputs=record(inputs),agreement=record(plan),assets=pinned,asset_origins=assets,orca=record(ORCA),implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},tasks=[],compute_budget=None,wall_time_limit=None,execution_status='prepared_only_requires_owner_integration_and_gate',energy_scope='finite_embedded_electronic_component_only',full_hybrid_qualified=False,classification=None,affinity=None)
 if trah_scout:m['solver']='TRAH'
 if origin_metals is not None:m['origin_metals']=origin_metals
 if seed_records is not None:m['initial_guess']='MORead'
 for label in ('A','B'):
  for z in ('La','Dy'):
   if trah_scout and (z,label)!=('Dy','A'):continue
   if origin_metals is not None and (label!='A' or z not in origin_metals):continue
   e=c['configurations'][label]['endpoints'][z];d=out/f'{z}_{label}';d.mkdir();shutil.copyfile(verify(e['xyz']),d/'core.xyz');shutil.copyfile(verify(c['configurations'][label]['pointcharges']),d/'environment.pc');b,a=asset_text(z,pinned[z]['basis'],pinned[z]['aux']);(d/'endpoint.inp').write_text(input_text(e['charge'],b,a,trah_scout));st=state(d/'core.xyz',z,e['charge'])
   t=dict(task_id=d.name,metal=z,configuration=label,charge=e['charge'],multiplicity=1,physical_multiplicity=e['multiplicity'],electronic_state=st,xyz=record(d/'core.xyz'),pointcharges=record(d/'environment.pc'),core_mapping=c['core_mapping'][label],boundary_mapping=c['boundary_mapping'],input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient')
   if seed_records is not None:
    if t['task_id'] not in seed_records:raise InvalidArtifact('required seed missing')
    t['orbital_seed']=record(seed_records[t['task_id']]);s=check_seed(t['orbital_seed'],t,METHOD,m['orca'],pinned)
    shutil.copyfile(verify(s['staged_gbw']),d/'initial.gbw');t['initial_gbw']=record(d/'initial.gbw')
    (d/'endpoint.inp').write_text(input_text(e['charge'],b,a,trah_scout,True));t['input']=record(d/'endpoint.inp')
   t['cache_key']=cache_key(dict(task=t,protocol=protocol,method=METHOD,inputs=m['inputs'],assets=pinned,orca=m['orca'],implementation=pins));m['tasks'].append(t)
 if seed_records is not None and set(seed_records)!={t['task_id'] for t in m['tasks']}:raise InvalidArtifact('extra seed tasks')
 write_new(out/'manifest.json',m);return validate(out/'manifest.json')

def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='op',required=True);a=sp.add_parser('prepare')
 for key in ('inputs','plan','output','la-basis','dy-basis','la-aux','dy-aux'):a.add_argument('--'+key,required=True)
 a.add_argument('--trah-scout',action='store_true')
 a.add_argument('--origin-metals',nargs='+',choices=['La','Dy'],help='Prepare only the declared A origins; distinct protocol, no automatic reuse or submission')
 a.add_argument('--seed-records',help='JSON mapping every task ID to a staged SEED.json path; explicit MORead only')
 a.add_argument('--mpi-ranks',required=True,type=int);a.add_argument('--workers',required=True,type=int);a=sp.add_parser('dry-run');a.add_argument('--manifest',required=True);a=p.parse_args()
 if a.op=='dry-run':result=validate(a.manifest)
 else:result=prepare(a.inputs,a.plan,a.output,{z:{'basis':record(getattr(a,z.lower()+'_basis')),'aux':record(getattr(a,z.lower()+'_aux'))} for z in ('La','Dy')},a.mpi_ranks,a.workers,a.trah_scout,a.origin_metals,read_json(a.seed_records) if a.seed_records else None)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
