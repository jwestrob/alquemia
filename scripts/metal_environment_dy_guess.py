"""Two declared native initial guesses on one real isolated Dy core; research only."""
import argparse,json,os,shutil
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,cache_key
from affordable_workflow import dry_run,execute
from metal_environment_reference import ORCA
from metal_environment_electronic_state import METHOD,PBE0_PROFILE,PBE0_METHOD,describe,input_text,parse
PROTOCOL='nikasha_isolated_Dy_native_guess_diagnostic_v1'

def validate(path):
 m=read_json(path)
 profile=m.get('method_profile'); expected=PROTOCOL if profile is None else 'nikasha_isolated_Dy_pbe0_capability_v1'
 guesses=m.get('declared_guesses',['PModel','HCore'])
 if profile not in (None,PBE0_PROFILE) or not guesses or len(set(guesses))!=len(guesses) or set(guesses)-{'HCore','PModel'}:raise InvalidArtifact('unsupported method/guess inventory')
 if m['protocol_id']!=expected or len(m['tasks'])!=len(guesses) or {t['scf_guess'] for t in m['tasks']}!=set(guesses):raise InvalidArtifact('declared guess inventory differs')
 verify(m['source_xyz'])
 for p in m['implementation'].values():verify(p)
 for t in m['tasks']:
  if t['xyz']['sha256']!=m['source_xyz']['sha256'] or t['metal']!='Dy' or t['charge']!=-1 or t['multiplicity']!=6:raise InvalidArtifact('frozen Dy state differs')
  describe(verify(t['xyz']),'Dy',-1,6)
  if verify(t['input']).read_text()!=input_text(-1,6,False,t['scf_guess'],profile):raise InvalidArtifact('guess or Hamiltonian differs')
  bare={k:v for k,v in t.items() if k!='cache_key'}
  if t['cache_key']!=cache_key({'task':bare,'source':m['source_xyz'],'protocol':m['protocol_id'],'implementation':m['implementation'],'orca':m['orca']}):raise InvalidArtifact('cache differs')
 return dry_run(path)

def prepare(source,plan,output,ranks,profile=None,guesses=None):
 guesses=guesses or ['PModel','HCore']; protocol=PROTOCOL if profile is None else 'nikasha_isolated_Dy_pbe0_capability_v1'
 describe(source,'Dy',-1,6);w=Path(output).resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir();pins={}
 for p in Path(__file__).resolve().parent.glob('*.py'):
  shutil.copyfile(p,impl/p.name);pins[p.name]=record(impl/p.name)
 shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
 for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
 m=dict(protocol_id=protocol,method=METHOD if profile is None else PBE0_METHOD,method_profile=profile,declared_guesses=guesses,source_xyz=record(source),agreement=record(plan),orca=record(ORCA),implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':ranks,'concurrent_tasks':len(guesses)},tasks=[],compute_budget=None,wall_time_limit=None)
 for guess in guesses:
  d=w/guess;d.mkdir();shutil.copyfile(source,d/'core.xyz');(d/'endpoint.inp').write_text(input_text(-1,6,False,guess,profile))
  t=dict(task_id=guess,metal='Dy',charge=-1,multiplicity=6,scf_guess=guess,method_profile=profile,xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient')
  t['cache_key']=cache_key({'task':t,'source':m['source_xyz'],'protocol':protocol,'implementation':pins,'orca':m['orca']});m['tasks'].append(t)
 write_new(w/'manifest.json',m);return validate(w/'manifest.json')

def collect(path):
 from ggr_sensitivity import executed
 validate(path);m,rows=executed(path)
 for t in m['tasks']:
  row=rows[t['task_id']]
  if row['status']!='complete':continue
  try:row.update(parse(t,verify(row['output']),t['engrad_path'],embedded=False))
  except (ValueError,OSError) as e:row.update(status='invalid',reason=str(e),energy_hartree=None)
 return dict(protocol_id=m['protocol_id'],manifest=record(path),rows=rows,status='complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete',classification=None,affinity=None,full_hybrid_qualified=False)

def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('prepare')
 for k in ('source','plan','output'):a.add_argument('--'+k,required=True)
 a.add_argument('--ranks',type=int,required=True);a.add_argument('--method-profile',choices=[PBE0_PROFILE]);a.add_argument('--guesses',nargs='+',choices=['PModel','HCore'])
 for op in ('dry-run','execute','collect'):
  a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
 a=p.parse_args()
 if a.op=='prepare':r=prepare(a.source,a.plan,a.output,a.ranks,a.method_profile,a.guesses)
 elif a.op=='dry-run':r=validate(a.manifest)
 elif a.op=='collect':r=collect(a.manifest)
 else:
  validate(a.manifest);m=read_json(a.manifest);n=m['execution_resources']['mpi_ranks']
  if int(os.environ.get('SLURM_NTASKS','0'))<len(m['tasks'])*n:raise InvalidArtifact('MPI slots insufficient')
  os.environ['METAL_ENV_WORKERS']=str(len(m['tasks']));r=execute(a.manifest)
 if a.op!='prepare' and a.output:write_new(a.output,r)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
