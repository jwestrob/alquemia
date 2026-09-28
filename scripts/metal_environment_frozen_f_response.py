"""Four declared physical displacements qualify the frozen-f origin gradient."""
import argparse,json,os,shutil
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,cache_key,xyz,HA_TO_KCAL
from affordable_workflow import dry_run,execute
import metal_environment_frozen_f as scout
PROTOCOL='nikasha_DyIII_frozen4f_directional_force_v1'
STEPS={'plus':.005,'minus':-.005,'half_plus':.0025,'half_minus':-.0025}
def direction(source):
 atoms=xyz(source);q=np.array([a[1:] for a in atoms]);metals=[i for i,a in enumerate(atoms) if a[0]=='Dy']
 if len(metals)!=1:raise InvalidArtifact('one Dy required')
 i=metals[0];donors=[j for j,a in enumerate(atoms) if a[0]=='O']
 if not donors:raise InvalidArtifact('real oxygen donor required')
 j=min(donors,key=lambda k:(np.linalg.norm(q[k]-q[i]),k));v=q[j]-q[i];v/=np.linalg.norm(v)
 return dict(metal_index=i,donor_index=j,unit_vector=v.tolist(),rule='nearest existing oxygen; tie by source atom index; metal motion only')
def validate(path):
 m=read_json(path)
 if m['protocol_id']!=PROTOCOL or {t['task_id'] for t in m['tasks']}!=set(STEPS):raise InvalidArtifact('four displacements required')
 for v in m['implementation'].values():verify(v)
 origin=verify(m['source_xyz']);atoms=xyz(origin);q=np.array([a[1:] for a in atoms]);d=direction(origin)
 if d!=m['direction']:raise InvalidArtifact('direction changed')
 for t in m['tasks']:
  a=xyz(verify(t['xyz']));expected=q.copy();expected[d['metal_index']]+=STEPS[t['task_id']]*np.array(d['unit_vector'])
  if [x[0] for x in a]!=[x[0] for x in atoms] or not np.allclose([x[1:] for x in a],expected,atol=6e-11,rtol=0):raise InvalidArtifact('physical displacement differs')
  if scout.state(verify(t['xyz']))!=m['electronic_state']:raise InvalidArtifact('state changed')
  if verify(t['input']).read_text()!=scout.input_text(verify(m['basis']).read_text(),verify(m['aux']).read_text()):raise InvalidArtifact('Hamiltonian changed')
  bare={k:v for k,v in t.items() if k!='cache_key'}
  if t['cache_key']!=cache_key(dict(task=bare,protocol=PROTOCOL,implementation=m['implementation'],origin=m['origin_collection'],basis=m['basis'],aux=m['aux'],state=m['electronic_state'],orca=m['orca'])):raise InvalidArtifact('cache mismatch')
 verify(m['origin_collection']);return dry_run(path)
def prepare(a):
 origin=scout.collect(a.origin_manifest)
 r=origin['rows']['Dy_origin']
 if r['status']!='complete' or 'gradient_kcal_mol_per_A' not in r:raise InvalidArtifact('origin energy/analytic gradient unavailable')
 om=read_json(a.origin_manifest);w=Path(a.output).resolve();w.mkdir(parents=True,exist_ok=False)
 write_new(w/'ORIGIN_COLLECTION.json',origin);impl=w/'implementation';impl.mkdir();pins={}
 for p in Path(__file__).resolve().parent.glob('*.py'):shutil.copyfile(p,impl/p.name);pins[p.name]=record(impl/p.name)
 shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
 for n in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[n]=record(impl/n)
 m={k:om[k] for k in ('method','electronic_state','source_xyz','basis','aux','orca')};m.update(protocol_id=PROTOCOL,agreement=record(a.plan),origin_collection=record(w/'ORIGIN_COLLECTION.json'),implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':a.ranks,'concurrent_tasks':a.workers},tasks=[],direction=direction(verify(om['source_xyz'])))
 atoms=xyz(verify(m['source_xyz']));q=np.array([x[1:] for x in atoms]);dv=m['direction']
 for name,h in STEPS.items():
  d=w/name;d.mkdir();pos=q.copy();pos[dv['metal_index']]+=h*np.array(dv['unit_vector'])
  (d/'core.xyz').write_text(str(len(atoms))+'\nmodeled metal displacement, not an observed structure\n'+''.join(f'{atom[0]} {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n' for atom,v in zip(atoms,pos)))
  (d/'endpoint.inp').write_text(scout.input_text(verify(m['basis']).read_text(),verify(m['aux']).read_text()))
  t=dict(task_id=name,metal='Dy',charge=-1,multiplicity=1,physical_multiplicity=6,xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient')
  t['cache_key']=cache_key(dict(task=t,protocol=PROTOCOL,implementation=pins,origin=m['origin_collection'],basis=m['basis'],aux=m['aux'],state=m['electronic_state'],orca=m['orca']));m['tasks'].append(t)
 write_new(w/'manifest.json',m);return validate(w/'manifest.json')
def collect(path):
 from ggr_sensitivity import executed
 validate(path);m,rows=executed(path)
 for t in m['tasks']:
  r=rows[t['task_id']]
  if r['status']=='complete':
   try:r.update(scout.parse(t,m))
   except (ValueError,OSError) as e:r.update(status='invalid',reason=str(e),energy_hartree=None)
 result=dict(protocol_id=PROTOCOL,rows=rows,force_consistency=None,classification=None,affinity=None)
 if all(r['status']=='complete' for r in rows.values()):
  o=read_json(verify(m['origin_collection']))['rows']['Dy_origin'];d=m['direction'];g=float(np.dot(o['gradient_kcal_mol_per_A'][d['metal_index']],d['unit_vector']))
  fd=(rows['plus']['energy_hartree']-rows['minus']['energy_hartree'])*HA_TO_KCAL/.01
  half=(rows['half_plus']['energy_hartree']-rows['half_minus']['energy_hartree'])*HA_TO_KCAL/.005
  result['force_consistency']=dict(analytic_kcal_mol_A=g,central_full=fd,central_half=half,full_error=abs(fd-g),half_error=abs(half-g),refinement_change=abs(fd-half),pass_declared_tolerance=abs(fd-g)<=.1 and abs(half-g)<=.1 and abs(fd-half)<=.05)
 return result
def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('prepare')
 for k in ('origin-manifest','plan','output'):a.add_argument('--'+k,required=True)
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
