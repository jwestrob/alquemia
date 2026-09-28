"""Separate Dy(III) 4f-in-core reference scout. Never changes physical spin metadata."""
import argparse,json,os,re,shutil
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,cache_key,xyz,energy,BOHR_TO_A,HA_TO_KCAL
from affordable_workflow import dry_run,execute
from affordable_response import read_engrad
from metal_environment_electronic_state import NUMBERS
from metal_environment_reference import ORCA
PROTOCOL='nikasha_DyIII_frozen4f_pbe0_force_scout_v1'
METHOD='PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad'
def state(path):
 atoms=xyz(path);symbols=[a[0] for a in atoms]
 metals=[x for x in symbols if x in ('La','Dy')]
 if len(metals)!=1 or set(symbols)-{'Dy','La','H','C','N','O'}:raise InvalidArtifact('one La/Dy and closed-shell CHNO ligands required')
 metal=metals[0];core=55 if metal=='Dy' else 46
 total=sum(NUMBERS[s] for s in symbols)+1;explicit=total-core
 if explicit%2:raise InvalidArtifact('effective closed-shell electron parity differs')
 return dict(metal=metal,oxidation_state_hypothesis=3,charge=-1,physical_multiplicity=6 if metal=='Dy' else 1,physical_f_occupation=9 if metal=='Dy' else 0,effective_multiplicity=1,ecp_core_electrons=core,all_electron_count=total,explicit_electrons=explicit,alpha_electrons=explicit//2,beta_electrons=explicit//2,representation='spin_free_frozen_4f9_valence_model' if metal=='Dy' else 'spin_free_4f0_valence_model',spin_orbit_included=False)
def input_text(basis,aux):
 return f'! {METHOD}\n%scf\n Guess PModel\nend\n%basis\n{basis}\n{aux}\nend\n* xyzfile -1 1 core.xyz\n'
def validate(path):
 m=read_json(path)
 if m['protocol_id']!=PROTOCOL or len(m['tasks'])!=1:raise InvalidArtifact('single scout required')
 for p in m['implementation'].values():verify(p)
 t=m['tasks'][0]
 if state(verify(t['xyz']))['metal']!='Dy' or state(verify(t['xyz']))!=m['electronic_state'] or t['xyz']['sha256']!=m['source_xyz']['sha256']:raise InvalidArtifact('state/source mismatch')
 if verify(t['input']).read_text()!=input_text(verify(m['basis']).read_text(),verify(m['aux']).read_text()):raise InvalidArtifact('input/basis differs')
 bare={k:v for k,v in t.items() if k!='cache_key'}
 if t['cache_key']!=cache_key(dict(task=bare,protocol=PROTOCOL,implementation=m['implementation'],orca=m['orca'],state=m['electronic_state'],basis=m['basis'],aux=m['aux'])):raise InvalidArtifact('cache identity mismatch')
 return dry_run(path)
def prepare(a):
 w=Path(a.output).resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir();pins={}
 for p in Path(__file__).resolve().parent.glob('*.py'):shutil.copyfile(p,impl/p.name);pins[p.name]=record(impl/p.name)
 shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
 for n in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[n]=record(impl/n)
 m=dict(protocol_id=PROTOCOL,method=METHOD,electronic_state=state(a.source),source_xyz=record(a.source),basis=record(a.basis),aux=record(a.aux),agreement=record(a.plan),orca=record(ORCA),implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':a.ranks,'concurrent_tasks':1},tasks=[],compute_budget=None,wall_time_limit=None)
 d=w/'Dy_origin';d.mkdir();shutil.copyfile(a.source,d/'core.xyz');(d/'endpoint.inp').write_text(input_text(Path(a.basis).read_text(),Path(a.aux).read_text()))
 t=dict(task_id='Dy_origin',metal='Dy',charge=-1,multiplicity=1,physical_multiplicity=6,xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient')
 t['cache_key']=cache_key(dict(task=t,protocol=PROTOCOL,implementation=pins,orca=m['orca'],state=m['electronic_state'],basis=m['basis'],aux=m['aux']));m['tasks']=[t];write_new(w/'manifest.json',m);return validate(w/'manifest.json')
def parse(t,m):
 text=Path(t['output_path']).read_text();e=energy(t['output_path']);st=m['electronic_state']
 checks=[(r'Program Version\s+(6\.1\.1)\b','6.1.1'),(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)','-1'),(r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)','1'),(r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',str(st['explicit_electrons'])),(r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)','RHF')]
 for pattern,want in checks:
  if re.findall(pattern,text)!=[want]:raise InvalidArtifact('executed state/version mismatch: '+pattern)
 ecps=re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',text)
 if len(ecps)!=1 or ecps[0][0]!=st['metal'] or ecps[0][2]!=str(st['ecp_core_electrons']):raise InvalidArtifact('actual element/core differs')
 for pattern in [r'ORCA SCF GRADIENT CALCULATION',r'DISPERSION GRADIENT',r'ECP gradient\s+\(SHARK\)\s+\.{2,}\s+done',r'CARTESIAN GRADIENT',r'DFTD4']:
  if not re.search(pattern,text,re.I):raise InvalidArtifact('missing analytic component '+pattern)
 if re.search(r'gCP correction\s+[-+0-9.]|^\s*CPCM SOLVATION MODEL',text,re.M):raise InvalidArtifact('undeclared energy term')
 r=read_engrad(t['engrad_path']);atoms=xyz(verify(t['xyz']))
 if r['atomic_numbers'].tolist()!=[NUMBERS[x[0]] for x in atoms] or abs(r['energy_Ha']-e)>1e-8:raise InvalidArtifact('gradient inventory/energy mismatch')
 if not np.allclose(r['coordinates_bohr']*BOHR_TO_A,[x[1:] for x in atoms],atol=1e-6,rtol=0) or not np.isfinite(r['gradient_Ha_per_bohr']).all():raise InvalidArtifact('gradient coordinates/values invalid')
 return dict(energy_hartree=e,electronic_state=st,gradient_kcal_mol_per_A=(r['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A).tolist(),engrad=record(t['engrad_path']),executed_ecp=ecps,quantity='gradient_not_force',force_consistency_qualified=False,physical_state_qualified=False)
def collect(path):
 from ggr_sensitivity import executed
 validate(path);m,rows=executed(path)
 for t in m['tasks']:
  r=rows[t['task_id']]
  if r['status']=='complete':
   try:r.update(parse(t,m))
   except (ValueError,OSError) as e:r.update(status='invalid',reason=str(e),energy_hartree=None)
 return dict(protocol_id=PROTOCOL,manifest=record(path),collector_implementation=record(__file__),rows=rows,classification=None,affinity=None)
def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('prepare')
 for k in ('source','basis','aux','plan','output'):a.add_argument('--'+k,required=True)
 a.add_argument('--ranks',type=int,required=True)
 for op in ('dry-run','execute','collect'):
  a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
 a=p.parse_args()
 if a.op=='prepare':r=prepare(a)
 elif a.op=='dry-run':r=validate(a.manifest)
 elif a.op=='collect':r=collect(a.manifest)
 else:
  validate(a.manifest);os.environ['METAL_ENV_WORKERS']='1';r=execute(a.manifest)
 if a.op!='prepare' and a.output:write_new(a.output,r)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
