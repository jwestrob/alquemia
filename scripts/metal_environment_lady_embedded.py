"""Owned execution/collection bridge for prepared frozen-core embedded response."""
import argparse,json,os,re
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,HA_TO_KCAL
from affordable_workflow import execute
from metal_environment_lady_embedded_prepare import validate,PROTOCOL
from metal_environment_frozen_f import parse
from metal_environment_lanm_reference import torsion_tangent
from metal_environment_reference import read_pcgrad

def collect(path):
 from ggr_sensitivity import executed
 validate(path);m,rows=executed(path);c=read_json(verify(m['inputs']))
 for t in m['tasks']:
  r=rows[t['task_id']]
  if r['status']!='complete':continue
  try:
   text=verify(r['output']).read_text()
   if m.get('initial_guess')=='MORead':
    actual='\n'.join(line for line in text.splitlines() if not re.match(r'\s*\|\s*\d+>',line))
    if 'initial.gbw' not in actual or 'INITIAL GUESS: MODEL POTENTIAL' in actual:raise InvalidArtifact('explicit orbital-read execution evidence missing or fallback used')
    r['orbital_initialization']=dict(seed=t['orbital_seed'],gbw=t['initial_gbw'],target_energy_reused=False)
   if m.get('solver')=='TRAH' and not re.search(r'\((?:TRAH|NR) MAcro\)',text):raise InvalidArtifact('actual TRAH iteration evidence missing')
   n=int(verify(t['pointcharges']).read_text().splitlines()[0])
   counts=re.findall(r'Reading point charge file\s+\.{2,}\s+ok\s+\((\d+) point charges\)',text)
   if not counts or any(int(x)!=n for x in counts) or 'environment.pc' not in text:raise InvalidArtifact('executed permanent field differs')
   r.update(parse(t,{'electronic_state':t['electronic_state']}));tan=torsion_tangent(c,t['configuration'])
   r['torsion_gradient_kcal_mol_per_radian']=float(np.sum(np.array(r['gradient_kcal_mol_per_A'])*tan))
   p=Path(t['output_path']).with_name('endpoint.runtime.pcgrad');g=read_pcgrad(p,n)
   r.update(pointcharge_gradient=record(p),pointcharge_gradient_units='Hartree/bohr',pointcharge_gradient_max_abs=float(abs(g).max()))
  except (ValueError,OSError) as e:r.update(status='invalid',reason=str(e),energy_hartree=None)
 works={}
 for metal in ('La','Dy'):
  a,b=[rows.get(metal+'_'+x,{'status':'unavailable'}) for x in ('A','B')];works[metal]=(b['energy_hartree']-a['energy_hartree'])*HA_TO_KCAL if a['status']==b['status']=='complete' else None
 return dict(protocol_id=m['protocol_id'],manifest=record(path),collector=record(__file__),rows=rows,per_metal_perturbation_work_kcal_mol=works,Dy_minus_La_perturbation_work_kcal_mol=works['Dy']-works['La'] if all(x is not None for x in works.values()) else None,source_id=m['source_id'],energy_scope='finite embedded electronic component; not full hybrid',classification=None,affinity=None,force_consistency_qualified=False)

def main():
 p=argparse.ArgumentParser();p.add_argument('op',choices=['collect','execute','dry-run']);p.add_argument('--manifest',required=True);p.add_argument('--output');p.add_argument('--force-gate');a=p.parse_args()
 if a.op=='collect':r=collect(a.manifest)
 elif a.op=='dry-run':r=validate(a.manifest)
 else:
  validate(a.manifest)
  if not a.force_gate or not (read_json(a.force_gate).get('force_consistency') or {}).get('pass_declared_tolerance'):raise InvalidArtifact('real local force gate must pass')
  m=read_json(a.manifest);layout=m['execution_resources']
  if int(os.environ.get('SLURM_NTASKS','0'))<layout['mpi_ranks']*layout['concurrent_tasks']:raise InvalidArtifact('MPI slots insufficient')
  if 'metal_environment_lady_embedded.py' not in m['implementation'] or record(__file__)!=m['implementation']['metal_environment_lady_embedded.py']:raise InvalidArtifact('execute from new pinned implementation snapshot')
  os.environ['METAL_ENV_WORKERS']=str(layout['concurrent_tasks'])
  write_new(Path(a.manifest).parent/'ADMISSION.json',dict(force_gate=record(a.force_gate),manifest=record(a.manifest),executor=record(__file__)))
  r=execute(a.manifest)
 if a.output:write_new(a.output,r)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
