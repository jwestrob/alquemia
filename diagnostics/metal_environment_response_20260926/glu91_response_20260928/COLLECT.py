"""Collect the declared real origin/minus/plus common pool, without molecular calls."""
import argparse,math,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,HA_TO_KCAL,InvalidArtifact
from metal_environment_force_assembly import assemble

def collect(config,output):
 c=read_json(config);out=Path(output);out.mkdir(parents=True,exist_ok=False)
 motion=read_json(c['motion_inputs']);t=read_json(verify(motion['tangents']))
 if record(c['origin_inputs'])!=t['origin_inputs']:raise InvalidArtifact('Origin mismatch')
 rows={};missing=[]
 for metal in ('La','Dy'):
  rows[metal]={}
  for state,configuration in [('origin','A'),('minus','A'),('plus','B')]:
   origin=state=='origin';inputs=c['origin_inputs'] if origin else c['motion_inputs']
   cp=Path(c['origin_collections'][metal] if origin else c['motion_collection']);task=f'{metal}_{configuration}'
   cl=Path(c['origin_classical_dir'] if origin else c['motion_classical_dir'])/f'{metal}_{configuration}.json'
   if not cp.is_file() or not cl.is_file():missing.append(f'{metal}/{state}');rows[metal][state]=None;continue
   col=read_json(cp)
   if col['rows'].get(task,{}).get('status')!='complete':missing.append(f'{metal}/{state}');rows[metal][state]=None;continue
   a=assemble(inputs,configuration,metal,cp,task,cl);classical=read_json(cl)
   if classical['inputs']['sha256']!=record(inputs)['sha256']:raise InvalidArtifact('Classical preparation mismatch')
   energies=classical['components_kcal_mol']
   if abs(sum(energies.values())-classical['total_kcal_mol'])>1e-8:raise InvalidArtifact('Component sum mismatch')
   ap=out/f'{metal}_{state}_force.json';write_new(ap,a)
   ids=[v['id'] for v in a['physical_atoms']];ix={s:i for i,s in enumerate(t['physical_source_ids'])}
   if set(ids)!=set(ix):raise InvalidArtifact('Physical source inventory mismatch')
   tangent=t['states']['origin' if origin else configuration];q=np.array(tangent['raw_physical_A_per_radian'])[[ix[s] for s in ids]]
   slopes={k:float((np.asarray(g)*q).sum()) for k,g in [('electronic',a['electronic']['gradient_kcal_mol_A']),('finite',a['total_gradient_kcal_mol_A'])]}
   eh=col['rows'][task]['energy_hartree'];rows[metal][state]=dict(energy_hartree=eh,classical_kcal_mol=classical['total_kcal_mol'],gradient_kcal_mol_radian=slopes,assembly=record(ap),collection=record(cp),classical=record(cl))
 result=dict(implementation=record(__file__),configuration=record(config),rows=rows,missing=missing,pool=None,derivatives=None,affinity=None,full_hybrid_qualified=False)
 if not missing:
  pools={};derivatives={};norm=t['states']['origin']['physical_norm_A_per_radian'];h=math.pi/180
  for scope in ('electronic','finite'):
   works={};derivatives[scope]={}
   for metal in ('La','Dy'):
    r=rows[metal];o=r['origin'];works[metal]={s:(v['energy_hartree']-o['energy_hartree'])*HA_TO_KCAL+(v['classical_kcal_mol']-o['classical_kcal_mol'] if scope=='finite' else 0.) for s,v in r.items()}
    fd=(works[metal]['plus']-works[metal]['minus'])/(2*h);analytic=o['gradient_kcal_mol_radian'][scope];error=abs(fd-analytic)/norm
    derivatives[scope][metal]=dict(central_kcal_mol_radian=fd,analytic_kcal_mol_radian=analytic,normalized_absolute_error_kcal_mol_A=error,diagnostic_tolerance_kcal_mol_A=.1,within_diagnostic_tolerance=bool(error<=.1),refinement_tested=False)
   selected={m:min(works[m],key=works[m].get) for m in works};delta={m:works[m][selected[m]] for m in works}
   pools[scope]=dict(work_kcal_mol=works,mathematical_selected=selected,relaxation_work_kcal_mol=delta,Dy_minus_La_pool_change_kcal_mol=delta['Dy']-delta['La'],interpretation='Finite candidate selection, not optimized minima or populations')
  result.update(pool=pools,derivatives=derivatives)
 write_new(out/'RESULT.json',result);return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=collect(a.config,a.output);print({'missing':r['missing'],'pool':r['pool']})
