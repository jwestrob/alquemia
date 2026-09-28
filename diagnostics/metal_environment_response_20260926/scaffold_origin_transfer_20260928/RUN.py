"""Selected fixed-cap Hans classical ledger; finite true-source phi checks."""
from pathlib import Path
import sys,os,json,time,collections,re,resource
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify
sys.path.insert(0,str(Path(__file__).parents[1]/"coupled_scaffold_20260928"))
from result_io import write_new
from affordable_response import cap_jacobians
from metal_environment_component_checks import build_system,GROUPS
from metal_environment_mechanics import lj_from_amber,graph_distances
P=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Hans8DQ2'
IONBASE=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm')
def ion(z):
 p=IONBASE/('frcmod.ions1lm_126_tip3p' if z=='Na' else 'frcmod.ions234lm_126_tip3p');key={'La':'La3+','Dy':'Dy3+','Na':'Na+','Nd':'Nd3+'}[z];matches=[l.split() for l in p.read_text().splitlines() if l.split() and l.split()[0]==key and len(l.split())>=3];assert len(matches)==1;r,e=map(float,matches[0][1:3]);return lj_from_amber(r,e),{'file':record(p),'line':' '.join(matches[0]),'Rmin_over2_A':r,'epsilon_kcal_mol':e}

def prepare(source,out):
 P=ROOT/"workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2"/source
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('one shared CPU allocation required')
 out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();c=read_json(P/'INPUTS.json');ex=read_json(verify(c['export']));a=read_json(verify(c['physical_source']));parent=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());native=read_json(verify(ex['atoms']));nb=next(f for f in parent.getForces() if isinstance(f,mm.NonbondedForce));bmap=read_json(verify(c['boundary_mapping']));selected=set(bmap['all_selected_QM_source_ids']);caps=read_json(verify(c['core_mapping']['A']));env=read_json(verify(c['environment_atoms']));envby={v['id']:v for v in env};particles=[];qset=set();pars={};ionpins={}
 for z in sorted({'La','Dy'}|{v['element'] for v in c['spectator_occupancy']}):pars[z],ionpins[z]=ion(z)
 assert [v["id"] for v in a]==[v["id"] for v in native]
 assert parent.getNumParticles()==len(a)
 assert selected <= {v["id"] for v in a}
 for i,v in enumerate(a):
  q,s,e=nb.getParticleParameters(i);isq=v['id'] in selected;qtilde=0. if isq else envby[v['id']]['charge_e'] if v['id'] in envby else 0.;p=dict(v,region='QM' if isq else 'MM',charge_e=qtilde,native_charge_e=q.value_in_unit(unit.elementary_charge),sigma_nm=s.value_in_unit(unit.nanometer),epsilon_kJ_mol=e.value_in_unit(unit.kilojoule_per_mole));particles.append(p)
  assert abs(p["native_charge_e"]-native[i]["charge_e"])<1e-10
  if isq:qset.add(i)
 for v in c['spectator_occupancy']:particles.append(dict(v,region='MM',representation='frozen formal charge plus explicit pure12-6 TIP3P IOD LJ',**pars[v['element']]))
 target=caps[0];particles.append({'id':'target_EF3','element':'La','kind':'metal','region':'QM','charge_e':0.,'xyz_A':target['xyz_A'],'endpoint_lj':{z:pars[z] for z in ['La','Dy']}});qset.add(len(particles)-1);assert len(particles)==len(a)+4 and len(qset)==len(selected)+1;ids={v['id']:i for i,v in enumerate(particles)};omitted_to_cap={ids[v['omitted_source_id']]:v for v in caps if v['kind']=='cap'};bonded=[];decisions=[];counts=collections.Counter()
 for f in parent.getForces():
  typ=type(f).__name__
  if typ in ['NonbondedForce','CMMotionRemover']:continue
  get,num,k=('getBondParameters','getNumBonds',2) if typ=='HarmonicBondForce' else ('getAngleParameters','getNumAngles',3) if typ=='HarmonicAngleForce' else ('getTorsionParameters','getNumTorsions',4) if typ=='PeriodicTorsionForce' else ('getTorsionParameters','getNumTorsions',8)
  for j in range(getattr(f,num)()):
   raw=getattr(f,get)(j);support=set(map(int,raw[1:] if typ=='CMAPTorsionForce' else raw[:k]));qm=support&qset;mmset=support-qset
   if not mmset:keep=False;reason='wholly_QM'
   elif not qm:keep=True;reason='wholly_MM'
   elif typ=='HarmonicBondForce':keep=True;reason='cross_radial_stiffness'
   elif all(i in omitted_to_cap and ids[omitted_to_cap[i]['retained_source_id']] in support for i in mmset):keep=False;reason='cross_represented_by_actual_caps'
   else:keep=True;reason='cross_unrepresented_exterior'
   counts[typ+':'+reason]+=1
   if keep:bonded.append({'force':typ,'index':j})
   if qm and mmset:decisions.append({'force':typ,'index':j,'source_atoms':[particles[i]['id'] for i in sorted(support)],'keep':keep,'reason':reason})
 assert counts['HarmonicBondForce:cross_radial_stiffness']==4;assert len(omitted_to_cap)==4
 adjacency={i:set() for i in range(len(a))}
 for u,v in read_json(verify(ex['bonds'])):i,j=ids[u],ids[v];adjacency[i].add(j);adjacency[j].add(i)
 exceptions=[]
 for k in range(nb.getNumExceptions()):
  i,j,q,s,e=nb.getExceptionParameters(k);i,j=int(i),int(j);distance=graph_distances(adjacency,i)[j];scale=0. if distance<=2 else 1/1.2 if distance==3 else None;assert scale is not None;origq=particles[i]['native_charge_e']*particles[j]['native_charge_e']*scale;assert abs(origq-q.value_in_unit(unit.elementary_charge**2))<1e-9;region='QM' if i in qset and j in qset else 'MM' if i not in qset and j not in qset else 'cross';exceptions.append({'atoms':[i,j],'region':region,'chargeprod_e2':scale*particles[i]['charge_e']*particles[j]['charge_e'],'sigma_nm':s.value_in_unit(unit.nanometer),'epsilon_kJ_mol':e.value_in_unit(unit.kilojoule_per_mole)})
 artifacts={'particles.json':particles,'exceptions.json':exceptions,'bonded_terms.json':bonded};ledger={'native_system':ex['system']};assert abs(sum(v['charge_e'] for v in particles)-sum(v['charge_e'] for v in env))<1e-9
 for name,value in artifacts.items():write_new(out/name,value)
 write_new(out/'cross_term_decisions.json',decisions);write_new(out/'LEDGER.json',{'inputs':record(P/'INPUTS.json'),'parent_export':c['export'],'native_system':ex['system'],'ion_parameters':ionpins,'counts':dict(counts),'artifacts':{name:record(out/name) for name in artifacts},'implementation':record(__file__),'reused_builder':record(ROOT/'scripts/metal_environment_component_checks.py'),'full_hybrid_energy':'unavailable_electronic_component_not_evaluated','boundary_policy':'fixed_cap_radial_bonds_retained_substitutable_angular_terms_omitted_v1'})

 assert len(ids)==len(particles)
 for v in env:
  assert np.array_equal(particles[ids[v['id']]]['xyz_A'],v['xyz_A'])
  assert abs(particles[ids[v['id']]]['charge_e']-v['charge_e'])<1e-12
 for v in caps:
  if v['kind']=='metal':pos=np.array(particles[-1]['xyz_A'])
  elif v['kind']=='source':pos=np.array(particles[ids[v['source_id']]]['xyz_A'])
  else:
   r=np.array(particles[ids[v['retained_source_id']]]['xyz_A']);m=np.array(particles[ids[v['omitted_source_id']]]['xyz_A']);pos=r+v['length_A']*(m-r)/np.linalg.norm(m-r)
  assert np.max(abs(pos-np.array(v['xyz_A'])))<1e-10
 return c,ledger,artifacts,particles

def main():
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();results={};calls=0
 for source in ['Hans8FNR','Mex8FNS']:
  d=out/source;c,ledger,artifacts,particles=prepare(source,d);x=np.array([v['xyz_A'] for v in particles]);rr={}
  site=read_json(verify(c['core_mapping']['A']))[0]['source_site'];actual_target=f"{site['chain']}/{site['resid']}/{site['insertion_code']}/{site['name']}"
  physical_ids=[v['id'] if v['id']!='target_EF3' else actual_target for v in particles]
  for metal in ['La','Dy']:
   system,description=build_system(ledger,artifacts,metal);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));ctx.setPositions(x*.1);ee={};gradients={}
   for group,name in enumerate(GROUPS):
    st=ctx.getState(getEnergy=True,getForces=True,groups={group});calls+=1;ee[name]=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);gradients[name]=-np.asarray(st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))
   result={'source_ids':physical_ids,'coordinates_A':x,'quantity':'gradient','units':'kcal/mol/angstrom','components_kcal_mol':ee,'component_gradients_kcal_mol_A':gradients,'gradient_kcal_mol_A':sum(gradients.values()),'total_kcal_mol':sum(ee.values()),'metal':metal,'target_source_id':actual_target,'QM_MM_Coulomb_included':False,'C4_induction_included':False,'ledger':record(d/'LEDGER.json'),'inputs':record(ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2'/source/'INPUTS.json'),'spectators':c['spectator_occupancy'],'builder_description':description}
   write_new(d/f'{metal}_A.json',result);rr[metal]=result;del ctx;del integ
  checks={name:bool(abs(rr['La']['components_kcal_mol'][name]-rr['Dy']['components_kcal_mol'][name])<1e-9 and np.max(abs(rr['La']['component_gradients_kcal_mol_A'][name]-rr['Dy']['component_gradients_kcal_mol_A'][name]))<1e-9) for name in GROUPS[:3]}
  results[source]={'endpoints':{m:record(d/f'{m}_A.json') for m in rr},'common_component_checks':checks,'Dy_minus_La_classical_kcal_mol':rr['Dy']['total_kcal_mol']-rr['La']['total_kcal_mol'],'components_Dy_minus_La_kcal_mol':{k:rr['Dy']['components_kcal_mol'][k]-rr['La']['components_kcal_mol'][k] for k in GROUPS}}
 result={'results':results,'component_queries':calls,'QM_calls':0,'optimization_steps':0,'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':time.monotonic()-start,'peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'plan':record(Path(__file__).with_name('PLAN.md')),'implementation':record(__file__),'full_hybrid_qualified':False}
 write_new(out/'RESULT.json',result);print((out/'RESULT.json').read_text(),end='')
if __name__=='__main__':main()
