"""Selected fixed-cap Hans classical ledger; finite true-source phi checks."""
from pathlib import Path
import sys,os,json,time,collections,re,resource
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
from affordable_response import cap_jacobians
from metal_environment_component_checks import build_system,GROUPS
from metal_environment_mechanics import lj_from_amber,graph_distances
P=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Hans8DQ2'
IONBASE=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm')
def ion(z):
 p=IONBASE/('frcmod.ions1lm_126_tip3p' if z=='Na' else 'frcmod.ions234lm_126_tip3p');key={'La':'La3+','Dy':'Dy3+','Na':'Na+'}[z];matches=[l.split() for l in p.read_text().splitlines() if l.split() and l.split()[0]==key and len(l.split())>=3];assert len(matches)==1;r,e=map(float,matches[0][1:3]);return lj_from_amber(r,e),{'file':record(p),'line':' '.join(matches[0]),'Rmin_over2_A':r,'epsilon_kcal_mol':e}

def main():
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('one shared CPU allocation required')
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();c=read_json(P/'INPUTS.json');ex=read_json(verify(c['export']));a=read_json(verify(c['physical_source']));parent=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());native=read_json(verify(ex['atoms']));nb=next(f for f in parent.getForces() if isinstance(f,mm.NonbondedForce));bmap=read_json(verify(c['boundary_mapping']));selected=set(bmap['all_selected_QM_source_ids']);caps=read_json(verify(c['core_mapping']['A']));env=read_json(verify(c['environment_atoms']));envby={v['id']:v for v in env};particles=[];qset=set();pars={};ionpins={}
 for z in ['La','Dy','Na']:pars[z],ionpins[z]=ion(z)
 for i,v in enumerate(a):
  q,s,e=nb.getParticleParameters(i);isq=v['id'] in selected;qtilde=0. if isq else envby[v['id']]['charge_e'] if v['id'] in envby else 0.;p=dict(v,region='QM' if isq else 'MM',charge_e=qtilde,native_charge_e=q.value_in_unit(unit.elementary_charge),sigma_nm=s.value_in_unit(unit.nanometer),epsilon_kJ_mol=e.value_in_unit(unit.kilojoule_per_mole));particles.append(p)
  if isq:qset.add(i)
 for v in c['spectator_occupancy']:particles.append(dict(v,region='MM',**pars[v['element']]))
 target=caps[0];particles.append({'id':'target_EF3','element':'La','kind':'metal','region':'QM','charge_e':0.,'xyz_A':target['xyz_A'],'endpoint_lj':{z:pars[z] for z in ['La','Dy']}});qset.add(len(particles)-1);assert len(particles)==1891 and len(qset)==191;ids={v['id']:i for i,v in enumerate(particles)};omitted_to_cap={ids[v['omitted_source_id']]:v for v in caps if v['kind']=='cap'};bonded=[];decisions=[];counts=collections.Counter()
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
 assert counts['HarmonicBondForce:cross_radial_stiffness']==4;assert counts['HarmonicAngleForce:cross_represented_by_actual_caps']==8;assert counts['PeriodicTorsionForce:cross_represented_by_actual_caps']==11
 adjacency={i:set() for i in range(len(a))}
 for u,v in read_json(verify(ex['bonds'])):i,j=ids[u],ids[v];adjacency[i].add(j);adjacency[j].add(i)
 exceptions=[]
 for k in range(nb.getNumExceptions()):
  i,j,q,s,e=nb.getExceptionParameters(k);i,j=int(i),int(j);distance=graph_distances(adjacency,i)[j];scale=0. if distance<=2 else 1/1.2 if distance==3 else None;assert scale is not None;origq=particles[i]['native_charge_e']*particles[j]['native_charge_e']*scale;assert abs(origq-q.value_in_unit(unit.elementary_charge**2))<1e-9;region='QM' if i in qset and j in qset else 'MM' if i not in qset and j not in qset else 'cross';exceptions.append({'atoms':[i,j],'region':region,'chargeprod_e2':scale*particles[i]['charge_e']*particles[j]['charge_e'],'sigma_nm':s.value_in_unit(unit.nanometer),'epsilon_kJ_mol':e.value_in_unit(unit.kilojoule_per_mole)})
 artifacts={'particles.json':particles,'exceptions.json':exceptions,'bonded_terms.json':bonded};ledger={'native_system':ex['system']};assert abs(sum(v['charge_e'] for v in particles)-7)<1e-10
 for name,value in artifacts.items():write_new(out/name,value)
 write_new(out/'cross_term_decisions.json',decisions);write_new(out/'LEDGER.json',{'inputs':record(P/'INPUTS.json'),'parent_export':c['export'],'native_system':ex['system'],'ion_parameters':ionpins,'counts':dict(counts),'artifacts':{name:record(out/name) for name in artifacts},'implementation':record(__file__),'reused_builder':record(ROOT/'scripts/metal_environment_component_checks.py'),'full_hybrid_energy':'unavailable_electronic_component_not_evaluated','boundary_policy':'fixed_cap_radial_bonds_retained_substitutable_angular_terms_omitted_v1'})
 x=np.array([v['xyz_A'] for v in particles]);ni,ci=ids['A/83//N'],ids['A/83//CA'];graph={i:set(v) for i,v in adjacency.items()};graph[ni].remove(ci);graph[ci].remove(ni);moving={ci};todo=[ci]
 while todo:
  i=todo.pop()
  for j in graph[i]-moving:moving.add(j);todo.append(j)
 assert len(moving)==814 and ni not in moving;mi=sorted(moving);axis=x[ci]-x[ni];axis/=np.linalg.norm(axis);pivot=x[ci];tangent=np.zeros_like(x);tangent[mi]=np.cross(axis,x[mi]-pivot);h=np.deg2rad(.05);pool={}
 for label,angle in [('origin',0.),('minus',-h),('plus',h)]:
  xx=x.copy();v=x[mi]-pivot;xx[mi]=pivot+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+np.outer(v@axis,axis)*(1-np.cos(angle));pool[label]=xx
 cap_checks=[];capxyz={label:[] for label in pool};pcindices=[ids[v['id']] for v in env]
 for label,xx in pool.items():
  qm=copy_map=[]
  for v in caps:
   if v['kind']=='metal':pos=xx[-1]
   elif v['kind']=='source':pos=xx[ids[v['source_id']]]
   else:
    i,j=ids[v['retained_source_id']],ids[v['omitted_source_id']];u=xx[j]-xx[i];pos=xx[i]+v['length_A']*u/np.linalg.norm(u);capxyz[label].append(pos.tolist())
   qm.append({'id':v['id'],'xyz_A':pos.tolist(),'kind':v['kind']})
  write_new(out/f'{label}_physical.json',{'source_ids':[v['id'] for v in particles],'xyz_A':xx.tolist()});write_new(out/f'{label}_core_map.json',qm);write_new(out/f'{label}_field.json',[dict(v,xyz_A=xx[ids[v['id']]].tolist()) for v in env])
 for index,v in enumerate(a for a in caps if a['kind']=='cap'):
  i,j=ids[v['retained_source_id']],ids[v['omitted_source_id']];xx,yy=x[i],x[j];length=v['length_A'];ja,jb=cap_jacobians(xx,yy,length);u=(yy-xx)/np.linalg.norm(yy-xx);f=lambda q,m:q+length*(m-q)/np.linalg.norm(m-q);jerr=0.;eps=1e-5
  for which,J in [(0,ja),(1,jb)]:
   for k in range(3):
    d=np.eye(3)[k]*eps;fd=(f(xx+d,yy)-f(xx-d,yy))/(2*eps) if which==0 else (f(xx,yy+d)-f(xx,yy-d))/(2*eps);jerr=max(jerr,float(np.max(abs(fd-J[:,k]))))
  analytic=ja@tangent[i]+jb@tangent[j];fd=(np.array(capxyz['plus'][index])-np.array(capxyz['minus'][index]))/(2*h);cap_checks.append({'id':v['id'],'jacobian_residual':jerr,'radial_null':float(np.linalg.norm(jb@u)),'phi_tangent_residual_A_per_rad':float(np.max(abs(fd-analytic))),'passed':jerr<1e-8 and np.linalg.norm(jb@u)<1e-12 and np.max(abs(fd-analytic))<1e-5})
 radial=[];bf=next(f for f in parent.getForces() if isinstance(f,mm.HarmonicBondForce))
 for row in decisions:
  if row['reason']!='cross_radial_stiffness':continue
  i,j,r0,k=bf.getBondParameters(row['index']);r=np.linalg.norm(x[int(i)]-x[int(j)]);eq=r0.value_in_unit(unit.angstrom);kk=k.value_in_unit(unit.kilocalorie_per_mole/unit.angstrom**2);dr=1e-4;e=lambda t:.5*kk*(t-eq)**2;numeric=(e(r+dr)-2*e(r)+e(r-dr))/dr**2;err=abs(numeric-kk);radial.append({'atoms':row['source_atoms'],'native_k_kcal_mol_A2':kk,'curvature_residual':err,'passed':err<=1e-6+1e-8*abs(kk)})
 bondres=max(abs(np.linalg.norm(xx[i]-xx[j])-np.linalg.norm(x[i]-x[j])) for xx in pool.values() for i,vv in adjacency.items() for j in vv);af=next(f for f in parent.getForces() if isinstance(f,mm.HarmonicAngleForce));angleres=0.
 def angle(xx,i,j,k):u=xx[i]-xx[j];v=xx[k]-xx[j];return np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v),-1,1))
 for index in range(af.getNumAngles()):
  i,j,k,*_=af.getAngleParameters(index);i,j,k=int(i),int(j),int(k)
  for xx in pool.values():angleres=max(angleres,float(abs(angle(xx,i,j,k)-angle(x,i,j,k))))
 results={};calls=0
 for metal in ['La','Dy']:
  system,description=build_system(ledger,artifacts,metal);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));rr={}
  for label,xx in pool.items():
   ctx.setPositions(xx*.1);ee={};der={}
   for group,name in enumerate(GROUPS):
    st=ctx.getState(getEnergy=True,getForces=True,groups={group});calls+=1;ee[name]=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);grad=-np.asarray(st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom));der[name]=float(np.sum(grad*tangent)) if label=='origin' else None
   rr[label]={'components_kcal_mol':ee,'total_kcal_mol':sum(ee.values()),'origin_directional_derivative_kcal_mol_rad':der}
  checks=[]
  for name in GROUPS:
   fd=(rr['plus']['components_kcal_mol'][name]-rr['minus']['components_kcal_mol'][name])/(2*h);an=rr['origin']['origin_directional_derivative_kcal_mol_rad'][name];tol=.02+2e-5*abs(an);checks.append({'component':name,'analytic':an,'central_difference':fd,'residual':fd-an,'tolerance':tol,'passed':abs(fd-an)<=tol})
  results[metal]={'configurations':rr,'checks':checks,'builder_description':description};del ctx;del integ
 result={'status':'passed' if all(v['passed'] for v in cap_checks+radial+[v for r in results.values() for v in r['checks']]) and bondres<1e-10 and angleres<1e-10 else 'failed','ledger':record(out/'LEDGER.json'),'plan':record(Path(__file__).with_name('PLAN.md')),'implementation':record(__file__),'job_id':os.environ['SLURM_JOB_ID'],'results':results,'cap_checks':cap_checks,'radial_checks':radial,'max_bond_length_change_A':bondres,'max_angle_change_rad':angleres,'moving_atoms':len(moving),'moving_real_QM':len(moving&qset),'moving_MM':len(moving-qset),'moving_field_rows':sum(i in moving for i in pcindices),'max_field_displacement_A':float(np.linalg.norm(pool['plus'][pcindices]-x[pcindices],axis=1).max()),'actual_component_energy_force_calls':calls,'wall_seconds':time.monotonic()-start,'peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'QM_calls':0,'optimization_steps':0,'full_hybrid_qualified':False};result=json.loads(json.dumps(result,default=lambda value:value.item() if isinstance(value,np.generic) else str(value)));write_new(out/'RESULT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
