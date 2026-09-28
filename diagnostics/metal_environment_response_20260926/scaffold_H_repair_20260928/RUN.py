"""Bounded real-source C-bound-H geometry preparation; no electronic calls."""
from pathlib import Path
import sys,os,time,resource,json
import numpy as np
import openmm as mm
from openmm import unit
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new

def run(source,out):
 p=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1'/source;c=read_json(p/'INPUTS.json');ex=read_json(verify(c['export']));atoms=read_json(verify(ex['atoms']));bonds=read_json(verify(ex['bonds']));ids={v['id']:i for i,v in enumerate(atoms)};adj={i:set() for i in range(len(atoms))}
 for a,b in bonds:i,j=ids[a],ids[b];adj[i].add(j);adj[j].add(i)
 free=np.array([i for i,v in enumerate(atoms) if v['element']=='H' and len(adj[i])==1 and atoms[next(iter(adj[i]))]['element']=='C']);fixed=sorted(set(range(len(atoms)))-set(free));x0=np.array([a['xyz_A'] for a in atoms]);system=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());nb=next(v for v in system.getForces() if isinstance(v,mm.NonbondedForce))
 for i in range(nb.getNumParticles()):q,s,e=nb.getParticleParameters(i);nb.setParticleParameters(i,0.,s,e)
 for i in range(nb.getNumExceptions()):a,b,q,s,e=nb.getExceptionParameters(i);nb.setExceptionParameters(i,a,b,0.,s,e)
 nb.setUseDispersionCorrection(False);integrator=mm.VerletIntegrator(.001);ctx=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'));calls=0;start=time.monotonic()
 def evaluate(y):
  nonlocal calls
  x=x0.copy();x[free]=np.asarray(y).reshape(-1,3);ctx.setPositions(x*.1);st=ctx.getState(getEnergy=True,getForces=True);calls+=1;e=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);g=-st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom);return e,np.asarray(g),x
 def objective(y):e,g,x=evaluate(y);return e,g[free].ravel()
 before,gb,_=evaluate(x0[free]);initial=x0[free].ravel();res=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=list(zip(initial-2.2,initial+2.2)),options={'maxiter':1000,'maxfun':5000,'maxls':20,'gtol':.05,'ftol':1e-12});after,ga,x=evaluate(res.x)
 af=next(v for v in system.getForces() if isinstance(v,mm.HarmonicAngleForce));bf=next(v for v in system.getForces() if isinstance(v,mm.HarmonicBondForce))
 def angle_audit(xx):
  total=0.;maxdev=0.;rows=[]
  for i in range(af.getNumAngles()):
   a,b,d,t,k=af.getAngleParameters(i);a,b,d=int(a),int(b),int(d);u=xx[a]-xx[b];v=xx[d]-xx[b];theta=np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1));dev=theta-t.value_in_unit(unit.radian);energy=.5*k.value_in_unit(unit.kilojoule_per_mole/unit.radian**2)*dev**2/4.184
   if any(atoms[j]['element']=='H' for j in (a,b,d)):total+=energy
   if any(j in freeset for j in (a,b,d)):maxdev=max(maxdev,abs(float(np.degrees(dev))));rows.append({'ids':[atoms[j]['id'] for j in (a,b,d)],'deviation_deg':float(np.degrees(dev)),'energy_kcal_mol':float(energy)})
  return {'H_angle_energy_kcal_mol':total,'max_freeH_angle_deviation_deg':maxdev,'largest_freeH_angles':sorted(rows,key=lambda v:abs(v['deviation_deg']),reverse=True)[:10]}
 freeset=set(free.tolist());ab=angle_audit(x0);aa=angle_audit(x);bondmax=0.
 for i in range(bf.getNumBonds()):
  a,b,r,k=bf.getBondParameters(i);a,b=int(a),int(b)
  if a in freeset or b in freeset:bondmax=max(bondmax,abs(np.linalg.norm(x[a]-x[b])-r.value_in_unit(unit.angstrom)))
 newclashes=[];mind=float('inf')
 for i in free:
  excluded={int(i)}|adj[i]|set().union(*(adj[j] for j in adj[i]));others=np.array(sorted(set(range(len(atoms)))-excluded));d=np.linalg.norm(x[others]-x[i],axis=1);mind=min(mind,float(d.min()))
  for j in others[d<.65]:
   if np.linalg.norm(x0[j]-x0[i])>=.65:newclashes.append([atoms[int(i)]['id'],atoms[int(j)]['id']])
 inversions=[]
 for i,neighbors in adj.items():
  if len(neighbors)!=4 or atoms[i]['element']=='H' or not (neighbors&freeset):continue
  nn=sorted(neighbors);v0=np.linalg.det((x0[nn[:3]]-x0[nn[3]]));v=np.linalg.det((x[nn[:3]]-x[nn[3]]))
  if abs(v0)>1e-6 and v*v0<0:inversions.append(atoms[i]['id'])
 displacement=np.linalg.norm(x-x0,axis=1);gradmax=float(np.max(abs(ga[free])));checks={'optimizer_success':bool(res.success),'freeH_gradient':gradmax<=.1,'fixed_exact':bool(np.array_equal(x[fixed],x0[fixed])),'C_H_lengths':bondmax<=.06,'H_angles_improve':aa['H_angle_energy_kcal_mol']<ab['H_angle_energy_kcal_mol'],'H_angles_bounded':aa['max_freeH_angle_deviation_deg']<=25,'displacement':float(displacement.max())<=2.2,'no_new_clashes':not newclashes,'no_inversions':not inversions};out.mkdir();(out/'system_objective.xml').write_text(mm.XmlSerializer.serialize(system));write_new(out/'coordinates_A.json',{'source_ids':[a['id'] for a in atoms],'xyz_A':x.tolist()});write_new(out/'forces.json',{'before_kcal_mol_A':(-gb).tolist(),'after_kcal_mol_A':(-ga).tolist()})
 result={'source_id':source,'inputs':record(p/'INPUTS.json'),'parent_export':c['export'],'objective_system':record(out/'system_objective.xml'),'coordinates':record(out/'coordinates_A.json'),'forces':record(out/'forces.json'),'free_H_source_ids':[atoms[int(i)]['id'] for i in free],'status':'admitted' if all(checks.values()) else 'failed','checks':checks,'optimizer':{'success':bool(res.success),'message':str(res.message),'iterations':int(res.nit),'evaluations':int(res.nfev)},'actual_energy_force_calls':calls,'wall_seconds':time.monotonic()-start,'before_energy_kcal_mol':before,'after_energy_kcal_mol':after,'before_angles':ab,'after_angles':aa,'max_freeH_gradient_kcal_mol_A':gradmax,'max_displacement_A':float(displacement.max()),'max_CH_length_error_A':bondmax,'minimum_freeH_nonbonded_distance_A':mind,'new_clashes':newclashes,'inversions':inversions,'missing_ion_terms':'explicitly absent from geometry-preparation surrogate','electronic_calls':0};write_new(out/'RESULT.json',result);return result

def main():
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('allocation required')
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);rows=[]
 for source in ['Hans8DQ2','Hans8FNR','Mex8FNS']:
  try:r=run(source,out/source);rows.append({'source_id':source,'status':r['status'],'result':record(out/source/'RESULT.json')});print(source,r['status'],r['checks'],flush=True)
  except Exception as e:
   import traceback;traceback.print_exc();rows.append({'source_id':source,'status':'error','reason':str(e)})
 write_new(out/'COLLECTION.json',{'rows':rows,'implementation':record(__file__),'plan':record(Path(__file__).with_name('PLAN.md')),'job_id':os.environ['SLURM_JOB_ID'],'allocated_CPUs':os.environ.get('SLURM_CPUS_ON_NODE'),'process_peakRSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'electronic_calls':0})
if __name__=='__main__':main()
