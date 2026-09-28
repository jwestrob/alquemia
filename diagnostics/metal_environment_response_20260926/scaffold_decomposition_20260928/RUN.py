"""Two real configurations, native OpenMM component accounting only."""
import json,sys,time,os,resource,collections
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
P=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2'

def main():
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('one shared CPU allocation required')
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();ex=read_json(P/'export_v1/EXPORT.json');atoms=read_json(verify(ex['atoms']));parent=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());config=read_json(P/'INPUTS.json');selected=set(read_json(verify(config['boundary_mapping']))['all_selected_QM_source_ids']);qset={i for i,a in enumerate(atoms) if a['id'] in selected};mset=set(range(len(atoms)))-qset;child=mm.System()
 for i in range(len(atoms)):child.addParticle(parent.getParticleMass(i))
 counts={};groupnames={};group=0;bondstrains=[]
 for f in parent.getForces():
  typ=type(f).__name__
  if typ in ('NonbondedForce','CMMotionRemover'):continue
  accessor,ngetter,adder,k=(('getBondParameters','getNumBonds','addBond',2) if typ=='HarmonicBondForce' else ('getAngleParameters','getNumAngles','addAngle',3) if typ=='HarmonicAngleForce' else ('getTorsionParameters','getNumTorsions','addTorsion',4) if typ=='PeriodicTorsionForce' else ('getTorsionParameters','getNumTorsions','addTorsion',8))
  fs={}
  for region in ('QM','cross','MM'):
   new=type(f)()
   if typ=='CMAPTorsionForce':
    for i in range(f.getNumMaps()):new.addMap(*f.getMapParameters(i))
   new.setForceGroup(group);groupnames[group]=typ+'_'+region;group+=1;fs[region]=new;counts[typ+'_'+region]=0
  for i in range(getattr(f,ngetter)()):
   pars=getattr(f,accessor)(i);ids=pars[1:] if typ=='CMAPTorsionForce' else pars[:k];n=sum(int(x) in qset for x in ids);region='QM' if n==len(ids) else 'MM' if n==0 else 'cross';getattr(fs[region],adder)(*pars);counts[typ+'_'+region]+=1
  for fnew in fs.values():child.addForce(fnew)
 nb=next(f for f in parent.getForces() if isinstance(f,mm.NonbondedForce));assert nb.getNonbondedMethod()==mm.NonbondedForce.NoCutoff
 for region,keep in [('QM',qset),('MM',mset)]:
  f=mm.NonbondedForce();f.setNonbondedMethod(f.NoCutoff);f.setUseDispersionCorrection(False);f.setForceGroup(group);groupnames[group]='nonbonded_'+region;group+=1
  for i in range(len(atoms)):
   q,s,e=nb.getParticleParameters(i);f.addParticle(q if i in keep else 0,s,e if i in keep else 0)
  for i in range(nb.getNumExceptions()):
   a,b,q,s,e=nb.getExceptionParameters(i);keep_pair=int(a) in keep and int(b) in keep;f.addException(a,b,q if keep_pair else 0,s,e if keep_pair else 0)
  child.addForce(f)
 expression='138.93545764438198*q1*q2/r+4*sqrt(epsilon1*epsilon2)*((0.5*(sigma1+sigma2)/r)^12-(0.5*(sigma1+sigma2)/r)^6)'
 cross=mm.CustomNonbondedForce(expression)
 for name in ('q','sigma','epsilon'):cross.addPerParticleParameter(name)
 cross.setNonbondedMethod(cross.NoCutoff);cross.setForceGroup(group);groupnames[group]='nonbonded_cross';cross.addInteractionGroup(qset,mset)
 for i in range(len(atoms)):
  q,s,e=nb.getParticleParameters(i);cross.addParticle([q.value_in_unit(unit.elementary_charge),s.value_in_unit(unit.nanometer),e.value_in_unit(unit.kilojoule_per_mole)])
 repl=mm.CustomBondForce('138.93545764438198*q/r+4*epsilon*((sigma/r)^12-(sigma/r)^6)')
 for name in ('q','sigma','epsilon'):repl.addPerBondParameter(name)
 repl.setForceGroup(group)
 for i in range(nb.getNumExceptions()):
  a,b,q,s,e=nb.getExceptionParameters(i);cross.addExclusion(a,b)
  if (int(a) in qset)!=(int(b) in qset):repl.addBond(a,b,[q.value_in_unit(unit.elementary_charge**2),s.value_in_unit(unit.nanometer),e.value_in_unit(unit.kilojoule_per_mole)])
 child.addForce(cross);child.addForce(repl)
 coords={label:np.array([a['xyz_A'] for a in atoms]) for label in ('A','B')};ids={a['id']:i for i,a in enumerate(atoms)}
 for label in coords:
  for a in read_json(verify(config['core_mapping'][label])):
   if a['kind']=='source':coords[label][ids[a['source_id']]]=a['xyz_A']
 results={};platform=mm.Platform.getPlatformByName('Reference');integs=[mm.VerletIntegrator(.001),mm.VerletIntegrator(.001)];contexts=[mm.Context(s,i,platform) for s,i in zip([parent,child],integs)]
 for label,x in coords.items():
  pc,cc=contexts;pc.setPositions(x*.1);cc.setPositions(x*.1);st=pc.getState(getEnergy=True,getForces=True);energy=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);forces=st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom);parts={name:cc.getState(getEnergy=True,groups={idx}).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole) for idx,name in groupnames.items()};total=sum(parts.values());res=total-energy;tol=1e-6+1e-10*abs(energy)
  ranked=[];bf=next(f for f in parent.getForces() if isinstance(f,mm.HarmonicBondForce))
  for i in range(bf.getNumBonds()):
   a,b,r0,k=bf.getBondParameters(i);r=np.linalg.norm(x[int(a)]-x[int(b)])*.1;de=.5*k.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2)*(r-r0.value_in_unit(unit.nanometer))**2/4.184;ranked.append({'source_atoms':[atoms[int(a)]['id'],atoms[int(b)]['id']],'distance_A':r*10,'equilibrium_A':r0.value_in_unit(unit.nanometer)*10,'energy_kcal_mol':de})
  results[label]={'parent_kcal_mol':energy,'components_kcal_mol':parts,'sum_components_kcal_mol':total,'residual_kcal_mol':res,'tolerance_kcal_mol':tol,'passed':abs(res)<tol,'max_force_kcal_mol_A':float(np.linalg.norm(forces,axis=1).max()),'max_force_atom':atoms[int(np.linalg.norm(forces,axis=1).argmax())]['id'],'largest_bond_strains':sorted(ranked,key=lambda d:d['energy_kcal_mol'],reverse=True)[:10]}
 result={'protocol':'Hans_exact_native_scaffold_decomposition_v1','inputs':record(P/'INPUTS.json'),'parent_export':record(P/'export_v1/EXPORT.json'),'implementation':record(__file__),'plan':record(Path(__file__).with_name('PLAN.md')),'platform':'Reference','openmm_version':mm.__version__,'job_id':os.environ['SLURM_JOB_ID'],'allocated_CPUs':os.environ.get('SLURM_CPUS_ON_NODE'),'wall_seconds':time.monotonic()-start,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'counts':counts,'results':results,'B_minus_A_kcal_mol':results['B']['parent_kcal_mol']-results['A']['parent_kcal_mol'],'component_changes_kcal_mol':{k:results['B']['components_kcal_mol'][k]-v for k,v in results['A']['components_kcal_mol'].items()},'ion_terms':'unavailable_no_ion_particles','hybrid_energy':None,'optimization_steps':0,'evaluated_configurations':2,'energy_state_queries':2+2*len(groupnames),'parent_force_queries':2,'status':'passed' if all(v['passed'] for v in results.values()) else 'failed'};write_new(out/'RESULT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
