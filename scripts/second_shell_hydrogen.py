"""Physical-H-only native OMOL preparation in frozen PQQ donor contexts."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from hydration_square import endpoint
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz

PROTOCOL='PQQ_native_OMOL_context_physical_H_preparation_v1'
SETTINGS={'optimizer':'SLSQP','maxiter':200,'ftol_eV':1e-9,'displacement_radius_A':.35,
          'projected_H_gradient_max_eV_A':.03,'boundary_tolerance_A':1e-5,
          'heavy_and_artificial_caps':'exactly_frozen','start':'single_archived_context',
          'proposal':'final_SLSQP_iterate_when_same_parent_and_domain_valid_regardless_convergence',
          'core_transfer':'physical_H_only_source_and_canonical_PQQ_H'}


def parents(rows,mobile):
    coords=np.array([a[1:] for a in rows]);heavy=[i for i,a in enumerate(rows) if a[0] not in ('H','D')]
    result=[]
    for i in mobile:
        j=min(heavy,key=lambda j:np.linalg.norm(coords[j]-coords[i]))
        distance=float(np.linalg.norm(coords[j]-coords[i]))
        if not .65<distance<1.5:raise InvalidArtifact('unsupported physical H parent distance')
        result.append(j)
    return result


def prepare(source,agreement,output):
    src=read_json(source);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    old=read_json(Path(source).parent/'mace_manifest.json');tasks=[]
    impl=out/'implementation';impl.mkdir();pins={}
    for path in Path(__file__).parent.glob('*.py'):
        q=impl/path.name;shutil.copyfile(path,q);pins[path.name]=record(q)
    for state in src['states']:
        if state['case'] not in ('1H4I','4MAE'):continue
        prep=read_json(verify(state['preparation']))
        for metal in ('Ca','La'):
            task=next(t for t in src['tasks'] if t['case']==state['case'] and t['metal']==metal)
            rows=xyz(verify(task['xyz']))
            mobile=[a['qm_index'] for a in prep['mapping']['source_to_qm']
                    if rows[a['qm_index']][0]=='H' and a['kind'] in ('source','opaque_cofactor')]
            core=state['source']['endpoints'][metal]
            t={'task_id':state['case']+'__'+metal,'case':state['case'],'metal':metal,'charge':task['charge'],
               'xyz':task['xyz'],'source_preparation':state['preparation'],'mobile':mobile,
               'hydrogen_parents':parents(rows,mobile),'core':core,'core_to_context':prep['core_to_context'],
               'state':check_atoms(rows,task['charge'])}
            t['cache_key']=cache_key(t | {'protocol':PROTOCOL,'settings':SETTINGS,'model':old['model'],'implementation':pins})
            tasks.append(t)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'source':record(source),'agreement':record(agreement),
       'tasks':tasks,'model':old['model'],'software':old['software'],'implementation':pins,
       'orca':src['orca'],'execution_policy':src['execution_policy'],'compute_budget':None}
    if len(tasks)!=4:raise InvalidArtifact('exact four PQQ endpoints required')
    write_new(out/'manifest.json',m);return {'status':'prepared','tasks':4,'mobile_H':{t['task_id']:len(t['mobile']) for t in tasks}}


def validate(manifest):
    m=read_json(manifest);verify(m['source']);verify(m['agreement'])
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or len(m['tasks'])!=4:raise InvalidArtifact('scope changed')
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        rows=xyz(verify(t['xyz']))
        if check_atoms(rows,t['charge'])!=t['state'] or parents(rows,t['mobile'])!=t['hydrogen_parents']:
            raise InvalidArtifact('source chemistry changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(payload | {'protocol':PROTOCOL,'settings':SETTINGS,'model':m['model'],'implementation':m['implementation']}):
            raise InvalidArtifact('cache changed')
    return {'status':'pass','tasks':4}


def chemical_checks(initial,coords,mobile,expected_parents):
    rows=[(a[0],*v) for a,v in zip(initial,coords)]
    if parents(rows,mobile)!=expected_parents:raise InvalidArtifact('physical H covalent parent changed')
    fixed=[i for i in range(len(initial)) if i not in mobile];start=np.array([a[1:] for a in initial])
    if not np.array_equal(coords[fixed],start[fixed]):raise InvalidArtifact('frozen heavy/cap coordinates moved')
    extent=np.linalg.norm(coords[mobile]-start[mobile],axis=1)
    if max(extent)>SETTINGS['displacement_radius_A']+1e-7:raise InvalidArtifact('H outside declared displacement domain')
    return extent


def execute(manifest):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from scipy.optimize import minimize
    from mace_omol import input_batch
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():raise InvalidArtifact('allocated GPU required')
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;start=time.monotonic()
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64)
    torch.cuda.reset_peak_memory_stats();calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64')
    model=calc.models[0]
    if type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs']:
        raise InvalidArtifact('native model differs')
    for parameter in model.parameters():parameter.requires_grad_(False)
    versions={k:p._version for k,p in model.named_parameters()};results=[]
    for t in m['tasks']:
        dest=root/'tasks'/t['task_id'];dest.mkdir(parents=True,exist_ok=False);began=time.monotonic()
        initial=xyz(verify(t['xyz']));start_coords=np.array([a[1:] for a in initial]);mobile=t['mobile'];n=len(mobile)
        atoms=Atoms([a[0] for a in initial],positions=start_coords,pbc=False);atoms.info.update(charge=t['charge'],spin=1);atoms.calc=calc
        batch=input_batch(calc,atoms,t['charge'],1);trace=[];last={}
        def objective(x):
            coords=start_coords.copy();coords[mobile]+=x.reshape(n,3);atoms.set_positions(coords)
            energy=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces(),dtype=float)
            gradient=-forces[mobile].reshape(-1)
            trace.append({'energy_eV':energy,'max_H_gradient_eV_A':float(np.max(np.abs(gradient))),
                          'max_H_displacement_A':float(np.max(np.linalg.norm(x.reshape(n,3),axis=1)))})
            last.update(energy=energy,forces=forces,coords=coords,x=x.copy(),gradient=gradient)
            return energy,gradient
        def constraint(x):return SETTINGS['displacement_radius_A']**2-(x.reshape(n,3)**2).sum(axis=1)
        def jacobian(x):
            j=np.zeros((n,3*n))
            for i in range(n):j[i,3*i:3*i+3]=-2*x[3*i:3*i+3]
            return j
        result=minimize(objective,np.zeros(3*n),jac=True,method='SLSQP',
              constraints=[{'type':'ineq','fun':constraint,'jac':jacobian}],
              options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['ftol_eV']})
        if not np.array_equal(result.x,last['x']):objective(result.x)
        extent=chemical_checks(initial,last['coords'],mobile,t['hydrogen_parents'])
        boundary=extent>=SETTINGS['displacement_radius_A']-SETTINGS['boundary_tolerance_A']
        g=last['gradient'].reshape(n,3).copy();u=result.x.reshape(n,3)
        for i in np.flatnonzero(boundary):
            radial=u[i]/np.linalg.norm(u[i]);g[i]-=min(float(g[i]@radial),0.)*radial
        projected=float(np.max(np.abs(g)))
        proposed=[(a[0],*v) for a,v in zip(initial,last['coords'])];xp=dest/'context_proposed.xyz';write_xyz(xp,proposed)
        core=xyz(verify(t['core']['xyz']));core_map={int(a):b for a,b in t['core_to_context'].items()}
        transferred=[]
        for i,j in core_map.items():
            if j in mobile:
                if core[i][0]!='H':raise InvalidArtifact('transferred non-H coordinate')
                core[i]=proposed[j];transferred.append(i)
        cp=dest/'core_proposed.xyz';write_xyz(cp,core)
        ca=Atoms([a[0] for a in core],positions=[a[1:] for a in core],pbc=False);ca.info.update(charge=t['core']['charge'],spin=1);ca.calc=calc
        core_batch=input_batch(calc,ca,t['core']['charge'],1);core_energy=float(ca.get_potential_energy());core_forces=np.asarray(ca.get_forces(),dtype=float)
        fp=dest/'context_forces_eV_A.npy';np.save(fp,last['forces']);cfp=dest/'core_forces_eV_A.npy';np.save(cfp,core_forces)
        r={'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'computed','cache_key':t['cache_key'],
           'manifest':record(manifest),'native_batch':batch,'core_native_batch':core_batch,
           'source_context_energy_eV':trace[0]['energy_eV'],'context_energy_eV':last['energy'],'core_energy_eV':core_energy,
           'context_xyz':record(xp),'core_xyz':record(cp),'context_forces':record(fp),'core_forces':record(cfp),
           'transferred_core_H_indices':transferred,'iterations':int(result.nit),'objective_evaluations':len(trace),
           'scipy_success':bool(result.success),'termination':str(result.message),'trace':trace,
           'maximum_H_displacement_A':float(max(extent)),'boundary_H_count':int(sum(boundary)),
           'projected_gradient_max_eV_A':projected,'constrained_stationary':projected<=SETTINGS['projected_H_gradient_max_eV_A'],
           'same_covalent_parents':True,'heavy_and_caps_exact':True,'wall_seconds':time.monotonic()-began,
           'slurm_job_id':os.environ['SLURM_JOB_ID'],'parameter_versions_unchanged':versions=={k:p._version for k,p in model.named_parameters()}}
        if not r['parameter_versions_unchanged']:raise InvalidArtifact('model parameters changed')
        write_new(dest/'result.json',r);results.append(record(dest/'result.json'))
        print(json.dumps({k:r[k] for k in ('task_id','iterations','objective_evaluations','boundary_H_count','constrained_stationary','wall_seconds')}),flush=True)
    receipt={'manifest':record(manifest),'status':'complete','results':results,'wall_seconds':time.monotonic()-start,
       'slurm_job_id':os.environ['SLURM_JOB_ID'],'allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),
       'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated(),
       'new_native_optimization_tasks':4,'new_native_core_evaluations':4}
    write_new(root/('collection_'+os.environ['SLURM_JOB_ID']+'.json'),receipt);return receipt


def prepare_dft(collection,output):
    col=read_json(collection);m=read_json(verify(col['manifest']));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    if col['status']!='complete' or len(col['results'])!=4:raise InvalidArtifact('complete four actual proposals required')
    impl=out/'implementation';shutil.copytree(Path(m['implementation']['second_shell_hydrogen.py']['path']).parent,impl)
    tasks=[]
    for pin in col['results']:
        result=read_json(verify(pin));t=next(t for t in m['tasks'] if t['task_id']==result['task_id'])
        if not result['same_covalent_parents'] or not result['heavy_and_caps_exact']:raise InvalidArtifact('invalid geometry proposal')
        d=out/'tasks'/t['task_id'];d.mkdir(parents=True);xp=d/'core.xyz';xp.write_bytes(verify(result['core_xyz']).read_bytes())
        ip=d/'endpoint.inp';ip.write_text('! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'+f"* xyzfile {t['core']['charge']} 1 core.xyz\n")
        tasks.append({'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'charge':t['core']['charge'],
           'multiplicity':1,'xyz':record(xp),'input':record(ip),'output_path':str(d/'endpoint.out'),
           'proposal_result':pin,'original_core':t['core']})
    manifest={'protocol_id':'native_r2scan3c_PQQ_context_prepared_physical_H_v1','tasks':tasks,'agreement':m['agreement'],
       'orca':m['orca'],'execution_policy':m['execution_policy'],'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},
       'source_collection':record(collection),'implementation':{p.name:record(p) for p in impl.glob('*.py')},
       'compute_budget':None,'reference':None,'calibrated_decision':None,'baseline_changed':False}
    write_new(out/'manifest.json',manifest)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect_dft(manifest,output):
    m=read_json(manifest);rows=[]
    for t in m['tasks']:
        op=Path(t['output_path']);r=endpoint(record(op),record(str(op)+'.execution.json'),t['xyz'],t['input'])
        rows.append({'case':t['case'],'metal':t['metal'],'endpoint':r,'before_energy_hartree':t['original_core']['energy_hartree'],
                     'proposal':read_json(verify(t['proposal_result']))})
    cases=[]
    for case in ('1H4I','4MAE'):
        pair={r['metal']:r for r in rows if r['case']==case};before=pair['Ca']['before_energy_hartree']-pair['La']['before_energy_hartree']
        after=pair['Ca']['endpoint']['energy_hartree']-pair['La']['endpoint']['energy_hartree']
        cases.append({'case':case,'before_R_hartree':before,'after_R_hartree':after,'delta_R_kcal_mol':(after-before)*HA_TO_KCAL,
            'MACE_after_R_kcal_mol':(pair['Ca']['proposal']['core_energy_eV']-pair['La']['proposal']['core_energy_eV'])*EV_TO_KCAL})
    result={'manifest':record(manifest),'status':'complete','rows':rows,'cases':cases,
        'PQQ_gap_before_kcal_mol':(cases[1]['before_R_hartree']-cases[0]['before_R_hartree'])*HA_TO_KCAL,
        'PQQ_gap_after_kcal_mol':(cases[1]['after_R_hartree']-cases[0]['after_R_hartree'])*HA_TO_KCAL,
        'MACE_gap_after_kcal_mol':cases[1]['MACE_after_R_kcal_mol']-cases[0]['MACE_after_R_kcal_mol'],
        'reference':None,'calibrated_decision':None,'baseline_changed':False}
    write_new(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare')
    for n in ('source','agreement','output'):q.add_argument('--'+n,required=True)
    for name in ('validate','execute'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
    q=s.add_parser('prepare-dft');q.add_argument('--collection',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('collect-dft');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op').replace('-','_');r=globals()[op](**a);print(json.dumps(r,indent=2))
