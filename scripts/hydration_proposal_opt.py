"""Exact rigid-water MACE proposals; native DFT remains the target evaluator."""
from __future__ import annotations
import argparse
import importlib.metadata
import os
from pathlib import Path
import resource
import shutil
import time
import traceback
import numpy as np
from affordable_common import HA_TO_KCAL, BOHR_TO_A, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import accepted_attempt, check_atoms
from hydration_mace import rotational_gradient

STAGE='hydration_orientation_optimization'
SETTINGS={'method':'BFGS','gtol_eV_radian':.001,'maxiter':200,'starts':['source','radial_away']}


def rotate_waters(initial, groups, parameters):
    """Return coordinates and SO(3) left Jacobians for physical water rotations."""
    from scipy.spatial.transform import Rotation
    coords=np.array(initial,copy=True);jacobians=[]
    for w,v in zip(groups,np.asarray(parameters).reshape(-1,3)):
        theta=float(np.linalg.norm(v));x,y,z=v
        skew=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
        if theta<1e-4:
            a=.5-theta**2/24+theta**4/720;b=1/6-theta**2/120+theta**4/5040
        else:a=(1-np.cos(theta))/theta**2;b=(theta-np.sin(theta))/theta**3
        jacobians.append(np.eye(3)+a*skew+b*(skew@skew))
        o=initial[w['oxygen_index']];ids=w['hydrogen_indices']
        coords[ids]=(initial[ids]-o)@Rotation.from_rotvec(v).as_matrix().T+o
    return coords,np.array(jacobians)


def prepare(proposal_collection, manifests, inventory, agreement, output, occupancy=False):
    import mace_omol as omol
    prior=read_json(proposal_collection)
    if not prior.get('checks',{}).get('proposal_utility_pass'):raise InvalidArtifact('proposal utility check did not pass')
    pm=read_json(verify(prior['manifest']))
    _,out,m=omol.common(inventory,verify(pm['software']),agreement,output,STAGE)
    for name in ('hydration_mace.py','hydration_proposal_opt.py'):
        target=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),target);m['implementation'][name]=record(target)
    tasks=[]
    for mp in manifests:
        parent=read_json(mp)
        keys=sorted({(t['case'],t['pattern'],t['metal']) for t in parent['tasks']
                     if not occupancy or t['variable_water_count']>0})
        for case,pattern,metal in keys:
            matches=[t for t in parent['tasks'] if (t['case'],t['pattern'],t['metal'])==(case,pattern,metal)]
            if {t['seed'] for t in matches}!=set(SETTINGS['starts']) or len(matches)!=2:
                raise InvalidArtifact('two predefined water orientation starts required')
            source=next(t for t in matches if t['seed']=='source')
            task={'task_id':case+'__'+(pattern+'__' if occupancy else '')+metal,'case_id':case,'metal':metal,'metal_index':0,
                'kind':'core','variant':'primary','energy_component':omol.COMPONENT,'energy_only':False,
                'charge':source['charge'],'spin_multiplicity':1,'xyz':source['xyz'],
                'state':check_atoms(xyz(verify(source['xyz'])),source['charge']),
                'water_groups':source['groups'],'mobile_indices':source['mobile_indices'],
                'starts':{t['seed']:t['xyz'] for t in matches},'water_orientation_optimization':True}
            if occupancy:task.update(pattern=pattern,variable_water_count=source['variable_water_count'])
            tasks.append(task)
    m.update(tasks=tasks,optimization=SETTINGS,source_manifests=[record(p) for p in manifests],
             proposal_collection=record(proposal_collection),evidence_use='consumed_water_orientation_development')
    if occupancy:m.update(occupancy_enumerated=True,protocol_id='mace_omol_rigid_water_occupancy_proposals_v1')
    return omol.seal(out,m)


def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest)
    if m['stage']!=STAGE or m['optimization']!=SETTINGS or m['model']!=omol.model(verify(m['software'])):
        raise InvalidArtifact('optimization settings/model changed')
    verify(m['agreement']);verify(m['proposal_collection'])
    for p in m['implementation'].values():verify(p)
    for p in m['source_manifests']:verify(p)
    if m.get('occupancy_enumerated'):
        expected={}
        for pin in m['source_manifests']:
            for t in read_json(verify(pin))['tasks']:
                if t['variable_water_count']>0:
                    expected.setdefault((t['case'],t['pattern'],t['metal']),{})[t['seed']]=t
        keys=[(t['case_id'],t['pattern'],t['metal']) for t in m['tasks']]
        if len(keys)!=len(set(keys)) or set(keys)!=set(expected):
            raise InvalidArtifact('missing/duplicate manifested occupancy proposals')
        for t in m['tasks']:
            sources=expected[(t['case_id'],t['pattern'],t['metal'])]
            s=sources['source']
            if (t['starts']!={k:v['xyz'] for k,v in sources.items()} or
                t['water_groups']!=s['groups'] or t['mobile_indices']!=s['mobile_indices'] or
                t['charge']!=s['charge'] or t['xyz']!=s['xyz']):
                raise InvalidArtifact('occupancy preparation changed')
            counterpart=(t['case_id'],t['pattern'],'La' if t['metal']=='Ca' else 'Ca')
            if counterpart not in expected:raise InvalidArtifact('unpaired occupancy state')
    elif {(t['case_id'],t['metal']) for t in m['tasks']}!={(c,z) for c in ('1F6S','6IP9') for z in ('Ca','La')} or len(m['tasks'])!=4:
        raise InvalidArtifact('four paired optimization tasks required')
    for t in m['tasks']:
        initial=xyz(verify(t['xyz']));fixed=[i for i in range(len(initial)) if i not in t['mobile_indices']]
        for pin in t['starts'].values():
            rows=xyz(verify(pin))
            if check_atoms(rows,t['charge'])!=t['state'] or any(rows[i]!=initial[i] for i in fixed):
                raise InvalidArtifact('seed composition/fixed geometry changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('optimization cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def worker(manifest,task_id,output,memory_mode):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from scipy.optimize import minimize
    from mace_omol import input_batch, UNAVAILABLE, COMPONENT
    from mace_hybrid import write_xyz
    m=read_json(manifest);t=next(t for t in m['tasks'] if t['task_id']==task_id);out=Path(output)
    start=time.monotonic();r={'task_id':task_id,'cache_key':t['cache_key'],'manifest':record(manifest),
        'status':'unavailable','energy_component':COMPONENT,'memory_mode':memory_mode,**UNAVAILABLE}
    try:
        if not os.environ.get('SLURM_JOB_ID') or memory_mode!='native' or not torch.cuda.is_available():
            raise InvalidArtifact('native allocated GPU required')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats();props=torch.cuda.get_device_properties(0)
        r.update(execution_device='cuda',device={'name':props.name,'total_memory_bytes':props.total_memory})
        calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64')
        model=calc.models[0]
        if type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs']:
            raise InvalidArtifact('model state embedding changed')
        for p in model.parameters():p.requires_grad_(False)
        versions={n:p._version for n,p in model.named_parameters()}
        r['model_load_seconds']=time.monotonic()-start
        groups=[w for w in t['water_groups'] if w['role']=='variable'];runs=[];forces_by_seed={}
        for seed in SETTINGS['starts']:
            rows=xyz(verify(t['starts'][seed]));initial=np.array([a[1:] for a in rows])
            atoms=Atoms([a[0] for a in rows],positions=initial,pbc=False);atoms.info.update(charge=t['charge'],spin=1);atoms.calc=calc
            r['input_state_check']=input_batch(calc,atoms,t['charge'],1)
            trace=[];last={}
            def objective(parameters):
                coords,jacobians=rotate_waters(initial,groups,parameters);atoms.set_positions(coords)
                energy=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces(),dtype=float)
                torques=rotational_gradient(coords,-forces,groups)
                gradient=np.einsum('wij,wi->wj',jacobians,torques).reshape(-1)
                if not np.isfinite(energy) or not np.isfinite(gradient).all():raise InvalidArtifact('nonfinite native objective')
                trace.append({'energy_eV':energy,'parameters_radian':parameters.tolist(),'gradient_eV_radian':gradient.tolist()})
                last.update(coords=coords,forces=forces,energy=energy,gradient=gradient,parameters=np.array(parameters,copy=True))
                return energy,gradient
            began=time.monotonic()
            result=minimize(objective,np.zeros(3*len(groups)),jac=True,method='BFGS',options={'gtol':SETTINGS['gtol_eV_radian'],'maxiter':SETTINGS['maxiter']})
            if not np.array_equal(result.x,last['parameters']):objective(result.x)
            converged=float(np.max(np.abs(last['gradient'])))<=SETTINGS['gtol_eV_radian']
            coords=last['coords'];fixed=[i for i in range(len(rows)) if i not in t['mobile_indices']]
            if not np.array_equal(coords[fixed],initial[fixed]):raise InvalidArtifact('exact frozen geometry changed')
            for w in groups:
                ids=w['indices'];before=initial[ids];after=coords[ids]
                if np.max(np.abs(np.linalg.norm(before[:,None]-before[None,:],axis=2)-np.linalg.norm(after[:,None]-after[None,:],axis=2)))>1e-10:
                    raise InvalidArtifact('rigid water geometry changed')
            xp=out/(seed+'.xyz');write_xyz(xp,[(a[0],*v) for a,v in zip(rows,coords)])
            fp=out/(seed+'_forces.npy');np.save(fp,last['forces']);forces_by_seed[seed]=record(fp)
            runs.append({'seed':seed,'energy_eV':last['energy'],'xyz':record(xp),'forces':record(fp),
                'converged':converged,'gradient_max_eV_radian':float(np.max(abs(last['gradient']))),
                'scipy_success':bool(result.success),'termination':str(result.message),'iterations':int(result.nit),
                'evaluations':len(trace),'trace':trace,'wall_seconds':time.monotonic()-began})
        best=min(runs,key=lambda run:run['energy_eV']);eligible=all(run['converged'] for run in runs)
        r.update(status='computed',energy_eV=best['energy_eV'],forces=best['forces'],charge_check=None,
                 force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy',
                 parameter_versions_unchanged=versions=={n:p._version for n,p in model.named_parameters()},
                 optimization_runs=runs,optimization_eligible=eligible,selected_seed=best['seed'] if eligible else None,
                 proposed_xyz=best['xyz'] if eligible else None,evaluation_seconds=sum(run['wall_seconds'] for run in runs),
                 energy_evaluation_count=sum(run['evaluations'] for run in runs))
    except Exception as exc:
        r.update(status='failed',reason=str(exc));traceback.print_exc()
    finally:
        r.update(wall_seconds=time.monotonic()-start,peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
            peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_available() else None,
            versions={k:importlib.metadata.version(k) for k in ('mace-torch','torch','ase','scipy')},
            slurm_job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
            allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'))
        write_new(out/'result.json',r)
    return r


def collect(manifest):
    validate(manifest);mp=Path(manifest);m=read_json(mp);rows=[]
    for t in m['tasks']:
        found=[(a,accepted_attempt(a,t,mp)) for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*'))]
        valid=[(a,r) for a,r in found if r is not None]
        if not valid:rows.append({'task_id':t['task_id'],'status':'unavailable'});continue
        a,r=valid[-1]
        rows.append({'task_id':t['task_id'],'status':'computed','optimization_eligible':r['optimization_eligible'],
                     'result':record(a/'result.json'),'receipt':record(a/'receipt.json'),'selected_seed':r['selected_seed'],
                     'proposed_xyz':r['proposed_xyz'],'energy_evaluation_count':r['energy_evaluation_count'],
                     'runs':[{k:v for k,v in run.items() if k!='trace'} for run in r['optimization_runs']]})
    return {'manifest':record(mp),'rows':rows,'status':'complete' if all(r['status']=='computed' for r in rows) else 'incomplete',
            'all_optimization_eligible':all(r.get('optimization_eligible',False) for r in rows),
            'DFT_adjudication':'not_run','baseline_changed':False,'occupancy_probabilities':None}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare');q.add_argument('--manifests',nargs='+',required=True)
    q.add_argument('--occupancy',action='store_true')
    for name in ('proposal-collection','inventory','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='prepare':result=prepare(**a)
    else:result=collect(a['manifest']);write_new(a['output'],result)
    print(json.dumps({k:result[k] for k in ('status','tasks','all_optimization_eligible') if k in result},indent=2))
