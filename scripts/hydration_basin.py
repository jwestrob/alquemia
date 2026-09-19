"""DFT-anchored coupled water response, using native MACE and existing runners."""
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
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_water_motion import checked_gradient, METHOD, CENTERS
from hydration_square import endpoint
from hydration_basin_coordinates import WaterCoordinates, projected_hessian, curvature_summary, MASSES
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms, write_xyz

STAGE='hydration_coupled_basin'
PROTOCOL='native_r2scan3c_cpcm_cartesian_anchored_mace_rigid_water_basin_v1'
SETTINGS={'hessian_steps':[.001,.0005], 'translation_limit_A':.20,'rotation_limit_radian':.35,
          'optimizer':'SLSQP','maxiter':200,'ftol_eV':1e-10,
          'minimum_gradient_max_kcal_scaled':.05,'boundary_margin':.001,
          'curvature_symmetry_absolute':.05,'curvature_symmetry_relative':.001,
          'curvature_refinement_absolute':.05,'curvature_refinement_relative':.01,
          'masses_amu':MASSES,'rotation_scale_A_per_radian':1.}


def source_centers(collection):
    c=read_json(collection);m=read_json(verify(c['manifest']))
    if c['status']!='complete':raise InvalidArtifact('complete source arrangement table required')
    rows={r['task_id']:r for r in c['rows']};centers={};model=None
    for t in m['tasks']+m['reused']:
        if not t['variable_water_count']:continue
        r=rows[t['task_id']];atoms=xyz(verify(t['xyz']))
        if verify(t['input']).read_text()!=METHOD+f"\n* xyzfile {t['charge']} 1 core.xyz\n":
            raise InvalidArtifact('target Hamiltonian changed')
        actual=endpoint(r['result']['output'],r['result']['receipt'],t['xyz'],t['input'])
        if actual!=r['result']:raise InvalidArtifact('DFT endpoint mismatch')
        checked_gradient(r['gradient'],actual,atoms)
        mr=read_json(verify(t['mace_result']));mm=read_json(verify(mr['manifest']))
        mt=next(v for v in mm['tasks'] if v['task_id']==mr['task_id'])
        if accepted_attempt(Path(t['mace_result']['path']).parent,mt,verify(mr['manifest']))!=mr:
            raise InvalidArtifact('unqualified center MACE artifact')
        run=next(v for v in mr['optimization_runs'] if v['seed']==t['selected_seed'])
        if xyz(verify(run['xyz']))!=atoms:raise InvalidArtifact('different MACE/DFT center')
        if model is not None and model!=mm['model']:raise InvalidArtifact('incompatible MACE centers')
        model=mm['model'];forces=np.load(verify(run['forces']),allow_pickle=False)
        if forces.shape!=(len(atoms),3) or not np.isfinite(forces).all():raise InvalidArtifact('invalid actual forces')
        centers[t['task_id']]={'xyz':t['xyz'],'groups':t['groups'],'charge':t['charge'],'metal':t['metal'],
            'case':t['case'],'pattern':t['pattern'],'water_count':t['variable_water_count'],
            'DFT_result':actual,'DFT_gradient':r['gradient'],'MACE_result':t['mace_result'],
            'MACE_forces':run['forces'],'MACE_energy_eV':run['energy_eV'],'MACE_manifest':mr['manifest']}
    if len(centers)!=20:raise InvalidArtifact('expected all20 nonempty metal states')
    return m,centers,mm['software'],model


def prepare(collection,inventory,agreement,output):
    import mace_omol as omol
    source,centers,software,model=source_centers(collection)
    _,out,m=omol.common(inventory,verify(software),agreement,output,STAGE)
    for name in ('hydration_basin.py','hydration_basin_coordinates.py','hydration_water_motion.py','hydration_square.py'):
        p=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),p);m['implementation'][name]=record(p)
    if m['model']!=model:raise InvalidArtifact('changed model')
    pp=out/'preparation.json'
    write_new(pp,{'centers':centers,'settings':SETTINGS,'source_collection':record(collection),
                  'agreement':record(agreement),'protocol_id':PROTOCOL})
    tasks=[]
    for name,c in sorted(centers.items()):
        tasks.append({'task_id':name,'case_id':c['case'],'metal':c['metal'],'metal_index':0,
            'kind':'core','variant':'primary','charge':c['charge'],'spin_multiplicity':1,
            'energy_component':omol.COMPONENT,'energy_only':False,'xyz':c['xyz'],
            'state':check_atoms(xyz(verify(c['xyz'])),c['charge']),'preparation':record(pp)})
    m.update(tasks=tasks,preparation=record(pp),settings=SETTINGS,protocol_id=PROTOCOL,
             evidence_use='consumed_water_occupancy_response_development')
    return omol.seal(out,m)


def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest);p=read_json(verify(m['preparation']))
    if m['stage']!=STAGE or m['settings']!=SETTINGS or p['settings']!=SETTINGS or m['model']!=omol.model(verify(m['software'])):
        raise InvalidArtifact('changed coupled water model')
    for pin in [m['agreement'],p['source_collection'],*m['implementation'].values()]:verify(pin)
    generation=p.get('generation',0)
    if generation==0:
        coverage=len(m['tasks'])==20
    else:
        coverage=generation in (1,2) and 0<len(m['tasks'])<=4 and set(p['centers'])<=set(CENTERS)
        verify(p['previous_validation'])
    if not coverage or len(m['tasks'])!=len(p['centers']) or {t['task_id'] for t in m['tasks']}!=set(p['centers']):
        raise InvalidArtifact('state coverage changed')
    for t in m['tasks']:
        c=p['centers'][t['task_id']];atoms=xyz(verify(t['xyz']))
        if t['xyz']!=c['xyz'] or t['charge']!=c['charge'] or t['preparation']!=m['preparation'] or t['spin_multiplicity']!=1:
            raise InvalidArtifact('changed physical state')
        WaterCoordinates(atoms,c['groups'])
        checked_gradient(c['DFT_gradient'],c['DFT_result'],atoms);verify(c['MACE_forces'])
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('changed implementation/preparation cache')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def worker(manifest,task_id,output,memory_mode):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from scipy.optimize import minimize
    from mace_omol import input_batch,UNAVAILABLE,COMPONENT
    m=read_json(manifest);t=next(v for v in m['tasks'] if v['task_id']==task_id)
    c=read_json(verify(m['preparation']))['centers'][task_id];out=Path(output);start=time.monotonic()
    r={'task_id':task_id,'cache_key':t['cache_key'],'manifest':record(manifest),'status':'unavailable',
       'memory_mode':memory_mode,'energy_component':COMPONENT,**UNAVAILABLE}
    try:
        if not os.environ.get('SLURM_JOB_ID') or memory_mode!='native' or not torch.cuda.is_available():
            raise InvalidArtifact('native allocated GPU required')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats();props=torch.cuda.get_device_properties(0)
        r.update(execution_device='cuda',device={'name':props.name,'total_memory_bytes':props.total_memory})
        calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64')
        model=calc.models[0]
        if type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs']:
            raise InvalidArtifact('native model embedding changed')
        for p in model.parameters():p.requires_grad_(False)
        versions={n:p._version for n,p in model.named_parameters()};r['model_load_seconds']=time.monotonic()-start
        rows=xyz(verify(c['xyz']));frame=WaterCoordinates(rows,c['groups'])
        atoms=Atoms(frame.symbols,positions=frame.initial,pbc=False);atoms.info.update(charge=t['charge'],spin=1);atoms.calc=calc
        r['input_state_check']=input_batch(calc,atoms,t['charge'],1)
        gdft=checked_gradient(c['DFT_gradient'],c['DFT_result'],rows)/EV_TO_KCAL
        fold=np.load(verify(c['MACE_forces']),allow_pickle=False)
        correction=gdft+fold;calls=0;trace=[]
        def cartesian(coords):
            nonlocal calls
            atoms.set_positions(coords);energy=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces(),dtype=float)
            calls+=1
            anchored_change=energy-c['MACE_energy_eV']+float(np.sum(correction*(coords-frame.initial)))
            gradient=-forces+correction
            if not np.isfinite(anchored_change) or not np.isfinite(gradient).all():raise InvalidArtifact('nonfinite potential')
            return anchored_change,gradient,energy,forces
        def objective(q,local=frame):
            e,g,_,_=cartesian(local.positions(q));return e,local.gradient(q,g)
        zero=np.zeros(frame.dimension);e0,g0,raw_energy,forces=cartesian(frame.initial)
        if abs(e0)*EV_TO_KCAL>1e-5 or np.max(abs(forces-fold))>1e-6:raise InvalidArtifact('center MACE reuse not reproducible')
        fp=out/'center_forces.npy';np.save(fp,forces)
        def curvature(local):
            def grad(q):return objective(q,local)[1]*EV_TO_KCAL
            coarse,fine=[projected_hessian(grad,local.dimension,step) for step in SETTINGS['hessian_steps']]
            s=curvature_summary(coarse,fine,local.mass_matrix())
            s.update(coarse_raw=coarse.tolist(),fine_raw=fine.tolist(),gradient=grad(np.zeros(local.dimension)).tolist())
            return s
        initial_curve=curvature(frame)
        def objective_trace(q):
            e,g=objective(q);trace.append({'q':q.tolist(),'change_kcal_mol':e*EV_TO_KCAL,
                                          'gradient_max_kcal_scaled':float(np.max(abs(g)))*EV_TO_KCAL})
            return e,g
        def limits(q):
            parts=q.reshape(-1,6)
            return np.r_[SETTINGS['translation_limit_A']**2-np.sum(parts[:,:3]**2,axis=1),
                         SETTINGS['rotation_limit_radian']**2-np.sum(parts[:,3:]**2,axis=1)]
        def limit_jac(q):
            parts=q.reshape(-1,6);n=len(parts);j=np.zeros((2*n,len(q)))
            for i in range(n):j[i,6*i:6*i+3]=-2*parts[i,:3];j[n+i,6*i+3:6*i+6]=-2*parts[i,3:]
            return j
        result=minimize(objective_trace,zero,jac=True,method=SETTINGS['optimizer'],
                        bounds=[(-limit,limit) for _ in frame.groups for limit in
                                ([SETTINGS['translation_limit_A']]*3+[SETTINGS['rotation_limit_radian']]*3)],
                        constraints=[{'type':'ineq','fun':limits,'jac':limit_jac}],
                        options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['ftol_eV']})
        q=result.x;coords=frame.positions(q);e,g=objective(q);parts=q.reshape(-1,6)
        translations=np.linalg.norm(parts[:,:3],axis=1);rotations=np.linalg.norm(parts[:,3:],axis=1)
        interior=bool(np.max(translations)<SETTINGS['translation_limit_A']-SETTINGS['boundary_margin'] and
                      np.max(rotations)<SETTINGS['rotation_limit_radian']-SETTINGS['boundary_margin'])
        xp=out/'proposed.xyz';write_xyz(xp,frame.rows(q))
        minimum_frame=WaterCoordinates(frame.rows(q),c['groups']);minimum_curve=curvature(minimum_frame)
        gradmax=float(np.max(abs(g)))*EV_TO_KCAL
        min_valid=bool(result.success and interior and minimum_curve['numerical_pass'] and
                       minimum_curve['positive_definite'] and gradmax<SETTINGS['minimum_gradient_max_kcal_scaled'])
        if not np.array_equal(coords[frame.fixed],frame.initial[frame.fixed]):raise InvalidArtifact('fixed atom moved')
        r.update(status='computed',energy_eV=raw_energy,forces=record(fp),charge_check=None,
            force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy',
            parameter_versions_unchanged=versions=={n:p._version for n,p in model.named_parameters()},
            initial_curvature=initial_curve,cartesian_anchor_correction_eV_A=correction.tolist(),
            minimum={'q':q.tolist(),'xyz':record(xp),'predicted_change_kcal_mol':e*EV_TO_KCAL,
                     'gradient_max_kcal_scaled':gradmax,'translations_A':translations.tolist(),'rotations_radian':rotations.tolist(),
                     'interior':interior,'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
                     'iterations':int(result.nit),'curvature':minimum_curve,'cheap_minimum_eligible':min_valid,'trace':trace},
            energy_evaluation_count=calls,evaluation_seconds=time.monotonic()-start-r['model_load_seconds'],
            response_model_status='native_coupled_validation_pending',occupancy_probabilities=None)
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
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        found=[(a,accepted_attempt(a,t,manifest)) for a in sorted((Path(manifest).parent/'execution'/t['task_id']).glob('attempt_*'))]
        found=[(a,r) for a,r in found if r is not None]
        if not found:rows.append({'task_id':t['task_id'],'status':'unavailable'});continue
        a,r=found[-1];v=r['minimum']
        rows.append({'task_id':t['task_id'],'status':'computed','result':record(a/'result.json'),
            'receipt':record(a/'receipt.json'),'initial_numerical_pass':r['initial_curvature']['numerical_pass'],
            'initial_negative_modes':r['initial_curvature']['negative_mode_count'],
            'minimum_numerical_pass':v['curvature']['numerical_pass'],'minimum_negative_modes':v['curvature']['negative_mode_count'],
            'minimum_interior':v['interior'],'minimum_eligible':v['cheap_minimum_eligible'],
            'minimum_gradient_max_kcal_scaled':v['gradient_max_kcal_scaled'],
            'predicted_relaxation_kcal_mol':v['predicted_change_kcal_mol'],'evaluations':r['energy_evaluation_count']})
    return {'manifest':record(manifest),'rows':rows,'status':'complete' if all(r['status']=='computed' for r in rows) else 'incomplete',
            'occupancy_probabilities':None,'response_model_status':'native_coupled_validation_pending','baseline_changed':False}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare')
    for name in ('collection','inventory','agreement','output'):q.add_argument('--'+name,required=True)
    q=sub.add_parser('collect')
    for name in ('manifest','output'):q.add_argument('--'+name,required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='collect':
        r=collect(a['manifest']);write_new(a['output'],r);print(json.dumps({'status':r['status'],'rows':len(r['rows'])}))
    else:print(json.dumps(prepare(**a),indent=2))
