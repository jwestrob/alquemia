"""Finite physical water configurational integrals, explicitly not occupancy."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import shutil
import time
import traceback
import numpy as np
from scipy.constants import R
from scipy.special import ndtri, logsumexp
from scipy.spatial.transform import Rotation
from scipy.stats import qmc, multivariate_normal
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_basin_coordinates import WaterCoordinates, left_jacobian
from hydration_water_motion import CENTERS, checked_gradient
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms, write_xyz

STAGE = 'water_basin_sampling'
PROTOCOL = 'native_r2scan3c_anchored_omol_finite_rigid_water_integral_v1'
SETTINGS = {'temperature_K':298.15, 'sobol_power':10, 'seeds':[190919,290919],
    'proposal_standard_deviation_multiplier':1.5,
    'domains_A_radian':[[.30,.50],[.45,.80],[.60,1.10]],
    'scramble_agreement_kcal_mol':.10, 'minimum_effective_samples':128,
    'DFT_validation_weighted_absolute_error_kcal_mol':.25,
    'DFT_validation_maximum_error_kcal_mol':.50,
    'representative_rules':['weighted_medoid','extent_q50','extent_q90','outer_highest_weight'],
    'measure':'product d3COM_A * sinc(|rotation|/2)^2 d3rotation_rad; common unit=1 A^3 per water',
    'occupancy_probabilities':None, 'absolute_entropy':None}


def coordinate_of(frame, positions):
    """Exact inverse on our SO(3) chart, retaining original hydrogen identities."""
    result=[]
    for w,c,m in zip(frame.groups,frame.centers,frame.masses):
        ids=w['indices'];target=np.asarray(positions)[ids]
        center=np.sum(target*m[:,None],axis=0)/m.sum()
        rotation,_=Rotation.align_vectors(target-center,frame.initial[ids]-c,weights=m)
        result.extend(np.r_[center-c,rotation.as_rotvec()])
    q=np.array(result)
    if not np.allclose(frame.positions(q),positions,atol=1e-9,rtol=0):
        raise InvalidArtifact('water shapes, fixed coordinates or rigid inverse changed')
    return q


def common_frame(ca,la):
    a=xyz(verify(ca['xyz']));b=xyz(verify(la['xyz']))
    if [r[0] for r in a[1:]]!=[r[0] for r in b[1:]] or ca['groups']!=la['groups']:
        raise InvalidArtifact('paired inventory/mapping differs')
    fa=WaterCoordinates(a,ca['groups']);fb=WaterCoordinates(b,la['groups'])
    if not np.array_equal(fa.initial[fa.fixed],fb.initial[fb.fixed]):
        raise InvalidArtifact('paired fixed physical context differs')
    aligned=fb.initial.copy();swaps=[]
    for w,c1,c2,m in zip(fa.groups,fa.centers,fb.centers,fa.masses):
        ids=w['indices'];hydrogen=[i for i in ids if fa.symbols[i]=='H']
        candidate=aligned.copy();candidate[hydrogen]=candidate[hydrogen[::-1]]
        direct=Rotation.align_vectors(aligned[ids]-c2,fa.initial[ids]-c1,weights=m)[0].magnitude()
        reverse=Rotation.align_vectors(candidate[ids]-c2,fa.initial[ids]-c1,weights=m)[0].magnitude()
        if reverse<direct:
            aligned=candidate;swaps.append(hydrogen)
    q=coordinate_of(fa,aligned)
    mid=fa.positions(q*.5)
    frame=WaterCoordinates([(z,*x) for z,x in zip(fa.symbols,mid)],ca['groups'])
    separation=np.linalg.norm(frame.centers[:,None]-frame.centers[None,:],axis=2)
    np.fill_diagonal(separation,np.inf)
    if np.min(separation)<=2*SETTINGS['domains_A_radian'][-1][0]:
        raise InvalidArtifact('distinct-water COM domains overlap; permutation accounting required')
    frame.paired_hydrogen_swaps=swaps
    return frame


def log_haar(q):
    v=np.asarray(q).reshape((-1,np.asarray(q).shape[-1]//6,6))[:,:,3:]
    theta=np.linalg.norm(v,axis=2)
    return np.sum(2*np.log(np.sinc(theta/(2*np.pi))),axis=1)


def domain_mask(q,translation,rotation):
    v=np.asarray(q).reshape((-1,np.asarray(q).shape[-1]//6,6))
    return np.all(np.linalg.norm(v[:,:,:3],axis=2)<=translation,axis=1)&np.all(np.linalg.norm(v[:,:,3:],axis=2)<=rotation,axis=1)


def extent(q,translation=.6,rotation=1.1):
    v=np.asarray(q).reshape((-1,np.asarray(q).shape[-1]//6,6))
    return np.maximum(np.max(np.linalg.norm(v[:,:,:3],axis=2),axis=1)/translation,
                      np.max(np.linalg.norm(v[:,:,3:],axis=2),axis=1)/rotation)


def draws(mean,covariance):
    pieces=[];seeds=[]
    chol=np.linalg.cholesky(covariance)
    for seed in SETTINGS['seeds']:
        u=qmc.Sobol(len(mean),scramble=True,seed=seed).random_base2(SETTINGS['sobol_power'])
        pieces.append(mean+ndtri(u)@chol.T);seeds.extend([seed]*len(u))
    q=np.concatenate(pieces)
    logp=multivariate_normal.logpdf(q,mean=mean,cov=covariance)
    return q,logp,np.array(seeds)


def integral(q,logp,seeds,changes):
    rt=R*SETTINGS['temperature_K']/4184;rows=[]
    for t,r in SETTINGS['domains_A_radian']:
        inside=domain_mask(q,t,r);logw=np.full(len(q),-np.inf)
        logw[inside]=-changes[inside]/rt+log_haar(q[inside])-logp[inside]
        if not np.isfinite(logw[inside]).all():raise InvalidArtifact('missing computed interior energy')
        total=float(logsumexp(logw));logz=total-np.log(len(q))
        free=-rt*logz;ess=float(np.exp(2*total-logsumexp(2*logw)))
        each=[float(-rt*(logsumexp(logw[seeds==seed])-np.log(np.sum(seeds==seed)))) for seed in SETTINGS['seeds']]
        near=extent(q,t,r)>.9
        tail=float(np.exp(logsumexp(logw[near])-total))
        halves=[]
        for n in (2**(SETTINGS['sobol_power']-1),2**SETTINGS['sobol_power']):
            selected=np.concatenate([np.flatnonzero(seeds==seed)[:n] for seed in SETTINGS['seeds']])
            halves.append(float(-rt*(logsumexp(logw[selected])-np.log(len(selected)))))
        rows.append({'translation_limit_A':t,'rotation_limit_radian':r,'conditional_F_config_kcal_mol':float(free),
            'log_integral':float(logz),'draws_in_domain':int(inside.sum()),'denominator_all_draws':len(q),
            'effective_samples':ess,'scramble_F_kcal_mol':each,'scramble_difference_kcal_mol':abs(each[0]-each[1]),
            'half_full_F_kcal_mol':halves,'outer_10percent_weight_fraction':tail,
            'numerical_gate_pass':ess>=SETTINGS['minimum_effective_samples'] and abs(each[0]-each[1])<=SETTINGS['scramble_agreement_kcal_mol']})
    return rows


def representatives(q,logp,changes):
    rt=R*SETTINGS['temperature_K']/4184
    inside=domain_mask(q,*SETTINGS['domains_A_radian'][-1]);idx=np.flatnonzero(inside)
    lw=-changes[idx]/rt+log_haar(q[idx])-logp[idx];weights=np.exp(lw-logsumexp(lw));ex=extent(q[idx])
    scale=np.tile([.6]*3+[1.1]*3,q.shape[1]//6);mean=weights@q[idx]
    orders=[np.argsort(np.sum(((q[idx]-mean)/scale)**2,axis=1))]
    sorted_idx=np.argsort(ex);cum=np.cumsum(weights[sorted_idx])
    for quantile in (.5,.9):
        pivot=ex[sorted_idx[min(int(np.searchsorted(cum,quantile)),len(idx)-1)]]
        orders.append(np.argsort(abs(ex-pivot)))
    outer=ex>=.9
    if not outer.any():raise InvalidArtifact('no outer-domain sample for prescribed validation')
    orders.append(np.flatnonzero(outer)[np.argsort(-weights[outer])])
    used=set();chosen=[]
    for name,order in zip(SETTINGS['representative_rules'],orders):
        pick=next((int(i) for i in order if int(idx[i]) not in used),None)
        if pick is None:raise InvalidArtifact('representative selection exhausted')
        actual=int(idx[pick]);used.add(actual)
        chosen.append({'kind':name,'sample_index':actual,'q':q[actual].tolist(),
            'normalized_importance_weight':float(weights[pick]),'extent':float(ex[pick]),
            'predicted_change_kcal_mol':float(changes[actual])})
    return chosen


def latest_centers(validations):
    result={};last_manifest=None
    for path in validations:
        v=read_json(path);p=read_json(verify(v['preparation']))
        dc=read_json(verify(v['DFT_collection']));mc=read_json(verify(v['MACE_collection']))
        dm=read_json(verify(dc['manifest']));mm=read_json(verify(mc['manifest']));last_manifest=mm
        dt={t['task_id']:t for t in dm['tasks']};dr={r['task_id']:r for r in dc['rows']};mr={r['task_id']:r for r in mc['rows']}
        if any(x['status']!='complete' for x in (v,dc,mc)):raise InvalidArtifact('partial anchor source')
        for row in v['rows']:
            if row['kind'] not in ('minimum','proposal'):continue
            name=row['center_id'];s=p['selection'][name];c=s['center'];t=dt[row['task_id']];d=dr[row['task_id']];a=mr[row['task_id']]
            curve=read_json(verify(s['basin_result']))['minimum']['curvature']
            if not curve['numerical_pass'] or not curve['positive_definite']:raise InvalidArtifact('invalid proposal covariance')
            checked_gradient(d['gradient'],d['result'],xyz(verify(t['xyz'])))
            result[name]={**c,'xyz':t['xyz'],'DFT_result':d['result'],'DFT_gradient':d['gradient'],
                'MACE_forces':a['forces'],'MACE_energy_eV':a['energy_eV'],'MACE_result':a['result'],
                'curvature':curve,'source_validation':record(path),'native_stationary':row['stationary_within_tolerance']}
    if set(result)!=set(CENTERS):raise InvalidArtifact('four consumed paired centers required')
    return result,last_manifest


def prepare(validation,agreement,output):
    import mace_omol as omol
    centers,source=latest_centers(validation)
    _,out,m=omol.common(verify(source['inventory']),verify(source['software']),agreement,output,STAGE)
    for name in ('water_basin_sampling.py','hydration_basin_coordinates.py','hydration_water_motion.py','hydration_square.py'):
        p=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),p);m['implementation'][name]=record(p)
    # Isolated dispatch extension: no shared source/default or frozen prior runner changes.
    dispatch=out/'implementation/mace_omol.py';text=dispatch.read_text()
    for name,args in [('validate','manifest'),('collect','manifest'),('worker','manifest, task_id, output, memory_mode')]:
        header=f'def {name}({args}):\n';insert=f"    if read_json(manifest).get('stage')=='{STAGE}':\n        from water_basin_sampling import {name} as local_operation\n        return local_operation({args})\n"
        if text.count(header)!=1:raise InvalidArtifact('runner dispatch shape changed')
        text=text.replace(header,header+insert)
    dispatch.write_text(text);m['implementation']['mace_omol.py']=record(dispatch)
    for case in ('1F6S__11','6IP9__110'):
        frame=common_frame(centers[case+'__Ca'],centers[case+'__La'])
        for metal in ('Ca','La'):
            name=case+'__'+metal;c=centers[name]
            endpoint_frame=WaterCoordinates(xyz(verify(c['xyz'])),c['groups'])
            midpoint_target=frame.initial.copy()
            if metal=='La':
                for ids in frame.paired_hydrogen_swaps:midpoint_target[ids]=midpoint_target[ids[::-1]]
            midpoint_q=coordinate_of(endpoint_frame,midpoint_target)
            endpoint_midpoint=endpoint_frame.positions(midpoint_q)
            frame_roundoff=float(np.max(abs(endpoint_midpoint-midpoint_target)))
            rows=[(metal if i==0 else z,*x) for i,(z,x) in enumerate(zip(frame.symbols,endpoint_midpoint))]
            fp=out/(name+'_common.xyz');write_xyz(fp,rows);f=WaterCoordinates(rows,c['groups'])
            q0=coordinate_of(f,np.array([a[1:] for a in xyz(verify(c['xyz']))]))
            b=np.eye(f.dimension)
            for i,q in enumerate(q0.reshape(-1,6)):b[6*i+3:6*i+6,6*i+3:6*i+6]=left_jacobian(q[3:])
            k=b.T@np.array(c['curvature']['symmetric_curvature'])@b
            covariance=(R*SETTINGS['temperature_K']/4184)*np.linalg.inv(k)*SETTINGS['proposal_standard_deviation_multiplier']**2
            q,logp,seeds=draws(q0,covariance);sp=out/(name+'_samples.npz')
            np.savez(sp,q=q,logp=logp,seeds=seeds)
            c.update(common_xyz=record(fp),anchor_q=q0.tolist(),proposal_covariance=covariance.tolist(),samples=record(sp),
                common_frame_shape_roundoff_max_A=frame_roundoff)
            c['common_physical_frame_H_permutations']=frame.paired_hydrogen_swaps if metal=='La' else []
    pp=out/'preparation.json';write_new(pp,{'centers':centers,'settings':SETTINGS,'validations':[record(p) for p in validation],
        'agreement':record(agreement),'protocol_id':PROTOCOL})
    tasks=[]
    for name,c in sorted(centers.items()):
        t={'task_id':name,'case_id':c['case'],'metal':c['metal'],'metal_index':0,'kind':'core','variant':'primary',
            'charge':c['charge'],'spin_multiplicity':1,'energy_component':omol.COMPONENT,'energy_only':False,
            'xyz':c['xyz'],'state':check_atoms(xyz(verify(c['xyz'])),c['charge']),'preparation':record(pp)}
        t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']});tasks.append(t)
    m.update(tasks=tasks,preparation=record(pp),settings=SETTINGS,protocol_id=PROTOCOL,evidence_use='consumed_finite_water_domain_development')
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);p=read_json(verify(m['preparation']))
    if m['stage']!=STAGE or m['settings']!=SETTINGS or p['settings']!=SETTINGS:raise InvalidArtifact('sampling model changed')
    for pin in [m['agreement'],*m['implementation'].values(),*p['validations']]:verify(pin)
    if {t['task_id'] for t in m['tasks']}!=set(CENTERS):raise InvalidArtifact('center coverage changed')
    for t in m['tasks']:
        c=p['centers'][t['task_id']]
        for key in ('xyz','common_xyz','samples','MACE_forces','DFT_gradient'):verify(c[key])
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('sampling cache changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def worker(manifest,task_id,output,memory_mode):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch,UNAVAILABLE,COMPONENT
    m=read_json(manifest);t=next(t for t in m['tasks'] if t['task_id']==task_id)
    c=read_json(verify(m['preparation']))['centers'][task_id];out=Path(output);start=time.monotonic()
    result={'task_id':task_id,'cache_key':t['cache_key'],'manifest':record(manifest),'status':'unavailable',
        'memory_mode':memory_mode,'energy_component':COMPONENT,**UNAVAILABLE}
    try:
        if not os.environ.get('SLURM_JOB_ID') or memory_mode!='native' or not torch.cuda.is_available():raise InvalidArtifact('allocated native GPU required')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64);torch.cuda.reset_peak_memory_stats()
        calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64');model=calc.models[0]
        for parameter in model.parameters():parameter.requires_grad_(False)
        versions={n:p._version for n,p in model.named_parameters()}
        frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups']);anchor=np.array([a[1:] for a in xyz(verify(c['xyz']))])
        atoms=Atoms(frame.symbols,positions=anchor,pbc=False);atoms.info.update(charge=c['charge'],spin=1);atoms.calc=calc
        result['input_state_check']=input_batch(calc,atoms,c['charge'],1)
        e0=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces());old=np.load(verify(c['MACE_forces']))
        if abs(e0-c['MACE_energy_eV'])*EV_TO_KCAL>1e-5 or np.max(abs(forces-old))>1e-6:raise InvalidArtifact('archived native MACE anchor not reproducible')
        fp=out/'anchor_forces.npy';np.save(fp,forces)
        correction=checked_gradient(c['DFT_gradient'],c['DFT_result'],xyz(verify(c['xyz'])))/EV_TO_KCAL+forces
        samples=np.load(verify(c['samples']));q=samples['q'];inside=domain_mask(q,*SETTINGS['domains_A_radian'][-1])
        changes=np.full(len(q),np.nan);energies=np.full(len(q),np.nan);eval_start=time.monotonic();calls=1
        for i in np.flatnonzero(inside):
            coordinates=frame.positions(q[i]);atoms.set_positions(coordinates)
            energy=float(atoms.get_potential_energy());energies[i]=energy
            changes[i]=(energy-e0+np.sum(correction*(coordinates-anchor)))*EV_TO_KCAL;calls+=1
            if calls%128==0:print(json.dumps({'task':task_id,'computed':calls,'seconds':time.monotonic()-eval_start}),flush=True)
        if not np.isfinite(changes[inside]).all():raise InvalidArtifact('nonfinite native energy')
        ap=out/'samples.npz';np.savez(ap,q=q,logp=samples['logp'],seeds=samples['seeds'],changes=changes,energies_eV=energies)
        result.update(status='computed',energy_eV=e0,forces=record(fp),charge_check=None,execution_device='cuda',
            force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy',
            parameter_versions_unchanged=versions=={n:p._version for n,p in model.named_parameters()},
            samples=record(ap),cartesian_anchor_correction_eV_A=correction.tolist(),
            integral=integral(q,samples['logp'],samples['seeds'],changes),representatives=representatives(q,samples['logp'],changes),
            candidate_draws=len(q),rejected_outside_largest_domain=int((~inside).sum()),actual_MACE_evaluations=calls,
            evaluation_seconds=time.monotonic()-eval_start,response_model_status='finite_domain_native_validation_pending',
            occupancy_probabilities=None)
    except Exception as exc:
        result.update(status='failed',reason=str(exc));traceback.print_exc()
    finally:
        result.update(wall_seconds=time.monotonic()-start,peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
            slurm_job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
            allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),versions={n:importlib.metadata.version(n) for n in ('torch','mace-torch','ase','scipy')})
        write_new(out/'result.json',result)
    return result


def collect(manifest):
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        candidates=[(a,accepted_attempt(a,t,manifest)) for a in sorted((Path(manifest).parent/'execution'/t['task_id']).glob('attempt_*'))]
        candidates=[(a,r) for a,r in candidates if r is not None]
        if not candidates:rows.append({'task_id':t['task_id'],'status':'unavailable'});continue
        a,r=candidates[-1];verify(r['samples'])
        rows.append({'task_id':t['task_id'],'status':'computed','result':record(a/'result.json'),'receipt':record(a/'receipt.json'),
            'integral':r['integral'],'actual_MACE_evaluations':r['actual_MACE_evaluations']})
    return {'status':'complete' if all(r['status']=='computed' for r in rows) else 'incomplete','manifest':record(manifest),
        'rows':rows,'occupancy_probabilities':None,'absolute_entropy':None,'baseline_changed':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--validation',action='append',required=True)
    q.add_argument('--agreement',required=True);q.add_argument('--output',required=True)
    for name in ('validate','collect'):
        q=sub.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('op')
    if op=='collect':
        output=args.pop('output');r=collect(**args);write_new(output,r);print(json.dumps({'status':r['status']}))
    else:print(json.dumps(globals()[op](**args),indent=2))
