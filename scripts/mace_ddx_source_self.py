"""Resolved PCM versus native GK on fixed, real full-protein source states."""
from __future__ import annotations
import argparse
import fcntl
import gc
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new

PROTOCOL='fixed_source_full_protein_ddPCM_GK_component_comparison_v1'
DIRECT_PROTOCOL=PROTOCOL+'_direct_source_v2'
CASES=['GGR_2FW0','GGR_2FVY']
GRIDS={'coarse':{'lmax':6,'n_lebedev':194},'primary':{'lmax':9,'n_lebedev':302},'refined':{'lmax':12,'n_lebedev':590}}
MODEL=dict(model='pcm',solvent_epsilon=78.3,solvent_kappa=0.,eta=.1,shift=0.,incore=False,
           maxiter=300,jacobi_n_diis=20,enable_fmm=True,fmm_multipole_lmax=12,fmm_local_lmax=12,
           n_proc=64,enable_force=False)
TOL=dict(solver=1e-10,potential_au=1e-8,psi=1e-12,refinement_kcal=.1,rigid_kcal=.05,
         repeat_kcal=1e-8,identity_kcal=1e-8,reciprocity_kcal=.05,energy_kcal=1e-8)


def prepare(config,output):
    from mace_density_gk_hybrid import accepted
    c=read_json(config);verify(c['plan']);r=read_json(verify(c['native_collection']));nm=read_json(verify(r['manifest']))
    g=read_json(verify(c['gk_accounting']));software=read_json(verify(c['software_receipt']))
    if c['protocol'] not in (PROTOCOL,DIRECT_PROTOCOL) or not(r['complete'] and r['numerical_pass'] and g['checks_pass']) or software['status']!='built_and_imported':raise InvalidArtifact('qualified native sources and solver build required')
    if c['protocol']==DIRECT_PROTOCOL:verify(c['recovery_plan'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for n in (Path(__file__).name,'affordable_common.py'):shutil.copyfile(Path(__file__).with_name(n),impl/n)
    m=dict(protocol=c['protocol'],config=record(config),plan=c['plan'],software_receipt=c['software_receipt'],python=c['python'],module=c['module'],
           implementation={p.name:record(p) for p in impl.glob('*.py')},model=MODEL,tolerances=TOL,grids=GRIDS,
           groups=[],new_DFT_calls=0,new_MACE_calls=0,requested_forward_solves=20,requested_model_setups=8,new_score=None,baseline_changed=False)
    if c['protocol']==DIRECT_PROTOCOL:m['source_potential_policy']='direct_Coulomb_fsum_v1'
    parent=Path(r['manifest']['path']).parent
    for case in CASES:
        for variant in ('primary','rigid'):
            data=None
            for metal in ('Ca','La'):
                t=next(t for t in nm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']==variant)
                item=accepted(t,parent/'tasks'/t['task_id'])
                if item is None:raise InvalidArtifact('actual native state unavailable')
                s=item[0];b=read_json(verify(t['boundary']));indices=np.array(t['frozen_indices'])-1
                coords=[a['xyz_A'] for a in s['parameters']['atoms']];radii=[a['radius_A'] for a in s['parameters']['atoms']]
                moments=np.array([a['global_'] for a in s['moments']])[indices]
                if np.any(moments[:,1:]!=0):raise InvalidArtifact('source is not monopolar')
                charge=np.zeros(len(coords));charge[indices]=moments[:,0]
                if abs(math.fsum(charge)-(-1 if metal=='Ca' else 0))>5e-5:raise InvalidArtifact('source charge inventory differs')
                now=dict(case_id=case,variant=variant,physical_ids=b['physical_ids'],coordinates_A=coords,radii_A=radii,source_indices=indices.tolist())
                if data is None:data={**now,'endpoints':{}}
                elif any(data[k]!=v for k,v in now.items()):raise InvalidArtifact('paired physical system differs')
                data['endpoints'][metal]=dict(charge_e=charge.tolist(),native_static_output=record(item[1]),native_task=t,
                      source_self_GK_primary_kcal=g['cases'][case]['endpoints'][metal]['source_self_kcal'])
            p=root/'sources'/f'{case}_{variant}.json';p.parent.mkdir(exist_ok=True);write_new(p,data)
            for grid in (GRIDS if variant=='primary' else ['primary']):
                states=['Ca','La']
                if variant==grid=='primary':
                    states.append('zero')
                    if case=='GGR_2FVY':states+=['Ca_repeat','La_repeat']
                group=dict(group_id=case+'_'+variant+'_'+grid,case_id=case,variant=variant,grid=grid,source=record(p),states=states)
                group['cache_key']=cache_key(dict(group=group,config=m['config'],model=MODEL,tolerances=TOL,grids=GRIDS,module=m['module'],implementation=m['implementation']))
                m['groups'].append(group)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol'] not in (PROTOCOL,DIRECT_PROTOCOL) or m['model']!=MODEL or m['tolerances']!=TOL or m['grids']!=GRIDS or len(m['groups'])!=8 or sum(len(g['states']) for g in m['groups'])!=20 or m['new_DFT_calls'] or m['new_MACE_calls']:raise InvalidArtifact('declared PCM scope differs')
    for pin in [m['config'],m['plan'],m['software_receipt'],m['python'],m['module'],*m['implementation'].values()]:verify(pin)
    c=read_json(verify(m['config']))
    if c['protocol']!=m['protocol']:raise InvalidArtifact('source evaluation protocol differs')
    for name in ('environment_configuration','dependency_inventory'):verify(c[name])
    if m['protocol']==DIRECT_PROTOCOL:
        verify(c['recovery_plan'])
        if m.get('source_potential_policy')!='direct_Coulomb_fsum_v1':raise InvalidArtifact('source potential policy differs')
    elif 'source_potential_policy' in m:raise InvalidArtifact('legacy source policy changed')
    expected={case+'_'+variant+'_'+grid for case in CASES for variant in ('primary','rigid') for grid in (GRIDS if variant=='primary' else ['primary'])}
    if {g['group_id'] for g in m['groups']}!=expected:raise InvalidArtifact('case/geometry/grid inventory differs')
    for group in m['groups']:
        payload={k:v for k,v in group.items() if k!='cache_key'}
        if group['cache_key']!=cache_key(dict(group=payload,config=m['config'],model=MODEL,tolerances=TOL,grids=GRIDS,module=m['module'],implementation=m['implementation'])):raise InvalidArtifact('PCM cache mismatch')
        data=read_json(verify(group['source']));n=len(data['physical_ids'])
        wanted=['Ca','La']
        if group['variant']==group['grid']=='primary':
            wanted.append('zero')
            if group['case_id']=='GGR_2FVY':wanted+=['Ca_repeat','La_repeat']
        if group['states']!=wanted or data['case_id']!=group['case_id'] or data['variant']!=group['variant']:raise InvalidArtifact('declared controls or source labels differ')
        if len(set(data['physical_ids']))!=n or len(data['coordinates_A'])!=n or len(data['radii_A'])!=n:raise InvalidArtifact('physical inventory malformed')
        if np.min(data['radii_A'])<=0 or not np.isfinite(data['coordinates_A']).all():raise InvalidArtifact('invalid cavity')
        for metal,e in data['endpoints'].items():
            s=read_json(verify(e['native_static_output']));t=e['native_task'];b=read_json(verify(t['boundary']));ix=np.array(data['source_indices'])
            if t['case_id']!=group['case_id'] or t['variant']!=group['variant'] or t['state']!=metal or b['physical_ids']!=data['physical_ids'] or [a['xyz_A'] for a in s['parameters']['atoms']]!=data['coordinates_A'] or [a['radius_A'] for a in s['parameters']['atoms']]!=data['radii_A'] or t['frozen_indices']!=list(ix+1):raise InvalidArtifact('source cavity changed')
            q=np.array(e['charge_e']);mask=np.ones(n,dtype=bool);mask[ix]=False
            if len(q)!=n or np.any(q[mask]!=0) or not np.array_equal(q[ix],np.array([a['global_'][0] for a in s['moments']])[ix]):raise InvalidArtifact('source/environment charges changed')
            if abs(math.fsum(q)-(-1 if metal=='Ca' else 0))>5e-5:raise InvalidArtifact('charge closure failed')
    return m


def potential_error(cavity,centres,q,phi):
    ix=np.flatnonzero(q);error=0.
    for start in range(0,cavity.shape[1],256):
        points=cavity[:,start:start+256].T
        if len(ix):
            distance=np.linalg.norm(points[:,None,:]-centres[:,ix].T[None,:,:],axis=2)
            if np.min(distance)<1e-8:raise InvalidArtifact('cavity point overlaps a source nucleus')
            expected=np.sum(q[ix][None,:]/distance,axis=1)
        else:expected=np.zeros(len(points))
        error=max(error,float(np.max(abs(expected-phi[start:start+len(points)]))))
    return error


def direct_potential(cavity,centres,q):
    """Analytical monopole source potential; no continuum solver is reimplemented."""
    ix=np.flatnonzero(q);phi=np.zeros(cavity.shape[1])
    if not len(ix):return phi
    for start in range(0,cavity.shape[1],256):
        points=cavity[:,start:start+256].T
        distance=np.linalg.norm(points[:,None,:]-centres[:,ix].T[None,:,:],axis=2)
        if np.min(distance)<1e-8:raise InvalidArtifact('cavity point overlaps source nucleus')
        values=q[ix][None,:]/distance
        phi[start:start+len(points)]=[math.fsum(row) for row in values]
    return phi


def execute(path,output):
    import pyddx
    m=validate(path);root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64:raise InvalidArtifact('declared solver allocation required')
    if record(pyddx.__file__)!=m['module'] or pyddx.__version__!='0.9.0' or np.__version__!='1.26.4':raise InvalidArtifact('unqualified solver module')
    groups={};start=time.monotonic();cpu=time.process_time()
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for group in m['groups']:
            d=root/'groups'/group['group_id'];d.mkdir(parents=True,exist_ok=True)
            if (d/'attempt_0001').exists():raise InvalidArtifact('retain prior attempt; explicit collection/recovery required')
            attempt=d/'attempt_0001';attempt.mkdir();data=read_json(verify(group['source']));setup=time.monotonic()
            result=dict(group_id=group['group_id'],cache_key=group['cache_key'],status='failed',rows={},source=group['source'],reciprocity_error_kcal=None)
            model=None;state=None;arrays={};psis={};solutions={}
            try:
                centres=np.asfortranarray(np.array(data['coordinates_A']).T/BOHR_TO_A);radii=np.array(data['radii_A'])/BOHR_TO_A
                kwargs={**MODEL,**GRIDS[group['grid']]};model=pyddx.Model(sphere_centres=centres,sphere_radii=radii,**kwargs)
                result.update(setup_wall_seconds=time.monotonic()-setup,n_cavity=model.n_cav,n_spheres=len(radii),n_basis=model.n_basis,
                              actual_parameters={k:v for k,v in model.input_parameters.items() if k not in ('sphere_centres','sphere_radii')})
                cavity=np.array(model.cavity);arrays['cavity_bohr']=cavity
                for label in group['states']:
                    metal=label.split('_')[0];q=np.zeros(len(radii)) if label=='zero' else np.array(data['endpoints'][metal]['charge_e'])
                    native_start=time.monotonic();row=dict(state=label,status='failed',forward_solve_started=False)
                    write_new(attempt/(label+'_admission.json'),dict(state=label,manifest=record(path),group=group['group_id'],time_unix=time.time(),slurm_job_id=os.environ['SLURM_JOB_ID']))
                    try:
                        multipoles=np.asfortranarray(q.reshape(1,-1)/np.sqrt(4*np.pi));phi=np.array(model.multipole_electrostatics(multipoles,derivative_order=0)['phi']);psi=np.array(model.multipole_psi(multipoles))
                        if m['protocol']==DIRECT_PROTOCOL:
                            arrays[label+'_native_FMM_phi']=phi.copy();row['native_FMM_potential_max_error_au']=potential_error(cavity,centres,q,phi)
                            phi=direct_potential(cavity,centres,q)
                        perr=potential_error(cavity,centres,q,phi);expected=np.zeros_like(psi);expected[0]=np.sqrt(4*np.pi)*q;serr=float(np.max(abs(psi-expected)))
                        row.update(potential_max_error_au=perr,psi_max_error=serr,charge_sum_e=math.fsum(q),preparation_wall_seconds=time.monotonic()-native_start)
                        if perr>TOL['potential_au'] or serr>TOL['psi']:raise InvalidArtifact('source potential/integral gate failed')
                        state=pyddx.State(model,np.asfortranarray(psi),phi);row['forward_solve_started']=True;solve_start=time.monotonic();state.solve(tol=TOL['solver'])
                        energy=float(state.energy());x=np.array(state.x);contract=.5*math.fsum((psi*x).ravel());err=abs(energy-contract)*HA_TO_KCAL
                        if not state.is_solved or not np.isfinite(energy) or not np.isfinite(x).all():raise InvalidArtifact('native PCM did not return a finite solved state')
                        arrays[label+'_phi']=phi;arrays[label+'_psi']=psi;arrays[label+'_x']=x;psis[label]=psi;solutions[label]=x
                        row.update(status='computed',energy_hartree=energy,energy_kcal_mol=energy*HA_TO_KCAL,contraction_error_kcal=err,
                                   contraction_pass=err<=TOL['energy_kcal'],passive_energy_pass=energy*HA_TO_KCAL<=TOL['energy_kcal'],
                                   iterations=state.x_n_iter,solve_wall_seconds=time.monotonic()-solve_start)
                    except Exception as exc:row['failure_reason']=str(exc)
                    row['wall_seconds']=time.monotonic()-native_start;result['rows'][label]=row;write_new(attempt/(label+'_result.json'),row)
                    print(group['group_id'],label,row['status'],flush=True);state=None
                if all(x in solutions for x in ('Ca','La')):
                    result['reciprocity_error_kcal']=.5*abs(math.fsum((psis['Ca']*solutions['La']).ravel())-math.fsum((psis['La']*solutions['Ca']).ravel()))*HA_TO_KCAL
                result['status']='complete' if all(r['status']=='computed' for r in result['rows'].values()) else 'incomplete'
            except Exception as exc:result['failure_reason']=str(exc)
            np.savez_compressed(attempt/'arrays.npz',**arrays);result['arrays']=record(attempt/'arrays.npz')
            result.update(wall_seconds=time.monotonic()-setup,peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,slurm_job_id=os.environ['SLURM_JOB_ID'])
            write_new(attempt/'result.json',result);groups[group['group_id']]=result
            state=None;model=None;arrays.clear();psis.clear();solutions.clear();gc.collect()
    r=dict(protocol=m['protocol'],manifest=record(path),groups=groups,complete=all(g['status']=='complete' for g in groups.values()),
           wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
           requested_forward_solves=20,actual_forward_solve_starts=sum(row['forward_solve_started'] for g in groups.values() for row in g['rows'].values()),
           new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,r);return r


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output);r={'status':'prepared','groups':len(r['groups']),'solves':r['requested_forward_solves']}
    elif a.command=='dry-run':r=validate(a.manifest);r={'status':'pass','groups':len(r['groups']),'solves':r['requested_forward_solves']}
    else:r=execute(a.manifest,a.output);r={k:v for k,v in r.items() if k!='groups'}
    print(json.dumps(r))

if __name__=='__main__':main()
