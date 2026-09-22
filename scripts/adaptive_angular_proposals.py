"""Four common angular donor coordinates, constrained native-MACE proposals."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import traceback
import numpy as np
from affordable_common import InvalidArtifact, paired, read_json, record, verify, write_new, xyz
from accommodation_nonlinear import WarmGPU
from accommodation_proposals import MACEProposal, PLM
from accommodation_torsion_profiles import geometry_check
from adaptive_force_diagnostic import endpoint, preview, project
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics

PROTOCOL = 'common_four_angular_native_OMOL_SLSQP_proposals_v1'
PILOT4 = ('1H4I', '4MAE', *PLM)
CASES = PILOT4  # Original import/API retained for the completed pilot.
POPULATIONS = ('pilot4', 'remaining26')
SETTINGS = {'optimizer': 'SLSQP', 'maxiter': 200, 'ftol_eV': 1e-9,
            'maximum_physical_heavy_displacement_A': .8, 'angular_bound_radian': .8,
            'final_physical_tolerance_A': 1e-7, 'angular_tolerance_radian': 1e-12,
            'origin_replay_tolerance_kcal_mol': .01, 'maximum_final_energy_increase_eV': 1e-7,
            'physical_boundary_margin_A': 1e-5, 'angular_boundary_margin_radian': .02,
            'selected_angular_modes': 4, 'starts_per_endpoint': 1,
            'selector': 'alternating_differential_individual_heavy_norm_MGS_v1',
            'trial_constraint_policy': 'infeasible_SLSQP_trials_recorded_final_candidate_must_be_feasible'}


def full_q(task, active):
    active = np.asarray(active, dtype=float)
    if active.shape != (4,) or not np.isfinite(active).all():
        raise InvalidArtifact('expected four finite angular primitives')
    q = np.zeros(task['mode_count']); q[task['active_indices']] = active
    return q


def constraints(kin, task, active):
    """Existing physical SLSQP inequality/Jacobian, restricted to selected columns."""
    q = full_q(task, active); p, _, j, _ = kin.evaluate(q)
    delta = p[kin.heavy] - kin.positions[kin.heavy]
    values = SETTINGS['maximum_physical_heavy_displacement_A']**2 - np.sum(delta*delta, axis=1)
    jac = -2*np.einsum('ij,mij->im', delta, j[task['active_indices']][:, kin.heavy])
    return values, jac


def physical_status(kin, task, active):
    q = full_q(task, active); p, _, _, _ = kin.evaluate(q)
    distance = np.linalg.norm(p[kin.heavy] - kin.positions[kin.heavy], axis=1)
    allowed = SETTINGS['maximum_physical_heavy_displacement_A']
    return {'maximum_heavy_displacement_A': float(distance.max()),
            'physical_feasible': bool(distance.max() <= allowed+SETTINGS['final_physical_tolerance_A']),
            'boundary_heavy_atoms': int(sum(distance >= allowed-SETTINGS['physical_boundary_margin_A'])),
            'angular_boundary_flag': bool(np.min(.8-np.abs(active)) < SETTINGS['angular_boundary_margin_radian'])}


def final_geometry(kin, task, active, symbols):
    active = np.asarray(active); q = full_q(task, active)
    if np.max(np.abs(active)) > .8+SETTINGS['angular_tolerance_radian']:
        raise InvalidArtifact('final angular ceiling violated')
    status = physical_status(kin, task, active)
    if not status['physical_feasible']: raise InvalidArtifact('final physical displacement bound violated')
    checks = geometry_check(kin, q, symbols)
    if not checks['pass']: raise InvalidArtifact('final source/cap/bond or overlap check failed')
    moving = {i for k in task['active_indices'] for i in kin.modes[k]['moving_indices']}
    fixed = [i for i in range(len(kin.positions)) if i not in moving]
    p = kin.evaluate(q)[0]
    if not np.allclose(p[fixed], kin.positions[fixed], atol=1e-12, rtol=0):
        raise InvalidArtifact('nonselected physical atoms moved')
    return status | {'mapping_checks': checks, 'all_other_physical_atoms_fixed': True,
                     'metal_water_PQQ_scaffold_modes_added': False}


def population_cases(parent, population='pilot4'):
    """Frozen source order, never score/label/availability-based selection."""
    ids = [c['case_id'] for c in parent['cases']]
    if len(ids) != 30 or len(set(ids)) != 30 or not set(PILOT4) <= set(ids):
        raise InvalidArtifact('expected original consumed30 source with pilot4')
    if population == 'pilot4':
        return list(PILOT4)
    if population == 'remaining26':
        return [cid for cid in ids if cid not in PILOT4]
    raise InvalidArtifact('unknown adaptive population: '+str(population))


def prepare(source, diagnostic, agreement, output, population='pilot4'):
    parent = read_json(source); diag = read_json(diagnostic)
    if record(source) not in [a['manifest'] for a in diag['archives']]:
        raise InvalidArtifact('diagnostic has another source')
    cases = population_cases(parent, population)
    tasks = []; choices = []
    for cid in cases:
        ts = [next(t for t in parent['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
        paired(verify(ts[1]['xyz']), verify(ts[0]['xyz']), ts[1]['charge'], ts[0]['charge'])
        if (read_json(verify(ts[0]['mapping'])) != read_json(verify(ts[1]['mapping'])) or
                ts[0]['source_preparation'] != ts[1]['source_preparation']):
            raise InvalidArtifact('common origin mapping/preparation differs')
        projections = [endpoint(source,parent,t,'origin')[1] for t in ts]
        if any(x is None for x in projections): raise InvalidArtifact('origin force unavailable')
        ca,la = projections
        choice = preview(ca['mapping']['modes'],ca['vectors'],ca['normalized'],la['normalized'])
        old = next(r for r in diag['pairs'] if r['case_id']==cid)
        if choice != old['selection_preview']: raise InvalidArtifact('frozen selector does not replay')
        ids = [r['id'] for r in choice['selected']]
        if len(ids)!=4 or len(set(ids))!=4: raise InvalidArtifact('four independent angular modes unavailable')
        choices.append({'case_id':cid,'selection':choice})
        for t in ts:
            if t['q0_status']!='available': raise InvalidArtifact('composite origin unavailable')
            kin = Kinematics(read_json(verify(t['mapping']))['context'])
            names = [m['id'] for m in kin.modes]
            new = copy.deepcopy(t)
            new.update(active_mode_ids=ids, active_indices=[names.index(name) for name in ids],
                       active_roles=['adaptive_physical_angular']*4,
                       original_active_mode_ids=t['active_mode_ids'], selector=choice)
            final_geometry(kin,new,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
            tasks.append(new)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        target=impl/p.name;shutil.copyfile(p,target);pins[p.name]=record(target)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'source_manifest':record(source),
       'diagnostic':record(diagnostic),'agreement':record(agreement),'tasks':tasks,'selections':choices,
       'cases':[c for cid in cases for c in parent['cases'] if c['case_id']==cid],
       'declared_case_ids':cases,'implementation':pins,
       **{k:parent[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable')},
       'resources':{'GPU_cpus':32,'GPUs':1,'GPU_host_mem_MiB':200000},
       'new_optimizer_starts':len(tasks),'new_DFT_calls':0,'new_GFN2_calls_in_this_adapter':0,
       'source_states_changed':False,'reference':None,'production_changed':False,
       'execution_requires_completed_shared_pool_pilot':True}
    # Keep default pilot manifests structurally compatible with the original.
    if population != 'pilot4':
        m.update(population=population, excluded_case_ids=list(PILOT4),
                 population_rule='original30_source_order_excluding_fixed_pilot4')
    path=out/'manifest.json';write_new(path,m);check=validate(path);write_new(out/'PREFLIGHT.json',check)
    return check


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS: raise InvalidArtifact('adaptive protocol/settings differ')
    parent=read_json(verify(m['source_manifest']));diag=read_json(verify(m['diagnostic']));verify(m['agreement'])
    cases=population_cases(parent,m.get('population','pilot4'))
    if (m['declared_case_ids']!=cases or len(m['tasks'])!=2*len(cases) or
            m['new_optimizer_starts']!=2*len(cases)):
        raise InvalidArtifact('declared population/endpoint scope differs')
    if m['source_manifest'] not in [a['manifest'] for a in diag['archives']]:
        raise InvalidArtifact('diagnostic has another source')
    if m['cases'] != [c for cid in cases for c in parent['cases'] if c['case_id']==cid]:
        raise InvalidArtifact('source case records changed')
    if (m.get('population','pilot4')=='remaining26' and
            (m.get('excluded_case_ids')!=list(PILOT4) or
             m.get('population_rule')!='original30_source_order_excluding_fixed_pilot4')):
        raise InvalidArtifact('fixed pilot exclusions differ')
    if {(t['case_id'],t['metal']) for t in m['tasks']}!={(c,z) for c in cases for z in ('Ca','La')}:
        raise InvalidArtifact('paired task membership differs')
    if any(m[k]!=parent[k] for k in ('model','software','orca','cpu_python','gpu_python')):
        raise InvalidArtifact('source electronic model or environment differs')
    for key in ('software','orca','cpu_executable','gpu_executable'):verify(m[key])
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        old=next(v for v in parent['tasks'] if v['task_id']==t['task_id'])
        for key in ('xyz','mapping','source_preparation','charge','multiplicity','q0'):
            if t[key]!=old[key]:raise InvalidArtifact('adaptive source state changed: '+key)
        kin=Kinematics(read_json(verify(t['mapping']))['context'])
        choice=next(r['selection_preview'] for r in diag['pairs'] if r['case_id']==t['case_id'])
        if t['selector']!=choice or t['active_mode_ids']!=[s['id'] for s in choice['selected']]:
            raise InvalidArtifact('frozen common selector differs')
        if ([kin.modes[i]['id'] for i in t['active_indices']]!=t['active_mode_ids'] or
                any(kin.modes[i]['kind'] not in ('sidechain_torsion','peptide_crankshaft') for i in t['active_indices'])):
            raise InvalidArtifact('unsupported or mismapped selected mode')
        final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
    for cid in cases:
        ca,la=[next(t for t in m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
        if ca['active_mode_ids']!=la['active_mode_ids'] or ca['active_indices']!=la['active_indices']:
            raise InvalidArtifact('endpoint active spaces differ')
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
    return {'status':'prepared_dry_run_pass','manifest':record(manifest),'cases':len(cases),'tasks':len(m['tasks']),
            'molecular_calls':0,'optimizer_executed':False,'new_DFT_calls':0}


def pool_gate(manifest, pool_collection):
    m=read_json(manifest);pool=read_json(pool_collection);pm=read_json(verify(pool['manifest']))
    if (pool['protocol_id']!='nikasha_common_geometry_native_OMOL_GFN2_ALPB_v1' or
        pm['source_manifest']!=m['source_manifest'] or pm['model']!=m['model'] or
        {c['case_id'] for c in pool['cases']}!=set(PILOT4) or
        len(pool['cases'])!=len(PILOT4) or
        any(c['pool']['status']!='available' for c in pool['cases'])):
        raise InvalidArtifact('complete compatible four-case shared-pool pilot required')
    return record(pool_collection)


def optimize(manifest, task, gpu):
    from scipy.optimize import minimize
    ev=MACEProposal(manifest,task,gpu);dest=ev.directory/'result.json'
    if dest.exists():
        old=read_json(dest)
        if old['manifest']!=record(manifest):raise InvalidArtifact('existing adaptive result source differs')
        return record(dest)
    began=time.monotonic();trace=[];accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],
            'manifest':record(manifest),'status':'unavailable','origin':None,'proposal':None,
            'optimizer':None,'unconstrained_minimum_claimed':False,'composite_gradient':None}
    def evaluate(q,purpose):
        status=physical_status(ev.kin,task,q)
        trial={'active_q_radian':np.asarray(q).tolist(),'purpose':purpose,**status,'MACE_complete':False}
        trace.append(trial)
        point=ev.evaluate(q,purpose);trial['MACE_complete']=True
        return point
    try:
        origin=evaluate(np.zeros(4),'origin');result['origin']=origin
        replay=(origin['MACE_eV']-task['q0']['components']['MACE_eV'])*EV_TO_KCAL
        result['origin_replay_kcal_mol']=replay
        if abs(replay)>SETTINGS['origin_replay_tolerance_kcal_mol']:
            raise InvalidArtifact('native origin replay failed')
        def objective(q):
            point=evaluate(q,'SLSQP_trial')
            return point['MACE_eV']-origin['MACE_eV'],np.asarray(point['gradient_kcal_mol_rad'])/EV_TO_KCAL
        opt=minimize(objective,np.zeros(4),jac=True,method='SLSQP',bounds=[(-.8,.8)]*4,
            constraints=[{'type':'ineq','fun':lambda q:constraints(ev.kin,task,q)[0],
                          'jac':lambda q:constraints(ev.kin,task,q)[1]}],
            callback=lambda q:accepted.append(np.asarray(q).tolist()),
            options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['ftol_eV']})
        result['optimizer']={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),
                             'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev)}
        # No clipping, coordinate projection, candidate substitution or rescue.
        check=final_geometry(ev.kin,task,opt.x,[a[0] for a in ev.atoms]);result['final_geometry']=check
        candidate=evaluate(opt.x,'final_candidate');result['proposal']=candidate
        if not opt.success:raise InvalidArtifact('SLSQP did not converge: '+str(opt.message))
        work=candidate['MACE_eV']-origin['MACE_eV']
        if work>SETTINGS['maximum_final_energy_increase_eV']:
            raise InvalidArtifact('native final energy increased')
        force=np.load(verify(candidate['forces']),allow_pickle=False)
        _,_,raw,normalized,lengths=project(ev.kin.data,candidate['full_q'],force)
        active=set(task['active_indices']);angular=[i for i,m in enumerate(ev.kin.modes) if m['unit']=='radian']
        result.update(status='proposal_available',MACE_proposal_work_kcal_mol=float(work*EV_TO_KCAL),
            boundary_flag=bool(check['boundary_heavy_atoms'] or check['angular_boundary_flag']),
            residual={'force_hamiltonian':'native_vacuum_OMOL',
                      'all_mode_ids':[m['id'] for m in ev.kin.modes],
                      'raw_gradient_kcal_per_unit':raw.tolist(),'normalized_gradient_kcal_mol_A':normalized.tolist(),
                      'selected_max_load':float(max(abs(normalized[i]) for i in active)),
                      'omitted_angular_max_load':float(max(abs(normalized[i]) for i in angular if i not in active))})
    except Exception as exc:
        result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(trials=trace,accepted_iterations=accepted,requests=ev.requests,
                  infeasible_trial_requests=sum(not t['physical_feasible'] for t in trace),
                  infeasible_completed_MACE_requests=sum(not t['physical_feasible'] and t['MACE_complete'] for t in trace),
                  maximum_trial_heavy_extent_A=max((t['maximum_heavy_displacement_A'] for t in trace),default=None),
                  wall_seconds=time.monotonic()-began,job_id=os.environ.get('SLURM_JOB_ID'))
    write_new(dest,result);return record(dest)


def execute(manifest,pool_collection):
    validate(manifest);pool=pool_gate(manifest,pool_collection)
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
            int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):
        raise InvalidArtifact('existing32CPU/200000MiB/GPU allocation required')
    m=read_json(manifest);root=Path(manifest).parent;lock=(root/'proposal.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    began=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for task in m['tasks']:
            pin=optimize(manifest,task,gpu);results.append(pin);r=read_json(verify(pin))
            print(json.dumps({'task_id':task['task_id'],'status':r['status'],
                              'infeasible_trials':r['infeasible_trial_requests']}),flush=True)
    except Exception as exc:error=str(exc)
    finally:
        try:
            if gpu is not None:gpu.close()
        except Exception as exc:error=(error+'; ' if error else '')+str(exc)
        elapsed=time.monotonic()-began
        receipt={'manifest':record(manifest),'shared_pool_pilot':pool,'results':results,'error':error,
                 'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,
                 'job_id':os.environ['SLURM_JOB_ID'],'GPU_command':gpu.command if gpu else None,
                 'new_DFT_calls':0,'new_GFN2_calls':0}
        write_new(root/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),receipt)
    if error:raise InvalidArtifact(error)
    return {'status':'proposal_phase_complete','endpoint_denominator':len(m['tasks']),'results':results,
            'common_pool_scoring':'required_next_stage_not_executed_here'}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;endpoints=[]
    for t in m['tasks']:
        path=root/'proposals'/t['task_id']/'result.json';r=read_json(path) if path.exists() else None
        if r and r['manifest']!=record(manifest):raise InvalidArtifact('candidate belongs to another manifest')
        ok=bool(r and r['status']=='proposal_available')
        row={'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
             'status':'candidate_available' if ok else 'unavailable',
             'reason':r.get('reason') if r else 'not_run','proposal_receipt':record(path) if r else None,
             'source_preparation':t['source_preparation'],'mapping':t['mapping'],
             'charge':t['charge'],'multiplicity':t['multiplicity'],'common_active_mode_ids':t['active_mode_ids'],
             'origin_xyz':t['xyz'],'origin_components':t['q0']['components'],
             'candidate':r['proposal'] if ok else None,'candidate_GFN2_vacuum_hartree':None,
             'candidate_GFN2_ALPB_hartree':None,'composite_candidate_energy':None,
             'residual':r.get('residual') if r else None,'boundary_flag':r.get('boundary_flag') if r else None}
        if ok:
            atoms=xyz(verify(r['proposal']['coordinate']));kin=Kinematics(read_json(verify(t['mapping']))['context'])
            q=np.asarray(r['proposal']['full_q']);inactive=set(range(len(q)))-set(t['active_indices'])
            if any(q[i]!=0 for i in inactive):raise InvalidArtifact('unselected candidate mode moved')
            final_geometry(kin,t,q[t['active_indices']],[a[0] for a in atoms])
            if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):
                raise InvalidArtifact('candidate coordinates do not reproduce physical map')
        endpoints.append(row)
    cases=[{'case_id':cid,'status':'ready_for_common_pool' if all(e['status']=='candidate_available' for e in endpoints if e['case_id']==cid) else 'unavailable',
            'selected_geometry':None,'score':None} for cid in m['declared_case_ids']]
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'source_manifest':m['source_manifest'],
            'cases':cases,'endpoints':endpoints,'case_denominator':len(cases),'endpoint_denominator':len(m['tasks']),
            'available_candidates':sum(e['status']=='candidate_available' for e in endpoints),
            'new_DFT_calls':0,'new_GFN2_calls_in_adapter':0,'common_pool_scoring_performed':False,
            'production_changed':False,'reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('cases','endpoints')}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='operation',required=True)
    q=sub.add_parser('prepare')
    for key in ('source','diagnostic','agreement','output'):q.add_argument('--'+key,required=True)
    q.add_argument('--population',choices=POPULATIONS,default='pilot4',
                   help='pilot4 preserves the original scope; remaining26 excludes exactly that fixed pilot')
    q=sub.add_parser('dry-run');q.add_argument('--manifest',required=True)
    q=sub.add_parser('execute');q.add_argument('--manifest',required=True);q.add_argument('--pool-collection',required=True)
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation');fn={'prepare':prepare,'dry-run':validate,'execute':execute,'collect':collect}[op]
    print(json.dumps(fn(**a),indent=2))


if __name__=='__main__':main()
