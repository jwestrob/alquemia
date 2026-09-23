"""Versioned SLSQP objective-unit scaling for the same four physical angles."""
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
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from scipy.optimize import minimize

from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, paired, read_json, record, verify, write_new, xyz
import adaptive_angular_proposals as angular
from adaptive_force_diagnostic import endpoint, project
from accommodation_proposals import MACEProposal
from accommodation_nonlinear import WarmGPU
from accommodation_torsion_profiles import geometry_check
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics

PROTOCOL='common_four_angular_native_OMOL_SLSQP_hartree_units_v2'
SCALE=EV_TO_KCAL/HA_TO_KCAL
PILOT=('1H4I','4MAE','mmol_1770-pqq-la_model','q9z4j7-pqq-la_model')
SETTINGS={**angular.SETTINGS,'optimizer_objective_unit':'hartree_equivalent',
          'objective_eV_multiplier':SCALE,'optimizer_ftol':angular.SETTINGS['ftol_eV']*SCALE,
          'numerical_change':'same_positive_scale_for_shifted_objective_gradient_and_ftol',
          'additional_starts_or_fallbacks':0}


def parents(paths):
    rows=[];seen=set();reference=None
    for path in paths:
        angular.validate(path);m=read_json(path)
        if reference is None:reference=m
        if any(m[k]!=reference[k] for k in ('source_manifest','diagnostic','model','software','orca','gpu_python','cpu_python')):
            raise InvalidArtifact('angular archives do not share exact source/scientific settings')
        for t in m['tasks']:
            if t['task_id'] in seen:raise InvalidArtifact('duplicate source endpoint')
            seen.add(t['task_id']);rows.append((Path(path).resolve(),m,t))
    original=read_json(verify(reference['source_manifest']));order=[c['case_id'] for c in original['cases']]
    if len(order)!=30 or {(t['case_id'],t['metal']) for _,_,t in rows}!={(c,z) for c in order for z in ('Ca','La')}:
        raise InvalidArtifact('all original30 paired sources must be explicit')
    return reference,original,rows


def origin_record(path,m,t):
    row,data=endpoint(path,m,t,'origin')
    if data is None:raise InvalidArtifact('actual origin force unavailable')
    rp=verify(row['proposal_receipt']);r=read_json(rp)['origin']
    if any(r['full_q']) or xyz(verify(r['coordinate']))!=xyz(verify(t['xyz'])):
        raise InvalidArtifact('origin geometry is not exact q0')
    return {'proposal_receipt':row['proposal_receipt'],'point':r}


class FirstTrial(Exception):pass


def first_trial(kin,t,origin):
    captured={};g=np.asarray(origin['point']['gradient_kcal_mol_rad'])/HA_TO_KCAL
    def objective(q):
        if np.any(q!=0):
            captured['q']=np.asarray(q).tolist()
            raise FirstTrial()
        return 0.,g.copy()  # exact energy difference and archived analytic gradient at q0 only
    try:
        opt=minimize(objective,np.zeros(4),jac=True,method='SLSQP',bounds=[(-.8,.8)]*4,
                     constraints=[{'type':'ineq','fun':lambda q:angular.constraints(kin,t,q)[0],
                                   'jac':lambda q:angular.constraints(kin,t,q)[1]}],
                     options={'maxiter':200,'ftol':SETTINGS['optimizer_ftol']})
        if 'q' not in captured:captured.update(q=np.zeros(4).tolist(),termination=str(opt.message))
    except FirstTrial:pass
    q=np.asarray(captured['q']);full=angular.full_q(t,q)
    checks=geometry_check(kin,full,[a[0] for a in xyz(verify(t['xyz']))])
    return {'active_q_radian':q.tolist(),'physical_status':angular.physical_status(kin,t,q),
            'geometry_checks':checks,'nonzero_objective_evaluations':0,
            'origin_provenance':origin['proposal_receipt']}


def probe(angular_manifest,agreement,output):
    _,_,sources=parents(angular_manifest);rows=[]
    for path,m,t in sources:
        o=origin_record(path,m,t);kin=Kinematics(read_json(verify(t['mapping']))['context'])
        rows.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
                     **first_trial(kin,t,o)})
    result={'protocol_id':PROTOCOL,'settings':SETTINGS,'agreement':record(agreement),
            'sources':[record(p) for p in angular_manifest],'rows':rows,'denominator':60,
            'geometry_guard_pass':sum(r['geometry_checks']['pass'] for r in rows),
            'inside_final_domain':sum(r['physical_status']['physical_feasible'] for r in rows),
            'new_molecular_calls':0,'nonzero_point_energy_or_gradient_returned':False,
            'scipy_version':scipy.__version__,'scipy_wrapper':record(scipy_slsqp.__file__),
            'scipy_kernel':record(scipy_kernel.__file__)}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','settings')}


def selected_cases(original,stage):
    order=[c['case_id'] for c in original['cases']]
    if stage=='all30':return order
    raise InvalidArtifact('unknown fixed population')


def prepare(angular_manifest,agreement,output,stage='all30'):
    parent,original,sources=parents(angular_manifest);cases=selected_cases(original,stage);tasks=[]
    for cid in cases:
        for z in ('Ca','La'):
            path,m,t=next(r for r in sources if (r[2]['case_id'],r[2]['metal'])==(cid,z))
            new=copy.deepcopy(t);new['prior_angular_manifest']=record(path);new['origin_reuse']=origin_record(path,m,t)
            tasks.append(new)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dst=impl/p.name;shutil.copyfile(p,dst);pins[p.name]=record(dst)
    m={**{k:parent[k] for k in ('source_manifest','diagnostic','model','software','orca','cpu_python','gpu_python',
                              'cpu_executable','gpu_executable','resources')},
       'protocol_id':PROTOCOL,'settings':SETTINGS,'stage':stage,'agreement':record(agreement),
       'prior_angular_manifests':[record(p) for p in angular_manifest],
       'cases':[next(c for c in original['cases'] if c['case_id']==cid) for cid in cases],
       'declared_case_ids':cases,'tasks':tasks,'implementation':pins,
       'optimizer_software':{'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),
                             'kernel':record(scipy_kernel.__file__)},
       'pilot_task_ids':[c+'__'+z for c in PILOT for z in ('Ca','La')],
       'new_optimizer_starts':len(tasks),'new_DFT_calls':0,'new_GFN2_calls_in_this_adapter':0,
       'source_states_changed':False,'reference':None,'production_changed':False}
    write_new(out/'manifest.json',m);r=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',r);return r


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS:raise InvalidArtifact('numerical version changed')
    original=read_json(verify(m['source_manifest']));cases=selected_cases(original,m['stage'])
    if m['declared_case_ids']!=cases or [(t['case_id'],t['metal']) for t in m['tasks']]!=[(c,z) for c in cases for z in ('Ca','La')]:
        raise InvalidArtifact('fixed paired population changed')
    if m['pilot_task_ids']!=[c+'__'+z for c in PILOT for z in ('Ca','La')]:
        raise InvalidArtifact('fixed pilot task selector changed')
    verify(m['agreement']);verify(m['diagnostic'])
    for k in ('software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('installed SciPy version differs')
    archives={a['path']:read_json(verify(a)) for a in m['prior_angular_manifests']}
    for t in m['tasks']:
        old=next(o for o in archives[t['prior_angular_manifest']['path']]['tasks'] if o['task_id']==t['task_id'])
        if {k:v for k,v in t.items() if k not in ('prior_angular_manifest','origin_reuse')}!=old:
            raise InvalidArtifact('four-angle task source/physics changed')
        origin=t['origin_reuse'];r=read_json(verify(origin['proposal_receipt']))
        if r['origin']!=origin['point'] or r['manifest']!=t['prior_angular_manifest']:
            raise InvalidArtifact('origin replay pin changed')
        native=read_json(verify(origin['point']['MACE']));req=read_json(verify(native['request']))
        if (native['status']!='complete' or native['model']!=m['model'] or
                req['charge']!=t['charge'] or req['multiplicity']!=t['multiplicity'] or
                req['xyz']!=origin['point']['coordinate'] or native['forces']!=origin['point']['forces']):
            raise InvalidArtifact('cached origin method/state differs')
        verify(native['forces'])
        if xyz(verify(req['xyz']))!=xyz(verify(t['xyz'])) or any(origin['point']['full_q']):
            raise InvalidArtifact('cached origin coordinates differ')
        kin=Kinematics(read_json(verify(t['mapping']))['context'])
        angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
    for cid in cases:
        ca,la=[next(t for t in m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
        if ca['active_mode_ids']!=la['active_mode_ids'] or ca['selector']!=la['selector']:
            raise InvalidArtifact('common source selector differs')
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
    return {'status':'prepared_dry_run_pass','manifest':record(manifest),'cases':len(cases),'tasks':len(m['tasks']),
            'reusable_exact_origins':len(m['tasks']),'new_molecular_calls':0}


class CompletionProposal(MACEProposal):
    def evaluate(self,q,purpose):
        q=np.asarray(q,dtype=float)
        if q.shape!=(4,) or np.any(q!=0):return super().evaluate(q,purpose)
        key=cache_key({'manifest':self.manifest,'task':self.task['task_id'],'q':q.tolist()})
        d=self.directory/'evaluations'/key;p=d/'result.json'
        if not p.exists():
            d.mkdir(parents=True,exist_ok=False)
            r={**self.task['origin_reuse']['point'],'manifest':self.manifest,'cache_key':key,
               'reused_scientific_origin':self.task['origin_reuse']['proposal_receipt']}
            write_new(p,r)
        r=read_json(p)
        if r['manifest']!=self.manifest or r['cache_key']!=key or r['status']!='complete':
            raise InvalidArtifact('origin reuse cache differs')
        self.requests.append({'purpose':purpose,'result':record(p),'reused':True,'scientific_origin_reused':True})
        return r


def optimize(manifest,task,gpu):
    ev=CompletionProposal(manifest,task,gpu);dest=ev.directory/'result.json'
    if dest.exists():
        if read_json(dest)['manifest']!=record(manifest):raise InvalidArtifact('result manifest differs')
        return record(dest)
    began=time.monotonic();trace=[];accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],'manifest':record(manifest),
            'status':'unavailable','origin':None,'proposal':None,'optimizer':None,
            'unconstrained_minimum_claimed':False,'composite_gradient':None,'optimizer_objective_scale':SCALE}
    def evaluate(q,purpose):
        trial={'active_q_radian':np.asarray(q).tolist(),'purpose':purpose,
               **angular.physical_status(ev.kin,task,q),'MACE_complete':False};trace.append(trial)
        point=ev.evaluate(q,purpose);trial['MACE_complete']=True;return point
    try:
        origin=evaluate(np.zeros(4),'origin');result['origin']=origin
        replay=(origin['MACE_eV']-task['q0']['components']['MACE_eV'])*EV_TO_KCAL
        result['origin_replay_kcal_mol']=replay
        if abs(replay)>SETTINGS['origin_replay_tolerance_kcal_mol']:raise InvalidArtifact('native origin replay failed')
        def objective(q):
            point=evaluate(q,'SLSQP_trial')
            return ((point['MACE_eV']-origin['MACE_eV'])*SCALE,
                    np.asarray(point['gradient_kcal_mol_rad'])/HA_TO_KCAL)
        opt=minimize(objective,np.zeros(4),jac=True,method='SLSQP',bounds=[(-.8,.8)]*4,
                     constraints=[{'type':'ineq','fun':lambda q:angular.constraints(ev.kin,task,q)[0],
                                   'jac':lambda q:angular.constraints(ev.kin,task,q)[1]}],
                     callback=lambda q:accepted.append(np.asarray(q).tolist()),
                     options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['optimizer_ftol']})
        result['optimizer']={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),
                             'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev),
                             'ftol_hartree_equivalent':SETTINGS['optimizer_ftol']}
        check=angular.final_geometry(ev.kin,task,opt.x,[a[0] for a in ev.atoms]);result['final_geometry']=check
        candidate=evaluate(opt.x,'final_candidate');result['proposal']=candidate
        if not opt.success:raise InvalidArtifact('SLSQP did not converge: '+str(opt.message))
        work=candidate['MACE_eV']-origin['MACE_eV']
        if work>SETTINGS['maximum_final_energy_increase_eV']:raise InvalidArtifact('native final energy increased')
        forces=np.load(verify(candidate['forces']),allow_pickle=False)
        _,_,raw,normed,_=project(ev.kin.data,candidate['full_q'],forces)
        active=set(task['active_indices']);angles=[i for i,v in enumerate(ev.kin.modes) if v['unit']=='radian']
        result.update(status='proposal_available',MACE_proposal_work_kcal_mol=float(work*EV_TO_KCAL),
                      boundary_flag=bool(check['boundary_heavy_atoms'] or check['angular_boundary_flag']),
                      residual={'force_hamiltonian':'native_vacuum_OMOL','all_mode_ids':[v['id'] for v in ev.kin.modes],
                                'raw_gradient_kcal_per_unit':raw.tolist(),'normalized_gradient_kcal_mol_A':normed.tolist(),
                                'selected_max_load':float(max(abs(normed[i]) for i in active)),
                                'omitted_angular_max_load':float(max(abs(normed[i]) for i in angles if i not in active))})
    except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(trials=trace,accepted_iterations=accepted,requests=ev.requests,
                  infeasible_trial_requests=sum(not t['physical_feasible'] for t in trace),
                  infeasible_completed_MACE_requests=sum(not t['physical_feasible'] and t['MACE_complete'] for t in trace),
                  maximum_trial_heavy_extent_A=max((t['maximum_heavy_displacement_A'] for t in trace),default=None),
                  wall_seconds=time.monotonic()-began,job_id=os.environ.get('SLURM_JOB_ID'))
    write_new(dest,result);return record(dest)


def execute(manifest,pool_collection,phase):
    validate(manifest);pool=angular.pool_gate(manifest,pool_collection)
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
            int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):raise InvalidArtifact('existing32CPU/200000MiB/GPU allocation required')
    m=read_json(manifest);root=Path(manifest).parent;lock=(root/'proposal.lock').open('a+')
    if m['stage']!='all30' or phase not in ('pilot8','complete60'):
        raise InvalidArtifact('execution requires the full original30 manifest and explicit phase')
    chosen=set(m['pilot_task_ids']) if phase=='pilot8' else {t['task_id'] for t in m['tasks']}
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);began=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for t in m['tasks']:
            if t['task_id'] not in chosen:continue
            pin=optimize(manifest,t,gpu);results.append(pin)
            print(json.dumps({'task_id':t['task_id'],'status':read_json(verify(pin))['status']}),flush=True)
    except Exception as exc:error=str(exc)
    finally:
        try:
            if gpu is not None:gpu.close()
        except Exception as exc:error=(error+'; ' if error else '')+str(exc)
        elapsed=time.monotonic()-began
        write_new(root/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),
                  {'manifest':record(manifest),'shared_pool_pilot':pool,'results':results,'error':error,
                   'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,
                   'job_id':os.environ['SLURM_JOB_ID'],'GPU_command':gpu.command if gpu else None,
                   'phase':phase,'selected_task_ids':sorted(chosen),'new_DFT_calls':0,'new_GFN2_calls':0})
    if error:raise InvalidArtifact(error)
    return {'status':'proposal_phase_complete','endpoint_denominator':len(m['tasks']),
            'phase_task_count':len(chosen),'results':results}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;rows=[]
    for t in m['tasks']:
        path=root/'proposals'/t['task_id']/'result.json';r=read_json(path) if path.exists() else None
        if r and r['manifest']!=record(manifest):raise InvalidArtifact('candidate source differs')
        ok=bool(r and r['status']=='proposal_available')
        if ok:
            kin=Kinematics(read_json(verify(t['mapping']))['context']);q=np.asarray(r['proposal']['full_q'])
            if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('unselected mode moved')
            atoms=xyz(verify(r['proposal']['coordinate']));angular.final_geometry(kin,t,q[t['active_indices']],[a[0] for a in atoms])
            if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate mapping differs')
        rows.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
                     'status':'candidate_available' if ok else 'unavailable','reason':r.get('reason') if r else 'not_run',
                     'proposal_receipt':record(path) if r else None,'source_preparation':t['source_preparation'],
                     'mapping':t['mapping'],'charge':t['charge'],'multiplicity':t['multiplicity'],
                     'common_active_mode_ids':t['active_mode_ids'],'origin_xyz':t['xyz'],'origin_components':t['q0']['components'],
                     'candidate':r['proposal'] if ok else None,'candidate_GFN2_vacuum_hartree':None,
                     'candidate_GFN2_ALPB_hartree':None,'composite_candidate_energy':None,
                     'residual':r.get('residual') if r else None,'boundary_flag':r.get('boundary_flag') if r else None})
    cases=[{'case_id':cid,'status':'ready_for_common_pool' if all(r['status']=='candidate_available' for r in rows if r['case_id']==cid) else 'unavailable',
            'selected_geometry':None,'score':None} for cid in m['declared_case_ids']]
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'source_manifest':m['source_manifest'],
            'cases':cases,'endpoints':rows,'case_denominator':len(cases),'endpoint_denominator':len(rows),
            'available_candidates':sum(r['status']=='candidate_available' for r in rows),
            'new_DFT_calls':0,'new_GFN2_calls_in_adapter':0,'common_pool_scoring_performed':False,
            'production_changed':False,'reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('cases','endpoints')}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for name in ('probe','prepare'):
        q=s.add_parser(name);q.add_argument('--angular-manifest',action='append',required=True)
        q.add_argument('--agreement',required=True);q.add_argument('--output',required=True)
        if name=='prepare':q.add_argument('--stage',choices=('all30',),default='all30')
    q=s.add_parser('dry-run');q.add_argument('--manifest',required=True)
    q=s.add_parser('execute');q.add_argument('--manifest',required=True);q.add_argument('--pool-collection',required=True)
    q.add_argument('--phase',choices=('pilot8','complete60'),required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op');fn={'probe':probe,'prepare':prepare,'dry-run':validate,'execute':execute,'collect':collect}[op]
    print(json.dumps(fn(**a),indent=2))


if __name__=='__main__':main()
