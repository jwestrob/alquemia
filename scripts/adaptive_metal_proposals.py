"""Common metal translations plus the frozen four angular primitives; opt-in."""
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

from affordable_common import InvalidArtifact, cache_key, paired, read_json, record, verify, write_new, xyz
import adaptive_angular_proposals as angular
from adaptive_force_diagnostic import endpoint, project
from accommodation_nonlinear import WarmGPU
from accommodation_proposals import MACEProposal
from accommodation_torsion_profiles import geometry_check
from mace_hybrid import EV_TO_KCAL, write_xyz
from mace_site_kinematics import Kinematics

PROTOCOL = 'common_metal3_four_angular_native_OMOL_SLSQP_proposals_v1'
METAL_IDS = ['metal_x', 'metal_y', 'metal_z']
UNITS = ['angstrom']*3+['radian']*4
GRADIENT_UNITS = ['eV/angstrom']*3+['eV/radian']*4
SETTINGS = {**angular.SETTINGS, 'metal_translation_components':3,
            'translation_component_bound_A':.8,
            'coordinate_order':'metal_xyz_then_unchanged_q0_selected_four_angles',
            'objective_gradient_units':GRADIENT_UNITS}


def full_q(task, active):
    active=np.asarray(active,dtype=float)
    if active.shape!=(7,) or not np.isfinite(active).all():
        raise InvalidArtifact('expected three finite translations and four finite angles')
    q=np.zeros(task['mode_count']);q[task['active_indices']]=active
    return q


def constraints(kin, task, active):
    p,_,j,_=kin.evaluate(full_q(task,active))
    delta=p[kin.heavy]-kin.positions[kin.heavy]
    return (.8**2-np.sum(delta*delta,axis=1),
            -2*np.einsum('ij,mij->im',delta,j[task['active_indices']][:,kin.heavy]))


def physical_status(kin,task,active):
    active=np.asarray(active);p=kin.evaluate(full_q(task,active))[0]
    distance=np.linalg.norm(p[kin.heavy]-kin.positions[kin.heavy],axis=1)
    metal=float(np.linalg.norm(p[kin.metal]-kin.positions[kin.metal]))
    return {'maximum_heavy_displacement_A':float(distance.max()),'metal_displacement_A':metal,
            'physical_feasible':bool(distance.max()<=.8+SETTINGS['final_physical_tolerance_A']),
            'boundary_heavy_atoms':int(sum(distance>=.8-SETTINGS['physical_boundary_margin_A'])),
            'metal_boundary_flag':metal>=.8-SETTINGS['physical_boundary_margin_A'],
            'angular_boundary_flag':bool(np.min(.8-np.abs(active[3:]))<SETTINGS['angular_boundary_margin_radian'])}


def final_geometry(kin,task,active,symbols):
    active=np.asarray(active);q=full_q(task,active)
    if np.max(np.abs(active[:3]))>.8+1e-12 or np.max(np.abs(active[3:]))>.8+1e-12:
        raise InvalidArtifact('final component bound violated (angstrom/radian)')
    status=physical_status(kin,task,active)
    if not status['physical_feasible']:raise InvalidArtifact('final physical displacement bound violated')
    check=geometry_check(kin,q,symbols)
    if not check['pass']:raise InvalidArtifact('final source/cap/bond or overlap check failed')
    moving={i for k in task['active_indices'] for i in kin.modes[k]['moving_indices']}
    fixed=[i for i in range(len(kin.positions)) if i not in moving]
    p=kin.evaluate(q)[0]
    if not np.allclose(p[fixed],kin.positions[fixed],atol=1e-12,rtol=0):
        raise InvalidArtifact('nonselected physical atom moved')
    return status|{'mapping_checks':check,'all_other_physical_atoms_fixed':True,
                   'metal_modes_added':True,'water_PQQ_scaffold_modes_added':False}


def projected_gradient_eV(kin,task,active,forces):
    """Energy derivatives; translation eV/Å, torsion eV/rad, cap chain rules included."""
    _,coords,_,jac=kin.evaluate(full_q(task,active));forces=np.asarray(forces,dtype=float)
    if forces.shape!=coords.shape or not np.isfinite(forces).all():
        raise InvalidArtifact('invalid actual Cartesian force array')
    return -np.einsum('mij,ij->m',jac[task['active_indices']],forces)


def residual_inventory(source):
    m=read_json(source);rows=[]
    for t in m['tasks']:
        r,data=endpoint(source,m,t,'proposal')
        if data is None or r['proposal_status']!='proposal_available':
            raise InvalidArtifact('completed angular candidate/forces unavailable')
        metal=[next(v for v in r['modes'] if v['id']==name) for name in METAL_IDS]
        rows.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
                     'proposal_receipt':r['proposal_receipt'],'forces':r['forces'],'coordinate':r['coordinate'],
                     'gradient_kcal_mol_A':[v['raw_gradient_kcal_per_unit'] for v in metal],
                     'metal_gradient_norm_kcal_mol_A':r['metal_translation_load_norm'],
                     'selected_angular_max_load_kcal_mol_A':r['prior_active_max_load'],
                     'omitted_angular_max_load_kcal_mol_A':r['omitted_angular_max_load'],
                     'angular_pilot_boundary_flag':r['boundary_flag']})
    return {'source_manifest':record(source),'rows':rows,'force_hamiltonian':'native_vacuum_OMOL',
            'comparison_geometry':'each_endpoint_own_final_candidate_not_matched_differential',
            'new_model_calls':0}


def prepare(angular_manifest,agreement,output):
    angular.validate(angular_manifest);parent=read_json(angular_manifest)
    if parent['declared_case_ids']!=list(angular.PILOT4) or len(parent['tasks'])!=8:
        raise InvalidArtifact('only the original four consumed contexts are declared')
    residual=residual_inventory(angular_manifest);tasks=[]
    for old in parent['tasks']:
        kin=Kinematics(read_json(verify(old['mapping']))['context']);names=[m['id'] for m in kin.modes]
        ids=METAL_IDS+old['active_mode_ids'];t=copy.deepcopy(old)
        t.update(active_mode_ids=ids,active_indices=[names.index(n) for n in ids],
                 active_roles=['metal_translation']*3+['adaptive_physical_angular']*4,
                 active_coordinate_units=UNITS,objective_gradient_units=GRADIENT_UNITS,
                 angular_parent_mode_ids=old['active_mode_ids'])
        final_geometry(kin,t,np.zeros(7),[a[0] for a in xyz(verify(t['xyz']))]);tasks.append(t)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    write_new(out/'RESIDUALS.json',residual)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        target=impl/p.name;shutil.copyfile(p,target);pins[p.name]=record(target)
    m={**{k:parent[k] for k in ('source_manifest','diagnostic','cases','declared_case_ids',
                              'model','software','orca','cpu_python','gpu_python','cpu_executable',
                              'gpu_executable','resources')},
       'protocol_id':PROTOCOL,'settings':SETTINGS,'angular_manifest':record(angular_manifest),
       'agreement':record(agreement),'residual_inventory':record(out/'RESIDUALS.json'),
       'implementation':pins,'tasks':tasks,'new_optimizer_starts':8,'new_DFT_calls':0,
       'new_GFN2_calls_in_this_adapter':0,'source_states_changed':False,'reference':None,
       'production_changed':False,'initial_delivery':'prepared_only_no_scientific_execution'}
    write_new(out/'manifest.json',m);r=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',r)
    return r


def validate(manifest):
    m=read_json(manifest);parent=read_json(verify(m['angular_manifest']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS:
        raise InvalidArtifact('joint metal/angular method differs')
    if (m['declared_case_ids']!=list(angular.PILOT4) or len(m['tasks'])!=8 or
            [(t['case_id'],t['metal']) for t in m['tasks']]!=[(t['case_id'],t['metal']) for t in parent['tasks']]):
        raise InvalidArtifact('fixed four-case eight-endpoint scope differs')
    verify(m['agreement']);residual=read_json(verify(m['residual_inventory']))
    if residual['source_manifest']!=m['angular_manifest'] or len(residual['rows'])!=8:
        raise InvalidArtifact('residual source changed')
    for r in residual['rows']:
        for k in ('proposal_receipt','forces','coordinate'):verify(r[k])
    for k in ('source_manifest','diagnostic','cases','model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable'):
        if m[k]!=parent[k]:raise InvalidArtifact('source configuration changed: '+k)
    for k in ('software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    for t,old in zip(m['tasks'],parent['tasks']):
        for k in ('xyz','mapping','source_preparation','charge','multiplicity','q0','selector','mode_count'):
            if t[k]!=old[k]:raise InvalidArtifact('source origin/selector changed: '+k)
        if (t['active_mode_ids']!=METAL_IDS+old['active_mode_ids'] or
                t['angular_parent_mode_ids']!=old['active_mode_ids'] or
                t['active_coordinate_units']!=UNITS or t['objective_gradient_units']!=GRADIENT_UNITS):
            raise InvalidArtifact('common coordinates or explicit units differ')
        kin=Kinematics(read_json(verify(t['mapping']))['context'])
        modes=[kin.modes[i] for i in t['active_indices']]
        if [x['id'] for x in modes]!=t['active_mode_ids'] or [x['unit'] for x in modes]!=UNITS:
            raise InvalidArtifact('actual kinematic IDs/units differ')
        for i,mode in enumerate(modes[:3]):
            if (mode['kind']!='metal_translation' or mode['moving_indices']!=[kin.metal] or
                    not np.array_equal(mode['axis'],np.eye(3)[i])):
                raise InvalidArtifact('unsupported physical metal block')
        final_geometry(kin,t,np.zeros(7),[a[0] for a in xyz(verify(t['xyz']))])
    for cid in m['declared_case_ids']:
        ca,la=[next(t for t in m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
        if ca['active_mode_ids']!=la['active_mode_ids'] or ca['active_indices']!=la['active_indices']:
            raise InvalidArtifact('metal pair active spaces differ')
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if read_json(verify(la['mapping']))!=read_json(verify(ca['mapping'])):
            raise InvalidArtifact('metal pair physical mappings differ')
    return {'status':'prepared_dry_run_pass','manifest':record(manifest),'cases':4,'tasks':8,
            'active_coordinate_units':UNITS,'molecular_calls':0,'optimizer_executed':False}


class MixedProposal(MACEProposal):
    """Existing native worker/cache/geometry contract with correctly named mixed units."""
    def evaluate(self,active,purpose):
        active=np.asarray(active,dtype=float);q=full_q(self.task,active)
        if np.max(np.abs(active))>.8+1e-12:raise InvalidArtifact('trial component bounds exceeded')
        _,coords,_,_=self.kin.evaluate(q)
        key=cache_key({'manifest':self.manifest,'task':self.task['task_id'],'q':active.tolist(),'units':UNITS})
        directory=self.directory/'evaluations'/key;path=directory/'result.json'
        if path.exists():
            r=read_json(path)
            if r['cache_key']!=key or r['manifest']!=self.manifest:raise InvalidArtifact('MACE cache differs')
            self.requests.append({'purpose':purpose,'result':record(path),'reused':True})
            if r['status']!='complete':raise InvalidArtifact('previous failed evaluation retained')
            return r
        directory.mkdir(parents=True,exist_ok=False)
        r={'cache_key':key,'manifest':self.manifest,'task_id':self.task['task_id'],
           'active_q':active.tolist(),'active_coordinate_units':UNITS,'full_q':q.tolist(),
           'full_q_units':[v['unit'] for v in self.kin.modes],'status':'failed','MACE':None}
        try:
            check=geometry_check(self.kin,q,[a[0] for a in self.atoms]);r['geometry_checks']=check
            if not check['pass']:raise InvalidArtifact('unsupported physical proposal geometry')
            xp=directory/'context.xyz'
            if np.array_equal(active,np.zeros(7)):shutil.copyfile(verify(self.task['xyz']),xp)
            else:write_xyz(xp,[(a[0],*p) for a,p in zip(self.atoms,coords)])
            request={'manifest':self.manifest,'task_id':self.task['task_id'],'xyz':record(xp),
                     'charge':self.task['charge'],'multiplicity':self.task['multiplicity']}
            rp=directory/'mace_request.json';write_new(rp,request)
            native,pin=self.gpu.evaluate(record(rp));r['MACE']=pin
            if native['status']!='complete':raise InvalidArtifact('native MACE unavailable: '+native.get('reason','unknown'))
            forces=np.load(verify(native['forces']),allow_pickle=False)
            g=projected_gradient_eV(self.kin,self.task,active,forces)
            r.update(status='complete',coordinate=record(xp),MACE_eV=native['energy_eV'],forces=native['forces'],
                     gradient_eV_per_coordinate=g.tolist(),gradient_units=GRADIENT_UNITS)
        except Exception as exc:r.update(reason=str(exc),traceback=traceback.format_exc())
        write_new(path,r);self.requests.append({'purpose':purpose,'result':record(path),'reused':False})
        if r['status']!='complete':raise InvalidArtifact(r['reason'])
        return r


def optimize(manifest,task,gpu):
    from scipy.optimize import minimize
    ev=MixedProposal(manifest,task,gpu);dest=ev.directory/'result.json'
    if dest.exists():
        if read_json(dest)['manifest']!=record(manifest):raise InvalidArtifact('prior candidate manifest differs')
        return record(dest)
    began=time.monotonic();trace=[];accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],'manifest':record(manifest),
            'status':'unavailable','origin':None,'proposal':None,'optimizer':None,
            'unconstrained_minimum_claimed':False,'composite_gradient':None,'active_coordinate_units':UNITS}
    def evaluate(q,purpose):
        trial={'active_q':np.asarray(q).tolist(),'active_coordinate_units':UNITS,'purpose':purpose,
               **physical_status(ev.kin,task,q),'MACE_complete':False};trace.append(trial)
        point=ev.evaluate(q,purpose);trial['MACE_complete']=True;return point
    try:
        origin=evaluate(np.zeros(7),'origin');result['origin']=origin
        replay=(origin['MACE_eV']-task['q0']['components']['MACE_eV'])*EV_TO_KCAL
        result['origin_replay_kcal_mol']=replay
        if abs(replay)>SETTINGS['origin_replay_tolerance_kcal_mol']:raise InvalidArtifact('native origin replay failed')
        def objective(q):
            point=evaluate(q,'SLSQP_trial')
            return point['MACE_eV']-origin['MACE_eV'],np.asarray(point['gradient_eV_per_coordinate'])
        opt=minimize(objective,np.zeros(7),jac=True,method='SLSQP',bounds=[(-.8,.8)]*7,
                     constraints=[{'type':'ineq','fun':lambda q:constraints(ev.kin,task,q)[0],
                                   'jac':lambda q:constraints(ev.kin,task,q)[1]}],
                     callback=lambda q:accepted.append(np.asarray(q).tolist()),
                     options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['ftol_eV']})
        result['optimizer']={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),
                             'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev)}
        check=final_geometry(ev.kin,task,opt.x,[a[0] for a in ev.atoms]);result['final_geometry']=check
        candidate=evaluate(opt.x,'final_candidate');result['proposal']=candidate
        if not opt.success:raise InvalidArtifact('SLSQP did not converge: '+str(opt.message))
        work=candidate['MACE_eV']-origin['MACE_eV']
        if work>SETTINGS['maximum_final_energy_increase_eV']:raise InvalidArtifact('native final energy increased')
        forces=np.load(verify(candidate['forces']),allow_pickle=False)
        _,_,raw,normed,lengths=project(ev.kin.data,candidate['full_q'],forces)
        selected=set(task['active_indices']);angles=[i for i,v in enumerate(ev.kin.modes) if v['unit']=='radian']
        result.update(status='proposal_available',MACE_proposal_work_kcal_mol=float(work*EV_TO_KCAL),
                      boundary_flag=bool(check['boundary_heavy_atoms'] or check['angular_boundary_flag']),
                      residual={'force_hamiltonian':'native_vacuum_OMOL','all_mode_ids':[v['id'] for v in ev.kin.modes],
                                'all_mode_units':[v['unit'] for v in ev.kin.modes],
                                'raw_gradient_kcal_per_unit':raw.tolist(),'normalized_gradient_kcal_mol_A':normed.tolist(),
                                'metal_gradient_norm_kcal_mol_A':float(np.linalg.norm(raw[task['active_indices'][:3]])),
                                'selected_normalized_max_load':float(max(abs(normed[i]) for i in selected)),
                                'omitted_angular_max_load':float(max(abs(normed[i]) for i in angles if i not in selected))})
    except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(trials=trace,accepted_iterations=accepted,requests=ev.requests,
                  infeasible_trial_requests=sum(not t['physical_feasible'] for t in trace),
                  infeasible_completed_MACE_requests=sum(not t['physical_feasible'] and t['MACE_complete'] for t in trace),
                  maximum_trial_heavy_extent_A=max((t['maximum_heavy_displacement_A'] for t in trace),default=None),
                  wall_seconds=time.monotonic()-began,job_id=os.environ.get('SLURM_JOB_ID'))
    write_new(dest,result);return record(dest)


def execute(manifest,pool_collection):
    validate(manifest);pool=angular.pool_gate(manifest,pool_collection)
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
            int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):
        raise InvalidArtifact('existing32CPU/200000MiB/GPU allocation required')
    m=read_json(manifest);root=Path(manifest).parent;lock=(root/'proposal.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);began=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for t in m['tasks']:
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
                   'new_DFT_calls':0,'new_GFN2_calls':0})
    if error:raise InvalidArtifact(error)
    return {'status':'proposal_phase_complete','endpoint_denominator':8,'results':results,'scores_unavailable':True}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;rows=[]
    for t in m['tasks']:
        path=root/'proposals'/t['task_id']/'result.json';r=read_json(path) if path.exists() else None
        if r and r['manifest']!=record(manifest):raise InvalidArtifact('candidate manifest differs')
        ok=bool(r and r['status']=='proposal_available')
        if ok:
            kin=Kinematics(read_json(verify(t['mapping']))['context']);point=r['proposal'];q=np.asarray(point['full_q'])
            if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):
                raise InvalidArtifact('unselected mode moved')
            atoms=xyz(verify(point['coordinate']));final_geometry(kin,t,q[t['active_indices']],[a[0] for a in atoms])
            if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):
                raise InvalidArtifact('candidate XYZ differs from mapping')
        rows.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
                     'status':'candidate_available' if ok else 'unavailable','reason':r.get('reason') if r else 'not_run',
                     'proposal_receipt':record(path) if r else None,'source_preparation':t['source_preparation'],
                     'mapping':t['mapping'],'charge':t['charge'],'multiplicity':t['multiplicity'],
                     'common_active_mode_ids':t['active_mode_ids'],'active_coordinate_units':UNITS,
                     'origin_xyz':t['xyz'],'origin_components':t['q0']['components'],'candidate':r['proposal'] if ok else None,
                     'candidate_GFN2_vacuum_hartree':None,'candidate_GFN2_ALPB_hartree':None,
                     'composite_candidate_energy':None,'residual':r.get('residual') if r else None,
                     'boundary_flag':r.get('boundary_flag') if r else None})
    cases=[{'case_id':cid,'status':'ready_for_common_pool' if all(r['status']=='candidate_available' for r in rows if r['case_id']==cid) else 'unavailable',
            'selected_geometry':None,'score':None} for cid in m['declared_case_ids']]
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'source_manifest':m['source_manifest'],
            'cases':cases,'endpoints':rows,'case_denominator':4,'endpoint_denominator':8,
            'available_candidates':sum(r['status']=='candidate_available' for r in rows),
            'new_DFT_calls':0,'new_GFN2_calls_in_adapter':0,'common_pool_scoring_performed':False,
            'production_changed':False,'reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('cases','endpoints')}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for k in ('angular_manifest','agreement','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q=s.add_parser('dry-run');q.add_argument('--manifest',required=True)
    q=s.add_parser('execute');q.add_argument('--manifest',required=True);q.add_argument('--pool-collection',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation');fn={'prepare':prepare,'dry-run':validate,'execute':execute,'collect':collect}[op]
    print(json.dumps(fn(**a),indent=2))


if __name__=='__main__':main()
