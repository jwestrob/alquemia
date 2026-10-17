#!/usr/bin/env python3
"""Finite matched fixed-FF LBT continuation; numerical work only inside Slurm.
--stage preflight|scout|remaining --config CONFIG --output FRESH
Uses the frozen ten-mode mechanics/search and the qualified native FF solvent.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import time
import traceback
from types import SimpleNamespace

PROTOCOL = 'lbt3_matched_fixedFF_response_v1'
STATES = (('7CCO','La'),('7CCO','Dy'),('7CCN','La'),('7CCN','Dy'))
START_IDS = dict(zip(STATES, ('geometry_11','geometry_13','geometry_15','geometry_17')))
BASIS_SHA = 'b8cd5cb3015584eceada929662f7945fc7ec878607cf5ab2a7c0c494ba1f07d5'
MATCHED_SHA = 'c9856d48a0fd32dde3df620d097eb304eb9ba81c0c4e9afd76465abb44f79452'
EV, HA, BOHR = 23.06054783061903, 627.509474, .529177210903
COST_KEYS = ('new_composite_evaluations','MACE_calls_attempted','ALPB_calls_attempted',
             'MACE_call_seconds','ALPB_call_seconds','model_load_seconds','cache_hits')


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p=Path(record['path']).resolve(strict=True)
    require(pin(p)['sha256']==record['sha256'],'Changed pin: '+str(p)); return p


def read(record): return json.loads(verify(record).read_text())


def write(path,value):
    with Path(path).open('x') as f: json.dump(value,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n')


def module(name,record):
    spec=importlib.util.spec_from_file_location(name,verify(record)); obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj); return obj


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def recipe(world,source,metal):
    c=world['base_cells'][(source['id'],metal)]
    return {'protocol_id':PROTOCOL,'physical_charge':-2,'physical_multiplicity':1 if metal=='La' else 6,
        'target_index':255,'physical_ids':[a['id'] for a in source['atoms']], 'symbols':c['symbols'],
        'model':world['mcfg']['model'],'MACE_versions':world['mcfg']['expected_versions'],
        'native_ALPB_binary':world['binary'],'native_ALPB_build':world['build'],
        'fixed_FF_charges_e':c['q'],'field_and_MM_ledger':'absent by whole-peptide partition',
        'solvent_definition':'native ALPB water fixed ff19SB plus formal3 target; no SCC'}


def key_for(world,source,metal,x):
    return digest({'recipe':recipe(world,source,metal),'positions_A':np.asarray(x,dtype=float).tolist()})


def relative(a,b):
    mace=(a['MACE_eV']-b['MACE_eV'])*EV
    solvent=(a['fixed_ALPB_hartree']-b['fixed_ALPB_hartree'])*HA
    return {'MACE_kcal_mol':mace,'fixed_ALPB_kcal_mol':solvent,'composite_kcal_mol':math.fsum((mace,solvent))}


def raw_components(world,source,metal,point,need_gradient):
    """Admit only the pinned MACE and native fixed-q solvent, never old composites."""
    h=world['h']; raw=read(point['MACE']); sol=read(point['ALPB']); h.pins_in(sol['artifacts'])
    h.state(raw['input_state_check'],256,-2,1 if metal=='La' else 6,255)
    require(math.isfinite(raw['energy_eV']),'Nonfinite MACE energy')
    start=read(sol['artifacts']['START.json']); c=world['base_cells'][(source['id'],metal)]
    sy,x=h.xyz(point['xyz']); ssy,sx=h.xyz(start['xyz'])
    require(sy==ssy==c['symbols'] and x==sx and start['physical_source']==point['physical_source'],
            'Fixed-q solvent physical input differs')
    require([float(v) for v in verify(start['charges']).read_text().split()]==c['q']
            and start['command']==[str(verify(world['binary'])),start['xyz']['path'],start['charges']['path']],
            'Native solvent command/charges differ')
    execution=read(sol['artifacts']['EXECUTION.json'])
    require(execution['returncode']==0 and execution['command']==start['command']
            and abs(sol['charge_e']+2)<=1e-5 and sol['rounding_bounds_hartree']==[0.,0.]
            and set(sol['components_hartree'])==set(h.PARTS)
            and sol['total_hartree']==math.fsum(sol['components_hartree'].values())
            and sol['gradient_units']=='Hartree/Bohr; archived only','Native solvent admission differs')
    require(all(math.isfinite(v) for v in sol['components_hartree'].values()),'Nonfinite ALPB energy')
    native={'MACE_eV':raw['energy_eV'],'field_Eh':0.,'MM_kcal_mol':0.,
            'ALPB_Eh':sol['total_hartree'],'ALPB_parts_Eh':sol['components_hartree']}
    gradient=None
    if need_gradient:
        require(raw['force_convention']=='negative derivative of the same native vacuum energy','MACE force convention differs')
        forces=np.asarray(h.vectors(raw['forces_eV_A'],256),dtype=float)
        ag=np.asarray(h.vectors([[float(v) for v in line.split()] for line in verify(sol['artifacts']['gradient.dat']).read_text().splitlines()],256),dtype=float)
        gradient=-forces*EV+ag*(HA/BOHR)
        require(gradient.shape==(256,3) and np.isfinite(gradient).all(),'Invalid matched gradient')
    return native,gradient


def admit_point(world,source,metal,record,require_gradient=False):
    point=read(record)
    require(point['status']=='complete' and point['source_id']==source['id'] and point['metal']==metal
            and point['recipe']==recipe(world,source,metal),'Point state/recipe differs')
    x=world['pool'].check_coordinates(world['kinworld'],source,metal,point['xyz'],point['physical_source'],point['positions'])
    require(point['cache_key']==key_for(world,source,metal,x),'Point cache identity differs')
    native,gradient=raw_components(world,source,metal,point,require_gradient)
    require(point['native_energies']==native and point['components']=={'MACE_eV':native['MACE_eV'],'fixed_ALPB_hartree':native['ALPB_Eh']},'Point raw energy replay differs')
    if require_gradient:
        saved=np.load(verify(point['gradient_kcal_mol_A']),allow_pickle=False)
        require(saved.shape==(256,3) and np.isfinite(saved).all() and np.array_equal(saved,gradient),'Point matched raw-gradient replay differs')
    return point,x,native


def archive_saved(world,source,metal,c,row,d):
    """Reconstruct saved matched energy/forces with zero molecular calls."""
    d.mkdir(); ev=row['MACE_evidence']; raw=read(ev['raw_MACE'])
    if ev['kind']=='origin':
        raw=next(r for r in raw['rows'] if r['task_id']==c['task']['task_id'])
    write(d/'MACE.json',raw); np.save(d/'POSITIONS.npy',np.asarray(c['x'],dtype=float),allow_pickle=False)
    p={'status':'complete','source_id':source['id'],'metal':metal,'recipe':recipe(world,source,metal),
       'cache_key':key_for(world,source,metal,c['x']),'positions':pin(d/'POSITIONS.npy'),
       'physical_source':ev['physical_source'],'xyz':ev['xyz'],'MACE':pin(d/'MACE.json'),'ALPB':row['ALPB'],
       'matched_parent':world['config']['matched_pool_result'],'matched_candidate_id':row['candidate_id'],
       'archived_MACE_evidence':ev,'new_MACE_calls':0,'new_ALPB_calls':0}
    native,gradient=raw_components(world,source,metal,p,True)
    require(native==row['native_energies'],'Saved matched native components differ')
    p.update(native_energies=native,components={'MACE_eV':native['MACE_eV'],'fixed_ALPB_hartree':native['ALPB_Eh']},
             gradient_convention='positive derivative kcal/mol/A; -MACE force + fixed-q ALPB gradient')
    np.save(d/'GRADIENT.npy',gradient,allow_pickle=False); p['gradient_kcal_mol_A']=pin(d/'GRADIENT.npy')
    write(d/'POINT.json',p); record=pin(d/'POINT.json'); admit_point(world,source,metal,record,True)
    return p,record,gradient


def load_inputs(cfg,out):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Response adapter differs'); verify(cfg['plan'])
    matched=read(cfg['matched_pool_result']); mpcfg=read(matched['config'])
    verify(matched['implementation']);verify(mpcfg['implementation'])
    require(matched['status']=='complete' and matched['slurm_job_id']=='12720'
            and matched['implementation']['sha256']==mpcfg['implementation']['sha256']==MATCHED_SHA
            and len(matched['rows'])==36 and all(r['status']=='complete' for r in matched['rows']), 'Complete frozen matched pool required')
    mp=module('matched_response_pool_parent',matched['implementation']); d=out/'matched_parent_admission'; d.mkdir()
    h,cells,_,scfg,build,binary,old=mp.inputs(mpcfg,d)
    basis_parent=read(old['basis_scout_result']); bcfg=read(basis_parent['config'])
    require(basis_parent['implementation']['sha256']==BASIS_SHA and basis_parent['implementation']==bcfg['implementation'],'Frozen basis helper differs')
    basis=module('matched_response_basis',basis_parent['implementation'])
    (_,range_module,pool,base,kinworld,settings,old_starts,_,_,_,lbfgs,old_cfg)=basis.inputs(verify(basis_parent['config']))
    global np,FALLBACK
    import numpy as np
    FALLBACK=basis.FALLBACK
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
    for k in ('LD_PRELOAD','PYTHONHOME','PYTHONPATH'):env.pop(k,None)
    for lib in build['linked_libraries']:verify(lib)
    by={(r['candidate_id'],r['metal']):r for r in matched['rows']}
    cb={(r['candidate_id'],r['metal']):c for c,r in cells}
    world={'config':cfg,'matched_pool':matched,'old_pool':old,'h':h,'binary':binary,'build':build,'env':env,
        'basis':basis,'base':base,'range':range_module,'pool':pool,'kinworld':kinworld,'sources':kinworld['sources'],
        'settings':settings,'mcfg':kinworld['mcfg'],'native_helper':kinworld['fd'],'lbfgs':lbfgs,'starts':{},
        'base_cells':{(sid,m):cb[('geometry_00' if sid=='7CCO' else 'geometry_01',m)] for sid,m in STATES},
        'matched_cells':cb,'basis_implementation':basis_parent['implementation'],'matched_implementation':matched['implementation']}
    checks=[]; saved=out/'matched_saved_starts';saved.mkdir()
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS
    for sid,metal in STATES:
        source=world['sources'][sid];old_start=old_starts[(sid,metal,'candidate')];cid=START_IDS[(sid,metal)]
        slot=next(r for r in old['searches'] if r['source_id']==sid and r['metal']==metal and r['arm']=='candidate')
        search_old=read(slot['search']); row=by[(cid,metal)]; c=cb[(cid,metal)]
        require(slot['status']=='admitted' and slot['candidate_id']==cid and slot['point']==search_old['candidate']
                and row['MACE_evidence']['evidence_record']==slot['point'] and search_old['candidate_admitted']
                and search_old['final_geometry']['pass'] and search_old['source_rebased'] is False
                and search_old['active_indices']==old_start['active'] and search_old['coordinate_scaling_T']==old_start['metric']['coordinate_scaling_T'],
                'Exact ten-mode terminal lineage differs')
        q=np.asarray(search_old['full_q'],dtype=float);x=np.asarray(c['x'],dtype=float)
        require(q.shape==(64,) and np.isfinite(q).all(),'Invalid old terminal coordinates')
        active=old_start['active'];require(len(active)==10 and np.count_nonzero(q[[i for i in range(64) if i not in active]])==0,'Mode inventory differs')
        replay,_,jac,_=source['kin'].evaluate(q)
        require(np.allclose(replay,x,rtol=0,atol=1e-12),'Terminal source-map replay differs')
        guard=range_module.range_guard(base,source,active,Geometry,GEOMETRY_SETTINGS,metal)
        origin_guard=guard.check(np.zeros(64),source['x'],source['x']);start_guard=guard.check(q,x,x)
        require(origin_guard['pass'] and start_guard['pass'],'Original/warm geometry admission failed')
        origin_id='geometry_00' if sid=='7CCO' else 'geometry_01'
        origin,op,_=archive_saved(world,source,metal,cb[(origin_id,metal)],by[(origin_id,metal)],saved/(sid+'_'+metal+'_origin'))
        point,pp,gradient=archive_saved(world,source,metal,c,row,saved/(sid+'_'+metal+'_start'))
        start=dict(old_start,search=search_old,search_pin=slot['search'],point=point,point_pin=pp,
            origin=origin,origin_pin=op,x=x,q=q,jac=jac,gradient=gradient,guard=guard,origin_guard=origin_guard,
            start_guard=start_guard,cache={origin['cache_key']:op,point['cache_key']:pp},records=[op,pp],
            old_objective_terminal=slot['point'])
        world['starts'][(sid,metal)]=start
        checks.append({'source_id':sid,'metal':metal,'pass':True,'candidate_id':cid,'previous_search':slot['search'],
            'old_objective_terminal':slot['point'],'matched_start':pp,'matched_origin':op,'initialization_molecular_calls':0,
            'old_termination':search_old['termination'],'active_indices':active,'active_names':search_old['active_names'],
            'coordinate_scaling_T':start['transform'].tolist(),'origin_guard':origin_guard,'start_guard':start_guard,
            'matched_work_from_origin':relative(point['components'],origin['components']),
            'matched_force_diagnostic':basis.projected_diagnostic(jac,gradient,active,source['kin'].modes),
            'old_composite_gradient_reused':False,'fresh_history':True})
    world['preflight']={'status':'passed','declared_states':4,'molecular_calls':0,'model_loaded':False,'rows':checks}
    return world


class MatchedEvaluator:
    def __init__(self,world,source,metal,out,start=None,native=None):
        self.world,self.source,self.metal,self.out=world,source,metal,Path(out)
        self.native=native;self.cache=dict(start['cache']) if start else {};self.records=list(start['records']) if start else []
        self.requests=[];self.cost={k:0 for k in COST_KEYS}

    def key(self,positions):return key_for(self.world,self.source,self.metal,positions)

    def evaluate(self,positions,compute_force=True):
        x=np.asarray(positions,dtype=float);require(x.shape==(256,3) and np.isfinite(x).all(),'Invalid whole LBT positions')
        key=self.key(x)
        if key in self.cache:
            require(self.cache[key] is not None,'Previously failed exact state; no retry')
            p,_,_=admit_point(self.world,self.source,self.metal,self.cache[key],compute_force)
            self.cost['cache_hits']+=1;return p,self.cache[key],True
        self.cost['new_composite_evaluations']+=1
        d=self.out/('evaluation_%03d'%self.cost['new_composite_evaluations']);d.mkdir()
        physical=[dict(a,xyz_A=v.tolist()) for a,v in zip(self.source['atoms'],x)]
        write(d/'PHYSICAL.json',physical);np.save(d/'POSITIONS.npy',x,allow_pickle=False)
        c0=self.world['base_cells'][(self.source['id'],self.metal)];sy=c0['symbols']
        (d/'whole.xyz').write_text('256\nMatched fixed-FF LBT response\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(sy,x)))
        task=dict(c0['task'],task_id=self.source['id']+'_'+self.metal+'_'+d.name,
            physical_source=pin(d/'PHYSICAL.json'),xyz=pin(d/'whole.xyz'),origin_task_id=c0['task']['task_id'])
        write(d/'TASK.json',task)
        self.world['pool'].check_coordinates(self.world['kinworld'],self.source,self.metal,task['xyz'],task['physical_source'],pin(d/'POSITIONS.npy'))
        request={'status':'unavailable','source_id':self.source['id'],'metal':self.metal,'cache_key':key,
            'recipe':recipe(self.world,self.source,self.metal),'positions':pin(d/'POSITIONS.npy'),
            'physical_source':task['physical_source'],'xyz':task['xyz'],'task':pin(d/'TASK.json'),
            'compute_force':compute_force,'new_MACE_calls':0,'new_ALPB_calls':0}
        self.cache[key]=None;tick=time.monotonic()
        try:
            if self.native is None:
                t=time.monotonic();self.native=self.world['native_helper'].NativeMACE(self.world['mcfg']);self.cost['model_load_seconds']+=time.monotonic()-t
            before=self.native.calls;t=time.monotonic()
            try:raw=self.native.evaluate(task,sy,x,compute_force=compute_force)
            finally:
                request['new_MACE_calls']=self.native.calls-before;self.cost['MACE_calls_attempted']+=request['new_MACE_calls'];self.cost['MACE_call_seconds']+=time.monotonic()-t
            write(d/'MACE.json',raw)
            for lib in self.world['build']['linked_libraries']:verify(lib)
            c=dict(c0,task=task,x=x.tolist(),physical=physical)
            t=time.monotonic();self.cost['ALPB_calls_attempted']+=1;request['new_ALPB_calls']=1
            try:self.world['h'].solvent(self.world['binary'],c,d/'ALPB',self.world['env'])
            finally:self.cost['ALPB_call_seconds']+=time.monotonic()-t
            point=dict(request,status='complete',MACE=pin(d/'MACE.json'),ALPB=pin(d/'ALPB'/'RESULT.json'))
            native,gradient=raw_components(self.world,self.source,self.metal,point,compute_force)
            point.update(native_energies=native,components={'MACE_eV':native['MACE_eV'],'fixed_ALPB_hartree':native['ALPB_Eh']})
            if compute_force:
                np.save(d/'GRADIENT.npy',gradient,allow_pickle=False)
                point.update(gradient_kcal_mol_A=pin(d/'GRADIENT.npy'),gradient_convention='positive derivative kcal/mol/A; -MACE force + fixed-q ALPB gradient')
            point['elapsed_seconds']=time.monotonic()-tick;write(d/'POINT.json',point);record=pin(d/'POINT.json')
            admit_point(self.world,self.source,self.metal,record,compute_force)
            self.cache[key]=record;self.records.append(record);request.update(status='complete',point=record)
            return point,record,False
        except Exception as exc:
            request.update(failure=str(exc),traceback=traceback.format_exc());raise
        finally:
            request['elapsed_seconds']=time.monotonic()-tick;write(d/'REQUEST.json',request)
            self.requests.append(dict(request,record=pin(d/'REQUEST.json')))


# Search copied from frozen b8cd5cb3; only objective adapter, receipt fields,
# three-value evaluator return and final warm-start energy gate differ.
def projected_diagnostic(jac,gradient,active,modes):
    return BASIS.projected_diagnostic(jac,gradient,active,modes)


def feasible_coordinates(start,q,positions,jac,gradient):
    return BASIS.feasible_coordinates(start,q,positions,jac,gradient)


def search(base, start, arm, evaluator, out, lbfgs_direction):
    source, active, transform = start['source'], start['active'], start['transform']
    kin, guard = source['kin'], start['guard']; n = len(active)
    q0 = start['q'].copy(); q = q0.copy(); v = np.zeros(n)
    positions, jac, gradient = (start[key].copy() for key in ('x', 'jac', 'gradient'))
    point, point_pin = start['point'], start['point_pin']; origin = start['origin']; history = []
    gv = transform.T@np.einsum('kij,ij->k', jac[active], gradient)
    current_guard = start['start_guard']
    result = {'protocol_id': PROTOCOL, 'source_id': source['id'], 'metal': base.METAL, 'arm': arm,
        'status': 'unavailable', 'candidate_admitted': False, 'active_indices': active,
        'active_names': [kin.modes[k]['name'] for k in active], **start['metric'],
        'coordinate_offset_q64': q0.tolist(), 'coordinate_rule': 'q_active=q_start_active+T_arm@v; v initially zero',
        'previous_coordinate_scaling_T': start['search']['coordinate_scaling_T'],
        'authoritative_origin': start['search']['authoritative_origin'], 'original_point': start['origin_pin'],
        'previous_search': start['search_pin'], 'continuation_start': start['point_pin'],
        'old_objective_terminal': start['old_objective_terminal'],
        'descent_reference': 'exact matched-objective warm start; original work is reported without a sign gate',
        'origin_guard': start['origin_guard'], 'continuation_start_guard': current_guard,
        'initial_history_length': 0, 'initial_new_evaluations': 0, 'initial_offset_v_A': v.tolist(),
        'source_rebased': False, 'trace': [], 'accepted_steps': 0, 'accepted_fallback_steps': 0,
        'fallback_rule': FALLBACK, 'stationary_minimum_claimed': False}
    try:
        while True:
            diagnostic = projected_diagnostic(jac, gradient, active, kin.modes)
            row = {'accepted_step': result['accepted_steps'], 'offset_v_A': v.tolist(), 'full_q': q.tolist(), 'point': point_pin,
                'work_from_origin': base.relative(point['components'], origin['components']),
                'work_from_continuation_start': base.relative(point['components'], start['point']['components']),
                'geometry': current_guard, 'force_diagnostic': diagnostic}
            result['trace'].append(row)
            with (out/'TRACE.jsonl').open('a') as handle:
                handle.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
            if diagnostic['low_active_span_force']:
                termination = 'low_active_span_force'; break
            if result['accepted_steps'] >= 12:
                termination = 'accepted_step_limit'; break
            if evaluator.cost['new_composite_evaluations'] >= 24:
                termination = 'evaluation_limit'; break
            direction = lbfgs_direction(gv, history)
            if float(gv@direction) >= -1e-14:
                history = []; direction = -gv
            tangent = np.einsum('k,kij->ij', transform@direction, jac[active])
            speed = float(np.linalg.norm(tangent, axis=1).max())
            if speed <= 1e-14:
                termination = 'zero_active_direction'; break
            direction *= min(1., .05/speed); slope = float(gv@direction)
            require(slope < 0, 'Own-energy direction is not descent')
            row['trials'] = []; accepted = None; used_fallback = False
            for multiplier in base.SETTINGS['trial_multipliers']:
                trial_v = v+multiplier*direction; trial_q = q0.copy(); trial_q[active] = q0[active]+transform@trial_v
                trial_positions, _, trial_jac, _ = kin.evaluate(trial_q)
                check = guard.check(trial_q, trial_positions, positions)
                trial = {'multiplier': multiplier, 'offset_v_A': trial_v.tolist(), 'full_q': trial_q.tolist(),
                    'geometry': check, 'molecular_evaluation_requested': False, 'accepted': False}
                row['trials'].append(trial)
                if not check['pass']:
                    trial['reason'] = 'geometry_or_bound'; continue
                if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(trial_positions) not in evaluator.cache:
                    trial['reason'] = 'evaluation_limit'; break
                trial['molecular_evaluation_requested'] = True  # cache-energy trials also prohibit fallback entry
                next_point, next_pin, _ = evaluator.evaluate(trial_positions)
                change = base.relative(next_point['components'], point['components'])
                trial.update(point=next_pin, work_from_current=change, armijo_rhs_kcal_mol=1e-4*multiplier*slope)
                if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                    trial['accepted'] = True
                    accepted = (trial_v, trial_q, trial_positions, trial_jac, next_point, next_pin, check); break
                trial['reason'] = 'own_composite_Armijo'
            fallback_allowed = (accepted is None and len(row['trials']) == 4
                and all(t.get('reason') == 'geometry_or_bound' and not t['molecular_evaluation_requested'] for t in row['trials']))
            row['fallback'] = {'entered': fallback_allowed, 'evaluated_candidates': []}
            if fallback_allowed:
                candidates, proposals = feasible_coordinates(start, q, positions, jac, gradient)
                row['fallback'].update(coordinate_proposals=proposals,
                    ranked_mode_indices=[c['record']['mode_index'] for c in candidates])
                for candidate in candidates[:2]:
                    if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(candidate['positions']) not in evaluator.cache:
                        row['fallback']['stop'] = 'evaluation_limit'; break
                    trial = {'mode_index': candidate['record']['mode_index'], 'accepted': False,
                        'predicted_own_decrease_kcal_mol': -candidate['slope'],
                        'armijo_rhs_kcal_mol': 1e-4*candidate['slope'], 'molecular_evaluation_requested': True}
                    row['fallback']['evaluated_candidates'].append(trial)
                    next_point, next_pin, _ = evaluator.evaluate(candidate['positions'])
                    change = base.relative(next_point['components'], point['components'])
                    trial.update(point=next_pin, work_from_current=change)
                    if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                        trial['accepted'] = True; used_fallback = True
                        accepted = (candidate['v'], candidate['q'], candidate['positions'], candidate['jac'],
                                    next_point, next_pin, candidate['guard']); break
                    trial['reason'] = 'own_composite_Armijo'
            write(out/('ITERATION_%02d.json' % result['accepted_steps']), row)
            if accepted is None:
                termination = 'evaluation_limit' if evaluator.cost['new_composite_evaluations'] >= 24 else 'no_admissible_decreasing_trial'; break
            new_v, q, positions, jac, point, point_pin, current_guard = accepted
            gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
            new_gv = transform.T@np.einsum('kij,ij->k', jac[active], gradient)
            if used_fallback:
                history = []; result['accepted_fallback_steps'] += 1
            else:
                displacement, difference = new_v-v, new_gv-gv; sy = float(displacement@difference)
                if sy > 1e-10*np.linalg.norm(displacement)*np.linalg.norm(difference) and sy > 1e-16:
                    history.append((displacement, difference, 1./sy)); history = history[-8:]
            v, gv = new_v, new_gv; result['accepted_steps'] += 1
        final_check = guard.check(q, positions, positions); require(final_check['pass'], 'Terminal physical gate failed')
        work = base.relative(point['components'], origin['components'])
        incremental = base.relative(point['components'], start['point']['components'])
        require(incremental['composite_kcal_mol'] <= 1e-8, 'Matched warm-start energy increased')
        result.update(status='finite_proposal_available', candidate_admitted=True, termination=termination, candidate=point_pin,
            full_q=q.tolist(), offset_v_A=v.tolist(), final_geometry=final_check, final_force_diagnostic=diagnostic,
            work_from_origin=work, work_from_continuation_start=incremental, preserved_boundary_flags=final_check['boundary_flags'])
    except Exception as exc:
        result.update(status='unavailable', candidate_admitted=False, termination='backend_or_numerical_failure',
            failure=str(exc), traceback=traceback.format_exc(), last_archived_valid_point=point_pin)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','scout','remaining'),required=True)
    for key in ('config','output'):parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args();require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    result={'protocol_id':PROTOCOL,'stage':args.stage,'status':'unavailable','implementation':pin(__file__),
        'config':pin(args.config),'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),
        'preflight':{'status':'not_passed'},'rows':[],'points':[],'new_initialization_molecular_calls':0,
        'new_GFN2_calls':0,'new_MM_calls':0,'new_DFT_calls':0,'source_rebased':False,
        'pool_comparison_performed':False,'stationary_minimum_claimed':False}
    write(out/'START.json',result)
    warm=None
    try:
        cfg=json.loads(args.config.read_text());world=load_inputs(cfg,out)
        global BASIS
        BASIS=world['basis']
        result.update(plan=cfg['plan'],matched_pool_result=cfg['matched_pool_result'],settings=world['settings'],
            fallback_rule=FALLBACK,preflight=world['preflight'],basis_implementation=world['basis_implementation'],
            matched_implementation=world['matched_implementation'],native_ALPB_binary=world['binary'],
            model=world['mcfg']['model'],versions=world['mcfg']['expected_versions'])
        write(out/'PREFLIGHT.json',world['preflight']);result['preflight_result']=pin(out/'PREFLIGHT.json')
        if args.stage=='preflight':result['status']='complete'
        else:
            require(int(os.environ['SLURM_CPUS_PER_TASK'])>=1 and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'), 'Allocated GPU required')
            if args.stage=='remaining':
                scout=read(cfg['scout_result']);scfg=read(scout['config'])
                verify(scout['implementation']);verify(scfg['implementation'])
                require(scout['status']=='complete' and scout['stage']=='scout' and scout['protocol_id']==PROTOCOL
                    and scout['continuation_execution_and_physical_admission'] and scout['preflight']['status']=='passed'
                    and scout['implementation']['sha256']==scfg['implementation']['sha256']==cfg['implementation']['sha256']
                    and scout['plan']==cfg['plan'] and scout['matched_pool_result']==cfg['matched_pool_result']
                    and scout['settings']==world['settings'] and scout['fallback_rule']==FALLBACK
                    and len(scout['rows'])==1,'Admitted unchanged scout required')
                sr=scout['rows'][0];sp=read(sr['search'])
                require((sr['source_id'],sr['metal'])==STATES[0] and sr['status']=='complete' and sr['candidate_admitted']
                    and read(sr['result'])=={k:v for k,v in sr.items() if k!='result'}
                    and sp['candidate_admitted'] and sp['final_geometry']['pass']
                    and sp['work_from_continuation_start']['composite_kcal_mol']<=1e-8,'Scout execution/physical gate failed')
                result['scout_result']=cfg['scout_result']
            states=STATES[:1] if args.stage=='scout' else STATES[1:]
            result['declared_tasks']=[{'source_id':s,'metal':m} for s,m in states]
            for sid,metal in states:
                d=out/(sid+'_'+metal);d.mkdir();start=world['starts'][(sid,metal)]
                evaluator=MatchedEvaluator(world,start['source'],metal,d,start,native=warm)
                prior=warm.calls if warm is not None else 0
                scope=SimpleNamespace(METAL=metal,SETTINGS=world['base'].SETTINGS,relative=relative)
                row={'source_id':sid,'metal':metal,'status':'unavailable','candidate_admitted':False,
                    'physical_multiplicity':1 if metal=='La' else 6,'previous_search':start['search_pin'],
                    'continuation_start':start['point_pin'],'old_objective_terminal':start['old_objective_terminal']}
                try:
                    proposal=search(scope,start,'matched_fixedFF',evaluator,d,world['lbfgs'])
                    write(d/'SEARCH.json',proposal)
                    row.update(status='complete' if proposal['candidate_admitted'] else 'unavailable',
                        candidate_admitted=proposal['candidate_admitted'],search=pin(d/'SEARCH.json'))
                except Exception as exc:row.update(failure=str(exc),traceback=traceback.format_exc())
                warm=evaluator.native;row.update(actual_cost=evaluator.cost,points=evaluator.records,requests=evaluator.requests)
                result['points'].extend(evaluator.records)
                if warm is not None:
                    row['native_runtime']={'versions':warm.versions,'device':warm.device,'worker_prior_forwards':prior,
                        'worker_cumulative_forwards':warm.calls,'actual_model_forwards_this_search':warm.calls-prior,
                        'torch_threads':warm.torch.get_num_threads()}
                    require(warm.calls-prior==evaluator.cost['MACE_calls_attempted'],'Native forward accounting differs')
                write(d/'RESULT.json',row);row['result']=pin(d/'RESULT.json');result['rows'].append(row)
                print(json.dumps({'source_id':sid,'metal':metal,'status':row['status']}),flush=True)
            result['status']='complete' if all(r['candidate_admitted'] for r in result['rows']) else 'unavailable'
            result['continuation_execution_and_physical_admission']=result['status']=='complete'
        verify(result['config']);verify(cfg['plan'])
    except Exception as exc:result.update(status='execution_failure',failure=str(exc),traceback=traceback.format_exc())
    result['actual_cost']={k:sum(r.get('actual_cost',{}).get(k,0) for r in result['rows']) for k in COST_KEYS}
    limit=1 if args.stage=='scout' else 3 if args.stage=='remaining' else 0
    require(all(result['actual_cost'][k]<=24*limit for k in ('new_composite_evaluations','MACE_calls_attempted','ALPB_calls_attempted')), 'Declared molecular ceiling exceeded')
    if warm is not None:
        result['native_runtime']={'actual_model_forwards':warm.calls,'model_loaded_once':True,
            'peak_cuda_allocated_bytes':warm.torch.cuda.max_memory_allocated()}
    result['elapsed_seconds']=time.monotonic()-tick
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
    result['application_allocated_GPU_seconds']=result['elapsed_seconds'] if args.stage!='preflight' else 0.
    result['peak_host_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation']={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_CPU','CUDA_VISIBLE_DEVICES')}
    write(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}),flush=True)
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
