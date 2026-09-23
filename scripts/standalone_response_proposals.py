"""Finite composite-gradient proposals on six consumed, exactly pinned PQQ sources."""
from __future__ import annotations
import argparse
import concurrent.futures as futures
import copy
import fcntl
import json
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from scipy.optimize import minimize
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,read_json,record,verify,write_new,xyz
import standalone_xtb as xtb
from accommodation_nonlinear import WarmGPU,project_components,relative_components
from accommodation_torsion_profiles import geometry_check
import adaptive_angular_proposals as angular
from adaptive_origin_recovery import snapshot,allocation
from mace_hybrid import write_xyz
from mace_site_kinematics import Kinematics
from nikasha_pool import choose_rows
from nikasha_pool_compare import call

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL='standalone_xtb671_six_source_bounded_composite_response_v1'
CASES=('q9z4j7-pqq-la_model','c5axv8-pqq-la_model','q88jh5-pqq-la_model',
'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1',
'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4',
'a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4')
BASE=('origin','adaptive_Ca','adaptive_La')
SETTINGS={'maxiter':20,'maximum_distinct_evaluations':40,'ftol_Ha_equivalent':.001/HA_TO_KCAL,
'angle_bound_radian':.8,'maximum_heavy_displacement_A':.8,'start_replay_tolerance_kcal_mol':.1,
'start_rule':'own_composite_minimum_origin_adaptive_Ca_adaptive_La_tie_order',
'proposal_rule':'lowest_energy_complete_feasible_evaluation_earliest_tie_on_success_or_declared_limit',
'failure_policy':'molecular_geometry_gradient_or_other_optimizer_failure_invalidates_endpoint_no_retry',
'objective':'own_composite_energy_shifted_from_start_in_Hartree_equivalent',
'accuracy':'0.02','maximum_new_standalone_calls':1000,'maximum_new_MACE_calls':480,
'endpoint_workers':2,'solvent_workers_per_endpoint':2,'solvent_threads':8}


def paths():
    w=ROOT/'workspaces'
    return {'canonical_native':w/'adaptive_completion_20260922/original30_pool_v1/final_collection.json',
    'fold_native':w/'adaptive_completion_20260922/primary225_pool_v1/final_collection.json',
    'canonical_proposals':w/'adaptive_completion_20260922/original30_v1/manifest.json',
    'fold_proposals':w/'adaptive_completion_20260922/primary225_v2_sharded/manifest.json',
    'external_reference':w/'standalone_xtb_reference_20260923/run_v1/REFERENCE_v1.json',
    'external_canonical':w/'standalone_xtb_reference_20260923/run_v1/collection_v1.json',
    'external_pilot':w/'standalone_xtb_20260923/run_v1/collection_v2.json',
    'external_folds':w/'standalone_xtb_transfer_20260923/run_v1/collection_v1.json'}


def low_task(cell,medium,directory,grad):
    d=Path(directory);d.mkdir(parents=True,exist_ok=False);shutil.copyfile(verify(cell['xyz']),d/'core.xyz');(d/'xcontrol').write_text(xtb.control())
    return {'task_id':d.name,'kind':'pool','cell_id':d.name,'accuracy':'0.02','case_id':cell['case_id'],
      'candidate':cell['candidate'],'metal':cell['metal'],'medium':medium,'offset_radian':None,
      'charge':cell['charge'],'multiplicity':1,'gradient_requested':grad,'xyz':record(d/'core.xyz'),
      'input':record(d/'xcontrol'),'directory':str(d.resolve())}


def external_reuse(row,cell,medium,backend):
    rec=read_json(verify(row['receipt']));mp=verify(rec['manifest']);m=read_json(mp)
    t=next(t for t in m['tasks'] if t['task_id']==rec['task_id'])
    if (row['status']!='complete' or t['accuracy']!='0.02' or t['medium']!=medium or t['charge']!=cell['charge'] or
        t['multiplicity']!=1 or xyz(verify(t['xyz']))!=xyz(verify(cell['xyz'])) or
        any(m[k]!=backend[k] for k in ('executable','installed_parameter','package')) or
        verify(m['parameter']).read_bytes()!=verify(backend['parameter']).read_bytes()):raise InvalidArtifact('standalone baseline reuse differs')
    actual=xtb.parse_task(mp,t)
    if actual['status']!='complete' or actual['energy_hartree']!=row['energy_hartree']:raise InvalidArtifact('actual standalone row no longer replays')
    return {'status':'complete','row':actual,'task':t,'manifest':record(mp)}


def prepare(agreement,output):
    src={k:record(p) for k,p in paths().items()};data={k:read_json(verify(p)) for k,p in src.items()}
    bm=read_json(verify(data['external_pilot']['manifest']));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    params=out/'parameters';params.mkdir();shutil.copyfile(verify(bm['parameter']),params/'param_gfn2-xtb.txt')
    method=read_json(verify(data['canonical_native']['manifest']))
    m={k:method[k] for k in ('model','software','gpu_python','cpu_python')}
    m.update({k:bm[k] for k in ('executable','installed_parameter','package')});m['parameter']=record(params/'param_gfn2-xtb.txt')
    external={}
    for name in ('external_canonical','external_pilot','external_folds'):
        for r in data[name]['rows']:
            if r['accuracy']=='0.02' and r['kind']=='pool':external.setdefault((r['case_id'],r['candidate'],r['metal'],r['medium']),r)
    tasks=[];cases=[];missing=[];force_count=0
    for cid in CASES:
        canonical='__conditioned_' not in cid;npool=data['canonical_native' if canonical else 'fold_native'];proposals=data['canonical_proposals' if canonical else 'fold_proposals']
        npth=verify(src['canonical_proposals' if canonical else 'fold_proposals']);source=next(c for c in npool['cases'] if c['case_id']==cid)
        pair=[copy.deepcopy(next(t for t in proposals['tasks'] if t['case_id']==cid and t['metal']==z)) for z in ('Ca','La')]
        if pair[0]['active_mode_ids']!=pair[1]['active_mode_ids'] or len(pair[0]['active_indices'])!=4:raise InvalidArtifact('existing common four modes differ')
        cells={z:{} for z in ('Ca','La')};candidates=[]
        for q in BASE:
            for t in pair:
                z=t['metal'];old=source['matrix'][z][q];nr=read_json(verify(old['MACE']));rq=read_json(verify(nr['request']))
                if (old['status']!='complete' or nr['status']!='complete' or nr['model']!=m['model'] or nr['energy_eV']!=old['components']['MACE_eV'] or
                    rq['charge']!=t['charge'] or rq['multiplicity']!=1 or xyz(verify(rq['xyz']))!=xyz(verify(old['xyz']))):raise InvalidArtifact('native source state differs')
                verify(nr['forces']);force_count+=1
                if q=='origin':full=[0.]*t['mode_count']
                else:full=read_json(npth.parent/'proposals'/(cid+'__'+q.split('_')[-1])/'result.json')['proposal']['full_q']
                kin=Kinematics(read_json(verify(t['mapping']))['context']);atoms=xyz(verify(old['xyz']))
                if not np.allclose(kin.evaluate(full)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('source candidate map differs')
                angular.final_geometry(kin,t,np.array(full)[t['active_indices']],[a[0] for a in atoms])
                cell={'case_id':cid,'candidate':q,'metal':z,'xyz':old['xyz'],'charge':t['charge'],'multiplicity':1,
                      'native':old['MACE'],'forces':nr['forces'],'MACE_eV':nr['energy_eV'],'full_q':full,'low':{}}
                for medium in ('vacuum','alpb'):
                    row=external.get((cid,q,z,medium))
                    if row:cell['low'][medium]=external_reuse(row,cell,medium,m)
                    else:
                        tid='__'.join((cid,q,z,medium));lt=low_task(cell,medium,out/'baseline/tasks'/tid,False)
                        missing.append(lt);cell['low'][medium]={'status':'pending','task_id':tid}
                cells[z][q]=cell
            a,b=[xyz(verify(cells[z][q]['xyz'])) for z in ('Ca','La')]
            if a[1:]!=b[1:] or a[0][1:]!=b[0][1:]:raise InvalidArtifact('candidate rows do not share geometry')
            candidates.append({'id':q,'xyz':cells['Ca'][q]['xyz']})
        cases.append({'case_id':cid,'known_class_for_report_only':source['old_result']['expected_class'],
            'biological_group':cid.split('__conditioned_')[0],'native_static_R':source['old_result']['R0'],'native_adaptive':source['pool'],
            'cells':cells,'candidates':candidates})
        tasks.extend(pair)
    if len(missing)!=16 or force_count!=36:raise InvalidArtifact('expected baseline/force reuse counts differ')
    m.update(protocol_id=PROTOCOL,settings=SETTINGS,backend_settings=xtb.SETTINGS,sources=src,agreement=record(agreement),
       implementation=snapshot(Path(__file__).parent,out/'implementation'),cases=cases,tasks=tasks,baseline_tasks=missing,
       existing_baseline_cells=56,missing_baseline_cells=16,verified_native_force_cells=36,production_changed=False,
       optimizer_software={'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),'kernel':record(scipy_kernel.__file__)})
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['backend_settings']!=xtb.SETTINGS or [c['case_id'] for c in m['cases']]!=list(CASES):raise InvalidArtifact('fixed scientific scope differs')
    if len(m['tasks'])!=12 or len(m['baseline_tasks'])!=16 or m['verified_native_force_cells']!=36:raise InvalidArtifact('finite task population differs')
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer version changed')
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    for p in m['sources'].values():verify(p)
    for p in m['implementation'].values():verify(p)
    for k in ('agreement','software','executable','parameter','installed_parameter','package'):verify(m[k])
    if verify(m['parameter']).read_bytes()!=verify(m['installed_parameter']).read_bytes():raise InvalidArtifact('parameter changed')
    if {p.name for p in verify(m['parameter']).parent.iterdir()}!={'param_gfn2-xtb.txt'}:raise InvalidArtifact('parameter override')
    for c in m['cases']:
        pair=[t for t in m['tasks'] if t['case_id']==c['case_id']]
        if [t['metal'] for t in pair]!=['Ca','La'] or pair[0]['active_mode_ids']!=pair[1]['active_mode_ids']:raise InvalidArtifact('common basis/state changed')
        for t in pair:
            verify(t['mapping']);verify(t['source_preparation'])
            for cell in c['cells'][t['metal']].values():
                nr=read_json(verify(cell['native']));rq=read_json(verify(nr['request']));force=np.load(verify(cell['forces']),allow_pickle=False)
                if nr['model']!=m['model'] or nr['forces']!=cell['forces'] or nr['energy_eV']!=cell['MACE_eV'] or rq['charge']!=t['charge'] or xyz(verify(rq['xyz']))!=xyz(verify(cell['xyz'])) or force.shape!=(len(xyz(verify(cell['xyz']))),3):raise InvalidArtifact('actual force/source changed')
                if cell['charge']!=t['charge'] or cell['multiplicity']!=1:raise InvalidArtifact('state changed')
    return {'status':'validated','manifest':record(manifest),'sources':6,'searches':12,'baseline_reuses':56,'new_baseline_calls':16,
            'maximum_search_calls':960,'maximum_cross_calls':24,'maximum_new_standalone':1000,'maximum_new_MACE':480,'new_calls_in_validation':0}


def run_low(manifest,m,t):
    rp=Path(t['directory'])/'execution.json'
    if not rp.exists():xtb.execute_cell(manifest,m,t)
    r=xtb.parse_task(manifest,t)
    return r


def components(native,low):
    if any(low[s]['status']!='complete' for s in ('vacuum','alpb')):raise InvalidArtifact('required standalone cell failed: '+str({s:low[s].get('reason') for s in low if low[s]['status']!='complete'}))
    return {'MACE_eV':native,'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':low['alpb']['energy_hartree']}


class LimitReached(Exception):pass
class SharedGPU:
    def __init__(self,manifest):self.gpu=WarmGPU(manifest);self.lock=threading.Lock()
    def evaluate(self,pin):
        with self.lock:return self.gpu.evaluate(pin)
    def close(self):self.gpu.close()


class Evaluator:
    def __init__(self,manifest,m,t,start,gpu):
        self.path=Path(manifest);self.pin=record(manifest);self.m=m;self.t=t;self.start=start;self.gpu=gpu
        self.kin=Kinematics(read_json(verify(t['mapping']))['context']);self.atoms=xyz(verify(t['xyz']))
        self.root=self.path.parent/'searches'/t['task_id'];self.root.mkdir(parents=True,exist_ok=True)
        self.points=[];self.cache={}
    def evaluate(self,q):
        q=np.asarray(q,dtype=float);full=angular.full_q(self.t,q);key=cache_key({'manifest':self.pin,'task_id':self.t['task_id'],'q':q.tolist()})
        if key in self.cache:return self.cache[key]
        if len(self.points)>=40:raise LimitReached()
        if not np.isfinite(q).all() or max(abs(q))>.8+1e-12:raise InvalidArtifact('angle_bound_violation')
        d=self.root/'evaluations'/key;d.mkdir(parents=True,exist_ok=False)
        point={'status':'unavailable','active_q_radian':q.tolist(),'full_q':full.tolist(),'key':key,'new_MACE_calls':0,'new_standalone_calls':0}
        self.points.append(point)
        try:
            checks=geometry_check(self.kin,full,[a[0] for a in self.atoms])
            if not checks['pass']:raise InvalidArtifact('geometry_trial_failure')
            coord=self.kin.evaluate(full)[1];xp=d/'context.xyz'
            initial=np.array_equal(full,np.asarray(self.start['full_q']))
            if initial:shutil.copyfile(verify(self.start['xyz']),xp)
            else:write_xyz(xp,[(a[0],*p) for a,p in zip(self.atoms,coord)])
            if initial:native=read_json(verify(self.start['native']));nrecord=self.start['native']
            else:
                req=d/'request.json';write_new(req,{'manifest':self.pin,'task_id':self.t['task_id'],'xyz':record(xp),'charge':self.t['charge'],'multiplicity':1})
                point['new_MACE_calls']=1;native,nrecord=self.gpu.evaluate(record(req))
            point.update(native=nrecord,coordinate=record(xp))
            if native['status']!='complete':raise InvalidArtifact('native_MACE_failed')
            cell={**self.t,'xyz':record(xp),'candidate':'response'};lows=[low_task(cell,s,d/s,True) for s in ('vacuum','alpb')]
            lm={k:self.m[k] for k in ('executable','parameter','installed_parameter','package')};lm.update(tasks=lows,parent_manifest=self.pin)
            lp=d/'low_manifest.json';write_new(lp,lm);point['new_standalone_calls']=2
            with futures.ThreadPoolExecutor(max_workers=2) as ex:
                results=list(ex.map(lambda t:run_low(lp,lm,t),lows))
            low={r['medium']:r for r in results};point.update(low=low,low_manifest=record(lp));comp=components(native['energy_eV'],low)
            if any(not low[s]['gradient'] for s in low):raise InvalidArtifact('analytic_gradient_unavailable')
            if initial:
                for s,f in (('vacuum','GFN2_vacuum_hartree'),('alpb','GFN2_ALPB_hartree')):
                    if abs(comp[f]-self.start['components'][f])*HA_TO_KCAL>.1:raise InvalidArtifact('start_energy_replay_failed')
            jac=self.kin.evaluate(full)[3][self.t['active_indices']]
            ng,sg,g=project_components(jac,np.load(verify(native['forces']),allow_pickle=False),low['vacuum']['gradient']['kcal_mol_A'],low['alpb']['gradient']['kcal_mol_A'])
            physical=angular.physical_status(self.kin,self.t,q)
            point.update(status='complete',components=comp,coordinate=record(xp),native=nrecord,low=low,
                gradient_kcal_mol_rad=g.tolist(),native_gradient_kcal_mol_rad=ng.tolist(),solvent_gradient_kcal_mol_rad=sg.tolist(),
                physical=physical,geometry_checks=checks,eligible=physical['physical_feasible'])
        except Exception as exc:point['reason']=str(exc)
        write_new(d/'result.json',point);point['receipt']=record(d/'result.json');self.cache[key]=point
        if point['status']!='complete':raise InvalidArtifact(point['reason'])
        return point


def search(manifest,m,t,start,gpu):
    ev=Evaluator(manifest,m,t,start,gpu);rp=ev.root/'result.json'
    if rp.exists():raise InvalidArtifact('completed or failed search must not rerun')
    tick=time.monotonic();r={'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],'manifest':record(manifest),'status':'unavailable','reason':None,'candidate':None,'termination':None}
    try:
        first=ev.evaluate(np.asarray(start['full_q'])[t['active_indices']]);r['start']=first
        def fun(q):
            p=ev.evaluate(q);return relative_components(p['components'],first['components'])['composite_kcal_mol']/HA_TO_KCAL,np.array(p['gradient_kcal_mol_rad'])/HA_TO_KCAL
        try:
            opt=minimize(fun,np.asarray(start['full_q'])[t['active_indices']],jac=True,method='SLSQP',bounds=[(-.8,.8)]*4,
                constraints=[{'type':'ineq','fun':lambda q:angular.constraints(ev.kin,t,q)[0],'jac':lambda q:angular.constraints(ev.kin,t,q)[1]}],
                options={'maxiter':20,'ftol':SETTINGS['ftol_Ha_equivalent']})
            r['optimizer']={'success':bool(opt.success),'status':int(opt.status),'message':str(opt.message),'iterations':int(opt.nit),'function_evaluations':int(opt.nfev)}
            if not opt.success and opt.status!=9:raise InvalidArtifact('optimizer_failure: '+str(opt.message))
            r['termination']='converged' if opt.success else 'iteration_limit'
        except LimitReached:r['termination']='evaluation_limit'
        eligible=[p for p in ev.points if p.get('status')=='complete' and p.get('eligible')]
        best=min(eligible,key=lambda p:relative_components(p['components'],first['components'])['composite_kcal_mol'])
        angular.final_geometry(ev.kin,t,np.array(best['active_q_radian']),[a[0] for a in ev.atoms])
        r.update(status='finite_proposal_available',candidate=best,work_from_start=relative_components(best['components'],first['components']),
            new_geometry=xyz(verify(best['coordinate']))!=xyz(verify(start['xyz'])))
    except Exception as exc:
        reason=str(exc);kind=('geometry_failure' if 'geometry' in reason or 'angle_bound' in reason else 'SCF_or_gradient_failure' if 'standalone' in reason or 'gradient' in reason or 'replay' in reason else 'optimizer_failure' if 'optimizer' in reason else 'native_or_execution_failure')
        r.update(reason=reason,termination=kind)
    r.update(points=ev.points,wall_seconds=time.monotonic()-tick,new_MACE_calls=sum(x['new_MACE_calls'] for x in ev.points),
             new_standalone_calls=sum(x['new_standalone_calls'] for x in ev.points),stationary_minimum_claimed=False)
    write_new(rp,r);print(json.dumps({'task':t['task_id'],'status':r['status'],'termination':r['termination'],'points':len(ev.points),'reason':r['reason']}),flush=True);return r


def execute(manifest):
    validate(manifest);allocation();m=read_json(manifest);root=Path(manifest).parent;starttime=time.monotonic();gpu=None;error=None
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            with futures.ThreadPoolExecutor(max_workers=4) as ex:baseline=list(ex.map(lambda t:run_low(manifest,m,t),m['baseline_tasks']))
            write_new(root/'baseline_collection.json',{'manifest':record(manifest),'rows':baseline})
            pending={r['task_id']:r for r in baseline};cases=copy.deepcopy(m['cases']);starts={};searches=[]
            for c in cases:
                matrix={z:{} for z in ('Ca','La')}
                try:
                    for z in matrix:
                        for q,cell in c['cells'][z].items():
                            low={s:(v['row'] if v['status']=='complete' else pending[v['task_id']]) for s,v in cell['low'].items()}
                            cell['executed_low']=low
                            try:
                                comp=components(cell['MACE_eV'],low);cell['components']=comp
                                matrix[z][q]={'status':'complete','components':comp}
                            except InvalidArtifact as exc:matrix[z][q]={'status':'unavailable','reason':str(exc)}
                    prior=choose_rows(matrix,BASE);c.update(base_matrix=matrix,base_pool=prior)
                    if prior['status']!='available':raise InvalidArtifact('required starting three-candidate matrix incomplete')
                    c.update(status='prepared')
                    for z in ('Ca','La'):
                        name=prior['rows'][z]['mathematical_candidate'];t=next(t for t in m['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
                        starts[t['task_id']]=c['cells'][z][name]
                except Exception as exc:c.update(status='unavailable',reason='baseline_missing: '+str(exc))
            write_new(root/'STARTS.json',{'manifest':record(manifest),'starts':starts,'cases':cases})
            gpu=SharedGPU(manifest)
            with futures.ThreadPoolExecutor(max_workers=2) as ex:
                pending={ex.submit(search,manifest,m,t,starts[t['task_id']],gpu):t for t in m['tasks'] if t['task_id'] in starts}
                for f in futures.as_completed(pending):searches.append(f.result())
            write_new(root/'search_collection.json',{'manifest':record(manifest),'rows':searches})
            by={(r['case_id'],r['metal']):r for r in searches};rows=[]
            for c in cases:
                row={**c,'extended_pool':None,'extension_status':'unavailable'};matrix=None
                try:
                    if c['status']!='prepared':raise InvalidArtifact(c['reason'])
                    matrix=copy.deepcopy(c['base_matrix']);candidates=list(BASE);geometries={q:xyz(verify(c['cells']['Ca'][q]['xyz'])) for q in BASE};aliases={}
                    for maker in ('Ca','La'):
                        r=by[c['case_id'],maker]
                        if r['status']!='finite_proposal_available':raise InvalidArtifact('required search unavailable: '+maker)
                        p=r['candidate'];name='response_'+maker;atoms=xyz(verify(p['coordinate']));normalized=[('Ca',*atoms[0][1:]),*atoms[1:]]
                        duplicate=next((q for q,a in geometries.items() if a==normalized),None)
                        if duplicate:aliases[name]=duplicate;continue
                        geometries[name]=normalized;candidates.append(name)
                        matrix[maker][name]={'status':'complete','components':p['components'],'search_receipt':record(root/'searches'/(c['case_id']+'__'+maker)/'result.json')}
                        other='La' if maker=='Ca' else 'Ca';t=next(t for t in m['tasks'] if t['case_id']==c['case_id'] and t['metal']==other)
                        d=root/'cross'/(c['case_id']+'__'+name+'__'+other);d.mkdir(parents=True);xp=d/'context.xyz';write_xyz(xp,[(other,*atoms[0][1:]),*atoms[1:]])
                        req=d/'request.json';write_new(req,{'manifest':record(manifest),'task_id':t['task_id'],'xyz':record(xp),'charge':t['charge'],'multiplicity':1})
                        nr,pin=gpu.evaluate(record(req))
                        if nr['status']!='complete':raise InvalidArtifact('cross_MACE_failed')
                        cell={**t,'candidate':name,'xyz':record(xp)};lt=[low_task(cell,s,d/s,False) for s in ('vacuum','alpb')]
                        lm={k:m[k] for k in ('executable','parameter','installed_parameter','package')};lm.update(tasks=lt,parent_manifest=record(manifest));lp=d/'manifest.json';write_new(lp,lm)
                        with futures.ThreadPoolExecutor(max_workers=2) as ex:low=list(ex.map(lambda t:run_low(lp,lm,t),lt))
                        comp=components(nr['energy_eV'],{r['medium']:r for r in low});matrix[other][name]={'status':'complete','components':comp,'native':pin,'low':low}
                    row.update(extension_status='available',extended_pool=choose_rows(matrix,candidates),extended_matrix=matrix,candidate_ids=candidates,aliases=aliases)
                except Exception as exc:
                    row['extension_reason']=str(exc)
                    if matrix is not None:row['attempted_extended_matrix']=matrix
                rows.append(row)
            write_new(root/'collection.json',{'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'denominator':6,'available':sum(r['extension_status']=='available' for r in rows),'searches':record(root/'search_collection.json')})
        except Exception as exc:error=str(exc)
        finally:
            try:
                if gpu:gpu.close()
            except Exception as exc:error=(error+'; ' if error else '')+str(exc)
            elapsed=time.monotonic()-starttime;write_new(root/'EXECUTION.json',{'manifest':record(manifest),'error':error,'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'requested_GPU_seconds':elapsed})
    if error:raise InvalidArtifact(error)
    return {'collection':record(root/'collection.json')}


def report(collection,output):
    r=read_json(collection);m=read_json(verify(r['manifest']));ref=read_json(verify(m['sources']['external_reference']));bands=ref['variants']['static']['bands'];rows=[]
    search_rows=read_json(verify(r['searches']))['rows']
    def value(p):return p['operational']['composite_R_model_kcal_mol'] if p and p.get('status')=='available' and p.get('operational') else None
    for c in r['rows']:
        values={'static':value(choose_rows(c['base_matrix'],('origin',))) if c.get('base_matrix') else None,
                'prior_three':value(c.get('base_pool')),'solvent_response':value(c.get('extended_pool'))}
        response={}
        for z in ('Ca','La'):
            s=next((x for x in search_rows if x['case_id']==c['case_id'] and x['metal']==z),None)
            if s is None:response[z]={'status':'not_run'};continue
            response[z]={k:s.get(k) for k in ('status','termination','reason','work_from_start','new_MACE_calls','new_standalone_calls','optimizer')}
            response[z]['distinct_evaluations']=len(s['points'])
            response[z]['active_mode_ids']=next(t['active_mode_ids'] for t in m['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
            for stage in ('start','candidate'):
                point=s.get(stage)
                response[z][stage]=({k:point.get(k) for k in ('active_q_radian','gradient_kcal_mol_rad','native_gradient_kcal_mol_rad','solvent_gradient_kcal_mol_rad','physical')} if point else None)
        rows.append({'case_id':c['case_id'],'known_class_for_report_only':c['known_class_for_report_only'],'biological_group':c['biological_group'],'values':values,
           'static_band_calls_transferred_for_changed_models':{k:call(v,bands) for k,v in values.items()},'native_static_R':c['native_static_R'],'native_adaptive':c['native_adaptive'],
           'extension_status':c['extension_status'],'reason':c.get('extension_reason') or c.get('reason'),'response_components':response})
    gaps={}
    for method in ('static','prior_three','solvent_response'):
        if any(x['values'][method] is None for x in rows):gaps[method]=None;continue
        gaps[method]=min(x['values'][method] for x in rows if x['known_class_for_report_only']=='La')-max(x['values'][method] for x in rows if x['known_class_for_report_only']=='Ca')
    result={'protocol_id':PROTOCOL,'collection':record(collection),'rows':rows,'raw_class_gaps_model_kcal_mol':gaps,'case_denominator':6,'biological_groups':5,'own_response_reference':None,'source_selection':'inspected_development_failures','new_calibration':False,'production_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,names in {'prepare':('agreement','output'),'validate':('manifest',),'execute':('manifest',),'report':('collection','output')}.items():
        a=sub.add_parser(op)
        for n in names:a.add_argument('--'+n,required=True,type=Path)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
if __name__=='__main__':main()
