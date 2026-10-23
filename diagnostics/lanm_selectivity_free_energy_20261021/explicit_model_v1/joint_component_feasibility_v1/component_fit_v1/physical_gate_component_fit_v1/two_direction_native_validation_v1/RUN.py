#!/usr/bin/env python3
"""Two fixed native witnesses; no fit, solver, direction selection or promotion."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import socket
import sys
import time

HERE=Path(__file__).resolve().parent
LIMITS={'new_parameter_vectors':2,'new_native_programs':112,'workers':12,'solver_calls':0,'numerical_threads':1}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed authority '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def fingerprint(cfg):return digest({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation')})

def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['adapter']==cfg['implementation'],'Frozen runner')
    require(cfg['limits']==LIMITS and cfg['heldout_authorized'] is False,'Exact two-point scope')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical/native thread')
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful login')
    require(cfg['hostname']==socket.gethostname() and pin(sys.executable)==cfg['interpreter'],'Host/interpreter')
    inv=read(cfg['inventory']);inp=read(inv['inputs']);s=inp['sources']
    for p in s.values():verify(p)
    require(inv['status']=='prepared_two_direction_metadata_no_execution' and inv['implementation']==cfg['inventory_implementation'] and inv['plan']==cfg['plan'],'Exact inventory/plan')
    require(cfg['reuse_inventory']==cfg['inventory'],'Inherited immutable evidence alias')
    plan=read(cfg['plan']);panel=read(plan['panel_result']);affine=read(plan['affine_result']);root=read(s['panel_root_review']);audit=read(root['saved_review'])
    require(root['status']=='admitted_finite_final_local_attribution' and root['result']==plan['panel_result'] and root['affine_result']==plan['affine_result'] and root['next_plan']==cfg['plan'],'Root actual panel admission')
    require(audit['status']=='passed_final_local_attribution_saved_review' and audit['result']==plan['panel_result'] and audit['affine_result']==plan['affine_result'],'Actual saved panel review')
    require(panel['status']=='completed_final_local_attribution_panel' and panel['error'] is None,'Terminal measured derivatives')
    require(len(inv['rows'])==2 and [r['id'] for r in inv['rows']]==['B_half','D_half'] and inv['potential_new_programs']<=112,'Exact fixed pair')
    for row,declared in zip(inv['rows'],plan['rows']):
        require(all(row[k]==v for k,v in declared.items()),'Literal preregistered point')
        pred=next(x for x in affine['cases'][row['case']]['predictions'] if x['fraction']==row['fraction'])
        require(pred['parameters']==row['parameters'] and pred['key']==row['key']==digest(row['parameters']) and pred['literal_original_box_pass'],'Exact in-domain half point')
    center=read(s['selected']);require(root['center']==s['selected'] and digest(center['parameters'])==inv['center_key']==plan['center_key'],'Exact immutable center')
    parent=read(cfg['parent_configuration']);require(parent['implementation']==cfg['parent_evaluator']==s['evaluator'],'Frozen scientific evaluator')
    require(cfg['parent_configuration']==s['parent_configuration'] and cfg['core']==s['core'] and cfg['targets']==s['targets'],'Scientific input pins')
    core=module(cfg['core'],'two_direction_core');core.clarabel.DefaultSolver=lambda *a,**k:(_ for _ in ()).throw(ValueError('No solver authorized'))
    base=module(cfg['parent_evaluator'],'two_direction_evaluator');old=module(cfg['original_adapter'],'two_direction_parent_adapter')
    require(old.J.pin(old.J.__file__)==cfg['journal_implementation'],'Admitted journal')
    env=read(s['environment']);require(str(Path(sys.prefix).resolve())==env['prefix'],'Pinned environment')
    import importlib.metadata
    for name,version in env['versions'].items():require(importlib.metadata.version(name)==version,'Pinned package '+name)
    for p in env['files']:verify(p)
    return old,base,core,parent,inv,center,s

def evaluator_class(old,base):
    Parent=old.evaluator_class(base)
    class WitnessEvaluator(Parent):
        def __init__(self,cfg,path,mode,inv,center,sources):
            require(mode in ('prepare','run') and not Path(path).exists(),'Fresh finite cache')
            self.batch=copy.deepcopy(cfg);self.mode=mode;self.batch_fingerprint=fingerprint(cfg)
            self.allowed={r['key']:{'parameters':r['parameters'],'role':r['role'],'id':r['id']} for r in inv['rows']}
            self.complete_vectors={};self.imported_vector_keys={inv['center_key']}
            self.executor_identity={'mode':'login','job_id':None,'hostname':socket.gethostname(),'controller_pid':os.getpid(),'run_id':cfg['run_id'],'implementation':cfg['implementation'],'parent_implementation':cfg['parent_evaluator']}
            base.Evaluator.__init__(self,read(cfg['parent_configuration']),path,root_release=cfg.get('root_release'))
            self.program_attempts=[];self.cell_attempts={};self.started_cells=set()
            key=inv['center_key'];self.registry=copy.deepcopy(inv['registry'])
            self.loaded_saved_vectors={};self.saved_metal_index=copy.deepcopy(inv['saved_metal_index'])
            self.imported_cell_keys=set();self._import_saved(key,center)
            require(len(self.imported_cell_keys)==46,'Exact46 selected imported cells')
            self.ledger.update(new_policy=cfg['implementation'],limits=LIMITS,configuration_fingerprint=self.batch_fingerprint,
                imported_origin_vectors=1,imported_derivative_vectors=0,additional_pilot_vectors=0,available_saved_vectors=len(self.registry),
                initial_imported_unique_cells=46,native_returncodes_observed=0)
            self.persistence=old.J.Journal(self.cache_dir/'EVENTS.jsonl',self.ledger,self.lock,
                {'configuration_fingerprint':self.batch_fingerprint,'implementation':cfg['implementation'],'root_release':cfg.get('root_release')})
            self.persistence.materialize(self.ledger_path,expected=self.ledger)
            for row in inv['rows']:
                require(self.required_programs(row['parameters'])==row['potential_new_programs'],'Actual cache identities equal metadata budget')
                for metal,match in row['metal_cache_matches'].items():
                    key=match['matching_saved_key']
                    if key is not None:
                        value=self._import_saved(key)
                        for cid in match['cells']:
                            token,identity,_,meta=self.cache_identity(row['parameters'],cid);record=value['cells'][cid]
                            require(token==record['cache_key'] and identity==record['cache_identity'] and meta==record['request_provenance']['render'],'Exact reused scientific identity/render')
            self.expected_programs={row['key']:row['potential_new_programs'] for row in inv['rows']}

        def _gate(self):
            require(self.mode=='run' and not self.stop_event.is_set(),'No native work in preparation or after failure')
            require(self.ledger['native_program_invocations']<112,'Hard112 invocation-boundary ceiling')
            require(all(os.environ.get(k)=='1' for k in THREADS) and socket.gethostname()==self.batch['hostname'],'One-thread native declared host')
            require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'No fabricated Slurm metadata')
            if self._gated:return
            c=self.config
            for p in c['executables'].values():verify(p)
            build=read(c['build_receipt']);qual=read(c['qualification_review']);verify(c['source_manifest']);verify(qual['implementation']);verify(qual['result'])
            require(build['status']=='complete' and build['manifest']==c['source_manifest'] and all(build['executables'][k]==v for k,v in c['executables'].items()),'Exact admitted native build')
            require(qual['status']=='passed_local_response_native_qualification' and all(qual[k]==c[k] for k in ('source_manifest','build_receipt','executables')),'Actual native qualification')
            require(all(qual[k] for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Inherited numerical gates')
            check_release(self.batch);self._gated=True

        def authorize(self,*args,**kwargs):raise ValueError('Allowlist frozen by exact final inventory')

        def evaluate(self,values):
            key=digest(values)
            require(key in self.allowed and self.allowed[key]['parameters']==values,'Only two preregistered half vectors')
            require(key not in self.ledger['vector_keys'] and len(self.ledger['vector_keys'])<2,'No retry or third vector')
            needed=self.required_programs(values);require(needed==self.expected_programs[key] and self.ledger['native_program_invocations']+needed<=112,'Exact verified whole-vector budget')
            before=self.ledger['native_program_invocations'];answer=super().evaluate(values)
            require(self.ledger['native_program_invocations']-before==needed and len(answer['cells'])==46,'One complete fixed native witness')
            return answer

        def snapshot(self):
            with self.lock:self.persistence.materialize(self.ledger_path,expected=self.ledger)
            s=base.Evaluator.snapshot(self)
            s.update(max_parameter_vectors=2,max_native_programs=112,active_limits=LIMITS,
                new_parameter_vectors=len(self.ledger['vector_keys']),native_returncodes_observed=self.ledger['native_returncodes_observed'],
                imported_unique_cells=len(self.imported_cell_keys),imported_origin_vectors=1,imported_derivative_vectors=0,additional_pilot_vectors=0,
                executor=copy.deepcopy(self.executor_identity),configuration_fingerprint=self.batch_fingerprint)
            return s
    return WitnessEvaluator

def check_release(cfg):
    review=read(cfg['source_review']);prep=read(cfg['preparation']);release=read(cfg['root_release']);fp=fingerprint(cfg)
    require(review['status']=='passed_two_direction_native_validation_source_review' and review['configuration_fingerprint']==fp,'Independent source contract')
    require(prep['status']=='prepared_two_direction_zero_calls' and prep['configuration_fingerprint']==fp and prep['native_programs']==prep['solver_calls']==0,'Exact zero-call preparation')
    require(release['status']=='released_two_fixed_native_witnesses' and release['configuration_fingerprint']==fp and release['source_review']==cfg['source_review'] and release['preparation']==cfg['preparation'],'Exact root release')
    for k in ('implementation','core','inventory','plan','limits','output_directory','parent_configuration','original_adapter'):
        require(review[k]==cfg[k] and release[k]==cfg[k],'Exact source/release pin '+k)

def summarize(answer,prediction,core,fixed,targets):
    metrics=core.assess(answer,fixed);chi=answer['response']['chi_states_a0_cubed'];C=answer['response']['C_dipoles_Debye']
    increments={m:chi[m+'_trimer']-chi[m+'_dimer'] for m in ('La','Dy')};increments['Dy_minus_La']=increments['Dy']-increments['La']
    reference={m:targets['chi_a0_cubed'][m+'_trimer']-targets['chi_a0_cubed'][m+'_dimer'] for m in ('La','Dy')};reference['Dy_minus_La']=reference['Dy']-reference['La']
    return {'metrics':metrics,'predicted':prediction,'actual_minus_predicted_component_errors_kcal_mol':{k:v-prediction['component_errors_kcal_mol'][k] for k,v in metrics['component_errors_kcal_mol'].items()},
        'actual_minus_predicted_q':metrics['q']-prediction['actual_component_q'],
        'signed_chi_model_a0_cubed':chi,'signed_chi_errors_a0_cubed':{k:v-targets['chi_a0_cubed'][k] for k,v in chi.items()},
        'chi_conditional_increments_a0_cubed':{k:{'model':v,'reference':reference[k],'error':v-reference[k]} for k,v in increments.items()},
        'W9_C_Dy_minus_La_Debye':{'model':[a-b for a,b in zip(C['aquo9_central_Dy'],C['aquo9_central_La'])],'reference':[a-b for a,b in zip(targets['C_Debye']['aquo9_central_Dy'],targets['C_Debye']['aquo9_central_La'])]},
        'physical_failures_are_diagnostic':True,'model_admitted':False,'selected_candidate':None}


def run(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg=json.loads(Path(config_path).read_text());out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    old=None;ev=None;result={'status':'started','mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),
        'configuration_fingerprint':fingerprint(cfg),'native_programs':0,'solver_calls':0,'error':None,'vectors':{},
        'limits':LIMITS,'selected_candidate':None,'reference_calls':0,'heldout_reads':0,'job_id':None,
        'rows':{k:{'status':'not_attempted'} for k in ('B_half','D_half')}}
    def persist():
        if ev is not None:result['evaluator']=ev.snapshot();result['native_programs']=ev.ledger['native_program_invocations']
        if old is not None:old.save(out/'RESULT.json',result)
        else:(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    try:
        old,base,core,parent,inv,center,sources=context(cfg)
        if mode=='run':require(str(out)==cfg['output_directory'],'Exact immutable run output');check_release(cfg)
        ev=evaluator_class(old,base)(cfg,out/'cache',mode,inv,center,sources)
        result.update(center=sources['selected'],center_key=inv['center_key'],inventory=cfg['inventory'],plan=cfg['plan'],source_review=cfg.get('source_review'),root_release=cfg.get('root_release'))
        declaration=read(cfg['campaign_declaration']);origin=read(declaration['start']['artifact'])
        fixed=core.make_fixed(origin,ev.m,read(sources['decomposition']),ev.domain,declaration)
        result['center_metrics']=core.assess(center,fixed);require(result['center_metrics']['all25_original_physical_gates'],'Admitted final center25')
        affine=read(read(cfg['plan'])['affine_result']);targets=read(cfg['targets']);persist()
        if mode=='prepare':
            result.update(status='prepared_two_direction_zero_calls',preparation_checks={'two_exact_half_vectors':True,'46_selected_cells_exact':True,'verified_program_budget':sum(ev.expected_programs.values()),'no_solver_or_native':True})
        else:
            (out/'vectors').mkdir()
            for row in inv['rows']:
                record={'status':'started','parameters':row['parameters'],'id':row['id']};result['vectors'][row['key']]=record;result['rows'][row['id']]=record;persist()
                try:
                    answer=ev.evaluate(row['parameters']);p=out/'vectors'/(row['key']+'.json');old.save(p,answer)
                    record.update(status='complete',artifact=pin(p));persist()
                    prediction=next(x for x in affine['cases'][row['case']]['predictions'] if x['fraction']==row['fraction'])
                    record['comparison']=summarize(answer,prediction,core,fixed,targets);persist()
                except Exception as exc:record.update(status='failed',error=repr(exc));raise
            require(len(result['vectors'])==2 and result['native_programs']==inv['potential_new_programs']<=112,'Exact two fixed endpoints and actual cache budget')
            result['status']='completed_two_direction_native_validation'
    except Exception as exc:
        if ev is not None:ev._stop(exc)
        result.update(status='failed_two_direction_native_validation' if mode=='run' else 'failed_preparation',error=repr(exc))
    finally:
        now=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,
            'child_user_CPU_seconds':now.ru_utime-child.ru_utime,'child_system_CPU_seconds':now.ru_stime-child.ru_stime,
            'native_workers':12,'numerical_threads':1,'scheduler_allocated_CPU_seconds':0,'execution_mode':'login','inherited_costs_not_added':True};persist()
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'native_programs':result['native_programs'],'solver_calls':0,'error':result['error']}))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    raise SystemExit(bool(run(a.mode,a.config,a.out)['error']))
