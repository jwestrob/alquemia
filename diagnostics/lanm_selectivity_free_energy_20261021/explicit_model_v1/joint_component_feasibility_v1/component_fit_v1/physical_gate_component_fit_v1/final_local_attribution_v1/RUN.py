#!/usr/bin/env python3
"""Fixed24 derivative batch followed by four descriptive affine cases; no fit loop."""
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
LIMITS={'new_parameter_vectors':24,'new_native_programs':672,'workers':12,'solver_calls':4,'numerical_threads':1}
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
    require(cfg['implementation']==pin(__file__) and cfg['solver']==pin(HERE/'SOLVE.py'),'Frozen implementation')
    for name in ('specification','plan_ready','root_plan_review'):verify(cfg[name])
    require(cfg['limits']==LIMITS and cfg['heldout_authorized'] is False,'Exact finite scope')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical/native thread')
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful login execution')
    require(cfg['hostname']==socket.gethostname() and pin(sys.executable)==cfg['interpreter'],'Host/interpreter')
    inv=read(cfg['inventory']);inp=read(inv['inputs']);s=inp['sources']
    require(inv['implementation']==cfg['inventory_implementation'] and inv['status']=='prepared_final_local_probe_inventory_no_execution','Fixed inventory producer')
    for p in s.values():verify(p)
    require(cfg['settings']==s['settings'] and cfg['targets']==s['targets'] and cfg['environment']==s['environment'],'Fixed solver/reference authorities')
    require(cfg['reuse_inventory']==cfg['inventory'] and cfg['adapter']==cfg['implementation'],'Inherited immutable-evidence aliases')
    for k in ('parent_configuration','core','evaluator'):
        require(cfg[k if k!='evaluator' else 'parent_evaluator']==s[k],'Scientific input '+k)
    require(inv['potential_new_programs_after_exact_metadata_reuse']==672 and inv['full_vector_matches']==0,'Exact24 missing halves')
    fit=read(s['fit_result']);audit=read(s['fit_review']);center=read(s['selected'])
    require(fit['status']=='bounded_physical_gate_component_fit_complete' and fit['error'] is None,'Terminal fit')
    require(audit['status']=='passed_physical_gate_component_fit_actual_review' and audit['result']==s['fit_result'] and audit['selected']==s['selected'] and audit['selected_native_admission_passed'],'Terminal actual admission')
    require(fit['selected']==s['selected'] and fit['selected_key']==inv['center_key']==digest(center['parameters']),'Exact final state')
    require(fit['selected_metrics']['all25_original_physical_gates'],'Center25 physical gates')
    parent=read(cfg['parent_configuration']);require(parent['implementation']==cfg['parent_evaluator'],'Frozen scientific evaluator')
    core=module(cfg['core'],'final_attribution_core');base=module(cfg['parent_evaluator'],'final_attribution_evaluator')
    old=module(cfg['original_adapter'],'final_attribution_parent_adapter')
    require(old.J.pin(old.J.__file__)==cfg['journal_implementation'],'Admitted journal')
    require([r['parameters'] for r in inv['probe_rows']]==core.probe_points(center['parameters'],read(parent['domain'])),'Exact final24 coordinates')
    require(len(inv['probe_rows'])==24 and len({r['key'] for r in inv['probe_rows']})==24,'Unique24')
    for row in inv['probe_rows']:
        require(digest(row['parameters'])==row['key'] and row['potential_new_programs']==28,'One-metal28 per probe')
        matches=[x['matching_saved_vector'] for x in row['metal_cache_matches'].values()]
        require(sum(x is not None for x in matches)==1 and all(x is None or (x['key']==inv['center_key'] and x['artifact']==s['selected']) for x in matches),'Only selected unchanged-metal reuse')
    env=read(s['environment']);require(str(Path(sys.prefix).resolve())==env['prefix'],'Pinned existing numerical environment')
    import importlib.metadata
    for name,version in env['versions'].items():require(importlib.metadata.version(name)==version,'Pinned package '+name)
    for p in env['files']:verify(p)
    return old,base,core,parent,inv,center,s

def evaluator_class(old,base):
    Parent=old.evaluator_class(base)
    class FinalEvaluator(Parent):
        def __init__(self,cfg,path,mode,inv,center,sources):
            require(mode in ('prepare','run') and not Path(path).exists(),'Fresh finite cache')
            self.batch=copy.deepcopy(cfg);self.mode=mode;self.batch_fingerprint=fingerprint(cfg)
            self.allowed={r['key']:{'parameters':r['parameters'],'role':'final_probe','ordinal':r['ordinal']} for r in inv['probe_rows']}
            self.complete_vectors={};self.imported_vector_keys={inv['center_key']}
            self.executor_identity={'mode':'login','job_id':None,'hostname':socket.gethostname(),'controller_pid':os.getpid(),'run_id':cfg['run_id'],'implementation':cfg['implementation'],'parent_implementation':cfg['parent_evaluator']}
            base.Evaluator.__init__(self,read(cfg['parent_configuration']),path,root_release=cfg.get('root_release'))
            self.program_attempts=[];self.cell_attempts={};self.started_cells=set()
            key=inv['center_key'];self.registry={key:{'parameters':center['parameters'],'artifact':sources['selected'],'kind':'vector','source':sources['fit_result']}}
            self.loaded_saved_vectors={};self.saved_metal_index={m:{digest(base.effective(center['parameters'],m,self.domain)):key} for m in ('La','Dy')}
            self.imported_cell_keys=set();self._import_saved(key,center)
            require(len(self.imported_cell_keys)==46,'Exact46 selected imported cells')
            self.ledger.update(new_policy=cfg['implementation'],limits=LIMITS,configuration_fingerprint=self.batch_fingerprint,
                imported_origin_vectors=1,imported_derivative_vectors=0,additional_pilot_vectors=0,available_saved_vectors=1,
                initial_imported_unique_cells=46,native_returncodes_observed=0)
            self.persistence=old.J.Journal(self.cache_dir/'EVENTS.jsonl',self.ledger,self.lock,
                {'configuration_fingerprint':self.batch_fingerprint,'implementation':cfg['implementation'],'root_release':cfg.get('root_release')})
            self.persistence.materialize(self.ledger_path,expected=self.ledger)
            for row in inv['probe_rows']:
                require(self.required_programs(row['parameters'])==28,'Verified unchanged-metal cache requires exactly28 programs')
                for cid,record in center['cells'].items():
                    metal=self.m['rows'][cid]['metal'];match=row['metal_cache_matches'][metal]['matching_saved_vector']
                    if match is not None:
                        token,identity,_,meta=self.cache_identity(row['parameters'],cid)
                        require(token==record['cache_key'] and identity==record['cache_identity'] and token in self.cache,'Exact unchanged-metal scientific cache identity')
                        require(meta==record['request_provenance']['render'],'Exact unchanged-metal render')

        def _gate(self):
            require(self.mode=='run' and not self.stop_event.is_set(),'No native work in preparation or after failure')
            require(self.ledger['native_program_invocations']<672,'Hard672 invocation-boundary ceiling')
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
            require(key in self.allowed and self.allowed[key]['parameters']==values,'Only exact24 final derivative vectors')
            require(key not in self.ledger['vector_keys'] and len(self.ledger['vector_keys'])<24,'No retry or25th vector')
            require(self.required_programs(values)==28 and self.ledger['native_program_invocations']+28<=672,'Exact verified whole-vector budget')
            before=self.ledger['native_program_invocations'];answer=super().evaluate(values)
            require(self.ledger['native_program_invocations']-before==28 and len(answer['cells'])==46,'One admitted28-program derivative')
            return answer

        def snapshot(self):
            with self.lock:self.persistence.materialize(self.ledger_path,expected=self.ledger)
            s=base.Evaluator.snapshot(self)
            s.update(max_parameter_vectors=24,max_native_programs=672,active_limits=LIMITS,
                new_parameter_vectors=len(self.ledger['vector_keys']),native_returncodes_observed=self.ledger['native_returncodes_observed'],
                imported_unique_cells=46,imported_origin_vectors=1,imported_derivative_vectors=0,additional_pilot_vectors=0,
                executor=copy.deepcopy(self.executor_identity),configuration_fingerprint=self.batch_fingerprint)
            return s
    return FinalEvaluator

def check_release(cfg):
    review=read(cfg['source_review']);prep=read(cfg['preparation']);release=read(cfg['root_release']);fp=fingerprint(cfg)
    require(review['status']=='passed_final_local_attribution_source_review' and review['configuration_fingerprint']==fp,'Independent source contract')
    require(prep['status']=='prepared_final_local_attribution_zero_calls' and prep['configuration_fingerprint']==fp and prep['native_programs']==prep['solver_calls']==0,'Exact zero-call preparation')
    require(release['status']=='released_one_final_local_attribution_panel' and release['configuration_fingerprint']==fp and release['source_review']==cfg['source_review'] and release['preparation']==cfg['preparation'],'Exact root release')
    for k in ('implementation','solver','inventory','limits','output_directory','parent_configuration','original_adapter'):
        require(review[k]==cfg[k] and release[k]==cfg[k],'Exact source/release pin '+k)

def run(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg=json.loads(Path(config_path).read_text());out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    old=None;ev=None;result={'status':'started','mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),
        'configuration_fingerprint':fingerprint(cfg),'native_programs':0,'solver_calls':0,'error':None,'vectors':{},
        'limits':LIMITS,'selected_candidate':None,'native_trials':0,'reference_calls':0,'heldout_reads':0,'job_id':None}
    def persist():
        if ev is not None:result['evaluator']=ev.snapshot();result['native_programs']=ev.ledger['native_program_invocations']
        if old is not None:old.save(out/'RESULT.json',result)
        else:(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    try:
        old,base,core,parent,inv,center,sources=context(cfg)
        if mode=='run':require(str(out)==cfg['output_directory'],'Exact immutable run output');check_release(cfg)
        ev=evaluator_class(old,base)(cfg,out/'cache',mode,inv,center,sources)
        result.update(center=sources['selected'],center_key=inv['center_key'],inventory=cfg['inventory'],source_review=cfg.get('source_review'),root_release=cfg.get('root_release'))
        declaration=read(cfg['campaign_declaration']);origin=read(declaration['start']['artifact'])
        fixed=core.make_fixed(origin,ev.m,read(sources['decomposition']),ev.domain,declaration)
        metrics=core.assess(center,fixed);require(metrics['all25_original_physical_gates'],'Fixed final center25')
        result['center_metrics']=metrics;persist()
        if mode=='prepare':
            result.update(status='prepared_final_local_attribution_zero_calls',preparation_checks={'24_exact_probes':True,'46_selected_cells_exact':True,'672_after_verified_reuse':True,'no_solver_or_native':True})
        else:
            (out/'vectors').mkdir();probes=[]
            for row in inv['probe_rows']:
                record={'status':'started','parameters':row['parameters'],'ordinal':row['ordinal']};result['vectors'][row['key']]=record;persist()
                try:
                    answer=ev.evaluate(row['parameters']);p=out/'vectors'/(row['key']+'.json');old.save(p,answer)
                    record.update(status='complete',artifact=pin(p));probes.append(answer);persist()
                except Exception as exc:record.update(status='failed',error=repr(exc));raise
            require(len(probes)==24 and result['native_programs']==672,'Complete finite derivative panel before solves')
            solver=module(cfg['solver'],'final_attribution_solver')
            result['affine_diagnostics']=solver.run_cases(core,center,probes,fixed,
                dict(cfg,configuration_fingerprint_for_solver=fingerprint(cfg)),out/'affine_cases')
            result['solver_calls']=result['affine_diagnostics']['solver_calls']
            require(result['solver_calls']==4,'Four attempted affine cases')
            result['status']='completed_final_local_attribution_panel'
    except Exception as exc:
        if ev is not None:ev._stop(exc)
        result.update(status='failed_final_local_attribution_panel' if mode=='run' else 'failed_preparation',error=repr(exc))
    finally:
        # The solver independently saves every attempt/return before reporting.
        ap=out/'affine_cases/RESULT.json'
        if ap.exists():result['affine_result']=pin(ap);result['solver_calls']=json.loads(ap.read_text())['solver_calls']
        now=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,
            'child_user_CPU_seconds':now.ru_utime-child.ru_utime,'child_system_CPU_seconds':now.ru_stime-child.ru_stime,
            'native_workers':12,'numerical_threads':1,'scheduler_allocated_CPU_seconds':0,'execution_mode':'login','inherited_costs_not_added':True}
        persist()
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'native_programs':result['native_programs'],'solver_calls':result['solver_calls'],'error':result['error']}))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    raise SystemExit(bool(run(a.mode,a.config,a.out)['error']))
