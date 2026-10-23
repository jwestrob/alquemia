#!/usr/bin/env python3
"""Exactly one frozen consumed-data native vector; unchanged scientific evaluator."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import socket
import sys
import threading
import time

HERE=Path(__file__).resolve().parent
LIMITS={'new_parameter_vectors':1,'native_program_attempts':56,'workers':12,'aggregate_login_thread_budget':16}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')


def require(ok,msg):
    if not ok:raise ValueError(msg)


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])==p,'Changed pin '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def save(p,x):
    p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)


def module(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def fingerprint(cfg):return digest({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation')})


def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Batch source/specification')
    require(cfg['limits']==LIMITS and cfg['reserved_outcomes_authorized'] is False,'Finite declared scope')
    for k in ('declaration','parent_evaluator','parent_configuration','half_review','failed_proposal','origin','fit_review','root_fit_review','decomposition','component_diagnostic'):
        verify(cfg[k])
    parent=read(cfg['parent_configuration']);require(parent['implementation']==cfg['parent_evaluator'],'Unchanged parent evaluator')
    declaration=read(cfg['declaration'])
    base=module(cfg['parent_evaluator'],'frozen_joint_evaluator')
    cert=read(cfg['half_review']);failed=read(cfg['failed_proposal']);origin=read(cfg['origin']);fitreview=read(cfg['fit_review']);root=read(cfg['root_fit_review'])
    require(cert['status']=='passed_saved_half_linear_feasibility_review' and cert['strict_linear_feasibility_verified'] is True and cert['all25_exact_threshold_pass'] is True,'Independent strict half certificate')
    require(cert['original_result']==cfg['failed_proposal'] and cert['origin_artifact']==cfg['origin'] and cert['original_protocol_passed'] is False,'Certificate original failure/origin joins')
    verify(cert['implementation']);require(failed['status']=='failed_numerical_local_proposal' and failed['solver']['success'] is False and failed['solver_calls']==1,'Original failed solve preserved')
    require(fitreview['status']=='passed_terminal_joint_fit_review_continuation' and root['status']=='admitted_completed_frozen_fit_provenance_and_numerical_review' and root['independent_review']==cfg['fit_review'],'Completed fit admission')
    require(fitreview['fit_result']==root['fit_result']==failed['sources']['artifacts']['fit'],'Exact original fit lineage')
    fit=read(root['fit_result']);require(fit['vectors'][cert['origin']]['artifact']==cfg['origin'] and failed['origin']==cert['origin'] and origin['status']=='complete','Immutable admitted origin')
    values=base.checked(cert['physical_half_parameters'],read(parent['domain']))
    require(values==cfg['candidate']['parameters'] and digest(values)==cfg['candidate']['key']==cert['half_parameter_token'] and cfg['candidate']['fraction']==.5,'One exact certified half vector')
    require(values==[x+w*(d/2) for x,w,d in zip(failed['origin_parameters'],[p['upper']-p['lower'] for p in read(parent['domain'])['parameters']],failed['normalized_direction'])],'Literal half arithmetic')
    require(declaration['origin_artifact']==cfg['origin'] and declaration['half_certificate']==cfg['half_review'] and declaration['failed_proposal']==cfg['failed_proposal'] and declaration['physical_parameters']==values and declaration['parameter_token']==cfg['candidate']['key'],'Exact development declaration')
    require(declaration['original_protocol_passed'] is False and declaration['execution']['launcher']=='authorized_login' and declaration['execution']['vector_cap']==1 and declaration['execution']['native_attempt_cap']==56 and declaration['execution']['native_workers']==12 and declaration['execution']['maximum_owned_aggregate_threads']==16,'Declared scope/resources')
    require(len(origin['cells'])==46 and all(v['status']=='complete' for v in origin['cells'].values()),'Admitted origin46 denominator')
    fitcfg=read(fit['configuration'])
    require(fitcfg['evaluator_config']==cfg['parent_configuration'] and fitcfg['evaluator']==cfg['parent_evaluator'],'Actual original scientific evaluator configuration')
    dec=read(cfg['decomposition'])
    require(dec['status']=='complete_saved_pair_collective_diagnostic' and dec['sources']['evaluator_manifest']==parent['manifest'] and dec['sources']['targets']==parent['target_inventory'],'Exact literal component/static authority')
    require(failed['sources']['artifacts']['diagnostic']==cfg['component_diagnostic'],'Exact saved linear prediction authority')
    require(cfg['executor']=={'mode':'login','hostname':socket.gethostname(),'workers':12,'numerical_threads':1,'aggregate_thread_budget':16},'Declared truthful login envelope')
    require(pin(sys.executable)==cfg['interpreter'],'Exact existing interpreter')
    require(all(os.environ.get(k)=='1' for k in THREADS) and os.environ.get('BLIS_NUM_THREADS','1')=='1','Single native/numerical threads')
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Login execution cannot impersonate Slurm')
    return base,parent,cert,origin


class NotAttempted(RuntimeError):pass


def evaluator_class(base):
    class FiniteEvaluator(base.Evaluator):
        def __init__(self,cfg,cache_dir,mode):
            require(mode in ('prepare','run'),'Finite mode')
            self.batch=copy.deepcopy(cfg);self.batch_fingerprint=fingerprint(cfg);self.mode=mode
            self.allowed_key=cfg['candidate']['key'];self.allowed_values=list(cfg['candidate']['parameters'])
            self.executor_identity={'mode':'login','hostname':socket.gethostname(),'controller_pid':os.getpid(),'run_id':cfg['run_id'],'job_id':None,'implementation':pin(__file__),'parent_implementation':cfg['parent_evaluator']}
            require(not Path(cache_dir).exists(),'Fresh batch cache only; no old or resumed ledger')
            super().__init__(read(cfg['parent_configuration']),cache_dir,root_release=cfg.get('root_release'))
            self.started_cells=set();self.program_attempts=[];self.cell_attempts={}
            self.ledger.update(batch_fingerprint=self.batch_fingerprint,batch_limits=LIMITS,executor=self.executor_identity,native_returncodes_observed=0)
            save(self.ledger_path,self.ledger)

        def allowed(self,values):
            checked=self.checked_values(values)
            require(checked==self.allowed_values and digest(checked)==self.allowed_key,'Only exact fixed certified candidate is authorized')
            return checked

        def _stop(self,exc):
            self.stop_event.set()
            with self.lock:
                self.ledger['stopped']=True
                self.ledger.setdefault('failure',repr(exc))
                save(self.ledger_path,self.ledger)

        def _gate(self):
            require(self.mode=='run','Preparation forbids native calls')
            require(not self.stop_event.is_set(),'Global native stop active')
            require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'No fabricated scheduler identity')
            require(all(os.environ.get(k)=='1' for k in THREADS) and os.environ.get('BLIS_NUM_THREADS','1')=='1','One native numerical thread')
            require(socket.gethostname()==self.batch['executor']['hostname'],'Declared host')
            if self._gated:return
            cfg=self.config
            require(self.root_release is not None and all(isinstance(cfg.get(k),dict) for k in ('source_manifest','build_receipt','executables','qualification_review')),'Actual reviewed native bindings unavailable')
            require(set(cfg['executables'])=={'analyze','testgrad'},'Exact native executables')
            for q in cfg['executables'].values():verify(q)
            build=read(cfg['build_receipt']);qualification=read(cfg['qualification_review'])
            verify(qualification['implementation']);verify(qualification['result']);verify(cfg['source_manifest'])
            require(build['status']=='complete' and build['manifest']==cfg['source_manifest'] and all(build['executables'][k]==v for k,v in cfg['executables'].items()),'Completed build/source/executable joins')
            require(qualification['status']=='passed_local_response_native_qualification' and all(qualification[k]==cfg[k] for k in ('source_manifest','build_receipt','executables')),'Actual native qualification joins')
            require(all(qualification[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'All four native qualification flags')
            review=read(self.batch['source_review']);prep=read(self.batch['preparation']);release=read(self.root_release)
            require(review['status']=='passed_fixed_half_native_source_review' and review['implementation']==pin(__file__) and review['configuration_fingerprint']==self.batch_fingerprint,'Independent finite-source review')
            require(prep['status']=='prepared_fixed_half_native_zero_calls' and prep['configuration_fingerprint']==self.batch_fingerprint and prep['native_programs']==0,'Exact zero-call preparation')
            require(release['status']=='released_one_fixed_half_native_pilot' and release['implementation']==pin(__file__) and release['configuration_fingerprint']==self.batch_fingerprint,'Exact finite root release')
            for k in ('parent_configuration','half_review','origin','declaration','source_review','preparation','candidate','limits','executor'):
                require(release[k]==self.batch[k],'Root release binding '+k)
            require(release['output_directory']==self.batch['output_directory'] and release['cache_directory']==str(self.cache_dir) and release['reserved_outcomes_authorized'] is False,'Exact fresh output/cache and consumed scope')
            require(self.batch['output_directory']==str(self.cache_dir.parent),'Cache lies in declared new output')
            self._gated=True

        def _invoke(self,program,folder,record):
            with self.slots:
                with self.lock:
                    if self.stop_event.is_set():raise NotAttempted('Global stop before invocation')
                    self._gate();require(self.ledger['native_program_invocations']<56,'Finite56 native-attempt ceiling')
                    attempt={'number':len(self.program_attempts)+1,'cell_id':record['cell_id'],'program':program,'folder':str(folder),'status':'invoke_boundary_entered'}
                    self.program_attempts.append(attempt);self.ledger['native_program_invocations']+=1
                    self.cell_attempts[record['cell_id']]=self.cell_attempts.get(record['cell_id'],0)+1
                    self.ledger['native_attempts']=copy.deepcopy(self.program_attempts);save(self.ledger_path,self.ledger)
                args=['system.xyz','EM'] if program=='analyze' else ['system.xyz','Y','N']
                before=len(record['receipts'])
                try:
                    self.h['native_helper'].invoke([self.config['executables'][program]['path']]+args,folder,folder/(program+'.out'),record['receipts'])
                    require(len(record['receipts'])==before+1 and record['receipts'][-1].get('returncode')==0,'Exactly one successful native receipt')
                    with self.lock:attempt['status']='returned_success'
                except Exception as exc:
                    self.stop_event.set()
                    with self.lock:attempt.update(status='failed_native_attempt',error=repr(exc))
                    self._stop(exc);raise
                finally:
                    with self.lock:
                        for receipt in record['receipts'][before:]:receipt['executor']=copy.deepcopy(self.executor_identity)
                        attempt['receipt_count']=len(record['receipts'])-before
                        attempt['returncode_observed']=any('returncode' in r for r in record['receipts'][before:])
                        self.ledger['native_returncodes_observed']=sum(x.get('returncode_observed',False) for x in self.program_attempts)
                        self.ledger['native_attempts']=copy.deepcopy(self.program_attempts);save(self.ledger_path,self.ledger)

        def _cell(self,values,cid):
            values=self.allowed(values)
            if self.stop_event.is_set():raise NotAttempted('Untouched cell after global stop')
            cachekey,identity,key,meta=self.cache_identity(values,cid);row=self.m['rows'][cid]
            with self.lock:lock=self.cell_locks.setdefault(cachekey,threading.Lock())
            with lock:
                if self.stop_event.is_set():raise NotAttempted('Untouched cell after global stop')
                require(cachekey not in self.cache,'No repeated cell request in one-vector batch')
                require(base.effective(values,row['metal'],self.domain)!=base.effective(self.initial_vector,row['metal'],self.domain),'Both metals changed; no implicit baseline substitution')
                directory=self.cache_dir/'cells'/cachekey;rp=directory/'RESULT.json'
                require(not directory.exists(),'Fresh per-cell output only; failed/existing cells are never retried')
                with self.lock:self._gate()
                directory.mkdir(parents=True,exist_ok=False)
                record={'status':'started','cache_key':cachekey,'cache_identity':identity,'cell_id':cid,'first_requested_vector':values[:],
                        'inherited_native':False,'implementation':self.config['implementation'],'render':meta,'receipts':[],
                        'program_inputs':{},'job_id':None,'executor':copy.deepcopy(self.executor_identity),'root_release':self.root_release}
                with self.lock:self.started_cells.add(cid)
                try:
                    for program in (('analyze','testgrad') if row['needs_gradient'] else ('analyze',)):
                        if self.stop_event.is_set():raise NotAttempted('Stop before next native program')
                        folder=directory/program;folder.mkdir();text=key if program=='analyze' else self.render(values,cid,exports=False)[0]
                        for k,n in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/n)
                        (folder/'system.key').write_text(text)
                        record['program_inputs'][program]={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                        save(rp,record);self._invoke(program,folder,record)
                        parsed=self._parse(record,row);record.update(parsed);save(rp,record)
                    require(not row['needs_gradient'] or 'model' in record,'All ten native analytic-force rows')
                    record['status']='complete';self.cache[cachekey]=record
                except Exception as exc:
                    self._stop(exc);record.update(status='not_attempted' if isinstance(exc,NotAttempted) and not record['receipts'] else 'failed',error=repr(exc));raise
                finally:save(rp,record)
                result=copy.deepcopy(record);result['request_provenance']={'vector':values[:],'cache_key':cachekey,'cache_identity':identity,'render':meta,'executor':copy.deepcopy(self.executor_identity)}
                return result

        def evaluate(self,values):
            values=self.allowed(values)
            with self.evaluation_lock:
                with self.lock:
                    require(not self.ledger['vector_keys'] and not self.stop_event.is_set(),'Exactly one attempt; no second vector or retry')
                    self._gate();self.ledger['vector_keys'].append(self.allowed_key)
                    request={'values':values[:],'vector_key':self.allowed_key,'status':'started','cells':{cid:{'status':'not_attempted'} for cid in self.m['rows']}}
                    self.ledger['requests'].append(request);save(self.ledger_path,self.ledger)
                cells={};first_error=None;stopped_error=None
                try:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
                        futures={pool.submit(self._cell,values,cid):cid for cid in self.m['rows']}
                        for future in concurrent.futures.as_completed(futures):
                            cid=futures[future]
                            try:
                                value=future.result();cells[cid]=value;state={'status':'complete','cache_key':value['cache_key'],'inherited_native':False}
                            except NotAttempted as exc:
                                attempts=self.cell_attempts.get(cid,0)
                                state={'status':'failed_partial' if attempts else 'not_attempted','error':repr(exc),'native_attempts':attempts,'prepared_without_invocation':cid in self.started_cells and attempts==0}
                                if stopped_error is None:stopped_error=exc
                            except Exception as exc:
                                self._stop(exc);state={'status':'failed','error':repr(exc)}
                                if first_error is None:first_error=exc
                            with self.lock:request['cells'][cid]=state;save(self.ledger_path,self.ledger)
                    if first_error is not None:raise first_error
                    if stopped_error is not None:raise stopped_error
                    answer=self._assemble(cells)
                    require(len(cells)==46 and self.ledger['native_program_invocations']==56 and self.ledger['native_returncodes_observed']==56,'Exact46 cells/56 completed native programs')
                    with self.lock:request.update(status='complete',objective=answer['score']['objective']);save(self.ledger_path,self.ledger)
                    return answer
                except Exception as exc:
                    self._stop(exc)
                    with self.lock:request.update(status='failed',error=repr(exc));save(self.ledger_path,self.ledger)
                    raise

        def snapshot(self):
            answer=super().snapshot()
            answer.update(executor=copy.deepcopy(self.executor_identity),batch_fingerprint=self.batch_fingerprint,
                          parent_limits={'max_parameter_vectors':673,'max_native_programs':24000,'workers':12},
                          max_parameter_vectors=1,max_native_programs=56,active_limits=LIMITS,
                          native_count_definition='Invocation attempts at unchanged helper boundary; pre-copy/guard failures excluded; observed returncodes separately retained.',
                          native_returncodes_observed=self.ledger['native_returncodes_observed'])
            return answer
    return FiniteEvaluator


def component_report(cfg,parent,answer,origin):
    dec=read(cfg['decomposition']);manifest=read(parent['manifest']);static={k:v['energy_kcal_mol'] for k,v in manifest['static'].items()}
    def one(cells):
        energies=dict(static,**{k:v['energy'] for k,v in cells.items()});parents={}
        for name,row in dec['parents'].items():
            parents[name]={k:math.fsum(c*energies[cid] for cid,c in terms.items()) for k,terms in row['model_terms'].items()}
        paired={}
        for n in (8,9):
            reference=dec['paired_Dy_minus_La']['W'+str(n)]['reference_total_kcal_mol']
            model={k:parents[f'aquo{n}_central_Dy'][k]-parents[f'aquo{n}_central_La'][k] for k in ('C','P','B')}
            paired['W'+str(n)]={'model':model,'reference':reference,'error':{k:model[k]-reference[k] for k in model}}
        return {'parents':parents,'paired':paired,'J_P':sum(p['error']['P']**2 for p in paired.values())/2,'J_B':sum(p['error']['B']**2 for p in paired.values())/2}
    old,new=one(origin['cells']),one(answer['cells'])
    diagnostic=read(cfg['component_diagnostic']);cert=read(cfg['half_review']);proposal=read(cfg['failed_proposal'])
    row=diagnostic['origins'][-1];require(row['key']==cert['origin'],'Saved prediction origin')
    e0=diagnostic['vectors'][cert['origin']]['component_errors_kcal_mol'];direction=cert['normalized_half_direction']
    predicted={name:e0[name]+math.fsum(a*b for a,b in zip(j,direction)) for name,j in zip(diagnostic['component_order'],row['Jacobian_components'])}
    actual={f'{k}{n}':new['paired']['W'+str(n)]['error'][k] for n in (8,9) for k in ('P','B')}
    oldblocks=origin['score']['block_objectives'];newblocks=answer['score']['block_objectives']
    return {'definition':'Literal paired Dy-La C/P/B; no absolute-parent total-energy target.', 'origin':old,'candidate':new,
            'predicted_component_errors_kcal_mol':predicted,'actual_component_errors_kcal_mol':actual,
            'actual_minus_predicted_component_errors_kcal_mol':{k:actual[k]-predicted[k] for k in actual},
            'predicted_half_blocks':proposal['half']['blocks'],'origin_original_block_scores':oldblocks,'candidate_original_block_scores':newblocks,
            'original_block_changes':{k:newblocks[k]-oldblocks[k] for k in oldblocks},
            'strict_component_score_improvement':{k:new[k]<old[k] for k in ('J_P','J_B')},
            'original_block_exact_nonincrease':{k:newblocks[k]<=oldblocks[k] for k in oldblocks},
            'original25_pass':answer['score']['original_EF_checks_passed'],'physical_threshold_added':False}


def main(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg=json.loads(Path(config_path).read_text());out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),'configuration_fingerprint':fingerprint(cfg),'error':None,'job_id':None,'native_programs':0,'reserved_outcomes_opened':False,'physical_accuracy_qualified':False}
    ev=None
    try:
        base,parent,cert,origin=context(cfg)
        require(mode=='prepare' or str(out)==cfg['output_directory'],'Exact released output directory')
        cls=evaluator_class(base);ev=cls(cfg,out/'cache',mode)
        require(set(ev.m['rows'])==set(origin['cells']) and len(ev.m['rows'])==46 and sum(r['needs_gradient'] for r in ev.m['rows'].values())==10 and sum(r['needs_export'] for r in ev.m['rows'].values())==4,'Exact finite row/program/export denominator')
        for metal in ('La','Dy'):require(base.effective(cfg['candidate']['parameters'],metal,ev.domain)!=base.effective(ev.initial_vector,metal,ev.domain),'Both candidate metals changed')
        # Origin reuse is a saved immutable result, never a disk-cache _parse replay.
        old=ev._assemble(copy.deepcopy(origin['cells']))
        require(old['model']==origin['model'] and old['response']==origin['response'],'Exact immutable origin model/response')
        require(max(abs(a-b) for a,b in zip(old['score']['weighted_residuals'],origin['score']['weighted_residuals']))<=1e-10,'Immutable origin score arithmetic')
        result.update(parent_configuration=cfg['parent_configuration'],origin=cfg['origin'],half_review=cfg['half_review'],candidate=cfg['candidate'],executor=ev.executor_identity,
                      immutable_origin_original_score=origin['score'],origin_reassembly_does_not_replace_stored_score=True,limits=LIMITS)
        if mode=='prepare':
            require(cfg.get('root_release') is None,'Preparation has no execution release')
            # Verify every exact input and rendered key, without starting cells.
            cacheids={}
            for cid in ev.m['rows']:
                ck,identity,key,meta=ev.cache_identity(cfg['candidate']['parameters'],cid)
                cacheids[cid]={'cache_key':ck,'cache_identity':identity,'render':meta,'key_sha256':hashlib.sha256(key.encode()).hexdigest()}
            result.update(status='prepared_fixed_half_native_zero_calls',prepared_cells=cacheids,native_programs=0)
        else:
            result['answer']=ev.evaluate(cfg['candidate']['parameters'])
            result['components']=component_report(cfg,parent,result['answer'],origin)
            result.update(status='complete_one_fixed_half_native_pilot',original25_pass=result['answer']['score']['original_EF_checks_passed'],native_programs=ev.snapshot()['native_program_invocations'])
    except Exception as exc:result.update(status='failed_fixed_half_native_pilot',error=repr(exc))
    if ev is not None:result['snapshot']=ev.snapshot();result['native_programs']=result['snapshot']['native_program_invocations']
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'native_child_user_CPU_seconds':after.ru_utime-child.ru_utime,'native_child_system_CPU_seconds':after.ru_stime-child.ru_stime,'scheduler_allocated_CPU_seconds':0,'native_worker_limit':12,'controller_numerical_threads':1,'aggregate_login_budget':16,'per_receipt_child_CPU_may_overlap_do_not_sum':True}
    save(out/'RESULT.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    r=main(a.mode,a.config,a.out);print(json.dumps({k:r.get(k) for k in ('status','error','native_programs','cost')},indent=2))
    raise SystemExit(0 if r['status'] in ('prepared_fixed_half_native_zero_calls','complete_one_fixed_half_native_pilot') else 1)
