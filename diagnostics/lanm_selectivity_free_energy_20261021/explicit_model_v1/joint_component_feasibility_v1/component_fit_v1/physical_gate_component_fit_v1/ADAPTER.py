#!/usr/bin/env python3
"""New bounded execution policy; inherited scientific rendering/parsing/assembly."""
import concurrent.futures
import copy
import hashlib
import importlib.util
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import shutil
import socket
import threading

_rs=importlib.util.spec_from_file_location("physical_fit_reuse",Path(__file__).with_name("REUSE.py"))
R=importlib.util.module_from_spec(_rs);_rs.loader.exec_module(R)

_js=importlib.util.spec_from_file_location("component_ledger_journal",Path(__file__).with_name("JOURNAL.py"))
J=importlib.util.module_from_spec(_js);_js.loader.exec_module(J)

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
LIMITS={'outer_steps':8,'new_parameter_vectors':200,'new_native_programs':6496,'reused_origin_vectors':1,'reused_initial_derivative_vectors':24,'workers':12,'numerical_threads_per_process':1,'aggregate_login_threads':16,'native_timeout_seconds':300}

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def save(p,x):
    J.atomic_json(p,x)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def fingerprint(c):return digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation')})
def cell_digest(r):return digest({k:v for k,v in r.items() if k!='request_provenance'})


def context(cfg):
    require(cfg['adapter']==pin(__file__) and cfg['limits']==LIMITS,'New policy source/limits')
    require(cfg['journal_implementation']==pin(Path(__file__).with_name('JOURNAL.py')) and cfg['persistence_policy']==J.POLICY,'Versioned journal source/policy')
    for k in ('implementation','core','declaration','execution_context','parent_evaluator','parent_configuration','reuse_inventory','decomposition','interpreter'):verify(cfg[k])
    require(pin(sys.executable)==cfg['interpreter'],'Pinned interpreter')
    env=read(cfg['solver_environment'])
    require(str(Path(sys.prefix).resolve())==env['prefix'],'Pinned Clarabel environment')
    for name,version in env['versions'].items():require(importlib.metadata.version(name)==version,'Pinned numerical package '+name)
    for value in env['files']:verify(value)
    require(all(os.environ.get(k)=='1' for k in THREADS) and os.environ.get('BLIS_NUM_THREADS','1')=='1','One numerical thread')
    declaration=read(cfg['declaration']);require(declaration['limits']==LIMITS,'Frozen declaration limits')
    execution=read(cfg['execution_context']);require(execution['campaign_declaration']==cfg['declaration'] and execution['execution']['mode']=='login' and execution['execution']['workers']==12 and execution['execution']['hostname']==cfg['hostname'] and execution['execution']['aggregate_login_limit']==16 and execution['execution']['job_id'] is None,'Declared truthful login envelope')
    parent=read(cfg['parent_configuration']);require(parent['implementation']==cfg['parent_evaluator'],'Unchanged scientific evaluator')
    inv=read(cfg['reuse_inventory']);require(inv['declaration']==cfg['declaration'] and inv['execution_context']==cfg['execution_context'] and inv['parent_configuration']==cfg['parent_configuration'],'Immutable import inventory lineage')
    require(cfg['reuse_implementation']==pin(Path(__file__).with_name('REUSE.py')),'Frozen reuse reader')
    registry=R.build(inv['authorities'],cfg['parent_configuration'],cfg['parent_evaluator'])
    require(registry==inv['registry'],'Exact immutable metadata registry')
    origin=declaration['start'];old=read(inv['authorities']['original_inventory'])
    require(origin['key']==inv['origin_key']==old['origin_key'] and origin['artifact']==old['vectors'][old['origin_key']],'Exact origin24')
    require(inv['derivative_keys']==old['derivative_keys'] and len(inv['derivative_keys'])==24,'Exact saved24 stencil')
    keys=[inv['origin_key'],*inv['derivative_keys'],inv['pilot_key'],inv['full_pilot_key']]
    require(len(set(keys))==27,'One origin,24 probes,two pilots')
    vectors={key:R.load(registry[key]) for key in keys}
    require(inv['pilot_key']==old['pilot_key'] and inv['full_pilot_key']==R.read(inv['authorities']['full_pilot_root'])['candidate']['key'],'Exact pilot identities')
    require(cfg['reserved_outcomes_authorized'] is False,'Consumed scope only')
    return module(cfg['parent_evaluator'],'component_parent'),parent,inv,vectors


class NotAttempted(RuntimeError):pass
class BudgetStop(RuntimeError):pass


def evaluator_class(base):
    class CampaignEvaluator(base.Evaluator):
        def __init__(self,cfg,cache_dir,mode,inv,vectors):
            require(mode in ('prepare','run'),'Mode');require(not Path(cache_dir).exists(),'Fresh cache only; never read old mutable ledger')
            self.batch=copy.deepcopy(cfg);self.mode=mode;self.batch_fingerprint=fingerprint(cfg)
            self.allowed={};self.complete_vectors=copy.deepcopy(vectors);self.imported_vector_keys=set(vectors)
            self.executor_identity={'mode':'login','job_id':None,'hostname':socket.gethostname(),'controller_pid':os.getpid(),'implementation':cfg['adapter'],'driver':cfg['implementation'],'parent_implementation':cfg['parent_evaluator']}
            super().__init__(read(cfg['parent_configuration']),cache_dir,root_release=cfg.get('root_release'))
            self.program_attempts=[];self.cell_attempts={};self.started_cells=set()
            self.registry=inv['registry'];self.imported_vector_keys=set(self.registry)
            self.saved_metal_index={m:{} for m in ('La','Dy')};self.loaded_saved_vectors={}
            for key,entry in sorted(self.registry.items()):
                for metal in self.saved_metal_index:
                    effective=digest(base.effective(entry['parameters'],metal,self.domain))
                    self.saved_metal_index[metal].setdefault(effective,key)
            self.imported_cell_keys=set()
            for key,v in vectors.items():self._import_saved(key,v)
            self.ledger.update(new_policy=cfg['adapter'],limits=LIMITS,configuration_fingerprint=self.batch_fingerprint,imported_origin_vectors=1,imported_derivative_vectors=24,additional_pilot_vectors=2,available_saved_vectors=len(self.registry),initial_imported_unique_cells=len(self.imported_cell_keys),native_returncodes_observed=0)
            self.persistence=J.Journal(self.cache_dir/'EVENTS.jsonl',self.ledger,self.lock,{'configuration_fingerprint':self.batch_fingerprint,'implementation':cfg['adapter'],'journal_implementation':cfg['journal_implementation'],'root_release':cfg.get('root_release')})
            self.persistence.materialize(self.ledger_path,expected=self.ledger)

        def _import_saved(self,key,value=None):
            # The existing RLock protects lazy loading. No historical ledger/replay.
            with self.lock:
                if key in self.loaded_saved_vectors:return self.loaded_saved_vectors[key]
                value=R.load(self.registry[key]) if value is None else copy.deepcopy(value)
                require(value['parameters']==self.registry[key]['parameters'] and digest(value['parameters'])==key,'Exact saved vector')
                require(set(value['cells'])==set(self.m['rows']),'All46 saved cells')
                for cid,record in value['cells'].items():
                    token,identity,_,meta=self.cache_identity(value['parameters'],cid)
                    require(record['status']=='complete' and record['cache_key']==token and record['cache_identity']==identity,'Saved scientific cell identity')
                    rp=record['request_provenance'];require(rp['cache_key']==token and rp['cache_identity']==identity and rp['render']==meta and rp['vector']==value['parameters'],'Saved exact request render')
                    # Multiple admitted executions may have different receipt paths.
                    # Deterministic first authority wins; never compare derived LAPACK bytes.
                    if token not in self.cache:self.cache[token]=copy.deepcopy(record)
                    self.imported_cell_keys.add(token)
                self.loaded_saved_vectors[key]=value
                return value

        def _saved_metal(self,values,metal):
            return self.saved_metal_index[metal].get(digest(base.effective(values,metal,self.domain)))

        def authorize(self,values,role,step):
            values=self.checked_values(values);key=digest(values)
            require(role in ('probe','trial') and 1<=step<=8,'Only finite declared driver roles')
            self.allowed[key]={'parameters':values,'role':role,'outer_step':step}
            return key

        def _stop(self,exc):
            self.stop_event.set()
            with self.lock:
                self.ledger['stopped']=True;self.ledger.setdefault('failure',repr(exc))
                if hasattr(self,'persistence'):
                    self.persistence.append('stop',{'error':repr(exc)})
                    self.persistence.materialize(self.ledger_path)
                else:save(self.ledger_path,self.ledger)

        def _gate(self):
            require(self.mode=='run','Preparation forbids native calls');require(not self.stop_event.is_set(),'Global stop')
            require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful login execution; no fabricated Slurm identity')
            require(socket.gethostname()==self.batch['hostname'],'Declared login host')
            require(all(os.environ.get(k)=='1' for k in THREADS) and os.environ.get('BLIS_NUM_THREADS','1')=='1','One numerical thread')
            if self._gated:return
            c=self.config
            for p in c['executables'].values():verify(p)
            build=read(c['build_receipt']);q=read(c['qualification_review']);verify(c['source_manifest']);verify(q['implementation']);verify(q['result'])
            require(build['status']=='complete' and build['manifest']==c['source_manifest'] and all(build['executables'][k]==v for k,v in c['executables'].items()),'Unchanged native build')
            require(q['status']=='passed_local_response_native_qualification' and all(q[k]==c[k] for k in ('source_manifest','build_receipt','executables')),'Actual native qualification')
            require(all(q[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'All native admission gates')
            review=read(self.batch['source_review']);prep=read(self.batch['preparation']);release=read(self.root_release)
            require(review['status']=='passed_physical_gate_component_fit_source_review' and review['configuration_fingerprint']==self.batch_fingerprint,'Independent campaign source review')
            require(prep['status']=='prepared_physical_gate_component_fit_zero_native' and prep['configuration_fingerprint']==self.batch_fingerprint and prep['native_programs']==0 and prep['real_proposal_calls']==0,'Actual saved preparation')
            require(release['status']=='released_bounded_physical_gate_component_fit' and release['configuration_fingerprint']==self.batch_fingerprint,'Exact campaign release')
            for k in ('implementation','adapter','core','declaration','execution_context','parent_configuration','reuse_inventory','reuse_implementation','journal_implementation','persistence_policy','solver_environment','saved_initial_proposal','source_review','preparation','limits','output_directory'):
                require(release[k]==self.batch[k],'Release binding '+k)
            require(review['implementation']==self.batch['implementation'] and review['adapter']==self.batch['adapter'] and review['core']==self.batch['core'],'Independent exact three-source review')
            require(release['cache_directory']==str(self.cache_dir) and self.cache_dir.parent==Path(self.batch['output_directory']) and release['reserved_outcomes_authorized'] is False,'Fresh exact execution scope')
            self._gated=True

        def _invoke(self,program,folder,record):
            with self.slots:
                with self.lock:
                    if self.stop_event.is_set():raise NotAttempted('Global stop before native invocation')
                    self._gate();require(self.ledger['native_program_invocations']<6496,'Hard6496 actual-attempt cap')
                    attempt={'number':len(self.program_attempts)+1,'cache_key':record['cache_key'],'cell_id':record['cell_id'],'program':program,'folder':str(folder),'status':'invoke_boundary_entered'}
                    self.program_attempts.append(attempt);self.ledger['native_program_invocations']+=1
                    self.cell_attempts[record['cache_key']]=self.cell_attempts.get(record['cache_key'],0)+1
                    self.ledger['native_attempts']=self.program_attempts
                    self.persistence.append('attempt_intent',{'attempt':attempt,'request_index':len(self.ledger['requests'])-1,'program_inputs':record['program_inputs'][program]})
                args=['system.xyz','EM'] if program=='analyze' else ['system.xyz','Y','N'];before=len(record['receipts'])
                try:
                    self.h['native_helper'].invoke([self.config['executables'][program]['path']]+args,folder,folder/(program+'.out'),record['receipts'])
                    require(len(record['receipts'])==before+1 and record['receipts'][-1].get('returncode')==0,'One successful native receipt')
                    with self.lock:attempt['status']='returned_success'
                except Exception as exc:
                    self.stop_event.set()
                    with self.lock:attempt.update(status='failed_native_attempt',error=repr(exc))
                    self._stop(exc);raise
                finally:
                    with self.lock:
                        for receipt in record['receipts'][before:]:receipt['executor']=copy.deepcopy(self.executor_identity)
                        old_observed=attempt.get('returncode_observed',False)
                        attempt['returncode_observed']=any('returncode' in r for r in record['receipts'][before:])
                        self.ledger['native_returncodes_observed']+=int(attempt['returncode_observed'])-int(old_observed)
                        self.ledger['native_attempts']=self.program_attempts
                        self.persistence.append('attempt_return',{'attempt':attempt,'receipt_delta':record['receipts'][before:],'program_inputs':record['program_inputs'][program]})

        def _cell(self,values,cid):
            if self.stop_event.is_set():raise NotAttempted('Queued cell untouched after stop')
            token,identity,key,meta=self.cache_identity(values,cid);row=self.m['rows'][cid]
            with self.lock:lock=self.cell_locks.setdefault(token,threading.Lock())
            with lock:
                if self.stop_event.is_set():raise NotAttempted('Queued cell untouched after stop')
                saved_key=self._saved_metal(values,row['metal'])
                if token not in self.cache and saved_key is not None:self._import_saved(saved_key)
                if token in self.cache:record=self.cache[token]
                elif base.effective(values,row['metal'],self.domain)==base.effective(self.initial_vector,row['metal'],self.domain):
                    record=copy.deepcopy(self.m['baseline_cells'][cid]);record.update(cache_identity=identity,request_render=meta,cache_key=token);self.cache[token]=record
                else:
                    directory=self.cache_dir/'cells'/token;rp=directory/'RESULT.json'
                    require(not directory.exists(),'No failed-cell retry or mutable disk-cache import')
                    with self.lock:self._gate()
                    directory.mkdir(parents=True,exist_ok=False)
                    record={'status':'started','cache_key':token,'cache_identity':identity,'cell_id':cid,'first_requested_vector':values[:],
                            'inherited_native':False,'implementation':self.config['implementation'],'render':meta,'receipts':[],'program_inputs':{},
                            'job_id':None,'executor':copy.deepcopy(self.executor_identity),'root_release':self.root_release}
                    with self.lock:self.started_cells.add(token)
                    try:
                        for program in (('analyze','testgrad') if row['needs_gradient'] else ('analyze',)):
                            if self.stop_event.is_set():raise NotAttempted('Stop before next native program')
                            folder=directory/program;folder.mkdir();text=key if program=='analyze' else self.render(values,cid,exports=False)[0]
                            for k,n in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/n)
                            (folder/'system.key').write_text(text)
                            record['program_inputs'][program]={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                            save(rp,record);self._invoke(program,folder,record)
                            parsed=self._parse(record,row);record.update(parsed);save(rp,record)
                        require(not row['needs_gradient'] or 'model' in record,'Ten admitted analytic-force rows')
                        record['status']='complete';self.cache[token]=record
                    except Exception as exc:
                        self._stop(exc);record.update(status='not_attempted' if isinstance(exc,NotAttempted) and not record['receipts'] else 'failed',error=repr(exc));raise
                    finally:save(rp,record)
                result=copy.deepcopy(record);result['request_provenance']={'vector':values[:],'cache_key':token,'cache_identity':identity,'render':meta,'executor':copy.deepcopy(self.executor_identity)}
                return result

        def required_programs(self,values):
            total=0
            for cid,r in self.m['rows'].items():
                token,_,_,_=self.cache_identity(values,cid)
                if token not in self.cache and self._saved_metal(values,r['metal']) is None and base.effective(values,r['metal'],self.domain)!=base.effective(self.initial_vector,r['metal'],self.domain):total+=1+int(r['needs_gradient'])
            return total

        def evaluate(self,values):
            values=self.checked_values(values);key=digest(values)
            require(key in self.allowed and self.allowed[key]['parameters']==values,'Driver must authorize explicit probe/trial role')
            with self.evaluation_lock:
                self._gate()
                if key in self.complete_vectors:return copy.deepcopy(self.complete_vectors[key])
                if key in self.registry:
                    value=self._import_saved(key);self.complete_vectors[key]=copy.deepcopy(value);return copy.deepcopy(value)
                needed=self.required_programs(values)
                with self.lock:
                    if len(self.ledger['vector_keys'])>=200 or self.ledger['native_program_invocations']+needed>6496:raise BudgetStop('Whole-vector hard budget unavailable')
                    require(key not in self.ledger['vector_keys'],'No failed-vector retry')
                    self.ledger['vector_keys'].append(key)
                    request={'values':values[:],'vector_key':key,'role':copy.deepcopy(self.allowed[key]),'status':'started','required_new_programs':needed,'cells':{cid:{'status':'not_attempted'} for cid in self.m['rows']}}
                    self.ledger['requests'].append(request)
                    request_index=len(self.ledger['requests'])-1
                    self.persistence.append('request_start',{'index':request_index,'request':request})
                cells={};first_error=None
                try:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
                        futures={pool.submit(self._cell,values,cid):cid for cid in self.m['rows']}
                        for future in concurrent.futures.as_completed(futures):
                            cid=futures[future];token,_,_,_=self.cache_identity(values,cid)
                            try:
                                value=future.result();cells[cid]=value;state={'status':'complete','cache_key':token,'immutable_import':token in self.imported_cell_keys}
                            except Exception as exc:
                                attempts=self.cell_attempts.get(token,0)
                                state={'status':'failed_partial' if attempts else ('not_attempted' if self.stop_event.is_set() else 'failed_before_native'),'error':repr(exc),'native_attempts':attempts}
                                if first_error is None or isinstance(first_error,NotAttempted):first_error=exc
                                self._stop(exc)
                            with self.lock:
                                request['cells'][cid]=state
                                cell_path=self.cache_dir/'cells'/token/'RESULT.json'
                                evidence={'result':pin(cell_path)} if cell_path.exists() else {'cache_key':token,'parent_configuration':self.batch['parent_configuration'],'reuse_inventory':self.batch['reuse_inventory'],'immutable_import':token in self.imported_cell_keys}
                                self.persistence.append('request_cell',{'index':request_index,'cell_id':cid,'state':state,'evidence':evidence})
                    if first_error is not None:raise first_error
                    answer=self._assemble(cells);answer.update(status='complete',parameters=values[:])
                    with self.lock:
                        request.update(status='complete',objective=answer['score']['objective'])
                        self.persistence.append('request_update',{'index':request_index,'updates':{'status':'complete','objective':answer['score']['objective']}})
                        self.persistence.materialize(self.ledger_path,expected=self.ledger)
                    self.complete_vectors[key]=copy.deepcopy(answer);return answer
                except Exception as exc:
                    self._stop(exc)
                    with self.lock:
                        request.update(status='failed',error=repr(exc))
                        self.persistence.append('request_update',{'index':request_index,'updates':{'status':'failed','error':repr(exc)}})
                        self.persistence.materialize(self.ledger_path,expected=self.ledger)
                    raise

        def snapshot(self):
            with self.lock:self.persistence.materialize(self.ledger_path,expected=self.ledger)
            s=super().snapshot()
            with self.lock:s.update(max_parameter_vectors=200,max_native_programs=6496,active_limits=LIMITS,new_parameter_vectors=len(self.ledger['vector_keys']),native_returncodes_observed=self.ledger['native_returncodes_observed'],imported_unique_cells=len(self.imported_cell_keys),imported_origin_vectors=1,imported_derivative_vectors=24,additional_pilot_vectors=2,available_saved_vectors=len(self.registry),loaded_saved_vectors=len(self.loaded_saved_vectors),executor=copy.deepcopy(self.executor_identity),configuration_fingerprint=self.batch_fingerprint)
            return s
    return CampaignEvaluator
