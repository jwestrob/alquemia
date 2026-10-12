"""Owned source-boundary drain. Only the explicit drain operation sends signals."""
import argparse, fcntl, hashlib, json, os, signal, socket, time, traceback
from pathlib import Path

def read(p): return json.loads(Path(p).read_text())
def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    if pin(p['path'])!=p: raise ValueError('changed artifact '+p['path'])
    return Path(p['path'])
def write(p,d):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); temp=p.with_name(p.name+'.tmp.'+str(os.getpid()))
    with temp.open('x') as f: json.dump(d,f,indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    try: os.link(temp,p)
    finally: temp.unlink()
def process(pid):
    p=Path('/proc')/str(pid)
    try:
        a=(p/'stat').read_text().rsplit(')',1)[1].split(); status={s.split(':',1)[0]:s.split(':',1)[1].strip() for s in (p/'status').read_text().splitlines() if ':' in s}
        return {'pid':pid,'start_ticks':a[19],'state':a[0],'ppid':int(a[1]),'pgid':int(a[2]),
                'user_cpu_ticks':int(a[11]),'system_cpu_ticks':int(a[12]),'clock_ticks_per_second':os.sysconf('SC_CLK_TCK'),
                'process_cpu_seconds':(int(a[11])+int(a[12]))/os.sysconf('SC_CLK_TCK'),
                'VmHWM':status.get('VmHWM'),'VmRSS':status.get('VmRSS'),'uid':int(status['Uid'].split()[0]),
                'command':(p/'cmdline').read_bytes().decode().strip('\0').split('\0'),'observed_epoch':time.time()}
    except (FileNotFoundError,ProcessLookupError): return None

def check(t,missing=False):
    a=process(t['pid'])
    if a is None or a['state']=='Z':
        if missing: return None
        raise RuntimeError('owned process no longer live: '+str(t['pid']))
    if (a['start_ticks']!=t['start_ticks'] or a['command']!=t['command'] or a['pgid']!=t['pgid'] or a['uid']!=os.geteuid()):
        raise RuntimeError('process ownership mismatch: '+str(t['pid']))
    return a

def source_state(batch,m,lane):
    active=[]; done=[]; optimizers=[]; dispatched=[]; started=[]
    for c in m['cases']:
        if c['lane']!=lane: continue
        d=batch/'sources'/c['case_id']; terminal=d/'SOURCE_RESULT.json'
        if terminal.exists(): done.append(c['case_id']); continue
        if not d.exists(): continue
        active.append(c['case_id'])
        for p in d.glob('proposals/*/OPTIMIZER_STARTED.json'):
            if not (p.parent/'result.json').exists() and not (p.parent/'INTERRUPTED.json').exists(): optimizers.append({'task_id':p.parent.name,'case_id':c['case_id'],'path':str(p)})
        for p in list(d.glob('origins/*/DISPATCH_STARTED.json'))+list(d.glob('cross_pool/*/DISPATCH_STARTED.json')):
            if not (p.parent/'mace_result.json').exists(): dispatched.append({'path':str(p),'case_id':c['case_id'],'cpu_started':(p.parent/'CPU_CALL_STARTED.json').exists()})
        for p in d.glob('**/CPU_CALL_STARTED.json'):
            if not (p.parent/'mace_result.json').exists(): started.append({'path':str(p),'case_id':c['case_id']})
    return {'active_sources':active,'terminal_sources':done,'incomplete_optimizers':optimizers,'incomplete_dispatches':dispatched,'incomplete_started_calls':started}

def action(s):
    if s['incomplete_optimizers']: return 'finish_source'
    if s['incomplete_started_calls']: return 'wait_call'
    if s['incomplete_dispatches']: return 'flush_dispatch'
    return 'safe'

class Supervisor:
    def __init__(self,batch,output,signals=False):
        self.batch=Path(batch).resolve(); self.output=Path(output).resolve(); self.signals=signals; self.stopped={}; self.events=[]
        self.m=read(self.batch/'manifest.json'); self.manifest=pin(self.batch/'manifest.json'); launch=read(self.batch/'LAUNCH.json')
        if launch['manifest']!=self.manifest: raise ValueError('launch manifest differs')
        verify(self.m['implementation']); self.attempt=self.batch/'attempts/001'
        self.coordinator={'pid':launch['worker_pid'],'start_ticks':launch['worker_start_ticks'],'pgid':launch['worker_pid'],'command':launch['worker_command']}
        self.lanes=[]
        for i in (0,1):
            p=read(self.attempt/('lane'+str(i))/'LANE_START.json'); w=read(self.attempt/('lane'+str(i))/'CPU_WORKER_START.json')
            self.lanes.append({'lane':i,'parent':p,'model':w,'phase':'running','target':None})
        self.targets=[self.coordinator]+[t[k] for t in self.lanes for k in ('parent','model')]
    def event(self,name,**data):
        r={'event':name,'epoch':time.time(),**data}; self.events.append(r)
        if self.signals:
            p=self.output.with_name('DRAIN_EVENTS.jsonl'); p.parent.mkdir(parents=True,exist_ok=True)
            with p.open('a') as f: f.write(json.dumps(r)+'\n'); f.flush(); os.fsync(f.fileno())
        print(json.dumps(r),flush=True)
    def send(self,t,sig):
        if not self.signals: raise RuntimeError('signals disabled in read-only operation')
        a=check(t,missing=True)
        if a is None: return
        os.kill(t['pid'],sig); self.event('signal',pid=t['pid'],signal=signal.Signals(sig).name)
        if sig==signal.SIGSTOP: self.stopped[t['pid']]=t
        if sig==signal.SIGCONT: self.stopped.pop(t['pid'],None)
    def stop(self,t):
        self.send(t,signal.SIGSTOP)
        while True:
            a=check(t,missing=True)
            if a is None or a['state'] in ('T','t'): return
            time.sleep(.01)
    def preflight(self):
        records=[check(t) for t in self.targets]
        for lane in self.lanes:
            if check(lane['parent'])['ppid']!=self.coordinator['pid'] or check(lane['model'])['ppid']!=lane['parent']['pid']: raise RuntimeError('owned parent-child link differs')
        return {'status':'read_only_preflight_pass','manifest':self.manifest,'host':socket.gethostname(),'implementation':pin(__file__),
                'processes':records,'lanes':[{'lane':l['lane'],**source_state(self.batch,self.m,l['lane'])}for l in self.lanes],
                'signals_sent':0,'molecular_calls':0}
    def thaw(self):
        errors=[]
        # Parents first, then coordinator. Both parents are safe before any TERM.
        ts=[t for pid,t in self.stopped.items() if pid!=self.coordinator['pid']]
        if self.coordinator['pid'] in self.stopped: ts.append(self.coordinator)
        for t in ts:
            try: self.send(t,signal.SIGCONT)
            except Exception as e: errors.append(str(e))
        return errors
    def run(self):
        self.output.parent.mkdir(parents=True,exist_ok=True)
        if self.output.exists(): raise RuntimeError('drain receipt already exists')
        lock=self.output.with_name('DRAIN.lock').open('a+'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        before=self.preflight(); write(self.output.with_name('DRAIN_STARTED.json'),before)
        old={s:signal.getsignal(s) for s in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP)}
        def interrupted(s,f): raise KeyboardInterrupt('supervisor signal '+str(s))
        for s in old: signal.signal(s,interrupted)
        error=None; accounting=[]; status='failed'; final=None
        try:
            self.stop(self.coordinator)
            for l in self.lanes:
                s=source_state(self.batch,self.m,l['lane']); l['target']=s['active_sources'][0] if s['active_sources'] else None
                self.event('source_boundary_target',lane=l['lane'],case_id=l['target'])
            while not all(l['phase']=='safe' for l in self.lanes):
                for l in self.lanes:
                    if l['phase']=='safe': continue
                    parent=check(l['parent'],missing=True)
                    if parent is None:
                        if check(l['model'],missing=True) is not None: raise RuntimeError('parent exited with model still live')
                        if action(source_state(self.batch,self.m,l['lane']))!='safe': raise RuntimeError('parent exited with incomplete work')
                        l['phase']='safe'; continue
                    s=source_state(self.batch,self.m,l['lane'])
                    if l['phase']=='running':
                        if l['target'] is not None and l['target'] not in s['terminal_sources']: continue
                        if l['target'] is not None:
                            r=read(self.batch/'sources'/l['target']/'SOURCE_RESULT.json')
                            if r['batch_manifest']!=self.manifest: raise ValueError('terminal source manifest mismatch')
                        self.stop(l['parent']); l['phase']='paused'; s=source_state(self.batch,self.m,l['lane'])
                    if l['phase']=='flushing':
                        # Stop again once that dispatched request actually begins or ends.
                        if not s['incomplete_started_calls'] and s['incomplete_dispatches'] and not s['incomplete_optimizers']: continue
                        self.stop(l['parent']); l['phase']='paused'; s=source_state(self.batch,self.m,l['lane'])
                    decision=action(s)
                    if decision=='finish_source':
                        l['target']=s['incomplete_optimizers'][0]['case_id']; l['phase']='running'
                        self.event('raced_optimizer_preserved',lane=l['lane'],case_id=l['target']); self.send(l['parent'],signal.SIGCONT)
                    elif decision=='flush_dispatch':
                        l['phase']='flushing'; self.event('allow_pending_dispatch',lane=l['lane'],dispatches=s['incomplete_dispatches']); self.send(l['parent'],signal.SIGCONT)
                    elif decision=='safe':
                        l['phase']='safe'; self.event('lane_safe_paused',lane=l['lane'],state=s)
                    # Otherwise parent stays STOP while its separate model finishes.
                time.sleep(.05)
            # Both parent processes are paused with no unfinished optimizer/call.
            for l in self.lanes:
                s=source_state(self.batch,self.m,l['lane'])
                if action(s)!='safe': raise RuntimeError('safe lane changed before shutdown')
                accounting.append({'lane':l['lane'],'parent':check(l['parent'],missing=True),'model':check(l['model'],missing=True),'state':s})
            write(self.output.with_name('DRAIN_PRE_TERMINATION.json'),{'manifest':self.manifest,'accounting':accounting,'host':socket.gethostname(),'epoch':time.time()})
            for l in self.lanes: self.send(l['parent'],signal.SIGTERM)
            for l in self.lanes: self.send(l['parent'],signal.SIGCONT)
            while any(check(l['parent'],missing=True) is not None for l in self.lanes): time.sleep(.1)
            self.send(self.coordinator,signal.SIGCONT)
            while any(check(t,missing=True) is not None for t in self.targets): time.sleep(.1)
            final=self.batch/'FINAL_COLLECTION.json'
            if not final.exists() or read(final)['manifest']!=self.manifest: raise RuntimeError('coordinator final collection missing/mismatched')
            if any(action(source_state(self.batch,self.m,l['lane']))!='safe' for l in self.lanes): raise RuntimeError('incomplete scientific work after drain')
            status='drained'
        except BaseException as exc:
            error=f'{type(exc).__name__}: {exc}'; self.event('exception',error=error,traceback=traceback.format_exc())
        finally:
            cleanup=self.thaw()
            for s,h in old.items(): signal.signal(s,h)
            live=[check(t,missing=True) for t in self.targets]; live=[p for p in live if p is not None]
            stopped=[p for p in live if p['state'] in ('T','t')]
            incomplete=[o['task_id'] for l in self.lanes for o in source_state(self.batch,self.m,l['lane'])['incomplete_optimizers']]
            receipt={'status':status,'manifest':self.manifest,'host':socket.gethostname(),'implementation':pin(__file__),'error':error,'cleanup_errors':cleanup,
                     'live_owned_processes':live,'stopped_owned_processes':stopped,'incomplete_optimizer_task_ids':incomplete,
                     'pre_termination_accounting':accounting,'events':self.events,'batch_collection':pin(final) if final and final.exists() else None,
                     'molecular_calls_by_supervisor':0,'frozen_code_and_parameters_changed':False,
                     'operational_stop_not_scientific_failure':status=='drained'}
            if live or stopped or incomplete or cleanup: receipt['status']='failed'
            write(self.output,receipt); lock.close()
        return receipt

def selftest():
    # Pure decision tests; no target process is signalled.
    base={'incomplete_optimizers':[],'incomplete_dispatches':[],'incomplete_started_calls':[]}
    assert action(base)=='safe'
    assert action({**base,'incomplete_started_calls':[{}]})=='wait_call'
    assert action({**base,'incomplete_dispatches':[{}]})=='flush_dispatch'
    assert action({**base,'incomplete_optimizers':[{}],'incomplete_started_calls':[{}]})=='finish_source'
    own=process(os.getpid()); bad={**own,'start_ticks':str(int(own['start_ticks'])+1)}
    try: check(bad)
    except RuntimeError: pass
    else: raise AssertionError('PID reuse identity not rejected')
    # Exercise emergency thaw ordering using a recorder, with zero os.kill calls.
    s=object.__new__(Supervisor); s.coordinator={'pid':10}; s.stopped={10:s.coordinator,20:{'pid':20},30:{'pid':30}}; sent=[]
    s.send=lambda t,sig:sent.append((t['pid'],sig))
    assert not s.thaw() and sent==[(20,signal.SIGCONT),(30,signal.SIGCONT),(10,signal.SIGCONT)]
    return {'status':'six_file_only_checks_pass','signals_sent':0,'molecular_calls':0}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('operation',choices=('preflight','drain','selftest')); p.add_argument('--batch');p.add_argument('--output')
    a=p.parse_args()
    if a.operation=='selftest': result=selftest()
    else:
        s=Supervisor(a.batch,a.output,signals=a.operation=='drain'); result=s.preflight() if a.operation=='preflight' else s.run()
    print(json.dumps(result,indent=2),flush=True)
    if result['status']=='failed': raise SystemExit(1)
