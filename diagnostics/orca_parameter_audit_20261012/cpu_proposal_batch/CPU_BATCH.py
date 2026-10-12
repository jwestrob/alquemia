"""Finite two-lane CPU transport for unchanged historical OMOL proposal kernels."""
from __future__ import annotations
import argparse, copy, fcntl, hashlib, json, os, signal, subprocess, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/orca_parameter_audit_20261012/cpu_proposal_batch'
ARCHIVE=ROOT/'workspaces/motion_envelope_20260923/pilot34_searches_v1'
PREP=ROOT/'workspaces/orca_parameter_audit_20261012/ca_conditioned_preparation_v3'
SCOUT=ROOT/'workspaces/orca_parameter_audit_20261012/cpu_proposal_scout_v1'
ORDER=['origin','adaptive_Ca','adaptive_La']; METALS=('Ca','La')
PREFIX='CPU_OMOL_JSON '; EV_TO_KCAL=23.06054783061903
THREAD_KEYS=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']

def pin(path):
    p=Path(path).resolve(); h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p: raise ValueError('changed pin: '+p['path'])
    return Path(p['path'])
def read(path): return json.loads(Path(path).read_text())
def write(path,data):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+'.tmp.'+str(os.getpid()))
    with temp.open('x') as f: json.dump(data,f,indent=2,allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    try: os.link(temp,p)
    finally: temp.unlink()
def write_or_check(path,data):
    if Path(path).exists():
        if read(path)!=data: raise ValueError('existing artifact differs: '+str(path))
    else: write(path,data)
def xyz(path):
    a=Path(path).read_text().splitlines(); r=[s.split() for s in a[2:] if s.strip()]
    if len(r)!=int(a[0]): raise ValueError('invalid XYZ')
    return [(v[0],*map(float,v[1:4])) for v in r]
def parent_env():
    if any(os.environ.get(k)!='1' for k in THREAD_KEYS): raise ValueError('parent numerical thread variables must all equal1 before imports')
def modules(m):
    sys.path.insert(0,m['source_snapshot'])
def identity(pid):
    try:
        s=Path('/proc',str(pid),'stat').read_text(); return {'pid':pid,'start_ticks':s.rsplit(')',1)[1].split()[19]}
    except (FileNotFoundError,ProcessLookupError): return None
def live(record): return identity(record['pid'])=={k:record[k] for k in ('pid','start_ticks')}
def stopped(signum,frame): raise KeyboardInterrupt('owned batch received signal '+str(signum))
def handlers():
    old={s:signal.getsignal(s) for s in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT)}
    for s in old: signal.signal(s,stopped)
    return old

def prepare(output):
    import shutil
    out=Path(output).resolve(); out.mkdir(parents=True,exist_ok=False)
    source=read(ARCHIVE/'manifest.json'); inputs=read(PREP/'ORIGIN_INPUTS.json'); scout=read(SCOUT/'manifest.json')
    scope=ROOT/'workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json'
    cases=[]
    for t in inputs['tasks']:
        if t['metal']=='Ca':
            cases.append({'case_id':t['case_id'],'source_case_id':t['source_case_id'],'selection_id':t['selection_id'],
                          'lane':len(cases)%2,'prepared_admission':t['prepared_admission'],
                          'state_signature':t['state_signature'],'source_preparation':t['context_preparation']})
    fixtures=[f for f in scout['fixtures'] if f['id'] in ('00_Ca_origin','01_La_origin')]
    ft=[]
    for f in fixtures:
        t=next(t for t in scout['tasks'] if t['task_id']==f['task_id'])
        ft.append({k:t[k] for k in ('task_id','case_id','metal','xyz','charge','multiplicity')})
    shutil.copyfile(__file__,out/'CPU_BATCH.py'); shutil.copyfile(HERE/'PLAN.md',out/'PLAN.md')
    m={'protocol_id':'corrected_context_106_Ca_sources_CPU_OMOL_original_four_mode_v1',
       'implementation':pin(out/'CPU_BATCH.py'),'plan':pin(out/'PLAN.md'),'scientific_scope':pin(scope),
       'origin_inputs':pin(PREP/'ORIGIN_INPUTS.json'),'preparation':pin(PREP/'PREPARATION.json'),
       'proposal_contract':pin(HERE.parent/'ca_conditioned_preparation/PROPOSAL_CONTRACT.json'),
       'source_manifest':pin(ARCHIVE/'manifest.json'),'source_modules':source['implementation'],
       'source_snapshot':str(ARCHIVE/'implementation'),'scout_result':pin(SCOUT/'FINAL.json'),'scout_manifest':pin(SCOUT/'manifest.json'),
       'model':source['model'],'settings':source['settings'],'pool_settings':read(HERE.parent/'ca_conditioned_preparation/PROPOSAL_CONTRACT.json')['pool_settings'],
       'optimizer_software':source['optimizer_software'],'optimizer_python':source['cpu_python'],
       'model_python':source['gpu_python'],'calculator_source':scout['calculator_source'],
       'optimizer_executable':pin(source['cpu_python']),'model_executable':pin(source['gpu_python']),
       'tasks':inputs['tasks']+ft,'origin_task_count':len(inputs['tasks']),'cases':cases,'fixtures':fixtures,
       'threads_per_model':7,'lanes':2,'parent_threads':1,'maximum_active_numerical_threads':14,
       'gates':{'energy_kcal_mol':.01,'force_max_eV_A':.001},
       'force_gate_precedent':source['implementation']['mace_omol.py'],
       'finite_scope':{'sources':106,'origin_calls':212,'optimizer_starts':212,'maximum_cross_calls':212,
                       'qualification_calls_per_attempt':4,'additional_call_or_wall_cap':None,'GFN2_calls':0,'DFT_calls':0},
       'restart_policy':'reuse_terminal_sources_and_complete_cells; interrupted_started_search_unavailable_no_restart',
       'authorization':'Jacob2026-10-12<=16aggregate login numerical threads; root allocates two sequential-parent/model lanes7 each.'}
    write(out/'manifest.json',m)
    result=validate(out/'manifest.json'); write(out/'DRY_RUN.json',result); print(json.dumps(result),flush=True)

def validate(manifest):
    parent_env(); m=read(manifest)
    for k in ('implementation','plan','scientific_scope','origin_inputs','preparation','proposal_contract','source_manifest',
              'scout_result','scout_manifest','calculator_source','optimizer_executable','model_executable'): verify(m[k])
    for p in m['source_modules'].values(): verify(p)
    verify(m['model']['checkpoint']); modules(m)
    import numpy as np, scipy
    import scipy.optimize._slsqplib as kernel
    import scipy.optimize._slsqp_py as wrapper
    from mace_site_kinematics import Kinematics
    from accommodation_torsion_profiles import geometry_check
    from affordable_common import paired
    from slsqp_precision import engine
    source=read(verify(m['source_manifest'])); inputs=read(verify(m['origin_inputs'])); contract=read(verify(m['proposal_contract'])); scout=read(verify(m['scout_manifest']))
    if m['gates']!={'energy_kcal_mol':.01,'force_max_eV_A':.001} or m['fixtures']!=[f for f in scout['fixtures'] if f['id'] in ('00_Ca_origin','01_La_origin')]: raise ValueError('frozen qualification differs')
    if m['source_modules']!=source['implementation'] or m['source_snapshot']!=str(ARCHIVE/'implementation') or m['optimizer_python']!=source['cpu_python'] or m['model_python']!=source['gpu_python']: raise ValueError('frozen runtime differs')
    if (pin(kernel.__file__)!=m['optimizer_software']['kernel'] or pin(wrapper.__file__)!=m['optimizer_software']['wrapper']
        or scipy.__version__!=m['optimizer_software']['version']): raise ValueError('original optimizer binary differs')
    if (m['model']!=source['model'] or m['model']!=inputs['model'] or m['settings']!=engine()[0].completion.SETTINGS
        or m['settings']!=contract['optimizer_settings'] or m['pool_settings']!=contract['pool_settings']): raise ValueError('method differs')
    if (len(m['cases'])!=106 or m['origin_task_count']!=212 or m['tasks'][:212]!=inputs['tasks'] or len(m['tasks'])!=214
        or len({c['case_id'] for c in m['cases']})!=106 or m['threads_per_model']!=7 or m['lanes']!=2
        or m['maximum_active_numerical_threads']!=14): raise ValueError('finite scope/resource differs')
    expected=[r['case_id'] for r in read(verify(m['preparation']))['rows'] if r['status']=='prepared' and not r.get('prepared_artifact_reused')]
    if set(c['source_case_id'] for c in m['cases'])!=set(expected) or [c['case_id'] for c in m['cases']]!=[t['case_id'] for t in inputs['tasks'] if t['metal']=='Ca']: raise ValueError('exact prepared106 population/order differs')
    checks=[]
    for i,c in enumerate(m['cases']):
        ts=[t for t in m['tasks'][:212] if t['case_id']==c['case_id']]
        if [t['metal'] for t in ts]!=list(METALS) or c['lane']!=i%2: raise ValueError('paired lane/source differs')
        ca,la=ts; paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if read(verify(ca['mapping']))!=read(verify(la['mapping'])): raise ValueError('paired map differs')
        for t in ts:
            for k in ('mapping','state_signature','context_preparation','prepared_admission'): verify(t[k])
            admission=read(verify(t['prepared_admission'])); ep=admission['representations']['context']['endpoints'][t['metal']]
            if any(t[k]!=ep[k] for k in ('xyz','charge','multiplicity')): raise ValueError('prepared state differs')
            kin=Kinematics(read(verify(t['mapping']))['context']); zero=np.zeros(len(kin.modes)); a=xyz(verify(t['xyz']))
            delta=float(np.max(np.abs(kin.evaluate(zero)[1]-np.asarray([v[1:] for v in a]))))
            check=geometry_check(kin,zero,[v[0] for v in a])
            if delta>1e-12 or not check['pass'] or t['multiplicity']!=1: raise ValueError('physical origin failed')
        checks.append({'case_id':c['case_id'],'atom_count':len(a),'modes':len(kin.modes),'q0_max_error_A':delta,'pass':True})
    for f in m['fixtures']:
        for k in ('xyz','forces','native_receipt'): verify(f[k])
    return {'status':'dry_run_pass','manifest':pin(manifest),'new_molecular_calls':0,'model_loaded':False,
            'source_checks':checks,'sources':106,'origin_cells':212,'optimizer_starts':212,'maximum_cross_cells':212,
            'qualification_cells_per_attempt':4,'threads_per_model':7,'maximum_active_numerical_threads':14,
            'selector_status':'deferred_until_actual_paired_origin_forces','all_historical_module_pins_verified':len(m['source_modules'])}

def authorize_request(manifest,m,request):
    """Allow only pinned root tasks or actual same-source child manifests."""
    caller=read(verify(request['manifest'])); rootpin=pin(manifest)
    t=next(t for t in caller['tasks'] if t['task_id']==request['task_id'])
    if request['manifest']==rootpin:
        if t not in m['tasks'] or request['xyz']!=t['xyz']: raise ValueError('root request differs')
    else:
        if caller['batch_manifest']!=rootpin or caller['model']!=m['model']: raise ValueError('child manifest not bound to batch')
        base=next(t0 for t0 in m['tasks'][:m['origin_task_count']] if (t0['case_id'],t0['metal'])==(t['case_id'],t['metal']))
        if any(t[k]!=base[k] for k in ('charge','multiplicity')): raise ValueError('child state differs')
        if caller['kind']=='proposals':
            if caller['settings']!=m['settings'] or t['xyz']!=base['xyz'] or t['mapping']!=base['mapping']: raise ValueError('optimizer method differs')
        elif caller['kind']=='cross_pool':
            if request['xyz']!=t['xyz']: raise ValueError('cross geometry differs')
        else: raise ValueError('unknown child manifest kind')
    if (request['charge'],request['multiplicity'])!=(t['charge'],t['multiplicity']): raise ValueError('request state differs')
    a=xyz(verify(request['xyz'])); origin=xyz(verify(t['xyz']))
    if [v[0] for v in a]!=[v[0] for v in origin] or a[0][0]!=t['metal']: raise ValueError('atom order/metal differs')
    return t,a

def worker(manifest,attempt,lane):
    m=read(manifest); verify(m['implementation']); verify(m['calculator_source']); verify(m['model']['checkpoint']); modules(m)
    import resource, importlib.metadata
    import numpy as np, torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    out=Path(attempt)/('lane'+str(lane)); began=time.monotonic(); cpu=time.process_time(); attempts=forwards=failures=0
    torch.set_num_threads(7); torch.set_num_interop_threads(1); torch.set_default_dtype(torch.float64)
    calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cpu',default_dtype='float64'); model=calc.models[0]
    if type(model).__name__!='ScaleShiftMACE' or list(model.heads)!=['omol'] or dict(model.embedding_specs)!=m['model']['embedding_specs'] or calc.energy_units_to_eV!=1.: raise ValueError('native model differs')
    for p in model.parameters(): p.requires_grad_(False)
    versions={n:p._version for n,p in model.named_parameters()}
    def before_forward(*args):
        nonlocal forwards
        forwards+=1
    hook=model.register_forward_pre_hook(before_forward); load=time.monotonic()-began
    software={n:importlib.metadata.version(n) for n in ('torch','mace-torch','numpy','ase','e3nn')}
    print(PREFIX+json.dumps({'status':'ready','load_seconds':load,'threads':torch.get_num_threads(),'versions':software}),flush=True)
    try:
        for line in sys.stdin:
            p=json.loads(line)['request']; rp=verify(p); req=read(rp); attempts+=1
            result={'request':p,'manifest':req['manifest'],'worker_manifest':pin(manifest),'task_id':req['task_id'],'status':'failed',
                    'energy_eV':None,'forces':None,'model_call_started':False,'device':'cpu','model':m['model']}
            tick=time.monotonic(); before=forwards; cpu0=time.process_time()
            try:
                t,rows=authorize_request(manifest,m,req)
                if t.get('case_id') in {c['case_id'] for c in m['cases']} and next(c['lane'] for c in m['cases'] if c['case_id']==t['case_id'])!=lane: raise ValueError('wrong worker lane')
                atoms=Atoms([a[0] for a in rows],positions=[a[1:] for a in rows],pbc=False)
                atoms.info.update(charge=t['charge'],spin=t['multiplicity']); atoms.calc=calc
                state=input_batch(calc,atoms,t['charge'],t['multiplicity'])
                if state['graph_count']!=1: raise ValueError('graph count differs')
                calc.reset(); write(rp.parent/'CPU_CALL_STARTED.json',{'request':p,'worker':identity(os.getpid()),'time_unix':time.time()})
                result['model_call_started']=True
                energy=float(atoms.get_potential_energy()); force=np.asarray(atoms.get_forces(),dtype=float)
                if not np.isfinite(energy) or force.shape!=(len(rows),3) or not np.isfinite(force).all(): raise ValueError('nonfinite energy/force')
                if versions!={n:q._version for n,q in model.named_parameters()}: raise ValueError('model weights mutated')
                fp=rp.parent/'forces_eV_A.npy'; np.save(fp,force,allow_pickle=False)
                result.update(status='complete',energy_eV=energy,forces=pin(fp),state=state,parameter_versions_unchanged=True,
                              force_definition='negative Cartesian total vacuum OMOL energy gradient, eV/angstrom')
            except Exception as exc: failures+=1; result.update(reason=str(exc),traceback=traceback.format_exc())
            result.update(wall_seconds=time.monotonic()-tick,process_cpu_seconds=time.process_time()-cpu0,
                          actual_model_forwards=forwards-before,threads=torch.get_num_threads(),versions=software,
                          peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scheduler_allocation=None)
            dest=rp.parent/'mace_result.json'; write(dest,result)
            print(PREFIX+json.dumps({'result':pin(dest)}),flush=True)
    finally:
        hook.remove(); write(out/'CPU_WORKER_SUMMARY.json',{'manifest':pin(manifest),'attempted_requests':attempts,'actual_model_forwards':forwards,
             'failed_requests':failures,'model_load_seconds':load,'wall_seconds':time.monotonic()-began,'process_cpu_seconds':time.process_time()-cpu,
             'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':torch.get_num_threads(),'versions':software,'device':'cpu','allocated_gpu_seconds':0})

class WarmCPU:
    def __init__(self,manifest,m,attempt,lane):
        self.calls=[]; self.log=(Path(attempt)/('lane'+str(lane))/'CPU_WORKER.log').open('x')
        env=os.environ.copy()
        for k in THREAD_KEYS: env[k]='1'
        for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OMP_THREAD_LIMIT'): env[k]='7'
        cmd=[m['model_python'],'-u',str(verify(m['implementation'])),'worker','--manifest',str(manifest),'--attempt',str(attempt),'--lane',str(lane)]
        self.process=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,start_new_session=True,env=env)
        write(Path(attempt)/('lane'+str(lane))/'CPU_WORKER_START.json',{**identity(self.process.pid),'pgid':os.getpgid(self.process.pid),'command':cmd,'parent_pid':os.getpid(),'thread_environment':{k:env[k] for k in env if 'THREAD' in k}})
        try:
            self.ready=self.receive()
            if self.ready.get('status')!='ready': raise RuntimeError('CPU worker did not initialize')
        except BaseException: self.abort(); raise
    def receive(self):
        while True:
            line=self.process.stdout.readline()
            if not line: raise RuntimeError('CPU worker ended: '+str(self.process.poll()))
            self.log.write(line); self.log.flush()
            if line.startswith(PREFIX): return json.loads(line[len(PREFIX):])
    def evaluate(self,p):
        self.process.stdin.write(json.dumps({'request':p})+'\n'); self.process.stdin.flush()
        rpin=self.receive()['result']; r=read(verify(rpin))
        if r['request']!=p: raise ValueError('worker returned another request')
        self.calls.append(rpin)
        print(json.dumps({'task':r['task_id'],'status':r['status'],'wall_seconds':r['wall_seconds'],'forwards':r['actual_model_forwards']}),flush=True)
        return r,rpin
    def abort(self):
        if self.process.poll() is None:
            os.killpg(self.process.pid,signal.SIGTERM)
            try: self.process.wait(timeout=10)
            except subprocess.TimeoutExpired: os.killpg(self.process.pid,signal.SIGKILL); self.process.wait()
        if not self.log.closed: self.log.close()
    def close(self):
        if self.process.poll() is None:
            self.process.stdin.close()
            for line in self.process.stdout: self.log.write(line)
            code=self.process.wait()
        else: code=self.process.returncode
        self.log.close()
        if code: raise RuntimeError('worker exit '+str(code))

def cell(directory,request,transport,model):
    """Finished calls reused; previously begun incomplete calls never repeated."""
    d=Path(directory); rp=d/'mace_request.json'; result=d/'mace_result.json'; d.mkdir(parents=True,exist_ok=True)
    write_or_check(rp,request)
    if result.exists():
        r=read(result)
        if r['request']!=pin(rp) or r['model']!=model: raise ValueError('cached call identity differs')
        if r['status']=='complete': verify(r['forces'])
        return r,pin(result)
    unavailable=d/'INTERRUPTED.json'
    if unavailable.exists(): return read(unavailable),pin(unavailable)
    if (d/'DISPATCH_STARTED.json').exists() or (d/'CPU_CALL_STARTED.json').exists():
        r={'status':'unavailable','reason':'previously_started_call_without_terminal_receipt_not_repeated','request':pin(rp),'model':model,'energy_eV':None,'forces':None}
        write(unavailable,r); return r,pin(unavailable)
    write(d/'DISPATCH_STARTED.json',{'request':pin(rp),'parent':identity(os.getpid()),'time_unix':time.time()})
    return transport.evaluate(pin(rp))

def proposal(manifest,t,transport,completion):
    d=Path(manifest).parent/'proposals'/t['task_id']; d.mkdir(parents=True,exist_ok=True)
    dest=d/'result.json'; interruption=d/'INTERRUPTED.json'
    if dest.exists():
        r=read(dest)
        if r['manifest']!=pin(manifest): raise ValueError('proposal manifest differs')
        return pin(dest)
    if interruption.exists(): return pin(interruption)
    if (d/'OPTIMIZER_STARTED.json').exists():
        write(interruption,{'status':'unavailable','task_id':t['task_id'],'manifest':pin(manifest),'origin':None,'proposal':None,
                            'reason':'interrupted_started_search_preserved_no_q0_restart','partial_evaluations_retained':True})
        return pin(interruption)
    write(d/'OPTIMIZER_STARTED.json',{'manifest':pin(manifest),'task_id':t['task_id'],'parent':identity(os.getpid()),'time_unix':time.time()})
    return completion.optimize(manifest,t,transport)

def native_cell(t,r,rpin,coordinate,origin):
    if r['status']!='complete': return {'status':'unavailable','reason':r.get('reason',r['status']),'MACE':rpin,'xyz':coordinate,'charge':t['charge'],'multiplicity':t['multiplicity']}
    verify(coordinate); verify(r['forces'])
    req=read(verify(r['request']))
    if req['xyz']!=coordinate or (req['charge'],req['multiplicity'])!=(t['charge'],t['multiplicity']): raise ValueError('native reused cell differs')
    return {'status':'complete','xyz':coordinate,'MACE':rpin,'energy_eV':r['energy_eV'],'components':{'MACE_eV':r['energy_eV']},
            'forces':r['forces'],'charge':t['charge'],'multiplicity':t['multiplicity'],'candidate_origin':origin,
            'composite_status':'unavailable_not_evaluated'}

def run_source(manifest,m,c,transport):
    import motion_envelope_run as builder
    import slsqp_precision
    from nikasha_pool import same_geometry
    from mace_hybrid import write_xyz
    root=Path(manifest).parent; d=root/'sources'/c['case_id']; d.mkdir(parents=True,exist_ok=True)
    terminal=d/'SOURCE_RESULT.json'
    if terminal.exists():
        r=read(terminal)
        if r['batch_manifest']!=pin(manifest): raise ValueError('terminal source belongs to another manifest')
        return pin(terminal)
    ts=[t for t in m['tasks'][:212] if t['case_id']==c['case_id']]; prepared=read(verify(c['prepared_admission']))
    row={'case_id':c['case_id'],'source_case_id':c['source_case_id'],'selection_id':c['selection_id'],'source':prepared,
         'mapping':{t['metal']:t['mapping'] for t in ts},'source_preparation':c['source_preparation'],'native_endpoints':{},'solvent_endpoints':None}
    result={**c,'expected_class':prepared['source']['expected_class'],'biological_group':prepared['source']['biological_group'],'batch_manifest':pin(manifest),'source':{'origin_row':row},'prepared_admission':c['prepared_admission'],
            'candidates':{},'aliases':{},'matrix':{z:{} for z in METALS},'proposals':{},'status':'unavailable',
            'pool':{'status':'unavailable','reason':'composite_solvent_not_evaluated'},'composite_status':'unavailable_not_evaluated'}
    try:
        for t in ts:
            z=t['metal']; req={'manifest':pin(manifest),'task_id':t['task_id'],'xyz':t['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
            r,rp=cell(d/'origins'/z,req,transport,m['model'])
            result['matrix'][z]['origin']=native_cell(t,r,rp,t['xyz'],'origin')
            row['native_endpoints'][z]={'xyz':t['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity'],'status':r['status'],
                                      'native_MACE_receipt':rp,'native_MACE_energy_eV':r.get('energy_eV')}
        result['candidates']['origin']={'status':'available','xyz':ts[0]['xyz'],'origin':'prepared_source'}
        result['aliases']['origin']={'representative':'origin'}
        if not all(ep['status']=='complete' for ep in row['native_endpoints'].values()): raise ValueError('actual paired origin force unavailable')
        tasks,selector=builder.build_tasks(row,m['model']); row['status']='native_complete'
        pm={'kind':'proposals','batch_manifest':pin(manifest),'model':m['model'],'settings':m['settings'],'tasks':tasks,
            'source':row,'selector':selector,'prepared_admission':c['prepared_admission']}
        pp=d/'proposal_manifest.json'; write_or_check(pp,pm); result['proposal_manifest']=pin(pp); result['selection']=selector
        completion=slsqp_precision.engine()[0].completion
        for t in tasks:
            rp=proposal(pp,t,transport,completion); r=read(verify(rp)); name='adaptive_'+t['metal']; result['proposals'][t['metal']]=rp
            if r['status']=='proposal_available':
                point=r['proposal']; result['candidates'][name]={'status':'available','xyz':point['coordinate'],'origin':t['metal'],'proposal_receipt':rp,'boundary_flag':r['boundary_flag']}
                native=read(verify(point['MACE'])); result['matrix'][t['metal']][name]=native_cell(t,native,point['MACE'],point['coordinate'],name)
            else: result['candidates'][name]={'status':'unavailable','reason':r.get('reason',r['status']),'proposal_receipt':rp}
            if transport.process.poll() is not None: raise RuntimeError('model worker stopped')
        for name in ORDER:
            candidate=result['candidates'][name]
            if candidate['status']!='available': continue
            representative=next((earlier for earlier in ORDER[:ORDER.index(name)] if result['candidates'][earlier]['status']=='available' and same_geometry(xyz(verify(candidate['xyz'])),xyz(verify(result['candidates'][earlier]['xyz'])))),name)
            result['aliases'][name]={'representative':representative,'coordinate_tolerance_A':1e-12}
            for t in ts:
                z=t['metal']
                if name in result['matrix'][z]: continue
                if representative!=name and result['matrix'][z].get(representative,{}).get('status')=='complete':
                    result['matrix'][z][name]={**result['matrix'][z][representative],'alias_of':representative,'candidate_origin':name}; continue
                rows=xyz(verify(candidate['xyz'])); rows[0]=(z,*rows[0][1:]); cd=d/'cross_pool'/(name+'__'+z); cd.mkdir(parents=True,exist_ok=True); xp=cd/'context.xyz'
                if not xp.exists(): write_xyz(xp,rows)
                if not same_geometry(xyz(xp),rows) or xyz(xp)[0][0]!=z: raise ValueError('cross coordinates differ')
                ct={'task_id':c['case_id']+'__'+name+'__'+z,'case_id':c['case_id'],'metal':z,'xyz':pin(xp),'charge':t['charge'],'multiplicity':t['multiplicity']}
                cm={'kind':'cross_pool','batch_manifest':pin(manifest),'model':m['model'],'tasks':[ct],'candidate_source':candidate}
                cp=cd/'manifest.json'; write_or_check(cp,cm)
                req={'manifest':pin(cp),'task_id':ct['task_id'],'xyz':ct['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
                r,rp=cell(cd,req,transport,m['model']); result['matrix'][z][name]=native_cell(t,r,rp,ct['xyz'],name)
        result['status']='native_matrix_complete' if all(result['matrix'][z].get(n,{}).get('status')=='complete' for z in METALS for n in ORDER) else 'unavailable'
    except Exception as exc:
        result['reason']=str(exc); result['traceback']=traceback.format_exc()
        if transport.process.poll() is not None: raise
    for name in ORDER:
        result['candidates'].setdefault(name,{'status':'unavailable','reason':result.get('reason','proposal unavailable')})
        for z in METALS: result['matrix'][z].setdefault(name,{'status':'unavailable','reason':result['candidates'][name].get('reason',result.get('reason','cell unavailable'))})
    result['native_pool_status']='available' if result['status']=='native_matrix_complete' else 'unavailable'
    write(terminal,result); print(json.dumps({'source':c['case_id'],'status':result['status']}),flush=True); return pin(terminal)

def lane_run(manifest,attempt,lane):
    import resource
    parent_env(); m=read(manifest); modules(m); out=Path(attempt)/('lane'+str(lane)); old=handlers(); transport=None
    import numpy as np
    started=time.monotonic(); cpu=time.process_time(); error=None; results=[]
    try:
        transport=WarmCPU(manifest,m,attempt,lane); checks=[]
        for f in m['fixtures']:
            t=next(t for t in m['tasks'] if t['task_id']==f['task_id']); td=out/'qualification'/f['id']
            req={'manifest':pin(manifest),'task_id':t['task_id'],'xyz':f['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
            r,rp=cell(td,req,transport,m['model']); check={'fixture':f,'receipt':rp,'pass':False}
            if r['status']=='complete':
                de=(r['energy_eV']-f['energy_eV'])*EV_TO_KCAL; df=float(np.max(np.abs(np.load(verify(r['forces']),allow_pickle=False)-np.load(verify(f['forces']),allow_pickle=False))))
                check.update(energy_error_kcal_mol=de,force_max_error_eV_A=df,**{'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']})
            checks.append(check)
        write(out/'QUALIFICATION.json',{'checks':checks,'pass':all(c['pass'] for c in checks),'threads':7})
        if not all(c['pass'] for c in checks): raise RuntimeError('7thread qualification failed before new sources')
        while not (Path(attempt)/'QUALIFICATION_RELEASE.json').exists(): time.sleep(.2)
        for c in m['cases']:
            if c['lane']==lane: results.append(run_source(manifest,m,c,transport))
    except BaseException as exc: error=f'{type(exc).__name__}: {exc}'; traceback.print_exc()
    finally:
        if transport:
            try:
                if error: transport.abort()
                else: transport.close()
            except BaseException as exc: transport.abort(); error=(error+'; ' if error else '')+str(exc)
        for s,h in old.items(): signal.signal(s,h)
        write(out/'LANE_RESULT.json',{'status':'complete' if error is None else 'failed','error':error,'sources':results,
              'wall_seconds':time.monotonic()-started,'parent_process_cpu_seconds':time.process_time()-cpu,
              'parent_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'model_receipts':transport.calls if transport else []})
    if error: raise SystemExit(1)

def terminate_owned(process):
    if process.poll() is None:
        os.killpg(process.pid,signal.SIGTERM)
        try: process.wait(timeout=20)
        except subprocess.TimeoutExpired: os.killpg(process.pid,signal.SIGKILL); process.wait()

def collect(manifest,attempt,error):
    m=read(manifest); root=Path(manifest).parent; rows=[]; pending=[]
    for c in m['cases']:
        p=root/'sources'/c['case_id']/'SOURCE_RESULT.json'
        if p.exists(): rows.append(read(p))
        else: pending.append(c['case_id'])
    receipts=list(root.glob('sources/*/origins/*/mace_result.json'))+list(root.glob('sources/*/cross_pool/*/mace_result.json'))+list(root.glob('sources/*/proposals/*/evaluations/*/mace_result.json'))+list(root.glob('attempts/*/lane*/qualification/*/mace_result.json'))
    rr=[read(p) for p in receipts]; started=list(root.glob('sources/*/proposals/*/OPTIMIZER_STARTED.json'))
    interrupted_calls=list(root.glob('sources/**/CPU_CALL_STARTED.json'))
    r={'protocol_id':m['protocol_id'],'manifest':pin(manifest),'attempt':str(attempt),'status':'complete' if not pending else 'interrupted_or_unstarted',
       'failure':error,'terminal':not pending,'cases':rows,'pending_source_ids':pending,'denominator':106,'terminal_sources':len(rows),
       'native_matrix_complete':sum(r['status']=='native_matrix_complete' for r in rows),'native_matrix_unavailable':sum(r['status']!='native_matrix_complete' for r in rows),
       'costs':{'recorded_model_requests':len(rr),'recorded_model_forwards':sum(r['actual_model_forwards'] for r in rr),
               'failed_model_requests':sum(r['status']!='complete' for r in rr),'optimizer_starts':len(started),
               'recorded_inference_wall_seconds_sum':sum(r['wall_seconds'] for r in rr),
               'recorded_inference_process_cpu_seconds':sum(r['process_cpu_seconds'] for r in rr),
               'model_calls_started_without_terminal_receipt':[str(p) for p in interrupted_calls if not (p.parent/'mace_result.json').exists()],
               'cost_note':'Interrupted forwards without terminal receipt have unknown completion/cost; worker and lane summaries retain import/other CPU cost.'},
       'source_results':[pin(root/'sources'/c['case_id']/'SOURCE_RESULT.json') for c in rows],
       'new_GFN2_calls':0,'new_DFT_calls':0,'production_changed':False,'composite_status':'unavailable_not_evaluated'}
    dest=Path(attempt)/'FINAL_COLLECTION.json'; write(dest,r)
    # This is the only mutable index; every prior complete snapshot stays in attempts/.
    tmp=root/('FINAL_COLLECTION.json.tmp.'+str(os.getpid())); tmp.write_bytes(dest.read_bytes()); os.replace(tmp,root/'FINAL_COLLECTION.json')
    return r

def execute(manifest):
    validate(manifest); m=read(manifest); root=Path(manifest).parent; lock=(root/'RUN.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    for p in root.glob('attempts/*/lane*/CPU_WORKER_START.json'):
        if live(read(p)): raise RuntimeError('prior owned model worker still alive; inspect before restart: '+str(p))
    for p in root.glob('attempts/*/lane*/LANE_START.json'):
        if live(read(p)): raise RuntimeError('prior owned parent still alive; inspect before restart: '+str(p))
    attempt=root/'attempts'/('%03d'%(len(list((root/'attempts').glob('*')))+1)); attempt.mkdir(parents=True)
    if all((root/'sources'/c['case_id']/'SOURCE_RESULT.json').exists() for c in m['cases']):
        write(attempt/'RUN_RESULT.json',{'manifest':pin(manifest),'status':'terminal_recollection_only','new_molecular_calls':0})
        result=collect(manifest,attempt,None)
        print(json.dumps({'status':result['status'],'new_molecular_calls':0,'collection':pin(root/'FINAL_COLLECTION.json')})); return
    old=handlers(); processes=[]; logs=[]; error=None; began=time.monotonic()
    write(attempt/'RUN_STARTED.json',{'manifest':pin(manifest),'parent':identity(os.getpid()),'time_unix':time.time(),'maximum_active_numerical_threads':14})
    try:
        for lane in (0,1):
            ld=attempt/('lane'+str(lane)); ld.mkdir(); log=(ld/'PARENT.log').open('x'); logs.append(log)
            cmd=[m['optimizer_python'],'-u',str(verify(m['implementation'])),'lane','--manifest',str(manifest),'--attempt',str(attempt),'--lane',str(lane)]
            p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True); processes.append(p)
            write(ld/'LANE_START.json',{**identity(p.pid),'pgid':os.getpgid(p.pid),'command':cmd})
        while True:
            q=[attempt/('lane'+str(i))/'QUALIFICATION.json' for i in (0,1)]
            if all(p.exists() for p in q):
                if not all(read(p)['pass'] for p in q): raise RuntimeError('qualification failed; no new sources admitted')
                write(attempt/'QUALIFICATION_RELEASE.json',{'checks':[pin(p) for p in q],'all_four_pass':True}); break
            if any(p.poll() is not None for p in processes): raise RuntimeError('worker failed before qualification barrier')
            time.sleep(.2)
        while any(p.poll() is None for p in processes):
            if any(p.poll() not in (None,0) for p in processes): raise RuntimeError('lane failed; preserving partial work')
            time.sleep(.5)
        if any(p.returncode!=0 for p in processes): raise RuntimeError('lane failed')
    except BaseException as exc: error=f'{type(exc).__name__}: {exc}'; traceback.print_exc()
    finally:
        for p in processes: terminate_owned(p)
        # A lane hard-killed during cleanup can leave its separate model session.
        for p in attempt.glob('lane*/CPU_WORKER_START.json'):
            owned=read(p)
            if live(owned):
                try: os.killpg(owned['pgid'],signal.SIGTERM)
                except ProcessLookupError: pass
        for f in logs: f.close()
        for s,h in old.items(): signal.signal(s,h)
        write(attempt/'RUN_RESULT.json',{'manifest':pin(manifest),'failure':error,'wall_seconds':time.monotonic()-began,'lane_exit_codes':[p.returncode for p in processes]})
        result=collect(manifest,attempt,error)
    print(json.dumps({k:result[k] for k in ('status','terminal_sources','native_matrix_complete','native_matrix_unavailable','failure')}),flush=True)
    if error: raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('operation',choices=('prepare','validate','execute','lane','worker'))
    p.add_argument('--output'); p.add_argument('--manifest'); p.add_argument('--attempt'); p.add_argument('--lane',type=int); a=p.parse_args()
    if a.operation=='prepare': prepare(a.output)
    elif a.operation=='validate': print(json.dumps(validate(a.manifest)),flush=True)
    elif a.operation=='execute': execute(a.manifest)
    elif a.operation=='lane': lane_run(a.manifest,a.attempt,a.lane)
    else: worker(a.manifest,a.attempt,a.lane)
