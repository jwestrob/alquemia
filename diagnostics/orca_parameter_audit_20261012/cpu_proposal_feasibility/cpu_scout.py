"""Isolated two-process CPU parity scout; historical proposal kernels untouched."""
from __future__ import annotations
import argparse, copy, hashlib, json, os, signal, subprocess, sys, time, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
ARCHIVE=ROOT/'workspaces/motion_envelope_20260923/pilot34_searches_v1'
CASE='a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1'
EV_TO_KCAL=23.06054783061903
PREFIX='CPU_OMOL_JSON '

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])['sha256']!=p['sha256']:raise ValueError('artifact changed: '+p['path'])
    return Path(p['path'])
def read(path):return json.loads(Path(path).read_text())
def write(path,data):
    with Path(path).open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
def xyz(path):
    lines=Path(path).read_text().splitlines();rows=[s.split() for s in lines[2:] if s.strip()]
    if len(rows)!=int(lines[0]):raise ValueError('invalid XYZ')
    return [(r[0],*map(float,r[1:4])) for r in rows]
def modules(manifest):
    for p in manifest['source_modules'].values():verify(p)
    sys.path.insert(0,str(ARCHIVE/'implementation'))

def prepare(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    source=read(ARCHIVE/'manifest.json');audit=read(HERE/'RESULT.json')
    ts=[copy.deepcopy(t) for t in source['tasks'] if t['source_case_id']==CASE]
    if len(ts)!=2 or [t['metal'] for t in ts]!=['Ca','La']:raise ValueError('paired fixture absent')
    fixtures={t['metal']:next(f for f in audit['force_proposal_fixtures'] if f['task_id']==t['task_id']) for t in ts}
    refs=[]
    for metal,stage in [('Ca','origin'),('La','origin'),('Ca','proposal'),('La','proposal'),('Ca','origin')]:
        f=fixtures[metal];point=f['points'][stage]
        refs.append({'id':f'{len(refs):02d}_{metal}_{stage}','task_id':f['task_id'],'metal':metal,'stage':stage,
                    'xyz':point['coordinate'],'forces':point['forces'],'energy_eV':point['MACE_eV'],'native_receipt':point['MACE']})
    for f in refs:
        for key in ['xyz','forces','native_receipt']:verify(f[key])
    verify(source['model']['checkpoint'])
    names=['slsqp_precision.py','adaptive_completion.py','adaptive_angular_proposals.py','accommodation_proposals.py',
           'accommodation_nonlinear.py','mace_site_kinematics.py','adaptive_force_diagnostic.py','mace_omol.py',
           'mace_hybrid.py','affordable_common.py','accommodation_torsion_profiles.py','union_adaptive.py']
    source_modules={n:source['implementation'][n] for n in names}
    for p in source_modules.values():verify(p)
    worker_root=Path(source['gpu_python']).parent.parent
    m={'protocol_id':'isolated_native_OMOL_CPU_force_two_proposal_parity_20261012_v1',
       'implementation':pin(__file__),'plan':pin(HERE/'SCOUT_PLAN.md'),'source_manifest':pin(ARCHIVE/'manifest.json'),
       'feasibility_inventory':pin(HERE/'RESULT.json'),'model':source['model'],'tasks':ts,'fixtures':refs,
       'archived_fixtures':fixtures,'settings':source['settings'],'optimizer_software':source['optimizer_software'],
       'source_modules':source_modules,'source_snapshot':str(ARCHIVE/'implementation'),
       'optimizer_python':source['cpu_python'],'model_python':source['gpu_python'],
       'calculator_source':pin(worker_root/'lib/python3.11/site-packages/mace/calculators/mace.py'),
       'threads':8,'parent_BLAS_threads':1,'sequential_parent_child_numerics':True,
       'gates':{'energy_kcal_mol':.01,'force_max_eV_A':.001,'proposal_energy_kcal_mol':.001},
       'force_gate_precedent':source_modules['mace_omol.py'],'proposal_gate_precedent':source_modules['slsqp_precision.py'],
       'finite_scope':{'force_replays':5,'optimizer_starts':2,'additional_call_or_wall_cap':None,'new_solvent_or_DFT_calls':0},
       'authorization':'Jacob2026-10-12<=16 aggregate login threads; root assigns8 total active threads to this sequential scout.'}
    write(out/'manifest.json',m)
    print(json.dumps({'manifest':pin(out/'manifest.json'),'prepared':True,'new_molecular_calls':0}))

def worker(manifest):
    m=read(manifest);verify(m['implementation']);verify(m['calculator_source']);verify(m['model']['checkpoint']);modules(m)
    import resource, importlib.metadata
    import numpy as np
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    out=Path(manifest).parent;began=time.monotonic();cpu=time.process_time();attempts=0;forwards=0;failures=0
    torch.set_num_threads(m['threads']);torch.set_num_interop_threads(1);torch.set_default_dtype(torch.float64)
    calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cpu',default_dtype='float64');model=calc.models[0]
    if (type(model).__name__!='ScaleShiftMACE' or list(model.heads)!=['omol'] or
        dict(model.embedding_specs)!=m['model']['embedding_specs'] or calc.energy_units_to_eV!=1.):raise ValueError('native model differs')
    for p in model.parameters():p.requires_grad_(False)
    parameter_versions={n:p._version for n,p in model.named_parameters()}
    def before_forward(*args):
        nonlocal forwards
        forwards+=1
    hook=model.register_forward_pre_hook(before_forward)
    load=time.monotonic()-began
    software={n:importlib.metadata.version(n) for n in ['torch','mace-torch','numpy','ase','e3nn']}
    print(PREFIX+json.dumps({'status':'ready','load_seconds':load,'threads':torch.get_num_threads(),'versions':software}),flush=True)
    try:
        for line in sys.stdin:
            p=json.loads(line)['request'];rp=verify(p);req=read(rp)
            result={'request':p,'manifest':pin(manifest),'task_id':req['task_id'],'status':'failed',
                    'energy_eV':None,'forces':None,'model_call_started':False,'device':'cpu','model':m['model']}
            tick=time.monotonic();before=forwards;cpu0=time.process_time();attempts+=1
            try:
                if req['manifest']!=pin(manifest):raise ValueError('request manifest differs')
                t=next(t for t in m['tasks'] if t['task_id']==req['task_id'])
                if (req['charge'],req['multiplicity'])!=(t['charge'],t['multiplicity']):raise ValueError('state mismatch')
                rows=xyz(verify(req['xyz']));atoms=Atoms([a[0] for a in rows],positions=[a[1:] for a in rows],pbc=False)
                atoms.info.update(charge=t['charge'],spin=t['multiplicity']);atoms.calc=calc
                state=input_batch(calc,atoms,t['charge'],t['multiplicity'])
                if state['graph_count']!=1 or rows[0][0]!=t['metal']:raise ValueError('graph or metal state mismatch')
                # Clear ASE's result cache, retaining the loaded model. Every requested
                # replay is a real force call; energy/forces within a request share one.
                calc.reset()
                write(rp.parent/'CPU_CALL_STARTED.json',{'request':p,'pid':os.getpid(),'time_unix':time.time()})
                result['model_call_started']=True
                energy=float(atoms.get_potential_energy());force=np.asarray(atoms.get_forces(),dtype=float)
                if not np.isfinite(energy) or force.shape!=(len(rows),3) or not np.isfinite(force).all():raise ValueError('nonfinite energy/force')
                if parameter_versions!={n:q._version for n,q in model.named_parameters()}:raise ValueError('model weights mutated')
                fp=rp.parent/'forces_eV_A.npy';np.save(fp,force,allow_pickle=False)
                result.update(status='complete',energy_eV=energy,forces=pin(fp),state=state,
                    force_definition='negative Cartesian total vacuum OMOL energy gradient, eV/angstrom',parameter_versions_unchanged=True)
            except Exception as exc:
                failures+=1;result.update(reason=str(exc),traceback=traceback.format_exc())
            result.update(wall_seconds=time.monotonic()-tick,process_cpu_seconds=time.process_time()-cpu0,
               actual_model_forwards=forwards-before,threads=torch.get_num_threads(),versions=software,
               peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scheduler_allocation=None)
            dest=rp.parent/'mace_result.json';write(dest,result)
            print(PREFIX+json.dumps({'result':pin(dest)}),flush=True)
    finally:
        hook.remove()
        write(out/'CPU_WORKER_SUMMARY.json',{'manifest':pin(manifest),'attempted_requests':attempts,'actual_model_forwards':forwards,
            'failed_requests':failures,'model_load_seconds':load,'wall_seconds':time.monotonic()-began,
            'process_cpu_seconds':time.process_time()-cpu,'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'threads':torch.get_num_threads(),'versions':software,'device':'cpu','allocated_gpu_seconds':0})

class WarmCPU:
    def __init__(self,manifest,m):
        self.manifest=Path(manifest);self.log=(self.manifest.parent/'CPU_WORKER.log').open('x');self.calls=[]
        env=os.environ.copy()
        for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OMP_THREAD_LIMIT']:env[key]=str(m['threads'])
        for key in ['OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
        env.pop('SLURM_JOB_ID',None)
        command=[m['model_python'],'-u',str(verify(m['implementation'])),'worker','--manifest',str(self.manifest)]
        self.process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                      text=True,bufsize=1,start_new_session=True,env=env)
        write(self.manifest.parent/'CPU_WORKER_START.json',{'pid':self.process.pid,'pgid':os.getpgid(self.process.pid),
              'command':command,'parent_pid':os.getpid(),'thread_environment':{k:env[k] for k in env if 'THREAD' in k}})
        try:self.ready=self.receive()
        except BaseException:self.abort();raise
        if self.ready.get('status')!='ready':self.abort();raise RuntimeError('CPU worker did not initialize')
    def receive(self):
        while True:
            line=self.process.stdout.readline()
            if not line:raise RuntimeError('CPU worker ended: '+str(self.process.poll()))
            self.log.write(line);self.log.flush()
            if line.startswith(PREFIX):return json.loads(line[len(PREFIX):])
    def evaluate(self,p):
        self.process.stdin.write(json.dumps({'request':p})+'\n');self.process.stdin.flush()
        response=self.receive();r=read(verify(response['result']))
        if r['request']!=p:raise ValueError('worker response belongs to another request')
        self.calls.append(response['result'])
        print(json.dumps({'call':len(self.calls),'task':r['task_id'],'status':r['status'],'wall_seconds':r['wall_seconds'],
                          'forwards':r['actual_model_forwards']}),flush=True)
        return r,response['result']
    def abort(self):
        if self.process.poll() is None:
            os.killpg(self.process.pid,signal.SIGTERM)
            try:self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(self.process.pid,signal.SIGKILL);self.process.wait()
        if not self.log.closed:self.log.close()
    def close(self):
        if self.process.poll() is None:
            self.process.stdin.close()
            for line in self.process.stdout:self.log.write(line)
            code=self.process.wait()
        else:code=self.process.returncode
        self.log.close()
        if code:raise RuntimeError('worker exit '+str(code))

def execute(manifest):
    import resource
    m=read(manifest);out=Path(manifest).parent;verify(m['implementation']);verify(m['plan']);modules(m)
    if any(os.environ.get(k)!='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']):
        raise ValueError('optimizer parent requires one numerical thread before imports')
    import numpy as np
    import scipy
    import scipy.optimize._slsqplib as kernel
    import scipy.optimize._slsqp_py as wrapper
    import slsqp_precision
    from adaptive_force_diagnostic import project,preview
    from mace_site_kinematics import Kinematics
    if (pin(kernel.__file__)!=m['optimizer_software']['kernel'] or pin(wrapper.__file__)!=m['optimizer_software']['wrapper']
       or scipy.__version__!=m['optimizer_software']['version']):raise ValueError('original optimizer binary required')
    started=time.monotonic();cpu0=time.process_time();gpu=None;replays=[];proposals=[];failure=None;status='failed'
    write(out/'RUN_STARTED.json',{'pid':os.getpid(),'manifest':pin(manifest),'time_unix':time.time(),
          'optimizer_python_invocation':sys.executable,'numpy_version':np.__version__,'scipy_version':scipy.__version__,
          'threads_total_reserved':m['threads'],'active_parent_child_overlap':False})
    def stopped(signum,frame):raise KeyboardInterrupt('owned scout received signal '+str(signum))
    old={s:signal.getsignal(s) for s in [signal.SIGTERM,signal.SIGHUP,signal.SIGINT]}
    for s in old:signal.signal(s,stopped)
    try:
        gpu=WarmCPU(manifest,m);origin_forces={};first=None
        for f in m['fixtures']:
            t=next(t for t in m['tasks'] if t['task_id']==f['task_id']);td=out/'replays'/f['id'];td.mkdir(parents=True)
            rp=td/'mace_request.json';write(rp,{'manifest':pin(manifest),'task_id':t['task_id'],'xyz':f['xyz'],
                                             'charge':t['charge'],'multiplicity':t['multiplicity']})
            r,rpin=gpu.evaluate(pin(rp))
            if r['status']!='complete':raise RuntimeError('CPU fixture failed: '+str(r.get('reason')))
            force=np.load(verify(r['forces']),allow_pickle=False);arch=np.load(verify(f['forces']),allow_pickle=False)
            de=(r['energy_eV']-f['energy_eV'])*EV_TO_KCAL;df=float(np.max(np.abs(force-arch)))
            check={'id':f['id'],'receipt':rpin,'energy_error_kcal_mol':de,'force_max_error_eV_A':df,
                   'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']}
            if f['stage']=='origin' and f['metal'] not in origin_forces:origin_forces[f['metal']]=force
            if first is None:first=(r['energy_eV'],force.copy())
            if f['id'].startswith('04'):
                check.update(repeat_energy_error_kcal_mol=(r['energy_eV']-first[0])*EV_TO_KCAL,
                             repeat_force_max_error_eV_A=float(np.max(np.abs(force-first[1]))))
                check['pass'] &= abs(check['repeat_energy_error_kcal_mol'])<=m['gates']['energy_kcal_mol'] and check['repeat_force_max_error_eV_A']<=m['gates']['force_max_eV_A']
            replays.append(check);write(td/'CHECK.json',check)
            if not check['pass']:raise RuntimeError('CPU energy/force parity gate failed')
        projections=[]
        for t in m['tasks']:
            kin=Kinematics(read(verify(t['mapping']))['context'])
            projections.append(project(kin.data,np.zeros(len(kin.modes)),origin_forces[t['metal']]))
        selected=preview(kin.modes,projections[0][1],projections[0][3],projections[1][3]);expected=m['tasks'][0]['selector']
        def categorical(s):
            x=copy.deepcopy(s)
            for v in x['selected']:
                for k in ['Ca_kcal_mol_A','La_kcal_mol_A','differential_kcal_mol_A','independent_geometric_fraction']:v.pop(k)
            return x
        selection={'fresh':selected,'archived':expected,'ordered_modes_and_categories_match':categorical(selected)==categorical(expected)}
        write(out/'PHASE1.json',{'replays':replays,'selection':selection,'pass':all(x['pass'] for x in replays) and selection['ordered_modes_and_categories_match']})
        if not selection['ordered_modes_and_categories_match']:raise RuntimeError('paired mode selection differs')
        # Private module copies set the historical ftol1e-8 without editing files.
        completion=slsqp_precision.engine()[0].completion
        if completion.SETTINGS!=m['settings']:raise ValueError('exact precision settings differ')
        for t in m['tasks']:
            rp=completion.optimize(manifest,t,gpu);r=read(verify(rp));arch=m['archived_fixtures'][t['metal']]
            row={'task_id':t['task_id'],'receipt':rp,'status':r['status'],'archived_receipt':arch['proposal_receipt']}
            if r['status']=='proposal_available':
                a=np.array([x[1:] for x in xyz(verify(r['proposal']['coordinate']))]);b=np.array([x[1:] for x in xyz(verify(arch['points']['proposal']['coordinate']))])
                de=(r['proposal']['MACE_eV']-arch['points']['proposal']['MACE_eV'])*EV_TO_KCAL
                row.update(energy_error_kcal_mol=de,maximum_coordinate_difference_A=float(np.linalg.norm(a-b,axis=1).max()),
                   full_q_difference=(np.array(r['proposal']['full_q'])-np.array(arch['points']['proposal']['full_q'])).tolist(),
                   boundary_flag=r['boundary_flag'],boundary_matches=r['boundary_flag']==arch['historical_cost']['boundary_flag'],
                   final_geometry=r['final_geometry'],optimizer=r['optimizer'],residual=r['residual'])
                row['pass']=abs(de)<=m['gates']['proposal_energy_kcal_mol'] and row['boundary_matches'] and r['final_geometry']['physical_feasible']
            else:row.update(reason=r.get('reason'),pass_=False);row['pass']=False
            proposals.append(row);write(out/(t['metal']+'_PROPOSAL_CHECK.json'),row)
            if not row['pass']:raise RuntimeError('CPU proposal parity/admission gate failed')
        status='complete'
    except BaseException as exc:
        failure=f'{type(exc).__name__}: {exc}';traceback.print_exc()
    finally:
        if gpu:
            try:
                if isinstance(failure,str) and failure.startswith('KeyboardInterrupt'):gpu.abort()
                else:gpu.close()
            except BaseException as exc:gpu.abort();failure=(failure+'; ' if failure else '')+str(exc);status='failed'
        for s,h in old.items():signal.signal(s,h)
        summary=read(out/'CPU_WORKER_SUMMARY.json') if (out/'CPU_WORKER_SUMMARY.json').exists() else None
        write(out/'FINAL.json',{'status':status,'failure':failure,'manifest':pin(manifest),'replays':replays,'proposals':proposals,
           'actual_requests':len(gpu.calls) if gpu else 0,'worker_summary':summary,
           'wall_seconds':time.monotonic()-started,'parent_process_cpu_seconds':time.process_time()-cpu0,
           'parent_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
           'model_receipts':gpu.calls if gpu else [],'all_scientific_gates_pass':status=='complete',
           'new_solvent_calls':0,'new_DFT_calls':0,'new_optimizer_starts':len(proposals),
           'source_scientific_method_unchanged':True,'full107_campaign_launched':False})
    print(json.dumps({'status':status,'failure':failure,'actual_requests':len(gpu.calls) if gpu else 0}),flush=True)
    if status!='complete':raise SystemExit(1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('operation',choices=['prepare','execute','worker'])
    parser.add_argument('--output');parser.add_argument('--manifest');args=parser.parse_args()
    if args.operation=='prepare':prepare(args.output)
    elif args.operation=='worker':worker(args.manifest)
    else:execute(args.manifest)
