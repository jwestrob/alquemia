#!/usr/bin/env python3
"""Four real CPU/GPU comparisons through the unchanged historical WarmGPU."""
import argparse,copy,fcntl,hashlib,importlib.util,json,os,shutil,socket,subprocess,sys,time,traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/gpu_reuse'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
SOURCE='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3'
METALS=('Ca','La');GEOMETRIES=('amide_flip','relaxed_flip_La');EV=23.06054783061903
PROTOCOL='existing_WarmGPU_four_real_CPU_points_compatibility_v1'
def require(v,why):
    if not v:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def freeze(output):
    pp=pin(WORK/'asn_canonical_relaxation_v1/execution/PREPARATION.json');p=read(pp);pm=read(p['manifest'])
    require(p['status']=='complete' and p['prepared_sources']==25,'Actual complete canonical geometry required')
    row,=[r for r in p['rows'] if r['source_case_id']==SOURCE];require(read(row['regression'])['status']=='passed' and row['executed_result_alias']['status']=='complete','Actual canonical reuse regression missing')
    old=read(WORK/'asn_rotamer_scoring_v1/CONFIG.json');wm=read(old['worker_manifest'])
    precedent=pin(ROOT/'workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/manifest.json');cm=read(precedent)
    historical=pin(ROOT/'workspaces/motion_envelope_20260923/pilot34_searches_v1/manifest.json');hm=read(historical)
    require(wm['model']==cm['model']==hm['model'] and cm['gates']=={'energy_kcal_mol':.01,'force_max_eV_A':.001},'Qualified model/gates differ')
    gp=hm['implementation']['accommodation_nonlinear.py'];require(gp['sha256']=='49f493c1d343ef14b645b653072e544c44e0e4739a04690ac5346f6be12adbaa','Historical GPU implementation differs')
    tasks=[];baseline={}
    for g in GEOMETRIES:
        for z in METALS:
            cell=row['executed_result_alias']['matrix'][z][g];r=read(cell['receipt']);req=read(r['request']);tid=SOURCE+'__'+g+'__'+z+'__GPU_replay'
            require(cell['status']=='complete' and r['status']=='complete' and req['xyz']==cell['xyz'] and r['model']==wm['model'],'Actual CPU point unavailable')
            check(cell['xyz']);check(cell['forces']);tasks.append({'task_id':tid,'case_id':SOURCE,'source_case_id':SOURCE,'candidate':g,'metal':z,'xyz':cell['xyz'],
                'charge':req['charge'],'multiplicity':req['multiplicity'],'atom_count':204});baseline[tid]=cell
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for n in ('SCOUT.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/n,out/n)
    gpu={'protocol_id':PROTOCOL,'implementation':hm['implementation'],'gpu_python':hm['gpu_python'],'gpu_executable':hm['gpu_executable'],
         'model':wm['model'],'tasks':tasks,'source_canonical_preparation':pp,'new_requests':4,'no_optimizer':True}
    write(out/'gpu_manifest.json',gpu)
    c={'protocol':PROTOCOL,'implementation':pin(out/'SCOUT.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(out/'PLAN.md'),'preparation':pp,'preparation_manifest':p['manifest'],
       'source_row':row,'completed_CPU_collection':pm['completed_search_collection'],'completed_CPU_execution':pm['completed_search_execution'],
       'gpu_manifest':pin(out/'gpu_manifest.json'),'historical_GPU_manifest':historical,'gpu_implementation':gp,'parent_python':old['python'],
       'calculator_source':cm['calculator_source'],'cpu_admission_helper':old['implementation'],'CPU_qualification':old['qualification'],'CPU_worker_summary':old['qualified_worker_summary'],
       'qualification_precedent':precedent,'gates':cm['gates'],'model':wm['model'],'baseline':baseline,'tasks':tasks,'new_GPU_requests':4,'new_CPU_requests':0,
       'kcal_per_eV':EV,'output':str(out/'execution')}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

PROBE="""import importlib.metadata,json,sys,os,torch
import mace.calculators.mace as calculator
print('GPU_RUNTIME_JSON '+json.dumps({'executable':sys.executable,'prefix':sys.prefix,'versions':{n:importlib.metadata.version(n) for n in ('torch','mace-torch','numpy','ase','e3nn')},'torch_version':torch.__version__,'CUDA_version':torch.version.cuda,'CUDA_available':torch.cuda.is_available(),'visible_GPU_count':torch.cuda.device_count(),'device':torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,'CUDA_VISIBLE_DEVICES':os.environ.get('CUDA_VISIBLE_DEVICES'),'calculator_source':calculator.__file__,'new_model_calls':0}))
"""
def preflight(c,np):
    require(c['protocol']==PROTOCOL and c['new_GPU_requests']==4 and c['new_CPU_requests']==0 and c['kcal_per_eV']==EV,'Finite scope changed')
    for k in ('wrapper','plan','preparation','preparation_manifest','completed_CPU_collection','completed_CPU_execution','historical_GPU_manifest','calculator_source'):check(c[k])
    require(read(c['completed_CPU_execution'])['status']=='complete','CPU execution not complete')
    require(pin(sys.executable)==c['parent_python'],'Parent interpreter differs')
    p=read(c['preparation']);row,=[r for r in p['rows'] if r['source_case_id']==SOURCE];require(row==c['source_row'] and read(row['regression'])['status']=='passed','Canonical source/alias changed')
    old=read(c['qualification_precedent']);require(c['gates']==old['gates']=={'energy_kcal_mol':.01,'force_max_eV_A':.001},'Existing parity gates changed')
    gm=read(c['gpu_manifest']);hm=read(c['historical_GPU_manifest']);require(gm['tasks']==c['tasks'] and len(c['tasks'])==4 and gm['model']==hm['model']==c['model']==old['model'],'GPU model/task binding differs')
    require(gm['implementation']==hm['implementation'] and gm['gpu_python']==hm['gpu_python']==old['model_python'] and pin(gm['gpu_python'])==gm['gpu_executable']==old['model_executable'],'Historical isolated GPU invocation differs')
    for v in gm['implementation'].values():check(v)
    check(c['model']['checkpoint']);check(c['gpu_implementation']);require(gm['implementation']['accommodation_nonlinear.py']==c['gpu_implementation'],'Wrong GPU server')
    summary=read(c['CPU_worker_summary']);qualification=read(c['CPU_qualification']);require(qualification['pass'],'Existing CPU qualification failed')
    scorer=load('_GPU_CPU_admission',c['cpu_admission_helper']);state=read(row['original_row']['state_signature']);coordinates={};forces={}
    for t in c['tasks']:
        z=t['metal'];g=t['candidate'];cell=c['baseline'][t['task_id']];saved=row['executed_result_alias']['matrix'][z][g]
        require(cell==saved and t['xyz']==cell['xyz'] and (t['charge'],t['multiplicity'],t['atom_count'])==(state['endpoint_charges'][z],1,204),'Exact CPU source/state join differs')
        result=read(cell['receipt']);request=read(result['request']);scorer.admitted_result(result,request,t,c['model'],summary['versions'])
        require(result['energy_eV']==cell['energy_eV'] and result['forces']==cell['forces'],'CPU energy/force pin differs')
        force=np.load(check(cell['forces']),allow_pickle=False);require(force.shape==(204,3) and np.isfinite(force).all(),'CPU forces invalid');forces[t['task_id']]=force
        atoms=[line.split() for line in check(t['xyz']).read_text().splitlines()[2:]];require(len(atoms)==204 and atoms[0][0]==z,'Wrong CPU point atom inventory')
        coordinates[(g,z)]=[(a[0],*map(float,a[1:])) for a in atoms]
    require([(t['candidate'],t['metal']) for t in c['tasks']]==[(g,z) for g in GEOMETRIES for z in METALS],'Four-point order/identity changed')
    for g in GEOMETRIES:require(coordinates[(g,'Ca')][1:]==coordinates[(g,'La')][1:] and coordinates[(g,'Ca')][0][1:]==coordinates[(g,'La')][0][1:],'Matched physical geometry differs')
    command=['srun','--exact','--overlap','--ntasks=1','--cpus-per-task=1','--gres=gpu:1',gm['gpu_python'],'-c',PROBE]
    probe=subprocess.run(command,text=True,capture_output=True);record={'argv':command,'returncode':probe.returncode,'stdout':probe.stdout,'stderr':probe.stderr}
    write(Path(c['output'])/'RUNTIME_PROBE.json',record);require(probe.returncode==0,'GPU runtime probe failed')
    payload,=[line.removeprefix('GPU_RUNTIME_JSON ') for line in probe.stdout.splitlines() if line.startswith('GPU_RUNTIME_JSON ')];runtime=json.loads(payload)
    require(runtime['CUDA_available'] and runtime['visible_GPU_count']==1 and runtime['versions']==summary['versions'],'Actual CUDA/runtime compatibility failed')
    require(Path(runtime['prefix']).resolve()==Path(gm['gpu_python']).parent.parent.resolve() and pin(runtime['calculator_source'])==c['calculator_source'],'Actual GPU environment/calculator differs')
    return {'status':'passed','runtime':runtime,'actual_GPU_worker_configured_threads':1,'threads_evidence':'Pinned historical gpu_server torch.set_num_threads(1), not a newly added native receipt field',
        'baseline_cells_checked':4,'gates':c['gates'],'model':c['model'],'new_model_calls':0},forces,scorer

def run(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_JOB_NUM_NODES','0'))==1,'Actual one-node Slurm allocation required before numerical imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(pin(__file__)==c['implementation'],'Wrong frozen scout')
    out=Path(c['output']);start=time.monotonic();cpu=time.process_time();gpu=None;results={};pf=None
    allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_JOB_PARTITION','SLURM_JOB_GPUS','SLURM_GPUS_ON_NODE','SLURM_CPUS_ON_NODE','SLURM_JOB_CPUS_PER_NODE','SLURM_MEM_PER_NODE','SLURM_MEM_PER_CPU','SLURM_CPUS_PER_GPU')}
    receipt={'status':'failed','config':cp,'allocation':allocation,'host':socket.gethostname(),'new_CPU_calls':0,'new_optimizer_calls':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Parent thread limits differ')
            gm=read(c['gpu_manifest']);sys.path.insert(0,str(Path(c['gpu_implementation']['path']).parent))
            import numpy as np
            pf,forces,scorer=preflight(c,np);write(out/'PREFLIGHT.json',pf)
            native=load('_historical_WarmGPU',c['gpu_implementation']);gpu=native.WarmGPU(check(c['gpu_manifest']));receipt['GPU_command']=gpu.command
            for t in c['tasks']:
                d=out/'cells'/t['task_id'];d.mkdir(parents=True,exist_ok=False);req={'manifest':c['gpu_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
                rp=d/'request.json';write(rp,req);write(d/'DISPATCH_STARTED.json',{'request':pin(rp),'time_unix':time.time()})
                v,vp=gpu.evaluate(pin(rp));row={'task_id':t['task_id'],'candidate':t['candidate'],'metal':t['metal'],'status':'unavailable','GPU_receipt':vp,'CPU':c['baseline'][t['task_id']],'xyz':t['xyz']}
                try:
                    require(v['request']==pin(rp) and v['model_call_started'],'GPU request not actually evaluated');scorer.physical_result(v,req,t,c['model'],v['model'])
                    force=np.load(check(v['forces']),allow_pickle=False);require(force.shape==(204,3) and np.isfinite(force).all(),'GPU force array invalid')
                    de=(v['energy_eV']-row['CPU']['energy_eV'])*EV;df=float(np.max(np.abs(force-forces[t['task_id']])))
                    checks={'energy_pass':abs(de)<=c['gates']['energy_kcal_mol'],'force_pass':df<=c['gates']['force_max_eV_A']}
                    row.update(status='complete',GPU_energy_eV=v['energy_eV'],GPU_forces=v['forces'],GPU_minus_CPU_energy_kcal_mol=de,maximum_force_component_difference_eV_A=df,
                        checks=checks,pass_all=all(checks.values()),GPU_wall_seconds=v['wall_seconds'],CPU_wall_seconds=read(row['CPU']['receipt'])['wall_seconds'])
                except Exception as exc:row.update(reason=str(exc),pass_all=False)
                results[t['task_id']]=row
            gpu.close();gpu=None;summary=Path(c['gpu_manifest']['path']).parent/('gpu_summary_'+os.environ['SLURM_JOB_ID']+'.json');s=read(summary)
            receipt['native_GPU_summary']=pin(summary)
            require(s['manifest']==c['gpu_manifest'] and s['new_MACE_calls']==4 and s['failed_requests']==0,'Native four-request completion unavailable')
            receipt.update(status='complete' if len(results)==4 and all(r['pass_all'] for r in results.values()) else 'parity_failed',native_GPU_summary=pin(summary))
        except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc())
        finally:
            try:
                if gpu is not None:
                    if gpu.process.poll() is None:gpu.process.terminate()
                    try:gpu.process.wait(timeout=10)
                    except subprocess.TimeoutExpired:gpu.process.kill();gpu.process.wait()
                    if not gpu.log.closed:gpu.log.close()
            except BaseException as exc:receipt.update(status='failed',cleanup_error=str(exc))
            try:
                rows=[results.get(t['task_id'],{'task_id':t['task_id'],'candidate':t['candidate'],'metal':t['metal'],'status':'unavailable','pass_all':False,'reason':'Not completed','CPU':c['baseline'][t['task_id']],'xyz':t['xyz']}) for t in c['tasks']]
                write(out/'RESULT.json',{'protocol':PROTOCOL,'config':cp,'status':receipt['status'],'all_four_pass':receipt['status']=='complete','rows':rows,'gates':c['gates'],
                    'new_GPU_requests_maximum':4,'new_CPU_requests':0,'model':c['model'],'runtime_admitted':pf is not None,'native_GPU_summary':receipt.get('native_GPU_summary'),
                    'limitation':'Compatibility of four fixed real points on the actual GPU only; no optimizer/classification/transport promotion'})
            except BaseException as exc:receipt.update(status='failed',reporting_error=str(exc),reporting_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(wall_seconds=elapsed,parent_CPU_seconds=time.process_time()-cpu)
                if allocation['SLURM_CPUS_ON_NODE']:receipt['allocated_CPU_seconds']=int(allocation['SLURM_CPUS_ON_NODE'].split('(')[0])*elapsed
                receipt['allocated_GPU_seconds']=elapsed;write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);f=s.add_parser('freeze');f.add_argument('--output',required=True)
    r=s.add_parser('run');r.add_argument('--config',required=True);r.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='freeze':freeze(a.output)
    else:sys.exit(run(a.config,a.config_sha))
