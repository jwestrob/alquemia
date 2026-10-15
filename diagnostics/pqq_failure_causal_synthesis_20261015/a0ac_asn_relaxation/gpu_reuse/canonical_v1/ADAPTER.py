#!/usr/bin/env python3
"""Request/admission bridge; import the qualified canonical scientific runner unchanged."""
import argparse,copy,fcntl,hashlib,importlib.util,json,os,shutil,signal,subprocess,sys,time,traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/gpu_reuse/canonical_v1'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
KIND='qualified_GPU_request_bridge_v1'
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
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def bind(output):
    """Metadata and immutable hashes only; no scientific imports."""
    base=pin(WORK/'asn_canonical_search_remaining_v1/CONFIG.json');c=read(base);m=read(c['worker_manifest'])
    qp=pin(WORK/'asn_gpu_reuse_root_v2/execution/RESULT.json');q=read(qp);qc=read(q['config']);qe=pin(WORK/'asn_gpu_reuse_root_v2/execution/EXECUTION.json');e=read(qe)
    require(q['all_four_pass'] and q['status']==e['status']=='complete' and all(r['pass_all'] for r in q['rows']) and len(q['rows'])==4,'Actual four-point qualification required')
    require(c['stage']=='remaining' and len(c['rows'])==23 and len(m['tasks'])==46 and read(c['scout_collection'])['stage_complete'],'Exact completed scout/remaining scope required')
    require(not Path(c['output']).exists(),'CPU fallback has started; do not duplicate canonical work')
    gm=read(qc['gpu_manifest']);require(gm['model']==m['model'],'Qualified GPU model differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for n in ('ADAPTER.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/n,out/n)
    gm.update(tasks=copy.deepcopy(m['tasks']),source_worker_manifest=c['worker_manifest'],protocol_id=KIND)
    gm.pop('new_requests',None);gm.pop('no_optimizer',None);write(out/'gpu_manifest.json',gm)
    c.update(output=str(out/'execution'),wrapper=pin(out/'RUN.sbatch'),execution_adapter=pin(out/'ADAPTER.py'),GPU_plan=pin(out/'PLAN.md'),
        base_CPU_config=base,GPU_manifest=pin(out/'gpu_manifest.json'),GPU_qualification=qp,GPU_qualification_execution=qe,
        GPU_qualification_config=q['config'],GPU_launcher=qc['implementation'],GPU_server=qc['gpu_implementation'],
        GPU_qualified_runtime=pin(WORK/'asn_gpu_reuse_root_v2/execution/PREFLIGHT.json'),GPU_transport=KIND)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def gpu_admission(result,request,task,model,versions,np,canonical,b):
    """Public receipt check for downstream reuse; native facts stay in the raw receipt."""
    require(result.get('transport_kind')==KIND,'Unknown GPU adapter receipt')
    c=read(result['adapter_config']);require(c['GPU_transport']==KIND,'Unknown transport configuration')
    raw=read(result['native_GPU_receipt']);bridge=read(result['bridge_request']);rt=read(result['GPU_runtime'])
    require(rt['status']=='passed' and rt['config']==result['adapter_config'] and rt['qualification']==c['GPU_qualification'],'GPU runtime qualification differs')
    q=read(c['GPU_qualification']);require(q['all_four_pass'] and read(c['GPU_qualification_execution'])['status']=='complete','GPU parity not qualified')
    require(rt['runtime']['versions']==versions and raw['device']==rt['runtime']['device']=='NVIDIA RTX A5000','Actual GPU runtime/device differs')
    original=read(result['request']);require(original==request and raw['request']==result['bridge_request'],'Request provenance differs')
    wm=read(c['worker_manifest']);t,atoms=b.authorize_request(check(c['worker_manifest']),wm,request)
    root,=[x for x in wm['tasks'] if (x['case_id'],x['metal'])==(t['case_id'],t['metal'])]
    expected={'manifest':c['GPU_manifest'],'task_id':root['task_id'],**{k:request[k] for k in ('xyz','charge','multiplicity')}}
    require(bridge==expected and result['task_id']==request['task_id'] and result['manifest']==request['manifest'],'Flat-manifest authorization changed coordinates/state')
    require(read(c['GPU_manifest'])['tasks']==wm['tasks'] and len(atoms)==task['atom_count'],'GPU task/count mismatch')
    for k,v in raw.items():
        if k!='request':require(result[k]==v,'Adapter altered native result: '+k)
    require(raw['status']=='complete' and raw['model_call_started'],'Native GPU model call unavailable')
    canonical.physical_result(raw,bridge,task,model,raw['model'])
    f=np.load(check(raw['forces']),allow_pickle=False);require(f.shape==(task['atom_count'],3) and np.isfinite(f).all(),'GPU forces invalid')

class Bridge:
    def __init__(self,c,cp,m,b,native,launcher,runtime):
        self.c=c;self.cp=cp;self.m=m;self.b=b;self.runtime=runtime;self.calls=[];self.launcher=launcher
        self.gpu=launcher.direct_transport(native,check(c['GPU_manifest']))
    def evaluate(self,p):
        req=read(p);t,atoms=self.b.authorize_request(check(self.c['worker_manifest']),self.m,req)
        base,=[v for v in self.m['tasks'] if (v['case_id'],v['metal'])==(t['case_id'],t['metal'])]
        require(len(atoms)==t['atom_count'],'Request count differs')
        d=check(p).parent/'GPU_NATIVE';d.mkdir(exist_ok=False)
        request={'manifest':self.c['GPU_manifest'],'task_id':base['task_id'],**{k:req[k] for k in ('xyz','charge','multiplicity')}}
        rp=d/'request.json';write(rp,request)
        write(d/'DISPATCH_STARTED.json',{'original_request':p,'bridge_request':pin(rp),'adapter_config':self.cp})
        raw,rawpin=self.gpu.evaluate(pin(rp))
        value={**raw,'request':p,'manifest':req['manifest'],'task_id':req['task_id'],'transport_kind':KIND,
            'adapter_config':self.cp,'native_GPU_receipt':rawpin,'bridge_request':pin(rp),'GPU_runtime':self.runtime}
        dest=check(p).parent/'mace_result.json';write(dest,value);self.calls.append(pin(dest))
        print(json.dumps({'task':req['task_id'],'status':raw['status'],'GPU_wall_seconds':raw['wall_seconds']}),flush=True)
        return value,pin(dest)
    def close(self):self.gpu.close()
    def abort(self):self.launcher.stop_owned(self.gpu)

def runtime(c,cp,launcher,out):
    gm=read(c['GPU_manifest']);q=read(c['GPU_qualification']);qc=read(c['GPU_qualification_config']);qualified=read(c['GPU_qualified_runtime'])
    require(q['all_four_pass'] and read(c['GPU_qualification_execution'])['status']=='complete','GPU qualification failed')
    for k in ('GPU_plan','wrapper','GPU_launcher','GPU_server'):check(c[k])
    prior=read(qc['gpu_manifest'])
    for k in ('implementation','gpu_python','gpu_executable','model'):require(gm[k]==prior[k],'Qualified GPU field changed: '+k)
    require(pin(gm['gpu_python'])==gm['gpu_executable'] and gm['implementation']['accommodation_nonlinear.py']==c['GPU_server'],'GPU executable differs')
    for p in gm['implementation'].values():check(p)
    check(gm['model']['checkpoint'])
    cmd=[gm['gpu_python'],'-c',launcher.PROBE];v=subprocess.run(cmd,text=True,capture_output=True)
    write(out/'GPU_RUNTIME_PROBE.json',{'argv':cmd,'returncode':v.returncode,'stdout':v.stdout,'stderr':v.stderr});require(v.returncode==0,'GPU probe failed')
    payload,=[x.removeprefix('GPU_RUNTIME_JSON ') for x in v.stdout.splitlines() if x.startswith('GPU_RUNTIME_JSON ')];rt=json.loads(payload)
    require(rt['CUDA_available'] and rt['visible_GPU_count']==1 and rt['device']==qualified['runtime']['device']=='NVIDIA RTX A5000','Unqualified GPU device')
    for k in ('versions','prefix','executable','calculator_source','torch_version','CUDA_version'):require(rt[k]==qualified['runtime'][k],'GPU runtime changed: '+k)
    require(pin(rt['calculator_source'])==qc['calculator_source'] and rt['cgroup']==Path('/proc/self/cgroup').read_text() and rt['CUDA_VISIBLE_DEVICES']==os.environ['CUDA_VISIBLE_DEVICES'],'GPU process scope/calculator differs')
    rp=out/'GPU_RUNTIME.json';write(rp,{'status':'passed','config':cp,'qualification':c['GPU_qualification'],'runtime':rt,'new_model_calls':0});return pin(rp)

def run(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_JOB_NUM_NODES')=='1','Actual one-node Slurm required before scientific imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(pin(__file__)==c['execution_adapter'] and c['GPU_transport']==KIND,'Wrong adapter')
    base=read(c['base_CPU_config'])
    for k,v in base.items():
        if k not in ('output','wrapper'):require(c[k]==v,'Canonical scientific configuration changed: '+k)
    require(not Path(base['output']).exists(),'CPU fallback started; do not duplicate')
    require(os.environ.get('SLURM_GPUS_ON_NODE')=='1' and bool(os.environ.get('SLURM_JOB_GPUS')) and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles') and len(os.environ['CUDA_VISIBLE_DEVICES'].split(','))==1,'Actual one-GPU grant required')
    require(int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])>=1,'Actual CPU grant required')
    out=Path(c['output']);start=time.monotonic();cpu=time.process_time();transport=None;completion=None;original=None;handlers={};pf=None
    allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_JOB_PARTITION','SLURM_JOB_GPUS','SLURM_GPUS_ON_NODE','SLURM_CPUS_ON_NODE','CUDA_VISIBLE_DEVICES')}
    receipt={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'allocation':allocation,'stage':c['stage'],'new_SCC_calls':0,'new_DFT_calls':0,'actual_model_forwards':None}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        canonical=load('_GPU_canonical_unchanged',c['implementation']);m=read(c['worker_manifest']);b=load('_GPU_canonical_batch',m['implementation'])
        scorer=load('_GPU_canonical_scorer',c['old_scorer']);scorer.physical_result=canonical.physical_result;old_admission=scorer.admitted_result
        try:
            b.parent_env();b.modules(m)
            import numpy as np
            pf,completion=canonical.preflight(c,b,scorer,np);write(out/'PREFLIGHT.json',pf)
            launcher=load('_GPU_qualified_launcher',c['GPU_launcher']);rt=runtime(c,cp,launcher,out)
            def admit(v,r,t,model,versions):
                if v.get('transport_kind')==KIND:return gpu_admission(v,r,t,model,versions,np,canonical,b)
                return old_admission(v,r,t,model,versions)
            scorer.admitted_result=admit
            helper=load('_GPU_canonical_contacts',c['generic_contact_helper']);original,guarded=canonical.guarded_class(completion,helper,np);completion.CompletionProposal=guarded
            handlers=b.handlers();native=load('_GPU_native_server',c['GPU_server']);transport=Bridge(c,cp,m,b,native,launcher,rt)
            receipt.update(GPU_command=transport.gpu.command,worker_cgroup=transport.gpu.worker_cgroup)
            partial=[]
            for row in c['rows']:
                r=canonical.source_run(c,row,m,b,scorer,transport,completion,np,pf,out)
                partial.append({'source_case_id':row['source_case_id'],'result':pin(out/'sources'/row['source_case_id']/'SOURCE_RESULT.json'),'status':r['status']})
                write(out/'progress'/f'{len(partial):02d}.json',{'completed_source_attempts':partial.copy()})
            requests=len(transport.calls);transport.close();transport=None
            sp=Path(c['GPU_manifest']['path']).parent/('gpu_summary_'+os.environ['SLURM_JOB_ID']+'.json');s=read(sp);receipt['native_GPU_summary']=pin(sp)
            require(s['manifest']==c['GPU_manifest'] and s['new_MACE_calls']==requests and s['failed_requests']==0,'Native request completion differs')
            receipt.update(status='complete',native_model_requests=s['new_MACE_calls'],native_failed_requests=s['failed_requests'])
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
        finally:
            try:
                if transport:transport.abort()
            except BaseException as e:receipt.update(status='failed',cleanup_error=str(e))
            try:
                if completion is not None and original is not None:completion.CompletionProposal=original
                scorer.admitted_result=old_admission
                for sig,h in handlers.items():signal.signal(sig,h)
                collection=canonical.collect(c,out,scorer,pf);collection['GPU_transport']=KIND;write(out/'FINAL_COLLECTION.json',collection)
                if not collection['stage_complete']:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,parent_CPU_seconds=time.process_time()-cpu,allocated_GPU_seconds=elapsed,
                    allocated_CPU_seconds=int(os.environ['SLURM_CPUS_ON_NODE'].split('(')[0])*elapsed)
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('bind');a.add_argument('--output',required=True)
    r=s.add_parser('run');r.add_argument('--config',required=True);r.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='bind':bind(a.output)
    else:sys.exit(run(a.config,a.config_sha))
