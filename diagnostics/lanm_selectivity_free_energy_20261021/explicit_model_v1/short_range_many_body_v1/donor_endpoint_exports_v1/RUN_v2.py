"""Eight baseline exports. Root release and actual standard_expand Slurm required."""
import argparse, concurrent.futures, hashlib, importlib.util, json, multiprocessing, os, platform, re, subprocess, sys, time
from pathlib import Path

_WORKER_THREADS_INITIALIZED = False

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
    assert pin(p['path'])==p;return json.loads(Path(p['path']).read_text())
def save(p,x):
    p=Path(p);temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');temp.replace(p)

def one(manifest_pin,id,folder):
    start,cpu=time.monotonic(),time.process_time();m=read(manifest_pin);cfg=m['methods'].copy()
    for p in cfg.values():assert pin(p['path'])==p
    parent=Path(cfg['implementation']['path']);sys.path.insert(0,str(parent.parent))
    spec=importlib.util.spec_from_file_location('donor_baseline_parent',parent);h=importlib.util.module_from_spec(spec);sys.modules[spec.name]=h;spec.loader.exec_module(h)
    assert pin(h.PROTOCOL.__file__)==cfg['protocol']
    assert pin(sys.modules['PREPARE'].__file__)==cfg['metadata_helper']
    global _WORKER_THREADS_INITIALIZED
    if not _WORKER_THREADS_INITIALIZED:
        h.torch.set_num_threads(1);h.torch.set_num_interop_threads(1)
        _WORKER_THREADS_INITIALIZED = True
    cfg.update(units=m['units'],timeout_seconds=m['timeout_seconds'])
    kernel=h.module(cfg['kernel']);combiner=h.module(cfg['combiner']);slot=m['states'][id]
    r=h.run_state(cfg,m,slot,kernel,combiner,Path(folder))
    if r['status']=='native_force_state_complete':
        # Parse only the already-produced export; no second solve or native call.
        b=h.PROTOCOL.baseline(r['baseline_state']['path'],slot,read(slot['fixture']))
        fixture={'schema':'short_range_many_body_donor_endpoint_baseline_v1','id':id,
                 'baseline':m['shared_baseline'],'geometry':slot['geometry'],'inputs':slot['inputs'],
                 'compatible':True,'baseline_observables':{'energy_kcal_mol':r['energy_kcal_mol'],
                 'gradient_kcal_mol_A':r['gradient_kcal_mol_A'],'forces_kcal_mol_A':r['forces_kcal_mol_A']},
                 'native_operator':{'A0_Aminus3':b['A'],'Ed_e_per_A2':b['Ed'],'Ep_e_per_A2':b['Ep'],
                 'alpha_A3':[s['alpha'] for s in b['sites']],'site_atom_indices':[s['atom']-1 for s in b['sites']],
                 'K_kcal_A_per_mol_e2':b['K'],'permanent_mu_eA':b['permanent_mu_eA'],'permanent_mu_origin_A':[0.,0.,0.],
                 'nonpolar_energy_kcal':r['other_energy_kcal_mol'],'u_zero_dense_eA':r['u'],'v_zero_dense_eA':r['v'],
                 'raw_operator_retained':True,'coordinate_derivatives':None},
                 'native_qualification':{'producer':pin(Path(folder)/'RESULT.json'),'baseline_state':r['baseline_state'],
                 'force_bridge_result':r['force_bridge_result'],'symmetry_projection':r['symmetry_projection'],
                 'minimum_eigenvalue':r['minimum_corrected_eigenvalue'],'solve_residual':r['solve_residual'],
                 'scope':'fresh zero-correction dense baseline state; not an archived iterative solve or full derivative provider'},
                 'historical_evidence':m['evidence'][id],'consumed_development':True,
                 'new_quantum_calls':0,'model_coefficients_fitted':False}
        save(Path(folder)/'FIXTURE.json',fixture)
    return {'result':pin(Path(folder)/'RESULT.json'),'fixture':pin(Path(folder)/'FIXTURE.json') if r['status']=='native_force_state_complete' else None,
            'status':r['status'],'error':r['error'],'native_program_attempts':r['native_program_attempts'],
            'native_program_invocations':r['native_program_invocations'],'native_cost':r['cost'],
            'worker_cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}}

def main():
    a=argparse.ArgumentParser();a.add_argument('--manifest',type=Path,required=True);a.add_argument('--release',type=Path,required=True);a.add_argument('--out',type=Path,required=True);args=a.parse_args()
    start,cpu=time.monotonic(),time.process_time();mp=pin(args.manifest);m=read(mp);release=json.loads(args.release.read_text())
    assert m['runtime_source']==pin(__file__) and release['implementation']==pin(__file__) and release['manifest']==mp
    assert release['status']=='approved_eight_saved_donor_baseline_exports' and release['output']==str(args.out.resolve())
    assert release['maximum_native_programs']==8 and release['workers']==4 and release['qos']=='standard_expand'
    assert len(m['ordered_ids'])==len(set(m['ordered_ids']))==len(m['states'])==8
    assert all(s['direction'] in ('radial','rotation') and s['coefficient_mode']=='zero' and s['parameters']=={'pairs':{},'triples':[]} for s in m['states'].values())
    job=os.environ['SLURM_JOB_ID'];assert job.isdigit() and int(os.environ['SLURM_CPUS_PER_TASK'])>=8
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):assert os.environ[key]=='1'
    command=['scontrol','show','job','-o',job];allocation=subprocess.run(command,text=True,capture_output=True,check=True)
    assert re.search(r'\bQOS=standard_expand\b',allocation.stdout)
    args.out.mkdir(parents=True,exist_ok=False)
    r={'status':'started','implementation':pin(__file__),'manifest':mp,'root_release':pin(args.release),'rows':{},
       'allocation':{'job_id':job,'allocated_CPUs':int(os.environ['SLURM_CPUS_PER_TASK']),'host':platform.node(),'pid':os.getpid(),
       'inspection_argv':command,'scontrol_stdout':allocation.stdout},'error':None,'native_program_attempts':0,'native_program_invocations':0}
    save(args.out/'RESULT.json',r)
    try:
        # Two finite waves: stop before the second wave if the first has a failure.
        with concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
            for offset in (0,4):
                jobs={pool.submit(one,mp,id,str(args.out/'states'/id)):id for id in m['ordered_ids'][offset:offset+4]}
                for future in concurrent.futures.as_completed(jobs):
                    id=jobs[future]
                    try:row=future.result()
                    except Exception as exc:
                        row={'status':'failed_adapter','error':repr(exc),'result':None,'fixture':None,'native_program_attempts':None,'native_program_invocations':None}
                        rp=args.out/'states'/id/'RESULT.json'
                        if rp.exists():
                            actual=json.loads(rp.read_text());row.update(result=pin(rp),native_program_attempts=actual['native_program_attempts'],native_program_invocations=actual['native_program_invocations'])
                    r['rows'][id]=row;save(args.out/'RESULT.json',r)
                if any(v['status']!='native_force_state_complete' for v in r['rows'].values()):r['error']='Endpoint failed; completed/current wave retained and unstarted wave untouched';break
        for key in ('native_program_attempts','native_program_invocations'):
            values=[v[key] for v in r['rows'].values()];r[key]=None if any(v is None for v in values) else sum(values)
        r['status']='eight_saved_donor_baseline_exports_complete' if r['error'] is None and len(r['rows'])==8 else 'failed'
    except Exception as exc:r['status']='failed';r['error']=repr(exc)
    r['cost']={'parent_wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,
               'worker_and_native_costs':'per-row; native child costs are separate from worker controller CPU; scheduler accounting retained by root'}
    save(args.out/'RESULT.json',r);print(json.dumps({'status':r['status'],'rows':len(r['rows']),'error':r['error']}));return 0 if r['status'].endswith('_complete') else 1

if __name__=='__main__':raise SystemExit(main())
