#!/usr/bin/env python3
"""One Slurm native contraction at archived u/v; no response solve or retry."""
import argparse, hashlib, importlib.util, json, math, os, platform, resource
import selectors, shutil, subprocess, sys, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('changed pinned input: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):
    q=Path(str(p)+'.tmp');q.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');q.replace(p)
def load(p):
    spec=importlib.util.spec_from_file_location('archived_dual_protocol',verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def flat(a):return [x for row in a for x in row]
def maximum(a):return max(abs(x) for x in a)
def child_cpu():
    x=resource.getrusage(resource.RUSAGE_CHILDREN);return x.ru_utime+x.ru_stime

def execute(config_path,release_path,out):
    start=time.monotonic();cpu=time.process_time();initial_child=child_cpu();proc=None;native_cpu_start=None
    out.mkdir(parents=True,exist_ok=False)
    r={'status':'started','implementation':pin(__file__),'configuration':pin(config_path),'release':pin(release_path),
       'native_program_attempts':0,'native_program_invocations':0,'response_solves':0,'kernel_evaluations':0,
       'error':None,'host':platform.node(),'controller_pid':os.getpid()}
    save(out/'RESULT.json',r)
    try:
        cfg=json.loads(config_path.read_text());release=json.loads(release_path.read_text())
        assert cfg['implementation']==pin(__file__) and cfg['maximum_native_programs']==1 and cfg['threads']==1 and cfg['timeout_seconds']==300
        assert release['status']=='approved_single_archived_dual_force_replay' and release['configuration']==pin(config_path)
        assert release['implementation']==cfg['implementation'] and release['output']==str(out.resolve())
        assert os.environ.get('SLURM_JOB_ID') and int(os.environ['SLURM_CPUS_PER_TASK'])==1
        job=os.environ['SLURM_JOB_ID'];p=subprocess.run(['/usr/local/bin/scontrol','show','job','-o',job],capture_output=True,text=True,check=True)
        (out/'ALLOCATION.txt').write_text(p.stdout);allocation=dict(v.split('=',1) for v in p.stdout.split() if '=' in v)
        assert allocation['QOS']=='standard_expand' and allocation['Partition']=='standard'
        assert allocation['NumNodes']=='1' and allocation['NumCPUs']=='1' and allocation['NumTasks']=='1'
        r['allocation']={'job_id':job,'source':pin(out/'ALLOCATION.txt'),'fields':allocation}
        P=load(cfg['protocol']);dense=read(cfg['dense_state']);fixture=read(cfg['fixture']);panel=read(cfg['panel'])
        slot=panel['states'][cfg['state_id']];op=fixture['native_operator'];build=read(cfg['build']);ready=read(cfg['bridge_ready'])
        assert slot['case_id']=='matched_dry_acetate_mono_central_La' and slot['direction'] is None and slot['coefficient_mode']=='zero'
        assert slot['parameters']=={'pairs':{},'triples':[]} and slot==dense['state'] and slot['fixture']==cfg['fixture']
        assert dense['status']=='failed' and dense['error']=="ValueError('available zero regression channel failed')"
        assert dense['returncode']==0 and dense['native_program_invocations']==1 and not dense['zero_regression']['gradient']['passed']
        assert all(row['passed'] for row in dense['zero_regression']['available_channels'].values())
        assert build['binary']==cfg['native_binary']==dense['native_binary'] and ready['root_build']==cfg['build']
        assert ready['sources']['BRIDGE.f90']==cfg['bridge_source'] and build['recipe']==ready['sources']['BUILD_RECIPE.json']
        assert cfg['native_binary']['sha256']=='68a47932f106d7d40bedaac1ddc034f6a314df182f8db07f318b1c9742f05e6a'
        verify(cfg['bridge_source']);verify(cfg['native_binary']);verify(dense['force_bridge_result']);verify(dense['supplied_uv'])
        for p in (fixture['native_qualification']['producer'],fixture['native_qualification']['review']):verify(p)
        u=op['u_native_eA'];v=op['v_native_eA'];n=len(slot['geometry']['xyz_A']);dim=3*n
        assert len(u)==len(v)==dim and all(math.isfinite(x) for x in u+v)
        # These are the fixture's admitted raw native vectors, not dense solutions.
        selected=[s for s in fixture['native_qualification']['snapshots'] if s['source']==op['snapshot']]
        assert len(selected)==1; snapshot=selected[0]
        assert [row[2] for row in snapshot['vectors']]==u and [row[3] for row in snapshot['vectors']]==v
        verify(snapshot['source'])
        gp=dense['zero_regression']['gradient']['source'];assert gp==fixture['baseline_observables']['native_outputs']['gradient_output']
        tokens=[];old_gradient=[]
        for line in verify(gp).read_text().splitlines():
            t=line.split()
            if len(t)>=5 and t[0]=='Anlyt' and t[1].isdigit():
                assert int(t[1])==len(old_gradient)+1;tokens.extend(t[2:5]);old_gradient.append([float(x.replace('D','E')) for x in t[2:5]])
        assert len(old_gradient)==n and old_gradient==[[-x for x in row] for row in fixture['baseline_observables']['forces_kcal_mol_A']]
        old_baseline=P.baseline(verify(dense['baseline_state']),slot,fixture)
        for name,p in slot['inputs'].items():shutil.copyfile(verify(p),out/{'xyz':'system.xyz','key':'system.key','parameters':'amoeba09.prm'}[name])
        r.update(state=slot,fixture=cfg['fixture'],dense_state=cfg['dense_state'],native_binary=cfg['native_binary'],
                 archived_dual_authority=fixture['native_qualification'],archived_gradient_source=gp,
                 original_failed_regression=dense['zero_regression'],u=u,v=v)
        env=os.environ.copy();env.pop('LD_PRELOAD',None);env.pop('LD_LIBRARY_PATH',None)
        env.update({k:'1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')})
        argv=[str(verify(cfg['native_binary'])),'system.xyz'];r.update(argv=argv,native_program_attempts=1);save(out/'RESULT.json',r)
        native_cpu_start=child_cpu();deadline=time.monotonic()+cfg['timeout_seconds']
        with (out/'stderr.txt').open('wb') as err,(out/'stdout.txt').open('wb') as log:
            proc=subprocess.Popen(argv,cwd=out,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,bufsize=0,start_new_session=True)
            r.update(native_program_invocations=1,native_pid=proc.pid);save(out/'RESULT.json',r)
            pending=b'';handshake=False
            with selectors.DefaultSelector() as selector:
                selector.register(proc.stdout,selectors.EVENT_READ)
                while not handshake:
                    wait=deadline-time.monotonic()
                    if wait<=0 or not selector.select(wait):raise TimeoutError('native baseline handshake timeout')
                    chunk=os.read(proc.stdout.fileno(),65536)
                    if not chunk:raise RuntimeError('native exited before SRMB_READY')
                    log.write(chunk);log.flush();pending+=chunk
                    while b'\n' in pending:
                        line,pending=pending.split(b'\n',1)
                        if line.strip()==b'SRMB_READY':handshake=True
                        if line.startswith(b'SRMB_ERROR'):raise RuntimeError(line.decode())
            b=P.baseline(out/'baseline_state.dat',slot,fixture)
            assert b==old_baseline,'fresh baseline differs from saved dense-state baseline'
            uv=out/'archived_uv.dat';P.write_uv(uv,b,u,v);r.update(baseline_state=pin(out/'baseline_state.dat'),archived_uv=pin(uv),baseline_exactly_reused=True);save(out/'RESULT.json',r)
            proc.stdin.write((str(uv.resolve())+'\n').encode());proc.stdin.flush();proc.stdin.close();proc.stdin=None
            tail,_=proc.communicate(timeout=max(.001,deadline-time.monotonic()));log.write(tail);log.flush();r['returncode']=proc.returncode
            if proc.returncode!=0:raise RuntimeError('native contraction failed')
        native=P.result(out/'force_bridge_result.dat',b,u,v);g=native['baseline_gradient'];gd=dense['gradient_kcal_mol_A'];go=old_gradient
        assert dense['local_delta_V_kcal_mol']==0 and all(x==0 for x in flat(dense['local_fixed_uv_gradient']))
        replay_error=[a-bb for a,bb in zip(flat(g),flat(go))]
        dense_shift=[a-bb for a,bb in zip(flat(gd),flat(g))]
        dense_error=[a-bb for a,bb in zip(flat(gd),flat(go))]
        closure=[a-bb-cc for a,bb,cc in zip(dense_error,dense_shift,replay_error)]
        # Use the failed run's original literal bound, without widening or retuning.
        limit=dense['zero_regression']['gradient']['limit'];assert limit==5.060311087690934e-09
        checks={'archived_dual_vs_archived_print':{'maximum_absolute_error':maximum(replay_error),'original_limit':limit,'passed':maximum(replay_error)<=limit},
                'saved_dense_vs_replay':{'maximum_absolute_difference':maximum(dense_shift),'full_flat_difference':dense_shift},
                'saved_dense_vs_archived_print':{'maximum_absolute_error':maximum(dense_error),'original_limit':limit,'passed':maximum(dense_error)<=limit},
                'difference_decomposition':{'maximum_absolute_closure_error':maximum(closure)},
                'units':'kcal/mol/angstrom'}
        r.update(status='archived_dual_replay_passed' if checks['archived_dual_vs_archived_print']['passed'] else 'archived_dual_replay_failed_gradient_gate',
                 native_result=pin(out/'force_bridge_result.dat'),native_readback=native,checks=checks,archived_gradient_tokens=tokens,
                 replay_gradient=g,saved_dense_gradient=gd,archived_print_gradient=go,force_scope='fixed_archived_dual_native_contraction',
                 counts=native['calls'],original_failed_state_preserved=True,completed_panel_states_repeated=0,
                 interpretation='Isolates contraction change from changing u/v; does not readmit the failed dense-to-iterative regression or qualify missing FD states.')
    except Exception as exc:
        r.update(status='failed',error=repr(exc))
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(proc.pid,15)
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,9);proc.wait()
            except Exception as stop:r['cleanup_error']=repr(stop)
        if proc is not None:r['returncode']=proc.returncode
    r['cost']={'wall_seconds':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-cpu,
               'all_child_CPU_seconds':child_cpu()-initial_child,'native_child_CPU_seconds':None if native_cpu_start is None else child_cpu()-native_cpu_start,
               'threads':1,'host':platform.node(),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),
               'peak_controller_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               'peak_child_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    save(out/'RESULT.json',r);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=execute(a.config,a.release,a.out);print(json.dumps({'status':r['status'],'error':r['error']}));sys.exit(0 if r['status']=='archived_dual_replay_passed' else 1)
