#!/usr/bin/env python3
"""Finite native force bridge harness. Execution requires exact root release."""
import argparse,hashlib,importlib.util,json,math,os,resource,selectors,shutil,platform,subprocess,sys,time
from pathlib import Path
from decimal import Decimal
BOOT_WALL=time.monotonic()
BOOT_CPU=time.process_time()
import torch
from PREPARE import pin,read,verify,save
import PROTOCOL

def module(p):
    spec=importlib.util.spec_from_file_location('force_bridge_'+p['sha256'][:12],verify(p));m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
def tensor(x):return torch.tensor(x,dtype=torch.float64)
def arithmetic(dim,*values):return 128*dim*torch.finfo(torch.float64).eps*max(1.,*[abs(float(x)) for x in values])
def params(mod,slot):
    p=slot['parameters'];return mod.Parameters(pairs={k:mod.Pair(**v) for k,v in p['pairs'].items()},triples={tuple(r['kinds']):mod.Triple(**r['coefficients']) for r in p['triples']})
def fp(cfg):
    body={k:v for k,v in cfg.items() if k!='root_release'};return hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def cpu_children():
    r=resource.getrusage(resource.RUSAGE_CHILDREN);return r.ru_utime+r.ru_stime
def half_token(s):return float(Decimal(5).scaleb(Decimal(s.replace('D','E')).as_tuple().exponent-1))
def zero_regression(slot,fixture,b,out,E,G):
    if slot['coefficient_mode']!='zero' or slot['direction'] is not None:return None
    op=fixture['native_operator'];q=fixture['native_qualification'];n=b['n'];dim=3*n;bounds=q['snapshots'];u=out['u'];v=out['v']
    ub=max(x['numerical_bounds']['native_vs_dense_u'] for x in bounds);vb=max(x['numerical_bounds']['native_vs_dense_v'] for x in bounds)
    pe=q['polarization_energy'];oldE=fixture['baseline_observables']['energy_kcal_mol'];eprint=5e-9 # inherited analyze total energy is printed to eight decimal places
    elim=eprint+pe['print_half_unit']+.5*b['K']*float(torch.linalg.vector_norm(tensor(op['Ep_e_per_A2'])))*ub+arithmetic(dim,E,oldE)
    rows={'u':{'error':float(torch.linalg.vector_norm(u-tensor(op['u_native_eA']))),'limit':ub},'v':{'error':float(torch.linalg.vector_norm(v-tensor(op['v_native_eA']))),'limit':vb},'energy':{'error':abs(E-oldE),'limit':elim}}
    obs=fixture['baseline_observables'];oldforce=obs.get('forces_kcal_mol_A')
    if oldforce is None:gradient={'status':'unavailable','reason':'No archived full gradient for this exact water fragment; FD check remains required.'}
    else:
        gp=obs['native_outputs']['gradient_output'];text=verify(gp).read_text();tokens=[]
        for line in text.splitlines():
            t=line.split()
            if len(t)>=5 and t[0]=='Anlyt' and t[1].isdigit():tokens.extend(t[2:5])
        # Parent testgrad parser authority may use a different row tag; do not invent precision.
        if len(tokens)!=dim:raise ValueError('archived gradient printing grammar needs source-bound adapter')
        pbound=max(half_token(s) for s in tokens);old=-tensor(oldforce);err=float(torch.max(torch.abs(G-old)));limit=pbound+arithmetic(dim,*G.flatten(),*old.flatten());gradient={'status':'checked','error':err,'limit':limit,'passed':err<=limit,'source':gp,'CG_force_error_propagation_available':False}
    for row in rows.values():row['passed']=row['error']<=row['limit']
    return {'available_channels':rows,'gradient':gradient,'passed':all(x['passed'] for x in rows.values()) and gradient.get('passed',True),'missing_regression_channels':['gradient'] if oldforce is None else []}

def run_state(cfg,panel,slot,mod,combiner,folder):
    folder.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();ccpu=cpu_children();fixture=read(slot['fixture']);proc=None
    r={'status':'started','id':slot['id'],'state':slot,'fixture':slot['fixture'],'kernel':cfg['kernel'],'native_binary':cfg['native_binary'],'native_program_invocations':0,'native_program_attempts':0,'error':None}
    for name,src in slot['inputs'].items():shutil.copyfile(verify(src),folder/({'xyz':'system.xyz','key':'system.key','parameters':'amoeba09.prm'}[name]))
    save(folder/'RESULT.json',r)
    try:
        env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
        argv=[str(verify(cfg['native_binary'])),'system.xyz'];r['argv']=argv;r['native_program_attempts']=1;save(folder/'RESULT.json',r)
        with (folder/'stderr.txt').open('w') as err,(folder/'stdout.txt').open('w') as log:
            proc=subprocess.Popen(argv,cwd=folder,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,bufsize=0,start_new_session=True)
            r['native_program_invocations']=1;r['pid']=proc.pid;save(folder/'RESULT.json',r)
            selector=selectors.DefaultSelector();selector.register(proc.stdout,selectors.EVENT_READ);deadline=time.monotonic()+cfg['timeout_seconds'];ready=False;buffer=b''
            while not ready:
                wait=deadline-time.monotonic()
                if wait<=0 or not selector.select(wait):raise TimeoutError('native baseline handshake timeout')
                chunk=os.read(proc.stdout.fileno(),65536)
                if not chunk:raise RuntimeError('native exited before SRMB_READY')
                log.write(chunk.decode('utf-8'));log.flush();buffer+=chunk
                while b'\n' in buffer:
                    line,buffer=buffer.split(b'\n',1)
                    if line.strip()==b'SRMB_READY':ready=True
                    if line.startswith(b'SRMB_ERROR'):raise RuntimeError(line.decode())
            selector.close()
            b=PROTOCOL.baseline(folder/'baseline_state.dat',slot,fixture);r['baseline_state']=pin(folder/'baseline_state.dat');n=b['n'];dim=3*n;xyz=tensor(slot['geometry']['xyz_A']).requires_grad_();Araw=tensor(b['A']);ed=tensor(b['Ed']);ep=tensor(b['Ep'])
            assert b['debye']==cfg['units']['debye_per_eA'] and b['bohr']==cfg['units']['bohr_A']
            bound=arithmetic(dim,*Araw.flatten(),*ed,*ep);symerr=float(torch.max(torch.abs(Araw-Araw.T)));assert abs(symerr-b['reported_symmetry_error'])<=bound and b['reported_dual_operator_error']<=bound
            A,receipt=mod.symmetrize_native_operator(Araw,max_absolute_asymmetry=bound,admission_authority={'fresh_baseline':r['baseline_state'],'phase':panel['phase']})
            g=slot['geometry'];geometry=mod.Geometry(g['metal_atom'],tuple(g['atomic_numbers']),tuple(mod.Ligand.from_dict(l) for l in g['ligands']));parameters=params(mod,slot);base=mod.Baseline(A,ed,ep,tuple(s['atom']-1 for s in b['sites']),b['K'],None,tensor(b['permanent_mu_eA']),authority=panel['shared_baseline']['sha256'])
            out=mod.evaluate(xyz,geometry,parameters,base);u=out['u'].detach();v=out['v'].detach();delta=out['delta'];uvpath=folder/'supplied_uv.dat';PROTOCOL.write_uv(uvpath,b,u.tolist(),v.tolist());r['supplied_uv']=pin(uvpath)
            proc.stdin.write((str(uvpath.resolve())+'\n').encode());proc.stdin.flush();proc.stdin.close();proc.stdin=None
            remaining=deadline-time.monotonic();tail,_=proc.communicate(timeout=max(.001,remaining));log.write(tail.decode('utf-8'));log.flush();r['returncode']=proc.returncode
            if proc.returncode!=0:raise RuntimeError('native fixed-gradient process failed')
        native=PROTOCOL.result(folder/'force_bridge_result.dat',b,u.tolist(),v.tolist());r['force_bridge_result']=pin(folder/'force_bridge_result.dat');rawpolar=-.5*b['K']*float(u@ep);assert abs(native['energy']['RAW_REDUCED_POLARIZATION']-rawpolar)<=arithmetic(dim,rawpolar)
        raw=native['energy'];assert abs(raw['OTHER_ENERGY']+raw['RAW_REDUCED_POLARIZATION']-raw['RAW_NATIVE_ENERGY'])<=arithmetic(dim,*raw.values())
        combined=combiner.combine(mod,xyz,geometry,parameters,out,{'other_energy':raw['OTHER_ENERGY'],'baseline_gradient':native['baseline_gradient'],'supplied_u':native['u'],'supplied_v':native['v'],'site_atom_indices':base.site_atom_indices},b['K']);E=combined['energy_kcal'];G=combined['gradient_kcal_per_A'].detach();localG=combined['local_explicit_gradient_kcal_per_A'].detach();assert torch.isfinite(G).all() and math.isfinite(E)
        regression=zero_regression(slot,fixture,b,out,E,G)
        r.update(status='native_force_state_complete',energy_kcal_mol=E,gradient_kcal_mol_A=G.tolist(),forces_kcal_mol_A=(-G).tolist(),native_gradient=native['baseline_gradient'],local_fixed_uv_gradient=localG.tolist(),combination_source=cfg['combiner'],local_delta_V_kcal_mol=float(delta['V'].detach()),corrected_polarization_energy_kcal_mol=float(out['polarization_energy'].detach()),other_energy_kcal_mol=raw['OTHER_ENERGY'],raw_native=native,site_order=b['sites'],direct_RHS=b['Ed'],energy_RHS=b['Ep'],u=u.tolist(),v=v.tolist(),native_flags=b['native_flags'],symmetry_projection=receipt,minimum_corrected_eigenvalue=float(out['eigenvalues'][0].detach()),solve_residual=float(out['solve_residual'].detach()),zero_regression=regression,force_scope='native_full_baseline_derivative_plus_local_fixed_dual_derivative',counts=native['calls'])
        if regression and not regression['passed']:raise ValueError('available zero regression channel failed')
    except Exception as exc:
        r['status']='failed';r['error']=repr(exc)
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid,15)
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,9);proc.wait()
        if proc is not None:r['returncode']=proc.returncode
    r['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'native_child_CPU_seconds':cpu_children()-ccpu,'native_threads':1,'controller_threads':1,'host':platform.node(),'peak_process_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_child_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss};save(folder/'RESULT.json',r);return r

def gate(cfg,panel,out):
    release=read(cfg['root_release'])
    assert release['status']=='approved_native_force_bridge_finite_run' and release['configuration_fingerprint']==fp(cfg) and release['output']==str(out.resolve()) and release['mode']==cfg['mode']
    for key in ('implementation','protocol','combiner','source_review','build','native_binary','panel'):assert release[key]==cfg[key]
    assert cfg['implementation']==pin(__file__) and cfg['protocol']==pin(Path(__file__).with_name('PROTOCOL.py')) and cfg['metadata_helper']==pin(Path(__file__).with_name('PREPARE.py'))
    review=read(cfg['source_review']);assert review['status']=='passed_native_force_qualification_source_review'
    for key in ('implementation','protocol','combiner','panel','bridge_ready'):assert review[key]==cfg[key]
    bridge=read(cfg['bridge_ready']);build=read(cfg['build'])
    assert bridge['root_build']==cfg['build'] and bridge['kernel']==cfg['kernel']==panel['kernel'] and bridge['phase']==panel['phase']
    assert bridge['sources']['COMBINE.py']==cfg['combiner'] and build['binary']==cfg['native_binary'] and build['recipe']==bridge['sources']['BUILD_RECIPE.json']
    assert build['status']=='built_isolated_native_force_bridge' and all(c['returncode']==0 for c in build['commands'])
    verify(cfg['native_binary'])
    assert cfg['mode'] in ('scout','remaining') and cfg['threads']==1 and cfg['maximum_native_programs']==(6 if cfg['mode']=='scout' else 96) and cfg['timeout_seconds']==300
    assert len(panel['states'])==102 and len(panel['scout_ids'])==6 and len(panel['remaining_ids'])==96 and set(panel['scout_ids']).isdisjoint(panel['remaining_ids'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();cfg=json.loads(a.config.read_text());panel=read(cfg['panel']);gate(cfg,panel,a.out)
    torch.set_num_threads(1);torch.set_num_interop_threads(1);a.out.mkdir(parents=True,exist_ok=False);t=time.monotonic();cpu=time.process_time();mod=module(cfg['kernel']);combiner=module(cfg['combiner']);ids=panel['scout_ids'] if cfg['mode']=='scout' else panel['remaining_ids'];reuse={}
    if cfg['mode']=='remaining':
        old=read(cfg['scout_result']);admission=read(cfg['scout_admission'])
        assert old['status']=='native_force_scout_complete' and old['panel']==cfg['panel'] and old['native_binary']==cfg['native_binary'] and old['kernel']==cfg['kernel'] and set(old['rows'])==set(panel['scout_ids'])
        assert admission['status']=='approved_native_force_scout_for_remaining_panel' and admission['result']==cfg['scout_result'] and admission['panel']==cfg['panel'] and admission['native_binary']==cfg['native_binary'];reuse=old['rows']
    result={'status':'started','implementation':pin(__file__),'configuration':pin(a.config),'configuration_fingerprint':fp(cfg),'panel':cfg['panel'],'native_binary':cfg['native_binary'],'kernel':cfg['kernel'],'combiner':cfg['combiner'],'build':cfg['build'],'source_review':cfg['source_review'],'bridge_ready':cfg['bridge_ready'],'root_release':cfg['root_release'],'mode':cfg['mode'],'rows':{},'reused_scout_rows':reuse,'native_program_invocations':0,'native_program_attempts':0,'error':None};save(a.out/'RESULT.json',result)
    try:
        for id in ids:
            assert result['native_program_attempts']<cfg['maximum_native_programs'];r=run_state(cfg,panel,panel['states'][id],mod,combiner,a.out/'states'/id)
            result['native_program_invocations']+=r['native_program_invocations'];result['native_program_attempts']+=r['native_program_attempts'];result['rows'][id]=pin(a.out/'states'/id/'RESULT.json');save(a.out/'RESULT.json',result)
            if r['status']!='native_force_state_complete':result['error']={'state':id,'reason':r['error']};break
    except Exception as exc:result['error']={'harness':repr(exc)}
    result['status']=('native_force_scout_complete' if cfg['mode']=='scout' else 'native_force_remaining_complete') if result['error'] is None else 'failed'
    result['cost']={'wall_seconds':time.monotonic()-BOOT_WALL,'controller_CPU_seconds':time.process_time()-BOOT_CPU,'run_wall_seconds_excluding_imports':time.monotonic()-t,'run_controller_CPU_seconds_excluding_imports':time.process_time()-cpu,'native_child_CPU_seconds':sum(read(p)['cost']['native_child_CPU_seconds'] for p in result['rows'].values()),'threads':1,'host':platform.node(),'peak_process_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_child_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    save(a.out/'RESULT.json',result);print(json.dumps({'status':result['status'],'new_native_programs':result['native_program_invocations'],'error':result['error']}));return 0 if result['error'] is None else 1
if __name__=='__main__':raise SystemExit(main())
