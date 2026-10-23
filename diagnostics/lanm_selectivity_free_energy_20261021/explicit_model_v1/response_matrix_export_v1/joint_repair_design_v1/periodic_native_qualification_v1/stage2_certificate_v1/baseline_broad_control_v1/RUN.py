"""Exactly one bound initializer or gradient, followed only by saved readback."""
import argparse
import math
import shutil
from BINDING import *

def sites(inputs):
    xyz=verify(inputs['xyz']).read_text().splitlines();key=verify(inputs['key']).read_text();prm=verify(inputs['parameters']).read_text()
    require(len(xyz)==11960 and int(xyz[0].split()[0])==11958,'literal atom count')
    cell=list(map(float,xyz[1].split()));require(len(cell)==6 and cell[3:]==[90.]*3 and min(cell[:3])>=44,'orthogonal declared cell')
    alpha={}
    for text in (prm,key):
        for line in text.splitlines():
            v=line.split()
            if v and v[0].lower()=='polarize':alpha[int(v[1])]=float(v[2])
    output=[]
    for i,line in enumerate(xyz[2:],1):
        v=line.split();t=int(v[5]);r=list(map(float,v[2:5]));a=alpha[t]
        require(int(v[0])==i and t in (900,402,403,405,414) and a>=1e-8 and all(math.isfinite(x) for x in [a,*r]),'native input site')
        output.append({'atom':i,'type':t,'xyz':r,'alpha':a})
    return output,cell

def inspect_initialization(out,row,m):
    reader=load(m['helpers']['snapshot'],'stage2_snapshot');old=load(m['helpers']['pme_checks'],'stage2_pme_checks')
    with (out/'native/pme_norm_export_001.dat').open() as f:snapshot=reader.read(f)
    expected,cell=sites(row['inputs']);site_checks=old.check_sites(snapshot,expected,cell)
    norms=reader.replay(snapshot);comparison=old.norm_comparison(snapshot,norms)
    require(comparison['passed'] and norms['q_negative_count']==0 and norms['q_min']>=0,'PME bounds/readback')
    return {'site_checks':site_checks,'summary_comparison':comparison,'operator_norms':norms,'e_app':None,'application_error_bound_established':False}

def run(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'commands':{},'native_attempts':0,'automatic_retry':False,'error':None}
    save(out/'RESULT.json',result)
    try:
        cfg,m,rel,alloc=source_context(config,release,out,__file__,'native',1)
        kind=cfg['operation'];require(kind in ('initialize','gradient'),'declared native kind')
        state_id=cfg['state_id'];row=m['rows'][state_id];build=native_build(cfg,m)
        require(rel['operation']==kind and rel['state_id']==state_id and rel['maximum_native_programs']==1,'root finite row release')
        require(kind=='gradient' or (row['role']=='ACTIVE' and state_id!='central'),'no repeat central initializer')
        cert=None
        if kind=='gradient' and row['role']=='ACTIVE':cert=certificate(cfg,m)
        executable=build['executables']['initialize' if kind=='initialize' else ('gradient_old' if row['role']=='BYPASS_OLD_ABSENT' else 'gradient_clone')]
        native=out/'native';native.mkdir()
        for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]:
            shutil.copyfile(verify(row['inputs'][key]),native/name);require(pin(native/name)['sha256']==row['inputs'][key]['sha256'],'staged input bytes')
        binding=pin(config)['sha256'];role='INITIALIZE' if kind=='initialize' else row['role']
        (native/'stage2_ticket.dat').write_text('STAGE2_BOUND_INPUT_V1\n'+role+'\n'+state_id+'\n'+binding+'\nEND\n')
        if cert is not None:shutil.copyfile(verify(cert['expected_state']),native/'expected_state.dat')
        result.update(manifest=cfg['manifest'],state_id=state_id,role=row['role'],source_inputs=row['source_inputs'],inputs=row['inputs'],input_binding=binding,build=cfg['build'],source_manifest=m['native_sources'],allocation=alloc,certificate=cfg.get('certificate') if cert is not None else None,certificate_admission=cfg.get('certificate_admission'),expected_state=None if cert is None else cert['expected_state'],ticket=pin(native/'stage2_ticket.dat'))
        verify(executable)
        loader=load(m['helpers']['slurm1'],'stage2_native_loader')
        ld=subprocess.run([m['loader_command']['path'],executable['path']],check=True,capture_output=True,text=True,timeout=30)
        require(not ld.stderr.strip(),'native loader stderr')
        result['native_runtime']=loader.native_loader(ld.stdout,m['native_runtime'])
        result['native_attempts']=1;save(out/'RESULT.json',result)
        invoke(result,out,'native',[executable['path'],'system.xyz'],native,3600,3*1024**3,3600)
        for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]:require(pin(native/name)['sha256']==row['inputs'][key]['sha256'],'native mutated input')
        if kind=='initialize':
            receipt=(native/'stage2_initialization_receipt.dat').read_text().splitlines()
            require(receipt==['STAGE2_INITIALIZATION_V1','COUNTS direct=1 mutual=0 solve=0 energy=0 gradient=0 factor=0','STATE_EXPORTED_NOT_CERTIFIED','END'],'complete initialization receipt')
            result.update(status='stage2_initialization_complete',state_output=pin(native/'certificate_actual_state.dat'),pme_export=pin(native/'pme_norm_export_001.dat'),initialization_receipt=pin(native/'stage2_initialization_receipt.dat'),counts={'native_programs':1,'direct_initializations':1,'mutual_applications':0,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0})
            # Same-state replay only; it neither solves nor factors the operator.
            set_limits(3*1024**3,1800)
            def replay_timeout(signum,frame):raise TimeoutError('saved PME readback exceeded 1800 seconds')
            signal.signal(signal.SIGALRM,replay_timeout);signal.alarm(1800)
            try:result['PME_checks']=inspect_initialization(out,row,m)
            finally:signal.alarm(0)
        else:
            parser=load(m['collector'],'stage2_saved_gradient_parser')
            parsed=parser.parse_gradient(native/'stage2_gradient.dat',n=11958)
            require(parsed is not None,'complete native gradient grammar')
            if cert is not None:
                actual=pin(native/'certificate_actual_state.dat')
                require(actual['sha256']==cert['expected_state']['sha256'],'exact native state before solve')
                result['actual_state']=actual;result['exact_state_match']=True
            result.update(status='stage2_gradient_complete',gradient_output=pin(native/'stage2_gradient.dat'),counts={'native_programs':1,'gradient_calls':1,'energy_only_calls':0,'extra_operator_applications':0})
    except BaseException as exc:result.update(status='stage2_native_failed',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'controller_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'scheduler_cost_separate':True}
    save(out/'RESULT.json',result);return 0 if result['status'] in ('stage2_initialization_complete','stage2_gradient_complete') else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(run(a.config,a.release,a.out))
