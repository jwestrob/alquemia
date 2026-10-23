"""One explicitly released initializer, parity program or gradient; no retries."""
import argparse
import shutil
from COMMON import *
from METADATA import sites,physical_metadata,key_settings

def inspect_initialization(out,row,m):
    reader=load(m['helpers']['snapshot'],'v2_snapshot');old=load(m['helpers']['pme_checks'],'v2_pme_checks')
    with (out/'native/pme_norm_export_001.dat').open() as f:snapshot=reader.read(f)
    expected,cell=sites(row['inputs']);site_checks=old.check_sites(snapshot,expected,cell)
    norms=reader.replay(snapshot);comparison=old.norm_comparison(snapshot,norms)
    require(comparison['passed'] and norms['q_negative_count']==0 and norms['q_min']>=0,'PME readback')
    state=identity_state(m);r=state.Rows(out/'native/certificate_actual_state.dat')
    try:
        header=state.read_header(r,{'cell':cell[:3],**{k:m['model'][k] for k in ('settings','coefficients','radial')}})
        require(header['ewald']==snapshot['ewald'] and header['volume']==snapshot['volume'],'Same initialized PME scalars')
        require(header['trans']==snapshot['trans'].tolist(),'Same initialized PME transform')
    finally:r.stream.close()
    return {'site_checks':site_checks,'summary_comparison':comparison,'operator_norms':norms,
            'V2_header':header,'full_body_admission_deferred_to_certificate':True,
            'e_app':None,'application_error_bound_established':False}

def bypass_inputs(row,m):
    center=m['rows']['central'];require(row['inputs']['xyz']['sha256']==center['inputs']['xyz']['sha256'] and row['inputs']['parameters']['sha256']==center['inputs']['parameters']['sha256'],'Same bypass geometry/parameters')
    def split(text):
        kept=[];response=[]
        for line in text.splitlines(keepends=True):
            words=line.split('#',1)[0].split()
            if words and words[0].lower() in ('polar-local-response','polar-radial-drive'):response.append(words)
            else:kept.append(line)
        return ''.join(kept),response
    text=verify(row['inputs']['key']).read_text();base,records=split(text)
    expected,_=split(verify(center['inputs']['key']).read_text());require(base==expected,'Only response-record bypass delta')
    key_settings(text)
    if row['role'] in ('BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT'):require(not records,'Absent response records')
    else:
        wanted={('polar-local-response',900,402):4,('polar-local-response',901,402):4,
                ('polar-radial-drive',900,402):1,('polar-radial-drive',901,402):1}
        seen=set()
        for v in records:
            key=(v[0].lower(),int(v[1]),int(v[2]));require(key in wanted and key not in seen,'Zero bypass record identity');seen.add(key)
            require(len(v)==3+wanted[key] and all(float(x)==0. for x in v[3:]),'Literal zero bypass')
        require(seen==set(wanted),'Both old/radial zero slots')

def run(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'commands':{},'native_attempts':0,'automatic_retry':False,'error':None}
    save(out/'RESULT.json',result)
    try:
        cfg,m,rel,alloc=source_context(config,release,out,__file__,'native',1)
        kind=cfg['operation'];require(kind in ('initialize','parity','gradient'),'Native operation')
        sid=cfg['state_id'];row=m['rows'][sid];build=native_build(cfg,m)
        require(rel['operation']==kind and rel['state_id']==sid and rel['maximum_native_programs']==1,'Exact one-native release')
        needs=policy.phase_needs(kind,sid,row['role'])
        cert=certificate(cfg,m) if needs['own_certificate'] else None
        if row['role']=='ACTIVE':physical_metadata(row,m)
        else:bypass_inputs(row,m)
        label='initialize' if kind=='initialize' else ('parity' if kind=='parity' else ('gradient_old' if row['role']=='BYPASS_OLD_ABSENT' else 'gradient_clone'))
        executable=build['executables'][label];native=out/'native';native.mkdir()
        for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]:
            shutil.copyfile(verify(row['inputs'][k]),native/name);require(pin(native/name)['sha256']==row['inputs'][k]['sha256'],'Staged input bytes')
        binding=pin(config)['sha256']
        if kind!='parity':(native/'stage2_ticket.dat').write_text(policy.ticket(kind,row['role'],sid,binding))
        if cert is not None:shutil.copyfile(verify(cert['expected_state']),native/'expected_state.dat')
        if kind=='parity':
            require(cfg['probe_source']==m['probe_source'] and cfg['probe_data']==m['probe_data'],'Exact ten-probe binding')
            verify(cfg['probe_source']);shutil.copyfile(verify(cfg['probe_data']),native/'probes.dat')
        result.update(manifest=cfg['manifest'],state_id=sid,role=row['role'],source_inputs=row['source_inputs'],inputs=row['inputs'],input_binding=binding,build=cfg['build'],source_manifest=m['native_sources'],allocation=alloc,
                      certificate=cfg.get('certificate') if cert is not None else None,certificate_admission=cfg.get('certificate_admission'),expected_state=None if cert is None else cert['expected_state'],
                      campaign_identity=policy.campaign_identity(m,cfg['build']),ticket=pin(native/'stage2_ticket.dat') if kind!='parity' else None)
        verify(executable);loader=load(m['helpers']['slurm1'],'v2_native_loader')
        ld=subprocess.run([m['loader_command']['path'],executable['path']],check=True,capture_output=True,text=True,timeout=30)
        require(not ld.stderr.strip(),'Native loader stderr');result['native_runtime']=loader.native_loader(ld.stdout,m['native_runtime'])
        result['native_attempts']=1;save(out/'RESULT.json',result)
        invoke(result,out,'native',[executable['path'],'system.xyz'],native,3600,3*1024**3,3600)
        for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]:require(pin(native/name)['sha256']==row['inputs'][k]['sha256'],'Native mutated input')
        if kind=='initialize':
            receipt=(native/'stage2_initialization_receipt.dat').read_text().splitlines()
            require(receipt==['STAGE2_INITIALIZATION_V1','COUNTS direct=1 mutual=0 solve=0 energy=0 gradient=0 factor=0','STATE_EXPORTED_NOT_CERTIFIED','END'],'Initialization count receipt')
            state_pin=pin(native/'certificate_actual_state.dat');pme_pin=pin(native/'pme_norm_export_001.dat')
            result.update(state_output=state_pin,pme_export=pme_pin,initialization_receipt=pin(native/'stage2_initialization_receipt.dat'),
                identity=policy.full_identity(m,cfg['build'],sid,state_pin,pme_pin),counts={'native_programs':1,'direct_initializations':1,'mutual_applications':0,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0})
            set_limits(3*1024**3,1800)
            def timeout(signum,frame):raise TimeoutError('Saved readback exceeded1800seconds')
            signal.signal(signal.SIGALRM,timeout);signal.alarm(1800)
            try:result['PME_checks']=inspect_initialization(out,row,m)
            finally:signal.alarm(0)
            result['status']='stage2_initialization_complete'
        elif kind=='parity':
            actual=pin(native/'parity_state.dat');require(actual['sha256']==cert['expected_state']['sha256'],'Exact V2 state before probes')
            require(pin(native/'pme_norm_export_001.dat')['sha256']==cert['pme_export']['sha256'],'Same initialized PME')
            replay=load(m['execution_sources']['PARITY_REPLAY.py'],'v2_parity_replay')
            set_limits(3*1024**3,1800)
            def timeout(signum,frame):raise TimeoutError('Saved parity replay exceeded1800seconds')
            signal.signal(signal.SIGALRM,timeout);signal.alarm(1800)
            try:record=replay.replay(cfg,m,cert,native)
            finally:signal.alarm(0)
            save(out/'REPLAY.json',record)
            require(record['status']=='passed_finite_complete_periodic_operator_parity','Actual parity checks')
            result.update(status='periodic_v2_parity_complete',identity=cert['identity'],actual_state=actual,replay=pin(out/'REPLAY.json'),
                          counts={'native_programs':1,'direct_initializations':1,'complete_operator':30,'separate_reciprocal':30,'total_reciprocal':60,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0})
        else:
            parser=load(m['execution_sources']['GRADIENT_PARSE.py'],'v2_gradient_parser')
            parsed=parser.parse_gradient(native/'stage2_gradient.dat',row['role'],m,n=11958)
            if cert is not None:
                actual=pin(native/'certificate_actual_state.dat');require(actual['sha256']==cert['expected_state']['sha256'],'Exact state before solve')
                require(parsed['radial']==m['model']['radial'],'Observed both radial slots')
                require(parsed['shared_h']['passed'] and parsed['polarization_contraction']['passed'],'Actual shared-h/dual contraction gates')
                result.update(actual_state=actual,exact_state_match=True,identity=cert['identity'])
            result.update(status='stage2_gradient_complete',gradient_output=pin(native/'stage2_gradient.dat'),counts={'native_programs':1,'gradient_calls':1,'energy_only_calls':0,'extra_operator_applications':0})
    except BaseException as exc:result.update(status='periodic_v2_native_failed',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'controller_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'scheduler_cost_separate':True}
    save(out/'RESULT.json',result)
    return 0 if result['status'] in ('stage2_initialization_complete','stage2_gradient_complete','periodic_v2_parity_complete') else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(run(a.config,a.release,a.out))
