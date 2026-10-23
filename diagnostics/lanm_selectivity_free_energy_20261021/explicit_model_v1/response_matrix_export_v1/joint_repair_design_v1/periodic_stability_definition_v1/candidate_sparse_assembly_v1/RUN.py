"""Pinned candidate-specific real operator; prepare is metadata-only."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import signal
import socket
import subprocess
import sys
import time
import traceback

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def pin(p):
    p=Path(p).resolve(); h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return dict(path=str(p),sha256=h.hexdigest())

def verify(p):
    assert pin(p['path'])==p,'changed pin: '+p['path']
    return Path(p['path'])

def read(p): return json.loads(verify(p).read_text())

def write(p,obj):
    p=Path(p); tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n'); tmp.replace(p)

def load(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p)); m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def unity(scales):
    assert set(scales)=={'mutual-'+s+'-scale' for s in ('11','12','13','14')}
    assert all(x==1. for x in scales.values()),'non-unity mutual scaling unsupported'

def weight(r):
    assert math.isfinite(r) and r>=0
    if r>=21.: return 0.
    fzero=1./((1.+5.+25./3.)*math.exp(-5.))**2
    t=2.*r;p=1.+t+t*t/3.; f=fzero*(p*math.exp(-t))**2
    if r>20.:
        u=r-20.;f*=((1.-u)**3*(1.+3.*u+6.*u*u))
    return f

def environment(weights,c,d):
    s=0.;q=0.
    for f in weights:
        assert math.isfinite(f) and f>=0
        q=q+2.*f*s;s=s+f
    eta=q/(1.+s*s);g=2.*s/(1.+s);sigma=g*(c+d*eta)
    return dict(s=s,cross_sum=q,eta=eta,g=g,sigma=sigma,normalized_beta=math.exp(sigma)-1.)

def metadata(inputs):
    for p in inputs['sources'].values(): verify(p)
    t=inputs['template'];comp=read(t['composition_manifest'])
    assert comp['status']=='prepared_only_not_executable'
    assert comp['candidate']['selected']==t['candidate']
    for p in comp['candidate'].values(): verify(p)
    verify(comp['configuration']);verify(comp['source_review'])
    f=next(x['prepared'] for x in comp['files'] if x['metal']=='La')
    for p in f.values(): verify(p)
    assert f['xyz']==t['geometry'] and t['metal']=='La'
    result=read(t['pme_result']);root=read(t['pme_root_review']);ind=read(t['pme_bound_review'])
    assert result['status']=='full_box_operator_initialization_complete' and result['error'] is None
    assert root['status']=='admitted_full_box_operator_initialization_nominal_data'
    assert root['execution']==t['pme_result'] and root['independent_actual_review']==t['pme_bound_review']
    assert root['composition']==t['composition_manifest']
    cell=result['cells']['La'];assert cell['export']==t['pme_export']
    assert cell['site_checks']['passed'] and cell['summary_comparison']['passed']
    for p in cell['inputs'].values(): verify(p)
    assert cell['inputs']['system.xyz']['sha256']==f['xyz']['sha256']
    assert cell['inputs']['amoebabio09.prm']['sha256']==f['parameters']['sha256']
    original=verify(f['key']).read_text().splitlines()
    stripped=[x for x in original if not x.lower().startswith('polar-local-response ')]
    nativekey=verify(cell['inputs']['system.key']).read_text().splitlines()
    assert nativekey[:len(stripped)]==stripped
    settings={}
    for line in nativekey[len(stripped):]:
        a=line.split('#',1)[0].split()
        if a:assert a[0].lower() not in settings;settings[a[0].lower()]=a[1:]
    for k,v in {'ewald':[],'ewald-cutoff':['9.0'],'ewald-alpha':['0.4'],'pewald-alpha':['0.4'],
                'pme-grid':['64','64','64'],'pme-order':['5'],'ppme-order':['5'],'fft-package':['FFTW']}.items():
        assert settings[k]==v
    assert t['U_R']==cell['operator_norms']['U_R_nominal']
    assert t['delta_odd']==cell['operator_norms']['evenness_perturbation_nominal']
    assert t['e_app'] is None and t['application_error_bound_established'] is False
    read(t['nominal_contract']);verify(t['pme_export'])
    polar={};scales={};pairs={};response={}
    for source in ('parameters','key'):
        seen=set()
        for line in verify(f[source]).read_text().splitlines():
            a=line.split('#',1)[0].split()
            if not a: continue
            name=a[0].lower()
            if name=='polarize' and int(a[1]) in (402,403,405,414,900,901):
                k=(name,int(a[1]));assert k not in seen;seen.add(k)
                polar[str(k[1])]=dict(alpha=float(a[2]),thole=float(a[3]))
            elif name.startswith('mutual-') and name.endswith('-scale'):
                assert name not in seen;seen.add(name);scales[name]=float(a[1])
            elif name=='polpair':
                k=tuple(sorted(map(int,a[1:3])));assert (name,k) not in seen;seen.add((name,k));pairs[k]=float(a[3])
            elif name=='polar-local-response':
                k=tuple(map(int,a[1:3]));assert k not in response;response[k]=list(map(float,a[3:]))
            elif name=='polar-donor-response': raise AssertionError('donor response is outside scope')
    unity(scales)
    assert set(polar)=={'402','403','405','414','900','901'}
    assert all(math.isfinite(v) for x in polar.values() for v in x.values())
    assert all(x['alpha']>0 and x['thole']>=0 for x in polar.values())
    assert set(pairs)=={(402,900),(402,901)} and set(response)=={(900,402),(901,402)}
    for m,ty in [('La',900),('Dy',901)]:
        effective=comp['effective_a_b_c_d_R_epsilon_POLPAIR'][m]
        assert response[(ty,402)]==effective[:4] and pairs[(402,ty)]==effective[6]
    lines=verify(f['xyz']).read_text().splitlines()
    assert len(lines)==11960 and int(lines[0].split()[0])==11958
    assert list(map(float,lines[1].split()))==[49.324]*3+[90.]*3
    rows=[]
    for i,line in enumerate(lines[2:]):
        a=line.split(); assert int(a[0])==i+1
        rows.append(dict(type=int(a[5]),xyz=list(map(float,a[2:5])),bonds=list(map(int,a[6:]))))
    assert Counter(x['type'] for x in rows)=={900:1,402:3980,403:7960,405:7,414:10}
    for i,x in enumerate(rows,1):
        assert all(math.isfinite(v) for v in x['xyz'])
        if x['type']==402:
            assert len(x['bonds'])==2 and len(set(x['bonds']))==2
            assert all(1<=j<=len(rows) and rows[j-1]['type']==403 and rows[j-1]['bonds']==[i] for j in x['bonds'])
        elif x['type']==403:
            assert len(x['bonds'])==1 and rows[x['bonds'][0]-1]['type']==402 and i in rows[x['bonds'][0]-1]['bonds']
        else: assert not x['bonds']
    # The admitted raw initialization fixes the actual positive-alpha DOF order.
    count=0;grid=None;ewald=None
    with verify(t['pme_export']).open() as raw:
        for line in raw:
            a=line.split()
            if not a:continue
            if a[0]=='GRID':grid=list(map(int,a[1:]))
            elif a[0]=='EWALD_VOLUME':ewald=float(a[1])
            elif a[0]=='SITE':
                count+=1; assert int(a[1])==count
                x=rows[count-1]; assert int(a[2])==x['type']
                assert float(a[6])==polar[str(x['type'])]['alpha'] and list(map(float,a[7:10]))==x['xyz']
    assert count==11958 and grid==[1,64,64,64,5,11958,11958,262144] and ewald==.4
    graph=read(inputs['graph_result'])
    assert graph['status']=='completed_single_saved_graph_symbolic_check'
    for k in ('indices','indptr'): assert graph['assembly'][k+'_sha256']==t[k+'_payload_sha256']
    assert graph['assembly']['info']['lower_entries']==t['nnz'] and t['shape']==[35874,35874]
    return dict(status='prepared_candidate_sparse_assembly_metadata',candidate=t['candidate'],geometry=t['geometry'],
        composition_manifest=t['composition_manifest'],type_parameters=polar,mutual_scales=scales,
        pair_overrides=[dict(types=list(k),gamma=v) for k,v in sorted(pairs.items())],
        response_coefficients=response[(900,402)],metal_site_index=next(i for i,x in enumerate(rows) if x['type']==900),
        cell_A=[49.324]*3,cutoff_A=9.,ewald_inverse_A=.4,sites=count,water_count=3980,
        pme_site_map_exact=True,H_evaluated=False,all_alpha_positive=True,source_inputs=inputs['sources'])

def prepare(inputs_path,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.perf_counter();cpu=time.process_time()
    inp=json.loads(inputs_path.read_text());meta=metadata(inp);write(out/'METADATA.json',meta)
    # Six declared artificial controls, no real geometry/assembly/backend call.
    checks=[abs(weight(2.5)-1.)<1e-15,weight(21.)==weight(22.)==0.]
    r=20.5;t=2*r;base=(1+t+t*t/3)**2*math.exp(-2*t)/((1+5+25/3)*math.exp(-5))**2
    checks.append(abs(weight(r)/base-.5)<1e-14)
    e=environment([1.,2.],.3,.4);checks.append(e['cross_sum']==4. and e['eta']==.4)
    checks.append(environment([1.,2.],0.,0.)['normalized_beta']==0.)
    bad={'mutual-'+s+'-scale':1. for s in ('11','12','13','14')};bad['mutual-12-scale']=.5
    try:unity(bad);checks.append(False)
    except AssertionError:checks.append(True)
    assert all(checks)
    write(out/'RESULT.json',dict(status='passed_candidate_sparse_assembly_preparation',implementation=pin(__file__),inputs=pin(inputs_path),
        metadata=pin(out/'METADATA.json'),synthetic_checks=checks,actual_assembly_calls=0,numeric_factorizations=0,native_calls=0,
        cost=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1)))

def context(config,release_path,out):
    cfg=json.loads(config.read_text());rel=json.loads(release_path.read_text())
    assert cfg['implementation']==pin(__file__) and rel['implementation']==pin(__file__)
    assert rel['status']=='released_one_terminal_candidate_La_real_assembly'
    assert rel['configuration']==pin(config) and Path(rel['output']).resolve()==out.resolve()
    review=read(rel['source_review']);assert review['status']=='passed_candidate_sparse_assembly_source_review'
    assert review['implementation']==pin(__file__) and review['configuration']==pin(config)
    assert rel['resources']==cfg['resources']
    r=cfg['resources'];assert r==dict(hostname=socket.gethostname(),address_space_bytes=8*1024**3,assembly_bytes=4*1024**3,wall_seconds=180,CPU_seconds=180,numerical_threads=1)
    assert all(os.environ.get(k)=='1' for k in THREADS) and not os.environ.get('SLURM_JOB_ID')
    env=read(cfg['environment']);assert Path(sys.prefix).resolve()==Path(env['prefix'])
    for p in env['runtime_files']:verify(p)
    prep=read(cfg['preparation']);assert prep['status']=='passed_candidate_sparse_assembly_preparation'
    assert prep['implementation']==pin(__file__) and prep['inputs']==cfg['inputs'] and prep['metadata']==cfg['metadata']
    return cfg,read(cfg['inputs'])

def worker(config,release_path,out):
    caps=json.loads(config.read_text())['resources'];limits={}
    for kind,k in [(resource.RLIMIT_AS,'address_space_bytes'),(resource.RLIMIT_CPU,'CPU_seconds')]:
        old=resource.getrlimit(kind);value=caps[k]
        assert all(v==resource.RLIM_INFINITY or value<=v for v in old)
        resource.setrlimit(kind,(value,value));limits[k]=list(resource.getrlimit(kind))
    cfg,inp=context(config,release_path,out);dest=out/'worker';dest.mkdir(exist_ok=False)
    tick=time.perf_counter();cpu=time.process_time()
    result=dict(status='started',phase='metadata',implementation=pin(__file__),configuration=pin(config),release=pin(release_path),
        assembly_attempts=0,actual_assembly_calls=0,numeric_factorizations=0,symbolic_analyses=0,native_calls=0,
        actual_limits=limits,periodic_stability_admitted=False,native_operator_parity_admitted=False,error=None)
    def persist():
        result['cost']=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(dest/'RESULT.json',result)
    persist()
    try:
        import numpy as np
        meta=metadata(inp);assert meta==read(cfg['metadata'])
        reader=load(inp['sources']['reader'],'candidate_frozen_reader')
        xyz,types,alpha,thole,groups=reader.read_fixture(dict(name='La',metal_type=900,xyz=meta['geometry']),meta)
        m=meta['metal_site_index'];weights=[];radii=[]
        for i,ty in enumerate(types):
            if ty!=402:continue
            # Native ANINT rounds ties away from zero. Ties are outside support.
            vec=[]
            for axis,L in enumerate(meta['cell_A']):
                d=float(xyz[i,axis]-xyz[m,axis]);n=math.copysign(math.floor(abs(d/L)+.5),d/L);vec.append(d-L*n)
            r=math.sqrt(sum(v*v for v in vec));radii.append(r);weights.append(weight(r))
        a,b,c,d=meta['response_coefficients'];h=environment(weights,c,d)
        assert all(math.isfinite(v) for v in h.values())
        h.update(metal_site_index=m,intrinsic_alpha=float(alpha[m]),native_beta=h['normalized_beta']/float(alpha[m]),
                 coefficients_a_b_c_d=[a,b,c,d],water_count=len(weights),supported_water_count=sum(r<21 for r in radii),
                 weight_sum_native_order=h['s'],compact_inner_A=20.,compact_outer_A=21.)
        np.save(dest/'water_weights.npy',np.asarray(weights,dtype=np.float64),allow_pickle=False)
        result['metal_stiffness']=h;result['water_weights']=pin(dest/'water_weights.npy')
        asm=load(inp['sources']['assembly'],'candidate_frozen_assembly')
        result['phase']='pair_query_and_assembly';result['assembly_attempts']=1;persist()
        lower,info=asm.assemble(xyz,types,alpha,thole,meta['cell_A'],cutoff=9.,ewald=.4,
            pair_overrides={tuple(x['types']):x['gamma'] for x in meta['pair_overrides']},stiffness={m:h['normalized_beta']},
            shift=0.,max_assembly_bytes=caps['assembly_bytes'],batch_pairs=4096,mutual_scales_all_one=True)
        result['actual_assembly_calls']=1;result['phase']='graph_and_array_export';persist()
        template=inp['template'];assert lower.shape==(35874,35874) and lower.nnz==16371009 and info['pairs']==1815015
        assert str(lower.data.dtype)=='float64' and str(lower.indices.dtype)==str(lower.indptr.dtype)=='int32'
        assert np.all(np.isfinite(lower.data)) and math.isfinite(info['unshifted_site_upper'])
        arrays={};hashes={}
        for k in ('data','indices','indptr'):
            array=getattr(lower,k);hashes[k]=hashlib.sha256(memoryview(np.ascontiguousarray(array)).cast('B')).hexdigest()
            if k!='data':assert hashes[k]==template[k+'_payload_sha256'],'changed graph '+k
            np.save(dest/(k+'.npy'),array,allow_pickle=False);arrays[k]=pin(dest/(k+'.npy'))
        result.update(status='completed_terminal_candidate_La_real_assembly',phase='complete',assembly=info,arrays=arrays,
                      payload_hashes=hashes,metal_coordinate_indices=list(range(3*m,3*m+3)))
        persist()
        manifest=dict(template);manifest.update(status='prepared_normalized_real_operator_for_single_static_factor',
            U_B=info['unshifted_site_upper'],arrays=arrays,metal_coordinate_indices=result['metal_coordinate_indices'],
            operator_producer=pin(__file__),operator_preparation_result=pin(dest/'RESULT.json'),
            metal_stiffness=h,source_preparation=cfg['preparation'])
        write(dest/'INPUT_MANIFEST.json',manifest)
        return 0
    except Exception as e:
        result['status']='failed_terminal_candidate_La_real_assembly'
        result['error']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc());persist();return 1

def supervise(config,release_path,out):
    cfg,_=context(config,release_path,out)
    available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
    assert available>=16*1024**3
    out.mkdir(parents=True,exist_ok=False);tick=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN)
    argv=[sys.executable,str(Path(__file__).resolve()),'worker','--config',str(config.resolve()),'--release',str(release_path.resolve()),'--out',str(out.resolve())]
    timeout=False
    with (out/'stdout.txt').open('x') as stdout,(out/'stderr.txt').open('x') as stderr:
        p=subprocess.Popen(argv,stdout=stdout,stderr=stderr,start_new_session=True)
        write(out/'START.json',dict(argv=argv,pid=p.pid,hostname=socket.gethostname(),scheduler_job_id=None,configuration=pin(config),release=pin(release_path),host_MemAvailable_bytes=available))
        try:code=p.wait(timeout=cfg['resources']['wall_seconds'])
        except subprocess.TimeoutExpired:
            timeout=True
            try:os.killpg(p.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            code=p.wait()
    u1=resource.getrusage(resource.RUSAGE_CHILDREN);rp=out/'worker/RESULT.json';mp=out/'worker/INPUT_MANIFEST.json'
    actual=json.loads(rp.read_text()) if rp.exists() else None
    ok=code==0 and actual is not None and actual['status']=='completed_terminal_candidate_La_real_assembly' and mp.exists()
    result=dict(status='completed_terminal_candidate_La_real_assembly' if ok else 'failed_terminal_candidate_La_real_assembly',
        implementation=pin(__file__),configuration=pin(config),release=pin(release_path),worker_result=pin(rp) if rp.exists() else None,
        input_manifest=pin(mp) if ok else None,exit_code=code,wall_timeout=timeout,automatic_retry=False,
        actual_assembly_calls_reported=None if actual is None else actual['actual_assembly_calls'],native_calls=0,numeric_factorizations=0,
        stdout=pin(out/'stdout.txt'),stderr=pin(out/'stderr.txt'),cost=dict(wall_seconds=time.perf_counter()-tick,
        child_user_CPU_seconds=u1.ru_utime-u0.ru_utime,child_system_CPU_seconds=u1.ru_stime-u0.ru_stime,child_max_RSS_KiB=u1.ru_maxrss,numerical_threads=1))
    write(out/'RESULT.json',result);print(json.dumps(result));return 0 if ok else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','run','worker']);p.add_argument('--inputs',type=Path)
    p.add_argument('--config',type=Path);p.add_argument('--release',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.mode=='prepare':prepare(a.inputs,a.out)
    else:raise SystemExit((supervise if a.mode=='run' else worker)(a.config,a.release,a.out))
