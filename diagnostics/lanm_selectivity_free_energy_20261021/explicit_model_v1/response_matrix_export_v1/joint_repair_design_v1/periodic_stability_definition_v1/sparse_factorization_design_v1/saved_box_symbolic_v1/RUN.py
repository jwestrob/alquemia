"""One released saved graph; supervisor is stdlib-only, worker is bounded."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import socket
import subprocess
import sys
import time
import traceback


def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def verify(p):
    assert pin(p['path'])==p,'Changed pin '+p['path']
    return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def load(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def metadata(cfg):
    handoff=read(cfg['handoff_inventory'])
    assert handoff['La_Dy_coordinates_topology_nonmetal_types_exactly_identical'] is True
    fixtures={x['metal']:x for x in handoff['fixtures']}
    normalized=[];models=[]
    for metal in ('La','Dy'):
        f=fixtures[metal];text=verify(f['geometry']).read_text().splitlines()
        assert int(text[0].split()[0])==11958 and len(text)==11960
        assert list(map(float,text[1].split()))==[49.324]*3+[90.]*3
        rows=[];counts={}
        for i,line in enumerate(text[2:]):
            t=line.split();ty=int(t[5]);assert int(t[0])==i+1
            counts[str(ty)]=counts.get(str(ty),0)+1
            rows.append([t[0],*t[2:5],'metal' if ty in (900,901) else t[5],*t[6:]])
        assert counts==f['type_counts'] and f['net_monopole_charge_e']=='0.00000'
        normalized.append(rows)
        polar={};scales={};pairs={}
        for authority in ('parameters','key'):
            seen=set()
            for line in verify(f[authority]).read_text().splitlines():
                tok=line.split('#',1)[0].split()
                if not tok:continue
                name=tok[0].lower()
                assert name not in ('polar-local-response','polar-donor-response')
                if name=='polarize' and int(tok[1]) in (402,403,405,414,900,901):
                    key=(name,int(tok[1]));assert key not in seen;seen.add(key)
                    polar[str(int(tok[1]))]=dict(alpha=float(tok[2]),thole=float(tok[3]))
                elif name in ('mutual-11-scale','mutual-12-scale','mutual-13-scale','mutual-14-scale'):
                    assert name not in seen;seen.add(name);scales[name]=float(tok[1])
                elif name=='polpair':
                    pair=tuple(sorted(map(int,tok[1:3])));key=(name,pair)
                    assert key not in seen;seen.add(key);pairs[pair]=float(tok[3])
        assert set(scales)=={'mutual-'+s+'-scale' for s in ('11','12','13','14')} and set(scales.values())=={1.}
        assert set(polar)=={'402','403','405','414','900','901'}
        assert all(v['alpha']>0 and v['thole']>=0 for v in polar.values())
        assert set(pairs)=={(402,900),(402,901)}
        models.append(dict(type_parameters=polar,mutual_scales=scales,pair_overrides=[dict(types=list(k),gamma=v) for k,v in sorted(pairs.items())]))
    assert normalized[0]==normalized[1] and models[0]==models[1]
    return dict(shared_graph=True,normalized_geometry_topology_sha256=hashlib.sha256(json.dumps(normalized[0],separators=(',',':')).encode()).hexdigest(),
                selected_fixture=fixtures['La'],other_fixture_geometry=fixtures['Dy']['geometry'],**models[0])


def context(config,release_path,out):
    cfg=json.loads(config.read_text());assert cfg['implementation']==pin(__file__)
    rel=json.loads(release_path.read_text())
    assert rel['status']=='released_one_saved_BIO_symbolic_graph'
    assert rel['configuration']==pin(config) and rel['implementation']==pin(__file__)
    assert Path(rel['output']).resolve()==out.resolve() and rel['resources']==cfg['resources']
    review=read(rel['source_review'])
    assert review['status']=='passed_saved_box_symbolic_source_review'
    assert review['implementation']==pin(__file__) and review['configuration']==pin(config)
    assert socket.gethostname()==cfg['resources']['hostname']
    assert cfg['resources']['address_space_bytes']==8*1024**3 and cfg['resources']['assembly_bytes']==4*1024**3
    assert cfg['resources']['wall_seconds']==cfg['resources']['CPU_seconds']==180
    assert cfg['cutoff_A']==9. and cfg['ewald_inverse_A']==.4 and cfg['shift']==0. and cfg['H']=='zero'
    assert cfg['batch_pairs']==4096
    env=read(cfg['environment']);assert Path(sys.prefix).resolve()==Path(env['prefix'])
    assert all(os.environ.get(k)=='1' for k in env['threads']['environment'])
    for p in env['runtime_files']:verify(p)
    sym=read(cfg['symbolic_ready']);symreview=read(cfg['symbolic_review'])
    assert symreview['status']=='passed_source_and_saved_tiny_symbolic_review' and symreview['ready']==cfg['symbolic_ready']
    assert read(sym['result'])['status']=='passed_tiny_symbolic_memory_controls'
    assembly=read(cfg['assembly_ready']);assert assembly['implementation']==cfg['assembly']
    verify(assembly['independent_source_review']);verify(cfg['reader']);verify(cfg['declaration'])
    oldreview=read(cfg['graph_reference_review']);assert oldreview['status']=='passed_real_bound_saved_result_review'
    old=read(cfg['graph_reference_result']);assert oldreview['result']==cfg['graph_reference_result']
    assert len(old['fixtures'])==2 and all(x['pairs']==1815015 for x in old['fixtures'])
    return cfg,rel


def worker(config,release_path,out):
    # Apply caps before importing numerical packages or enumerating pairs.
    preliminary=json.loads(config.read_text());caps=preliminary['resources']
    actual_limits={}
    for name,key in ((resource.RLIMIT_AS,'address_space_bytes'),(resource.RLIMIT_CPU,'CPU_seconds')):
        previous=resource.getrlimit(name);value=caps[key]
        if any(v!=resource.RLIM_INFINITY and value>v for v in previous):raise RuntimeError('would relax inherited resource limit')
        resource.setrlimit(name,(value,value));assert resource.getrlimit(name)==(value,value)
        actual_limits[key]=dict(previous=list(previous),actual=[value,value])
    cfg,release=context(config,release_path,out)
    directory=out/'worker';directory.mkdir(exist_ok=False)
    tick,cpu=time.perf_counter(),time.process_time()
    result=dict(status='started',phase='metadata',implementation=pin(__file__),configuration=pin(config),release=pin(release_path),
                actual_limits=actual_limits,symbolic_attempts=0,symbolic_analyses=0,numeric_factorizations=0,new_native_calls=0,
                numeric_stage_authorized=False,periodic_stability_admitted=False,error=None)
    def persist():
        result['cost']=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1,
                            max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        tmp=directory/'RESULT.tmp';tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(directory/'RESULT.json')
    persist()
    try:
        import numpy as np
        meta=metadata(cfg);assert meta==read(cfg['prepared_metadata'])
        f=meta['selected_fixture'];reader=load(cfg['reader'],'frozen_box_reader')
        xyz,types,alpha,thole,groups=reader.read_fixture(dict(name='La',metal_type=900,xyz=f['geometry']),meta)
        assert len(xyz)==11958 and np.all(np.isfinite(alpha)) and np.all(alpha>0)
        assembler=load(cfg['assembly'],'frozen_sparse_assembly')
        symbolic=load(read(cfg['symbolic_ready'])['implementation'],'frozen_symbolic')
        result['phase']='pair_query_and_assembly';persist()
        lower,info=assembler.assemble(xyz,types,alpha,thole,[49.324]*3,cutoff=9.,ewald=.4,
                      pair_overrides={tuple(x['types']):x['gamma'] for x in meta['pair_overrides']},
                      stiffness={},shift=0.,max_assembly_bytes=caps['assembly_bytes'],batch_pairs=4096,mutual_scales_all_one=True)
        assert info['pairs']==1815015 and lower.shape==(35874,35874) and lower.nnz==16371009
        def array_hash(a):return hashlib.sha256(memoryview(np.ascontiguousarray(a)).cast('B')).hexdigest()
        result['assembly']=dict(info=info,data_sha256=array_hash(lower.data),indices_sha256=array_hash(lower.indices),
                                indptr_sha256=array_hash(lower.indptr),index_dtype=str(lower.indices.dtype),
                                allocated_CSC_bytes=int(lower.data.nbytes+lower.indices.nbytes+lower.indptr.nbytes))
        result['phase']='symbolic_analysis';persist()
        headroom=max(0,caps['address_space_bytes']-symbolic.virtual_bytes())
        result['symbolic_attempts']=1;persist()
        factor,report=symbolic.analyze_only(lower,max_numeric_estimate_bytes=min(caps['assembly_bytes'],headroom),reserve_bytes=256*1024**2)
        result['symbolic_analyses']=1
        limit=int(np.iinfo(np.dtype(report['index_dtype'])).max)
        report['index_capacity_max']=limit;report['future_numeric_index_capacity_allowed']=report['symbolic_L_nnz']<=limit
        report['future_numeric_forecast_and_index_allowed']=report['numeric_memory_allowed'] and report['future_numeric_index_capacity_allowed']
        if factor is not None:assert not factor.is_numeric
        result['symbolic']=report;result['metadata']=cfg['prepared_metadata'];result['shared_geometry_graph_covers']=['La','Dy']
        result['status']='completed_single_saved_graph_symbolic_check';result['phase']='finished_no_numeric_stage'
        del factor,lower
    except Exception as e:
        result['status']='failed_single_saved_graph_symbolic_check'
        result['error']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc())
    persist();return 0 if result['status']=='completed_single_saved_graph_symbolic_check' else 1


def supervise(config,release_path,out):
    cfg,_=context(config,release_path,out)
    available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
    assert available>=16*1024**3,'insufficient declared host headroom'
    out.mkdir(parents=True,exist_ok=False)
    argv=[sys.executable,str(Path(__file__).resolve()),'worker','--config',str(config.resolve()),'--release',str(release_path.resolve()),'--out',str(out.resolve())]
    tick=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN);timed_out=False
    with (out/'stdout.txt').open('x') as stdout,(out/'stderr.txt').open('x') as stderr:
        proc=subprocess.Popen(argv,stdout=stdout,stderr=stderr,start_new_session=True)
        (out/'START.json').write_text(json.dumps(dict(argv=argv,worker_pid=proc.pid,hostname=socket.gethostname(),scheduler_job_id=None,
                                                      configuration=pin(config),release=pin(release_path),host_MemAvailable_bytes=available),indent=2)+'\n')
        try:code=proc.wait(timeout=cfg['resources']['wall_seconds'])
        except subprocess.TimeoutExpired:
            timed_out=True
            try:os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            code=proc.wait()
    u1=resource.getrusage(resource.RUSAGE_CHILDREN);rp=out/'worker/RESULT.json'
    actual=json.loads(rp.read_text()) if rp.exists() else None
    complete=code==0 and actual is not None and actual['status']=='completed_single_saved_graph_symbolic_check'
    result=dict(status='completed_single_saved_graph_symbolic_check' if complete else 'failed_single_saved_graph_symbolic_check',
                implementation=pin(__file__),configuration=pin(config),release=pin(release_path),worker_result=pin(rp) if rp.exists() else None,
                exit_code=code,wall_timeout=timed_out,stdout=pin(out/'stdout.txt'),stderr=pin(out/'stderr.txt'),
                worker_status=None if actual is None else actual['status'],symbolic_analyses_reported=None if actual is None else actual['symbolic_analyses'],
                numeric_factorizations=0,new_native_calls=0,automatic_retry=False,
                cost=dict(outer_wall_seconds=time.perf_counter()-tick,child_user_CPU_seconds=u1.ru_utime-u0.ru_utime,
                          child_system_CPU_seconds=u1.ru_stime-u0.ru_stime,child_max_RSS_KiB=u1.ru_maxrss,numerical_threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return 0 if complete else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    raise SystemExit((supervise if a.mode=='run' else worker)(a.config,a.release,a.out))
