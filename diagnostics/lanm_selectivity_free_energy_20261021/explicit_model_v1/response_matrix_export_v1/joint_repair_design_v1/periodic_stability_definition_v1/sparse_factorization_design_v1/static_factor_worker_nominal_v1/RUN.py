"""One nominal saved-operator factor attempt; no assembly or native engine."""
import argparse
import ctypes
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
import warnings

N = 35874
NNZ = 16371009
INDEX_SHA = 'f52075f7ef0db207af9a6f163e51f4b2303898ae5a99fb48a2c7bbacc6a526a1'
INDPTR_SHA = '095e3bacb7be69cd2a953545a84fc44250d3f708bfa3756fa60f1901021d7e25'


def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return dict(path=str(path),sha256=h.hexdigest())


def verify(p):
    if not isinstance(p,dict) or set(p)!= {'path','sha256'} or pin(p['path'])!=p:
        raise ValueError('missing or changed required pin')
    return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def load(p):
    spec=importlib.util.spec_from_file_location('qualified_static_backend',verify(p))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


NOMINAL_SHA = '2aa3ac434f3667302fb5cf70d3adda1444e942c6f839e052efef82305ec2e10d'


def nominal_parameters(inputs, review):
    # e_app is unknown, not zero. This worker implements only the nominal contract.
    for item in (inputs, review):
        assert item['e_app'] is None and item['application_error_bound_established'] is False
    for key in ('U_B','U_R','delta_odd'):
        assert inputs[key] == review[key]
        assert type(inputs[key]) in (int,float) and math.isfinite(inputs[key])
    ub,ur,odd=(float(inputs[key]) for key in ('U_B','U_R','delta_odd'))
    upper=ub+ur
    assert ur>=0 and odd>=0 and math.isfinite(upper) and upper>0
    tau=1e-6*upper+odd
    assert math.isfinite(tau) and tau>0
    return dict(U_B=ub,U_R=ur,delta_odd=odd,U_nominal=upper,tau=tau,
                e_app=None,application_error_bound_established=False,
                scope='nominal_saved_operator_only')


def nominal_provenance(cfg, inputs, review):
    assert cfg['nominal_contract']['sha256']==NOMINAL_SHA
    contract=read(cfg['nominal_contract'])['recommended_contract']
    assert contract['tau_nominal']=='1e-6*(U_B+U_R)+delta_odd'
    assert contract['e_app'] is None and contract['application_error_bound_established'] is False
    for key in ('nominal_contract','pme_result','pme_root_review','pme_bound_review'):
        assert inputs[key]==review[key]==cfg[key]
    assert inputs['metal']==review['metal']==cfg['metal']=='La'
    root=read(cfg['pme_root_review']);actual=read(cfg['pme_result']);independent=read(cfg['pme_bound_review'])
    assert root['status']=='admitted_full_box_operator_initialization_nominal_data'
    assert root['execution']==independent['result']==cfg['pme_result']
    assert root['independent_actual_review']==cfg['pme_bound_review']
    assert actual['status']=='full_box_operator_initialization_complete' and actual['error'] is None
    assert root['composition']==actual['composition']==cfg['composition_manifest']
    assert root['e_app'] is None and root['application_error_bound_established'] is False
    for key,source in (('result','campaign_result'),('configuration','campaign_configuration'),
                       ('selected','candidate'),('actual_review','terminal_fit_review')):
        assert actual['candidate'][key]==cfg[source]
    row=actual['cells']['La'];rr=root['cells']['La'];ir=independent['cells']['La']
    assert row['status']=='complete' and row['site_checks']['passed'] is True
    assert row['export']==rr['export']==ir['export']==cfg['pme_export']
    assert row['inputs']==rr['inputs'] and row['inputs']['system.xyz']['sha256']==inputs['geometry']['sha256']
    norms=row['operator_norms']
    assert norms==rr['operator_norms']==ir['operator_norms']
    assert norms['arithmetic']=='ordinary_double_not_outward'
    assert norms['q_negative_count']==0 and norms['q_min']>=0
    for key in ('q_min','q_max','q_abs_max','q_conjugate_difference_max','U_R_nominal','evenness_perturbation_nominal'):
        assert math.isfinite(norms[key]) and norms[key]>=0
    assert inputs['U_R']==norms['U_R_nominal'] and inputs['delta_odd']==norms['evenness_perturbation_nominal']
    for key in ('operator_producer','operator_preparation_result'):
        verify(inputs[key]);assert review[key]==inputs[key]
    nominal_parameters(inputs,review)


def context(config,release_path,out):
    cfg=json.loads(config.read_text())
    if cfg['status']!='bound_single_static_factor_configuration':
        raise ValueError('static configuration remains unbound')
    assert cfg['implementation']==pin(__file__)
    for key in ('campaign_result','campaign_configuration','candidate','terminal_fit_review',
                'composition_manifest','pme_export','pme_bound_review','input_manifest','input_review'):
        verify(cfg[key])
    campaign=read(cfg['campaign_result'])
    assert campaign['status']=='bounded_component_fit_complete' and campaign['error'] is None
    assert campaign['mode']=='run' and campaign['candidate_frozen'] is True
    assert campaign['configuration']==cfg['campaign_configuration'] and campaign['selected']==cfg['candidate']
    assert campaign['vectors'][campaign['selected_key']]['status']=='complete'
    assert campaign['vectors'][campaign['selected_key']]['artifact']==cfg['candidate']
    for key,status_key in [('terminal_fit_review','terminal_fit_review_pass_status'),('pme_bound_review','pme_review_pass_status')]:
        expected=cfg[status_key]
        assert isinstance(expected,str) and expected.startswith('passed_')
        assert read(cfg[key])['status']==expected
    inputs=read(cfg['input_manifest']);review=read(cfg['input_review'])
    assert inputs['status']=='prepared_normalized_real_operator_for_single_static_factor'
    assert review['status']=='accepted_candidate_PME_and_real_operator_for_one_static_factor'
    assert review['input_manifest']==cfg['input_manifest']
    for key in ('campaign_result','campaign_configuration','candidate','terminal_fit_review',
                'composition_manifest','pme_export','pme_bound_review'):
        assert review[key]==cfg[key]
    for key in ('candidate','composition_manifest','pme_export'):
        assert inputs[key]==cfg[key]
    nominal_provenance(cfg,inputs,review)
    assert all(review[k] is True for k in ('terminal_candidate_reviewed','initialized_PME_nominal_data_and_bound_reviewed',
               'real_space_operator_and_upper_bound_reviewed','symmetry_checks_passed','exact_active_map_and_graph_reviewed'))
    assert inputs['matrix']=='normalized_B0_plus_H_unshifted_lower' and inputs['shape']==[N,N] and inputs['nnz']==NNZ
    verify(inputs['geometry'])
    assert review['geometry']==inputs['geometry'] and inputs['cell_A']==[49.324]*3
    assert inputs['cutoff_A']==9. and inputs['ewald_inverse_A']==.4 and inputs['shift']==0.
    assert inputs['all_alpha_positive'] is True and inputs['mutual_scales_all_one'] is True
    assert inputs['U_B']==review['U_B'] and inputs['U_R']==review['U_R']
    assert math.isfinite(inputs['U_B']) and math.isfinite(inputs['U_R']) and inputs['U_R']>=0
    assert math.isfinite(inputs['U_B']+inputs['U_R']) and inputs['U_B']+inputs['U_R']>0
    assert math.isfinite(review['normalized_symmetry_max_absolute']) and 0<=review['normalized_symmetry_max_absolute']<=1e-10
    assert inputs['indices_payload_sha256']==INDEX_SHA and inputs['indptr_payload_sha256']==INDPTR_SHA
    graph=read(cfg['graph_authority'])
    assert graph['status']=='completed_single_saved_graph_symbolic_check'
    assert graph['assembly']['indices_sha256']==INDEX_SHA and graph['assembly']['indptr_sha256']==INDPTR_SHA
    assert graph['symbolic']['symbolic_L_nnz']==399831681
    mode=read(cfg['mode_qualification'])
    assert mode['status']=='passed_tiny_supernodal_numeric_path_with_explicit_converted_LL_storage'
    assert mode['backend']==cfg['backend'] and read(mode['result'])['status']=='passed_manufactured_supernodal_numeric_path_with_simplicial_storage'
    env=read(cfg['environment']);assert Path(sys.prefix).resolve()==Path(env['prefix'])
    for p in env['runtime_files']:verify(p)
    caps=cfg['resources']
    assert caps==dict(hostname='biotite',address_space_bytes=32*1024**3,wall_seconds=900,CPU_seconds=10800,
                     openblas_threads=12,openmp_threads=1,host_min_available_bytes=64*1024**3)
    assert socket.gethostname()==caps['hostname'] and not os.environ.get('SLURM_JOB_ID')
    for k in env['threads']['environment']:
        assert os.environ.get(k)==('12' if k=='OPENBLAS_NUM_THREADS' else '1')
    release=json.loads(release_path.read_text())
    assert release['status']=='released_one_bound_static_factor_attempt'
    assert release['implementation']==pin(__file__) and release['configuration']==pin(config)
    assert release['input_review']==cfg['input_review'] and release['resources']==caps
    assert Path(release['output']).resolve()==out.resolve()
    source_review=read(release['source_review'])
    assert source_review['status']=='passed_nominal_static_factor_worker_source_review'
    assert source_review['implementation']==pin(__file__) and source_review['configuration']==pin(config)
    return cfg,inputs


def worker(config,release_path,out):
    # Limits are installed before numerical imports or any matrix-array read.
    preliminary=json.loads(config.read_text());caps=preliminary['resources'];limits={}
    for key,which in [('address_space_bytes',resource.RLIMIT_AS),('CPU_seconds',resource.RLIMIT_CPU)]:
        old=resource.getrlimit(which);new=caps[key]
        if any(v!=resource.RLIM_INFINITY and new>v for v in old):raise RuntimeError('would relax inherited limit')
        resource.setrlimit(which,(new,new));limits[key]=dict(previous=list(old),actual=list(resource.getrlimit(which)))
    cfg,inputs=context(config,release_path,out)
    directory=out/'worker';directory.mkdir(exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time()
    result=dict(status='started',phase='imports',implementation=pin(__file__),configuration=pin(config),release=pin(release_path),
                input_manifest=cfg['input_manifest'],input_review=cfg['input_review'],limits=limits,
                nominal_contract=cfg['nominal_contract'],nominal_only=True,e_app=None,application_error_bound_established=False,
                symbolic_attempts=0,numeric_attempts=0,solve_attempts=0,new_native_calls=0,
                sufficient_numerical_test_passed=False,physical_accuracy_admitted=False,MD_admitted=False,error=None)
    def scalar(x):
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    def persist():
        result['cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,
                            peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,openblas_threads=12,openmp_threads=1)
        tmp=directory/'RESULT.tmp';tmp.write_text(json.dumps(result,indent=2,allow_nan=False,default=scalar)+'\n');tmp.replace(directory/'RESULT.json')
    def phase(name):
        key={'symbolic_analyses':'symbolic_attempts','numeric_factorizations':'numeric_attempts'}[name]
        result[key]+=1;assert result[key]<=1;result['phase']=name;persist()
    def solve_count():
        result['solve_attempts']+=1;assert result['solve_attempts']<=10
    try:
        import numpy as np
        from scipy.sparse import csc_array
        persist()
        assert json.loads(json.dumps({'numpy_bool':np.bool_(True)},default=scalar))=={'numpy_bool':True}
        blas=ctypes.CDLL(str(Path(sys.prefix)/'lib/libopenblasp-r0.3.34.so'));omp=ctypes.CDLL(str(Path(sys.prefix)/'lib/libgomp.so.1'))
        blas.openblas_get_num_threads.restype=omp.omp_get_max_threads.restype=ctypes.c_int
        result['actual_threads']=dict(openblas=blas.openblas_get_num_threads(),openmp=omp.omp_get_max_threads())
        assert result['actual_threads']==dict(openblas=12,openmp=1)
        backend=load(cfg['backend']);result['phase']='bound_array_read';persist()
        arrays={k:np.load(verify(inputs['arrays'][k]),mmap_mode='r',allow_pickle=False) for k in ('data','indices','indptr')}
        data,idx,ptr=(arrays[k] for k in ('data','indices','indptr'))
        assert data.dtype==np.float64 and data.shape==(NNZ,) and idx.dtype==np.int32 and idx.shape==(NNZ,)
        assert ptr.dtype==np.int32 and ptr.shape==(N+1,)
        assert hashlib.sha256(memoryview(idx).cast('B')).hexdigest()==INDEX_SHA
        assert hashlib.sha256(memoryview(ptr).cast('B')).hexdigest()==INDPTR_SHA
        k=csc_array((np.array(data,copy=True),idx,ptr),shape=(N,N),copy=False);backend.validate_lower(k)
        nominal=nominal_parameters(inputs,read(cfg['input_review']));tau=nominal['tau']
        # Canonical lower columns start at their explicit diagonal; pattern is unchanged.
        k.data[k.indptr[:-1]]-=tau;backend.validate_lower(k)
        result['operator']=dict(shape=[N,N],nnz=NNZ,**nominal,
                                shifted_data_sha256=hashlib.sha256(memoryview(k.data).cast('B')).hexdigest(),
                                interpretation='Nominal saved-operator sufficient numerical test only; no floating application bound, trajectory or energy/force admission. Failure does not establish full-PME instability.')
        f,error,ww,pre=backend.factor_once(k,on_phase=phase)
        result['pre_numeric']=pre;result['post_numeric']=backend.metadata(f);result['backend_exception']=error;result['backend_warnings']=ww;persist()
        if error is not None or ww:
            if error is not None and error['type']!='CholmodNotPositiveDefiniteError':
                result['status']='failed_static_numeric_attempt'
            else:
                result['status']='inconclusive_real_split_test' if error is not None else 'inconclusive_numeric_warning'
            result['phase']='discarded_factor_without_solve'
        else:
            coords=inputs['metal_coordinate_indices']
            assert len(coords)==3 and len(set(coords))==3 and all(isinstance(x,int) and 0<=x<N for x in coords)
            result['phase']='sparse_copy_and_sampled_residuals';persist()
            with warnings.catch_warnings(record=True) as observed:
                warnings.simplefilter('always')
                ar,lower,p=backend.audit_factor(f,k,backend.ten_probes(N,tuple(coords)),on_solve=solve_count)
            result['audit']=ar;result['audit_warnings']=[dict(type=type(x.message).__name__,message=str(x.message)) for x in observed]
            passed=ar['passed'] and not observed
            result['sufficient_numerical_test_passed']=bool(passed)
            result['status']='passed_nominal_static_sufficient_numerical_test' if passed else 'inconclusive_static_arithmetic_or_boundary'
            result['phase']='finished'
        del f
    except Exception as exc:
        result['status']='failed_static_numeric_attempt'
        result['error']=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    persist()
    return 0 if result['sufficient_numerical_test_passed'] else 2


def supervise(config,release_path,out):
    cfg,_=context(config,release_path,out)
    available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
    assert available>=cfg['resources']['host_min_available_bytes']
    out.mkdir(parents=True,exist_ok=False)
    argv=[sys.executable,str(Path(__file__).resolve()),'worker','--config',str(config.resolve()),'--release',str(release_path.resolve()),'--out',str(out.resolve())]
    tick=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);timeout=False
    with (out/'stdout.txt').open('x') as so,(out/'stderr.txt').open('x') as se:
        p=subprocess.Popen(argv,stdout=so,stderr=se,start_new_session=True)
        (out/'START.json').write_text(json.dumps(dict(argv=argv,worker_pid=p.pid,hostname=socket.gethostname(),job_id=None,
              configuration=pin(config),release=pin(release_path),host_available_bytes=available),indent=2)+'\n')
        try:code=p.wait(timeout=cfg['resources']['wall_seconds'])
        except subprocess.TimeoutExpired:
            timeout=True
            try:os.killpg(p.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            code=p.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN);rp=out/'worker/RESULT.json'
    actual=json.loads(rp.read_text()) if rp.exists() else None
    result=dict(status=actual['status'] if actual is not None and code in (0,2) else 'failed_or_timed_out_static_attempt',
                implementation=pin(__file__),configuration=pin(config),release=pin(release_path),worker_result=pin(rp) if rp.exists() else None,
                exit_code=code,wall_timeout=timeout,stdout=pin(out/'stdout.txt'),stderr=pin(out/'stderr.txt'),automatic_retry=False,new_native_calls=0,
                cost=dict(wall_seconds=time.monotonic()-tick,child_user_CPU_seconds=after.ru_utime-before.ru_utime,
                          child_system_CPU_seconds=after.ru_stime-before.ru_stime,child_peak_RSS_KiB=after.ru_maxrss))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return code if code>=0 else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    raise SystemExit((supervise if a.mode=='run' else worker)(a.config,a.release,a.out))
