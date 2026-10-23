"""One existing Hadamard4 fixture; preserve and reuse the completed tiny panel."""
import argparse
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import numpy as np
from scipy.sparse import csc_array


def pin(p):
    p=Path(p).resolve()
    return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def read(p):
    assert pin(p['path'])==p
    return json.loads(Path(p['path']).read_text())


def load(p):
    assert pin(p['path'])==p
    spec=importlib.util.spec_from_file_location('recovered_backend',p['path'])
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def json_scalar(value):
    if isinstance(value,np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def run(config,out):
    start,cpu=time.monotonic(),time.process_time()
    assert json.loads(json.dumps({'bool':np.bool_(True),'float':np.float64(1.)},default=json_scalar))=={'bool':True,'float':1.}
    cfg=json.loads(Path(config).read_text());assert pin(__file__)==cfg['implementation']
    dec=read(cfg['declaration']);old=read(cfg['preserved_result']);diag=read(cfg['diagnosis'])
    env=read(cfg['environment']);backend=load(cfg['backend'])
    assert Path(sys.prefix).resolve()==Path(env['prefix'])
    for p in env['runtime_files']:assert pin(p['path'])==p
    for p in [diag['source'],diag['source_archive'],diag['package_recipe'],diag['package_patch']]:assert pin(p['path'])==p
    assert diag['source_archive']['sha256']=='2f98b2c7d089a65d3eaa7533393966b9e695a429509cf0c4f5d9bca05a3f8687'
    assert all(os.environ.get(k)=='1' for k in env['threads']['environment'])
    blas=ctypes.CDLL(str(Path(sys.prefix)/'lib/libopenblasp-r0.3.34.so'));omp=ctypes.CDLL(str(Path(sys.prefix)/'lib/libgomp.so.1'))
    blas.openblas_get_num_threads.restype=omp.omp_get_max_threads.restype=ctypes.c_int
    assert blas.openblas_get_num_threads()==omp.omp_get_max_threads()==1
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    checks=[]
    result=dict(status='started',implementation=pin(__file__),backend=cfg['backend'],configuration=pin(config),declaration=cfg['declaration'],
                diagnosis=cfg['diagnosis'],preserved_result=cfg['preserved_result'],preserved_process=cfg['preserved_process'],
                original_result_status=old['status'],settings=backend.SETTINGS,checks=checks,
                symbolic_analyses=0,numeric_factorizations=0,solves=0,native_calls=0,saved_molecular_matrix_reads=0,
                error=None,threads=dict(openblas=blas.openblas_get_num_threads(),openmp=omp.omp_get_max_threads()))
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    def phase(name):
        result[name]+=1
        assert result[name]<=1
    def solve_count():
        result['solves']+=1
        assert result['solves']<=12
    def persist():
        result['cost']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
                            peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1)
        tmp=out/'RESULT.tmp';tmp.write_text(json.dumps(result,indent=2,allow_nan=False,default=json_scalar)+'\n');tmp.replace(out/'RESULT.json')
    persist()
    try:
        assert old['status']=='failed_manufactured_supernodal_mode_qualification' and old['error'] is None
        assert old['symbolic_analyses']==old['numeric_factorizations']==11 and old['solves']==62
        assert len(old['manufactured'])==11 and len(old['input_errors'])==5
        allowed=set();reused=[]
        for row in old['manufactured']:
            name=row['name'];ar=row.get('sparse_audit')
            if ar is not None:
                allowed.add(name+'/copy_flags_permutation_and_solves')
                b,a=ar['before_copy'],ar['after_copy']
                check(name+'/saved_expected_storage_copy_invariance',b==a and b['is_numeric'] and b['is_ll'] and b['is_lower'] and not b['is_super'])
                check(name+'/saved_arithmetic',ar['finite_positive_diagonal'] and all(x['passed'] and x['finite'] and
                      x['factor_action_error']<=ar['arithmetic_bound'] and x['solve_backward_error']<=ar['arithmetic_bound'] for x in ar['probes']))
                check(name+'/saved_dense_correspondence',row['tiny_dense_comparison']['whole_factor_backward_error']<=ar['arithmetic_bound'] and
                      all(x['K_action_difference']<=ar['arithmetic_bound'] and x['L_action_difference']<=ar['arithmetic_bound'] for x in row['tiny_dense_comparison']['action_correspondence']))
                if row['dense_classification']=='positive':allowed.add(name+'/positive_mode_result')
            else:
                check(name+'/saved_nonpositive_exception',row['exception']['type']=='CholmodNotPositiveDefiniteError' and not row['sufficient_test_admitted'])
            if row['dense_classification']=='numerically_unresolved':check(name+'/saved_boundary_unresolved',not row['sufficient_test_admitted'])
            reused.append(dict(name=name,dense_classification=row['dense_classification'],arithmetic_reused=ar is not None,
                               physical_full_operator_admission=False))
        failed={x['name'] for x in old['checks'] if not x['passed']}
        check('only_wrong_postnumeric_supernodal_guard_failed',failed==allowed==set(diag['failed_checks']))
        check('saved_five_malformed_inputs_rejected',all(x['exception'] is not None and x['numeric_attempts']==0 for x in old['input_errors']))
        corruption=old['manufactured'][0]['copy_corruption_control']
        check('saved_deliberate_copy_corruption_rejected',corruption['rejected'] and corruption['original_factor_unchanged'])
        result['reused_evidence']=reused
        assert all(x['passed'] for x in checks)
        prior=next(x for x in old['manufactured'] if x['name']==dec['fixture'])
        assert prior['name']=='hadamard4_spd' and prior['dense_classification']=='positive'
        k=np.array(prior['K'],dtype=float);assert k.shape==(4,4)
        lower=csc_array(np.tril(k));lower.setdiag(np.diag(k));lower.sort_indices()
        f,error,warnings,symbolic=backend.factor_once(lower,on_phase=phase)
        result['one_fixture']=dict(name=prior['name'],A=prior['A'],U=prior['U'],K=prior['K'],pre_numeric=symbolic,
                                  post_numeric=backend.metadata(f),exception=error,warnings=warnings)
        check('pre_numeric_supernodal_symbolic',symbolic['is_super'] and not symbolic['is_numeric'])
        check('one_numeric_call_without_exception',error is None and not warnings)
        assert error is None
        ar,l,p=backend.audit_factor(f,lower,backend.ten_probes(4),on_solve=solve_count)
        result['one_fixture']['sparse_audit']=ar
        check('expected_converted_simplicial_LL_and_sparse_arithmetic',ar['passed'])
        ld=l.toarray();dense_error=float(np.linalg.norm(k[np.ix_(p,p)]-ld@ld.T,np.inf))/max(1.,float(np.linalg.norm(k,np.inf)))
        result['one_fixture']['whole_factor_backward_error']=dense_error
        check('tiny_dense_factor_correspondence',dense_error<=ar['arithmetic_bound'])
        b=np.arange(1.,5.);solve_count();before=np.asarray(f.solve(b,system='A'))
        l.data[l.indptr[0]]*=1.1
        solve_count();after=np.asarray(f.solve(b,system='A'))
        independent=np.array_equal(before,after) and np.array_equal(f.get_factor(kind='LL',lower=True).toarray(),ld)
        result['one_fixture']['copied_L_independent']=bool(independent)
        check('nontrivial_copy_and_original_solve_independence',independent)
        check('one_fixture_cap',result['symbolic_analyses']==result['numeric_factorizations']==1 and result['solves']==12)
        result['combined_counts']=dict(original_matrices_reused=11,new_existing_fixture_rechecks=1,total_symbolic_analyses=12,total_numeric_factorizations=12,total_solves=74,new_native_calls=0)
        result['status']='passed_manufactured_supernodal_numeric_path_with_simplicial_storage' if all(x['passed'] for x in checks) else 'failed_supernodal_representation_recovery'
    except Exception as exc:
        result['status']='failed_supernodal_representation_recovery'
        result['error']=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    finally:persist()
    print(json.dumps(dict(status=result['status'],result=pin(out/'RESULT.json'),cost=result['cost'])))
    return 0 if result['status']=='passed_manufactured_supernodal_numeric_path_with_simplicial_storage' else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();raise SystemExit(run(a.config,a.out))
