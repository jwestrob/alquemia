"""The declared eleven manufactured supernodal controls; no molecular inputs."""
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
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def read_pin(p):
    assert pin(p['path']) == p
    return json.loads(Path(p['path']).read_text())


def load(p, name):
    assert pin(p['path']) == p
    spec = importlib.util.spec_from_file_location(name, p['path'])
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def run(config, out):
    start, cpu = time.monotonic(), time.process_time()
    cfg = json.loads(Path(config).read_text())
    assert pin(__file__) == cfg['implementation']
    dec = read_pin(cfg['declaration'])
    env = read_pin(cfg['environment'])
    old_result = read_pin(cfg['original_result'])
    old = load(cfg['original_source'], 'original_manufactured')
    backend = load(cfg['backend'], 'fixed_supernodal_backend')
    assert dec['maximum_dimension'] == 6 and dec['numeric_factorization_cap'] == 11
    assert env['status'] == 'passed_exact_isolated_environment_audit'
    assert Path(sys.prefix).resolve() == Path(env['prefix'])
    for p in env['runtime_files']:
        assert pin(p['path']) == p
    assert all(os.environ.get(k) == '1' for k in env['threads']['environment'])
    blas = ctypes.CDLL(str(Path(sys.prefix)/'lib/libopenblasp-r0.3.34.so'))
    omp = ctypes.CDLL(str(Path(sys.prefix)/'lib/libgomp.so.1'))
    blas.openblas_get_num_threads.restype = omp.omp_get_max_threads.restype = ctypes.c_int
    assert blas.openblas_get_num_threads() == omp.omp_get_max_threads() == 1
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    rows, checks = [], []
    result = dict(status='started', implementation=pin(__file__), backend=cfg['backend'],
                  configuration=pin(config), declaration=cfg['declaration'], environment=cfg['environment'],
                  original_source=cfg['original_source'], original_result=cfg['original_result'],
                  settings=backend.SETTINGS, manufactured=rows, input_errors=[], checks=checks,
                  symbolic_analyses=0, numeric_factorizations=0, solves=0, rank_updates=0,
                  molecular_matrix_reads=0, native_calls=0, new_installs=0, error=None,
                  threads=dict(openblas=blas.openblas_get_num_threads(), openmp=omp.omp_get_max_threads()))
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
    def count_phase(name):
        result[name] += 1
        assert result[name] <= dec['numeric_factorization_cap']
    def count_solve():
        result['solves'] += 1
    def persist():
        result['cost'] = dict(wall_seconds=time.monotonic()-start, CPU_seconds=time.process_time()-cpu,
                              peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1)
        tmp = out/'RESULT.tmp'
        tmp.write_text(json.dumps(result, indent=2, allow_nan=False,
                                  default=lambda x: x.item() if isinstance(x,np.generic) else str(x))+'\n')
        tmp.replace(out/'RESULT.json')
    persist()
    try:
        cases = old.manufactured()
        assert len(cases) == len(old_result['manufactured']) == 11
        for (name, a, u), prior in zip(cases, old_result['manufactured'], strict=True):
            assert name == prior['name'] and a.tolist() == prior['A'] and u == prior['U']
            assert len(a) <= 6
            a = old.guard(a); k = a-old.THRESHOLD*u*np.eye(len(a))
            dr, kr = old.dense(a), old.dense(k)
            margin = kr['eigenvalues'][0]
            resolution = 100*len(a)*old.EPS*max(1.,old.norm(a))
            expected = 'numerically_unresolved' if abs(margin)<=resolution else ('positive' if margin>0 else 'not_positive')
            lower = csc_array(np.tril(k))
            # The exact-zero shifted diagonal must remain explicitly represented.
            lower.setdiag(np.diag(k)); lower.sort_indices()
            f, error, warnings = backend.factor_once(lower, on_phase=count_phase)
            row = dict(name=name, A=a.tolist(), U=u, K=k.tolist(), dense=dr, shifted_dense=kr,
                       roundoff_resolution=resolution, dense_classification=expected,
                       exception=error, warnings=warnings, backend_metadata=backend.metadata(f),
                       sufficient_test_admitted=False, failed_factor_reused=False)
            rows.append(row)
            check(name+'/exact_original_fixture',a.tolist()==prior['A'] and u==prior['U'])
            check(name+'/dense_arithmetic',dr['passed'] and kr['passed'])
            check(name+'/supplied_upper_bound',dr['eigenvalues'][-1] <= u)
            if error is None:
                probes = backend.ten_probes(len(a))
                ar, l, p = backend.audit_factor(f, lower, probes, on_solve=count_solve)
                row['sparse_audit'] = ar
                ld = l.toarray()  # Only declared n<=6 manufactured matrices.
                dense_factor = old.norm(k[np.ix_(p,p)]-ld@ld.T)/max(1.,old.norm(k))
                correspond = []
                for pname, v in probes:
                    sparse_k = backend.action(lower,v)
                    sparse_l = np.asarray(l@(l.T@v[p]))
                    dense_l = ld@(ld.T@v[p])
                    scale=max(1.,old.norm(k))*old.norm(v)
                    correspond.append(dict(name=pname,
                        K_action_difference=old.norm(sparse_k-k@v)/scale,
                        L_action_difference=old.norm(sparse_l-dense_l)/scale))
                row['tiny_dense_comparison']=dict(whole_factor_backward_error=dense_factor,
                    norm_difference=abs(ar['K_infinity_norm']-old.norm(k)),action_correspondence=correspond)
                bound=ar['arithmetic_bound']
                check(name+'/sparse_dense_correspondence',dense_factor<=bound and
                      abs(ar['K_infinity_norm']-old.norm(k))<=bound*max(1.,old.norm(k)) and
                      all(z['K_action_difference']<=bound and z['L_action_difference']<=bound for z in correspond))
                check(name+'/copy_flags_permutation_and_solves',ar['passed'])
                if name=='diagonal6_spd':
                    before=l.toarray().copy(); b=probes[0][1]
                    count_solve(); xb=np.asarray(f.solve(b,system='A'))
                    col0=l.indptr[0]; assert l.indices[col0]==0
                    l.data[col0]*=1.1
                    errors=[backend.factor_action_error(lower,l,p,v,ar['K_infinity_norm'])[0] for _,v in probes]
                    count_solve(); xa=np.asarray(f.solve(b,system='A'))
                    fresh=f.get_factor(kind='LL',lower=True).toarray()
                    independent=np.array_equal(xb,xa) and np.array_equal(before,fresh) and not np.array_equal(before,l.toarray())
                    row['copy_corruption_control']=dict(relative_diagonal_change=.1,
                        maximum_action_error=max(errors),arithmetic_bound=bound,
                        rejected=bool(max(errors)>bound),original_factor_unchanged=bool(independent))
                    check('copied_diagonal_corruption_rejected',max(errors)>bound)
                    check('copied_L_independent_of_original_factor',independent)
                if expected=='positive':
                    row['sufficient_test_admitted']=bool(ar['passed'] and not warnings)
            if expected=='positive':
                check(name+'/positive_mode_result',error is None and row['sufficient_test_admitted'])
            elif expected=='not_positive':
                check(name+'/nonpositive_error_not_promoted',error is not None and
                      error['type']=='CholmodNotPositiveDefiniteError' and not row['sufficient_test_admitted'])
                row['interpretation']='Sufficient real-part test failed; no full-system instability inference.'
            else:
                check(name+'/boundary_unresolved',not row['sufficient_test_admitted'])
                row['interpretation']='Numerically unresolved boundary regardless of return status.'
            del f
            persist()
        bad=[('nonsquare',np.ones((2,3))),('asymmetric',np.array([[1,.01],[0,1]])),
             ('NaN',np.array([[np.nan,0],[0,1]])),('infinity',np.diag([np.inf,1])),('complex',np.eye(2,dtype=complex))]
        for name,a in bad:
            try:old.guard(a);error=None
            except ValueError as exc:error=old.exception(exc)
            result['input_errors'].append(dict(name=name,exception=error,numeric_attempts=0))
            check(name+'/original_guard_rejection',error is not None)
        check('exact_finite_inventory',len(rows)==11 and len(result['input_errors'])==5 and
              result['numeric_factorizations']==result['symbolic_analyses']==11 and result['rank_updates']==0)
        result['status']='passed_manufactured_supernodal_mode_qualification' if all(r['passed'] for r in checks) else 'failed_manufactured_supernodal_mode_qualification'
    except Exception as exc:
        result['status']='failed_manufactured_supernodal_mode_qualification'
        result['error']=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    finally:persist()
    print(json.dumps(dict(status=result['status'],result=pin(out/'RESULT.json'),cost=result['cost'])))
    return 0 if result['status']=='passed_manufactured_supernodal_mode_qualification' else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();raise SystemExit(run(args.config,args.out))
