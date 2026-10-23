"""Audit the exact isolated binary installation; no numeric factorization."""
import ctypes
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import platform
import sys
import time


def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    t0, c0 = time.perf_counter(), time.process_time()
    root = Path.cwd().resolve()
    rel = Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_stability_definition_v1/sparse_factorization_design_v1')
    d = root / 'diagnostics/lanm_selectivity_free_energy_20261021' / rel
    w = root / 'workspaces/lanm_selectivity_free_energy_20261021' / rel
    env = w / 'environment_v1'
    prefix = root / 'workspaces/_envs/lanm_cholmod_v1'
    assert Path(sys.prefix).resolve() == prefix
    thread_keys = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
                   'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
    assert all(os.environ.get(k) == '1' for k in thread_keys)
    plan_path = env / 'attempt2/EXACT_PLAN.json'
    packages = json.loads(plan_path.read_text())['packages']
    installed = {}
    for p in sorted((prefix / 'conda-meta').glob('*.json')):
        r = json.loads(p.read_text())
        installed[r['name']] = (r, p)
    assert len(installed) == len(packages) == 56
    records = []
    for wanted in packages:
        actual, metadata = installed[wanted['name']]
        fields = ('name', 'version', 'build', 'url', 'md5', 'sha256')
        assert all(actual[k] == wanted[k] for k in fields), wanted['name']
        records.append(dict(identity={k: actual[k] for k in fields}, metadata=pin(metadata)))
    import numpy
    import scipy
    from sksparse import cholmod
    openblas = prefix / 'lib/libopenblasp-r0.3.34.so'
    omp = prefix / 'lib/libgomp.so.1'
    b, g = ctypes.CDLL(str(openblas)), ctypes.CDLL(str(omp))
    b.openblas_get_num_threads.restype = ctypes.c_int
    g.omp_get_max_threads.restype = ctypes.c_int
    assert b.openblas_get_num_threads() == g.omp_get_max_threads() == 1
    files = [Path(sys.executable), Path(cholmod.__file__),
             Path(cholmod.__file__).with_name('cholmod.pxd'),
             prefix / 'include/suitesparse/cholmod.h', prefix / 'lib/libcholmod.so.5.3.1',
             openblas, omp]
    api = {}
    for name in ('factorize', 'get_factor', 'get_perm', 'solve', 'update', 'downdate', 'change_factor'):
        obj = getattr(cholmod.CholeskyFactor, name)
        api[name] = dict(signature=str(inspect.signature(obj)), documentation=obj.__doc__)
    assert not hasattr(cholmod.CholeskyFactor, 'minor')
    assert not hasattr(cholmod.CholeskyFactor, 'status')
    receipts = [env / 'resolve_v1/RECEIPT.json', env / 'attempt2/resolve_v1/RECEIPT.json',
                env / 'attempt2/install_v1/RECEIPT.json', env / 'attempt3/install_v1/RECEIPT.json']
    costs = [dict(receipt=pin(p), **json.loads(p.read_text())['cost']) for p in receipts]
    out = dict(status='passed_exact_isolated_environment_audit', implementation=pin(__file__),
               prefix=str(prefix), exact_plan=pin(plan_path), packages=records,
               runtime_files=[pin(p) for p in files],
               versions=dict(python=platform.python_version(), numpy=numpy.__version__,
                             scipy=scipy.__version__, scikit_sparse=importlib.metadata.version('scikit-sparse')),
               settings=dict(lower=True, order='amd', sym_kind='sym', supernodal_mode='simplicial',
                             initial_numeric_ldl=False, numeric_beta=0.0,
                             update_representation='explicit change_factor(kind="LDL"); inspect positive D before LL readback'),
               threads=dict(environment={k: os.environ[k] for k in thread_keys},
                            actual_openblas=b.openblas_get_num_threads(), actual_openmp=g.omp_get_max_threads()),
               api=api, public_minor_available=False, public_status_available=False,
               failure_contract='Record documented exceptions and partial public metadata. Positive admission additionally requires finite positive factors and explicit reconstruction/solve residuals. No private pointer access.',
               installation_attempt_costs=costs,
               summed_attempt_wall_seconds=sum(r['wall_seconds'] for r in costs),
               summed_child_CPU_seconds=sum(r['child_user_CPU_seconds'] + r['child_system_CPU_seconds'] for r in costs),
               native_calls=0, factorizations=0,
               cost=dict(wall_seconds=time.perf_counter()-t0, CPU_seconds=time.process_time()-c0, numerical_threads=1))
    p = env / 'ENVIRONMENT_AUDIT.json'
    with p.open('x') as f:
        json.dump(out, f, indent=2); f.write('\n')
    print(json.dumps(pin(p)))


if __name__ == '__main__':
    main()
