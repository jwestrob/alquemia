#!/usr/bin/env python3
"""Read-only origin-map arithmetic diagnosis; bind is metadata, run is Slurm-only."""
import argparse
from collections import Counter
import contextlib
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import platform
import struct
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
WORK = ROOT/'workspaces/lanm_coupled_mechanics_20261017'
HERE = Path(__file__).resolve().parent
PROTOCOL = 'lbt_hybrid_exact_origin_mapping_diagnosis_v1'
EXPECTED = [('7CCO', 'baseline'), ('7CCO', 'expanded'), ('7CCN', 'baseline'), ('7CCN', 'expanded')]
THREADS = {'OMP_NUM_THREADS': '12', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
           'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve(strict=True)
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def bind(out):
    out.mkdir(parents=True, exist_ok=False)
    config = pin(WORK/'lbt_hybrid_scout_v1/CONFIG.json')
    c = read(config)
    interpreter = ROOT/'workspaces/mace_hybrid_20260916/software_v1/venv/bin/python'
    cfg = {'protocol': PROTOCOL, 'implementation': pin(__file__), 'spec': pin(HERE/'SPEC.md'),
           'hybrid_config': config, 'failed_result': pin(WORK/'lbt_hybrid_scout_v1/execution/RESULT.json'),
           'passed_preflight': pin(WORK/'lbt_hybrid_preflight_v1/execution/RESULT.json'),
           'failed_launcher': pin(WORK/'lbt_hybrid_scout_v1/RUN.sbatch'),
           'passed_launcher': pin(WORK/'lbt_hybrid_preflight_v1/RUN.sbatch'),
           'interpreter_literal': str(interpreter), 'interpreter': pin(interpreter),
           'expected_host': 'node-224-2t-8gpu-1', 'expected_thread_environment': THREADS,
           'declared_source_regions': [list(x) for x in EXPECTED], 'repeats_per_origin': 2,
           'new_molecular_calls': 0, 'model_loads': 0, 'tolerance_changes': False}
    for label, dirname in [('failed', 'lbt_hybrid_scout_v1'), ('passed', 'lbt_hybrid_preflight_v1')]:
        # Watcher receipts can still update delivery fields. Freeze complete bytes now.
        raw = (WORK/dirname/'completion_event.json').read_bytes()
        target = out/(label+'_completion_snapshot.json')
        with target.open('xb') as stream:
            stream.write(raw)
        cfg[label+'_completion_snapshot'] = pin(target)
    for key in ('implementation', 'preparation_result', 'preparation_implementation',
                'preparation_config', 'gradient_helper', 'common_helper', 'GFN_backend'):
        verify(c[key])
        cfg['hybrid_'+key] = c[key]
    return write(out/'CONFIG.json', cfg)


def ordered_float_bits(value):
    bits = struct.unpack('>Q', struct.pack('>d', value))[0]
    return (0xffffffffffffffff-bits) if bits >> 63 else bits+0x8000000000000000


def leaf_differences(old, new, path='$'):
    """Enumerate every structural, type or unequal scalar leaf; no tolerance."""
    rows = []
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new)):
            p = path+'/'+str(key).replace('~', '~0').replace('/', '~1')
            if key not in old or key not in new:
                rows.append({'path': p, 'kind': 'missing_key', 'old_present': key in old,
                             'new_present': key in new, 'old': old.get(key), 'new': new.get(key)})
            else:
                rows.extend(leaf_differences(old[key], new[key], p))
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            rows.append({'path': path, 'kind': 'list_length', 'old_length': len(old), 'new_length': len(new)})
        for i in range(max(len(old), len(new))):
            if i >= len(old) or i >= len(new):
                rows.append({'path': path+'/'+str(i), 'kind': 'missing_item',
                             'old': old[i] if i < len(old) else None,
                             'new': new[i] if i < len(new) else None})
            else:
                rows.extend(leaf_differences(old[i], new[i], path+'/'+str(i)))
    elif type(old) is not type(new) or old != new:
        row = {'path': path, 'kind': 'scalar', 'old': old, 'new': new,
               'old_type': type(old).__name__, 'new_type': type(new).__name__}
        if type(old) in (int, float) and type(new) in (int, float):
            a, b = float(old), float(new)
            if math.isfinite(a) and math.isfinite(b):
                row.update(absolute_difference=abs(b-a), signed_difference=b-a,
                           old_hex=a.hex(), new_hex=b.hex(),
                           ulp_distance=abs(ordered_float_bits(a)-ordered_float_bits(b)))
        rows.append(row)
    return rows


def runtime(np, cfg):
    text = io.StringIO()
    with contextlib.redirect_stdout(text):
        np.show_config()
        if hasattr(np, 'show_runtime'):
            np.show_runtime()
    core = getattr(np, '_core', getattr(np, 'core', None))
    dispatch = getattr(core, '_multiarray_umath', None)
    return {'host': platform.node(), 'slurmd_nodename': os.environ.get('SLURMD_NODENAME'),
            'platform': platform.platform(), 'machine': platform.machine(),
            'python_version': sys.version, 'sys_executable_literal': sys.executable,
            'sys_executable': pin(sys.executable), 'sys_prefix': sys.prefix, 'sys_base_prefix': sys.base_prefix,
            'numpy_version': np.__version__, 'numpy_module': pin(np.__file__),
            'numpy_configuration': text.getvalue(),
            'numpy_CPU_features': getattr(dispatch, '__cpu_features__', None),
            'numpy_CPU_baseline': getattr(dispatch, '__cpu_baseline__', None),
            'numpy_CPU_dispatch': getattr(dispatch, '__cpu_dispatch__', None),
            'environment': {k: os.environ[k] for k in sorted(os.environ)
                            if k.startswith(('SLURM_', 'SLURMD_', 'OMP_', 'MKL_', 'OPENBLAS_', 'BLIS_', 'NUMEXPR_', 'NPY_'))},
            'cpuinfo': Path('/proc/cpuinfo').read_text() if Path('/proc/cpuinfo').exists() else None,
            'processor_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None}


def cap_primitives(np, partition, x, archived, regenerated):
    ids = {a['id']: i for i, a in enumerate(partition['physical_atoms'])}
    rows = []
    for c, a, b in zip(partition['caps'], archived['caps'], regenerated['caps']):
        ia, ib = ids[c['retained_source_id']], ids[c['omitted_source_id']]
        v = x[ib]-x[ia]
        values = [float(t) for t in v]
        norms = {'numpy_linalg_norm': float(np.linalg.norm(v)),
                 'numpy_sqrt_dot': float(np.sqrt(np.dot(v, v))),
                 'numpy_sqrt_sum': float(np.sqrt(np.sum(v*v))),
                 'python_sqrt_fsum': math.sqrt(math.fsum(t*t for t in values)),
                 'python_hypot': math.hypot(*values)}
        rows.append({'definition': c, 'retained_index': ia, 'omitted_index': ib,
                     'retained_xyz_A': x[ia].tolist(), 'omitted_xyz_A': x[ib].tolist(),
                     'vector_A': values, 'vector_hex': [v.hex() for v in values],
                     'norms_A': norms, 'norms_hex': {k: v.hex() for k, v in norms.items()},
                     'archived_cap': a, 'regenerated_cap': b})
    return rows


def run(cfg, out):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    require(cfg['protocol'] == PROTOCOL and cfg['declared_source_regions'] == [list(x) for x in EXPECTED]
            and cfg['repeats_per_origin'] == 2, 'Diagnostic declaration differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Diagnostic implementation differs')
    require(pin(sys.executable)['sha256'] == cfg['interpreter']['sha256'], 'Python binary differs')
    require(Path(sys.prefix).resolve() == Path(cfg['interpreter_literal']).parent.parent.resolve(), 'Python environment differs')
    host = os.environ.get('SLURMD_NODENAME', platform.node()).split('.')[0]
    require(host == cfg['expected_host'], 'Must diagnose the actual failed host')
    require(all(os.environ.get(k) == v for k, v in cfg['expected_thread_environment'].items()), 'Thread environment differs')
    for key, value in cfg.items():
        if isinstance(value, dict) and set(value) == {'path', 'sha256'}:
            verify(value)
    failed, passed = read(cfg['failed_result']), read(cfg['passed_preflight'])
    c = read(cfg['hybrid_config'])
    require(failed['status'] == 'execution_failure' and failed['slurm_job_id'] == '12975'
            and failed['failure'] == "ValueError('Exact origin mapping replay differs')", 'Actual failure differs')
    require(passed['status'] == 'complete' and passed['stage'] == 'preflight' and passed['slurm_job_id'] == '12972', 'Actual preflight differs')
    require(failed['config']['sha256'] == passed['config']['sha256'] == cfg['hybrid_config']['sha256'], 'Preflight/scout configs differ')
    for key in ('implementation', 'preparation_result', 'preparation_implementation',
                'preparation_config', 'gradient_helper', 'common_helper', 'GFN_backend'):
        require(c[key] == cfg['hybrid_'+key], 'Bound helper differs: '+key)
    prep = read(c['preparation_result'])
    require(prep['status'] == 'complete' and prep['slurm_job_id'] == '12953'
            and prep['new_molecular_calls'] == prep['new_MM_energy_calls'] == 0
            and [(r['source_id'], r['region']) for r in prep['rows']] == EXPECTED, 'Preparation denominator differs')
    initial_environment = {k: os.environ.get(k) for k in list(THREADS)+['VECLIB_MAXIMUM_THREADS']}
    import numpy as np
    sys.path.insert(0, str(verify(c['common_helper']).parent))
    helper = module('diagnostic_original_mapping', c['preparation_implementation'])
    gradient = module('diagnostic_original_cap_gradient', c['gradient_helper'])
    helper.np = np
    helper.cap_jacobians = gradient.cap_jacobians
    # The actual constructor imports this definition-only backend after NumPy.
    # Its module-level THREAD_KEYS loop sets OMP/BLAS variables to 1. Reproduce
    # that existing side effect without constructing or calling any calculator.
    module('diagnostic_original_GFN_backend', c['GFN_backend'])
    runtime_info = runtime(np, cfg)
    runtime_info['launch_thread_environment_before_imports'] = initial_environment
    runtime_info['GFN_backend_import'] = c['GFN_backend']
    rt = write(out/'RUNTIME.json', runtime_info)
    result = {'status': 'complete', 'protocol': PROTOCOL, 'runtime': rt, 'rows': [],
              'new_molecular_calls': 0, 'model_loads': 0, 'new_physical_states': 0,
              'tolerance_changes': False, 'diagnostic_is_admission': False}
    for row in prep['rows']:
        entry = {'source_id': row['source_id'], 'region': row['region'], 'status': 'unavailable',
                 'preparation': row['preparation']}
        result['rows'].append(entry)
        try:
            p = read(row['preparation'])
            partition, archived = read(p['partition']), read(p['mapping'])
            x = np.asarray([a['xyz_A'] for a in partition['physical_atoms']], float)
            require(x.shape == (256, 3) and np.isfinite(x).all(), 'Physical coordinates invalid')
            fresh = [helper.update_geometry(partition, x) for _ in range(2)]
            d = leaf_differences(archived, fresh[0])
            repeat = leaf_differences(fresh[0], fresh[1])
            name = row['source_id']+'_'+row['region']
            candidate_pins = [write(out/(name+'_mapping_repeat%d.json'%i), value) for i, value in enumerate(fresh)]
            differences = write(out/(name+'_DIFFERENCES.json'), {'archived_vs_recomputed': d, 'repeat_vs_repeat': repeat})
            primitive = write(out/(name+'_CAP_PRIMITIVES.json'), cap_primitives(np, partition, x, archived, fresh[0]))
            numeric = [r for r in d if 'absolute_difference' in r]
            xyz = [r for r in numeric if '/xyz_A/' in r['path']]
            jac = [r for r in numeric if '/jacobian_' in r['path']]
            entry.update(status='compared', partition=p['partition'], archived_mapping=p['mapping'],
                         regenerated_mappings=candidate_pins, differences=differences, cap_primitives=primitive,
                         original_dictionary_equality=(archived == fresh[0]), repeated_dictionary_equality=(fresh[0] == fresh[1]),
                         typed_leaf_difference_count=len(d), repeated_typed_leaf_difference_count=len(repeat),
                         difference_kind_counts=dict(Counter(v['kind'] for v in d)),
                         difference_top_level_counts=dict(Counter(v['path'].split('/')[1] for v in d)),
                         max_numeric_absolute_difference=max((v['absolute_difference'] for v in numeric), default=0.),
                         max_ulp_distance=max((v['ulp_distance'] for v in numeric), default=0),
                         max_coordinate_absolute_difference_A=max((v['absolute_difference'] for v in xyz), default=0.),
                         max_jacobian_absolute_difference=max((v['absolute_difference'] for v in jac), default=0.))
        except Exception as exc:
            entry.update(failure=repr(exc), traceback=traceback.format_exc())
            result['status'] = 'partial_diagnostic'
    result['compared_origins'] = sum(r['status'] == 'compared' for r in result['rows'])
    result['declared_origins'] = 4
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['bind', 'run'], required=True)
    parser.add_argument('--config')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output).resolve()
    if args.stage == 'bind':
        print(json.dumps(bind(out)))
        return
    require(args.config is not None, 'Run requires config')
    cp = pin(args.config)
    cfg = read(cp)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    try:
        result = run(cfg, out)
    except Exception as exc:
        result = {'status': 'execution_failure', 'failure': repr(exc), 'traceback': traceback.format_exc(),
                  'new_molecular_calls': 0, 'model_loads': 0}
    result.update(config=cp, implementation=pin(__file__), slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                  elapsed_seconds=time.monotonic()-started)
    write(out/'RESULT.json', result)
    if result['status'] != 'complete':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
