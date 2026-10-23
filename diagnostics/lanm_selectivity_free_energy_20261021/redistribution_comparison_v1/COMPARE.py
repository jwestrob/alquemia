#!/usr/bin/env python3
"""Saved, matched-parent comparison of quantum and model electric moments."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

BOHR = 0.529177210544
DEBYE = 4.803206802
CONVERSION = BOHR * DEBYE
GIDS = tuple(f'aquo{n}_central_{m}' for n in (8, 9) for m in ('La', 'Dy'))


def require(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verified(record):
    require(pin(record['path']) == record, 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verified(record).read_text())


def contains_pin(value, target):
    if value == target:
        return True
    if isinstance(value, dict):
        return any(contains_pin(v, target) for v in value.values())
    if isinstance(value, list):
        return any(contains_pin(v, target) for v in value)
    return False


def compare(cfg):
    require(cfg['implementation'] == pin(__file__), 'Comparison source changed')
    verified(cfg['declaration'])
    qm = read(cfg['quantum_result'])
    qr = read(cfg['quantum_review'])
    qa = read(cfg['quantum_admission'])
    require(qm['status'] == 'complete_saved_redistribution_diagnostic', 'Quantum incomplete')
    require(qr['status'] == 'passed_independent_saved_redistribution_review', 'Quantum review failed')
    require(qr['actual_result'] == cfg['quantum_result'] and qa['result'] == cfg['quantum_result']
            and qa['independent_actual_review'] == cfg['quantum_review'], 'Quantum joins differ')
    mr = read(cfg['model_result'])
    mc = read(cfg['model_config'])
    mi = read(cfg['model_independent_review'])
    require(mi['status'].startswith('passed') and contains_pin(mi, cfg['model_result']),
            'Model independent review not bound to actual result')
    manifest = read(mr['manifest'])
    spec = importlib.util.spec_from_file_location('admitted_model_moments', verified(mc['implementation']))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    replay = module.admit_result(mc, manifest, mr)
    collection = read(qr['collection'])
    common = read(collection['panel'])
    require(common['pair_panel'] == manifest['panel'], 'Parent panel differs between methods')
    panel = read(manifest['panel'])
    require(set(qm['parents']) == set(qr['parents']) == set(replay['model_cycles']['parents']) == set(GIDS),
            'Four-parent denominator differs')
    geometry_checks = {}
    parents = {}
    for gid in GIDS:
        raw = verified(manifest['rows'][gid]['input']['xyz']).read_text().splitlines()
        xyz = [line.split() for line in raw[1:] if line.strip()]
        g = panel['parent_geometries'][gid]
        require(int(raw[0].split()[0]) == len(xyz) == len(g['symbols']), 'Atom denominator differs')
        require([r[1].lower() for r in xyz] == [x.lower() for x in g['symbols']], 'Atom labels differ')
        geometry_error = max(abs(float(r[2+j]) - pos[j]) for r, pos in zip(xyz, g['positions_A']) for j in range(3))
        require(geometry_error <= 5e-10, 'Parent coordinates differ')
        geometry_checks[gid] = {'max_difference_A': geometry_error, 'bound_A': 5e-10}
        parents[gid] = {}
        for kind in ('C', 'P', 'B'):
            q = qm['parents'][gid]['vectors']['total_dipole_au'][kind]
            qfull = qr['parents'][gid]['vectors']['total_dipole_au'][kind]
            qbound = qm['parents'][gid]['dipole_roundoff_bounds_au']['total_dipole_au'][kind]
            require(len(q) == len(qfull) == 3 and qbound >= 0, 'Quantum vector/bound invalid')
            require(max(abs(a-b) for a,b in zip(q,qfull)) <= qbound, 'Full-precision quantum replay differs')
            model = replay['model_cycles']['parents'][gid][kind]
            x = model['induced_Debye']
            require(len(x) == 3 and model['induced_print_bound_Debye_per_component'] >= 0,
                    'Model vector/bound invalid')
            qd = [v * CONVERSION for v in q]
            qfd = [v * CONVERSION for v in qfull]
            difference = [a-b for a,b in zip(x,qd)]
            bound = model['induced_print_bound_Debye_per_component'] + CONVERSION*qbound
            require(all(math.isfinite(v) for v in x+qd+qfd+difference+[bound]), 'Nonfinite comparison')
            parents[gid][kind] = {
                'quantum_Debye': qd, 'quantum_full_precision_check_Debye': qfd,
                'model_Debye': x, 'model_minus_quantum_Debye': difference,
                'quantum_print_bound_Debye_per_component': CONVERSION*qbound,
                'model_print_bound_Debye_per_component': model['induced_print_bound_Debye_per_component'],
                'combined_print_bound_Debye_per_component': bound,
                'difference_exceeds_print_bounds': [abs(v)>bound for v in difference],
                'quantum_vector_norm_Debye': math.sqrt(math.fsum(v*v for v in qd)),
                'model_vector_norm_Debye': math.sqrt(math.fsum(v*v for v in x)),
                'difference_norm_Debye': math.sqrt(math.fsum(v*v for v in difference)),
            }
        for field in ('quantum_Debye', 'model_Debye', 'model_minus_quantum_Debye'):
            require(max(abs(parents[gid]['C'][field][j]-parents[gid]['P'][field][j]-parents[gid]['B'][field][j])
                        for j in range(3)) <= 1e-10, 'C/P/B arithmetic closure failed')
    return {'status':'complete_saved_quantum_model_dipole_comparison', 'inputs':cfg,
            'Debye_per_atomic_unit':CONVERSION, 'geometry_checks':geometry_checks, 'parents':parents,
            'new_native_calls':0, 'fitting_performed':False, 'physical_acceptance_threshold':None,
            'unique_polarization_or_charge_transfer_attribution':False,
            'limitations':['Print bounds exclude SCF/model/reference uncertainty.',
                           'W8 vector cancellation is symmetry dependent.',
                           'Dipole differences are not exchange energies or polarizabilities.']}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    out=Path(args.output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); cpu=time.process_time()
    try:
        result=compare(json.loads(Path(args.config).read_text()))
        filename='RESULT.json'
    except Exception as exc:
        result={'status':'failed_saved_comparison','error':repr(exc),'implementation':pin(__file__),'new_native_calls':0}
        filename='FAILURE.json'
    result['cost']={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'local_threads':1}
    (out/filename).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'output':str(out)}))
    raise SystemExit(0 if filename=='RESULT.json' else 1)
