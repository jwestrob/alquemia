#!/usr/bin/env python3
"""Independent arithmetic from already admitted saved rows; no native re-admission."""
import hashlib
import json
import math
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[4]
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/redistribution_comparison_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def load(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p
        p = p['path']
    return json.loads(Path(p).read_text())


def main():
    wall = time.monotonic(); cpu = time.process_time()
    cfg = load(W / 'CONFIG.json')
    q = load(cfg['quantum_result']); qr = load(cfg['quantum_review'])
    m = load(cfg['model_result']); mr = load(cfg['model_independent_review'])
    result = load(W / 'actual_v1/RESULT.json')
    assert mr['input_result'] == cfg['model_result']
    assert qr['actual_result'] == cfg['quantum_result']
    assert mr['status'].startswith('passed') and qr['status'].startswith('passed')
    assert result['inputs'] == cfg and result['new_native_calls'] == 0
    conversion = 4.803206802 * 0.529177210544
    assert result['Debye_per_atomic_unit'] == conversion
    checks = 0; max_error = 0.0; parents = {}

    def close(a, b):
        nonlocal checks, max_error
        assert len(a) == len(b) == 3
        e = max(abs(x-y) for x, y in zip(a, b))
        assert all(math.isfinite(v) for v in a+b) and e <= 1e-12
        checks += 1; max_error = max(max_error, e)

    for gid in ('aquo8_central_La', 'aquo8_central_Dy', 'aquo9_central_La', 'aquo9_central_Dy'):
        parents[gid] = {}
        for kind in ('C', 'P', 'B'):
            coeff = m['model_cycles']['parents'][gid][kind]['terms']
            mv = [math.fsum(c * math.fsum(v[j] for v in m['rows'][rid]['moments']['induced_Debye'])
                            for rid, c in coeff.items()) for j in range(3)]
            qterms = q['parents'][gid]['combinations'][kind]
            qv = [math.fsum(c * q['rows'][rid]['total_dipole_au'][j] for c, rid in qterms)
                  * conversion for j in range(3)]
            qb = math.fsum(abs(c) * q['rows'][rid]['dipole_roundoff_bounds_au']['total_dipole_au']
                           for c, rid in qterms) * conversion
            mb = math.fsum(abs(c) * len(m['rows'][rid]['moments']['induced_Debye']) * .5e-6
                           for rid, c in coeff.items())
            r = result['parents'][gid][kind]
            close(mv, r['model_Debye'])
            close(qv, r['quantum_Debye'])
            close(mv, mr['parents'][gid][kind]['induced_Debye'])
            close([x-y for x, y in zip(mv, qv)], r['model_minus_quantum_Debye'])
            assert abs(qb+mb-r['combined_print_bound_Debye_per_component']) <= 1e-15
            assert r['difference_exceeds_print_bounds'] == [abs(x-y)>qb+mb for x,y in zip(mv,qv)]
            checks += 2
            parents[gid][kind] = {'model_Debye': mv, 'quantum_Debye': qv,
                                  'combined_print_bound_Debye_per_component': qb+mb}
        for method in ('model_Debye', 'quantum_Debye'):
            close([parents[gid]['C'][method][j] - parents[gid]['P'][method][j] for j in range(3)],
                  parents[gid]['B'][method])
    paired = {method: [parents['aquo9_central_Dy']['B'][method][j] -
                       parents['aquo9_central_La']['B'][method][j] for j in range(3)]
              for method in ('model_Debye', 'quantum_Debye')}
    out = {'status': 'passed_independent_saved_comparison_arithmetic', 'implementation': pin(__file__),
           'comparison_result': pin(W/'actual_v1/RESULT.json'), 'inputs': cfg,
           'vector_or_bound_checks': checks, 'maximum_arithmetic_difference_Debye': max_error,
           'parents': parents, 'W9_Dy_minus_La_B': paired, 'new_native_calls': 0,
           'raw_admission_repeated': False,
           'limitations': ['Relies on pinned independent raw admissions; this is an arithmetic review.',
                          'Print bounds do not measure reference accuracy or SCF uncertainty.',
                          'W8 near-cancellation is not an accuracy qualification.'],
           'cost': {'wall_seconds': time.monotonic()-wall, 'CPU_seconds': time.process_time()-cpu,
                    'local_threads': 1}}
    path = W/'independent_arithmetic_review_v1/REVIEW.json'
    assert not path.exists()
    path.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'review':pin(path), 'checks':checks, 'paired':paired}))


if __name__ == '__main__':
    main()
