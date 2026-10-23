#!/usr/bin/env python3
"""Replay only the terminal two-cell mapping; independently derive Urey E/gradient."""
import datetime as dt
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_numerical_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    wall, cpu = time.monotonic(), time.process_time()
    out = WORK/'actual_review_v1'; out.mkdir(parents=True, exist_ok=False)
    checks = []
    def check(ok, message):
        if not ok:
            raise ValueError(message)
        checks.append(message)
    def verify(p):
        check(pin(p['path']) == p, 'Pin '+p['path'])
        return Path(p['path'])
    rp = WORK/'run_v1/execution/RESULT.json'
    result = json.loads(rp.read_text()); cfg = json.loads(verify(result['config']).read_text())
    release = json.loads(verify(cfg['root_release']).read_text())
    verify(release['independent_review']); verify(release['preflight'])
    check(release['implementation'] == cfg['implementation'], 'Exact released source')
    base = {k: v for k, v in cfg.items() if k != 'root_release'}
    check(hashlib.sha256(json.dumps(base, sort_keys=True).encode()).hexdigest() == release['config_sha256_without_release'], 'Exact released configuration')
    for p in cfg['authorities']:
        verify(p)
    for p in cfg['executables'].values():
        verify(p)
    verify(cfg['implementation']); helper_path = verify(cfg['native_helper'])
    spec = importlib.util.spec_from_file_location('strict_original_parser', helper_path)
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    check(result['status'] == 'mapping_controls_passed' and result['error'] is None and result['slurm_job'] == '17668', 'Terminal successful actual job')
    check(result['native_programs'] == 4 and len(result['rows']) == len(cfg['rows']) == 2, 'Exact finite program and cell denominator')
    check(not result['bulk_physics_qualified'] and not result['protein_qualified'], 'No bulk or protein qualification')
    evidence = []; receipt_count = 0; native_wall = native_cpu = 0.
    for actual_row, row in zip(result['rows'], cfg['rows']):
        gid = row['id']; folder = rp.parent/gid
        check(actual_row['id'] == gid and actual_row['status'] == 'passed', 'Complete row '+gid)
        for field, filename in [('xyz', 'system.xyz'), ('key', 'system.key'), ('parameters', row['parameter_filename'])]:
            verify(row[field])
            check(pin(folder/filename)['sha256'] == row[field]['sha256'], 'Actual input bytes '+gid+' '+field)
        for field in ('reference_xyz', 'reference_key', 'reference_parameters'):
            verify(row[field])
        for field in ('analyze_output', 'gradient_output'):
            verify(row['reference_native'][field])
        parsed = helper.parse(folder/'analyze.out', folder/'testgrad.out', {'atoms': 25, 'derivative_test': False})
        reference = helper.parse(row['reference_native']['analyze_output']['path'], row['reference_native']['gradient_output']['path'], {'atoms': 25, 'derivative_test': False})
        check(parsed == actual_row['native'], 'Exact raw BIO replay '+gid)
        check(reference == actual_row['reference'], 'Exact raw saved SI replay '+gid)
        check(len(actual_row['receipts']) == 2, 'Two actual receipts '+gid)
        for receipt, name, args in zip(actual_row['receipts'], ['analyze', 'testgrad'], [['system.xyz', 'E'], ['system.xyz', 'Y', 'N']]):
            verify(receipt['output']); receipt_count += 1
            check(receipt['argv'] == [cfg['executables'][name]['path']]+args and receipt['cwd'] == str(folder), 'Exact native invocation '+gid+' '+name)
            check(receipt['returncode'] == 0 and receipt['timeout_seconds'] == 300 and set(receipt['threads'].values()) == {'1'}, 'Native success and resources '+gid+' '+name)
            native_wall += receipt['elapsed_seconds']
            native_cpu += receipt['child_user_cpu_seconds']+receipt['child_system_cpu_seconds']
        atoms = [line.split() for line in (folder/'system.xyz').read_text().splitlines()[1:]]
        check(len(atoms) == 25, 'Full25 atom geometry '+gid)
        xyz = [[float(x) for x in a[2:5]] for a in atoms]
        k = -7.6; rsi = 1.5337; rbio = 1.5537; correction = 0.; grad = [[0.]*3 for a in atoms]; pairs = []
        for a in atoms:
            if int(a[5]) != 402:
                continue
            i, j = [int(x)-1 for x in a[6:]]; pairs.append((i, j))
            separation = [xyz[i][c]-xyz[j][c] for c in range(3)]
            distance = math.sqrt(sum(v*v for v in separation))
            # Expanded polynomial and separate native radial derivatives are
            # independent of the runner's difference-of-squares implementation.
            correction += k*(2*distance*(rsi-rbio)+rbio*rbio-rsi*rsi)
            derivative_si = 2*k*(distance-rsi)
            derivative_bio = 2*k*(distance-rbio)
            for c in range(3):
                delta = (derivative_bio-derivative_si)*separation[c]/distance
                grad[i][c] += delta; grad[j][c] -= delta
        check(len(pairs) == 8 and len(set(i for pair in pairs for i in pair)) == 16, 'Eight distinct Urey pairs '+gid)
        saved = actual_row['comparison']; old = saved['known_urey_difference']
        check(abs(correction-old['energy_kcal_mol']) < 1e-12 and max(abs(x-y) for a,b in zip(grad,old['gradient_kcal_mol_A']) for x,y in zip(a,b)) < 1e-12, 'Independent analytic correction agrees '+gid)
        de = parsed['total_energy_kcal_mol']-reference['total_energy_kcal_mol']
        error = de-correction
        residual = [[a-b-c for a,b,c in zip(aa,bb,cc)] for aa,bb,cc in zip(parsed['gradient_kcal_mol_A'], reference['gradient_kcal_mol_A'], grad)]
        max_gradient = max(abs(x) for rowg in residual for x in rowg)
        check(abs(error) <= 1e-4 and max_gradient <= 1e-3, 'Independent full E/gradient parity '+gid)
        check(abs(error-saved['energy_remainder_kcal_mol']) < 1e-12 and abs(max_gradient-saved['max_gradient_remainder_kcal_mol_A']) < 1e-12, 'Recorded numerical comparison replay '+gid)
        ac, rc = parsed['components'], reference['components']
        check(set(ac) == set(rc), 'Component denominator '+gid)
        components = {}
        for name in ac:
            check(ac[name]['interactions'] == rc[name]['interactions'], 'Component interaction counts '+gid+' '+name)
            delta = ac[name]['energy_kcal_mol']-rc[name]['energy_kcal_mol']
            components[name] = delta-(correction if name == 'Urey-Bradley' else 0.)
            check(abs(components[name]) <= 1e-4, 'Independent component parity '+gid+' '+name)
        evidence.append({'id':gid, 'BIO_minus_SI_energy':de, 'independent_Urey_energy':correction,
            'energy_remainder':error, 'max_gradient_remainder':max_gradient,
            'full25_gradient_remainder_kcal_mol_A':residual,
            'component_remainders':components, 'BIO_outputs':{k:parsed[k] for k in ('analyze_output','gradient_output')},
            'saved_SI_outputs':row['reference_native']})
    check(receipt_count == 4, 'All four receipts retained')
    receipt_path = WORK/'run_v1/completion_event.json'; blob = receipt_path.read_bytes()
    (out/'completion_event.snapshot.json').write_bytes(blob)
    scheduler = json.loads(blob)['records']['17668']
    check(scheduler['state']=='COMPLETED' and scheduler['exit_code']=='0:0', 'Actual scheduler completion')
    review = {'status':'actual_mapping_review_passed', 'job_id':'17668', 'utc':dt.datetime.now(dt.timezone.utc).isoformat(),
        'implementation':pin(__file__), 'actual_result':pin(rp), 'config':result['config'], 'root_release':cfg['root_release'],
        'source_review':release['independent_review'], 'completion_snapshot':pin(out/'completion_event.snapshot.json'),
        'checks':checks, 'check_count':len(checks), 'cells':evidence, 'native_programs_replayed':4,
        'saved_SI_cells_replayed':2, 'native_cost':{'scheduler_elapsed_seconds':int(scheduler['elapsed_seconds']),
        'allocated_CPUs':int(scheduler['allocated_cpus']), 'allocated_CPU_seconds':int(scheduler['elapsed_seconds'])*int(scheduler['allocated_cpus']),
        'summed_native_wall_seconds':native_wall, 'measured_child_CPU_seconds':native_cpu,
        'adapter_wall_seconds':result['elapsed_seconds']},
        'new_native_calls':0, 'repeated_controls':0, 'local_threads':1, 'review_wall_seconds':time.monotonic()-wall,
        'review_CPU_seconds':time.process_time()-cpu, 'finite_mapping_branch_closed':True,
        'bulk_physics_qualified':False, 'protein_qualified':False, 'reserved_donor_outcomes_read':False}
    (out/'REVIEW.json').write_text(json.dumps(review, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'review':pin(out/'REVIEW.json'), 'checks':len(checks), 'cost':review['native_cost'],
        'max_gradient_remainder':max(v['max_gradient_remainder'] for v in evidence)}))


if __name__ == '__main__':
    main()
