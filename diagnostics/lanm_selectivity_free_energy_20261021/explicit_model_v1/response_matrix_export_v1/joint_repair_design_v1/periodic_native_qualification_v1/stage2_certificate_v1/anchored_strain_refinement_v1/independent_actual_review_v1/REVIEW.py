"""Saved receipt/summary review only; no numerical helper or native imports."""
import hashlib
import json
import math
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
PANEL = HERE.parent
DIRECTIONS = tuple('strain_' + str(i) for i in range(3))
LABELS = ('1.25e-05', '6.25e-06')
STEPS = (0.0000125, 0.00000625)
IDS = tuple(d + '__' + h + '__' + s for d in DIRECTIONS for h in LABELS for s in ('minus', 'plus'))
COUNTS = dict(native_programs=1, gradient_calls=1, energy_only_calls=0, extra_operator_applications=0)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return dict(path=str(p), sha256=digest)


def verify(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'} and pin(item['path']) == item, 'Changed pin')
    return Path(item['path'])


def read(item):
    return json.loads(verify(item).read_text())


def run(config, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = dict(status='failed_anchored_strain_saved_review', implementation=pin(__file__), configuration=pin(config), error=None,
                  native_calls=0, solver_calls=0, arithmetic_replays=0, control_reruns=0, active_response_admission=False)
    try:
        cfg = json.loads(Path(config).read_text())
        require(cfg['implementation'] == pin(__file__), 'Frozen reviewer')
        frozen = read(cfg['inputs'])
        require(cfg['inputs'] == pin(HERE / 'INPUTS.json'), 'Frozen review inputs')
        for item in frozen.values():
            verify(item)
        source = read(frozen['source_review'])
        r = read(cfg['result'])
        c = read(r['configuration'])
        release = read(r['release'])
        m = read(c['manifest'])
        result.update(result=cfg['result'], producer_configuration=r['configuration'], producer_release=r['release'], source_review=frozen['source_review'])
        require(r['implementation'] == c['implementation'] == frozen['collector'], 'Frozen collector')
        require(c['status'] == 'bound_anchored_strain_collection' and c['inputs'] == frozen['collector_inputs'], 'Collector input contract')
        require(release['status'] == 'released_anchored_strain_collection' and release['implementation'] == frozen['collector'] and release['configuration'] == r['configuration'], 'Exact collection release')
        require(Path(release['output']).resolve() == Path(cfg['result']['path']).resolve().parent, 'Released collector output')
        require(r['source_review'] == release['source_review'] == frozen['source_review'] and source['status'] == 'passed_anchored_strain_refinement_source_review', 'Admitted source')
        require(c['manifest'] == r['manifest'] == source['manifest'] == frozen['manifest'], 'Prepared panel identity')
        require(r['complete'] and not r['failures'] and r['error'] is None, 'Complete collection required; preserve incomplete evidence')
        require(r['counts'] == dict(expected_new_gradients=12, admitted_new_gradients=12, reused_central_gradients=1, corrected_checks_completed=9, raw_checks_completed=9), 'Exact 12/1/9/9 denominators')
        require(tuple(r['rows']) == tuple(c['results']) == tuple(m['rows']) == IDS, 'Exact twelve endpoints')
        require(r['central_result'] == m['central_result'] and r['central_gradient'] == m['central_gradient'], 'Exact historical center reuse')
        center = read(r['central_result'])
        verify(r['central_gradient'])
        require(center['gradient_output'] == r['central_gradient'] and center['state_id'] == 'central', 'Central receipt/output identity')
        require(r['build'] == c['build'] == m['build'] and r['build_root_review'] == c['build_root_review'] == m['build_root_review'], 'Unchanged build')
        build = read(r['build'])
        verify(build['executables']['gradient_old'])
        read(r['build_root_review'])
        subs = read(c['submissions'])
        require(tuple(x['state_id'] for x in subs) == IDS and len({str(x['job_id']) for x in subs}) == 12, 'Twelve distinct allocations')
        for sid, sub in zip(IDS, subs):
            row = r['rows'][sid]
            require(row['result'] == c['results'][sid], 'Endpoint pin join')
            rec = read(row['result'])
            require(rec['status'] == 'stage2_gradient_complete' and rec['error'] is None and rec['state_id'] == sid and rec['role'] == 'BYPASS_OLD_ABSENT', 'Successful inactive endpoint')
            require(rec['implementation'] == source['implementation'] and rec['manifest'] == c['manifest'] and rec['build'] == c['build'] and rec['source_manifest'] == m['native_sources'], 'Native source lineage')
            require(rec['configuration'] == sub['configuration'] and rec['release'] == sub['release'] and rec['inputs'] == m['rows'][sid]['inputs'] and rec['source_inputs'] == m['rows'][sid]['source_inputs'], 'Submitted exact input identity')
            nc, nr = read(rec['configuration']), read(rec['release'])
            require(nc['state_id'] == nr['state_id'] == sid and nc['operation'] == nr['operation'] == 'gradient' and nr['maximum_native_programs'] == 1 and nr['configuration'] == rec['configuration'], 'One released gradient')
            require(nr['source_review'] == frozen['source_review'] and nr['status'] == 'released_anchored_strain_refinement', 'Endpoint source/release')
            require(Path(row['result']['path']).resolve() == Path(m['rows'][sid]['output']).resolve() / 'RESULT.json' and Path(nr['output']).resolve() == Path(m['rows'][sid]['output']).resolve(), 'Immutable output')
            require(rec['counts'] == COUNTS and rec['native_attempts'] == 1 and rec['certificate'] is None and rec['automatic_retry'] is False, 'One attempt, no extra operation')
            alloc = rec['allocation']
            require(str(alloc['job_id']) == str(sub['job_id']) == str(row['job_id']) and alloc['allocated_CPUs'] == alloc['nodes'] == alloc['tasks'] == 1, 'Actual one-CPU allocation')
            cmd = rec['commands']['native']
            require(cmd['status'] == 'complete' and cmd['returncode'] == 0 and cmd['argv'] == [build['executables']['gradient_old']['path'], 'system.xyz'], 'Unchanged actual binary')
            require(cmd['wall_seconds_limit'] == cmd['CPU_seconds_limit'] == 3600 and cmd['AS_bytes'] == 3 * 1024**3, 'Inherited child limits')
            require(row['commands'] == rec['commands'] and row['cost'] == rec['cost'] and row['gradient_output'] == rec['gradient_output'], 'Literal collector receipt/cost copies')
            for item in [*rec['inputs'].values(), cmd['stdout'], cmd['stderr']]:
                verify(item)
            with verify(row['gradient_output']).open() as f:
                header = [f.readline().strip() for _ in range(4)]
            require(header[0] == 'STAGE2_ONE_ANALYTIC_GRADIENT_V1' and header[2].split() == ['TOTAL_ENERGY', row['energy_token']] and header[3].split() == ['POLARIZATION_ENERGY', row['polarization_token']], 'Literal energy-token header join')
        require(tuple(r['directions']) == DIRECTIONS, 'All three axes')
        directions = {d['id']: d for d in read(m['directions'])['directions']}
        summary = {}
        for name in DIRECTIONS:
            d = r['directions'][name]
            raw, corrected = d['raw_virial_only_checks'], d['corrected_anchored_checks']
            require(d['original_historical_failures_preserved'] is True, 'Historical failure retained')
            for label, checks in [('raw', raw), ('corrected', corrected)]:
                require(tuple(s['half_step_A'] for s in checks['steps']) == STEPS, 'Declared two steps')
                require(checks['passed'] == (all(s['check']['passed'] for s in checks['steps']) and checks['two_step']['passed']), 'Three-check aggregation')
                for h, s in zip(LABELS, checks['steps']):
                    require(s['energy_minus_token'] == r['rows'][name+'__'+h+'__minus']['energy_token'] and s['energy_plus_token'] == r['rows'][name+'__'+h+'__plus']['energy_token'], 'Literal FD input tokens')
                    require(s['analytic_kcal_mol_A'] == float(d['projection'][label if label == 'raw' else 'corrected']) and s['analytic_print_bound'] == float(d['projection'][label+'_print']), 'Saved analytic/printing projection join')
                    require(all(math.isfinite(s[k]) and s[k] >= 0 for k in ('energy_print_bound', 'analytic_print_bound', 'printing_bound')), 'Finite nonnegative print bounds')
            require([s['FD_kcal_mol_A'] for s in raw['steps']] == [s['FD_kcal_mol_A'] for s in corrected['steps']], 'Same FD slopes for raw/corrected comparison')
            summary[name] = dict(anchor_A=directions[name]['anchor_A'], strain_matrix=directions[name]['strain_matrix'], L0_A=directions[name]['L0_A'], projection=d['projection'], raw=raw, corrected=corrected)
        corrected_pass = sum(s['check']['passed'] for d in summary.values() for s in d['corrected']['steps']) + sum(d['corrected']['two_step']['passed'] for d in summary.values())
        raw_pass = sum(s['check']['passed'] for d in summary.values() for s in d['raw']['steps']) + sum(d['raw']['two_step']['passed'] for d in summary.values())
        require(r['numerical_checks_passed'] == (corrected_pass == 9) and r['status'] == ('anchored_strain_checks_passed' if corrected_pass == 9 else 'anchored_strain_numerical_failed'), 'Numerical result distinct from completion')
        accounting = verify(r['scheduler_accounting'])
        require(r['scheduler_accounting'] == c['scheduler_accounting'], 'Scheduler source')
        expected_jobs = {str(s['job_id']) for s in subs}
        parents = {}
        for line in accounting.read_text().splitlines():
            if not line.strip():
                continue
            job, name, state, code, seconds, cpus, measured = line.split('|')
            require(job.split('.')[0] in expected_jobs, 'No unrelated allocation')
            if '.' not in job:
                require(job not in parents and state == 'COMPLETED' and code == '0:0' and int(cpus) == 1 and int(seconds) >= 0, 'Unique successful one-CPU parent')
                parents[job] = int(seconds) * int(cpus)
        sc = r['new_gradient_scheduler_cost']
        require(set(parents) == expected_jobs and sc['parent_jobs'] == sorted(parents) and sc['allocation_count'] == 12 and sc['allocated_CPU_seconds'] == sum(parents.values()) and sc['batch_or_extern_not_added'], 'Parent-only cost arithmetic')
        require(all(r[k] == 0 for k in ('new_native_calls', 'new_FFTs', 'new_factorizations', 'new_scheduler_queries')) and not r['physical_accuracy_qualified'] and not r['active_response_qualified'], 'Saved-only collection and limited admission scope')
        result.update(status='passed_anchored_strain_saved_review', numerical_checks_passed=r['numerical_checks_passed'], counts=dict(new_native_programs=12, reused_center=1, corrected_passed=corrected_pass, corrected_total=9, raw_passed=raw_pass, raw_total=9), directions=summary, scheduler_cost=sc, collector_cost=r['collector_cost'], reused_center=r['central_gradient'], native_binary=build['executables']['gradient_old'])
    except Exception as exc:
        result['error'] = repr(exc)
    result['review_cost'] = dict(wall_seconds=time.monotonic()-tick, CPU_seconds=time.process_time()-cpu, threads=1)
    with (out/'REVIEW.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(run(*sys.argv[1:]))
