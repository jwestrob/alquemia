"""Standard prepared-source dispatch: fast compatible PQQ or the existing baseline.

The explicit DFT reference path remains available. Unsupported fast preparation
fails without replacing its score with DFT or an older cache entry.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from affordable_common import InvalidArtifact, HA_TO_KCAL, classify_raw, read_json, record, verify, write_new

WORKFLOW = 'standard_discriminator_v1'
FAST_MODE = 'fast_PQQ_OMOL_GFN2_ALPB_v1'
SOURCE_PROTOCOL = 'fixed_core_PQQ_source_to_complete_context_v1'
CORE_PROTOCOL = 'pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
DEFAULT_RELEASE = Path(__file__).resolve().parents[1] / 'params/pqq_fast_v1.json'


def release(path):
    r = read_json(path)
    if r['workflow_id'] != WORKFLOW or r['fast_mode'] != FAST_MODE:
        raise InvalidArtifact('unsupported standard-mode release')
    for pin in r['artifacts'].values(): verify(pin)
    for name, pin in r['implementation'].items():
        verify(pin)
        if record(Path(__file__).parent / name)['sha256'] != pin['sha256']:
            raise InvalidArtifact('released implementation changed: ' + name)
    result = read_json(verify(r['artifacts']['release_result']))
    if (r['promoted'] is not True or result['promotion_gate'] is not True
            or result['calibration_correct'] != 25 or result['transfer_correct'] != 3
            or result['preparation_supported'] != 28):
        raise InvalidArtifact('fast PQQ release checks have not passed')
    return r


def choose_backend(request, mode='auto'):
    if mode not in ('auto', 'fast-pqq', 'dft-reference'):
        raise InvalidArtifact('unsupported standard mode')
    if request.get('workflow_id') == 'baseline_contextual_water_v1':
        if mode == 'fast-pqq':
            raise InvalidArtifact('fast mode requires the supported explicit PQQ source preparation')
        return 'baseline_water'
    if request.get('schema_version') != SOURCE_PROTOCOL:
        raise InvalidArtifact('unsupported input; provide a fixed-core PQQ source request or an existing baseline-water request')
    return 'dft_reference' if mode == 'dft-reference' else 'fast_pqq'


def prepare(request, output, mode='auto', release_path=DEFAULT_RELEASE):
    req = read_json(request); backend = choose_backend(req, mode)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    if backend == 'baseline_water':
        from baseline_water import prepare as existing_prepare
        result = existing_prepare(request, out / 'baseline')
        plan = {'workflow_id': WORKFLOW, 'backend': backend, 'mode': mode,
                'request': record(request), 'baseline_plan': result['plan'],
                'baseline_behavior_changed': False, 'next_stage': result['next_stage']}
    else:
        r = release(release_path)
        if req['config'] != r['source_configuration']:
            raise InvalidArtifact('source preparation/model configuration differs from the released method')
        ids = [c['case_id'] for c in req['cases']]
        if not ids or len(set(ids)) != len(ids): raise InvalidArtifact('empty/duplicate source case list')
        for case in req['cases']:
            verify(case['source_structure'])
        impl = out / 'implementation'; impl.mkdir(); pins = {}
        for p in Path(__file__).parent.glob('*.py'):
            dest = impl / p.name; shutil.copyfile(p, dest); pins[p.name] = record(dest)
        plan = {'workflow_id': WORKFLOW, 'backend': backend, 'mode': mode, 'request': record(request),
                'release': record(release_path), 'implementation': pins,
                'new_source_preparations': len(ids), 'new_MACE_calls': 2 * len(ids) if backend == 'fast_pqq' else 0,
                'new_GFN2_calls': 4 * len(ids) if backend == 'fast_pqq' else 0,
                'new_DFT_calls': 2 * len(ids) if backend == 'dft_reference' else 0,
                'source_domain_status': {c['case_id']: 'consumed_reference_geometry' if c['source_structure']['sha256'] in r['reference_source_sha256']
                                         else 'unvalidated_input_domain' for c in req['cases']},
                'fallback_policy': 'explicit_failure_no_score_substitution', 'next_stage': 'execute'}
    path = out / 'plan.json'; write_new(path, plan)
    return {'status': 'prepared', 'backend': backend, 'plan': record(path), 'next_stage': plan['next_stage']}


def checked_plan(plan):
    p = read_json(plan); req = read_json(verify(p['request']))
    if p['workflow_id'] != WORKFLOW or p['backend'] != choose_backend(req, p['mode']):
        raise InvalidArtifact('standard dispatch changed')
    if p['backend'] == 'baseline_water':
        from baseline_water import checked_plan as existing_checked
        existing_checked(verify(p['baseline_plan']))
    else:
        r = release(verify(p['release']))
        if req['config'] != r['source_configuration']: raise InvalidArtifact('source configuration changed')
        for pin in p['implementation'].values(): verify(pin)
        for c in req['cases']: verify(c['source_structure'])
    return p


def prepare_dft(preparation, output, release_path):
    """Exact original-core native reference recipe; no DFT execution here."""
    p = read_json(preparation); r = release(release_path)
    if p['supported'] != p['denominator']: raise InvalidArtifact('source preparation incomplete')
    baseline = read_json(verify(r['artifacts']['baseline_release']))
    refs = baseline['artifacts']; out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    tasks = []
    for c in p['cases']:
        for z in ('Ca', 'La'):
            e = c['core']['endpoints'][z]; td = out / 'tasks' / (c['case_id'] + '__' + z); td.mkdir(parents=True)
            ip = td / Path(e['input']['path']).name; xp = td / Path(e['xyz']['path']).name
            shutil.copyfile(verify(e['input']), ip); shutil.copyfile(verify(e['xyz']), xp)
            tasks.append({'task_id': td.name, 'case': c['case_id'], 'metal': z, 'charge': e['charge'],
                          'multiplicity': 1, 'xyz': record(xp), 'input': record(ip), 'output_path': str(ip.with_suffix('.out'))})
    manifest = {'protocol_id': CORE_PROTOCOL, 'source_preparation': record(preparation), 'release': record(release_path),
                'agreement': r['artifacts']['agreement'], 'orca': refs['orca'], 'tasks': tasks,
                'execution_policy': {'task_runner': refs['task_runner'], 'runtime_renderer': refs['runtime_renderer']},
                'execution_resources': {'mpi_ranks': 16, 'concurrent_tasks': 2},
                'reference': refs['canonical_release'], 'historical_records_preserved': True}
    path = out / 'manifest.json'; write_new(path, manifest)
    from affordable_workflow import dry_run
    return dry_run(path)


def execute(plan):
    if not os.environ.get('SLURM_JOB_ID'): raise InvalidArtifact('source-to-score execution requires an allocation')
    start = time.monotonic(); p = checked_plan(plan); out = Path(plan).parent
    if p['backend'] == 'baseline_water':
        # This retains the existing staged water/GPU/DFT resource policy unchanged.
        return {'status': 'existing_baseline_pipeline', 'plan': p['baseline_plan'],
                'next_stage': p['next_stage'], 'entrypoint': 'affordable_workflow.py baseline',
                'automatic_scientific_fallback': False}
    r = release(verify(p['release']))
    from pqq_fast_prepare import prepare as source_prepare
    source_prepare(verify(p['request']), out / 'source_preparation')
    preparation = out / 'source_preparation/preparation.json'
    if p['backend'] == 'fast_pqq':
        from pqq_fast_release import prepare_score
        prepare_score(preparation, verify(r['artifacts']['calibration']), verify(r['artifacts']['agreement']),
                      out / 'prepared_score', r['cpu_python'], r['gpu_python'])
        manifest = out / 'prepared_score/scoring/manifest.json'
        script = verify(read_json(manifest)['implementation']['compact_solvation_scanner.py'])
        subprocess.run([r['cpu_python'], str(script), 'execute', '--manifest', str(manifest)], check=True)
    else:
        prepare_dft(preparation, out / 'dft', verify(p['release']))
        manifest = out / 'dft/manifest.json'
        from affordable_workflow import execute as existing_execute
        existing_execute(manifest)
    elapsed = time.monotonic() - start
    receipt = {'plan': record(plan), 'backend': p['backend'], 'manifest': record(manifest),
               'source_preparation': record(preparation), 'source_to_score_seconds': elapsed,
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'allocated_CPUs': int(os.environ['SLURM_CPUS_ON_NODE']),
               'source_preparation_included': True, 'folding_included': False}
    write_new(out / 'execution.json', receipt)
    return {'status': 'executed', 'receipt': record(out / 'execution.json'), 'next_stage': 'collect'}


def collect(plan, output):
    p = checked_plan(plan); out = Path(plan).parent
    if p['backend'] == 'baseline_water':
        raise InvalidArtifact('collect the existing baseline DFT manifest with affordable_workflow.py baseline collect')
    r = release(verify(p['release']))
    if p['backend'] == 'fast_pqq':
        from compact_solvation_scanner import collect as existing_collect
        component_path = Path(output).with_name(Path(output).stem + '__components.json')
        result = existing_collect(out / 'prepared_score/scoring/manifest.json', component_path)
        result.update(standard_plan=record(plan), component_collection=record(component_path), fast_mode=FAST_MODE,
                      source_domain_status=p['source_domain_status'])
        for row in result['rows']:
            row['source_domain_status'] = p['source_domain_status'][row['case_id']]
        write_new(output, result)
        return {'status': 'complete' if result['available'] == result['case_denominator'] else 'partial',
                'collection': record(output), 'source_domain_status': p['source_domain_status']}
    from baseline_water import completed
    mp = out / 'dft/manifest.json'; m = read_json(mp); ref = read_json(verify(m['reference'])); rows = []
    for case in read_json(verify(p['request']))['cases']:
        cid = case['case_id']; eps = {z: completed(mp, cid + '__' + z) for z in ('Ca', 'La')}
        raw = (eps['Ca']['energy_hartree'] - eps['La']['energy_hartree']) * HA_TO_KCAL if all(eps.values()) else None
        rows.append({'case_id': cid, 'protocol_id': CORE_PROTOCOL, 'endpoints': eps,
                     'R_kcal_mol': raw, 'S_kcal_mol': raw - ref['aquo_reporting_gauge']['A_kcal_mol'] if raw is not None else None,
                     'decision': classify_raw(raw, ref, CORE_PROTOCOL) if raw is not None else 'unavailable',
                     'source_domain_status': p['source_domain_status'][cid]})
    result = {'workflow_id': WORKFLOW, 'backend': p['backend'], 'plan': record(plan), 'manifest': record(mp),
              'rows': rows, 'status': 'complete' if all(x['R_kcal_mol'] is not None for x in rows) else 'partial'}
    write_new(output, result); return {'status': result['status'], 'collection': record(output)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare')
    for name in ('request', 'output'): q.add_argument('--' + name, required=True)
    q.add_argument('--mode', choices=('auto', 'fast-pqq', 'dft-reference'), default='auto')
    q.add_argument('--release-path', default=str(DEFAULT_RELEASE))
    for name in ('validate', 'execute', 'collect'):
        q = sub.add_parser(name); q.add_argument('--plan', required=True)
        if name == 'collect': q.add_argument('--output', required=True)
    a = vars(p.parse_args(argv)); op = a.pop('op')
    result = checked_plan(**a) if op == 'validate' else globals()[op](**a)
    print(json.dumps(result if op != 'validate' else {'status': 'validated', 'backend': result['backend']}))


if __name__ == '__main__': main()
