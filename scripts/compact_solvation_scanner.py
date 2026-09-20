"""Opt-in prepared-context OMOL + native GFN2 solvent score and timing pilot.

Reuses the native OMOL worker and existing ORCA executor. No new Hamiltonian,
geometry generation, reference fitting or baseline/default modification.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms
import compact_solvation as solvent
from compact_solvation_compare import mix_pair
import mace_omol as omol

SCHEMA = 'compact_prepared_context_scanner_v1'
CASES = ('1H4I', 'q9z4j7-pqq-la_model', '1F6S', '1GLG')
RESOURCES = {'cpus': 32, 'gpus': 1, 'host_mem_MiB': 200000,
             'mpi_ranks': 8, 'concurrent_GFN2_tasks': 4}
TOLERANCES = {'native_MACE_endpoint_kcal_mol': .01,
              'composite_repeat_model_kcal_mol': .02}
CONTEXT_PROTOCOLS = ('native_OMOL_complete_polar_context_disulfide_v2', 'native_OMOL_second_shell_fixed_v1')


def pqq_decision_compatible(source):
    return (source['label_scope'] in ('canonical_PQQ_functional_class', 'PQQ_prediction')
            and source['representations']['context']['protocol_id'] == CONTEXT_PROTOCOLS[0])


def low_level_prepare(source, case, agreement, directory, implementation, orca):
    """Only assemble existing recipe/runner inputs; no archived energy is required."""
    tasks = []
    for metal, e in source['representations']['context']['endpoints'].items():
        for medium in ('vacuum', 'alpb'):
            tid = case + '__context__' + metal + '__' + medium + '__native'
            td = directory / 'tasks' / tid; td.mkdir(parents=True)
            xp = td / 'core.xyz'; shutil.copyfile(verify(e['xyz']), xp)
            ip = td / 'endpoint.inp'; ip.write_text(solvent.input_text(e['charge'], e['multiplicity'], medium, 'native'))
            task = {'task_id': tid, 'case_id': case, 'case': case, 'representation': 'context', 'metal': metal,
                    'medium': medium, 'solver': 'native', 'charge': e['charge'], 'multiplicity': e['multiplicity'],
                    'xyz': record(xp), 'input': record(ip), 'output_path': str(td / 'endpoint.out')}
            task['scientific_key'] = cache_key({'xyz_sha256': e['xyz']['sha256'], 'charge': e['charge'],
                'multiplicity': e['multiplicity'], 'medium': medium, 'solver': 'native', 'method': solvent.METHOD,
                'input_body': ip.read_text(), 'orca': orca})
            tasks.append(task)
    m = {'protocol_id': solvent.PROTOCOL, 'method_id': solvent.METHOD, 'agreement': record(agreement),
         'orca': orca, 'implementation': implementation, 'tasks': tasks,
         'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 4},
         'execution_policy': {'task_runner': implementation['run_orca_task_manifest.py'],
                              'runtime_renderer': implementation['render_orca_runtime_input.py']}}
    path = directory / 'manifest.json'; write_new(path, m)
    return record(path)


def decision(value, bands):
    if value is None:
        return None
    if value <= bands['Ca_supported_max_R_model_kcal_mol']:
        return 'Ca-supported'
    if value >= bands['La_supported_min_R_model_kcal_mol']:
        return 'La-supported'
    return 'inconclusive'


def task_key(task, manifest):
    payload = {k: v for k, v in task.items() if k != 'cache_key'}
    return cache_key({'task': payload, 'model': manifest['model'],
                      'software': manifest['software'], 'implementation': manifest['implementation']})


def prepare(inventory, comparison, agreement, output, cpu_python, gpu_python, prepared_pairs=False):
    inv = read_json(inventory)
    old = read_json(comparison)
    if (not prepared_pairs and old['inventory'] != record(inventory)) or old['solver'] != 'native':
        raise InvalidArtifact('calibration does not describe the source inventory/native solver')
    bands = old['calibration']['context']['bands']
    if not bands or bands['protocol_id'] != solvent.PROTOCOL:
        raise InvalidArtifact('published context calibration unavailable')
    if inv['model'] != omol.model(verify(inv['software'])):
        raise InvalidArtifact('source checkpoint/model differs')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    impl = out / 'implementation'; impl.mkdir()
    implementation = {}
    for source_script in Path(__file__).parent.glob('*.py'):
        dest = impl / source_script.name; shutil.copyfile(source_script, dest)
        implementation[source_script.name] = record(dest)
    orca = read_json(Path(__file__).resolve().parents[1] / 'params/baseline_water_v1.json')['artifacts']['orca']
    sources = {c['case_id']: c for c in inv['cases']}
    cases = list(sources) if prepared_pairs else list(CASES)
    if not cases or len(cases) != len(set(cases)) or any(not re.fullmatch(r'[a-zA-Z0-9_.-]+', c) for c in cases):
        raise InvalidArtifact('empty/duplicate/invalid prepared case identity')
    low = {}
    for case in cases:
        directory = out / 'cases' / case / 'solvent'
        low[case] = low_level_prepare(sources[case], case, agreement, directory, implementation, orca)
    mm = {'schema_version': SCHEMA, 'protocol_id': solvent.PROTOCOL,
          'stage': 'compact_context_fresh_scalar_timing',
          'inventory': record(inventory), 'agreement': record(agreement),
          'software': inv['software'], 'model': inv['model'],
          'implementation': implementation, 'tasks': []}
    for case in cases:
        for metal in ('Ca', 'La'):
            source = sources[case]['representations']['context']['endpoints'][metal]
            task = {'task_id': case + '__context__' + metal, 'case_id': case,
                    'kind': 'core', 'representation': 'context', 'metal': metal,
                    'metal_index': source['metal_index'], 'charge': source['charge'],
                    'spin_multiplicity': source['multiplicity'], 'xyz': source['xyz'],
                    'energy_component': omol.COMPONENT, 'energy_only': True,
                    'capture_native_readout': True, 'source_endpoint': source}
            task['cache_key'] = task_key(task, mm)
            mm['tasks'].append(task)
    mace_path = out / 'mace_manifest.json'
    write_new(mace_path, mm)
    manifest = {'schema_version': SCHEMA, 'protocol_id': solvent.PROTOCOL,
                'inventory': record(inventory), 'archived_comparison': record(comparison),
                'agreement': record(agreement), 'published_context_bands': bands,
                'implementation': implementation, 'mace_manifest': record(mace_path),
                'solvent_manifests': low, 'cases': cases, 'prepared_pairs': prepared_pairs,
                'cpu_python': record(cpu_python), 'gpu_python': record(gpu_python),
                'cpu_python_invocation': str(Path(cpu_python).absolute()),
                'gpu_python_invocation': str(Path(gpu_python).absolute()),
                'execution_resources': RESOURCES, 'repeat_tolerances': TOLERANCES,
                'expected_new_calls': {'native_MACE': 2 * len(cases), 'native_GFN2': 4 * len(cases)},
                'cache_reuse': False, 'new_geometry_generated': False,
                'baseline_changed': False, 'absolute_aquo_reference': None,
                'timing_scope': 'prepared-input score path; excludes folding/context preparation'}
    write_new(out / 'manifest.json', manifest)
    return validate(out / 'manifest.json')


def prepare_pairs(pairs, comparison, agreement, output, cpu_python, gpu_python):
    """Explicit context pairs use the inventory source schema, without old energies/receipts."""
    return prepare(pairs, comparison, agreement, output, cpu_python, gpu_python, True)


def validate(manifest):
    m = read_json(manifest)
    if (m['schema_version'] != SCHEMA or m['protocol_id'] != solvent.PROTOCOL
            or (not m['prepared_pairs'] and m['cases'] != list(CASES)) or m['execution_resources'] != RESOURCES
            or m['repeat_tolerances'] != TOLERANCES or m['cache_reuse']):
        raise InvalidArtifact('frozen scanner scope/resources differ')
    inv = read_json(verify(m['inventory']))
    old = read_json(verify(m['archived_comparison']))
    if ((not m['prepared_pairs'] and old['inventory'] != m['inventory'])
            or old['calibration']['context']['bands'] != m['published_context_bands']):
        raise InvalidArtifact('published bands/source changed')
    verify(m['agreement'])
    verify(m['cpu_python']); verify(m['gpu_python'])
    for kind in ('cpu', 'gpu'):
        if record(m[kind + '_python_invocation']) != m[kind + '_python']:
            raise InvalidArtifact('interpreter invocation resolves to changed executable')
    for p in m['implementation'].values():
        verify(p)
    mm = read_json(verify(m['mace_manifest']))
    if (mm['model'] != inv['model'] or mm['software'] != inv['software']
            or mm['implementation'] != m['implementation']
            or mm['model'] != omol.model(verify(mm['software']))):
        raise InvalidArtifact('native checkpoint/software differs')
    if ({(t['case_id'], t['metal']) for t in mm['tasks']} != {(c, z) for c in m['cases'] for z in ('Ca', 'La')}
            or len(mm['tasks']) != 2 * len(m['cases'])):
        raise InvalidArtifact('native endpoint coverage differs')
    sources = {c['case_id']: c for c in inv['cases']}
    for case in m['cases']:
        rep = sources[case]['representations']['context']; prep = read_json(verify(rep['preparation']))
        if (rep['protocol_id'] not in CONTEXT_PROTOCOLS or prep['core_source_coordinates_unchanged'] is not True
                or prep['water_inventory_unchanged'] is not True or set(rep['endpoints']) != {'Ca', 'La'}):
            raise InvalidArtifact('unsupported context preparation/boundary')
        if sources[case]['label_scope'] in ('canonical_PQQ_functional_class', 'PQQ_prediction') and not pqq_decision_compatible(sources[case]):
            raise InvalidArtifact('PQQ band request requires the canonical PQQ complete-context preparation protocol')
        ca, la = (xyz(verify(rep['endpoints'][z]['xyz'])) for z in ('Ca', 'La'))
        if (len(ca) != len(la) or ca[0][0] != 'Ca' or la[0][0] != 'La'
                or any(a[0] != b[0] for a, b in zip(ca[1:], la[1:]))
                or rep['endpoints']['La']['charge'] - rep['endpoints']['Ca']['charge'] != 1):
            raise InvalidArtifact('paired source composition/charge differs')
    for t in mm['tasks']:
        s = sources[t['case_id']]['representations']['context']['endpoints'][t['metal']]
        if (t['source_endpoint'] != s or t['xyz'] != s['xyz'] or t['charge'] != s['charge']
                or t['spin_multiplicity'] != s['multiplicity'] or t['metal_index'] != s['metal_index']
                or t['energy_component'] != omol.COMPONENT or t.get('energy_only') is not True
                or t.get('capture_native_readout') is not True or t['cache_key'] != task_key(t, mm)):
            raise InvalidArtifact('native state/coordinates/cache differs')
        check_atoms(xyz(verify(t['xyz'])), t['charge'])
        if t['spin_multiplicity'] != 1:
            raise InvalidArtifact('unsupported multiplicity')
        if s.get('native_MACE_receipt'):
            actual = read_json(verify(s['native_MACE_receipt']))
            previous = read_json(verify(s['source_manifest']))
            if actual['energy_eV'] != s['native_MACE_energy_eV'] or previous['model'] != mm['model']:
                raise InvalidArtifact('archived scalar/checkpoint differs')
    for case, pin in m['solvent_manifests'].items():
        path = verify(pin); sm = read_json(path)
        from affordable_workflow import dry_run
        dry_run(path)
        if (sm['method_id'] != solvent.METHOD or sm['execution_resources'] != {'mpi_ranks': 8, 'concurrent_tasks': 4}
                or len(sm['tasks']) != 4 or {(t['metal'], t['medium']) for t in sm['tasks']} !=
                    {(z, s) for z in ('Ca', 'La') for s in ('vacuum', 'alpb')}):
            raise InvalidArtifact('solvent task scope/reuse differs')
        for t in sm['tasks']:
            e = sources[case]['representations']['context']['endpoints'][t['metal']]
            body = solvent.input_text(e['charge'], e['multiplicity'], t['medium'], 'native')
            if (t['case_id'] != case or t['solver'] != 'native' or t['charge'] != e['charge']
                    or t['multiplicity'] != e['multiplicity'] or verify(t['input']).read_text() != body
                    or verify(t['xyz']).read_bytes() != verify(e['xyz']).read_bytes()
                    or t['scientific_key'] != cache_key({'xyz_sha256': e['xyz']['sha256'], 'charge': e['charge'],
                        'multiplicity': e['multiplicity'], 'medium': t['medium'], 'solver': 'native', 'method': solvent.METHOD,
                        'input_body': body, 'orca': sm['orca']})):
                raise InvalidArtifact('GFN2 prepared recipe/state differs')
    return {'status': 'validated', 'manifest': record(manifest), **m['expected_new_calls']}


def mace_stage(manifest, case):
    m = read_json(manifest); mp = verify(m['mace_manifest']); mm = read_json(mp)
    if case not in m['cases']:
        raise InvalidArtifact('unknown case')
    for task in mm['tasks']:
        if task['case_id'] != case:
            continue
        out = Path(manifest).parent / 'mace' / task['task_id']
        out.mkdir(parents=True, exist_ok=False)
        omol.worker(mp, task['task_id'], out, 'native')


def fresh_mace(manifest, task):
    m = read_json(manifest)
    path = Path(manifest).parent / 'mace' / task['task_id'] / 'result.json'
    if not path.exists():
        return {'status': 'unavailable', 'energy_eV': None, 'reason': 'no fresh result'}
    result = read_json(path)
    base = {'receipt': record(path), 'energy_eV': None, 'status': 'unavailable'}
    if (result.get('status') != 'computed' or result.get('manifest') != m['mace_manifest']
            or result.get('cache_key') != task['cache_key'] or not omol.accepted_state(result, task)):
        return {**base, 'reason': result.get('reason', 'native scalar/state acceptance failed')}
    old = task['source_endpoint'].get('native_MACE_energy_eV')
    delta = (result['energy_eV'] - old) * EV_TO_KCAL if old is not None else None
    return {**base, 'status': 'available', 'energy_eV': result['energy_eV'],
            'archived_energy_eV': old,
            'archived_difference_kcal_mol': delta,
            'repeat_pass': abs(delta) <= TOLERANCES['native_MACE_endpoint_kcal_mol'] if delta is not None else None,
            'analytic_gradient_status': result['analytic_gradient_status'],
            'model_load_seconds': result['model_load_seconds'],
            'evaluation_seconds': result['evaluation_seconds'],
            'wall_seconds': result['wall_seconds'], 'device': result['device'],
            'peak_host_RSS_KiB': result['peak_host_RSS_KiB'],
            'peak_cuda_allocated_bytes': result['peak_cuda_allocated_bytes']}


def low_level_collect(manifest, output):
    m = read_json(manifest); eps = {}; rows = []
    for task in m['tasks']:
        pin = solvent.completed(manifest, task['task_id'])
        result = {'status': 'unavailable', 'energy_hartree': None, 'task_id': task['task_id']}
        if pin:
            try:
                result = {'status': 'complete', **pin, **solvent.diagnostics(pin, task)}
            except InvalidArtifact as exc:
                result.update(reason=str(exc), attempted_receipt=pin)
        else:
            paths = [Path(task['output_path']), Path(task['output_path'] + '.execution.json')]
            result['available_artifacts'] = [record(p) for p in paths if p.exists()]
        eps[(task['metal'], task['medium'])] = result
    for metal in ('Ca', 'La'):
        pair = {s: eps[(metal, s)] for s in ('vacuum', 'alpb')}
        good = all(e['status'] == 'complete' for e in pair.values())
        vacuum, alpb = (pair[s]['energy_hartree'] for s in ('vacuum', 'alpb'))
        rows.append({'metal': metal, 'status': 'complete' if good else 'unavailable',
                     'endpoints': pair, 'vacuum_hartree': vacuum, 'alpb_hartree': alpb,
                     'delta_solv_hartree': alpb - vacuum if good else None})
    result = {'manifest': record(manifest), 'rows': rows}
    write_new(output, result)
    return result


def solvent_stage(manifest):
    from affordable_workflow import execute as existing_execute
    return existing_execute(manifest)


def collect(manifest, output):
    validate(manifest)
    m = read_json(manifest); out = Path(manifest).parent
    inv = read_json(verify(m['inventory'])); mm = read_json(verify(m['mace_manifest']))
    archived = read_json(verify(m['archived_comparison']))
    rows = []
    for case in m['cases']:
        source = next(c for c in inv['cases'] if c['case_id'] == case)
        original = (next(r for r in archived['rows'] if r['case_id'] == case and r['representation'] == 'context')
                    if not m['prepared_pairs'] else {})
        native = {t['metal']: fresh_mace(manifest, t) for t in mm['tasks'] if t['case_id'] == case}
        smp = verify(m['solvent_manifests'][case])
        # A fresh output destination makes manual partial/final collections immutable.
        cp = Path(output).parent / (Path(output).stem + '__' + case + '__GFN2.json')
        low = low_level_collect(smp, cp)
        eps = {r['metal']: r for r in low['rows']}
        row = {'case_id': case, 'representation': 'context', 'role': source['role'],
               'biological_group': source['biological_group'], 'label_scope': source['label_scope'],
               'expected_class': source['expected_class'], 'native_MACE': native,
               'low_level_collection': record(cp), 'low_level_endpoints': eps,
               'status': 'unavailable', 'composite_R_model_kcal_mol': None,
               'published_PQQ_decision': None, 'archived_composite_decision': original.get('composite_decision'),
               'archived_composite_R_model_kcal_mol': original.get('composite_R_model_kcal_mol'),
               'repeat_difference_model_kcal_mol': None, 'repeat_pass': None}
        if all(v['status'] == 'available' for v in native.values()) and all(v['status'] == 'complete' for v in eps.values()):
            score = mix_pair({z: native[z]['energy_eV'] for z in native},
                             {z: eps[z]['vacuum_hartree'] for z in eps},
                             {z: eps[z]['alpb_hartree'] for z in eps})
            delta = score['composite_R_model_kcal_mol'] - original['composite_R_model_kcal_mol'] if original else None
            row.update(score, status='available', repeat_difference_model_kcal_mol=delta,
                       repeat_pass=abs(delta) <= TOLERANCES['composite_repeat_model_kcal_mol'] if delta is not None else None,
                       low_level_repeat_differences_hartree={z: {medium: eps[z][medium + '_hartree'] - original['low_level_endpoints'][z][medium + '_hartree']
                           for medium in ('vacuum', 'alpb')} for z in eps} if original else None)
            if pqq_decision_compatible(source):
                row['published_PQQ_decision'] = decision(score['composite_R_model_kcal_mol'], m['published_context_bands'])
        timing = out / 'timing' / (case + '.json')
        row['timing'] = read_json(timing) if timing.exists() else None
        rows.append(row)
    index = {r['case_id']: r for r in rows}
    margin = (index['1F6S']['composite_R_model_kcal_mol'] - index['1GLG']['composite_R_model_kcal_mol']
              if all(c in index and index[c]['status'] == 'available' for c in ('1F6S', '1GLG')) else None)
    result = {'schema_version': SCHEMA, 'protocol_id': solvent.PROTOCOL,
              'manifest': record(manifest), 'rows': rows, 'case_denominator': len(m['cases']),
              'available': sum(r['status'] == 'available' for r in rows),
              'published_context_bands': m['published_context_bands'],
              'repeat_tolerances': TOLERANCES,
              'alpha_minus_GGR_model_kcal_mol': margin,
              'alpha_GGR_independent_biological_comparisons': 1 if margin is not None else 0,
              'numerical_qualification': 'see separate numerical qualification; this tests execution repeatability',
              'all_evidence_consumed': all(next(c for c in inv['cases'] if c['case_id'] == cid).get('all_evidence_consumed', False)
                                           for cid in m['cases']), 'absolute_aquo_reference': None,
              'timing_scope': m['timing_scope'], 'new_geometry_generated': False, 'baseline_changed': False}
    write_new(output, result)
    return result


def execute(manifest):
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('scanner execution requires an allocated GPU node')
    overall = time.monotonic(); m = read_json(manifest); out = Path(manifest).parent
    if (int(os.environ['SLURM_CPUS_ON_NODE']) != RESOURCES['cpus']
            or int(os.environ['SLURM_NTASKS']) != RESOURCES['cpus']
            or int(os.environ['SLURM_MEM_PER_NODE']) != RESOURCES['host_mem_MiB']):
        raise InvalidArtifact('allocation differs from frozen CPU/RAM request')
    write_new(out / 'execution_started.json', {'manifest': record(manifest),
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'time_unix': time.time()})
    validate(manifest)
    (out / 'timing').mkdir(exist_ok=False)
    script = verify(m['implementation']['compact_solvation_scanner.py'])
    for case in m['cases']:
        start = time.monotonic(); stages = {}
        commands = [('MACE', ['srun', '--exact', '--ntasks=1', '--cpus-per-task=32', '--gres=gpu:1',
                              m['gpu_python_invocation'], str(script), 'mace-stage', '--manifest', str(manifest), '--case', case])]
        smp = verify(m['solvent_manifests'][case]); sm = read_json(smp)
        commands.append(('GFN2', [m['cpu_python_invocation'], str(script),
                                 'solvent-stage', '--manifest', str(smp)]))
        for name, command in commands:
            began = time.monotonic()
            with (out / 'timing' / (case + '__' + name + '.log')).open('x') as log:
                run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
            stages[name] = {'wall_seconds': time.monotonic() - began, 'returncode': run.returncode,
                            'command': command, 'log': record(out / 'timing' / (case + '__' + name + '.log'))}
        elapsed = time.monotonic() - start
        receipt = {'case_id': case, 'stages': stages, 'fresh_endpoint_path_seconds': elapsed,
                   'allocated_core_seconds': elapsed * RESOURCES['cpus'], 'allocated_GPU_seconds': elapsed,
                   'slurm_job_id': os.environ['SLURM_JOB_ID'], 'allocated_resources': RESOURCES,
                   'score_collection_included': False}
        write_new(out / 'timing' / (case + '.json'), receipt)
        print(json.dumps({'case': case, 'fresh_endpoint_path_seconds': elapsed}), flush=True)
    result = collect(manifest, out / ('result_' + os.environ['SLURM_JOB_ID'] + '.json'))
    elapsed = time.monotonic() - overall
    write_new(out / 'execution_complete.json', {'manifest': record(manifest), 'available': result['available'],
              'prepared_input_score_path_seconds': elapsed, 'allocated_core_seconds': elapsed * RESOURCES['cpus'],
              'allocated_GPU_seconds': elapsed, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'collection_included': True, 'folding_or_context_preparation_included': False})
    return {'status': 'complete' if result['available'] == len(m['cases']) else 'partial', 'available': result['available'], 'wall_seconds': elapsed}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('prepare')
    for key in ('inventory', 'comparison', 'agreement', 'output', 'cpu-python', 'gpu-python'):
        a.add_argument('--' + key, required=True)
    a = sub.add_parser('prepare-pairs')
    for key in ('pairs', 'comparison', 'agreement', 'output', 'cpu-python', 'gpu-python'):
        a.add_argument('--' + key, required=True)
    for name in ('validate', 'execute', 'mace-stage', 'solvent-stage', 'collect'):
        a = sub.add_parser(name); a.add_argument('--manifest', required=True)
        if name == 'collect': a.add_argument('--output', required=True)
        if name == 'mace-stage': a.add_argument('--case', required=True)
    args = vars(p.parse_args()); cmd = args.pop('cmd').replace('-', '_')
    result = globals()[cmd](**args)
    print(json.dumps(result if cmd != 'collect' else {'available': result['available']}, sort_keys=True))


if __name__ == '__main__':
    main()
