#!/usr/bin/env python3
"""Four embedded GFN2 damping tests; no evaluator or response continuation."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

PROTOCOL = 'lbt_response_embedded_scc_four_damping_controls_v1'
CASES = (('failed_005', 'failed', .05), ('failed_010', 'failed', .10),
         ('preceding_005', 'preceding', .05), ('origin_005', 'origin', .05))
GATES = {'energy_hartree': 1e-7, 'gradient_component_kcal_mol_A': .005}
ARGS = ['core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0', '--etemp', '300',
        '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--grad',
        '--input', 'control.inp']
EMBEDDING = '$embedding\n interface=orca\n input=external.pc\n gradient=pcgrad\n$end\n'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def verify_tree(value):
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            verify(value)
        else:
            for v in value.values():
                verify_tree(v)
    elif isinstance(value, list):
        for v in value:
            verify_tree(v)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def controls_text(damping):
    return EMBEDDING + '$scc\n  broydamp=' + format(damping, '.2f') + '\n$end\n'


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['agreement_gates'] == GATES
            and cfg['primary_case'] == 'failed_005'
            and cfg['cases'] == [{'case_id': c, 'fixture': f, 'broydamp': d} for c, f, d in CASES],
            'Frozen diagnostic declaration differs')
    verify_tree(cfg)
    old = read(cfg['original_failed_response'])
    require(old['slurm_job_id'] == '13019' and old['status'] == 'unavailable'
            and (old['source_id'], old['region'], old['metal']) == ('7CCN', 'baseline', 'La')
            and old['termination'] == 'execution_failure'
            and old['continuation_execution_and_physical_admission'] is False
            and old['accepted_steps'] == 2 and old['distinct_evaluations'] == 3,
            'Original failed response must be retained')
    oldcfg, search = read(old['config']), read(old['search'])
    verify(old['implementation']); verify(oldcfg['implementation'])
    require(old['implementation']['sha256'] == oldcfg['implementation']['sha256']
            == '05b1c412feebfc72048405c8723d426168872325ce32d956732f5c898c028981'
            and search['status'] == 'unavailable' and search['termination'] == 'execution_failure'
            and search['accepted_steps'] == 2 and search['distinct_evaluations'] == 3
            and len(old['requests']) == 3 and old['requests'][-1] == cfg['fixtures']['failed']['point']
            and old['requests'][-2] == old['terminal_point']['point'] == cfg['fixtures']['preceding']['point'],
            'Failed-trial/preceding endpoint lineage differs')
    origins = read(cfg['original_origins_result'])
    verify(origins['implementation']); read(origins['config'])
    require(origins['status'] == 'complete' and origins['stage'] == 'origins'
            and origins['slurm_job_id'] == '12987'
            and len(origins['rows']) == 8 and all(r['status'] == 'complete' for r in origins['rows'])
            and [r['point'] for r in origins['rows'] if (r['source_id'], r['region'], r['metal'])
                 == ('7CCN', 'baseline', 'La')] == [cfg['fixtures']['origin']['point']],
            'Original accepted eight-origin control lineage differs')
    installed = read(cfg['installed_package'])
    require(installed['name'] == 'xtb' and installed['version'] == '6.7.1' and installed['build'] == 'h8876d29_4',
            'Installed native package lineage differs')
    paths = {r['_path']: r for r in installed['paths_data']['paths']}
    require(paths['bin/xtb']['sha256_in_prefix'] == cfg['runtime']['executable']['sha256']
            and paths['share/man/man7/xcontrol.7']['sha256_in_prefix'] == cfg['installed_xcontrol_manual']['sha256'],
            'Binary/manual do not match installed package')
    require(cfg['runtime']['executable']['sha256'] == '481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4'
            and cfg['runtime']['parameters']['sha256'] == 'f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545'
            and cfg['corrected_profile']['sha256'] == '635e111bac31e7c5755bf612c8e148ce73b401c110ff31887fd5e31b590fa85f',
            'Qualified binary/full-parameter lineage differs')
    manual = verify(cfg['installed_xcontrol_manual']).read_text()
    parser = verify(cfg['source_evidence']['set_module']).read_text()
    require('broydamp' in manual and 'Broyden convergence accelerator' in manual
            and "case('broydamp')" in parser and 'set%broydamp = ddum' in parser,
            'Installed documentation/source do not support requested parser key')
    # The installed manual is package-linked; the saved release source documents
    # the legacy algorithm. Actual output must separately prove parser uptake.
    backend = module('scc_diagnostic_backend', cfg['backend'])
    sys.path.insert(0, str(verify(cfg['common_helper']).parent))
    gradients = module('scc_diagnostic_gradients', cfg['gradient_helper'])
    admission = module('scc_diagnostic_admission', cfg['admission'])
    mapping = module('scc_diagnostic_mapping', cfg['mapping_helper'])
    import numpy as np
    mapping.np = np
    profile = read(cfg['corrected_profile'])
    require(all(cfg['runtime']['environment'][k] == '1' for k in
                ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')),
            'Native threads must remain one')
    fixtures = {}
    for name in ('failed', 'preceding', 'origin'):
        pins = cfg['fixtures'][name]
        point = read(pins['point'])
        native = read(pins['native_result'])
        verify_tree(native)
        task = read(native['task'])
        execution = read(native['execution'])
        require(point['component_receipts']['GFN_embedded'] == pins['native_result']
                and point['identity_sha256'] == old['identity_sha256'] == origins['identity_sha256']
                and (task['source_id'], task['region'], task['metal'], task['kind']) == ('7CCN', 'baseline', 'La', 'embedded')
                and task['atom_count'] == 172 and task['charge'] == -2 and task['uhf'] == 0
                and task['valence_electrons'] == 538 and task['multiplicity'] == task['physical_multiplicity'] == 1
                and task['pointcharge_count'] == 84 and task['field_charge_e'] == point['state']['field_charge_e']
                and task['runtime'] == cfg['runtime'] and task['full_corrected_profile'] == cfg['corrected_profile']
                and task['command'] == execution['command'] == [str(verify(cfg['runtime']['executable']))] + ARGS,
                'Archived embedded fixture state/command differs: ' + name)
        verify_tree(task)
        mp = read(point['mapping'])
        positions = read(point['positions'])
        sy, x = admission.xyz(task['xyz'])
        require(sy == ['La'] + [a['element'] for a in mp['core'][1:]]
                and np.array_equal(x, [a['xyz_A'] for a in mp['core']])
                and positions['source_ids'] == [a['id'] for a in mp['physical_atoms']]
                and np.array_equal(positions['positions_A'], [a['xyz_A'] for a in mp['physical_atoms']]),
                'Archived exact core/physical map differs')
        field = read(task['field'])
        require(field['rows'] == mp['field'] and len(mp['field']) == 84
                and field['parameter_source'] == cfg['runtime']['parameters']
                and field['potential_reference'] == 'finite hardness-damped potential zero at infinity; no PC self-energy'
                and verify(task['staged_files']['control.inp']).read_text() == EMBEDDING,
                'Archived embedding input/field map differs')
        expected_pc = '\n'.join(['84'] + [' '.join(format(v, '.17g') for v in
            (a['charge_e'], *a['xyz_A'], profile['element'][a['element']]['gam'])) for a in mp['field']]) + '\n'
        require(verify(task['staged_files']['external.pc']).read_text() == expected_pc
                and task['staged_files']['core.xyz'] == task['xyz']
                and task['staged_files']['parameters']['sha256'] == cfg['runtime']['parameters']['sha256'],
                'Archived PC coordinates/charges/hardness or input bytes differ')
        fixture = {'pins': pins, 'point': point, 'native': native, 'task': task, 'mapping': mp}
        directory = verify(native['task']).parent
        text = verify(native['artifacts']['stdout.txt']).read_text()
        values = re.findall(r'Broyden damping\s+([-+\d.Ee]+)', text)
        require(values and float(values[-1]) == .4, 'Archived control is not the accepted default damping')
        if name == 'failed':
            require(point['status'] == 'unavailable' and native['status'] == 'unavailable'
                    and execution['returncode'] == 1 and execution['owned_process_group_empty']
                    and 'convergence criteria cannot be satisfied within 500 iterations' in text,
                    'Expected failed SCC evidence absent')
            isolated = read(point['component_receipts']['GFN_isolated'])
            verify_tree(isolated)
            itask = read(isolated['task'])
            require(isolated['status'] == 'accepted' and itask['kind'] == 'isolated'
                    and itask['xyz']['sha256'] == task['xyz']['sha256']
                    and itask['runtime'] == task['runtime'], 'Accepted isolated companion differs')
        else:
            require(point['status'] == 'complete' and native['status'] == 'accepted', 'Passed control unavailable')
            replay = admission.admit(directory, task, execution, backend=backend, gradients=gradients, np=np)
            require(all(replay[k] == native[k] for k in ('energy_hartree', 'charges_e', 'scc_iterations',
                                                       'qm_gradient', 'pointcharge_gradient')),
                    'Saved control native replay differs')
            fixture['reference'] = replay
        fixtures[name] = fixture
    require(len({f['task']['xyz']['sha256'] for f in fixtures.values()}) == 3,
            'Controls must be distinct actual archived geometries')
    result['preflight'] = write(out / 'PREFLIGHT.json', {'status': 'passed', 'new_molecular_calls': 0,
        'old_failure_preserved': cfg['original_failed_response'], 'original_origins_result': cfg['original_origins_result'],
        'primary_case': cfg['primary_case'], 'installed_package': cfg['installed_package'],
        'manual_binary_package_link_verified': True, 'actual_printed_damping_required_after_call': True,
        'fixed_args': ARGS, 'cases': cfg['cases'], 'fixtures': cfg['fixtures'], 'agreement_gates': GATES,
        'controls': {c: controls_text(d) for c, _, d in CASES}, 'native_admission': cfg['admission']})
    return backend, gradients, admission, mapping, np, fixtures


def run_cell(cfg, out, row, helpers):
    backend, gradients, admission, _, np, fixtures = helpers
    fixture = fixtures[row['fixture']]
    d = out / row['case_id']
    d.mkdir()
    task = copy.deepcopy(fixture['task'])
    (d / 'core.xyz').write_bytes(verify(task['xyz']).read_bytes())
    (d / 'param_gfn2-xtb.txt').write_bytes(verify(cfg['runtime']['parameters']).read_bytes())
    (d / 'external.pc').write_bytes(verify(task['staged_files']['external.pc']).read_bytes())
    (d / 'control.inp').write_text(controls_text(row['broydamp']))
    staged = {'core.xyz': pin(d / 'core.xyz'), 'parameters': pin(d / 'param_gfn2-xtb.txt'),
              'external.pc': pin(d / 'external.pc'), 'control.inp': pin(d / 'control.inp')}
    require(staged['core.xyz']['sha256'] == task['xyz']['sha256']
            and staged['external.pc']['sha256'] == task['staged_files']['external.pc']['sha256']
            and staged['parameters']['sha256'] == cfg['runtime']['parameters']['sha256'], 'Staged scientific bytes differ')
    command = [str(verify(cfg['runtime']['executable']))] + ARGS
    task.update(xyz=staged['core.xyz'], command=command, staged_files=staged,
                original_task=fixture['native']['task'], solver_diagnostic_protocol=PROTOCOL,
                actual_external_pointcharges=staged['external.pc'],
                solver_control={'broydamp': row['broydamp'], 'control': staged['control.inp'],
                                'unchanged_embedding_control': EMBEDDING,
                                'native_adaptive_behavior_unchanged': True})
    row['task'] = write(d / 'TASK.json', task)
    env = backend.clean_env(staged['parameters']['path'])
    for k, v in cfg['runtime']['environment'].items():
        if v is None:
            env.pop(k, None)
        else:
            env[k] = v
    env['XTBPATH'] = str(d.resolve())
    verify_tree(cfg['runtime']); verify_tree(staged); verify(row['task'])
    row.update(status='unavailable_attempted', attempted=True)
    try:
        execution = backend.call(command, d, env)
        row['execution'] = pin(d / 'EXECUTION.json')
        verify_tree(staged)
        text = (d / 'stdout.txt').read_text()
        values = re.findall(r'Broyden damping\s+([-+\d.Ee]+)', text)
        require(values and float(values[-1]) == row['broydamp'], 'Native parser did not report requested damping')
        row['printed_initial_broydamp'] = float(values[-1])
        admitted = admission.admit(d, task, execution, backend=backend, gradients=gradients, np=np)
        row.update(status='accepted', native=admitted)
        native = dict(admitted, status='accepted', attempted=True, task=row['task'], execution=row['execution'],
                      solver_diagnostic_protocol=PROTOCOL, solver_control=task['solver_control'],
                      original_native_result=fixture['pins']['native_result'],
                      original_failure_not_relabelled=True,
                      artifacts={p.name: pin(p) for p in d.iterdir() if p.is_file()})
        row['native_receipt'] = write(d / 'NATIVE_RESULT.json', native)
    except Exception as exc:
        row.update(failure=repr(exc), traceback=traceback.format_exc())
        if (d / 'EXECUTION.json').is_file():
            row['execution'] = pin(d / 'EXECUTION.json')
    row['artifacts'] = {p.name: pin(p) for p in d.iterdir() if p.is_file()}
    row['receipt'] = write(d / 'RESULT.json', row)
    # Native lifecycle failure stops even other declared cells; SCC failure alone
    # retains the other independent, already-declared controls without retry.
    require('execution' in row and read(row['execution'])['owned_process_group_empty'],
            'Native lifecycle incomplete; stop scheduling')
    verify_tree(cfg['runtime']); verify_tree(staged); verify(row['task'])


def compare(a, b, mp, helpers):
    _, _, admission, mapping, np, _ = helpers
    ga = np.asarray(a['qm_gradient']['hartree_per_bohr']) * admission.HA_KCAL / admission.BOHR_A
    gb = np.asarray(b['qm_gradient']['hartree_per_bohr']) * admission.HA_KCAL / admission.BOHR_A
    require(ga.shape == gb.shape == (172, 3) and np.isfinite(ga).all() and np.isfinite(gb).all(), 'Full gradient unavailable')
    ca = np.asarray(a['pointcharge_gradient']['hartree_per_bohr']) * admission.HA_KCAL / admission.BOHR_A
    cb = np.asarray(b['pointcharge_gradient']['hartree_per_bohr']) * admission.HA_KCAL / admission.BOHR_A
    require(ca.shape == cb.shape == (84, 3) and np.isfinite(ca).all() and np.isfinite(cb).all(),
            'Full native PC gradient unavailable')
    pa, pb = mapping.scatter(mp, ga, ca), mapping.scatter(mp, gb, cb)
    require(pa.shape == pb.shape == (256, 3), 'Full physical gradient map differs')
    de = abs(a['energy_hartree'] - b['energy_hartree'])
    dg = float(np.max(np.abs(ga - gb)))
    dc = float(np.max(np.abs(ca - cb)))
    dp = float(np.max(np.abs(pa - pb)))
    qa, qb = np.asarray(a['charges_e']), np.asarray(b['charges_e'])
    require(qa.shape == qb.shape == (172,) and np.isfinite(qa).all() and np.isfinite(qb).all(), 'Charge vectors unavailable')
    return {'energy_absolute_difference_hartree': de, 'maximum_QM_gradient_component_difference_kcal_mol_A': dg,
            'maximum_PC_gradient_component_difference_kcal_mol_A': dc,
            'maximum_physical_gradient_component_difference_kcal_mol_A': dp,
            'gradient_quantity': 'dE/dR; force is negative gradient',
            'native_gradient_units': 'Hartree/bohr',
            'conversion': {'Hartree_to_kcal_mol': admission.HA_KCAL, 'bohr_to_angstrom': admission.BOHR_A,
                           'operation': 'array * Hartree_to_kcal_mol / bohr_to_angstrom'},
            'charge_vectors_e': {'a': qa.tolist(), 'b': qb.tolist()},
            'maximum_charge_difference_e': float(np.max(np.abs(qa - qb))),
            'charge_agreement_is_descriptive': True,
            'pass': de <= GATES['energy_hartree'] and max(dg, dc, dp) <= GATES['gradient_component_kcal_mol_A']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('preflight', 'run'), required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'stage': a.stage, 'status': 'execution_failure', 'config': pin(a.config),
              'implementation': pin(__file__), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'qualification_passed': False,
              'component_kind': 'embedded', 'primary_case': 'failed_005', 'optimizer_resumed': False,
              'declared_calls': 4, 'new_MACE_MM_ALPB_DFT_calls': 0,
              'cells': [{'case_id': c, 'fixture': f, 'broydamp': d, 'status': 'unavailable_not_attempted', 'attempted': False}
                        for c, f, d in CASES], 'checks': []}
    try:
        cfg = json.loads(a.config.read_text())
        helpers = prepare(cfg, out, result)
        if a.stage == 'preflight':
            result['status'] = 'complete'
        else:
            for row in result['cells']:
                run_cell(cfg, out, row, helpers)
            by_id = {r['case_id']: r for r in result['cells']}
            for name, left, right, fixture in (('target_branch_consistency', 'failed_005', 'failed_010', 'failed'),
                    ('preceding_control', 'preceding_005', None, 'preceding'), ('origin_control', 'origin_005', None, 'origin')):
                check = {'name': name, 'left': left, 'right': right or cfg['fixtures'][fixture]['native_result'], 'status': 'unavailable', 'pass': False}
                result['checks'].append(check)
                if by_id[left]['status'] == 'accepted' and (right is None or by_id[right]['status'] == 'accepted'):
                    ref = helpers[-1][fixture]['reference'] if right is None else by_id[right]['native']
                    check.update(compare(by_id[left]['native'], ref, helpers[-1][fixture]['mapping'], helpers), status='complete')
            result['qualification_passed'] = all(r['status'] == 'accepted' for r in result['cells']) and all(c['pass'] for c in result['checks'])
            result['status'] = 'complete' if result['qualification_passed'] else 'unavailable_native_or_agreement'
        verify_tree(cfg)
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    executions = [read(r['execution']) for r in result['cells'] if 'execution' in r]
    result.update(actual_attempted_GFN_calls=sum(r['attempted'] for r in result['cells']),
                  completed_processes=len(executions), native_child_CPU_seconds=sum(r['cpu_user_seconds'] + r['cpu_system_seconds'] for r in executions),
                  native_wall_seconds=sum(r['wall_seconds'] for r in executions), elapsed_seconds=time.monotonic() - tick,
                  allocation={k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK')})
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
