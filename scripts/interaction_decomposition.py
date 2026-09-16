"""Native, fixed-geometry EDA diagnostics; separate from production scoring."""
from __future__ import annotations

import argparse
import copy
import math
from pathlib import Path
import re
import shutil

import numpy as np

from affordable_common import (HA_TO_KCAL, InvalidArtifact, energy, read_json,
                               record, verify, write_new, xyz, digest)
from affordable_workflow import dry_run

NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'
EDA_LABELS = {
    'bond': 'Bond Energy', 'orbital': 'Orbital Energy',
    'electrostatic': 'Electrostatic Energy', 'pauli': 'Pauli Energy',
    'delta_xc': 'Delta E^0(XC)', 'dispersion': 'Delta Dispersion',
}
PROTOCOL = 'native_r2scan3c_asp303_interaction_eda_v1'


def prepare(states_manifest, endpoint_manifest, plan, agreement, output):
    """Prepare only the declared, consumed 1H4I Asp303 fragment comparison."""
    sm, qm = read_json(states_manifest), read_json(endpoint_manifest)
    if sm['endpoint_manifest'] != record(endpoint_manifest):
        raise InvalidArtifact('state and quantum manifests differ')
    source = {t['task_id']: t for t in qm['tasks']}
    output = Path(output).resolve()
    if output.exists():
        raise InvalidArtifact('preparation directory already exists')
    pending = []
    for metal in ('La', 'Ca'):
        keys = [f'1h4i_qm{size}_{metal}' for size in (33, 36)]
        states = [read_json(verify(sm['states'][key])) for key in keys]
        small, large = [s['core_atoms'] for s in states]
        if (len(small), len(large)) != (47, 54):
            raise InvalidArtifact('not the declared 47/54-atom comparison')
        ids = [[a['id'] for a in atoms] for atoms in (small, large)]
        if any(len(set(group)) != len(group) for group in ids) or not set(ids[0]) <= set(ids[1]):
            raise InvalidArtifact('fragment source IDs are not a unique subset')
        indices = [ids[1].index(key) for key in ids[0]]
        rest = [i for i, key in enumerate(ids[1]) if key not in ids[0]]
        if any(not (large[i]['id'].startswith('A/303/') or large[i]['id'] == 'cap/A:ASP303') for i in rest):
            raise InvalidArtifact('added atoms are not the declared capped Asp303')
        for a, i in zip(small, indices):
            b = large[i]
            if a['element'] != b['element'] or max(abs(x-y) for x, y in zip(a['xyz_A'], b['xyz_A'])) > 1e-10:
                raise InvalidArtifact('shared fragment coordinates changed')
        for key, atoms in zip(keys, (small, large)):
            coords = xyz(verify(source[key]['xyz']))
            if len(coords) != len(atoms) or any(row[0] != atom['element'] or
                max(abs(x-y) for x,y in zip(row[1:],atom['xyz_A'])) > 1e-10 for row, atom in zip(coords,atoms)):
                raise InvalidArtifact('quantum XYZ differs from physical mapping')
        charge1, charge = (-1, -2) if metal == 'La' else (-2, -3)
        if source[keys[0]]['charge'] != charge1 or source[keys[1]]['charge'] != charge:
            raise InvalidArtifact('source charge differs from declared fragments')
        text = ('! r2SCAN-3c NoAutostart DefGrid3 EDA\n%EDA\n'
                ' FRAG1 "r2SCAN-3c NoAutostart DefGrid3 PAL16"\n'
                ' FRAG2 "r2SCAN-3c NoAutostart DefGrid3 PAL16"\n'
                f' FRAG1_C {charge1}\n FRAG1_M 1\n FRAG2_C -1\n FRAG2_M 1\nend\n'
                f'* xyzfile {charge} 1 core.xyz\n%Frag\n Definition\n'
                '  1 {' + ' '.join(map(str, indices)) + '} end\n'
                '  2 {' + ' '.join(map(str, rest)) + '} end\n end\nend\n')
        pending.append((metal, keys, indices, rest, charge, text))
    output.mkdir(parents=True)
    tasks = []
    for metal, keys, indices, rest, charge, text in pending:
        d = output / metal; d.mkdir()
        (d / 'endpoint.inp').write_text(text)
        shutil.copyfile(verify(source[keys[1]]['xyz']), d / 'core.xyz')
        tasks.append({'task_id': '1h4i_asp303_eda_' + metal, 'case': '1h4i_asp303_eda',
            'metal': metal, 'charge': charge, 'multiplicity': 1,
            'input': record(d / 'endpoint.inp'), 'xyz': record(d / 'core.xyz'),
            'output_path': str(d / 'endpoint.out'), 'source_adduct_task': source[keys[1]],
            'source_fragment1_task': source[keys[0]], 'fragment1_indices': indices,
            'fragment2_indices': rest, 'source_states': [sm['states'][key] for key in keys]})
    m = {'schema_version': 'alquemia.interaction_decomposition.v1', 'protocol_id': PROTOCOL,
         'tasks': tasks, 'agreement': record(plan), 'parent_authorization': record(agreement),
         'source_states_manifest': record(states_manifest), 'source_quantum_manifest': record(endpoint_manifest),
         'implementation': record(__file__), 'orca': qm['orca'], 'execution_policy': qm['execution_policy'],
         'compute_budget': None, 'wall_time_limit': None, 'classification': None,
         'internal_scf_count': 'At least six; collect actual generated inputs/outputs.'}
    write_new(output / 'manifest.json', m)
    return dry_run(output / 'manifest.json')


def one_number(text, pattern):
    hits = re.findall(pattern, text, re.M)
    if len(hits) != 1:
        raise InvalidArtifact(f'expected one energy matching {pattern}; found {len(hits)}')
    value = float(hits[0].replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise InvalidArtifact('nonfinite energy component')
    return value


def parse_native_terms(path):
    text = Path(path).read_text()
    result = {'native_total_hartree': energy(path)}
    for name, label in [('scf', r'Total Energy\s*:\s*'),
                        ('dispersion', r'Dispersion correction\s+'),
                        ('gcp', r'gCP correction\s+')]:
        result[name + '_hartree'] = one_number(text, r'^\s*' + label + '(' + NUMBER + ')')
    result['native_closure_kcal_mol'] = HA_TO_KCAL * (result['native_total_hartree'] -
        sum(result[key + '_hartree'] for key in ('scf', 'dispersion', 'gcp')))
    if abs(result['native_closure_kcal_mol']) > .01:
        raise InvalidArtifact('native SCF/D4/gCP accounting does not close')
    return result


def parse_eda(path):
    text = Path(path).read_text()
    native = parse_native_terms(path)
    components = {key: one_number(text, r'^\s*' + re.escape(label) + r'\s+(' + NUMBER + ')')
                  for key, label in EDA_LABELS.items()}
    closure = components['bond'] - sum(value for key, value in components.items() if key != 'bond')
    if abs(closure * HA_TO_KCAL) > .01:
        raise InvalidArtifact('EDA components, including separate delta-XC, do not close')
    return {'native': native, 'eda_hartree': components,
            'eda_kcal_mol': {key: value * HA_TO_KCAL for key, value in components.items()},
            'eda_component_closure_kcal_mol': closure * HA_TO_KCAL}


def inspect_fragment_input(path, task, fragment):
    """Retain native ghost basis centers separately from physical nuclei."""
    text = Path(path).read_text()
    decl = re.findall(r'^\s*\*xyz\s+(-?\d+)\s+(\d+)\s*$', text, re.M)
    expected_charge = task['charge'] + 1 if fragment == 1 else -1
    if decl != [(str(expected_charge), '1')]:
        raise InvalidArtifact('native fragment charge/multiplicity mismatch')
    matches = re.findall(r'^\s*([A-Z][a-z]?)\s*(:?)\((\d+)\)\s+(' + NUMBER +
                         r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', text, re.M)
    source = xyz(verify(task['xyz']))
    if len(matches) != len(source):
        raise InvalidArtifact('native fragment coordinate count differs')
    real_indices = task[f'fragment{fragment}_indices']
    for i, (element, ghost, group, *_coords) in enumerate(matches):
        if element != source[i][0] or bool(ghost) != (i not in real_indices):
            raise InvalidArtifact('native fragment nuclei/ghost mapping differs')
        if int(group) != (1 if i in task['fragment1_indices'] else 2):
            raise InvalidArtifact('native fragment grouping differs')
    delta = np.array([[float(x) for x in row[3:]] for row in matches]) - np.array([row[1:] for row in source])
    # ORCA emits seven decimals after applying a common origin shift. This
    # tolerance checks serialization (not a configurable physical displacement).
    shift = delta.mean(axis=0)
    shape_error = float(np.max(np.abs(delta - shift)))
    if shape_error > 2e-6 or float(np.max(np.abs(delta))) > 2e-6:
        raise InvalidArtifact('native fragment geometry changed beyond serialization precision')
    keyword = re.findall(r'^!(.+)$', text, re.M)
    if keyword != ['r2SCAN-3c NoAutostart DefGrid3 PAL16']:
        raise InvalidArtifact('native fragment method changed')
    return {'input': record(path), 'physical_atom_count': len(real_indices),
            'ghost_atom_count': len(source) - len(real_indices),
            'common_origin_shift_A': shift.tolist(), 'max_shape_serialization_error_A': shape_error,
            'charge': expected_charge, 'multiplicity': 1,
            'basis_reference': 'full_adduct_basis_including_other_fragment_ghost_centers'}


def inventory_outputs(output):
    inventory = []
    for op in sorted(output.parent.glob('*.out')):
        text = op.read_text()
        ip = op.with_suffix('.inp') if op != output else output.parent / 'endpoint.runtime.inp'
        row = {'output': record(op), 'input': record(ip) if ip.exists() else None,
               'normal_termination': 'ORCA TERMINATED NORMALLY' in text,
               'converged_scf_markers': len(re.findall(r'SCF CONVERGED AFTER', text)),
               'nonconvergence_reported': bool(re.search(r'SCF NOT CONVERGED|SCF failed to converge', text, re.I)),
               'final_energy_markers': len(re.findall(r'FINAL SINGLE POINT ENERGY', text)),
               'role': 'adduct' if op == output else 'atomic_reference' if '_atom' in op.name else 'fragment'}
        try:
            row['accepted_native_energy_hartree'] = energy(op)
        except InvalidArtifact as exc:
            row.update(accepted_native_energy_hartree=None, energy_unavailable_reason=str(exc))
        inventory.append(row)
    return inventory


def collect(manifest):
    from global_electrostatic import validate_vacuum_output
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    m, normalized = load_manifest_tasks(Path(manifest))
    if m['protocol_id'] != PROTOCOL or sorted(t['task_id'] for t in m['tasks']) != [
            '1h4i_asp303_eda_Ca', '1h4i_asp303_eda_La']:
        raise InvalidArtifact('not the declared interaction protocol/task set')
    verify(m['agreement']); verify(m['parent_authorization'])
    rows = []
    for task, normal in zip(m['tasks'], normalized):
        row = {'task_id': task['task_id'], 'metal': task['metal'], 'status': 'unavailable'}
        output = Path(task['output_path'])
        receipt = Path(str(output) + '.execution.json')
        row['scf_inventory'] = inventory_outputs(output)
        if output.exists(): row['output'] = record(output)
        if receipt.exists():
            row['execution_receipt'] = record(receipt)
            execution = read_json(receipt)
            if execution.get('manifest') == record(manifest) and execution.get('task_id') == task['task_id']:
                row['execution_finished'] = bool(execution.get('finished_at_utc'))
                row['execution_returncode'] = execution.get('returncode')
        # A completed fragment remains useful provenance even if the later EDA
        # property calculation fails. It is never substituted for the adduct.
        try:
            states = [read_json(verify(item)) for item in task['source_states']]
            for state, name in zip(states, ('source_fragment1_task', 'source_adduct_task')):
                if verify(state['endpoint_output']) != Path(task[name]['output_path']):
                    raise InvalidArtifact('archived reference output mapping differs')
            f1 = output.parent / 'endpoint.runtime_frag1.out'
            inspect_fragment_input(f1.with_suffix('.inp'), task, 1)
            offset = (energy(f1) - energy(verify(states[0]['endpoint_output']))) * HA_TO_KCAL
            row['completed_fragment1_reference_check'] = {
                'output': record(f1), 'archive': states[0]['endpoint_output'],
                'offset_kcal_mol': offset, 'tolerance_kcal_mol': .01, 'passed': abs(offset) <= .01}
        except (ValueError, OSError, KeyError) as exc:
            row['completed_fragment1_reference_check'] = {'status': 'unavailable', 'reason': str(exc)}
        try:
            template = verify(task['input']).read_text()
            if re.findall(r'^!.*$', template, re.M) != ['! r2SCAN-3c NoAutostart DefGrid3 EDA'] or re.search(
                    r'%basis|%pointcharges|CPCM|NumGrad|\bOpt\b', template, re.I):
                raise InvalidArtifact('unsupported interaction Hamiltonian')
            if not _completed_attempt_is_valid(receipt, output, manifest_sha256=digest(manifest), task=normal,
                    runner_identity=m['execution_policy']['task_runner'],
                    runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                raise InvalidArtifact('incomplete or invalid native EDA receipt')
            if 'offset_kcal_mol' not in row['completed_fragment1_reference_check']:
                raise InvalidArtifact('completed fragment/reference provenance unavailable')
            validate_vacuum_output(task['source_adduct_task'], output)
            result = parse_eda(output)
            fragments = {}
            atoms = xyz(verify(task['xyz']))
            for number in (1, 2):
                p = output.parent / f'endpoint.runtime_frag{number}'
                info = inspect_fragment_input(Path(str(p) + '.inp'), task, number)
                op = Path(str(p) + '.out')
                text = op.read_text()
                charge = one_number(text, r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)')
                mult = one_number(text, r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)')
                nel = one_number(text, r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)')
                z = {'H':1, 'C':6, 'N':7, 'O':8, 'Ca':20, 'La':57}
                real = [atoms[i][0] for i in task[f'fragment{number}_indices']]
                expected_nel = sum(z[e] - (46 if e == 'La' else 0) for e in real) - info['charge']
                if charge != info['charge'] or mult != 1 or nel != expected_nel:
                    raise InvalidArtifact('fragment actual electronic state differs')
                fragments[str(number)] = {'mapping': info, 'energies': parse_native_terms(op),
                                          'output': record(op), 'explicit_electrons': nel}
            a, b, c = result['native'], fragments['1']['energies'], fragments['2']['energies']
            interaction = {name: a[name] - b[name] - c[name] for name in
                           ('native_total_hartree', 'scf_hartree', 'dispersion_hartree', 'gcp_hartree')}
            composite_closure = (interaction['native_total_hartree'] -
                                 result['eda_hartree']['bond'] - interaction['gcp_hartree']) * HA_TO_KCAL
            references = {name: energy(task[field]['output_path']) for name, field in
                          [('adduct', 'source_adduct_task'), ('fragment1', 'source_fragment1_task')]}
            offsets = {'adduct_kcal_mol': (a['native_total_hartree'] - references['adduct']) * HA_TO_KCAL,
                       'fragment1_kcal_mol': (b['native_total_hartree'] - references['fragment1']) * HA_TO_KCAL}
            row.update(status='computed', output=record(output), execution_receipt=record(receipt),
                       result=result, fragments=fragments, interaction=interaction,
                       composite_closure_kcal_mol=composite_closure, archived_energies_hartree=references,
                       archived_offsets=offsets,
                       composite_accounting_passed=abs(composite_closure) <= .01,
                       archived_reference_check_passed=all(abs(x) <= .01 for x in offsets.values()))
        except (ValueError, OSError, KeyError) as exc:
            row['reason'] = str(exc)
        rows.append(row)
    pair = {r['metal']: r for r in rows}
    checks, comparison = {}, None
    if set(pair) == {'Ca', 'La'} and all(r['status'] == 'computed' for r in rows):
        ca, la = pair['Ca'], pair['La']
        fragment2_shift = (ca['fragments']['2']['energies']['native_total_hartree'] -
                           la['fragments']['2']['energies']['native_total_hartree']) * HA_TO_KCAL
        checks = {'composite_accounting': all(r['composite_accounting_passed'] for r in rows),
                  'archived_reference_agreement': all(r['archived_reference_check_passed'] for r in rows),
                  'identical_fragment2_energy': abs(fragment2_shift) <= .01}
        comparison = {'Ca_minus_La_components_kcal_mol': {
            key: (ca['result']['eda_hartree'][key] - la['result']['eda_hartree'][key]) * HA_TO_KCAL for key in EDA_LABELS},
            'Ca_minus_La_native_interaction_kcal_mol':
                (ca['interaction']['native_total_hartree'] - la['interaction']['native_total_hartree']) * HA_TO_KCAL,
            'Ca_minus_La_delta_gcp_kcal_mol':
                (ca['interaction']['gcp_hartree'] - la['interaction']['gcp_hartree']) * HA_TO_KCAL,
            'fragment2_Ca_minus_La_kcal_mol': fragment2_shift}
    if comparison is not None:
        status = 'complete_checks_passed' if all(checks.values()) else 'complete_checks_failed'
    elif all(row.get('execution_finished') for row in rows):
        status = 'failed_native_execution' if any(row.get('execution_returncode') != 0 for row in rows) else 'complete_unscorable'
    else:
        status = 'incomplete'
    return {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'implementation': record(__file__),
            'status': status, 'rows': rows, 'checks': checks, 'comparison': comparison,
            'attribution_to_original_partition_supported': bool(checks) and all(checks.values()),
            'affinity_score': None, 'classification': None,
            'claim': 'Native interaction decomposition. Failed reference checks prohibit attributing these components to the original partition discrepancy.'}


def report(result_path, output):
    result = read_json(result_path)
    lines = ['# Native Asp303 interaction diagnostic', '', f"Status: **{result['status']}**.", '',
             'Baseline/default unchanged. No affinity score or calibrated classification.', '',
             'The native fragments include partner ghost basis functions. Reference checks',
             'must pass before attributing these components to the original partition error.', '',
             '| Metal | Completed core-fragment reference shift (kcal/mol) | 0.01 check |',
             '|---|---:|---|']
    for row in result['rows']:
        check = row['completed_fragment1_reference_check']
        offset = check.get('offset_kcal_mol')
        lines.append(f"| {row['metal']} | {offset if offset is not None else 'unavailable'} | {check.get('passed', 'unavailable')} |")
    lines += ['', '## Native execution', '']
    for row in result['rows']:
        inv = row['scf_inventory']
        lines.append(f"- {row['metal']}: {row['status']}; {sum(x['converged_scf_markers'] for x in inv)} converged SCF markers in {len(inv)} native output files, including atomic-reference work.")
        if row.get('reason'): lines.append('  Collection limitation: ' + row['reason'])
    if result['comparison'] is not None:
        lines += ['', '## Native Ca-minus-La components', '', '| Component | kcal/mol |', '|---|---:|']
        for key,value in result['comparison']['Ca_minus_La_components_kcal_mol'].items():
            lines.append(f'| {key} | {value:.12f} |')
        lines += ['', 'The orbital term mixes polarization and charge transfer. The printed',
                  'Pauli term excludes the separately retained delta-XC. gCP is audited',
                  'separately against adduct-minus-fragments native energies.', '',
                  f"Checks: {result['checks']}"]
    else:
        lines += ['', 'No complete two-metal decomposition is available. Missing components stay unavailable.']
    lines += ['', 'Result record: ' + str(Path(result_path).resolve()), '']
    with Path(output).open('x') as f:
        f.write('\n'.join(lines))
    return {'status': 'report_written'}


def prepare_order_retry(manifest, output):
    """Recover ORCA's required geometry-before-fragments ordering only."""
    m = copy.deepcopy(read_json(manifest))
    output = Path(output).resolve()
    if output.exists():
        raise InvalidArtifact('retry directory already exists')
    prepared = []
    for task in m['tasks']:
        source = verify(task['input']).read_text()
        old_output = Path(task['output_path'])
        if 'fragments (in COORDS) should come after the geometry input' not in old_output.read_text():
            raise InvalidArtifact('not the recognized input-order failure')
        if 'FINAL SINGLE POINT ENERGY' in old_output.read_text():
            raise InvalidArtifact('unexpected scientific result in input-order failure')
        match = re.fullmatch(r'(.+?)(%Frag\n.+?end\n)(\* xyzfile -?\d+ 1 core\.xyz\n)', source, re.S)
        if not match:
            raise InvalidArtifact('unexpected fragment input syntax')
        corrected = match[1] + match[3] + match[2]
        prepared.append((task, corrected, verify(task['xyz']), record(old_output)))
    output.mkdir(parents=True)
    for task, corrected, xyz, failed in prepared:
        d = output / task['metal']; d.mkdir()
        (d / 'endpoint.inp').write_text(corrected)
        shutil.copyfile(xyz, d / 'core.xyz')
        task.update(input=record(d / 'endpoint.inp'), xyz=record(d / 'core.xyz'),
                    output_path=str(d / 'endpoint.out'), failed_attempt=failed)
    m.update(retry_of=record(manifest), implementation=record(__file__),
             retry_reason='Geometry before %Frag as required by native ORCA parser; identical method, fragments, charges and coordinates.')
    write_new(output / 'manifest.json', m)
    return dry_run(output / 'manifest.json')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='operation', required=True)
    prep = sub.add_parser('prepare')
    for name in ('states-manifest', 'endpoint-manifest', 'plan', 'agreement', 'output'):
        prep.add_argument('--' + name, required=True)
    retry = sub.add_parser('prepare-order-retry')
    retry.add_argument('--manifest', required=True)
    retry.add_argument('--output', required=True)
    coll = sub.add_parser('collect')
    coll.add_argument('--manifest', required=True)
    coll.add_argument('--output', required=True)
    rep = sub.add_parser('report')
    rep.add_argument('--result', required=True)
    rep.add_argument('--output', required=True)
    args = p.parse_args()
    if args.operation == 'prepare':
        result = prepare(args.states_manifest, args.endpoint_manifest, args.plan, args.agreement, args.output)
    elif args.operation == 'prepare-order-retry':
        result = prepare_order_retry(args.manifest, args.output)
    elif args.operation == 'collect':
        result = collect(args.manifest)
        write_new(args.output, result)
    else:
        result = report(args.result, args.output)
    print(result['status'])


if __name__ == '__main__':
    main()
