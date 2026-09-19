"""Native CPCM DFT on all four charge-changing compact PQQ contexts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from hydration_square import endpoint
from mace_hybrid import check_atoms

CASES = ('a0acd6b9f2-pqq-la_model', 'a8r3s4-pqq-la_model', 'q88jh5-pqq-la_model', 'q9z4j7-pqq-la_model')
PROTOCOL = 'native_r2scan3c_CPCM_complete_context_charge_PQQ_v1'
HEADER = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3'


def prepare(source, agreement, output):
    source = Path(source).resolve(); old = read_json(source)
    config = read_json(verify(read_json(verify(old['inventory']))['static_config']))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dest = impl/p.name; shutil.copyfile(p, dest); pins[p.name] = record(dest)
    chosen = [c for c in old['cases'] if read_json(verify(c['preparation']))['added_formal_charge'] == -1]
    if tuple(c['source']['case'] for c in chosen) != CASES:
        raise InvalidArtifact('all-four charge-changing context selection differs')
    tasks = []
    for c in chosen:
        s = c['source']; prep = read_json(verify(c['preparation']))
        for metal in ('Ca', 'La'):
            original = s['endpoints'][metal]
            source_lines = [line.strip() for line in verify(original['input']).read_text().splitlines()
                            if line.strip() and not line.strip().startswith('#')]
            if ([line for line in source_lines if line.startswith('!')] != [HEADER]
                    or any(not line.startswith(('!', '%maxcore ', '* xyzfile ')) for line in source_lines)):
                raise InvalidArtifact('archived native single-point recipe differs')
            task = next(t for t in old['tasks'] if t['case_id'] == s['case'] and t['metal'] == metal)
            d = out/'tasks'/task['task_id']; d.mkdir(parents=True)
            xp = d/'core.xyz'; shutil.copyfile(verify(task['xyz']), xp)
            ip = d/'endpoint.inp'; ip.write_text(HEADER+'\n'+f"* xyzfile {task['charge']} 1 core.xyz\n")
            tasks.append({'task_id': task['task_id'], 'case': s['case'], 'metal': metal,
                          'charge': task['charge'], 'multiplicity': 1, 'atom_count': prep['new_atom_count'],
                          'xyz': record(xp), 'input': record(ip), 'output_path': str(d/'endpoint.out'),
                          'preparation': c['preparation'], 'source_OMOL_task': task})
    manifest = {'protocol_id': PROTOCOL, 'source_context_manifest': record(source), 'agreement': record(agreement),
                'states': chosen, 'tasks': tasks, 'implementation': pins, 'orca': config['orca'],
                'execution_policy': config['execution_policy'], 'execution_resources': {'mpi_ranks': 16, 'concurrent_tasks': 4},
                'reference': None, 'calibrated_decision': None, 'compute_budget': None, 'new_DFT_endpoints': 8,
                'baseline_changed': False, 'source_MACE_result': record(source.parent/'result_1202474.json'),
                'neutral_DFT_manifest': record(source.parents[2]/'second_shell_20260919/prepared_v2/manifest.json'),
                'neutral_DFT_collection': record(source.parents[2]/'second_shell_20260919/prepared_v2/collection_1202082.json')}
    write_new(out/'manifest.json', manifest)
    return validate(out/'manifest.json')


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or tuple(s['source']['case'] for s in m['states']) != CASES or len(m['tasks']) != 8:
        raise InvalidArtifact('finite native DFT scope differs')
    for pin in m['implementation'].values(): verify(pin)
    for key in ('source_context_manifest', 'source_MACE_result', 'neutral_DFT_manifest', 'neutral_DFT_collection'): verify(m[key])
    for s in m['states']:
        prep = read_json(verify(s['preparation']))
        if prep['added_formal_charge'] != -1: raise InvalidArtifact('selected charge effect changed')
        for metal in ('Ca', 'La'):
            original = s['source']['endpoints'][metal]
            e = endpoint(original['output'], original['receipt'], original['xyz'], original['input'])
            if e['energy_hartree'] != original['energy_hartree']: raise InvalidArtifact('actual original-core energy differs')
            t = next(t for t in m['tasks'] if t['case'] == s['source']['case'] and t['metal'] == metal)
            if t['xyz']['sha256'] != t['source_OMOL_task']['xyz']['sha256'] or t['charge'] != original['charge']-1:
                raise InvalidArtifact('context geometry/charge differs from actual native OMOL input')
            rows = xyz(verify(t['xyz'])); check_atoms(rows, t['charge'])
            for i, j in prep['core_to_context'].items():
                if rows[j] != xyz(verify(original['xyz']))[int(i)]: raise InvalidArtifact('original core coordinates changed')
    from affordable_workflow import dry_run
    return dry_run(manifest)


def dielectric(pin):
    text = verify(pin).read_text()
    terms = re.findall(r'^CPCM Dielectric\s*:\s*([-+0-9.]+)\s+Eh', text, re.M)
    if len(terms) != 1: raise InvalidArtifact('expected one native CPCM dielectric diagnostic')
    return float(terms[0])


def compare_state(case, core, expanded, mace):
    raw = {}; dielectric_components = {}
    for metal in ('Ca', 'La'):
        c, e = core[metal], expanded[metal]
        raw[metal] = {'core_energy_hartree': c['energy_hartree'], 'context_energy_hartree': e['energy_hartree'],
                      'context_minus_core_hartree': e['energy_hartree']-c['energy_hartree']}
        dc, de = dielectric(c['output']), dielectric(e['output'])
        dielectric_components[metal] = {'core_hartree': dc, 'context_hartree': de, 'delta_hartree': de-dc}
    core_r = core['Ca']['energy_hartree']-core['La']['energy_hartree']
    expanded_r = expanded['Ca']['energy_hartree']-expanded['La']['energy_hartree']
    return {'case': case, 'endpoints': {'core': core, 'context': expanded}, 'endpoint_changes': raw,
            'core_R_hartree': core_r, 'context_R_hartree': expanded_r,
            'delta_R_kcal_mol': (expanded_r-core_r)*HA_TO_KCAL,
            'MACE_core_R_model_kcal_mol': mace['core_R_model_kcal_mol'],
            'MACE_context_R_model_kcal_mol': mace['context_R_model_kcal_mol'],
            'MACE_delta_R_model_kcal_mol': mace['context_R_model_kcal_mol']-mace['core_R_model_kcal_mol'],
            'CPCM_dielectric_diagnostic': dielectric_components,
            'CPCM_dielectric_delta_contrast_kcal_mol': (dielectric_components['Ca']['delta_hartree']-dielectric_components['La']['delta_hartree'])*HA_TO_KCAL,
            'CPCM_term_already_in_total_energy': True, 'unique_physical_component_attribution': None}


def collect(manifest, output):
    m = read_json(manifest); rows = []; mr = read_json(verify(m['source_MACE_result']))
    for s in m['states']:
        eps = {}; core = s['source']['endpoints']; case = s['source']['case']
        for metal in ('Ca', 'La'):
            t = next(t for t in m['tasks'] if t['case'] == case and t['metal'] == metal)
            eps[metal] = endpoint(record(t['output_path']), record(t['output_path']+'.execution.json'), t['xyz'], t['input'])
        row = compare_state(case, core, eps, next(r for r in mr['rows'] if r['case'] == case))
        row.update(added_formal_charge=-1, expected_class=s['source']['expected_class'], evaluation_role='consumed_development')
        rows.append(row)
    old = read_json(verify(m['neutral_DFT_collection'])); om = read_json(verify(m['neutral_DFT_manifest']))
    cfg = read_json(verify(om['configuration']))
    for case in ('1H4I', '4MAE'):
        core = next(s['endpoints'] for s in cfg['cases'] if s['case'] == case)
        expanded = next(r['endpoints'] for r in old['rows'] if r['case'] == case)
        row = compare_state(case, core, expanded, next(r for r in mr['rows'] if r['case'] == case))
        row.update(added_formal_charge=0, expected_class='Ca' if case=='1H4I' else 'La', evaluation_role='reused_neutral_context_control')
        rows.append(row)
    result = {'manifest': record(manifest), 'rows': rows, 'new_DFT_endpoints': 8, 'new_MACE_endpoints': 0,
              'baseline_changed': False, 'calibrated_bands': None, 'reference': None,
              'interpretation': 'Matched context effects; composition/charge/cavity vary together, no unique cause inferred.'}
    write_new(output, result); return result


def main():
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('prepare')
    for key in ('source', 'agreement', 'output'): q.add_argument('--'+key, required=True)
    q = sub.add_parser('validate'); q.add_argument('--manifest', required=True)
    q = sub.add_parser('collect'); q.add_argument('--manifest', required=True); q.add_argument('--output', required=True)
    a = vars(p.parse_args()); operation = a.pop('command'); print(json.dumps(globals()[operation](**a), indent=2))


if __name__ == '__main__': main()
