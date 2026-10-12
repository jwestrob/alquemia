"""Independent raw-artifact audit; no model calls or scoring-runner imports."""
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3'
OUT = Path(__file__).resolve().parent
ORDER = ['origin', 'adaptive_Ca', 'adaptive_La']
EV, HA = 23.06054783061903, 627.509474


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(path):
    if isinstance(path, dict):
        assert pin(path['path']) == path
        path = path['path']
    return json.loads(Path(path).read_text())


collection = read(BASE / 'CANONICAL_COLLECTION.json')
manifest = read(collection['manifest'])
tasks = {t['id']: t for t in manifest['tasks']}
cases = [c for c in manifest['cases'] if c['role'] == 'canonical_calibration']
rows = {r['case_id']: r for r in collection['rows']}
reference = read(manifest['historical_reference'])
assert {c['case_id'] for c in cases} == {r['case_id'] for r in reference['rows'] if r['role'] == 'canonical_calibration'}
assert len(cases) == 25
cost = {reuse: {backend: {'calls': 0, 'wall_seconds_sum': 0., 'cpu_seconds_sum': 0., 'peak_process_RSS_KiB': 0} for backend in ('xtb', 'tblite')} for reuse in ('fresh', 'reused')}
absolute_errors, component_errors, population_errors = [], [], []
raw_cache = {}


def raw(task, result_pin):
    if task['id'] in raw_cache:
        return raw_cache[task['id']]
    result = read(result_pin)
    assert result['task']['xyz'] == task['xyz']
    assert all(result['task'][k] == task[k] for k in ('metal', 'charge', 'multiplicity', 'uhf', 'medium', 'valence_electrons', 'atom_count'))
    values, populations = {}, {}
    for backend in ('xtb', 'tblite'):
        stored = result['backends'][backend]
        assert stored['status'] == 'available'
        receipt = read(stored['execution'])
        assert receipt['returncode'] == 0 and receipt['owned_process_group_empty']
        folder = Path(stored['execution']['path']).parent
        assert pin(folder / 'core.xyz')['sha256'] == task['xyz']['sha256']
        assert pin(receipt['stdout']['path']) == receipt['stdout']
        assert pin(receipt['stderr']['path']) == receipt['stderr']
        start = read(folder / 'START.json')
        assert all(v == '1' for v in start['thread_environment'].values())
        assert start['command'] == receipt['command']
        text = Path(receipt['stdout']['path']).read_text()
        if backend == 'xtb':
            assert re.search(r'convergence criteria satisfied|SCC converged', text, re.I)
            assert 'FAILED TO CONVERGE' not in text.upper()
            energy = float(re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)', text)[-1])
            charges = [float(x) for x in (folder / 'charges').read_text().split()]
        else:
            science = read(folder / 'scientific.json')
            assert science['library_returned_without_error']
            assert abs(sum(science['orbital_occupations_e']) - task['valence_electrons']) <= 1e-5
            energy, charges = science['energy_hartree'], science['charges_e']
        assert math.isfinite(energy) and len(charges) == task['atom_count']
        assert all(math.isfinite(x) for x in charges)
        assert abs(sum(charges) - task['charge']) <= 1e-5
        assert energy == stored['energy_hartree'] and charges == stored['charges_e']
        values[backend], populations[backend] = energy, charges
        c = cost['reused' if task['reuse'] else 'fresh'][backend]
        c['calls'] += 1
        c['wall_seconds_sum'] += receipt['wall_seconds']
        c['cpu_seconds_sum'] += receipt['cpu_user_seconds'] + receipt['cpu_system_seconds']
        c['peak_process_RSS_KiB'] = max(c['peak_process_RSS_KiB'], receipt['peak_RSS_KiB'])
    population_error = max(abs(a-b) for a,b in zip(populations['xtb'], populations['tblite']))
    assert population_error <= 1e-5
    population_errors.append(population_error)
    absolute_errors.append(values['xtb'] - values['tblite'])
    raw_cache[task['id']] = values
    return values


independent_rows = []
for case in cases:
    row = rows[case['case_id']]
    signature = read(case['state_signature'])
    matrix = {'Ca': {}, 'La': {}}
    for metal in matrix:
        for candidate, cell in case['matrix'][metal].items():
            assert cell['charge'] == signature['endpoint_charges'][metal]
            assert cell['multiplicity'] == signature['multiplicities'][metal] == 1
            mace = read(cell['MACE'])
            assert mace['energy_eV'] == cell['MACE_eV']
            component = row['components'][metal + '__' + candidate]
            energies = {medium: raw(tasks[cell['inputs'][medium]], component['inputs'][medium]) for medium in ('vacuum', 'water')}
            solvent = {b: energies['water'][b] - energies['vacuum'][b] for b in ('xtb', 'tblite')}
            error = solvent['xtb'] - solvent['tblite']
            assert component['accepted'] and abs(error) <= 1e-6
            assert error == component['component_error_hartree']
            component_errors.append(error)
            matrix[metal][candidate] = (mace['energy_eV'], solvent['xtb'])
    chosen, works = {}, {}
    for metal in matrix:
        zero = matrix[metal]['origin']
        works[metal] = {q: (e[0]-zero[0])*EV+(e[1]-zero[1])*HA for q,e in matrix[metal].items()}
        best = min(ORDER, key=lambda q: works[metal][q])
        chosen[metal] = {'origin': 'origin', 'mathematical': best, 'operational': best if works[metal][best] < -.1 else 'origin'}
    scores = {}
    for policy in ('operational', 'mathematical', 'origin'):
        ca, la = (matrix[z][chosen[z][policy]] for z in ('Ca', 'La'))
        scores[policy] = (ca[0]-la[0])*EV+(ca[1]-la[1])*HA
        assert abs(scores[policy]-row['pool'][policy]['R_model_kcal_mol']) < 1e-7
        for metal in matrix:
            assert chosen[metal][policy] == row['pool'][policy][metal]['candidate']
    independent_rows.append({'case_id': case['case_id'], 'protein': case['protein'], 'expected_class': case['expected_class'], 'scores_kcal_mol': scores, 'chosen': chosen, 'work_kcal_mol': works, 'accommodation_delta_R_kcal_mol': scores['operational']-scores['origin']})

calibration = {}
for policy in ('operational', 'mathematical', 'origin'):
    ca = max((r for r in independent_rows if r['expected_class'] == 'Ca'), key=lambda r:r['scores_kcal_mol'][policy])
    la = min((r for r in independent_rows if r['expected_class'] == 'La'), key=lambda r:r['scores_kcal_mol'][policy])
    gap = la['scores_kcal_mol'][policy] - ca['scores_kcal_mol'][policy]
    assert gap > .02 and abs(gap-collection['calibration'][policy]['gap_kcal_mol']) < 1e-7
    calibration[policy] = {'status': 'available', 'gap_kcal_mol': gap, 'bands': {'Ca_max':ca['scores_kcal_mol'][policy], 'La_min':la['scores_kcal_mol'][policy]}, 'Ca_boundary':ca['case_id'], 'La_boundary':la['case_id']}
start = read(BASE/'CANONICAL_EXECUTION_ATTEMPT_001.json')
finish = read(BASE/'CANONICAL_EXECUTION_ATTEMPT_001_FINISH.json')
assert not finish['exceptions'] and start['workers'] == 14
result = {'source_collection':pin(BASE/'CANONICAL_COLLECTION.json'), 'manifest':collection['manifest'], 'audit_script':pin(__file__), 'canonical_members':25, 'class_counts':dict(Counter(r['expected_class'] for r in independent_rows)), 'calibration':calibration, 'components_qualified':len(component_errors), 'raw_inputs_checked':len(absolute_errors), 'maximum_component_error_Eh':max(map(abs,component_errors)), 'maximum_population_error_e':max(population_errors), 'absolute_v1_gate_failures':sum(abs(e)>1e-6 for e in absolute_errors), 'maximum_absolute_backend_error_Eh':max(map(abs,absolute_errors)), 'cost':cost, 'canonical_execution_wall_seconds':finish['finished_epoch']-start['started_epoch'], 'workers':14, 'rows':independent_rows, 'scientific_calls_by_audit':0, 'interpretation':'Canonical training separation, not predictive validation; all absolute-v1 failures remain failures.'}
(OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
