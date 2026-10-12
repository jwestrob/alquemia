"""Independent raw-artifact audit; no model calls or scoring-runner imports."""
import hashlib
import json
import math
import re
import statistics
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


collection = read(BASE / 'FINAL_COLLECTION.json')
canonical_collection = read(BASE / 'CANONICAL_COLLECTION.json')
assert collection['calibration'] == canonical_collection['calibration']
manifest = read(collection['manifest'])
tasks = {t['id']: t for t in manifest['tasks']}
cases = manifest['cases']
rows = {r['case_id']: r for r in collection['rows']}
reference = read(manifest['historical_reference'])
assert {c['case_id'] for c in cases if c['role'] == 'canonical_calibration'} == {r['case_id'] for r in reference['rows'] if r['role'] == 'canonical_calibration'}
assert len(cases) == 135
cost = {reuse: {backend: {'calls': 0, 'wall_seconds_sum': 0., 'cpu_seconds_sum': 0., 'peak_process_RSS_KiB': 0} for backend in ('xtb', 'tblite')} for reuse in ('canonical', 'transfer', 'reused')}
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
        c = cost['reused' if task['reuse'] else task['stage']][backend]
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
    independent_rows.append({'case_id': case['case_id'], 'protein': case['protein'], 'expected_class': case['expected_class'], 'role': case['role'], 'source_case_id': case['source_case_id'], 'selection_id': case['selection_id'], 'biological_group': case['biological_group'], 'source_conditioning_metal': case['source_conditioning_metal'], 'scores_kcal_mol': scores, 'chosen': chosen, 'work_kcal_mol': works, 'accommodation_delta_R_kcal_mol': scores['operational']-scores['origin']})

calibration = {}
for policy in ('operational', 'mathematical', 'origin'):
    ca = max((r for r in independent_rows if r['expected_class'] == 'Ca' and r['role'] == 'canonical_calibration'), key=lambda r:r['scores_kcal_mol'][policy])
    la = min((r for r in independent_rows if r['expected_class'] == 'La' and r['role'] == 'canonical_calibration'), key=lambda r:r['scores_kcal_mol'][policy])
    gap = la['scores_kcal_mol'][policy] - ca['scores_kcal_mol'][policy]
    assert gap > .02 and abs(gap-collection['calibration'][policy]['gap_kcal_mol']) < 1e-7
    calibration[policy] = {'status': 'available', 'gap_kcal_mol': gap, 'bands': {'Ca_max':ca['scores_kcal_mol'][policy], 'La_min':la['scores_kcal_mol'][policy]}, 'Ca_boundary':ca['case_id'], 'La_boundary':la['case_id']}


def decision(value, bands, expected):
    if value is None:
        return {'decision': None, 'outcome': 'unavailable'}
    d = 'Ca' if value <= bands['Ca_max'] else 'La' if value >= bands['La_min'] else 'inconclusive'
    return {'decision': d, 'outcome': 'inconclusive' if d == 'inconclusive' else 'correct' if d == expected else 'wrong'}


policies = ('operational', 'mathematical', 'origin')
by_id = {r['case_id']: r for r in independent_rows}
for row in independent_rows:
    row['decisions'] = {p: decision(row['scores_kcal_mol'][p], calibration[p]['bands'], row['expected_class']) for p in policies}
    assert row['decisions'] == rows[row['case_id']]['decisions']
manifest_triples = {t['triple_id']: t for t in manifest['triples']}
inventory = read(manifest['historical_inventory'])
assert {t['triple_id'] for t in inventory['triples']} == {t['triple_id'] for t in collection['triples']}
triple_rows = []
for triple in collection['triples']:
    frozen = manifest_triples[triple['triple_id']]
    assert triple['members'] == frozen['members'] and triple['pool_case_ids'] == frozen['pool_case_ids']
    assert triple['prior_preparation_status'] == frozen['status']
    members = [by_id.get(cid) for cid in triple['pool_case_ids']]
    available = frozen['status'] == 'prepared' and all(r is not None for r in members)
    expected = next((r['expected_class'] for r in members if r is not None), None)
    if available:
        assert len(members) == 3 and len({r['expected_class'] for r in members}) == 1
    summary = {'triple_id': triple['triple_id'], 'protein': triple['protein'], 'status': 'available' if available else 'unavailable', 'median_R_kcal_mol': {}, 'spread_kcal_mol': {}, 'decisions': {}}
    for policy in policies:
        scores = [r['scores_kcal_mol'][policy] for r in members] if available else None
        median = statistics.median(scores) if scores is not None else None
        spread = max(scores)-min(scores) if scores is not None else None
        assert median == triple['median_R_kcal_mol'][policy]
        summary['median_R_kcal_mol'][policy] = median
        summary['spread_kcal_mol'][policy] = spread
        summary['decisions'][policy] = decision(median, calibration[policy]['bands'], expected)
        assert summary['decisions'][policy] == triple['decisions'][policy]
    summary['spread_change_accommodated_minus_origin_kcal_mol'] = summary['spread_kcal_mol']['operational']-summary['spread_kcal_mol']['origin'] if available else None
    triple_rows.append(summary)
assert len(triple_rows) == 100
primary_ids = {cid for t in collection['triples'] for cid in t['pool_case_ids'] if cid}
primary = [by_id[cid] for cid in sorted(primary_ids)]
assert len(primary) == 104
assert {(r['source_case_id'], r['selection_id']) for r in primary} == {(r['case_id'], r['selection_id']) for r in inventory['rows']}
assert len({r['source_case_id'] for r in primary}) == 98
counts = {}
for label, population in [('all135', independent_rows), ('canonical25', [r for r in independent_rows if r['role']=='canonical_calibration']), ('noncalibration110', [r for r in independent_rows if r['role']!='canonical_calibration']), ('crystals3', [r for r in independent_rows if r['role']=='crystal_transfer']), ('probes6', [r for r in independent_rows if r['role']=='separate_probe']), ('primary104_context_pairs', primary), ('triples100', triple_rows)]:
    counts[label] = {p: dict(Counter(r['decisions'][p]['outcome'] for r in population)) for p in policies}
paired = dict(Counter((r['decisions']['origin']['outcome']+' -> '+r['decisions']['operational']['outcome']) for r in independent_rows if r['role']!='canonical_calibration'))
valid_triples = [t for t in triple_rows if t['status']=='available']
spread = {'denominator':len(valid_triples), 'decreased':sum(t['spread_change_accommodated_minus_origin_kcal_mol'] < -1e-7 for t in valid_triples), 'increased':sum(t['spread_change_accommodated_minus_origin_kcal_mol'] > 1e-7 for t in valid_triples), 'same_within_1e7':sum(abs(t['spread_change_accommodated_minus_origin_kcal_mol']) <= 1e-7 for t in valid_triples), 'median_origin_range_kcal_mol':statistics.median(t['spread_kcal_mol']['origin'] for t in valid_triples), 'median_accommodated_range_kcal_mol':statistics.median(t['spread_kcal_mol']['operational'] for t in valid_triples), 'minimum_change_kcal_mol':min(t['spread_change_accommodated_minus_origin_kcal_mol'] for t in valid_triples), 'maximum_change_kcal_mol':max(t['spread_change_accommodated_minus_origin_kcal_mol'] for t in valid_triples), 'warning':'94 overlapping triples are correlated, not 94 biological observations'}
parent_groups = {}
for r in primary:
    parent_groups.setdefault(r['protein'], []).append(r)
protein_spread = []
for protein, members in sorted(parent_groups.items()):
    distinct = {r['source_case_id']:r for r in members}
    ranges = {p:max(r['scores_kcal_mol'][p] for r in members)-min(r['scores_kcal_mol'][p] for r in members) for p in policies}
    protein_spread.append({'protein':protein, 'source_context_pairs':len(members), 'unique_sources':len(distinct), 'ranges_kcal_mol':ranges, 'change_kcal_mol':ranges['operational']-ranges['origin']})
stages = {}
for stage in ('canonical','transfer'):
    start = read(BASE/(stage.upper()+'_EXECUTION_ATTEMPT_001.json'))
    finish = read(BASE/(stage.upper()+'_EXECUTION_ATTEMPT_001_FINISH.json'))
    assert not finish['exceptions'] and start['workers'] == 14
    stages[stage] = {'execution_wall_seconds':finish['finished_epoch']-start['started_epoch'], 'workers':14, 'fresh_inputs':start['pending_inputs']}
result = {'source_collection':pin(BASE/'FINAL_COLLECTION.json'), 'canonical_collection':pin(BASE/'CANONICAL_COLLECTION.json'), 'manifest':collection['manifest'], 'audit_script':pin(__file__), 'calibration_identical_to_canonical':True, 'calibration':calibration, 'counts':counts, 'paired_noncalibration_outcomes':paired, 'primary_context_pairs':104, 'primary_unique_structures':98, 'primary_biological_groups':len({r['biological_group'] for r in primary}), 'source_conditioning_counts':dict(Counter(r['source_conditioning_metal'] for r in independent_rows)), 'components_qualified':len(component_errors), 'raw_inputs_checked':len(absolute_errors), 'maximum_component_error_Eh':max(map(abs,component_errors)), 'maximum_population_error_e':max(population_errors), 'absolute_v1_gate_failures':sum(abs(e)>1e-6 for e in absolute_errors), 'maximum_absolute_backend_error_Eh':max(map(abs,absolute_errors)), 'cost':cost, 'stages':stages, 'triple_spread_summary':spread, 'protein_spread':protein_spread, 'rows':independent_rows, 'triples':triple_rows, 'scientific_calls_by_audit':0, 'interpretation':'Consumed development references; replica/triple correlations retained; all absolute-v1 failures remain failures.'}
(OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('rows','triples','protein_spread')},indent=2))
