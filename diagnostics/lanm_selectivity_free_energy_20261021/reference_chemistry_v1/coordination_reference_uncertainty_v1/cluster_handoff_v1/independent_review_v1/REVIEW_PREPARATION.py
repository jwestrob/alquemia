"""One saved-only literal-input review; does not import or execute the producer."""
import hashlib
import json
import time
from pathlib import Path

START = time.monotonic()
CPU = time.process_time()
BASE = Path(__file__).resolve().parents[6]
REL = Path('reference_chemistry_v1/coordination_reference_uncertainty_v1/cluster_handoff_v1')
D = BASE / 'diagnostics/lanm_selectivity_free_energy_20261021' / REL
W = BASE / 'workspaces/lanm_selectivity_free_energy_20261021' / REL

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def read(p):
    return json.loads(Path(p).read_text())

checks = {}
def check(name, value):
    assert value, name
    checks[name] = True

ready = read(D / 'READY.json')
check('frozen_READY', pin(D / 'READY.json')['sha256'] == '16c07dfa8a2c02dd30268eaa16e6de3b43b39a55f2220754f220e09d61e59fc8')
for name, authority in ready.items():
    if isinstance(authority, dict) and set(authority) == {'path', 'sha256'}:
        check('pin_' + name, pin(authority['path']) == authority)
p = read(ready['panel']['path'])
old = read(ready['parent_panel']['path'])
prep = read(ready['preparation']['path'])
controls = read(ready['saved_controls']['path'])
prior = [c for c in old['cells'] if c['stage'] == 'cluster']
check('exact_eight_ordered_ids', len(p['cells']) == len(prior) == 8 and [c['cell_id'] for c in p['cells']] == [c['cell_id'] for c in prior])
for new, before in zip(p['cells'], prior):
    cid = new['cell_id']
    check(cid + '_exact_scientific_definition', {k:v for k,v in new.items() if k != 'template_input'} == {k:v for k,v in before.items() if k != 'template_input'})
    for tag, c in [('old', before), ('new', new)]:
        check(cid + '_' + tag + '_input_pin', pin(c['template_input']['path']) == c['template_input'])
    a = Path(before['template_input']['path']).read_bytes()
    b = Path(new['template_input']['path']).read_bytes()
    check(cid + '_only_MPI_delta', a.count(b'%pal nprocs 4 end') == 1 and b == a.replace(b'%pal nprocs 4 end', b'%pal nprocs 24 end'))
    check(cid + '_original_input_join', p['original_inputs'][cid] == before['template_input'])
for name in ['assets', 'orca', 'native_modules', 'saved_complexes', 'tiers']:
    check('inherited_' + name, p[name] == old[name])
for name in ['parent_panel', 'qualified_controls', 'control_actual_review', 'control_root_review', 'failed_tight_controls', 'cost_accounting']:
    check('panel_' + name, p[name] == ready[name])
check('single_tighter_W9_Dy_scout', p['scout_cell_id'] == ready['scout_cell_id'] == 'cluster__tighter__aquo9_central_Dy')
check('original_scout_retained', p['original_scout_cell_id'] == old['scout_cell_id'] == 'cluster__tight__aquo9_central_Dy')
check('remaining_seven', ready['remaining_cell_ids'] == [c['cell_id'] for c in prior if c['cell_id'] != p['scout_cell_id']])
check('resources', p['ranks'] == ready['ranks'] == 24 and p['nodes'] == ready['nodes'] == 1 and ready['threads_per_rank'] == 1)
check('default_grants', ready['memory_override'] is None and ready['time_override'] is None)
check('preparation_complete_zero_native', prep['status'] == 'passed_cluster_handoff_preparation' and prep['declared'] == 8 and prep['new_native_calls'] == 0)
for name in ['implementation', 'panel', 'configuration']:
    check('preparation_' + name, prep[name] == ready[name])
check('saved_controls_read_not_rerun', controls['status'] == 'passed_cluster_handoff_saved_controls' and controls['check_count'] == len(controls['checks']) == 25 and all(v is True for v in controls['checks'].values()) and controls['new_native_calls'] == 0)
for name in ['implementation', 'panel', 'preparation', 'test_source']:
    check('controls_' + name, controls[name] == ready[name])
tight = read(ready['failed_tight_controls']['path'])
tighter = read(ready['qualified_controls']['path'])
check('old_Tight_failed_preserved', tight['status'] == 'local_controls_failed' and tight['complete'] == 18 and len(tight['checks']) == 35 and sum(tight['checks'].values()) == 32 and tight['all_controls_passed'] is False)
check('tighter18_passed', tighter['status'] == 'tighter_local_controls_passed' and tighter['complete'] == 18 and tighter['new_complete'] == 15 and tighter['reused_complete'] == 3 and len(tighter['checks']) == 35 and all(tighter['checks'].values()))
check('same_35_control_names', set(tighter['checks']) == set(tight['checks']))
receipt = {
    **ready['source_review_contract'],
    'ready': pin(D / 'READY.json'),
    'reviewer_implementation': pin(__file__),
    'declaration': ready['declaration'],
    'saved_controls': ready['saved_controls'],
    'checks': checks,
    'check_count': len(checks),
    'manual_source_review': {
        'original_deck_and_recovered_parser_unchanged': True,
        'actual_tighter18_and_failed_Tight18_authorities_separate': True,
        'exact_scout_then_remaining7_actual_review_and_cost_memory_scratch_gate': True,
        'resource_runtime_24_MPI_one_node_default_memory_MaxCore_floor_0_75_grant_per_rank': True,
        'immutable_output_paths_no_retry_and_raw_replay': True,
        'all8_required_before_signed_Dy_minus_La_W9_minus_W8_contrast': True,
        'original_0_25_kcal_cutoff_sensitivity_gate_retained': True,
        'CCSD_T1_and_PBE0_electronic_D4_quantities_kept_separate': True,
        'application_cost_not_claimed_full_scheduler_cost': True
    },
    'blocking_findings': [],
    'limitations': [
        'Source/preparation review only; no cluster native output or feasibility is admitted.',
        'Tight remains the lower-accuracy comparator with its original three failed controls.',
        'Small-control qualification and a cutoff-sensitivity check do not establish physical accuracy, cluster CP/B, or a unique electronic state.',
        'Remaining seven cells require actual scout admission and root cost, memory and scratch acceptance.'
    ],
    'author_suites_rerun': 0,
    'new_native_calls': 0,
    'threads': 1,
    'wall_seconds': time.monotonic() - START,
    'CPU_seconds': time.process_time() - CPU
}
out = W / 'independent_review_v1'
out.mkdir(parents=True, exist_ok=True)
target = out / 'REVIEW.json'
with target.open('x') as f:
    json.dump(receipt, f, indent=2); f.write('\n')
print(json.dumps({'review': pin(target), 'checks': len(checks)}))
