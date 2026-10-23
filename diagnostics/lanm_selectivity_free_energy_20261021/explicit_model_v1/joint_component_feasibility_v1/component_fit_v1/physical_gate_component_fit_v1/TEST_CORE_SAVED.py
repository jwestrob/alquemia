"""Focused saved-data policy replay; solver construction is forbidden."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent


def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def read(item):
    assert pin(item['path']) == item
    return json.loads(Path(item['path']).read_text())


def main(out):
    wall, cpu = time.monotonic(), time.process_time()
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('tested_physical_core', HERE/'CORE.py')
    core = importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
    calls = []
    def forbidden(*args, **kwargs):
        calls.append('solver_construction')
        raise AssertionError('No solver construction permitted in saved replay')
    core.clarabel.DefaultSolver = forbidden
    inputs = json.loads((HERE.parent/'conic_saved_comparison_v1/INPUTS.json').read_text())
    values = {k: read(v) for k, v in inputs['vectors'].items()}
    origin = values[inputs['origin_key']]
    declaration = read(inputs['campaign_declaration'])
    declaration['proposal'].update(method='Clarabel', settings=core.DEPENDENCIES['settings'])
    fixed = core.make_fixed(origin, read(inputs['manifest']), read(inputs['decomposition']),
                            read(inputs['domain']), declaration)
    probes = [values[k] for k in inputs['derivative_keys']]
    proposal = core.propose(origin, probes, fixed)
    saved = read(core.DEPENDENCIES['saved_initial_proposal'])
    checks = []
    def check(name, value):
        checks.append(dict(name=name, passed=bool(value)))
        assert value, name
    check('exact saved proposal reused without solver', proposal['solver_calls'] == 0 and
          proposal.get('imported_saved_proposal') == core.DEPENDENCIES['saved_initial_proposal'] and not calls)
    check('raw solution and status preserved', proposal['solution'] == saved['solution'] and
          proposal['solver_status'] == saved['solver_status'] and proposal['solver_success'] == saved['solver_success'])
    check('four ordered trial slots', [s['fraction'] for s in proposal['trial_slots']] == list(core.FRACTIONS))
    for slot, witness in zip(proposal['trial_slots'], saved['witnesses']):
        check('unchanged physical fraction '+str(slot['fraction']),
              slot['parameters'] == witness['constructed_source_certificate']['parameters'] and
              slot['key'] == witness['constructed_source_certificate']['key'])
        check('strict retained constraints '+str(slot['fraction']),
              slot['eligible'] and slot['certificate']['strict_retained_feasible'] and
              slot['certificate']['retained_scalar_margin_count'] == 847)
    full = proposal['trial_slots'][0]
    check('old ceiling failure is descriptive', not full['certificate']['fixed_old_ceiling_pass'] and full['eligible'])
    pilot_dir = Path(core.DEPENDENCIES['saved_initial_proposal']['path']).parent.parent/'nonlinear_check_v1'
    pilot = json.loads((pilot_dir/'run_v1/execution/RESULT.json').read_text())
    review = json.loads((pilot_dir/'run_v1/ROOT_ACTUAL_REVIEW.json').read_text())
    check('exact completed reviewed native pilot', review['result'] == pin(pilot_dir/'run_v1/execution/RESULT.json') and
          pilot['status'] == 'complete_one_ceiling_tradeoff_native_pilot' and pilot['native_programs'] == 56)
    vector = copy.deepcopy(pilot['answer'])
    # The finite pilot stores these two envelope fields outside its assembled answer.
    vector.update(parameters=pilot['candidate']['parameters'], status='complete')
    check('native vector parameter identity', vector['parameters'] == full['parameters'])
    measured = core.assess(vector, fixed)
    check('physical failure retained despite better q', measured['q'] < 1 and not measured['physical_constraints_pass'] and
          sum(g['passed'] for g in measured['original25'].values()) == 24)
    slots = copy.deepcopy(proposal['trial_slots'])
    for i, slot in enumerate(slots):
        if i:
            slot['eligible'] = False
        else:
            slot.update(execution_status='complete', native_admitted=True, vector=vector)
    selected = core.select(origin, slots, fixed)
    check('failed physical trial cannot replace incumbent', selected['selected_trial_index'] is None and
          selected['selected_key'] == inputs['origin_key'])
    check('exact selection epsilon retained', selected['strict_q_improvement_tolerance'] == 1e-10)
    result = dict(status='passed_saved_physical_component_policy_replay', implementation=pin(__file__),
                  core=pin(HERE/'CORE.py'), checks=checks, native_calls=0, solver_calls=0,
                  saved_first_proposal=core.DEPENDENCIES['saved_initial_proposal'],
                  first_proposal=proposal, failed_pilot_metrics=measured,
                  cost=dict(wall_seconds=time.monotonic()-wall, CPU_seconds=time.process_time()-cpu, threads=1))
    (out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], checks=len(checks), result=pin(out/'RESULT.json'), cost=result['cost'])))


if __name__ == '__main__':
    main(sys.argv[1])
