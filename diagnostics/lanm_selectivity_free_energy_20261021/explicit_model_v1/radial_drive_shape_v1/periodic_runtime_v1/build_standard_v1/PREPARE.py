"""Metadata-only scheduling recovery; never compiles, submits or runs a model."""
import copy
import datetime
import difflib
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
DECLARATION_SHA = '2deb466da81262f5337047f1ac6319605eb296e5ce1be0d4c278d1c3cecbabfd'
OLD_QOS = "'Account':'banfield','QOS':'standard_expand','Partition':'standard'"
NEW_QOS = "'Account':'banfield','QOS':'standard','Partition':'standard'"


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    assert pin(p['path']) == p, 'Changed authority'
    return json.loads(Path(p['path']).read_text())


def save(path, obj):
    with path.open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def main(out):
    tick, cpu = time.monotonic(), time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    dp = pin(HERE/'ROOT_DECLARATION.json')
    assert dp['sha256'] == DECLARATION_SHA
    d = read(dp)
    old_ready, old_inputs = read(d['old_ready']), read(d['old_inputs'])
    old_decl = read(d['old_declaration'])
    assert d['unchanged_build_implementation'] == old_ready['implementation'] == old_inputs['implementation']
    assert pin(d['unchanged_build_implementation']['path']) == d['unchanged_build_implementation']
    assert old_inputs['sources']['allocation12'] == d['old_allocation_parser']
    assert old_inputs['sources']['RUN.sbatch'] == d['old_wrapper']
    controls = read(old_ready['controls'])
    assert controls['status'] == 'passed_periodic_v2_build_source_controls' and controls['check_count'] == 30
    policy_pin = json.loads((HERE/'POLICY_AUTHORITY.json').read_text())
    policy = read(policy_pin)
    assert policy['status'] == 'read_only_scheduler_policy_snapshot' and policy['scheduler_mutations'] == 0
    for record in policy['commands'].values():
        assert record['returncode'] == 0
        assert pin(record['stdout']['path']) == record['stdout']
        assert pin(record['stderr']['path']) == record['stderr']

    original = Path(d['old_allocation_parser']['path']).read_text()
    assert pin(d['old_allocation_parser']['path']) == d['old_allocation_parser']
    assert original.count(OLD_QOS) == 1
    changed = original.replace(OLD_QOS, NEW_QOS, 1)
    assert changed.replace(NEW_QOS, OLD_QOS, 1) == original
    with (HERE/'ALLOCATION_STANDARD.py').open('x') as f:
        f.write(changed)
    patch = ''.join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True),
                                       fromfile=d['old_allocation_parser']['path'], tofile=str(HERE/'ALLOCATION_STANDARD.py')))
    with (HERE/'ALLOCATION_STANDARD.diff').open('x') as f:
        f.write(patch)
    wrapper = Path(d['old_wrapper']['path']).read_text()
    assert pin(d['old_wrapper']['path']) == d['old_wrapper']
    assert wrapper.count('#SBATCH --qos=standard_expand\n') == 1
    with (HERE/'RUN.sbatch').open('x') as f:
        f.write(wrapper.replace('#SBATCH --qos=standard_expand\n', '#SBATCH --qos=standard\n', 1))

    # Preserve the old schema required by unchanged BUILD.source_context.
    binding = copy.deepcopy(old_decl)
    binding['status'] = 'declared_candidate_independent_build_preparation_standard_recovery'
    binding['prospective_build']['qos'] = 'standard'
    binding['resource_recovery'] = {'root_declaration': dp, 'old_declaration': d['old_declaration'],
        'policy_snapshot': policy_pin, 'scheduler_defaults': 'standard default wall limit3days; default memory4000MiB/CPU; no flags override either',
        'old_job_resolution': 'Root-owned recheck and explicit resolution before any replacement submission; this preparation cannot mutate20739.'}
    save(out/'DECLARATION_BINDING.json', binding)
    x = copy.deepcopy(old_inputs)
    x['declaration'] = pin(out/'DECLARATION_BINDING.json')
    x['sources']['allocation12'] = pin(HERE/'ALLOCATION_STANDARD.py')
    x['sources']['RUN.sbatch'] = pin(HERE/'RUN.sbatch')
    x['sources']['standard_recovery_preparer'] = pin(__file__)
    x['sources']['standard_recovery_controls'] = pin(HERE/'TEST_POLICY.py')
    x['scheduling_recovery'] = {'declaration': dp, 'old_inputs': d['old_inputs'], 'old_ready': d['old_ready'],
        'old_controls_reused': old_ready['controls'], 'policy_snapshot': policy_pin}
    save(out/'INPUTS.json', x)
    cfg = read(old_ready['configuration_template'])
    cfg['inputs'] = pin(out/'INPUTS.json')
    cfg['scheduling_recovery_declaration'] = dp
    cfg['old_job_resolution_required_before_submission'] = True
    save(HERE/'CONFIG_TEMPLATE.json', cfg)
    release = read(old_ready['release_template'])
    release['scheduling_recovery_declaration'] = dp
    release['pending_old_job_resolution'] = None
    save(HERE/'RELEASE_TEMPLATE.json', release)
    result = {'status': 'prepared_standard_QOS_build_recovery_metadata', 'implementation': pin(__file__),
        'declaration': dp, 'inputs': pin(out/'INPUTS.json'), 'declaration_binding': pin(out/'DECLARATION_BINDING.json'),
        'allocation': pin(HERE/'ALLOCATION_STANDARD.py'), 'wrapper': pin(HERE/'RUN.sbatch'),
        'unchanged_BUILD': d['unchanged_build_implementation'], 'old_controls_reused': old_ready['controls'],
        'policy_snapshot': policy_pin, 'compiler_native_solver_calls': 0, 'scheduler_mutations': 0,
        'old_job_preserved': True, 'cost': {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'threads': 1}}
    save(out/'RESULT.json', result)
    print(json.dumps({'result': pin(out/'RESULT.json'), 'inputs': result['inputs']}))


if __name__ == '__main__':
    import sys
    main(Path(sys.argv[1]).resolve())
