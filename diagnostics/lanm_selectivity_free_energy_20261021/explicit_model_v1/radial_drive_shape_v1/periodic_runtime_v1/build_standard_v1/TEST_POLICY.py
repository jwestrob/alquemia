"""Focused manufactured checks for the one-literal scheduling recovery only."""
import ast
import copy
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path
import PREPARE as P


def main(prepared, out):
    tick, cpu = time.monotonic(), time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    result = {'status': 'started', 'implementation': P.pin(__file__), 'checks': checks,
              'compiler_native_solver_calls': 0, 'scheduler_queries_or_mutations': 0, 'error': None}
    def check(ok, label):
        assert ok, label
        checks.append(label)
    try:
        d = P.read(P.pin(P.HERE/'ROOT_DECLARATION.json'))
        x = P.read(P.pin(prepared/'INPUTS.json')); old = P.read(d['old_inputs'])
        old_text = Path(d['old_allocation_parser']['path']).read_text()
        new_text = Path(x['sources']['allocation12']['path']).read_text()
        check(new_text == old_text.replace(P.OLD_QOS, P.NEW_QOS, 1), 'allocation source changes exactly one expectedQOS literal')
        check(new_text.replace(P.NEW_QOS, P.OLD_QOS, 1) == old_text, 'inverse allocation patch restores original exact bytes')
        old_ast, new_ast = ast.parse(old_text), ast.parse(new_text)
        for before, after in zip(old_ast.body, new_ast.body):
            if isinstance(before, ast.FunctionDef) and before.name == 'parse_allocation':
                continue
            assert ast.dump(before) == ast.dump(after)
        check(len(old_ast.body) == len(new_ast.body), 'all non-parser AST nodes unchanged')
        spec = importlib.util.spec_from_file_location('manufactured_allocation_only', x['sources']['allocation12']['path'])
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        env = {'SLURM_JOB_ID': '123', 'SLURM_CPUS_PER_TASK': '12', 'SLURM_NTASKS': '1', 'SLURM_JOB_NUM_NODES': '1'}
        raw = ('JobId=123 JobState=RUNNING Account=banfield QOS=standard Partition=standard '
               'NumNodes=1 NumCPUs=12 NumTasks=1 CPUs/Task=12 Requeue=0 OverSubscribe=NO '
               'BatchHost=manufactured-node UserId=manufactured(601) AllocTRES=cpu=12,mem=48000M,node=1 MinMemoryCPU=4000M')
        rec = mod.parse_allocation(raw, env, 'manufactured-node', 601)
        check(rec['qos'] == 'standard' and rec['allocated_CPUs'] == 12 and rec['memory_grant_bytes'] == 48000*1024**2,
              'new parser accepts manufactured exact standard allocation and reports truthful QOS')
        try:
            mod.parse_allocation(raw.replace('QOS=standard ', 'QOS=standard_expand '), env, 'manufactured-node', 601)
        except AssertionError as exc:
            check(str(exc) == 'allocation QOS', 'new parser rejects old QOS rather than bypassing resource gate')
        else:
            raise AssertionError('Old QOS was accepted')
        check(P.pin(d['unchanged_build_implementation']['path']) == d['unchanged_build_implementation'] == x['implementation'],
              'BUILD.py exact original implementation retained')
        wrapper = Path(x['sources']['RUN.sbatch']['path']).read_text()
        old_wrapper = Path(d['old_wrapper']['path']).read_text()
        check(wrapper == old_wrapper.replace('#SBATCH --qos=standard_expand\n', '#SBATCH --qos=standard\n', 1),
              'wrapper changes exactly QOS and retains scheduler defaults')
        restored = copy.deepcopy(x)
        restored['declaration'] = old['declaration']
        for k in ['allocation12', 'RUN.sbatch']: restored['sources'][k] = old['sources'][k]
        for k in ['standard_recovery_preparer', 'standard_recovery_controls']: del restored['sources'][k]
        del restored['scheduling_recovery']
        check(restored == old, 'all other prepared inputs including handoff615sources tools libraries and roles unchanged')
        binding = P.read(x['declaration']); original_decl = P.read(d['old_declaration'])
        restored_binding = copy.deepcopy(binding)
        restored_binding['status'] = original_decl['status']; restored_binding['prospective_build']['qos'] = 'standard_expand'
        del restored_binding['resource_recovery']
        check(restored_binding == original_decl and binding['parents']['BUILD_HANDOFF.json'] == x['handoff']
              and binding['prospective_build']['archive_workers'] == 12,
              'compatible declaration retains original BUILD source_context required keys and recipe')
        cfg = json.loads((P.HERE/'CONFIG_TEMPLATE.json').read_text())
        rel = json.loads((P.HERE/'RELEASE_TEMPLATE.json').read_text())
        check(cfg['inputs'] == P.pin(prepared/'INPUTS.json') and cfg['implementation'] == rel['implementation'] == x['implementation']
              and cfg['candidate'] is None and cfg['composition'] is None and cfg['archive_workers'] == 12,
              'configuration binds actual new inputs and unchanged build without candidate')
        check(rel['status'] == 'unreleased_one_periodic_v2_build' and rel['configuration'] is None
              and rel['output'] is None and rel['build_source_review'] is None and rel['pending_old_job_resolution'] is None
              and [rel[k] for k in ['maximum_archive_builds','maximum_wrapper_compiles','maximum_wrapper_links','maximum_native_programs']] == [1,8,3,0],
              'release template remains unbound with exact original finite command limits')
        old_ready = P.read(d['old_ready']); old_controls = P.read(old_ready['controls'])
        check(old_controls['check_count'] == 30 and old_controls['status'] == 'passed_periodic_v2_build_source_controls',
              'thirty old build controls reused by pin without rerun')
        result.update(status='passed_standard_QOS_build_recovery_controls', preparation=P.pin(prepared/'RESULT.json'),
                      inputs=P.pin(prepared/'INPUTS.json'), old_controls_reused=old_ready['controls'])
    except BaseException as exc:
        result.update(status='failed_standard_QOS_build_recovery_controls', error=repr(exc))
    result.update(check_count=len(checks), cost={'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu,
        'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1})
    P.save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': P.pin(out/'RESULT.json'), 'checks': len(checks), 'error': result['error']}))
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()))
