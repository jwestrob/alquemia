#!/usr/bin/env python3
"""Metadata/AST receipt only; imports no producer and executes no scientific code."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parent

def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

def main():
    wall, cpu = time.monotonic(), time.process_time()
    checks, verified = [], {}
    def check(name, value):
        if not value:
            raise ValueError(name)
        checks.append({'name': name, 'passed': True})
    def verify(value):
        if isinstance(value, dict):
            if set(value) == {'path', 'sha256'}:
                if value['path'] not in verified:
                    check('pin:' + value['path'], pin(value['path']) == value)
                    verified[value['path']] = value
            else:
                for x in value.values(): verify(x)
        elif isinstance(value, list):
            for x in value: verify(x)
    def read(p):
        verify(p)
        return json.loads(Path(p['path']).read_text())
    ready = json.loads((BASE/'READY.json').read_text())
    check('frozen_READY', pin(BASE/'READY.json')['sha256'] == 'dc99de132cf0adc711140700c5492aa6a9fbc907ed22da8ddab2f4e6b63cdb3b')
    verify(ready)
    inputs = read(ready['inputs']); verify(inputs)
    cfg = read(ready['configuration_template'])
    controls = read(ready['controls']); prep = read(ready['preparation'])
    check('configuration_source_bindings', all(cfg[k] == ready[k] for k in ('implementation','field_adapter','inputs','declaration','limits')))
    check('configuration_unreleased', cfg['root_release'] is None and cfg['source_review'] is None and cfg['output'] is None)
    check('controls_bind_final_source', controls['implementation'] == ready['controls_source'] and controls['runner'] == ready['implementation'] and controls['field_adapter'] == ready['field_adapter'])
    check('saved_30_controls_pass', controls['status'] == 'passed_radial_field_recovery_controls' and controls['count'] == len(controls['checks']) == 30 and all(c['passed'] for c in controls['checks']))
    check('controls_zero_scientific_calls', all(controls[k] == 0 for k in ('native_calls','scientific_matrix_qualifications','optimizer_calls')) and controls['cost']['numerical_threads'] == 1)
    check('preflight_final_sources', prep['status'] == 'passed_radial_recovery_zero_call_preflight' and prep['implementation'] == ready['implementation'] and prep['inputs'] == ready['inputs'] and prep['configuration'] == ready['configuration_template'])
    check('preflight_scope', prep['reused_complete_rows'] == 42 and prep['saved_field_outputs'] == 4 and prep['untouched_native_slots'] == 3 and prep['new_native_calls'] == prep['scientific_matrix_qualifications'] == 0)
    check('saved_eight_header_checks', len(prep['saved_export_headers']) == 8 and all(x['maximum_error'] == 0.0 and x['native_external_tokens'][:2] == ['EXTERNAL','T'] for x in prep['saved_export_headers']))
    source = Path(ready['field_adapter']['path']).read_text()
    tree = ast.parse(source)
    guard = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'ZERO_GUARD' for t in n.targets))
    old_tree = ast.parse(Path(inputs['matrix_helper']['path']).read_text())
    old_fn = next(n for n in old_tree.body if isinstance(n, ast.FunctionDef) and n.name == 'parse_export')
    dump = lambda n: ast.dump(n, include_attributes=False)
    expected = ast.parse(guard).body[0]
    locations = [i for i, n in enumerate(old_fn.body) if dump(n) == dump(expected)]
    check('exact_single_inherited_zero_predicate', len(locations) == 1)
    scoped = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'scoped')
    text = ast.unparse(scoped)
    check('only_one_parser_statement_replacement', text.count('fn.body[hits[0]] =') == 1 and "ast.parse('check_external(lines[4],row,units)').body[0]" in text)
    replacement = copy.deepcopy(old_fn)
    replacement.body[locations[0]] = ast.parse('check_external(lines[4],row,units)').body[0]
    replacement.body[locations[0]] = copy.deepcopy(expected)
    check('inverse_parser_AST_exact', dump(replacement) == dump(old_fn))
    check('qualifier_code_object_inherited', 'types.FunctionType(matrix.qualify_exports.__code__, ns, matrix.qualify_exports.__name__, matrix.qualify_exports.__defaults__, matrix.qualify_exports.__closure__)' in text)
    parent = read(inputs['parent_result']); pc = read(inputs['parent_configuration'])
    check('exact_failed_parent_configuration', parent['configuration'] == pc and parent['error'] == "ValueError('Zero-field export required')" and parent['status'] == 'radial_drive_numerical_qualification_failed')
    check('parent_counts_46_42_132', parent['native_program_invocation_attempts'] == parent['native_programs_returned_zero'] == 46 and parent['native_programs_completed'] == 42 and parent['internal_FD_energies_observed'] == 132)
    check('parent_row_status_denominator', len(parent['rows']) == 49 and sum(r['status'] == 'complete' for r in parent['rows'].values()) == 42 and sum(r['status'] == 'failed' for r in parent['rows'].values()) == 4 and sum(r['status'] == 'not_attempted' for r in parent['rows'].values()) == 3)
    check('limits_unchanged_finite_recovery', ready['limits'] == {'new_native_programs':3,'combined_native_programs':49,'saved_complete_rows':42,'saved_field_reparses':4,'reused_internal_FD_energies':132,'new_internal_FD_energies':0,'native_workers':1,'allocated_CPUs':12,'program_timeout_seconds':120})
    wrapper = Path(ready['wrapper']['path']).read_text()
    check('scheduler_policy', all(s in wrapper for s in ('--partition=standard\n','--qos=standard_expand\n','--account=banfield\n','--ntasks=1\n','--cpus-per-task=12\n','--no-requeue\n')) and '--mem' not in wrapper and '--time' not in wrapper)
    review = dict(ready['source_review_contract'])
    review.update(reviewer_implementation=pin(__file__), ready=pin(BASE/'READY.json'), controls=ready['controls'], preparation=ready['preparation'], parent_result=inputs['parent_result'], original_failure_diagnosis=inputs['diagnosis'], checks=checks, unique_pins=len(verified), saved_author_controls_read=30, controls_rerun=False, source_findings=[
        'Only the exact zero-field parse predicate changes. Actual T flag, field sign, zero transverse components, exact prepared key, and native key/elefield conversion are required; raw exports are never rewritten.',
        'All remaining parser statements and the entire numerical qualifier code object are inherited. The original zero-field parser is unchanged for the three bypasses.',
        'Source preserves both E/M solves and dual RHS/solution/energy/residual gates. Native external-induced energy is already in the E-stage polarization component; no added energy correction or relaxed threshold is introduced.',
        'The 42 completed row objects and existing directional/operator checks are reused. Four saved field parses and both original field closures must pass before any new native call.',
        'Only the three previously untouched bypass slots can execute, serially and at most once. Intent precedes invocation; first error stops, preserves receipts and untouched remainder. Fresh output is required.',
        'Original failed result and its 288 allocated CPU-seconds remain separate from new process and scheduler costs. No sum of overlapping historical per-child deltas is used.'
    ], limitations=[
        'Source approval does not assert that the four field numerical qualifications or final three bypass programs will pass.',
        'The earlier source approval missed the inherited zero-field-only parser precondition; its receipt and the original failure remain immutable.',
        'This review does not reparse or numerically re-admit the 42 completed rows, rerun author controls, perform dense solves, or admit physical accuracy/cache equivalence.'
    ], native_calls=0, scientific_matrix_qualifications=0, optimizer_calls=0,
    cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
    (HERE/'REVIEW.json').write_text(json.dumps(review, indent=2, allow_nan=False)+'\n')
    (HERE/'REPORT.txt').write_text('PASS: narrow radial field-parser recovery source.\n\nExact nonzero field metadata is validated; all remaining parser and numerical gates remain inherited. Four saved field admissions and both field closures precede three serial untouched bypasses. The 42 completed rows, 132 FD energies, original failure and costs are preserved.\n\nRead final source and 30 saved author controls; performed metadata/AST joins only, with zero native, dense, optimizer, or control-suite execution. Actual recovery and conditional cache-bridge admission remain pending.\n')
    print(json.dumps({'review':pin(HERE/'REVIEW.json'),'report':pin(HERE/'REPORT.txt'),'checks':len(checks),'cost':review['cost']}))

if __name__ == '__main__':
    if (HERE/'REVIEW.json').exists(): raise SystemExit('Existing immutable review')
    try: main()
    except Exception as exc:
        (HERE/'FAILURE.json').write_text(json.dumps({'error':repr(exc),'native_calls':0,'scientific_matrix_qualifications':0})+'\n')
        raise
