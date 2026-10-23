#!/usr/bin/env python3
"""Independent saved-input review of the declared eight hydrated-amide SPs."""
import ast
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[4]
DIAG = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'
BRANCH = DIAG/'hydrated_amide_force_qualification_v1'
OUT = WORK/'hydrated_amide_force_independent_review_v1'
PINNED = {
    'REFERENCE.py': '7de9c17aa9198ba07c66c3f934390fc2a67d94a77860eca8ba4c0378b25ad677',
    'SPEC.md': 'f5de75cfa8b006a07e93267f6ebcdeef1d02e37a31d544872656698bc92112d4',
}
PANEL = WORK/'hydrated_amide_force_qualification_v1/prepared_v2/PANEL.json'
PANEL_SHA = '38b42c1020c47b12512076d59d66e299d4c8eb7f14754f8cfe5f71bec502e4cc'
checks = []
pins = {}

def check(value, label):
    checks.append({'check': label, 'passed': bool(value)})
    if not value:
        raise AssertionError(label)

def pin(path):
    path = Path(path)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

def verify(obj):
    if isinstance(obj, dict):
        if set(obj) == {'path', 'sha256'}:
            if obj['path'] not in pins:
                actual = pin(obj['path'])
                check(actual == obj, 'pin: '+obj['path'])
                pins[obj['path']] = actual
        else:
            for value in obj.values(): verify(value)
    elif isinstance(obj, list):
        for value in obj: verify(value)

def norm(v): return math.sqrt(sum(x*x for x in v))
def sub(a, b): return [x-y for x, y in zip(a, b)]

def canonical_input(text):
    # Only EnGrad and allocation controls can differ from the original deck.
    text = text.replace(' EnGrad', '')
    text = re.sub(r'^%maxcore \d+$', '%maxcore RESOURCE', text, flags=re.M)
    return re.sub(r'^%pal nprocs \d+ end$', '%pal nprocs RESOURCE end', text, flags=re.M)

def main():
    tick = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=False)
    for name, sha in PINNED.items(): check(pin(BRANCH/name)['sha256'] == sha, 'frozen '+name)
    check(pin(PANEL)['sha256'] == PANEL_SHA, 'frozen panel')
    p = json.loads(PANEL.read_text())
    ready = json.loads((BRANCH/'READY.json').read_text())
    verify(p); verify(ready)
    spec = importlib.util.spec_from_file_location('amide8_reviewed', BRANCH/'REFERENCE.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    def no_native(*args, **kwargs): raise AssertionError('Native invocation prohibited in independent review')
    m.subprocess.run = no_native
    check(p['maximum_native_calls'] == len(p['cells']) == 8, 'exact eight calls')
    check(p['ranks'] == 4 and p['steps_A'] == [.01, .005] and p['derivative_resolution_kcal_mol_A'] == .1, 'unchanged resource/step/resolution declaration')
    expected = {f'pbe0__{metal}__{side}_{step}__complex' for metal in ('La','Dy') for side in ('minus','plus') for step in ('010','005')}
    check({c['cell_id'] for c in p['cells']} == expected and len(set(ready['cells'])) == 8, 'exact slots and no center repeat')
    check(not p['physical_state_qualified'] and not p['method_accuracy_qualified'] and not p['validation_outcomes_opened'], 'scope limits preserved')
    metrics = []
    center_reviews = {}
    for metal in ('La', 'Dy'):
        o = p['central'][metal]
        g = o['original_cell']['geometry']
        original = Path(o['executed_input']['path']).read_text()
        saved = json.loads(Path(o['saved_result']['path']).read_text())
        raw = Path(o['output']['path']).read_text()
        ids = list(range(1,10)); u = g['radial_direction']
        check(g['ligand_indices'] == ids and g['metal_index'] == 0 and abs(norm(u)-1) < 1e-14, metal+' actual motion indices/unit direction')
        check(sum(u[j]*(g['positions_A'][1][j]-g['positions_A'][0][j]) for j in range(3)) > 0, metal+' positive motion is outward')
        derivative = -sum(saved['forces_kcal_mol_A'][i][j]*u[j] for i in ids for j in range(3))
        check(abs(derivative-o['analytic']['coordinate_derivative_kcal_mol_A']) < 1e-12, metal+' independent ligand-force projection/sign')
        strict = m.w.parse(raw, m.metadata(g), saved['returncode'])
        check(strict['status'] == 'complete' and strict['energy_hartree'] == saved['energy_hartree'], metal+' saved native central strict replay')
        check(strict == o['strict_native_replay'], metal+' exact retained native residual/stage evidence')
        check(m.parse(raw, m.metadata(g), 0)['status'] == 'failed_admission', metal+' SP adapter rejects actual EnGrad output')
        wrong = m.metadata(g); wrong['positions_A'] = [list(v) for v in g['positions_A']]; wrong['positions_A'][1][0] += .1
        check(m.w.parse(raw, wrong, 0)['status'] == 'failed_admission', metal+' wrong coordinates rejected')
        stages = strict['native_energy_stages']
        changed = json.loads(json.dumps(stages)); changed['has_analytic_gradient_stage'] = False
        comparison = m.w.compare_stages(stages, changed)
        check(comparison['expected_EnGrad_to_SP_difference'] and not comparison['changed_energy_stage_signature'] and not comparison['surface_equivalence_proven'], metal+' expected stage change does not assert surface equality')
        changed['recomputed_final_exchange'] = not changed['recomputed_final_exchange']
        check(m.w.compare_stages(stages, changed)['changed_energy_stage_signature'], metal+' actual exchange-stage change flagged')
        center_reviews[metal] = {'analytic_ligand_derivative_kcal_mol_A':derivative, 'native_checks':strict['checks'], 'printed_scf_criteria':strict['printed_scf_criteria'], 'SCF_mode':strict['native_scf_mode'], 'SCF_stop_signals':strict['native_scf_stop_signals'], 'basis_functions_by_center':strict['basis_functions_by_center']}
        for c in [x for x in p['cells'] if x['metal'] == metal]:
            lines = Path(c['source_geometry']['path']).read_text().splitlines(); entries = [line.split() for line in lines[2:]]
            x = [[float(v) for v in row[1:]] for row in entries]
            check(int(lines[0]) == len(entries) == 31 and [row[0] for row in entries] == g['symbols'] and x == c['positions_A'], c['cell_id']+' actual XYZ atom order/values')
            q = (-1 if 'minus' in c['shift_label'] else 1)*(.01 if c['shift_label'].endswith('010') else .005)
            check(c['q_A'] == q, c['cell_id']+' sign/step label')
            check(all(x[i] == g['positions_A'][i] for i in [0]+list(range(10,31))), c['cell_id']+' fixed metal/seven waters')
            err = max(abs(x[i][j]-g['positions_A'][i][j]-q*u[j]) for i in ids for j in range(3))
            check(err < 1e-12, c['cell_id']+' rigid ligand translation')
            internal = 0.
            for i in ids:
                for j in ids:
                    diff = abs(norm(sub(x[i],x[j]))-norm(sub(g['positions_A'][i],g['positions_A'][j])))
                    internal = max(internal,diff)
                    check(diff < 2e-12, c['cell_id']+f' internal distance {i}/{j}')
            deck = Path(c['template_input']['path']).read_text()
            check(canonical_input(deck) == canonical_input(original) and 'EnGrad' not in deck and '* xyzfile 3 1 geometry.xyz' in deck, c['cell_id']+' exact science deck except EnGrad/resources')
            check('%maxcore 3000\n%pal nprocs 4 end' in deck and 'NoAutostart' in deck and 'Guess PModel' in deck, c['cell_id']+' allocation template/no restart')
            check(c['charge'] == 3 and c['multiplicity'] == 1 and c['physical_multiplicity_of_parent'] == (1 if metal == 'La' else 6), c['cell_id']+' charge/effective/physical state')
            nuclear = [11]+[{'C':6,'N':7,'O':8,'H':1}[s] for s in g['symbols'][1:]]
            counts = [47]+[6 if s == 'H' else 31 for s in g['symbols'][1:]]
            check(c['expected_effective_nuclear_charges'] == nuclear and sum(nuclear)-3 == c['expected_electrons'] == 110, c['cell_id']+' independent electron count')
            check(c['expected_basis_functions_by_center'] == counts and sum(counts) == c['expected_basis_functions'] == 502, c['cell_id']+' independent AO count')
            check(not any(c['ghost_mask']) and c['expected_ecp_core'] == (46 if metal == 'La' else 55) and f'N_core {c["expected_ecp_core"]}' in deck, c['cell_id']+' real centers/ECP')
            metrics.append({'cell_id':c['cell_id'],'translation_error_A':err,'internal_distance_error_A':internal})
    for shift in ('minus_010','minus_005','plus_005','plus_010'):
        a,z = [next(c for c in p['cells'] if c['metal'] == metal and c['shift_label'] == shift) for metal in ('La','Dy')]
        check(a['positions_A'] == z['positions_A'] and a['symbols'][1:] == z['symbols'][1:], shift+' exact paired coordinates')
    # Pure analytic polynomial controls; these are not native or scientific outputs.
    rows = {}
    for metal,slope in (('La',3.),('Dy',8.)):
        for label,q in m.SHIFTS.items():
            electronic = (100+slope*q+4*q*q+5*q*q*q)/m.b.HA
            dispersion = (-2+.25*q)/m.b.HA
            rows[f'pbe0__{metal}__{label}__complex'] = {'electronic_energy_hartree':electronic,'D4_energy_hartree':dispersion,'energy_hartree':electronic+dispersion}
        for h in (.01,.005):
            d = m.derivatives(rows,metal,h)
            check(abs(d['electronic']-(slope+5*h*h)) < 1e-9 and abs(d['D4']-.25) < 1e-9 and abs(d['PBE0_D4']-(slope+.25+5*h*h)) < 1e-9, metal+f' signed FD/electronic/D4 polynomial h={h}')
    for h in (.01,.005):
        check(abs(m.derivatives(rows,'Dy',h)['PBE0_D4']-m.derivatives(rows,'La',h)['PBE0_D4']-5) < 1e-9, f'Dy-minus-La derivative sign h={h}')
    empty = m.collect(PANEL, [])
    check(empty['status'] == 'unavailable' and empty['admitted_cells'] == 0 and empty['expected_cells'] == 8 and len(empty['slots']) == 8 and 'per_metal' not in empty, 'missing cells retain eight slots/no derivative conclusion')
    source = (BRANCH/'REFERENCE.py').read_text(); ast.parse(source)
    check("two_step_resolved" in source and source.count("'sampled_analytic_consistency_at_both_steps'") == 2 and source.count("'prespecified_smaller_step_summary'") == 2, 'individual and paired summaries present')
    check("errors={h:v['PBE0_D4']-analytic for h,v in paired.items()}" in source and "changes={k:paired['0.005'][k]-paired['0.01'][k]" in source, 'paired analytic/two-step gates use paired values')
    launcher = Path(ready['launcher']['path']).read_text()
    check('#SBATCH --ntasks=4' in launcher and '#SBATCH --cpus-per-task=1' in launcher and '#SBATCH --no-requeue' in launcher and '#SBATCH --qos=standard_expand' in launcher and not re.search(r'^#SBATCH --(?:mem|time)',launcher,re.M), 'declared four-rank/default-resource/no-requeue launcher')
    result = {'status':'passed_independent_saved_review','reviewer':pin(__file__),'ready':pin(BRANCH/'READY.json'),'panel':pin(PANEL),'specification':pin(BRANCH/'SPEC.md'),'implementation':pin(BRANCH/'REFERENCE.py'),'checks':checks,'passed_checks':len(checks),'verified_unique_pins':list(pins.values()),'geometry_metrics':metrics,'central_reviews':center_reviews,'new_native_calls':0,'new_submissions':0,'local_threads':1,'elapsed_seconds':time.monotonic()-tick,'command':sys.argv,'no_launch_blocker_identified':True,'limits':['Eight native SPs remain unexecuted by this review.','Numerical derivative agreement would not establish frozen-f physical accuracy.','The original coarse-screen failure remains preserved.']}
    (OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'pins':len(pins),'elapsed_seconds':result['elapsed_seconds'],'result':pin(OUT/'RESULT.json')}))

if __name__ == '__main__': main()
