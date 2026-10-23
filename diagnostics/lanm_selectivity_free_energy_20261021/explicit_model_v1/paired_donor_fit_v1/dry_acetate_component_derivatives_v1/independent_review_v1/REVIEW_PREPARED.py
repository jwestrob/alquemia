#!/usr/bin/env python3
"""Independent saved input joins; no native invocation or author-suite rerun."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = next(p for p in HERE.parents if (p / 'diagnostics').is_dir() and (p / 'workspaces').is_dir())
OUT = ROOT / 'workspaces' / BASE.relative_to(ROOT / 'diagnostics') / 'independent_review_v1'
checks = {}
tick, cpu = time.monotonic(), time.process_time()

def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def check(name, value):
    checks[name] = bool(value)
    assert value, name

def verify(p, name):
    check(name, pin(p['path']) == p)
    return Path(p['path'])

def read(p, name):
    return json.loads(verify(p, name).read_text())

ready = json.loads((BASE / 'READY.json').read_text())
for k, v in ready.items():
    if isinstance(v, dict) and set(v) == {'path', 'sha256'}:
        verify(v, 'ready/' + k)
m = read(ready['manifest'], 'manifest')
cfg = read(ready['inputs'], 'inputs')
for k, v in cfg.items():
    if k != 'selected_fit_result' and isinstance(v, dict) and set(v) == {'path', 'sha256'}:
        verify(v, 'input/' + k)
old = read(cfg['selected_review'], 'selected_review')
check('selected_completed_authority', old['status'] == 'passed' and old['error'] is None and old['result'] == cfg['selected_fit_result'])
check('selected_eight_parameters', [x['selected'] for x in old['parameter_bounds']] == m['selected_parameters'] and len(m['selected_parameters']) == 8)
p = read(cfg['reference_panel'], 'reference_panel')
check('exact16_scope', set(m['rows']) == {c['cell_id'] for c in p['cells']} and len(m['rows']) == 16 and m['maximum_native_calls'] == 16 and m['workers'] == 1 and m['new_native_calls'] == 0)
check('two_steps_resolution', p['steps_A'] == [.01, .005] and p['derivative_resolution_kcal_mol_A'] == .1)
check('source_lineage', m['implementation'] == ready['implementation'] and m['inputs'] == ready['inputs'] and m['analyze'] == cfg['analyze'] and m['selected_review'] == cfg['selected_review'])
labels = ('Bond Stretching', 'Angle Bending', 'Stretch-Bend', 'Out-of-Plane Bend', 'Torsional Angle', 'Van der Waals', 'Atomic Multipoles', 'Polarization')
number = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?'
def dec(s): return Decimal(s.replace('D', 'E').replace('d', 'e'))
def dot(a, b): return sum(x*y for x,y in zip(a,b))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
central = {}
for metal in ('La', 'Dy'):
    c = m['central'][metal]; src = c['source']; tag = 'central/' + metal
    check(tag + '/selected_source', src == old['selected_sources'][c['geometry_id']])
    authority = read(src['authority'], tag + '/authority')
    check(tag + '/job17621', authority['job_id'] == '17621' and authority['status'] == 'complete' and not src['reused_baseline'])
    check(tag + '/effective_inputs', authority['input'] == src['effective_input'])
    for k,v in src['effective_input'].items(): verify(v, tag + '/input/' + k)
    for k,v in src['native_outputs'].items(): verify(v, tag + '/output/' + k)
    typed = [x.split() for x in Path(src['effective_input']['xyz']['path']).read_text().splitlines()[1:] if x.strip()]
    check(tag + '/original_typed_rows', typed == c['typed_rows'] and len(typed) == 8)
    check(tag + '/geometry', c['geometry'] == p['central'][metal]['original_cell']['geometry'])
    raw = Path(src['native_outputs']['analyze_output']['path']).read_text()
    tokens = re.findall(r'^\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+number+r')\s+(\d+)\s*$', raw.split('Energy Component Breakdown :')[1], re.M)
    energies = {name.strip(): dec(value) for name,value,n in tokens}
    total = dec(re.findall(r'Total Potential Energy\s*:\s*('+number+')', raw)[0])
    check(tag + '/components', len(tokens) == 8 and set(energies) == set(labels))
    check(tag + '/saved_tokens', all(energies[k] == dec(c['energy']['tokens'][k]) for k in labels) and total == dec(c['energy']['tokens']['Total']))
    remainder = total - sum(energies.values(), Decimal(0))
    bound = sum((Decimal(5).scaleb(x.as_tuple().exponent-1) for x in [total,*energies.values()]), Decimal(0))
    check(tag + '/decimal_closure', remainder == dec(c['energy']['component_sum_remainder_kcal_mol']) and bound == dec(c['energy']['component_sum_rounding_bound_kcal_mol']) and abs(remainder) <= bound)
    gradients = re.findall(r'^\s*Anlyt\s+(\d+)\s+('+number+')\s+('+number+')\s+('+number+')', Path(src['native_outputs']['gradient_output']['path']).read_text(), re.M)
    check(tag + '/gradient_inventory', [int(v[0]) for v in gradients] == list(range(1,9)))
    grad = [[float(dec(v)) for v in row[1:]] for row in gradients]
    check(tag + '/saved_force_sign', [[-v for v in row] for row in grad] == authority['forces'])
    frame = p['central'][metal]['analytic']['frame']; xyz = c['geometry']['positions_A']; derivatives = {}
    for direction in ('radial','rotation'):
        tangents = [frame['radial_unit'] if direction == 'radial' else [v/frame['rotation_length_A'] for v in cross(frame['rotation_unit'], [xyz[i][k]-frame['centroid_A'][k] for k in range(3)])] for i in frame['ligand_indices']]
        value = sum(dot(grad[i], t) for i,t in zip(frame['ligand_indices'], tangents))
        rounding = sum(abs(t[k])*float(Decimal(5).scaleb(dec(gradients[i][k+1]).as_tuple().exponent-1)) for i,t in zip(frame['ligand_indices'], tangents) for k in range(3))
        check(tag + '/derivative/' + direction, abs(value-c['analytic'][direction]) <= 1e-12)
        check(tag + '/gradient_print_bound/' + direction, abs(rounding-c['analytic_gradient_print_bound_kcal_mol_A'][direction]) <= 1e-18)
        derivatives[direction] = value
    central[metal] = dict(energy_kcal_mol=float(total), derivatives_kcal_mol_A=derivatives, component_rounding_bound_kcal_mol=str(bound))
for cell in p['cells']:
    cid = cell['cell_id']; row = m['rows'][cid]; c = m['central'][row['metal']]
    check(cid + '/metadata', all(row[k] == cell[k] for k in ('cell_id','metal','direction','shift_label','q_A','source_geometry')))
    ordinary = [x.split() for x in verify(cell['source_geometry'], cid+'/geometry_pin').read_text().splitlines()[2:] if x.strip()]
    typed = [x.split() for x in verify(row['input']['xyz'], cid+'/xyz_pin').read_text().splitlines()[1:] if x.strip()]
    check(cid + '/literal_coordinates', len(typed) == len(ordinary) == 8 and [r[2:5] for r in typed] == [r[1:4] for r in ordinary])
    check(cid + '/types_bonds_indices', [[*r[:2],*r[5:]] for r in typed] == [[*r[:2],*r[5:]] for r in c['typed_rows']])
    for key in ('key','parameters'):
        verify(row['input'][key], cid+'/'+key+'_pin')
        check(cid+'/'+key+'_unchanged', row['input'][key]['sha256'] == c['source']['effective_input'][key]['sha256'])
controls = read(ready['saved_controls'], 'controls')
check('author_controls_join', controls['implementation'] == ready['implementation'] and controls['manifest'] == ready['manifest'] and controls['check_count'] == 24 and all(controls['checks'].values()) and controls['actual_new_native_calls'] == 0)
prep = read(ready['preparation_result'], 'preparation')
check('preparation_join', prep['implementation'] == ready['implementation'] and prep['manifest'] == ready['manifest'] and prep['cells'] == 16 and prep['saved_central_raw_replays'] == 2 and prep['new_native_calls'] == 0)
result = dict(status='passed_independent_prepared_input_checks', implementation=ready['implementation'], manifest=ready['manifest'], reviewer_source=pin(__file__), checks=checks, check_count=len(checks), central_raw_observations=central, author_controls_rerun=False, new_native_calls=0, cost=dict(wall_seconds=time.monotonic()-tick, process_CPU_seconds=time.process_time()-cpu))
OUT.mkdir(parents=True, exist_ok=True)
(OUT/'CHECKS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'], checks=len(checks), receipt=pin(OUT/'CHECKS.json'), central=central),indent=2))
