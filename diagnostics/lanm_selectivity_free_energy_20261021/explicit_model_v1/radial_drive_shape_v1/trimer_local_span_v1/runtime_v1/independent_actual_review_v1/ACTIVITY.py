"""Saved-only dual/margin addendum; no matrix reconstruction or solver/native calls."""
import hashlib
import json
import math
from pathlib import Path
import time

wall, cpu = time.perf_counter(), time.process_time()
ROOT = Path(__file__).resolve().parents[7]
REL = Path('explicit_model_v1/radial_drive_shape_v1/trimer_local_span_v1/runtime_v1')
D = ROOT / 'diagnostics/lanm_selectivity_free_energy_20261021'
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021'
RUN = W / REL / 'solve_v1'
OUT = W / REL / 'independent_actual_review_v1/activity_v1'
OUT.mkdir(exist_ok=False)


def pin(path, expected=None):
    path = Path(path).resolve()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if expected is not None:
        assert digest == expected, str(path)
    return {'path': str(path), 'sha256': digest}


def read(path, expected=None):
    authority = pin(path, expected)
    return json.loads(Path(path).read_text()), authority


result, result_pin = read(RUN / 'execution/RESULT.json', 'a320ef808f64d2c5e593b7217c718e3bc26297280d501a6c2875f7c127c4756b')
layout, layout_pin = read(RUN / 'execution/ROW_LAYOUT.json', '838319fc726d7de9c999aaa29f6d8f53e3163ac02365ad8ed46741b49442c984')
review, review_pin = read(W / REL / 'independent_actual_review_v1/actual_v1/REVIEW.json', '29e07ff187a22b994182cc86c67fde37eb7bb6b7865ee086ef6b66bbf05ac5e0')
root_review, root_pin = read(RUN / 'ROOT_ACTUAL_REVIEW.json')
assert root_pin['sha256'].startswith('a4c35d31')
domain, domain_pin = read(D / 'explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/PARAMETER_DOMAIN.json', 'a9e53eb59749b45907343e72fe85290e13402c9a20326240861a1028b927799c')
math_pin = pin(D / 'explicit_model_v1/radial_drive_shape_v1/local_comparison_v2/MATH.py', '3ea3998f97a4cb6674863e4f0c27823fe3d465d7185529cd4d12d6d5a5381fbb')
parameters = domain['parameters'] + [dict(name=m + '_radial_drive', lower=-1., upper=1., units='dimensionless') for m in ('La', 'Dy')]
assert len(parameters) == 14
witness = result['witness']
assert result['solver_status'] == 'AlmostSolved' and not witness['exact_sign_feasible']
assert result['solver_calls'] == 1 and result['native_calls'] == 0

parameter_rows = []
for i, (definition, value, step) in enumerate(zip(parameters, witness['parameters'], witness['normalized_displacement'])):
    parameter_rows.append(dict(index=i, name=definition['name'], units=definition['units'],
        value=value, normalized_displacement=step, lower=definition['lower'], upper=definition['upper'],
        parameter_lower_distance=value-definition['lower'], parameter_upper_distance=definition['upper']-value,
        normalized_trust_lower_distance=.025+step, normalized_trust_upper_distance=.025-step))

z = result['solution']['z']
margins = result['independent_residuals']['primal_affine_cone_margins']
rows = []
for item in layout:
    name = item['name']
    row = dict(name=name, kind=item['kind'], start=item['start'], stop=item['stop'],
               dual_head=z[item['start']], saved_primal_cone_margin=margins[name])
    if name.startswith(('box/', 'trust/')):
        index = int(name.rsplit('/', 1)[1])
        row.update(parameter_index=index, parameter_name=parameters[index]['name'])
    rows.append(row)
ranked = sorted(rows, key=lambda row: row['dual_head'], reverse=True)
soc = [row for row in rows if row['kind'] == 'SOC']
trimer = [row for row in rows if row['name'].startswith('trimer/')]
selected_names = {row['name'] for row in ranked[:17]} | {'block/C_mu', 'block/chi', 'signed_W9_z/lower'}
compact = [row for row in ranked if row['name'] in selected_names]
res = result['independent_residuals']
diagnostics = {k: res[k] for k in ('primal_residual_inf', 'dual_residual_inf', 'primal_objective_t', 'dual_objective_t', 'signed_gap_t', 'complementarity_dot', 'rigorously_rounded_dual_certificate')}
doc = dict(status='saved_dual_margin_interpretation_only', implementation=pin(__file__),
    authorities=dict(result=result_pin, row_layout=layout_pin, independent_actual_review=review_pin,
                     root_actual_review=root_pin, original_domain=domain_pin, radial_domain_extension=math_pin),
    scope='Read saved dual coordinates, saved cone margins and exact domain metadata; no matrix, Jacobian, molecular or optimization replay.',
    solver_status=result['solver_status'], solver_success=result['solver_success'],
    t=witness['fraction'], q=witness['q'], literal_q_ceiling=witness['literal_q_ceiling'],
    exact_sign_feasible=witness['exact_sign_feasible'], all25_predicted=witness['all25_predicted'],
    negative_scalar_margins={k:v for k,v in witness['margins'].items() if v < 0},
    exact_parameter_mapping=parameter_rows, compact_activity=compact, all_SOC_activity=soc,
    trimer_dual_head_sum=math.fsum(row['dual_head'] for row in trimer),
    dual_convention='SOC head is the dual coordinate paired with its cone radius; scalar head is its single-side multiplier. These are in the saved normalized row units, not directly comparable physical-unit sensitivities.',
    activity_convention='Largest saved dual heads with their unaltered primal cone margins; no new activity/feasibility tolerance or sign clipping.',
    numerical_KKT=diagnostics, rigorous_global_bound=False, native_candidate_admission=False,
    interpretation='Returned point is crowded jointly by P/B caps, pair dipole norm, force rows, Dy_v lower box face, Dy_epsilon upper local trust face, signed W9 upper side and all six trimer groups. Collective dipole and susceptibility have substantial slack. No unique physical cause or rigorous optimum follows.',
    finite_test_closed=True, further_solve_or_fraction_sweep=False,
    cost=dict(wall_seconds=time.perf_counter()-wall, CPU_seconds=time.process_time()-cpu, threads=1,
              native_calls=0, solver_calls=0, matrix_replays=0))
(OUT/'ACTIVE_CONSTRAINTS.json').write_text(json.dumps(doc, indent=2, allow_nan=False)+'\n')
lines = [
    'Saved approximate constraint activity — finite LS test closed', '',
    'The AlmostSolved point proposes t=0.007435177318942966 (~0.744% common reduction),',
    'but violates the exact P/B q cap by 5.5569909784836113e-11. No positive feasible',
    'witness is admitted. All 25 predicted physical gates pass; this does not remove',
    'the separate q failure. No fraction sweep, solver retry or native trial follows.', '',
    'Largest saved dual heads and unchanged cone margins (normalized row units):',
    'constraint | dual head | primal cone margin']
for row in compact:
    name=row['name']
    if 'parameter_name' in row:
        name += ' ['+row['parameter_name']+']'
    lines.append(f"{name} | {row['dual_head']:.12g} | {row['saved_primal_cone_margin']:.12g}")
lines += ['',
    'This is a joint numerical boundary: both component caps, pair-dipole preservation,',
    'several force limits, the Dy_v lower domain face and the signed-W9 ceiling all',
    'have small margins and appreciable dual heads. Dy_epsilon reaches the upper',
    'local trust face (+0.02499998021895 of +0.025), not its parameter-domain upper',
    'bound. Dy_v=-0.9999999973905632 approaches its declared -1 domain boundary.',
    'C_mu and chi cone margins are 0.08042 and 0.12118 with dual heads ~2e-9;',
    'they do not look locally limiting at this returned point. All six new trimer',
    'group cones have small margins; their dual heads distribute the common-reduction',
    'objective across the groups. This does not establish a unique physical cause.', '',
    'SOC heads describe cone-radius sensitivities in the fixed normalized formulation.',
    'Their sizes are scaling-dependent; they are not comparable physical-unit costs.',
    'The raw minimization objective is -t. The saved primal-dual objective gap is',
    '3.538916917853255e-7; independent unscaled primal/dual infinity residuals are',
    '7.01126157132137e-9 / 1.1233623609108179e-7. These are numerical KKT diagnostics,',
    'not a rigorously feasible dual bound or proof of global model inadequacy.', '',
    'Exact zero-based parameter mapping; displacement is normalized by domain width:',
    'index | parameter | returned value | normalized displacement | physical domain']
for row in parameter_rows:
    lines.append(f"{row['index']} | {row['name']} | {row['value']:.17g} | {row['normalized_displacement']:.17g} | [{row['lower']:.17g}, {row['upper']:.17g}]")
lines += ['',
    'La_u/La_v and Dy_u/Dy_v are the declared transformed box coordinates, not native',
    'c/d coefficients: c=ln(4)*(u+v)/4, d=ln(4)*(u-v)/4. R is in angstrom;',
    'epsilon is in kcal/mol; the response and radial-drive coordinates are dimensionless.', '',
    'The original arithmetic review and root admission remain unchanged. Exact pins',
    'and this saved-only read cost are in ACTIVE_CONSTRAINTS.json. Zero new native',
    'calls, solver calls, matrix/Jacobian replay, or further diagnostics.']
(OUT/'REPORT.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'result':pin(OUT/'ACTIVE_CONSTRAINTS.json'),'report':pin(OUT/'REPORT.txt'),'cost':doc['cost']},indent=2))
