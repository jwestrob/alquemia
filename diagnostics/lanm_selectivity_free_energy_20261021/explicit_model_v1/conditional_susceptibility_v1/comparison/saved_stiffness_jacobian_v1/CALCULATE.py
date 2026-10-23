#!/usr/bin/env python3
"""Baseline oxygen-stiffness Jacobian from admitted saved field dipoles only."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

STATES = ('La_dimer', 'La_trimer', 'Dy_dimer', 'Dy_trimer')
PARAMETERS = ('a_La', 'b_La', 'c_La', 'a_Dy', 'b_Dy', 'c_Dy')
STEPS = {'h': .001, 'half': .0005}
TERMS = {'La': {'La_trimer': 1, 'La_dimer': -1}, 'Dy': {'Dy_trimer': 1, 'Dy_dimer': -1},
         'Dy_minus_La': {'Dy_trimer': 1, 'Dy_dimer': -1, 'La_trimer': -1, 'La_dimer': 1}}


def check(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(value):
    check(pin(value['path']) == value, 'Changed pinned input')
    return Path(value['path'])


def read(value):
    return json.loads(verify(value).read_text())


def dot(a, b):
    return math.fsum(x*y for x, y in zip(a, b))


def shape(r):
    t = 2*r
    return ((1+t+t*t/3)*math.exp(-t)/((1+5+25/3)*math.exp(-5)))**2


def metadata(row):
    lines = [s.split() for s in verify(row['input']['xyz']).read_text().splitlines()[1:] if s.strip()]
    check(len(lines) == row['atoms'] and [int(a[0]) for a in lines] == list(range(1, len(lines)+1)), 'XYZ denominator/order')
    check([int(a[5]) for a in lines] == row['types'], 'Frozen type order')
    xyz = [[float(v) for v in a[2:5]] for a in lines]; bonds = [set(int(v)-1 for v in a[6:]) for a in lines]
    records = {}; groups = {}; scales = {}; all_tokens = []
    for source in ('parameters', 'key'):
        for line in verify(row['input'][source]).read_text().splitlines():
            t = line.split('#', 1)[0].split()
            if not t: continue
            all_tokens.append(t)
            if t[0].lower() == 'polarize':
                typ = int(t[1]); records[typ] = {'alpha_A3': float(t[2]), 'source': source, 'record': line}
                groups[typ] = set(map(int, t[4:]))
            if t[0].lower() in ('polar-12-scale', 'polar-13-scale', 'polar-12-intra', 'polar-13-intra', 'direct-11-scale'):
                scales[t[0].lower()] = float(t[1])
    check(all(scales.get(k) == 0 for k in ('polar-12-scale', 'polar-13-scale', 'polar-12-intra', 'polar-13-intra', 'direct-11-scale')), 'Scoped equal direct/permanent exclusions')
    forbidden = {'polar-pair-response', 'exchpol', 'exchange-polar', 'chgpen', 'charge-flux', 'chgflux'}
    check(not any(t[0].lower() in forbidden for t in all_tokens), 'Nonbaseline response keyword')
    check(all(math.isfinite(records[t]['alpha_A3']) and records[t]['alpha_A3'] > 0 for t in row['types']), 'Positive complete alpha inventory')
    metal = [i for i, t in enumerate(row['types']) if t in (500, 501)]
    oxy = [i for i, t in enumerate(row['types']) if t == 36]
    check(metal == [0] and not bonds[0] and len(oxy) in (1, 2), 'Exact one-metal dimer/trimer inventory')
    covered = {0}; oxygen = []
    for i in oxy:
        hs = bonds[i]
        check(len(hs) == 2 and all(row['types'][j] == 37 and bonds[j] == {i} for j in hs), 'Intact complete-water graph')
        check(groups[36] == {37} and groups[37] == {36}, 'Exact water polarization groups')
        covered |= {i} | hs
        v = [x-y for x, y in zip(xyz[i], xyz[0])]; r = math.sqrt(dot(v, v))
        check(r > 0 and math.isfinite(r), 'Nonzero finite metal–oxygen distance')
        oxygen.append({'site': i, 'native_site': i+1, 'parent_atom_index': row['parent_atom_indices'][i],
                       'r_A': r, 'n': [x/r for x in v], 'f': shape(r), 'resolved_alpha': records[36]})
    check(covered == set(range(len(lines))), 'Unexpected graph/site')
    for o in oxygen:
        o['eta'] = math.fsum(p['f'] for p in oxygen if p['site'] != o['site'])
    return oxygen, records


def raw_sites(record, row):
    moments = record['moments']; verify(moments['analyze_output'])
    lines = [s.split() for s in verify(moments['induced_output']).read_text().splitlines() if s.strip()]
    check(int(lines[0][0]) == row['atoms'] and len(lines) == row['atoms']+1, 'Raw induced site denominator')
    vectors = []
    for i, t in enumerate(lines[1:]):
        check(len(t) == 5 and int(t[0]) == i+1 and t[1] == row['native_names'][i], 'Raw induced site/name order')
        v = list(map(float, t[2:])); check(all(math.isfinite(x) for x in v), 'Finite site dipoles')
        vectors.append(v)
    check(vectors == moments['induced_Debye'], 'Saved/raw site vectors differ')
    check(0 <= moments['induction']['maximum_printed_RMS_residual_Debye'] <= 1e-8, 'Original induction gate')
    return vectors


def singulars2(a, b, c):
    """Singular values of [[a,0],[b,c]], no fitting/inversion."""
    tau = a*a+b*b+c*c; det2 = (a*c)**2
    high = .5*(tau+math.sqrt(max(0., tau*tau-4*det2)))
    return [math.sqrt(high), math.sqrt(det2/high) if high else 0.]


def calculate(cfg):
    check(cfg['implementation'] == pin(__file__), 'Frozen calculation source')
    verify(cfg['declaration']); verify(cfg['design'])
    for value in cfg['native_sources'].values(): verify(value)
    model = read(cfg['model_result']); m = read(cfg['model_manifest']); review = read(cfg['model_review']); comp = read(cfg['comparison_result'])
    check(review['status'] == 'passed_independent_actual_model_susceptibility_review' and review['input_result'] == cfg['model_result'] and
          review['manifest'] == cfg['model_manifest'] and review['implementation'] == model['implementation'] and
          review['full_native_denominator'] == 19 and review['response_checks_passed'] == review['response_checks_declared'] == 36,
          'Exact model actual admission')
    verify(review['review_source']); verify(model['implementation'])
    check(model['status'] == 'conditional_model_susceptibility_complete' and model['manifest'] == cfg['model_manifest'] and model['error'] is None and
          model['controls_passed'] is True and model['summary']['numerical_checks_passed'] is True and model['native_program_invocations'] == 19,
          'Complete nineteen-call baseline')
    check(comp['status'] == 'conditional_susceptibility_comparison_complete' and comp['configuration']['model_result'] == cfg['model_result'] and
          comp['configuration']['model_review'] == cfg['model_review'] and comp['configuration']['model_manifest'] == cfg['model_manifest'], 'Admitted comparison join')
    units = m['conversion']; bohr = units['bohr_A']; debye = units['debye_per_eA']; kh = units['electric']/units['hartree_kcal_mol']
    coefficient = kh*bohr**2; g = kh/bohr
    check(units['MV_cm_per_field_au'] == units['hartree_kcal_mol']*units['elefield']/(units['electric']*bohr), 'Energy-conjugate field conversion')
    states = {}; shapes = {}
    for state in STATES:
        row = m['rows'][state]; oxygen, alphas = metadata(row); shapes[state] = oxygen; steps = {}
        for tag, field in STEPS.items():
            plus = model['rows'][state+'__plus_'+tag]; minus = model['rows'][state+'__minus_'+tag]
            check(plus['job']['field_au'] == [0., 0., field] and minus['job']['field_au'] == [0., 0., -field], 'Exact balanced field')
            for record in (plus, minus):
                check(record['status'] == 'complete' and record['job'] == m['jobs'][record['job']['id']], 'Native job identity')
            p, n = raw_sites(plus, row), raw_sites(minus, row)
            x = [[(p[i][j]-n[i][j])/(2*field*debye*bohr) for j in range(3)] for i in range(row['atoms'])]
            eps = .5e-6/(field*debye*bohr)
            chi = math.fsum(v[2] for v in x)
            check(abs(chi-model['summary']['states'][state]['steps'][tag]['chi_mu_column_a0_cubed'][2]) <= 1e-10, 'Site-sum conjugate moment closure')
            J = [0.]*6; bound = [0.]*6; offset = 0 if state.startswith('La') else 3; sites = []; sigma_upper = 0.; sigma_upper_error = 0.
            for o in oxygen:
                proj = dot(o['n'], x[o['site']]); perr = math.fsum(abs(v) for v in o['n'])*eps
                squared_error = 2*abs(proj)*perr+perr**2; weight = coefficient/o['resolved_alpha']['alpha_A3']
                db = -weight*o['f']*proj**2; dc = o['eta']*db
                eb = weight*o['f']*squared_error; ec = o['eta']*eb
                J[offset+1] += db; J[offset+2] += dc; bound[offset+1] += eb; bound[offset+2] += ec
                sigma_upper += .2*weight*proj**2; sigma_upper_error += .2*weight*squared_error
                sites.append(dict(o, response_au=x[o['site']], radial_projection=proj, radial_projection_print_bound=perr,
                                  derivative_b=db, derivative_c=dc, derivative_b_print_bound=eb, derivative_c_print_bound=ec))
            categories = {'metal': (500, 501), 'oxygen': (36,), 'hydrogen': (37,)}
            steps[tag] = {'x_site_au': x, 'site_component_print_bound_au': eps, 'chi_site_sum_a0_cubed': chi,
                          'site_category_z_response_a0_cubed': {k: math.fsum(x[i][2] for i,t in enumerate(row['types']) if t in ts) for k,ts in categories.items()},
                          'oxygen_terms': sites, 'J_a0_cubed_per_parameter': J, 'J_print_bound': bound,
                          'sigma_domain_first_order_abs_upper_bound': sigma_upper, 'sigma_domain_upper_bound_print_allowance': sigma_upper_error}
        diff = [[abs(steps['h']['x_site_au'][i][j]-steps['half']['x_site_au'][i][j]) for j in range(3)] for i in range(row['atoms'])]
        allowance = steps['h']['site_component_print_bound_au']+steps['half']['site_component_print_bound_au']
        states[state] = {'types': row['types'], 'parent_atom_indices': row['parent_atom_indices'], 'resolved_polarizabilities': alphas, 'steps': steps,
                         'max_site_response_step_difference': max(v for rowdiff in diff for v in rowdiff),
                         'step_difference_combined_print_bound': allowance,
                         'step_difference_within_print_bound': max(v for rowdiff in diff for v in rowdiff) <= allowance+1e-12}
    derived = {}; spectra = {}
    for name, terms in TERMS.items():
        derived[name] = {'terms': terms, 'steps': {}}
        for tag in STEPS:
            derived[name]['steps'][tag] = {
                'J_a0_cubed_per_parameter': [math.fsum(c*states[s]['steps'][tag]['J_a0_cubed_per_parameter'][j] for s,c in terms.items()) for j in range(6)],
                'J_print_bound': [math.fsum(abs(c)*states[s]['steps'][tag]['J_print_bound'][j] for s,c in terms.items()) for j in range(6)],
                'sigma_domain_first_order_abs_upper_bound': math.fsum(abs(c)*states[s]['steps'][tag]['sigma_domain_first_order_abs_upper_bound'] for s,c in terms.items()),
                'sigma_domain_upper_bound_print_allowance': math.fsum(abs(c)*states[s]['steps'][tag]['sigma_domain_upper_bound_print_allowance'] for s,c in terms.items())}
    for tag in STEPS:
        sv = []; sv4 = []
        for metal, offset in (('La', 0), ('Dy', 3)):
            a = states[metal+'_dimer']['steps'][tag]['J_a0_cubed_per_parameter'][offset+1]
            b,c = states[metal+'_trimer']['steps'][tag]['J_a0_cubed_per_parameter'][offset+1:offset+3]
            check(states[metal+'_dimer']['steps'][tag]['J_a0_cubed_per_parameter'][offset+2] == 0., 'Isolated-pair cooperative derivative must vanish')
            sv.extend(singulars2(a,b,c)); sv4.append(math.hypot(a,b))
        spectra[tag] = {}
        for family, values in (('six',sv), ('nested_four',sv4)):
            values = sorted(values, reverse=True); cutoff = values[0]*1e-6
            spectra[tag][family] = {'nonzero_candidate_singular_values': values, 'absolute_cutoff': cutoff,
                                    'rank_at_relative_1e_6': sum(v>cutoff for v in values), 'zero_drive_columns': True,
                                    'stiffness_column_condition': values[0]/values[-1] if values[-1]>cutoff else None}
    envelopes = []
    for r in comp['table']:
        q = r['quantity']; rec = states[q] if q in states else derived[q]; s = rec['steps'][r['step']]
        for family, cols in (('six',(1,2,4,5)), ('nested_four',(1,4))):
            box = .2*math.fsum(abs(s['J_a0_cubed_per_parameter'][j]) for j in cols)
            box_error = .2*math.fsum(s['J_print_bound'][j] for j in cols)
            upper = min(box, s['sigma_domain_first_order_abs_upper_bound'])
            printed_upper = min(box+box_error, s['sigma_domain_first_order_abs_upper_bound']+s['sigma_domain_upper_bound_print_allowance'])
            need = -r['model_minus_reference_a0_cubed']
            envelopes.append({'quantity': q, 'step': r['step'], 'estimator': r['estimator'], 'family': family,
                              'reference_minus_baseline_a0_cubed': need, 'conservative_first_order_abs_upper_bound': upper,
                              'upper_bound_including_site_print_allowance': printed_upper,
                              'target_shift_combined_print_bound': r['difference_print_bound_a0_cubed'],
                              'required_shift_exceeds_first_order_envelope': abs(need)>printed_upper+r['difference_print_bound_a0_cubed'],
                              'is_finite_parameter_feasibility_certificate': False})
    return {'status': 'saved_oxygen_stiffness_jacobian_complete', 'inputs': cfg, 'parameter_order': PARAMETERS,
            'unit_conversion': dict(units, exact_native_g=g, derivative_coefficient_with_alpha_A3=coefficient),
            'states': states, 'derived': derived, 'rank': spectra, 'linear_envelopes': envelopes,
            'assumptions': ['Frozen linear mutual response at fixed geometry; h/H are parameter corrections, not an applied field.',
                            'Ed=Ep for the exact intact-water exclusion graphs; both receive the same external field. No general dual-field extrapolation.',
                            'Six-decimal site rounding is propagated. Native solver residuals are retained in the admitted source, not converted to invented dipole error bounds.',
                            'One response column determines the baseline derivative but not the full response matrix or finite-parameter curvature. The envelopes are first order only.',
                            'Drive columns vanish for curvature in this fixed-geometry linear family; equilibrium dipole/energy constraints remain necessary for their identifiability.'],
            'coefficient_selection': None, 'optimization_performed': False, 'native_matrix_assembled': False, 'new_native_calls': 0}


def report(r):
    out = ['# Saved oxygen-stiffness Jacobian', '', 'Baseline derivatives only; no coefficient choice or finite-response feasibility claim.', '',
           '| State | Step | b derivative | c derivative | Total chi | Metal-site z contribution |', '|---|---|---:|---:|---:|---:|']
    for state, v in r['states'].items():
        offset = 0 if state.startswith('La') else 3
        for tag, s in v['steps'].items():
            J = s['J_a0_cubed_per_parameter']
            out.append(f"| {state} | {tag} | {J[offset+1]:.8g} | {J[offset+2]:.8g} | {s['chi_site_sum_a0_cubed']:.8g} | {s['site_category_z_response_a0_cubed']['metal']:.8g} |")
    out += ['', 'All derivative/rounding bounds, derived increment directions, both-family ranks and first-order domain envelopes are retained in RESULT.json.', '',
            'A needed shift outside an envelope exceeds this local linear leverage; it does not prove exclusion for finite signed stiffness. '
            'The field data contain only A⁻¹Bz, not A⁻¹ on every oxygen-radial direction or a stability certificate. '
            'Passing the old field gates does not turn print precision into a bound on solver/reference error.', '']
    return '\n'.join(out)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--config', type=Path, required=True); p.add_argument('--out', type=Path, required=True); args = p.parse_args()
    tick,cpu = time.monotonic(),time.process_time(); cfg = json.loads(args.config.read_text()); args.out.mkdir(parents=True, exist_ok=False)
    try:
        result = calculate(cfg); text = report(result); code = 0
    except Exception as exc:
        result = {'status': 'saved_oxygen_stiffness_jacobian_unavailable', 'error': repr(exc), 'new_native_calls': 0}; text = repr(exc)+'\n'; code = 1
    result['configuration'] = pin(args.config); result['cost'] = {'wall_seconds':time.monotonic()-tick, 'CPU_seconds':time.process_time()-cpu, 'local_threads':1}
    (args.out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n'); (args.out/'REPORT.md').write_text(text)
    print(json.dumps({'status':result['status'], 'result':pin(args.out/'RESULT.json')})); return code


if __name__ == '__main__':
    raise SystemExit(main())
