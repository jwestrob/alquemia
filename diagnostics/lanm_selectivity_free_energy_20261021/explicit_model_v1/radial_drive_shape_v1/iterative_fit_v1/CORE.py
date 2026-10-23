"""Radial14 campaign policy over frozen affine/scalar helpers; one solve per later step."""
import copy
import json
import math
import os
from pathlib import Path
import time
import numpy as np

FRACTIONS = (1., .5, .25, .125)
ELECTRICAL = ('pair_mu', 'C_mu', 'chi')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def probe_points(values, domain):
    require(len(values) == len(domain['parameters']) == 14, 'Fourteen coordinates')
    lower = np.array([p['lower'] for p in domain['parameters']]); upper = np.array([p['upper'] for p in domain['parameters']])
    center = np.asarray(values, dtype=float); width = upper-lower
    require(np.isfinite(center).all() and np.all(lower <= center) and np.all(center <= upper), 'Finite literal box')
    z = (center-lower)/width; points = []
    for j in range(14):
        for sign in (-1, 1):
            target = z.copy(); target[j] = np.clip(target[j]+sign*.005, 0, 1)
            value = center.copy(); value[j] = np.clip(center[j]+(target[j]-z[j])*width[j], lower[j], upper[j])
            if j >= 12 and center[j] == 0.:
                value[j] = sign*.01  # Preserve the admitted literal zero-origin radial stencil.
            points.append(value.tolist())
    return points


def restore_problem(saved):
    p = copy.deepcopy(saved)
    for k in ('center', 'lower', 'upper', 'width', 'z', 'bounds_lower', 'bounds_upper', 'Jacobian_original531', 'Jacobian_components'):
        p[k] = np.asarray(p[k], dtype=float)
    for b in [*p['blocks'].values(), p['energy'], *p['forces'].values(), *p['electrical_unweighted'].values()]:
        for k in ('r', 'J', 'weights'):
            if k in b: b[k] = np.asarray(b[k], dtype=float)
    p['signed_W9_z']['J_Debye'] = np.asarray(p['signed_W9_z']['J_Debye'], dtype=float)
    return p


class Policy:
    def __init__(self, core, attr, math_helper, solver_helper, fixed, inventory, targets, saved_affine, saved_problem, settings_pin):
        self.core, self.attr, self.M, self.S = core, attr, math_helper, solver_helper
        self.fixed, self.inv, self.targets = fixed, inventory, targets
        self.saved, self.initial_problem, self.settings_pin = saved_affine, saved_problem, settings_pin
        self.safety = inventory['fixed_safeguards']
        # Rebind only the isolated dimension-view's probe callback. Frozen modules remain unchanged.
        self.core.probe_points = probe_points
        self.core.linearize.__globals__['probe_points'] = probe_points
        require(saved_affine['settings'] == settings_pin, 'Exact saved solver settings')

    def assess(self, vector):
        r = self.core.assess(vector, self.fixed)
        electrical = {k:dict(score=r['original_block_scores'][k], ceiling=self.safety['actual_scores'][k],
                            passed=r['original_block_scores'][k] <= self.safety['actual_scores'][k]) for k in ELECTRICAL}
        C = vector['response']['C_dipoles_Debye']
        model = math.fsum([C['aquo9_central_Dy'][2], -C['aquo9_central_La'][2]])
        residual = math.fsum([C['aquo9_central_Dy'][2], -C['aquo9_central_La'][2], -self.safety['W9_reference_Debye']])
        signed = dict(model_Dy_minus_La_z_Debye=model, reference_Dy_minus_La_z_Debye=self.safety['W9_reference_Debye'],
                      model_minus_reference_Debye=residual, absolute_ceiling_Debye=self.safety['W9_absolute_ceiling_Debye'],
                      passed=abs(residual) <= self.safety['W9_absolute_ceiling_Debye'])
        passed = r['all25_original_physical_gates'] and all(v['passed'] for v in electrical.values()) and signed['passed']
        r.update(fixed_electrical_safeguards=electrical, signed_W9_z_safeguard=signed, all_declared_safeguards_pass=passed,
                 full_development_target=r['full_development_target'] and passed,
                 original25_physical_gate_meaning_unchanged=True, electrical_nonincrease_is_development_safeguard=True,
                 parameters=vector['parameters'], response=vector['response'])
        return r

    def problem(self, current, probes):
        p = self.M.problem(self.core, self.attr, current, probes, self.fixed, self.targets)
        for name in ELECTRICAL:
            p['blocks'][name]['denominator'] = self.safety['literal_squared_norm_denominators'][name]
        p['signed_W9_z']['absolute_ceiling_Debye'] = self.safety['W9_absolute_ceiling_Debye']
        require(p['signed_W9_z']['reference_Dy_minus_La_z_Debye'] == self.safety['W9_reference_Debye'], 'Fixed signed target')
        return p

    def fraction(self, p, point, fraction, success):
        cert = self.S.certificate(self.core, self.attr, p, point, fraction, success)
        tol = 1e-10*max(1., abs(p['raw_center_q']))
        # No 1e-8 feasibility relaxation. The inherited q-only cushion is preserved.
        eligible = bool(cert['retained_exact_sign_feasible'] and cert['all25_predicted'] and cert['signed_W9_safeguard_pass']
                        and cert['actual_component_q'] < p['raw_center_q']-tol)
        return dict(fraction=fraction, parameters=cert['parameters'], key=cert['key'], eligible=eligible,
                    prediction=cert, solver_success=success, execution_status='not_attempted')

    def propose(self, current, probes, step, out):
        require(1 <= step <= 8, 'Eight logical proposals maximum')
        require([v['parameters'] for v in probes] == probe_points(current['parameters'], self.fixed['domain']), 'Exact current28 stencil')
        base = self.core.conic; out = Path(out); out.mkdir(parents=True, exist_ok=False)
        result = dict(status='started', outer_step=step, origin_key=self.core.token(current['parameters']),
                      derivative_keys=[self.core.token(v['parameters']) for v in probes], solver_calls=0, solver_constructions=0,
                      error=None, trial_slots=[], native_calls=0, global_optimality_claim=False)
        save = lambda: base.save(out/'RESULT.json', result)
        try:
            if step == 1:
                require(result['origin_key'] == self.inv['origin_key'] and result['derivative_keys'] == self.inv['derivative_keys'], 'Exact saved first identities')
                p = restore_problem(self.initial_problem)
                require(p['center'].tolist() == current['parameters'] and p['spans'] == self.inv['actual_normalized_spans'], 'Literal saved initial problem')
                require(all(p['blocks'][k]['denominator'] == self.safety['literal_squared_norm_denominators'][k] for k in ELECTRICAL), 'Literal first electrical normalization')
                require(p['signed_W9_z']['absolute_ceiling_Debye'] == self.safety['W9_absolute_ceiling_Debye'], 'Literal first signed ceiling')
                result['imported_saved_proposal'] = self.inv['saved_affine_result']
            else:
                p = self.problem(current, probes)
            base.save(out/'LINEARIZATION.json', base.plain(p)); result['linearization'] = base.pin(out/'LINEARIZATION.json')
            model = self.S.assemble(self.core, self.attr, p, 'radial14')
            require(model['A'].shape == (1437, 15) and len(model['cones']) == 21, 'Fixed20 SOC/890 linear sides')
            np.savez(out/'CONIC_DATA.npz', A_data=model['A'].data, A_indices=model['A'].indices, A_indptr=model['A'].indptr, A_shape=model['A'].shape, b=model['b'], c=model['c'], columns=model['columns'])
            base.save(out/'ROW_LAYOUT.json', model['layout']); result.update(conic_data=base.pin(out/'CONIC_DATA.npz'), row_layout=base.pin(out/'ROW_LAYOUT.json'))
            settings = base.frozen_settings(self.settings_pin); result['effective_settings'] = base.attributes(settings); save()
            if step == 1:
                saved = np.load(base.verify(self.inv['saved_radial_conic']), allow_pickle=False)
                expected = dict(A_data=model['A'].data, A_indices=model['A'].indices, A_indptr=model['A'].indptr, A_shape=model['A'].shape, b=model['b'], c=model['c'], columns=model['columns'])
                require(all(np.array_equal(saved[k], v) for k, v in expected.items()) and model['layout'] == base.read(self.inv['saved_radial_layout']), 'Exact saved initial conic matrix/layout')
                old = self.saved['cases']['radial14']
                require(result['effective_settings'] == self.saved['effective_settings'], 'Exact initial effective settings')
                result.update({k:copy.deepcopy(old[k]) for k in ('solution', 'info', 'solver_status', 'solver_success', 'independent_residuals', 'dual_by_named_constraint')})
                result['inherited_solver_cost'] = old['cost']; result['new_solver_cost'] = dict(wall_seconds=0., CPU_seconds=0., threads=1)
            else:
                with (out/'ATTEMPT.json').open('x') as f:
                    json.dump(dict(step=step, origin=result['origin_key'], maximum_attempts=1, pid=os.getpid()), f); f.flush(); os.fsync(f.fileno())
                result.update(status='attempt_reserved', attempt=base.pin(out/'ATTEMPT.json'), solver_constructions=1); save()
                tick, cpu = time.monotonic(), time.process_time()
                solver = self.S.clarabel.DefaultSolver(model['P'], model['c'], model['A'], model['b'], model['cones'], settings)
                require(base.attributes(solver.get_settings()) == result['effective_settings'], 'All effective solver settings')
                solver.print_to_buffer(); result.update(status='solve_entered', solver_calls=1); save()
                sol = solver.solve(); info = solver.get_info(); (out/'SOLVER.log').write_text(solver.get_print_buffer())
                details = base.attributes(info); details['linsolver'] = base.attributes(info.linsolver)
                result.update(solution=base.plain(base.attributes(sol)), info=base.plain(details), solver_status=str(sol.status), solver_success=str(sol.status) == 'Solved',
                              new_solver_cost=dict(wall_seconds=time.monotonic()-tick, CPU_seconds=time.process_time()-cpu, threads=1),
                              independent_residuals=base.plain(base.diagnostics(model, sol)))
                z = np.asarray(sol.z)
                result['dual_by_named_constraint'] = {a['name']:dict(vector=z[a['start']:a['stop']].tolist(), minus_b_dot_z=float(-model['b'][a['start']:a['stop']]@z[a['start']:a['stop']])) for a in model['layout']}
            x = np.asarray(result['solution']['x'], dtype=float)
            if x.shape == (15,) and np.isfinite(x).all() and x[-1] >= 0:
                point = np.r_[x[:14], x[-1]**2]
                result['trial_slots'] = [self.fraction(p, point, f, result['solver_success']) for f in FRACTIONS]
                if step == 1:
                    for i, old in enumerate(self.saved['cases']['radial14']['predictions']):
                        require(result['trial_slots'][i]['prediction'] == old, 'Literal saved full/half certificate')
            else:
                result['prediction_rejection'] = 'Malformed/nonfinite/negative-t return; no fallback'
            result['status'] = 'completed_proposal'
        except Exception as exc:
            result.update(status='failed_proposal_no_retry', error=repr(exc)); raise
        finally:
            save()
        return result

    def select(self, current, slots):
        initial = self.assess(current); require(initial['all_declared_safeguards_pass'], 'Current safeguarded incumbent')
        require([s['fraction'] for s in slots] == list(FRACTIONS), 'All four ordered slots')
        require(all(not s['eligible'] or (s.get('execution_status') == 'complete' and s.get('native_admitted') is True) for s in slots), 'Complete eligible batch before selection')
        best, selected, reports = initial, None, []
        tol = 1e-10*max(1., abs(initial['q']))
        for i, slot in enumerate(slots):
            if not slot['eligible']:
                reports.append(dict(fraction=slot['fraction'], status='not_eligible')); continue
            v = slot['vector']; require(v['parameters'] == slot['parameters'] and self.core.token(v['parameters']) == slot['key'], 'Exact certified trial identity')
            metrics = self.assess(v); reports.append(dict(fraction=slot['fraction'], status='native_admitted', metrics=metrics))
            if metrics['all_declared_safeguards_pass'] and metrics['q'] < best['q']-tol:
                best, selected = metrics, i
        return dict(selected_trial_index=selected, selected_key=best['key'], selected_metrics=best, current_metrics=initial,
                    strict_q_improvement_tolerance=tol, trials=reports, current_wins_ties=True, earlier_fraction_wins_ties=True,
                    derivative_probes_are_candidates=False, global_convergence_claim=False)
