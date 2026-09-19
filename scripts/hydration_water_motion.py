"""Real water radial displacements: DFT energy/gradient and cheap-curvature check."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np

from affordable_common import (HA_TO_KCAL, BOHR_TO_A, InvalidArtifact, cache_key,
                               read_json, record, verify, write_new, xyz)
from affordable_response import read_engrad
from hydration_square import endpoint
from mace_hybrid import EV_TO_KCAL, Z, accepted_attempt, check_atoms, write_xyz

STAGE = 'hydration_water_motion'
PROTOCOL = 'native_r2scan3c_cpcm_anchored_mace_water_radial_motion_v1'
METHOD = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF EnGrad'
CENTERS = ('1F6S__11__Ca', '1F6S__11__La', '6IP9__110__Ca', '6IP9__110__La')
GRID = {'m050': -.05, 'm025': -.025, 'p025': .025, 'p050': .05}
TOL = {'energy_absolute_kcal_mol': .02, 'energy_relative': .25,
       'even_absolute_kcal_mol': .005, 'even_relative': .25,
       'refinement_absolute_kcal_mol': .005, 'refinement_relative': .05,
       'analytic_absolute_kcal_mol': .02, 'analytic_relative': .05,
       'gradient_absolute_kcal_mol_A': .4, 'gradient_relative': .25}


def radial_jacobian(atoms, groups):
    coords = np.array([a[1:] for a in atoms]); jac = np.zeros_like(coords)
    occupied = set()
    for water in groups:
        if water['role'] != 'variable':
            continue
        ids = water['indices']; delta = coords[water['oxygen_index']]-coords[0]
        if occupied.intersection(ids) or np.linalg.norm(delta) < 1.:
            raise InvalidArtifact('overlapping water or invalid metal-water direction')
        if sorted(atoms[i][0] for i in ids) != ['H', 'H', 'O']:
            raise InvalidArtifact('radial coordinate requires complete neutral water')
        occupied.update(ids); jac[ids] = delta/np.linalg.norm(delta)
    if not occupied:
        raise InvalidArtifact('no variable water for motion')
    return jac


def displaced(atoms, jac, h):
    return [(a[0], *(np.array(a[1:])+h*v)) for a, v in zip(atoms, jac)]


def checked_gradient(pin, result, atoms):
    g = read_engrad(verify(pin))
    if (abs(g['energy_Ha']-result['energy_hartree']) > 1e-7 or
            not np.array_equal(g['atomic_numbers'], [Z[a[0]] for a in atoms]) or
            not np.allclose(g['coordinates_bohr']*BOHR_TO_A,
                            [a[1:] for a in atoms], atol=2e-6, rtol=0)):
        raise InvalidArtifact('native gradient ordering/energy/coordinates differ')
    return g['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A


def centers_from(collection):
    c = read_json(collection); m = read_json(verify(c['manifest']))
    if c['status'] != 'complete':
        raise InvalidArtifact('completed occupancy table required')
    tasks = {t['task_id']: t for t in m['tasks']+m['reused']}
    rows = {r['task_id']: r for r in c['rows']}; centers = {}
    model = None
    for identity in CENTERS:
        t, r = tasks[identity], rows[identity]
        if verify(t['input']).read_text() != METHOD+f"\n* xyzfile {t['charge']} 1 core.xyz\n":
            raise InvalidArtifact('center DFT method differs')
        actual = endpoint(r['result']['output'], r['result']['receipt'], t['xyz'], t['input'])
        if actual != r['result']:
            raise InvalidArtifact('center energy differs from its actual output')
        atoms = xyz(verify(t['xyz'])); jac = radial_jacobian(atoms, t['groups'])
        gradient = checked_gradient(r['gradient'], actual, atoms)
        mr = read_json(verify(t['mace_result'])); mm = read_json(verify(mr['manifest']))
        mt = next(a for a in mm['tasks'] if a['task_id'] == mr['task_id'])
        if accepted_attempt(Path(t['mace_result']['path']).parent, mt, verify(mr['manifest'])) != mr:
            raise InvalidArtifact('center MACE receipt invalid')
        run = next(a for a in mr['optimization_runs'] if a['seed'] == t['selected_seed'])
        if not mr['optimization_eligible'] or xyz(verify(run['xyz'])) != atoms:
            raise InvalidArtifact('center MACE geometry/optimization differs')
        forces = np.load(verify(run['forces']), allow_pickle=False)
        if forces.shape != jac.shape or not np.isfinite(forces).all():
            raise InvalidArtifact('invalid center force array')
        if model is not None and model != mm['model']:
            raise InvalidArtifact('center model differs across endpoints')
        model = mm['model']
        centers[identity] = {'xyz': t['xyz'], 'groups': t['groups'], 'charge': t['charge'],
            'metal': t['metal'], 'case': t['case'], 'pattern': t['pattern'],
            'DFT_result': actual, 'DFT_gradient': r['gradient'], 'MACE_result': t['mace_result'],
            'MACE_manifest': mr['manifest'], 'MACE_forces': run['forces'], 'MACE_energy_eV': run['energy_eV'],
            'jacobian': jac.tolist(), 'DFT_g0_kcal_mol_A': float(np.sum(gradient*jac)),
            'MACE_g0_kcal_mol_A': -float(np.sum(forces*jac))*EV_TO_KCAL}
    for case, pattern in (('1F6S', '11'), ('6IP9', '110')):
        a, b = [centers[f'{case}__{pattern}__{z}'] for z in ('Ca', 'La')]
        if a['jacobian'] != b['jacobian'] or b['charge']-a['charge'] != 1 or a['groups'] != b['groups']:
            raise InvalidArtifact('paired physical coordinate/inventory differs')
    return m, centers, mm['software'], model


def prepare(collection, inventory, agreement, output):
    import mace_omol as omol
    source, centers, software, model = centers_from(collection)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    prep = {'protocol_id': PROTOCOL, 'centers': centers, 'source_collection': record(collection),
            'agreement': record(agreement), 'grid_A': GRID, 'tolerances': TOL,
            'coordinate': 'simultaneous rigid water translations along frozen metal-O radial vectors',
            'baseline_changed': False, 'response_model_status': 'response_model_not_validated'}
    write_new(out/'preparation.json', prep); pp = record(out/'preparation.json')
    _, mo, mm = omol.common(inventory, verify(software), agreement, out/'mace', STAGE)
    if mm['model'] != model:
        raise InvalidArtifact('changed MACE endpoint Hamiltonian')
    for name in ('hydration_water_motion.py', 'hydration_square.py'):
        path = mo/'implementation'/name; shutil.copyfile(Path(__file__).with_name(name), path)
        mm['implementation'][name] = record(path)
    do = out/'dft'; do.mkdir(); di = do/'implementation'; di.mkdir(); dpins = {}
    for name in ('hydration_water_motion.py', 'hydration_square.py', 'affordable_common.py',
                 'affordable_response.py', 'mace_hybrid.py'):
        path = di/name; shutil.copyfile(Path(__file__).with_name(name), path); dpins[name] = record(path)
    mts, dts = [], []
    for identity, c in centers.items():
        atoms = xyz(verify(c['xyz'])); jac = np.array(c['jacobian'])
        for label, h in GRID.items():
            tid = identity+'__'+label; directory = out/'geometries'/tid; directory.mkdir(parents=True)
            xp = directory/'core.xyz'; write_xyz(xp, displaced(atoms, jac, h))
            mt = {'task_id': tid, 'case_id': c['case'], 'metal': c['metal'], 'metal_index': 0,
                  'kind': 'core', 'variant': 'primary', 'charge': c['charge'], 'spin_multiplicity': 1,
                  'energy_component': omol.COMPONENT, 'energy_only': False,
                  'xyz': record(xp), 'state': check_atoms(xyz(xp), c['charge']),
                  'center_id': identity, 'amplitude_A': h, 'preparation': pp}
            mts.append(mt)
            if abs(h) != .05:
                continue
            dd = do/'tasks'/tid; dd.mkdir(parents=True)
            dx = dd/'core.xyz'; shutil.copyfile(xp, dx)
            ip = dd/'endpoint.inp'; ip.write_text(METHOD+f"\n* xyzfile {c['charge']} 1 core.xyz\n")
            dt = {'task_id': tid, 'case': c['case'], 'metal': c['metal'], 'charge': c['charge'],
                  'multiplicity': 1, 'input': record(ip), 'xyz': record(dx),
                  'output_path': str(dd/'endpoint.out'), 'engrad_path': str(dd/'endpoint.engrad'),
                  'center_id': identity, 'amplitude_A': h, 'preparation': pp}
            dt['cache_key'] = cache_key(dt | {'protocol_id': PROTOCOL}); dts.append(dt)
    mm.update(tasks=mts, preparation=pp, motion_tolerances=TOL,
              evidence_use='consumed_local_water_motion_development')
    mr = omol.seal(mo, mm)
    dm = {'protocol_id': PROTOCOL, 'tasks': dts, 'preparation': pp, 'agreement': record(agreement),
          'orca': source['orca'], 'execution_policy': source['execution_policy'],
          'execution_resources': {'mpi_ranks': 16, 'concurrent_tasks': 4}, 'implementation': dpins,
          'baseline_changed': False, 'new_endpoint_count': len(dts)}
    write_new(do/'manifest.json', dm)
    from affordable_workflow import dry_run
    return {'MACE': mr, 'DFT': dry_run(do/'manifest.json')}


def validate_mace(manifest):
    import mace_omol as omol
    m = read_json(manifest); p = read_json(verify(m['preparation']))
    if m['stage'] != STAGE or m['model'] != omol.model(verify(m['software'])) or m['motion_tolerances'] != TOL:
        raise InvalidArtifact('changed water-motion model or tolerances')
    if p['grid_A'] != GRID or p['tolerances'] != TOL or set(p['centers']) != set(CENTERS):
        raise InvalidArtifact('changed coordinate/center inventory')
    for pin in [m['agreement'], p['source_collection'], *m['implementation'].values()]:
        verify(pin)
    expected = {c+'__'+label for c in CENTERS for label in GRID}
    if len(m['tasks']) != len(expected) or {t['task_id'] for t in m['tasks']} != expected:
        raise InvalidArtifact('incomplete/duplicate water-motion task inventory')
    for t in m['tasks']:
        c = p['centers'][t['center_id']]; atoms = xyz(verify(c['xyz']))
        jac = radial_jacobian(atoms, c['groups']); label = t['task_id'].split('__')[-1]
        # Recomputed norms differ at machine precision between pinned NumPy builds.
        # Use the recorded coordinate for exact input validation; verify its
        # physical direction independently at a roundoff-only tolerance.
        actual = xyz(verify(t['xyz'])); wanted = displaced(atoms, np.array(c['jacobian']), GRID[label])
        if (not np.allclose(jac, c['jacobian'], atol=1e-14, rtol=0) or actual != wanted or t['amplitude_A'] != GRID[label] or
                t['charge'] != c['charge'] or t['spin_multiplicity'] != 1 or t['preparation'] != m['preparation'] or
                t['state'] != check_atoms(actual, c['charge'])):
            raise InvalidArtifact('changed physical displacement or electronic state')
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('water-motion cache mismatch')
    return {'status': 'pass', 'tasks': len(m['tasks']), 'manifest': record(manifest)}


def refresh_mace(manifest, output):
    """Fresh implementation snapshot after a startup failure; no scientific edits."""
    import mace_omol as omol
    validate_mace(manifest); old = read_json(manifest)
    _, out, m = omol.common(verify(old['inventory']), verify(old['software']),
                            verify(old['agreement']), output, STAGE)
    for name in ('hydration_water_motion.py', 'hydration_square.py'):
        path = out/'implementation'/name; shutil.copyfile(Path(__file__).with_name(name), path)
        m['implementation'][name] = record(path)
    if m['model'] != old['model']:
        raise InvalidArtifact('technical retry cannot change scientific model')
    m.update(tasks=[{k:v for k,v in t.items() if k!='cache_key'} for t in old['tasks']],
             preparation=old['preparation'], motion_tolerances=old['motion_tolerances'],
             evidence_use=old['evidence_use'], previous_manifest=record(manifest))
    return omol.seal(out, m)


def collect_mace(manifest):
    validate_mace(manifest); m = read_json(manifest); rows = []
    for t in m['tasks']:
        valid = [(a, accepted_attempt(a, t, manifest)) for a in
                 sorted((Path(manifest).parent/'execution'/t['task_id']).glob('attempt_*'))]
        valid = [(a, r) for a, r in valid if r is not None]
        if not valid:
            rows.append({'task_id': t['task_id'], 'status': 'unavailable'}); continue
        a, r = valid[-1]
        rows.append({'task_id': t['task_id'], 'status': 'computed', 'result': record(a/'result.json'),
                     'receipt': record(a/'receipt.json'), 'energy_eV': r['energy_eV'], 'forces': r['forces']})
    return {'manifest': record(manifest), 'rows': rows,
            'status': 'complete' if all(r['status']=='computed' for r in rows) else 'incomplete',
            'entropy_correction': None, 'response_model_status': 'response_model_not_validated'}


def collect_dft(manifest, output):
    m = read_json(manifest); p = read_json(verify(m['preparation'])); rows = []
    for t in m['tasks']:
        try:
            c = p['centers'][t['center_id']]; op = Path(t['output_path']); rp = Path(str(op)+'.execution.json')
            atoms = xyz(verify(t['xyz'])); original = xyz(verify(c['xyz']))
            if atoms != displaced(original, np.array(c['jacobian']), t['amplitude_A']):
                raise InvalidArtifact('changed DFT displacement')
            result = endpoint(record(op), record(rp), t['xyz'], t['input'])
            gp = read_json(rp)['artifacts']['engrad']; gradient = checked_gradient(gp, result, atoms)
            g = float(np.sum(gradient*np.array(c['jacobian'])))
            row = {'task_id': t['task_id'], 'status': 'complete', 'result': result,
                   'gradient': gp, 'projected_gradient_kcal_mol_A': g}
        except (OSError, ValueError, KeyError) as exc:
            row = {'task_id': t['task_id'], 'status': 'unavailable', 'failure': str(exc)}
        rows.append(row)
    result = {'manifest': record(manifest), 'rows': rows, 'baseline_changed': False,
              'status': 'complete' if all(r['status']=='complete' for r in rows) else 'incomplete'}
    write_new(output, result); return {'status': result['status'], 'completed': sum(r['status']=='complete' for r in rows)}


def compare(preparation, mace_collection, dft_collection, output):
    p = read_json(preparation); mc, dc = read_json(mace_collection), read_json(dft_collection)
    mm, dm = read_json(verify(mc['manifest'])), read_json(verify(dc['manifest']))
    if mm['preparation'] != record(preparation) or dm['preparation'] != record(preparation):
        raise InvalidArtifact('different physical preparation between methods')
    if mc != collect_mace(verify(mc['manifest'])) or dc['status'] != 'complete' or mc['status'] != 'complete':
        raise InvalidArtifact('complete actual results required')
    mr, dr = ({r['task_id']: r for r in collection['rows']} for collection in (mc, dc))
    dt = {t['task_id']: t for t in dm['tasks']}
    if set(dr) != set(dt) or len(dr) != len(dc['rows']):
        raise InvalidArtifact('missing/duplicate DFT endpoints')
    rows = []
    for identity, c in p['centers'].items():
        jac = np.array(c['jacobian']); mchanges, dchanges, mg, dg = {}, {}, {}, {}
        for label, h in GRID.items():
            tid = identity+'__'+label
            mchanges[label] = (mr[tid]['energy_eV']-c['MACE_energy_eV'])*EV_TO_KCAL
            forces = np.load(verify(mr[tid]['forces']), allow_pickle=False)
            mg[label] = -float(np.sum(forces*jac))*EV_TO_KCAL
            if abs(h) != .05:
                continue
            r, t = dr[tid], dt[tid]
            if endpoint(r['result']['output'], r['result']['receipt'], t['xyz'], t['input']) != r['result']:
                raise InvalidArtifact('DFT row differs from actual output')
            gradient = checked_gradient(r['gradient'], r['result'], xyz(verify(t['xyz'])))
            dg[label] = float(np.sum(gradient*jac))
            dchanges[label] = (r['result']['energy_hartree']-c['DFT_result']['energy_hartree'])*HA_TO_KCAL
        d_even = (dchanges['p050']+dchanges['m050'])/2
        m_even = (mchanges['p050']+mchanges['m050'])/2
        m_fine = (mchanges['p025']+mchanges['m025'])/2
        etol = max(TOL['energy_absolute_kcal_mol'], TOL['energy_relative']*abs(d_even))
        ctol = max(TOL['even_absolute_kcal_mol'], TOL['even_relative']*abs(d_even))
        rtol = max(TOL['refinement_absolute_kcal_mol'], TOL['refinement_relative']*abs(m_even))
        at = max(TOL['analytic_absolute_kcal_mol'], TOL['analytic_relative']*abs(.05*c['MACE_g0_kcal_mol_A']))
        analytic_residual = (mchanges['p050']-mchanges['m050'])/2-.05*c['MACE_g0_kcal_mol_A']
        predicted, errors, predicted_g, gradient_errors, gt = {}, {}, {}, {}, {}
        for label in ('m050', 'p050'):
            h = GRID[label]
            predicted[label] = h*c['DFT_g0_kcal_mol_A']+mchanges[label]-h*c['MACE_g0_kcal_mol_A']
            errors[label] = predicted[label]-dchanges[label]
            predicted_g[label] = c['DFT_g0_kcal_mol_A']+mg[label]-c['MACE_g0_kcal_mol_A']
            gradient_errors[label] = predicted_g[label]-dg[label]
            gt[label] = max(TOL['gradient_absolute_kcal_mol_A'], TOL['gradient_relative']*abs(dg[label]-c['DFT_g0_kcal_mol_A']))
        checks = {'anchored_energy': all(abs(e)<=etol for e in errors.values()),
                  'curvature_even': abs(m_even-d_even)<=ctol,
                  'cheap_refinement': abs(m_even-4*m_fine)<=rtol,
                  'MACE_analytic_gradient': abs(analytic_residual)<=at,
                  'anchored_gradient': all(abs(gradient_errors[label])<=gt[label] for label in gt)}
        rows.append({'center_id': identity, 'case': c['case'], 'metal': c['metal'],
                     'DFT_change_kcal_mol': dchanges, 'MACE_change_kcal_mol': mchanges,
                     'DFT_g0_kcal_mol_A': c['DFT_g0_kcal_mol_A'], 'MACE_g0_kcal_mol_A': c['MACE_g0_kcal_mol_A'],
                     'DFT_projected_gradients_kcal_mol_A': dg, 'MACE_projected_gradients_kcal_mol_A': mg,
                     'predicted_change_kcal_mol': predicted, 'energy_error_kcal_mol': errors,
                     'DFT_even_kcal_mol': d_even, 'MACE_even_kcal_mol': m_even, 'MACE_fine_even_kcal_mol': m_fine,
                     'DFT_secant_curvature_kcal_mol_A2': 2*d_even/.05**2,
                     'MACE_secant_curvature_kcal_mol_A2': 2*m_even/.05**2,
                     'energy_tolerance_kcal_mol': etol, 'even_tolerance_kcal_mol': ctol,
                     'refinement_tolerance_kcal_mol': rtol, 'MACE_analytic_residual_kcal_mol': analytic_residual,
                     'MACE_analytic_tolerance_kcal_mol': at, 'gradient_error_kcal_mol_A': gradient_errors,
                     'gradient_tolerance_kcal_mol_A': gt, 'checks': checks, 'directional_check_pass': all(checks.values())})
    contrasts = []
    for case in ('1F6S', '6IP9'):
        pair = {r['metal']: r for r in rows if r['case']==case}
        contrasts.append({'case': case, 'actual_delta_R_kcal_mol': {
            label: pair['Ca']['DFT_change_kcal_mol'][label]-pair['La']['DFT_change_kcal_mol'][label] for label in ('m050','p050')},
            'predicted_delta_R_kcal_mol': {
            label: pair['Ca']['predicted_change_kcal_mol'][label]-pair['La']['predicted_change_kcal_mol'][label] for label in ('m050','p050')}})
    result = {'protocol_id': PROTOCOL, 'preparation': record(preparation), 'MACE_collection': record(mace_collection),
              'DFT_collection': record(dft_collection), 'status': 'complete', 'rows': rows, 'contrasts': contrasts,
              'all_directional_checks_pass': all(r['directional_check_pass'] for r in rows), 'tolerances': TOL,
              'entropy_correction': None, 'relaxation_correction': None, 'occupancy_probabilities': None,
              'response_model_status': 'response_model_not_validated', 'baseline_changed': False,
              'implementation': record(__file__)}
    write_new(output, result)
    return {'status': 'complete', 'all_directional_checks_pass': result['all_directional_checks_pass'],
            'centers': [{'id': r['center_id'], 'checks': r['checks']} for r in rows]}


def report(result, output):
    """Export actual local energy comparisons; no additional scientific model."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r = read_json(result)
    if r['protocol_id'] != PROTOCOL or r['status'] != 'complete':
        raise InvalidArtifact('complete motion result required for report')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    fields = ('center_id', 'DFT_secant_curvature_kcal_mol_A2',
              'MACE_secant_curvature_kcal_mol_A2', 'directional_check_pass')
    with (out/'curvatures.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows({k:row[k] for k in fields} for row in r['rows'])
    table = ['| Center | DFT curvature | MACE curvature | Largest energy error | Passed |',
             '|---|---:|---:|---:|:---:|']
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), layout='constrained')
    for ax, row in zip(axes.flat, r['rows']):
        labels = ('m050', 'p050'); hs = [-.05, 0., .05]
        energies = [row['DFT_change_kcal_mol']['m050'], 0., row['DFT_change_kcal_mol']['p050']]
        model_h = [GRID[k] for k in GRID]
        model_e = [h*row['DFT_g0_kcal_mol_A']+row['MACE_change_kcal_mol'][k]
                   -h*row['MACE_g0_kcal_mol_A'] for k,h in GRID.items()]
        # The reused center is an exact anchor, not a straight segment between
        # the two nearest displaced samples.
        model_h.insert(2, 0.); model_e.insert(2, 0.)
        ax.plot(model_h, model_e, 'x--', color='#c34c22', label='DFT-anchored MACE')
        ax.scatter(hs, energies, color='#224b86', label='Native DFT', zorder=3)
        ax.axhline(0, color='0.8', linewidth=.5); ax.axvline(0, color='0.8', linewidth=.5)
        ax.set(title=row['center_id'], xlabel='Collective outward water motion (Å)',
               ylabel='Energy change (kcal/mol)')
        ax.ticklabel_format(axis='y', useOffset=False)
        error = max(abs(row['energy_error_kcal_mol'][k]) for k in labels)
        table.append(f"| {row['center_id']} | {row['DFT_secant_curvature_kcal_mol_A2']:.3f} | "
                     f"{row['MACE_secant_curvature_kcal_mol_A2']:.3f} | {error:.5f} | "
                     f"{'yes' if row['directional_check_pass'] else 'no'} |")
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('Local bound-water motion: fixed inventory and radial coordinate')
    fig.savefig(out/'energy_changes.svg'); fig.savefig(out/'energy_changes.png', dpi=170)
    plt.close(fig)
    (out/'TABLE.md').write_text('\n'.join(table)+'\n\nCurvature units: kcal/mol/Å²; energy errors: kcal/mol.\n'
        'Lines guide the eye between actual MACE evaluations. One direction does not establish basin entropy or occupancy.\n')
    write_new(out/'provenance.json', {'result': record(result), 'implementation': record(__file__)})
    return {'status': 'exported', 'directory': str(out)}


if __name__=='__main__':
    import json
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare')
    for name in ('collection', 'inventory', 'agreement', 'output'):q.add_argument('--'+name, required=True)
    q = sub.add_parser('collect-dft')
    for name in ('manifest', 'output'):q.add_argument('--'+name, required=True)
    q = sub.add_parser('refresh-mace')
    for name in ('manifest', 'output'):q.add_argument('--'+name, required=True)
    q = sub.add_parser('compare')
    for name in ('preparation', 'mace-collection', 'dft-collection', 'output'):q.add_argument('--'+name, required=True)
    q = sub.add_parser('report')
    for name in ('result', 'output'):q.add_argument('--'+name, required=True)
    a = vars(p.parse_args()); op = a.pop('op').replace('-', '_'); print(json.dumps(globals()[op](**a), indent=2))
