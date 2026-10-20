"""Full physical Cartesian L-BFGS and topology guards; no molecular engine."""
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

SETTINGS = dict(history=8, armijo=1e-4, backtracks=24, initial_trust_A=.025,
                minimum_trust_A=1e-6, maximum_trust_A=.10, rank_relative=1e-10,
                maximum_force=.10, RMS_atom_force=.03,
                bond_lower_factor=.65, bond_upper_factor=1.35,
                departure_distance_A=4.5)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def force_metrics(gradient):
    norms = np.linalg.norm(gradient, axis=1)
    return {'maximum_atom_kcal_mol_A': float(norms.max()),
            'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(norms**2))),
            'L2_kcal_mol_A': float(np.linalg.norm(gradient)),
            'converged': bool(norms.max() <= SETTINGS['maximum_force'] and
                              np.sqrt(np.mean(norms**2)) <= SETTINGS['RMS_atom_force'])}


def remove_rigid(x, vector):
    """Euclidean projection; no atomic force clipping or fixed physical atom."""
    center = x - x.mean(axis=0)
    columns = [np.tile(axis, (len(x), 1)).ravel() for axis in np.eye(3)]
    columns += [np.cross(axis, center).ravel() for axis in np.eye(3)]
    B = np.array(columns).T
    scales = np.linalg.norm(B, axis=0)
    require(np.all(scales > 1e-12), 'Degenerate global rigid-motion basis')
    B /= scales
    values, vectors = np.linalg.eigh(B.T @ B)
    keep = values > SETTINGS['rank_relative'] * values[-1]
    inverse = (vectors[:, keep] / values[keep]) @ vectors[:, keep].T
    flat = np.asarray(vector).ravel()
    removed = B @ (inverse @ (B.T @ flat))
    return (flat - removed).reshape(x.shape), {'rigid_rank': int(keep.sum()),
        'removed_L2': float(np.linalg.norm(removed)), 'normalized_Gram_eigenvalues': values.tolist()}


class PhysicalGuard:
    """Finite connectivity/stereochemistry guards, with no source-motion limit.

    equilibrium_A must come from the pinned native covalent force field.
    target_specs: index,label,site_donor_indices,coordinating_water_indices.
    Waters are specified by their physical oxygen indices. Labels are not forces.
    """
    def __init__(self, atoms, bond_pairs, equilibrium_A, source_x, geometry_module,
                 protein_indices, target_specs):
        self.atoms = atoms; self.ids = [a['id'] for a in atoms]
        self.source = np.asarray(source_x, dtype=float).copy()
        self.bonds = np.asarray(bond_pairs, dtype=int)
        self.equilibrium = np.asarray(equilibrium_A, dtype=float)
        self.protein = list(protein_indices); self.targets = target_specs
        self.geometry = geometry_module
        require(len(self.ids) == len(set(self.ids)) and self.source.shape == (len(atoms), 3)
                and np.isfinite(self.source).all(), 'Physical identity/coordinates invalid')
        require(self.bonds.ndim == 2 and self.bonds.shape[1] == 2 and
                self.equilibrium.shape == (len(self.bonds),) and np.all(self.equilibrium > 0)
                and np.isfinite(self.equilibrium).all(), 'Native bond references unavailable')
        require(np.all(self.bonds >= 0) and np.all(self.bonds < len(atoms)), 'Bond index invalid')
        inherited = np.linalg.norm(self.source[self.bonds[:, 1]]-self.source[self.bonds[:, 0]], axis=1)
        self.lower = np.minimum(SETTINGS['bond_lower_factor']*self.equilibrium, inherited)
        self.upper = np.maximum(SETTINGS['bond_upper_factor']*self.equilibrium, inherited)
        self.symbols = [a['element'] for a in atoms]
        graph = [set() for _ in atoms]
        for i, j in self.bonds:
            require(i != j and j not in graph[i], 'Duplicate/self bond')
            graph[i].add(j); graph[j].add(i)
        self.excluded = set()
        for i, neighbors in enumerate(graph):
            near = set(neighbors)
            for j in neighbors:
                near.update(graph[j])
            self.excluded.update(tuple(sorted((i, j))) for j in near if i != j)
        local = {global_index: k for k, global_index in enumerate(self.protein)}
        adapted = [dict(source_id=atoms[i]['id'], residue_id=atoms[i]['id'].rsplit('/', 1)[0],
                        name=atoms[i]['id'].rsplit('/', 1)[1], element=atoms[i]['element'],
                        resname=atoms[i]['resname']) for i in self.protein]
        pbonds = [[local[i], local[j]] for i, j in self.bonds if i in local and j in local]
        self.stereo = geometry_module.covalent_guards(adapted, pbonds, self.source[self.protein])
        self.acceptors = [i for i in self.protein if atoms[i]['element'] in ('N', 'O')]
        require(self.acceptors, 'No mapped peptide N/O atoms')
        for spec in target_specs:
            require(spec['index'] not in self.protein and spec['site_donor_indices']
                    and all(i in self.acceptors for i in spec['site_donor_indices']), 'Target/site mapping differs')

    def check(self, positions):
        x = np.asarray(positions, dtype=float)
        require(x.shape == self.source.shape and np.isfinite(x).all(), 'Invalid full physical coordinates')
        gs = self.geometry.SETTINGS
        lengths = np.linalg.norm(x[self.bonds[:, 1]]-x[self.bonds[:, 0]], axis=1)
        bad = np.flatnonzero((lengths < self.lower-1e-7) | (lengths > self.upper+1e-7))
        chiral = []; peptides = []; protein = x[self.protein]
        for row in self.stereo['chirality']:
            v = protein[row['neighbours']]; volume = float(np.linalg.det(v[:3]-v[3]))
            if volume*row['source_signed_volume_A3'] <= 0 or abs(volume) < gs['source_chirality_volume_min_A3']:
                chiral.append(self.protein[row['center']])
        for row in self.stereo['peptides']:
            try:
                angle = self.geometry.dihedral(protein[row['indices']])
            except self.geometry.InvalidArtifact:
                angle = float('nan')
            if not np.isfinite(angle) or np.cos(angle)*np.cos(row['source_omega_radian']) <= 0:
                peptides.append([self.protein[i] for i in row['indices']])
        overlaps = []
        for i, j in sorted(cKDTree(x).query_pairs(gs['new_overlap_heavy_A'])):
            if (i, j) in self.excluded:
                continue
            cutoff = gs['new_overlap_H_A'] if 'H' in (self.symbols[i], self.symbols[j]) else gs['new_overlap_heavy_A']
            current = float(np.linalg.norm(x[i]-x[j])); initial = float(np.linalg.norm(self.source[i]-self.source[j]))
            if current < min(cutoff, initial)-gs['overlap_roundoff_A']:
                overlaps.append({'indices': [i, j], 'ids': [self.ids[i], self.ids[j]],
                                 'distance_A': current, 'original_A': initial})
        states = []
        for spec in self.targets:
            i = spec['index']; limit = SETTINGS['departure_distance_A']
            nearest = float(np.linalg.norm(x[self.acceptors]-x[i], axis=1).min())
            site = float(np.linalg.norm(x[spec['site_donor_indices']]-x[i], axis=1).min())
            labels = []
            if nearest > limit:
                labels.append('departed_protein')
            elif site > limit:
                labels.append('departed_declared_site')
            water = [{'index': j, 'id': self.ids[j], 'distance_A': float(np.linalg.norm(x[j]-x[i]))}
                     for j in spec.get('coordinating_water_indices', [])]
            labels += ['lost_retained_water:' + a['id'] for a in water if a['distance_A'] > limit]
            states.append({'target': spec['label'], 'index': i, 'nearest_peptide_NO_A': nearest,
                           'nearest_declared_site_O_A': site, 'retained_waters': water, 'labels': labels})
        return {'pass': not len(bad) and not chiral and not peptides and not overlaps,
            'bond_violations': [{'indices': self.bonds[k].tolist(), 'length_A': float(lengths[k]),
                                 'lower_A': float(self.lower[k]), 'upper_A': float(self.upper[k])} for k in bad],
            'chirality_failures': chiral, 'peptide_basin_failures': peptides, 'new_severe_overlaps': overlaps,
            'state_compatible': not any(s['labels'] for s in states), 'state_labels': states,
            'maximum_source_displacement_A': float(np.linalg.norm(x-self.source, axis=1).max()),
            'source_displacement_is_not_a_gate': True, 'all_atoms_movable': len(x)}


def search(evaluator, guard, start, out, limits, lbfgs_direction, resume=None):
    """Evaluator contract: evaluate(X), difference(a,b), new_evaluations.

    Evaluations expose total_energy_kcal_mol, gradient_kcal_mol_A, positions_A,
    component_energies_kcal_mol, provenance. Checkpoint history is Cartesian.
    """
    out = Path(out); out.mkdir(exist_ok=False)
    require(set(limits) == {'accepted_steps', 'evaluations'} and
            all(type(v) is int and v > 0 for v in limits.values()), 'Invalid search limits')
    require(resume is not None, 'Only the exact step92 7CCN-La checkpoint may enter this completion attempt')
    if resume is not None:
        require(resume['settings'] == SETTINGS, 'Checkpoint optimizer settings differ')
        require(resume['limits'] == {'accepted_steps': 92, 'evaluations': 120}
                and limits == {'accepted_steps': 172, 'evaluations': 292}
                and hashlib.sha256(json.dumps(resume, sort_keys=True, separators=(',', ':'),
                    allow_nan=False).encode()).hexdigest() ==
                    '981ec084b57fa86599125619b9ecce7613c883597b23e38d3521be1eaa5c0da2',
                'Only declared172/292 completion from exact checkpoint317c2b24 is permitted')
        initial = evaluator.read_point(resume['initial_point'])
        require(np.array_equal(initial['positions_A'], start['positions_A']) and
                np.array_equal(initial['gradient_kcal_mol_A'], start['gradient_kcal_mol_A']) and
                initial['component_energies_kcal_mol'] == start['component_energies_kcal_mol'] and
                initial['total_energy_kcal_mol'] == start['total_energy_kcal_mol'],
                'Checkpoint original start/model/state differs')
        for key, bound in (('accepted_steps', 'accepted_steps'), ('distinct_evaluations', 'evaluations')):
            require(type(resume[key]) is int and 0 <= resume[key] <= resume['limits'][bound]
                    and resume[key] <= limits[bound], 'Checkpoint counters exceed declared budget')
        require(np.isfinite(resume['trust_A']) and SETTINGS['minimum_trust_A'] <= resume['trust_A']
                <= SETTINGS['maximum_trust_A'], 'Invalid checkpoint trust radius')
        require(len(resume['history']) <= SETTINGS['history'], 'Checkpoint history too long')
        for h in resume['history']:
            require(np.asarray(h['s']).shape == np.asarray(h['y']).shape == (np.asarray(start['positions_A']).size,)
                    and np.isfinite(h['s']).all() and np.isfinite(h['y']).all()
                    and np.isfinite(h['rho']) and h['rho'] > 0, 'Invalid checkpoint L-BFGS history')
    point = start if resume is None else evaluator.read_point(resume['point'])
    x = np.asarray(point['positions_A'], dtype=float)
    history = [] if resume is None else [(np.asarray(h['s']), np.asarray(h['y']), h['rho']) for h in resume['history']]
    trust = SETTINGS['initial_trust_A'] if resume is None else resume['trust_A']
    accepted = 0 if resume is None else resume['accepted_steps']
    prior_evaluations = 0 if resume is None else resume['distinct_evaluations']
    call_baseline = evaluator.new_evaluations
    trace = []; evaluated = []; termination = 'not_started'; failed = None
    prior_points = [] if resume is None else resume.get('evaluated_points', [])
    prior_requests = [] if resume is None else resume.get('evaluation_requests', [])
    def count():
        return prior_evaluations + evaluator.new_evaluations - call_baseline
    def checkpoint(final=False):
        value = {'point': point['provenance'], 'accepted_steps': accepted, 'distinct_evaluations': count(),
                 'trust_A': trust, 'history': [{'s': s.tolist(), 'y': y.tolist(), 'rho': rho} for s, y, rho in history],
                 'initial_point': start['provenance'], 'limits': limits, 'settings': SETTINGS,
                 'evaluated_points': prior_points+evaluated,
                 'evaluation_requests': prior_requests+list(getattr(evaluator, 'evaluation_requests', [])),
                 'trial_trace': trace, 'failure': failed}
        path = out / ('CHECKPOINT_FINAL.json' if final else 'CHECKPOINT_%03d.json' % accepted)
        with path.open('x') as handle:
            json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')
        return evaluator.pin(path)
    check = guard.check(x)
    require(check['pass'] and check['state_compatible'], 'Search start physically unavailable')
    checkpoint_pin = checkpoint()
    try:
        while True:
            g = np.asarray(point['gradient_kcal_mol_A'], dtype=float)
            metrics = force_metrics(g); pg, rigid = remove_rigid(x, g)
            entry = {'accepted_steps': accepted, 'point': point['provenance'], 'distinct_evaluations': count(),
                     'raw_full_force': metrics, 'rigid_removed_force': force_metrics(pg), 'rigid_projection': rigid,
                     'geometry': check, 'trust_A': trust, 'trials': []}
            trace.append(entry)
            if not check['state_compatible']:
                termination = 'state_change'; break
            if metrics['converged']:
                termination = 'full_force_converged'; break
            if accepted >= limits['accepted_steps']:
                termination = 'accepted_step_limit'; break
            if count() >= limits['evaluations']:
                termination = 'evaluation_limit'; break
            d = lbfgs_direction(pg.ravel(), history).reshape(x.shape)
            d, _ = remove_rigid(x, d)
            if float(np.sum(pg*d)) >= -1e-14:
                history = []; d = -pg
            maximum = float(np.linalg.norm(d, axis=1).max())
            if maximum < 1e-14:
                termination = 'zero_rigid_free_direction'; break
            d *= min(1., trust/maximum)
            slope = float(np.sum(g*d)); require(slope < 0, 'Non-descent Cartesian direction')
            next_state = None
            for trial_index in range(SETTINGS['backtracks']):
                if count() >= limits['evaluations']:
                    break
                multiplier = .5**trial_index; trial_x = x+multiplier*d
                record = {'multiplier': multiplier, 'evaluation_requested': False, 'accepted': False}
                entry['trials'].append(record)
                try:
                    trial_guard = guard.check(trial_x); record['geometry'] = trial_guard
                except ValueError as exc:
                    record['reason'] = 'geometry:' + str(exc); continue
                if not trial_guard['pass']:
                    record['reason'] = 'physical_guard'; continue
                record['evaluation_requested'] = True
                trial = evaluator.evaluate(trial_x)
                evaluated.append(trial['provenance']); record['point'] = trial['provenance']
                work = evaluator.difference(trial, point)
                record['work'] = work; record['Armijo_limit_kcal_mol'] = SETTINGS['armijo']*multiplier*slope
                if work['composite_kcal_mol'] <= record['Armijo_limit_kcal_mol']:
                    record['accepted'] = True
                    next_state = trial_x, trial, trial_guard, multiplier
                    break
                record['reason'] = 'Armijo_energy'
            with (out/'TRACE.jsonl').open('a') as handle:
                handle.write(json.dumps(entry, allow_nan=False)+'\n')
            if next_state is None:
                termination = 'evaluation_limit' if count() >= limits['evaluations'] else 'no_admissible_decreasing_trial'
                break
            nx, new, new_guard, multiplier = next_state
            ng, _ = remove_rigid(nx, np.asarray(new['gradient_kcal_mol_A']))
            s = (nx-x).ravel(); y = (ng-pg).ravel(); sy = float(s@y)
            if sy > 1e-10*np.linalg.norm(s)*np.linalg.norm(y) and sy > 1e-16:
                history.append((s, y, 1./sy)); history = history[-SETTINGS['history']:]
            ratio = evaluator.difference(new, point)['composite_kcal_mol']/(multiplier*slope)
            if multiplier == 1. and ratio > .75:
                trust = min(SETTINGS['maximum_trust_A'], 2*trust)
            elif ratio < .25:
                trust = max(SETTINGS['minimum_trust_A'], .5*trust)
            x, point, check = nx, new, new_guard; accepted += 1
            checkpoint_pin = checkpoint()
    except Exception as exc:
        termination = 'execution_failure'; failed = repr(exc)
    checkpoint_pin = checkpoint(final=True)
    final = force_metrics(np.asarray(point['gradient_kcal_mol_A']))
    physical = check['pass'] and check['state_compatible'] and failed is None
    return {'status': 'complete' if failed is None else 'unavailable', 'failure': failed,
            'termination': termination, 'converged': final['converged'] and physical,
            'physical_state_admitted': physical, 'terminal_point': point['provenance'],
            'start_point': start['provenance'], 'checkpoint': checkpoint_pin,
            'accepted_steps': accepted, 'distinct_evaluations': count(), 'final_geometry': check,
            'final_full_force': final, 'work_from_start': evaluator.difference(point, start),
            'trace': trace, 'evaluated_points': evaluated, 'settings': SETTINGS, 'limits': limits}
