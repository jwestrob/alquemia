"""Contained whole-chain La/Dy accommodation, using existing warm/native runners.

No absolute ion affinity or equilibrium occupancy is inferred. All energies are
conditional on the source composition; different sources are never geometry-pooled.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import shutil
import time

import numpy as np
from affordable_common import (InvalidArtifact, HA_TO_KCAL, cache_key, read_json,
                               record, verify, write_new, xyz)

PROTOCOL = 'nikasha_LanM_whole_chain_LaDy_occupancy_pool_v1'
EV_TO_KCAL = 23.06054783061903
SETTINGS = {
    'coordinate_box_A': .45, 'maxiter': 60, 'maxfun': 180,
    'ftol': 1e-10, 'gtol_eV_per_A': 1e-3,
    'selection_tolerance_kcal_mol': .1,
    'bond_ratio_limits': [.8, 1.2],
    'minimum_CA_oriented_volume_fraction': .2,
    'new_nonbonded_heavy_clash_A': 1.2,
    'adapter_energy_tolerance_eV': 1e-6,
    'adapter_force_tolerance_eV_per_A': 1e-6,
}
ORCA = '/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca'
METALS = ('La', 'Dy')


def states(manifest, selected=None):
    m = read_json(manifest)
    entries = m['states']
    out = []
    for r in entries:
        # Preparation records are pinned, not inferred from directory names.
        p = verify(r['state']) if 'state' in r else (
            verify(r) if 'path' in r else Path(verify(r['mapping'])).parent/'state.json')
        s = read_json(p)
        if 'state_id' in r and s != r:
            raise InvalidArtifact('inline and physical state records differ')
        if selected and s['state_id'] not in selected:
            continue
        if s['status'] != 'prepared' or s['n'] not in (2, 4):
            raise InvalidArtifact('unsupported prepared occupancy')
        mapping = read_json(verify(s['mapping']))
        if mapping['metal_indices'] != s['metal_indices']:
            raise InvalidArtifact('metal mapping differs')
        for key in ('source_preparation', 'source_water_inventory'):
            verify(s[key])
        a, b = (xyz(verify(s['endpoints'][z]['xyz'])) for z in METALS)
        if len(a) != s['n_atoms'] or len(b) != len(a) or len(mapping['atoms']) != len(a):
            raise InvalidArtifact('atom inventory differs')
        if not np.array_equal(np.array([x[1:] for x in a]), np.array([x[1:] for x in b])):
            raise InvalidArtifact('paired origins differ')
        for i, (aa, bb) in enumerate(zip(a, b)):
            if i in s['metal_indices']:
                if (aa[0], bb[0]) != METALS:
                    raise InvalidArtifact('metal identity differs')
            elif aa[0] != bb[0]:
                raise InvalidArtifact('nonmetal composition differs')
        for z in METALS:
            e = s['endpoints'][z]
            if e['charge'] != s['charge'] or e['native_effective_multiplicity'] != 1:
                raise InvalidArtifact('charge/native state differs')
            expected = 1 if z == 'La' else 1 + 5*s['n']
            if e['physical_multiplicity'] != expected or (e['all_electron_count']-expected+1) % 2:
                raise InvalidArtifact('physical multiplicity/parity differs')
        out.append((p, s, mapping))
    if not out or selected and {s['state_id'] for _, s, _ in out} != set(selected):
        raise InvalidArtifact('selected state absent')
    return out


def geometry_check(q0, q, mapping):
    """Fixed chemistry gates; no energy/label-dependent admission rules."""
    from scipy.spatial import cKDTree
    q0, q = np.asarray(q0), np.asarray(q)
    if q.shape != q0.shape or not np.isfinite(q).all():
        return {'admitted': False, 'reason': 'nonfinite_or_mismatched_coordinates'}
    delta = q-q0
    bad_bonds = []
    bonded = set()
    for b in mapping['bonds']:
        i, j = b['indices']; bonded.add(tuple(sorted((i, j))))
        d0 = np.linalg.norm(q0[i]-q0[j]); d = np.linalg.norm(q[i]-q[j])
        ratio = float(d/d0)
        if not SETTINGS['bond_ratio_limits'][0] <= ratio <= SETTINGS['bond_ratio_limits'][1]:
            bad_bonds.append({'indices': [i, j], 'origin_A': float(d0), 'candidate_A': float(d), 'ratio': ratio})
    residues = {}
    for a in mapping['atoms']:
        if a['kind'] != 'protein_source':
            continue
        src = a['source']; k = (src['chain'], src['resid'], src['insertion_code'])
        residues.setdefault(k, {})[src['name']] = a['index']
    bad_chiral = []
    for key, indices in residues.items():
        if not {'CA', 'N', 'C', 'CB'} <= indices.keys():
            continue  # Gly has no C-alpha stereocentre.
        ca, n, c, cb = [indices[k] for k in ('CA', 'N', 'C', 'CB')]
        def vol(p):
            return float(np.dot(p[n]-p[ca], np.cross(p[c]-p[ca], p[cb]-p[ca])))
        v0, v = vol(q0), vol(q)
        if abs(v0) < 1e-8 or v/v0 < SETTINGS['minimum_CA_oriented_volume_fraction']:
            bad_chiral.append({'residue': list(key), 'origin_volume_A3': v0, 'candidate_volume_A3': v})
    heavy = [a['index'] for a in mapping['atoms'] if a['element'] != 'H']
    def clashes(p):
        return {tuple(sorted((heavy[i], heavy[j]))) for i, j in
                cKDTree(p[heavy]).query_pairs(SETTINGS['new_nonbonded_heavy_clash_A'])
                if tuple(sorted((heavy[i], heavy[j]))) not in bonded}
    old_clash, new_clash = clashes(q0), clashes(q)
    excessive = float(np.abs(delta).max()) > SETTINGS['coordinate_box_A']+1e-9
    return {'admitted': not (bad_bonds or bad_chiral or new_clash-old_clash or excessive),
            'bad_bonds': bad_bonds, 'bad_CA_stereocentres': bad_chiral,
            'new_severe_heavy_clashes': [list(p) for p in sorted(new_clash-old_clash)],
            'origin_severe_heavy_clashes': [list(p) for p in sorted(old_clash)],
            'exceeds_box': excessive, 'max_atom_displacement_A': float(np.linalg.norm(delta, axis=1).max()),
            'RMS_atom_displacement_A': float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            'boundary_components': int(np.sum(np.abs(delta) >= SETTINGS['coordinate_box_A']-1e-6))}


def write_xyz(path, symbols, positions):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        f.write(f'{len(symbols)}\n{PROTOCOL}\n')
        for symbol, row in zip(symbols, positions):
            f.write(symbol+' '+ ' '.join(format(float(v), '.17g') for v in row)+'\n')
    return record(p)


def prepare(manifest, software, agreement, output, selected=None):
    source = states(manifest, selected)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir()
    pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dst = impl/p.name; shutil.copyfile(p, dst); pins[p.name] = record(dst)
    sw = read_json(software); verify(sw['checkpoint'])
    m = {'protocol_id': PROTOCOL, 'settings': SETTINGS, 'preparation': record(manifest),
         'states': [record(p) for p, _, _ in source], 'software': record(software),
         'model': sw['checkpoint'], 'agreement': record(agreement), 'implementation': pins,
         'production_changed': False, 'new_DFT': 0}
    write_new(out/'manifest.json', m)
    return {'manifest': record(out/'manifest.json'), 'state_count': len(source),
            'candidate_endpoint_ceiling': 6*len(source), 'native_scalar_ceiling': 12*len(source)}


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS:
        raise InvalidArtifact('protocol/settings mismatch')
    for k in ('preparation', 'software', 'model', 'agreement'):
        verify(m[k])
    for pin in m['implementation'].values():
        verify(pin)
    if record(__file__)['sha256'] != m['implementation'][Path(__file__).name]['sha256']:
        raise InvalidArtifact('use the pinned implementation')
    return m


class WarmMACE:
    def __init__(self, model_pin):
        import torch
        from mace.calculators import mace_omol
        if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():
            raise InvalidArtifact('GPU allocation required')
        torch.set_num_threads(min(8, int(os.environ['SLURM_CPUS_PER_TASK'])))
        torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats()
        start = time.monotonic()
        self.calc = mace_omol(model=str(verify(model_pin)), device='cuda', default_dtype='float64')
        self.model = self.calc.models[0]
        if self.calc.head != 'omol' or self.calc.energy_units_to_eV != 1.:
            raise InvalidArtifact('native OMOL head or units differ')
        for p in self.model.parameters():
            p.requires_grad_(False)
        self.adapted = False
        self.calls = 0
        self.load_seconds = time.monotonic()-start

    def adapt(self):
        from mace_omol_gradients import install
        install(self.model); self.adapted = True

    def evaluate(self, symbols, q, charge, mult, metal_indices, forces=True):
        import torch
        from ase import Atoms
        from mace_omol import input_batch
        from torch.utils.checkpoint import set_checkpoint_early_stop
        if self.adapted:
            for b in [*self.model.interactions, *self.model.products]:
                b._alquemia_forward_calls = []; b._alquemia_kernel_calls = 0
        a = Atoms(symbols, positions=q, pbc=False); a.info.update(charge=charge, spin=mult)
        batch = self.calc._atoms_to_batch(a)
        state = input_batch(self.calc, a, charge, mult, batch, metal_indices[0])
        state['metal_indices'] = metal_indices
        tick = time.monotonic(); self.calls += 1
        with set_checkpoint_early_stop(False), (contextlib.nullcontext() if forces else torch.no_grad()):
            v = self.model(batch.to_dict(), training=False, compute_force=forces,
                compute_virials=False, compute_stress=False, compute_displacement=False,
                compute_hessian=False, compute_edge_forces=False, compute_atomic_stresses=False)
        energy = float(v['energy'].detach().cpu().item())
        f = v['forces'].detach().cpu().numpy() if forces else None
        if not math.isfinite(energy) or forces and (f.shape != np.shape(q) or not np.isfinite(f).all()):
            raise InvalidArtifact('nonfinite MACE result')
        r = {'energy_eV': energy, 'wall_seconds': time.monotonic()-tick, 'input_state': state,
             'adapter': None, 'forces_eV_per_A': None if f is None else f.tolist()}
        if self.adapted:
            from mace_omol_gradients import receipt
            r['adapter'] = receipt(self.model, forces)
        return r


def run_mace(manifest, output, selected=None, capability_only=False):
    import torch
    from scipy.optimize import minimize
    m = validate(manifest); source = states(verify(m['preparation']), selected or
        [read_json(verify(p))['state_id'] for p in m['states']])
    if not {record(p)['sha256'] for p, _, _ in source} <= {p['sha256'] for p in m['states']}:
        raise InvalidArtifact('state outside finite manifest')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic(); calc = WarmMACE(m['model']); records = []
    for p, s, mapping in source:
        d = out/s['state_id']; d.mkdir()
        rows = xyz(verify(s['endpoints']['La']['xyz']))
        q0 = np.array([a[1:] for a in rows]); candidates = {'origin': q0}
        case = {'state': record(p), 'state_id': s['state_id'], 'proposals': {}, 'cells': [], 'status': 'incomplete'}
        # One real full-chain native/adapter force equivalence test before search.
        if not calc.adapted:
            try:
                ref = calc.evaluate([a[0] for a in rows], q0, s['charge'], 1, s['metal_indices'])
                write_new(d/'native_reference.json', ref)
            except Exception as exc:
                write_new(d/'native_reference_failure.json', {'reason': f'{type(exc).__name__}: {exc}',
                    'calls_attempted': calc.calls})
                # Preserve the failed native call; qualified archived exact adapter
                # still requires an actual local energy/force equivalence reference.
                raise
            calc.adapt()
            adapted = calc.evaluate([a[0] for a in rows], q0, s['charge'], 1, s['metal_indices'])
            write_new(d/'adapter_reference.json', adapted)
            de = abs(ref['energy_eV']-adapted['energy_eV'])
            df = float(np.max(np.abs(np.array(ref['forces_eV_per_A'])-adapted['forces_eV_per_A'])))
            qualification = {'energy_difference_eV': de, 'max_force_difference_eV_per_A': df,
                'passed': de <= SETTINGS['adapter_energy_tolerance_eV'] and df <= SETTINGS['adapter_force_tolerance_eV_per_A']}
            write_new(d/'adapter_qualification.json', qualification)
            if not qualification['passed']:
                raise InvalidArtifact('adapter qualification failed')
        if not capability_only:
            for metal in METALS:
                e = s['endpoints'][metal]; sym = [r[0] for r in xyz(verify(e['xyz']))]
                evaluations = []; store = {}
                def fun(x):
                    q = x.reshape(q0.shape)
                    key = q.tobytes()
                    if key not in store:
                        result = calc.evaluate(sym, q, s['charge'], e['physical_multiplicity'], s['metal_indices'])
                        store[key] = result
                        n = len(evaluations)
                        pin = write_xyz(d/metal/f'eval_{n:03d}.xyz', sym, q)
                        result.update(xyz=pin, evaluation=n)
                        write_new(d/metal/f'eval_{n:03d}.json', result)
                        evaluations.append({'energy_eV': result['energy_eV'], 'xyz': pin,
                                            'wall_seconds': result['wall_seconds']})
                    result = store[key]
                    return result['energy_eV']-evaluations[0]['energy_eV'], -np.array(result['forces_eV_per_A']).ravel()
                try:
                    qflat = q0.ravel()
                    answer = minimize(fun, qflat, method='L-BFGS-B', jac=True,
                        bounds=list(zip(qflat-SETTINGS['coordinate_box_A'], qflat+SETTINGS['coordinate_box_A'])),
                        options={k: SETTINGS[k] for k in ('maxiter', 'maxfun', 'ftol')} |
                                {'gtol': SETTINGS['gtol_eV_per_A']})
                    q = answer.x.reshape(q0.shape)
                    guard = geometry_check(q0, q, mapping)
                    cand = 'proposal_'+metal
                    case['proposals'][metal] = {'status': 'admitted' if guard['admitted'] else 'rejected_geometry',
                        'geometry': guard, 'optimizer_success': bool(answer.success), 'message': str(answer.message),
                        'iterations': int(answer.nit), 'function_evaluations': int(answer.nfev),
                        'vacuum_work_eV': float(answer.fun), 'actual_evaluations': evaluations,
                        'xyz': write_xyz(d/(cand+'.xyz'), sym, q)}
                    if guard['admitted']:
                        candidates[cand] = q
                except Exception as exc:
                    case['proposals'][metal] = {'status': 'failed', 'reason': f'{type(exc).__name__}: {exc}',
                                                'actual_evaluations': evaluations}
                write_new(d/(metal+'_proposal.json'), case['proposals'][metal])
                store.clear()
        # Exact-coordinate deduplication only. Candidate count is not a population.
        seen = {}; case['candidate_aliases'] = {}
        for name, q in candidates.items():
            key = q.tobytes()
            if key in seen:
                case['candidate_aliases'][name] = seen[key]; continue
            seen[key] = name
            for metal in METALS:
                e = s['endpoints'][metal]; sym = [r[0] for r in xyz(verify(e['xyz']))]
                pin = write_xyz(d/'candidates'/name/(metal+'.xyz'), sym, q)
                cell = {'candidate': name, 'metal': metal, 'xyz': pin, 'charge': s['charge'],
                        'physical_multiplicity': e['physical_multiplicity'], 'status': 'failed', 'energy_eV': None}
                try:
                    # Evaluate serialized coordinates so MACE and GFN2 see exactly the same geometry.
                    qq = np.array([a[1:] for a in xyz(verify(pin))])
                    val = calc.evaluate(sym, qq, s['charge'], e['physical_multiplicity'], s['metal_indices'], False)
                    cell.update(status='complete', **val)
                except Exception as exc:
                    cell['reason'] = f'{type(exc).__name__}: {exc}'
                case['cells'].append(cell)
                write_new(d/'candidates'/name/(metal+'.json'), cell)
        failed_search = any(v['status']=='failed' for v in case['proposals'].values())
        case['status'] = ('failed_search' if failed_search else
            ('complete' if all(c['status']=='complete' for c in case['cells']) else 'failed_cell'))
        case['search_status'] = ('not_requested' if capability_only else
            'unavailable' if failed_search else 'complete_finite_search')
        write_new(d/'result.json', case); records.append(case)
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'rows': records,
        'capability_only': capability_only, 'selected_states': [s['state_id'] for _, s, _ in source],
        'molecular_calls_attempted': calc.calls,
        'wall_seconds': time.monotonic()-start, 'model_load_seconds': calc.load_seconds,
        'peak_host_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'peak_cuda_allocated_bytes': torch.cuda.max_memory_allocated(),
        'device': torch.cuda.get_device_name(), 'job_id': os.environ['SLURM_JOB_ID'],
        'allocated_CPUs': int(os.environ['SLURM_CPUS_ON_NODE']),
        'versions': {k: importlib.metadata.version(k) for k in ('torch', 'mace-torch', 'ase', 'scipy')}}
    write_new(out/'result.json', result)
    return {k: v for k, v in result.items() if k != 'rows'}


def native_prepare(manifest, output, mace_result=None, selected=None, ranks=8, workers=1, reuse=None):
    from strict_native_pool import recipe
    from compact_solvation import completed
    m = validate(manifest)
    src = states(verify(m['preparation']), selected or [read_json(verify(p))['state_id'] for p in m['states']])
    if not {record(p)['sha256'] for p, _, _ in src} <= {p['sha256'] for p in m['states']}:
        raise InvalidArtifact('state outside finite manifest')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    mr = read_json(mace_result) if mace_result else None
    if mr and mr['manifest'] != record(manifest):
        raise InvalidArtifact('proposal manifest differs')
    if mr:
        validate_mace_rows(m, mr)
    prior = {}
    for rp in reuse or []:
        rm = read_json(rp)
        if rm['parent_manifest'] != record(manifest):
            raise InvalidArtifact('reuse protocol/inputs differ')
        for t in rm['tasks']:
            pin = completed(rp, t['task_id'])
            if pin:
                prior[t['scientific_key']] = pin
    tasks = []; reused = []
    for p, s, mapping in src:
        if mr:
            case = next(r for r in mr['rows'] if r['state_id']==s['state_id'])
            cells = case['cells']
        else:
            cells = [dict(candidate='origin', metal=z, xyz=s['endpoints'][z]['xyz']) for z in METALS]
        for cell in cells:
            for medium in ('vacuum', 'alpb'):
                text = recipe(s['charge'], medium, 'fresh')
                scientific_key = cache_key({'protocol': PROTOCOL, 'state': record(p),
                    'atoms': xyz(verify(cell['xyz'])), 'input_text': text, 'ranks': ranks})
                tid = '__'.join((s['state_id'], cell['candidate'], cell['metal'], medium))
                t = {'task_id': tid, 'state_id': s['state_id'], 'candidate': cell['candidate'],
                    'metal': cell['metal'], 'medium': medium, 'charge': s['charge'], 'multiplicity': 1,
                    'physical_multiplicity': s['endpoints'][cell['metal']]['physical_multiplicity'],
                    'scientific_key': scientific_key, 'gradient_requested': False}
                if scientific_key in prior:
                    reused.append(dict(t, actual=prior[scientific_key])); continue
                d = out/'tasks'/tid; d.mkdir(parents=True)
                shutil.copyfile(verify(cell['xyz']), d/'core.xyz')
                (d/'endpoint.inp').write_text(text)
                t.update(input=record(d/'endpoint.inp'), xyz=record(d/'core.xyz'), output_path=str(d/'endpoint.out'))
                tasks.append(t)
    nm = {'protocol_id': PROTOCOL, 'parent_manifest': record(manifest), 'agreement': m['agreement'],
        'tasks': tasks, 'reused': reused, 'orca': record(ORCA), 'mace_result': record(mace_result) if mace_result else None,
        'execution_resources': {'mpi_ranks': ranks, 'concurrent_tasks': workers},
        'execution_policy': {k: m['implementation'][n] for k, n in
            (('task_runner', 'run_orca_task_manifest.py'), ('runtime_renderer', 'render_orca_runtime_input.py'))}}
    write_new(out/'manifest.json', nm)
    from affordable_workflow import dry_run
    result = dry_run(out/'manifest.json'); write_new(out/'PREFLIGHT.json', result)
    return result


def choose(energies):
    best = min(energies, key=lambda k: (energies[k], k))
    selected = best if energies['origin']-energies[best] >= SETTINGS['selection_tolerance_kcal_mol'] else 'origin'
    return {'mathematical_minimum': best, 'selected': selected,
            'minimum_kcal_mol': energies[best], 'selected_kcal_mol': energies[selected],
            'accommodation_work_kcal_mol': energies[selected]-energies['origin']}


def validate_mace_rows(m, result):
    declared = {read_json(verify(p))['state_id']: p for p in m['states']}
    selected = result['selected_states']
    if not selected or len(set(selected)) != len(selected) or not set(selected) <= declared.keys():
        raise InvalidArtifact('MACE state selection differs')
    rows = result['rows']
    if len(rows) != len(selected) or {r['state_id'] for r in rows} != set(selected):
        raise InvalidArtifact('MACE selected state denominator differs')
    for row in rows:
        if row['state'] != declared[row['state_id']]:
            raise InvalidArtifact('MACE source state pin differs')
        s = read_json(verify(row['state']))
        required = {'origin'}
        if not result['capability_only']:
            if set(row['proposals']) != set(METALS):
                raise InvalidArtifact('missing search result')
            required |= {'proposal_'+z for z, v in row['proposals'].items() if v['status']=='admitted'}
        required -= row['candidate_aliases'].keys()
        expected = {(q, z) for q in required for z in METALS}
        cells = row['cells']
        if len(cells) != len(expected) or {(c['candidate'], c['metal']) for c in cells} != expected:
            raise InvalidArtifact('incomplete required geometry pool')
        for q in required:
            pair = [next(c for c in cells if c['candidate']==q and c['metal']==z) for z in METALS]
            a, b = [xyz(verify(c['xyz'])) for c in pair]
            if len(a) != s['n_atoms'] or len(b) != len(a) or any(x[1:] != y[1:] for x, y in zip(a, b)):
                raise InvalidArtifact('cross-metal candidate geometry differs')
            for i, (x, y) in enumerate(zip(a, b)):
                expected_symbols = METALS if i in s['metal_indices'] else (x[0], x[0])
                if (x[0], y[0]) != expected_symbols:
                    raise InvalidArtifact('cross-metal composition differs')
            for c in pair:
                if c['charge'] != s['charge'] or c['physical_multiplicity'] != s['endpoints'][c['metal']]['physical_multiplicity']:
                    raise InvalidArtifact('cross-score electronic state differs')


def collect(manifest, mace_result, native_manifests, output):
    from compact_solvation import completed, diagnostics
    from structure_informed_starts import scf_details
    m = validate(manifest); mr = read_json(mace_result)
    if mr['manifest'] != record(manifest):
        raise InvalidArtifact('MACE origin differs')
    validate_mace_rows(m, mr)
    scalars = {}; failures = []
    for npth in native_manifests:
        nm = read_json(npth)
        if nm['parent_manifest'] != record(manifest):
            raise InvalidArtifact('native parent differs')
        combined = [(t, npth) for t in nm['tasks']]
        for reuse in nm.get('reused', []):
            oldpath = verify(reuse['actual']['manifest']); old = read_json(oldpath)
            t = next(t for t in old['tasks'] if t['task_id']==reuse['actual']['task_id'])
            if (old['parent_manifest'] != record(manifest) or t['scientific_key'] != reuse['scientific_key']
                    or any(t[k] != reuse[k] for k in ('state_id', 'candidate', 'metal', 'medium'))):
                raise InvalidArtifact('native reused scientific identity differs')
            combined.append((t, oldpath))
        for t, actual_manifest in combined:
            key = (t['state_id'], t['candidate'], t['metal'], t['medium'])
            try:
                pin = completed(actual_manifest, t['task_id'])
                if not pin:
                    raise InvalidArtifact('missing/failed native receipt')
                audit = diagnostics(pin, t)
                scf = scf_details(verify(pin['output']).read_text())
                if scf['energy']['tolerance'] != 1e-10 or not scf['native_mixer_observed']:
                    raise InvalidArtifact('native tolerance/mixer differs')
                if audit['charge_sanity_status'] != 'pass':
                    raise InvalidArtifact('native charge sanity failed')
                scalars[key] = dict(pin, audit=audit, scf=scf)
            except Exception as exc:
                failures.append({'key': list(key), 'reason': str(exc)})
    rows = []
    for state_pin in m['states']:
        s = read_json(verify(state_pin))
        case = next((r for r in mr['rows'] if r['state_id']==s['state_id']), None)
        if case is None:
            rows.append({'state_id': s['state_id'], 'protein': s['protein'], 'source_id': s['source_id'],
                'occupied_sites': s['occupied_sites'], 'n': s['n'], 'status': 'unavailable',
                'reason': 'MACE_state_not_executed', 'selected_Dy_minus_La_kcal_mol': None,
                'static_Dy_minus_La_kcal_mol': None, 'selections': None, 'cells': []})
            continue
        s = read_json(verify(case['state'])); cells = []; energy_rows = {z: {} for z in METALS}
        for c in case['cells']:
            k = (s['state_id'], c['candidate'], c['metal'])
            row = {**c, 'composite_kcal_mol': None, 'solvent_transfer_kcal_mol': None}
            if c['status']=='complete' and (*k, 'alpb') in scalars and (*k, 'vacuum') in scalars:
                v, a = scalars[*k, 'vacuum'], scalars[*k, 'alpb']
                sol = (a['energy_hartree']-v['energy_hartree'])*HA_TO_KCAL
                total = c['energy_eV']*EV_TO_KCAL+sol
                row.update(composite_kcal_mol=total, solvent_transfer_kcal_mol=sol,
                           native={'alpb': a, 'vacuum': v})
                energy_rows[c['metal']][c['candidate']] = total
            else:
                row['status'] = 'unavailable'
            cells.append(row)
        required = {c['candidate'] for c in case['cells']}
        ok = case['status']=='complete' and all(set(energy_rows[z])==required for z in METALS)
        row = {'state_id': s['state_id'], 'protein': s['protein'], 'source_id': s['source_id'],
            'occupied_sites': s['occupied_sites'], 'n': s['n'], 'cells': cells, 'proposals': case['proposals'],
            'status': 'complete' if ok else 'unavailable', 'selected_Dy_minus_La_kcal_mol': None,
            'static_Dy_minus_La_kcal_mol': None, 'selections': None}
        if all('origin' in energy_rows[z] for z in METALS):
            row['static_Dy_minus_La_kcal_mol'] = energy_rows['Dy']['origin']-energy_rows['La']['origin']
        if ok:
            selections = {z: choose(energy_rows[z]) for z in METALS}
            row.update(selections=selections,
                selected_Dy_minus_La_kcal_mol=selections['Dy']['selected_kcal_mol']-selections['La']['selected_kcal_mol'],
                static_Dy_minus_La_kcal_mol=energy_rows['Dy']['origin']-energy_rows['La']['origin'])
        rows.append(row)
    comparisons = []
    for h in (r for r in rows if r['protein']=='Hans'):
        mex = next((r for r in rows if r['protein']=='Mex' and r['occupied_sites']==h['occupied_sites']), None)
        available = h['status']=='complete' and mex is not None and mex['status']=='complete'
        comparisons.append({'Hans_state': h['state_id'], 'Mex_state': mex['state_id'] if mex else None,
            'occupied_sites': h['occupied_sites'], 'status': 'complete' if available else 'unavailable',
            'static_relative_La_preference_kcal_mol': h['static_Dy_minus_La_kcal_mol']-mex['static_Dy_minus_La_kcal_mol'] if available else None,
            'accommodated_relative_La_preference_kcal_mol': h['selected_Dy_minus_La_kcal_mol']-mex['selected_Dy_minus_La_kcal_mol'] if available else None})
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'mace_result': record(mace_result),
        'native_manifests': [record(p) for p in native_manifests], 'rows': rows,
        'comparisons': comparisons, 'native_failures': failures,
        'absolute_affinity_available': False, 'occupancy_probabilities_available': False}
    write_new(output, result)
    return {'states': len(rows), 'available': sum(r['status']=='complete' for r in rows), 'comparisons': comparisons,
            'native_failure_count': len(failures)}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    pp = sub.add_parser('prepare'); pp.add_argument('--preparation', required=True); pp.add_argument('--software', required=True)
    pp.add_argument('--agreement', required=True); pp.add_argument('--output', required=True); pp.add_argument('--state', action='append')
    for name in ('mace', 'native-prepare', 'collect'):
        sp = sub.add_parser(name); sp.add_argument('--manifest', required=True); sp.add_argument('--output', required=True)
        if name != 'collect': sp.add_argument('--state', action='append')
        if name == 'mace': sp.add_argument('--capability-only', action='store_true')
        elif name == 'native-prepare':
            sp.add_argument('--mace-result'); sp.add_argument('--ranks', type=int, default=8)
            sp.add_argument('--workers', type=int, default=1); sp.add_argument('--reuse', action='append')
        else:
            sp.add_argument('--mace-result', required=True); sp.add_argument('--native-manifest', action='append', required=True)
    a = p.parse_args()
    if a.op == 'prepare': r = prepare(a.preparation, a.software, a.agreement, a.output, a.state)
    elif a.op == 'mace': r = run_mace(a.manifest, a.output, a.state, a.capability_only)
    elif a.op == 'native-prepare': r = native_prepare(a.manifest, a.output, a.mace_result, a.state, a.ranks, a.workers, a.reuse)
    else: r = collect(a.manifest, a.mace_result, a.native_manifest, a.output)
    print(json.dumps(r, indent=2))


if __name__ == '__main__':
    main()
