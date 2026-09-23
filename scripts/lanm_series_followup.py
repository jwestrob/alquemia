"""Bounded Hans/Mex La/Dy preparation; separate physical and f-in-core states."""
from __future__ import annotations
import argparse
import importlib.util
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz

PROTOCOL = 'LanM_Hans_Mex_LaDy_native_OMOL_fcore_GFN2_selfcontinuation_v1'
METALS = ('La', 'Dy')
SITES = ('EF1', 'EF2', 'EF3')
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'La': 57, 'Dy': 66}
PHYSICAL_MULTIPLICITY = {'La': 1, 'Dy': 6}
EFFECTIVE_MULTIPLICITY = {'La': 1, 'Dy': 1}


def recipe(charge, medium, continuation=False):
    if medium not in ('vacuum', 'alpb'):
        raise InvalidArtifact('unsupported medium')
    return ('! Native-GFN2-xTB' + ('' if continuation else ' NoAutostart')
            + (' ALPB(Water)' if medium == 'alpb' else '')
            + '\n%maxcore 2000\n%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            + '%scf\n MaxIter 500\n SmearTemp 300\n UseXTBMixer true\nend\n'
            + f'* xyzfile {charge} 1 core.xyz\n')


def state(rows, charge, metal, parameters):
    if rows[0][0] != metal or metal not in METALS or any(a[0] not in Z for a in rows):
        raise InvalidArtifact('unsupported atom identity/order')
    if sum(a[0] in METALS for a in rows) != 1:
        raise InvalidArtifact('exactly one selected lanthanide required')
    ne = sum(Z[a[0]] for a in rows) - charge
    mult = PHYSICAL_MULTIPLICITY[metal]
    valence = sum(sum(parameters['element'][a[0]]['refocc']) for a in rows) - charge
    if ne % 2 != (mult - 1) % 2 or valence % 2 != 0:
        raise InvalidArtifact('physical/effective electron parity inconsistent')
    p = parameters['element'][metal]
    if p['shells'] != ['5d', '6s', '6p'] or p['refocc'] != [1., 1., 1.]:
        raise InvalidArtifact('declared f-in-core valence representation differs')
    return {'physical_multiplicity': mult, 'effective_GFN_multiplicity': 1,
            'all_electron_count': ne, 'GFN_valence_electron_count': int(valence),
            'metal_reference_occupation': p['refocc'], 'metal_effective_shells': p['shells'],
            'oxidation_state': 3, 'charge': charge,
            'Dy_4f_treatment': 'physical_sextet_OMOL__f_in_core_GFN_singlet'}


def private_carver():
    """Add deposited Nd as an explicit geometric selector in an isolated module."""
    name = '_lanm_private_source_carver'
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name('carve_generic.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    module.METAL_ELEMENTS = module.METAL_ELEMENTS | {'ND'}
    return module


def heavy_atoms(source):
    import gemmi
    from coordination_policy import selected_residue_atoms
    st = gemmi.read_structure(str(source))
    if len(st) != 1:
        raise InvalidArtifact('one experimental source model required')
    result = {}
    for chain in st[0]:
        for residue in chain:
            for atom in selected_residue_atoms(residue):
                if atom.element.name in ('H', 'D'):
                    continue
                key = (chain.name, residue.seqid.num, residue.seqid.icode.strip(),
                       residue.name, atom.name.strip())
                if key in result:
                    raise InvalidArtifact('ambiguous source heavy atom')
                result[key] = (atom.element.name, *tuple(atom.pos))
    return result


def select_source_conformers(source, output):
    """Make PDBFixer use the same deposited conformer as the existing carver."""
    import gemmi
    from coordination_policy import selected_residue_atoms, atom_altloc
    st = gemmi.read_structure(str(source))
    selections = []
    for chain in st[0]:
        for residue in chain:
            chosen = [a.clone() for a in selected_residue_atoms(residue)]
            if len(chosen) != len(residue):
                selections.append({'chain': chain.name, 'resnum': residue.seqid.num,
                    'resname': residue.name, 'selected': [{'name': a.name, 'altloc': atom_altloc(a),
                    'occupancy': a.occ, 'xyz_A': list(a.pos)} for a in chosen]})
            while len(residue):
                del residue[0]
            for atom in chosen:
                atom.altloc = '\x00'
                residue.add_atom(atom)
    st.write_pdb(str(output))
    if heavy_atoms(source) != heavy_atoms(output):
        raise InvalidArtifact('selected-conformer serialization changes source atoms')
    return {'source': record(source), 'selected_source': record(output), 'selections': selections,
            'policy': 'existing_occupancy_aware_selected_residue_atoms',
            'implementation': record(Path(__file__).with_name('coordination_policy.py'))}


def prepare_sources(root, agreement, output):
    from affordable_peptide import repair
    from protonate_cif import protonate
    root = Path(root).resolve()
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    topology = Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml')
    hans = root / 'workspaces/benchmark_set_20260915/prepared/hans_lanm_v1'
    mex = root / 'workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif'
    report = {'protocol_id': PROTOCOL, 'agreement': record(agreement), 'root': str(root),
              'implementation': record(__file__), 'topology': record(topology), 'cases': [],
              'molecular_endpoint_calls': 0, 'protonation_pH': 5.,
              'Mex_source': record(mex), 'Hans_existing_preparation': record(hans / 'preparation_report.json')}
    for i, site in enumerate(SITES):
        p = hans / site / 'amide_v3/repair_manifest.json'
        m = read_json(p)
        report['cases'].append({'case_id': 'Hans_' + site, 'protein': 'Hans', 'site': site,
            'pdb': '8DQ2', 'conditioning_metal': 'La', 'chain': 'A', 'metal_resnum': 201 + i,
            'status': 'prepared', 'repair_manifest': record(p),
            'source': m['source_structure'], 'explicit_waters': m['explicit_water_inventory']})
    try:
        selected = out / 'mex_selected_source.pdb'
        report['Mex_conformer_selection'] = select_source_conformers(mex, selected)
        pp = out / 'mex_protonated.pdb'
        report['Mex_protonation'] = protonate(selected, pp, ph=5., add_missing_residues=False)
        before, after = heavy_atoms(mex), heavy_atoms(pp)
        missing = sorted(set(before) - set(after))
        shifts = [max(abs(a - b) for a, b in zip(before[k][1:], after[k][1:]))
                  for k in before if k in after]
        report['source_heavy_check'] = {'missing': missing, 'max_displacement_A': max(shifts, default=0.),
            'added': [list(k) for k in sorted(set(after) - set(before))]}
        if missing or max(shifts, default=0.) > 1e-8:
            raise InvalidArtifact('Mex deposited heavy atoms changed during preparation')
        carver = private_carver()
        for i, site in enumerate(SITES):
            row = {'case_id': 'Mex_' + site, 'protein': 'Mex', 'site': site, 'pdb': '8FNS',
                   'conditioning_metal': 'Nd', 'chain': 'A', 'metal_resnum': 201 + i,
                   'status': 'unsupported', 'source': record(pp)}
            try:
                d = out / ('Mex_' + site)
                stem = 'mex_lanm_' + site
                carver.carve(pp, d / 'selection_intermediate', stem, site_chain='A',
                    site_resnum=201+i, site_icode='', site_atom='ND', qm_inclusion_cut=3.3)
                mp = d / 'selection_intermediate' / (stem + '_carve_manifest.json')
                repaired = repair(mp, d / 'amide_v3', topology)
                # Added heavy atoms may exist at remote termini, never in a scored fragment.
                for atom in repaired['atom_graph']['source_to_qm']:
                    if atom['kind'] != 'source' or atom['source']['element'] in ('H', 'D'):
                        continue
                    s = atom['source']
                    k = (s['chain'], s['resnum'], s['insertion_code'], s['resname'], s['atom'])
                    if k not in before:
                        raise InvalidArtifact('selected heavy atom was missing from experimental source')
                row.update(status='prepared', repair_manifest=record(d/'amide_v3/repair_manifest.json'),
                    explicit_waters=repaired['explicit_water_inventory'])
            except Exception as exc:
                row.update(reason=f'{type(exc).__name__}: {exc}')
            report['cases'].append(row)
    except Exception as exc:
        report['Mex_preparation_error'] = f'{type(exc).__name__}: {exc}'
        represented = {r['case_id'] for r in report['cases']}
        for site in SITES:
            if 'Mex_' + site not in represented:
                report['cases'].append({'case_id': 'Mex_' + site, 'protein': 'Mex', 'site': site,
                    'status': 'unsupported', 'reason': report['Mex_preparation_error']})
    report['wall_seconds'] = time.monotonic() - tick
    report['supported'] = sum(r['status'] == 'prepared' for r in report['cases'])
    report['denominator'] = 6
    write_new(out / 'sources.json', report)
    return report


def prepare(sources, agreement, software, parameter_export, output):
    import mace_omol
    source = read_json(sources)
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    params = read_json(parameter_export)
    if params['meta']['name'] != 'GFN2-xTB':
        raise InvalidArtifact('native GFN2 parameter export required')
    model = mace_omol.model(software)
    inspection = read_json(verify(read_json(software)['inspection']))
    if not set(Z.values()) <= set(inspection['elements']):
        raise InvalidArtifact('checkpoint lacks required element')
    impl = out / 'implementation'
    impl.mkdir()
    pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target = impl / p.name
        shutil.copyfile(p, target)
        pins[p.name] = record(target)
    native = out / 'native_initial'
    native.mkdir()
    tasks, endpoints, cases = [], [], []
    for c in source['cases']:
        row = {**c, 'endpoints': {}}
        cases.append(row)
        for metal in METALS:
            e = {'endpoint_id': c['case_id'] + '__' + metal, 'case_id': c['case_id'],
                 'metal': metal, 'status': 'unavailable', 'energy_eV': None}
            endpoints.append(e)
            row['endpoints'][metal] = e['endpoint_id']
            if c['status'] != 'prepared':
                e['reason'] = c.get('reason', 'source preparation unavailable')
                continue
            repair = read_json(verify(c['repair_manifest']))
            old = repair['outputs']['La']
            atoms = xyz(verify(old['xyz']))
            atoms[0] = (metal, *atoms[0][1:])
            charge = old['charge']
            state_record = state(atoms, charge, metal, params)
            xp = out / 'endpoints' / e['endpoint_id'] / 'core.xyz'
            xp.parent.mkdir(parents=True)
            # Preserve every serialized coordinate token; only the metal label changes.
            lines = verify(old['xyz']).read_text().splitlines()
            lines[2] = metal + lines[2][len('La'):]
            xp.write_text('\n'.join(lines) + '\n')
            e.update(status='prepared', xyz=record(xp), charge=charge,
                     spin_multiplicity=PHYSICAL_MULTIPLICITY[metal],
                     state=state_record, source_preparation=c['repair_manifest'])
            e['cache_key'] = cache_key({'protocol': PROTOCOL, 'xyz': e['xyz'], 'state': state_record,
                                        'model': model, 'implementation': pins['lanm_series_followup.py']})
            for medium in ('vacuum', 'alpb'):
                tid = e['endpoint_id'] + '__' + medium
                d = native / 'tasks' / tid
                d.mkdir(parents=True)
                shutil.copyfile(xp, d/'core.xyz')
                (d/'endpoint.inp').write_text(recipe(charge, medium))
                tasks.append({'task_id': tid, 'endpoint_id': e['endpoint_id'], 'case_id': c['case_id'],
                    'metal': metal, 'medium': medium, 'charge': charge, 'multiplicity': 1,
                    'physical_multiplicity': PHYSICAL_MULTIPLICITY[metal], 'state': state_record,
                    'xyz': record(d/'core.xyz'), 'input': record(d/'endpoint.inp'),
                    'output_path': str(d/'endpoint.out')})
    release = read_json(Path(__file__).resolve().parents[1]/'params/baseline_water_v1.json')['artifacts']
    low = {'protocol_id': PROTOCOL, 'stage': 'initial', 'tasks': tasks,
        'orca': release['orca'], 'implementation': pins, 'parameter_export': record(parameter_export),
        'agreement': record(agreement),
        'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8},
        'execution_policy': {'task_runner': pins['run_orca_task_manifest.py'],
                             'runtime_renderer': pins['render_orca_runtime_input.py']}}
    write_new(native/'manifest.json', low)
    m = {'protocol_id': PROTOCOL, 'sources': record(sources), 'agreement': record(agreement),
         'software': record(software), 'model': model, 'parameter_export': record(parameter_export),
         'implementation': pins, 'cases': cases, 'endpoints': endpoints,
         'native_initial_manifest': record(native/'manifest.json'),
         'capability_first': {'MACE': 'Hans_EF1__Dy', 'GFN': 'Hans_EF1__Dy__vacuum'},
         'declared_counts': {'sites': 6, 'MACE': 12, 'GFN_cells': 24, 'GFN_calls': 48},
         'native_continuation_manifest': None, 'status': 'prepared_not_authorized_for_execution',
         'affinity_reference': None, 'absolute_affinity': None, 'baseline_changed': False}
    write_new(out/'manifest.json', m)
    result = validate(out/'manifest.json')
    write_new(out/'PREFLIGHT.json', result)
    return result


def validate(manifest):
    from run_orca_task_manifest import load_manifest_tasks
    m = read_json(manifest)
    verify(m['sources']); verify(m['agreement']); verify(m['software'])
    p = read_json(verify(m['parameter_export']))
    for pin in m['implementation'].values():
        verify(pin)
    if m['protocol_id'] != PROTOCOL or len(m['endpoints']) != 12:
        raise InvalidArtifact('scope differs')
    if {(c['protein'], c['site']) for c in m['cases']} != {(p,s) for p in ('Hans','Mex') for s in SITES}:
        raise InvalidArtifact('all ordered primary sites required')
    eps = {e['endpoint_id']: e for e in m['endpoints']}
    for c in m['cases']:
        pair = [eps[c['endpoints'][z]] for z in METALS]
        if c['status'] != 'prepared':
            if any(e['status'] != 'unavailable' for e in pair):
                raise InvalidArtifact('unsupported source became available')
            continue
        atoms = [xyz(verify(e['xyz'])) for e in pair]
        if atoms[0][1:] != atoms[1][1:] or atoms[0][0][1:] != atoms[1][0][1:] or pair[0]['charge'] != pair[1]['charge']:
            raise InvalidArtifact('metal pair changes nonmetal coordinates or charge')
        for e, rows in zip(pair, atoms):
            if state(rows, e['charge'], e['metal'], p) != e['state']:
                raise InvalidArtifact('electronic-state record differs')
    lm, tasks = load_manifest_tasks(verify(m['native_initial_manifest']))
    for t in lm['tasks']:
        e = eps[t['endpoint_id']]
        if verify(t['xyz']).read_bytes() != verify(e['xyz']).read_bytes() or verify(t['input']).read_text() != recipe(e['charge'], t['medium']):
            raise InvalidArtifact('native source/recipe mismatch')
    return {'status': 'prepared_not_executed', 'manifest': record(manifest),
        'supported_sites': sum(c['status'] == 'prepared' for c in m['cases']),
        'site_denominator': 6, 'prepared_MACE': sum(e['status'] == 'prepared' for e in m['endpoints']),
        'MACE_denominator': 12, 'initial_native_tasks': len(tasks),
        'continuation_native_tasks_pending_actual_seed_outputs': len(tasks),
        'new_endpoint_calls': 0, 'execution_gate': 'awaiting_scientific_choice_and_state_policy_agreement'}


def mace(manifest, authorization, output, task=None):
    """Native warm evaluator; physical Dy sextet, no shared Ca/La guard changes."""
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():
        raise InvalidArtifact('allocated GPU required')
    validate(manifest)
    auth = record(authorization)
    m = read_json(manifest)
    if record(__file__)['sha256'] != m['implementation']['lanm_series_followup.py']['sha256']:
        raise InvalidArtifact('execute the pinned adapter')
    chosen = set(task or [e['endpoint_id'] for e in m['endpoints'] if e['status']=='prepared'])
    known = {e['endpoint_id'] for e in m['endpoints'] if e['status']=='prepared'}
    if not chosen <= known or not chosen:
        raise InvalidArtifact('unknown/unprepared selected endpoint')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']))
    torch.set_default_dtype(torch.float64)
    torch.cuda.reset_peak_memory_stats()
    tick = time.monotonic()
    calc = mace_omol(model=str(verify(m['model']['checkpoint'])), device='cuda', default_dtype='float64')
    model = calc.models[0]
    if (type(model).__name__ != 'ScaleShiftMACE' or dict(model.embedding_specs) != m['model']['embedding_specs']
            or list(model.heads) != ['omol'] or calc.energy_units_to_eV != 1.):
        raise InvalidArtifact('native model/head/units differ')
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    model_load_seconds = time.monotonic() - tick
    rows = []
    ordered = sorted((e for e in m['endpoints'] if e['endpoint_id'] in chosen),
                     key=lambda e: (e['endpoint_id'] != m['capability_first']['MACE'], e['endpoint_id']))
    for e in ordered:
        start = time.monotonic()
        r = {'endpoint_id': e['endpoint_id'], 'xyz': e['xyz'], 'state': e['state'],
             'cache_key': e['cache_key'], 'energy_eV': None, 'status': 'failed', 'model_call_attempted': False}
        try:
            a = xyz(verify(e['xyz']))
            if not {Z[x[0]] for x in a} <= set(model.atomic_numbers.tolist()):
                raise InvalidArtifact('required element absent from loaded checkpoint')
            atoms = Atoms([x[0] for x in a], positions=[x[1:] for x in a], pbc=False)
            atoms.info.update(charge=e['charge'], spin=e['spin_multiplicity'])
            batch = calc._atoms_to_batch(atoms)
            r['input_state_check'] = input_batch(calc, atoms, e['charge'], e['spin_multiplicity'], batch, 0)
            with torch.no_grad():
                r['model_call_attempted'] = True
                v = model(batch.to_dict(), training=False, compute_force=False, compute_virials=False,
                    compute_stress=False, compute_displacement=False, compute_hessian=False,
                    compute_edge_forces=False, compute_atomic_stresses=False)
            if v['forces'] is not None:
                raise InvalidArtifact('unexpected force calculation')
            value = float(v['energy'].detach().cpu().item())
            if not math.isfinite(value):
                raise InvalidArtifact('nonfinite native energy')
            r.update(status='complete', energy_eV=value)
        except Exception as exc:
            r['reason'] = f'{type(exc).__name__}: {exc}'
        r['wall_seconds'] = time.monotonic() - start
        rows.append(r)
        write_new(out/(e['endpoint_id']+'.json'), r)
        if r['status'] != 'complete':
            break  # Unsupported Dy capability is never approximated away.
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'authorization': auth,
              'selected_endpoints': sorted(chosen), 'rows': rows,
              'molecular_calls_attempted': sum(r['model_call_attempted'] for r in rows),
              'model_load_seconds': model_load_seconds, 'wall_seconds': time.monotonic()-tick,
              'job_id': os.environ['SLURM_JOB_ID'], 'allocated_CPUs': os.environ['SLURM_CPUS_PER_TASK'],
              'peak_host_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'peak_cuda_allocated_bytes': torch.cuda.max_memory_allocated(),
              'device': torch.cuda.get_device_name(),
              'versions': {k: importlib.metadata.version(k) for k in ('torch','mace-torch','ase')}}
    write_new(out/'result.json', result)
    return result


def prepare_continuation(manifest, authorization, output, task=None):
    """Stage exact own-cell seeds after real successful initial calculations."""
    from compact_solvation import completed, diagnostics
    validate(manifest)
    m = read_json(manifest)
    initial_path = verify(m['native_initial_manifest'])
    initial = read_json(initial_path)
    chosen = set(task or [t['task_id'] for t in initial['tasks']])
    if not chosen or not chosen <= {t['task_id'] for t in initial['tasks']}:
        raise InvalidArtifact('unknown continuation selection')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    tasks, unavailable = [], []
    for t in initial['tasks']:
        if t['task_id'] not in chosen:
            continue
        original = completed(initial_path, t['task_id'])
        if original is None:
            unavailable.append({'task_id': t['task_id'], 'reason': 'initial cell unavailable'})
            continue
        audit = diagnostics(original, t)
        if audit['native_valence_electrons'] != t['state']['GFN_valence_electron_count']:
            raise InvalidArtifact('actual native effective electron count differs')
        if audit['charge_sanity_status'] != 'pass':
            unavailable.append({'task_id': t['task_id'], 'reason': 'initial charge sanity fails'})
            continue
        prior = verify(original['output']).parent
        if not all((prior / ('endpoint.runtime.'+suffix)).is_file() for suffix in ('xtbw','gbw')):
            unavailable.append({'task_id': t['task_id'], 'reason': 'initial matching restart files missing'})
            continue
        seed = record(prior/'endpoint.runtime.xtbw')
        gbw = record(prior/'endpoint.runtime.gbw')
        d = out/'tasks'/t['task_id']
        d.mkdir(parents=True)
        shutil.copyfile(verify(t['xyz']), d/'core.xyz')
        shutil.copyfile(verify(seed), d/'seed.immutable.xtbw')
        shutil.copyfile(verify(gbw), d/'seed.immutable.gbw')
        (d/'endpoint.inp').write_text(recipe(t['charge'], t['medium'], continuation=True))
        tasks.append({**t, 'xyz': record(d/'core.xyz'), 'input': record(d/'endpoint.inp'),
            'output_path': str(d/'endpoint.out'), 'original': original, 'initial_audit': audit,
            'immutable_seed': record(d/'seed.immutable.xtbw'),
            'immutable_gbw': record(d/'seed.immutable.gbw'),
            'seed_source': seed, 'gbw_source': gbw})
    continuation = {**initial, 'stage': 'self_continuation', 'parent': record(manifest),
        'authorization': record(authorization), 'tasks': tasks,
        'selected_cells': sorted(chosen), 'unavailable': unavailable}
    write_new(out/'manifest.json', continuation)
    from run_orca_task_manifest import load_manifest_tasks
    if tasks:
        load_manifest_tasks(out/'manifest.json')
    return {'manifest': record(out/'manifest.json'), 'prepared': len(tasks), 'unavailable': unavailable,
            'new_molecular_calls': 0}


def execute_native(manifest, authorization, task=None):
    import fcntl
    from run_orca_task_manifest import run_manifest
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['stage'] not in ('initial','self_continuation'):
        raise InvalidArtifact('wrong finite native manifest')
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('allocation required')
    if record(__file__)['sha256'] != m['implementation']['lanm_series_followup.py']['sha256']:
        raise InvalidArtifact('execute the pinned adapter')
    auth = record(authorization)
    out = Path(manifest).resolve().parent
    chosen = set(task or [t['task_id'] for t in m['tasks']])
    if not chosen or not chosen <= {t['task_id'] for t in m['tasks']}:
        raise InvalidArtifact('unknown native task')
    with (out/'execute.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        seeds = []
        for t in m['tasks']:
            if t['task_id'] not in chosen:
                continue
            verify(t['xyz']); verify(t['input'])
            if state(xyz(verify(t['xyz'])), t['charge'], t['metal'],
                     read_json(verify(m['parameter_export']))) != t['state']:
                raise InvalidArtifact('native effective state differs')
            if verify(t['input']).read_text() != recipe(t['charge'], t['medium'], m['stage']=='self_continuation'):
                raise InvalidArtifact('native recipe differs')
            if m['stage'] == 'self_continuation':
                d = Path(t['output_path']).parent
                # A completed continuation is immutable; never seed or rerun it.
                if Path(t['output_path']).exists():
                    raise InvalidArtifact('continuation already attempted; preserve it')
                for kind, suffix in (('immutable_seed','xtbw'), ('immutable_gbw','gbw')):
                    shutil.copyfile(verify(t[kind]), d/('endpoint.runtime.'+suffix))
                seeds.append({'task_id': t['task_id'], 'xtbw': record(d/'endpoint.runtime.xtbw'),
                              'gbw': record(d/'endpoint.runtime.gbw')})
        tick = time.monotonic()
        error = None
        results = None
        try:
            results = run_manifest(Path(manifest), orca_path=verify(m['orca']), workers=min(8,len(chosen)),
                                   nprocs=8, selected_tasks=chosen)
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
        result = {'manifest': record(manifest), 'authorization': auth, 'selected_cells': sorted(chosen),
            'seed_before': seeds, 'results': results, 'error': error,
            'wall_seconds': time.monotonic()-tick, 'job_id': os.environ['SLURM_JOB_ID'],
            'allocated_CPUs': os.environ.get('SLURM_CPUS_ON_NODE')}
        write_new(out/('execution_'+os.environ['SLURM_JOB_ID']+'.json'), result)
        if error:
            raise InvalidArtifact(error)
    return result


def balanced_difference(hans_la, hans_dy, mex_la, mex_dy):
    values = (hans_la, hans_dy, mex_la, mex_dy)
    if any(v is None for v in values):
        return None
    if not all(math.isfinite(v) for v in values):
        raise InvalidArtifact('nonfinite balanced exchange component')
    return (hans_dy-hans_la) - (mex_dy-mex_la)


def collect(manifest, output, mace_collection=None, continuation_manifest=None):
    from affordable_common import HA_TO_KCAL
    from compact_solvation import completed, diagnostics
    from mace_hybrid import EV_TO_KCAL
    from structure_informed_starts import scf_details
    validate(manifest)
    m = read_json(manifest)
    mp = record(manifest)
    mace_rows = {}
    for p in mace_collection or []:
        data = read_json(p)
        if data['manifest'] != mp or data['protocol_id'] != PROTOCOL:
            raise InvalidArtifact('MACE collection belongs to another preparation')
        for r in data['rows']:
            key = r['endpoint_id']
            if key in mace_rows:
                raise InvalidArtifact('duplicate MACE endpoint, not silently selected')
            mace_rows[key] = {**r, 'collection': record(p)}
    initial = verify(m['native_initial_manifest'])
    native_tasks = {t['task_id']: t for t in read_json(initial)['tasks']}
    continuations = {}
    for p in continuation_manifest or []:
        cm = read_json(p)
        if cm['parent'] != mp or cm['protocol_id'] != PROTOCOL or cm['stage'] != 'self_continuation':
            raise InvalidArtifact('continuation does not belong to prepared cells')
        for t in cm['tasks']:
            key = t['task_id']
            if key in continuations or key not in native_tasks:
                raise InvalidArtifact('duplicate/unknown continuation cell')
            original = native_tasks[key]
            for field in ('charge','multiplicity','metal','medium','state','endpoint_id'):
                if t[field] != original[field]:
                    raise InvalidArtifact('continuation state/identity differs')
            if verify(t['xyz']).read_bytes() != verify(original['xyz']).read_bytes():
                raise InvalidArtifact('continuation geometry differs')
            continuations[key] = (Path(p), t)
    rows = []
    for e in m['endpoints']:
        row = {'endpoint_id': e['endpoint_id'], 'case_id': e['case_id'], 'metal': e['metal'],
            'status': 'unavailable', 'native_MACE_eV': None, 'solvent_transfer_kcal_mol': None,
            'composite_kcal_mol': None, 'cells': {}}
        r = mace_rows.get(e['endpoint_id'])
        if r:
            row['MACE_receipt'] = r
            for key in ('xyz','state','cache_key'):
                if r[key] != e[key]:
                    raise InvalidArtifact('MACE source/state differs')
            if r['status'] == 'complete':
                check = r['input_state_check']
                if (check['charge'],check['spin_multiplicity'],check['head']) != (e['charge'],e['spin_multiplicity'],'omol'):
                    raise InvalidArtifact('actual physical MACE batch state differs')
                row['native_MACE_eV'] = r['energy_eV']
        for medium in ('vacuum','alpb'):
            tid = e['endpoint_id']+'__'+medium
            cell = {'status': 'unavailable', 'initial_hartree': None, 'continuation_hartree': None}
            row['cells'][medium] = cell
            if tid not in native_tasks:
                cell['reason'] = 'source preparation unsupported'
                continue
            t = native_tasks[tid]
            first = completed(initial, tid)
            if first:
                try:
                    first_audit = diagnostics(first,t)
                    cell.update(initial=first, initial_hartree=first['energy_hartree'], initial_audit=first_audit)
                    if tid not in continuations:
                        cell['reason'] = 'no actual self-continuation available'
                        continue
                    cp, ct = continuations[tid]
                    last = completed(cp,tid)
                    if last is None:
                        cell['reason'] = 'self-continuation failed or unrun'
                        continue
                    if ct['original'] != first:
                        raise InvalidArtifact('self-continuation original receipt differs')
                    audit = diagnostics(last,ct)
                    text = verify(last['output']).read_text()
                    if 'INITIAL GUESS: XTBRESTART' not in text:
                        raise InvalidArtifact('native self-continuation not positively activated')
                    if first_audit['parameter_export']['sha256'] != audit['parameter_export']['sha256']:
                        raise InvalidArtifact('parameters changed during continuation')
                    if audit['native_valence_electrons'] != e['state']['GFN_valence_electron_count']:
                        raise InvalidArtifact('actual f-in-core electron count differs')
                    if audit['charge_sanity_status'] != 'pass':
                        raise InvalidArtifact('native charge sanity fails')
                    cell.update(status='complete', continuation=last, continuation_hartree=last['energy_hartree'],
                        continuation_audit=audit, actual_restart_confirmed=True, scf_details=scf_details(text),
                        continuation_minus_initial_kcal_mol=(last['energy_hartree']-first['energy_hartree'])*HA_TO_KCAL)
                except (InvalidArtifact,KeyError,OSError,ValueError) as exc:
                    cell['reason'] = f'{type(exc).__name__}: {exc}'
            else:
                cell['reason'] = 'initial native cell failed or unrun'
        if all(c['status']=='complete' for c in row['cells'].values()):
            row['solvent_transfer_kcal_mol'] = (row['cells']['alpb']['continuation_hartree']
                                              - row['cells']['vacuum']['continuation_hartree'])*HA_TO_KCAL
            if row['native_MACE_eV'] is not None:
                row.update(status='complete', composite_kcal_mol=row['native_MACE_eV']*EV_TO_KCAL+row['solvent_transfer_kcal_mol'])
        rows.append(row)
    by_id = {r['endpoint_id']: r for r in rows}
    vector = []
    for site in SITES:
        quartet = [by_id[p+'_'+site+'__'+metal] for p in ('Hans','Mex') for metal in METALS]
        vector.append({'site': site, 'status': 'complete' if all(r['status']=='complete' for r in quartet) else 'unavailable',
            'D_composite_kcal_mol': balanced_difference(*(r['composite_kcal_mol'] for r in quartet)),
            'D_native_MACE_kcal_mol': balanced_difference(*(None if r['native_MACE_eV'] is None else r['native_MACE_eV']*EV_TO_KCAL for r in quartet)),
            'D_solvent_transfer_kcal_mol': balanced_difference(*(r['solvent_transfer_kcal_mol'] for r in quartet))})
    result = {'protocol_id': PROTOCOL, 'manifest': mp, 'endpoints': rows, 'ordered_vector': vector,
        'complete_endpoints': sum(r['status']=='complete' for r in rows), 'endpoint_denominator': 12,
        'new_calls_in_collection': 0, 'absolute_affinity': None, 'Kd': None, 'state_populations': None,
        'reference_cancellation': 'balanced_Hans_Mex_La_Dy_exchange',
        'positive_D_interpretation': 'Hans_more_La_selective_relative_to_Mex_in_conditional_electronic_descriptor',
        'electronic_continuity_qualification': 'not_established_by_one_self_continuation',
        'baseline_changed': False}
    write_new(output,result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='operation', required=True)
    q = sub.add_parser('prepare-sources')
    for key in ('root', 'agreement', 'output'):
        q.add_argument('--'+key, required=True)
    q = sub.add_parser('prepare')
    for key in ('sources', 'agreement', 'software', 'parameter-export', 'output'):
        q.add_argument('--'+key, required=True)
    q = sub.add_parser('validate')
    q.add_argument('--manifest', required=True)
    for op in ('mace','prepare-continuation','execute-native'):
        q = sub.add_parser(op)
        q.add_argument('--manifest', required=True)
        q.add_argument('--authorization', required=True)
        q.add_argument('--task', action='append')
        if op != 'execute-native':
            q.add_argument('--output', required=True)
    q = sub.add_parser('collect')
    q.add_argument('--manifest', required=True)
    q.add_argument('--output', required=True)
    q.add_argument('--mace-collection', action='append')
    q.add_argument('--continuation-manifest', action='append')
    args = vars(p.parse_args())
    operation = args.pop('operation').replace('-', '_')
    result = globals()[operation](**args)
    print(json.dumps({k:v for k,v in result.items() if k != 'cases'}, indent=2))


if __name__ == '__main__':
    main()
