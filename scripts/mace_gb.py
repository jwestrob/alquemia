"""Frozen learned-charge OBC-II descriptor; opt-in, with no calibrated class."""
from __future__ import annotations

import argparse
import copy
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

import numpy as np

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms, rotation, write_xyz

SCHEMA = 'alquemia.mace_gb.v1'
PROTOCOL = 'mace_frozen_monopole_obc2_v1'
EV_TO_KJ = 96.48533212331002
RADII = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'Ca': 1.8, 'La': 1.8}
SCALES = {'H': .85, 'C': .72, 'N': .79, 'O': .85, 'S': .96, 'Ca': .8, 'La': .8}
MODEL = {'solver': 'OpenMM_GBSAOBCForce_OBC2', 'radii_A': RADII, 'descreen_scales': SCALES,
         'solute_dielectric': 1., 'solvent_dielectric': 78.5, 'salt_M': 0., 'temperature_K': 298.15,
         'surface_coefficient': 0., 'periodic': False, 'cutoff': None, 'cuda_precision': 'double',
         'charge_representation': 'saved_MACE_density_coefficient_0_monopole',
         'dipolar_solvent_response': False, 'self_consistent_solvent_response': False,
         'energy_definition': 'GB_reaction_charging_only; direct_Coulomb_excluded',
         'force_definition': 'negative_gradient_of_GB_reaction_energy_at_fixed_monopoles',
         'metal_cavity_status': 'unvalidated_common_radius_and_generic_descreen_factor'}
TOL = {'total_charge_e': 1e-5, 'energy_kcal_mol': .01, 'force_max_eV_A': .001,
       'identity_energy_kcal_mol': 1e-6, 'identity_force_eV_A': 1e-6,
       'native_custom_energy_kcal_mol': .001}


def task_specs():
    specs = []
    for metal in ('La', 'Ca'):
        for solver in ('native', 'custom'):
            specs.append(('core', 'medium_core', metal, solver, 'Reference', 'primary'))
    for variant in ('primary', 'rotate', 'identity'):
        for model in ('medium', 'large'):
            for metal in ('La', 'Ca'):
                specs.append(('full', model, metal, 'native', 'CUDA', variant))
    for variant, platform in (('repeat', 'CUDA'), ('translate', 'CUDA'), ('reference', 'Reference')):
        specs.append(('full', 'medium', 'La', 'native', platform, variant))
    return specs


def source_data(collection_path, source_id):
    """Only admit archived successful calculations whose execution receipt agrees."""
    from mace_hybrid import accepted_attempt
    c = read_json(collection_path); mp = verify(c['manifest']); m = read_json(mp)
    t = next(t for t in m['tasks'] if t['task_id'] == source_id)
    r = c['rows'][source_id]
    candidates = [accepted_attempt(a, t, mp) for a in sorted((mp.parent/'execution'/source_id).glob('attempt_*'))]
    if r not in candidates or r['status'] != 'computed' or not r['charge_check']:
        raise InvalidArtifact('source MACE output has no matching successful receipt')
    atoms = xyz(verify(t['xyz'])); density = np.load(verify(r['density_coefficients']))
    if density.shape != (len(atoms), 4) or not np.isfinite(density).all():
        raise InvalidArtifact('invalid saved MACE density')
    if abs(float(density[:, 0].sum()) - t['charge']) > TOL['total_charge_e']:
        raise InvalidArtifact('saved monopole charge closure failed')
    return m, t, r, atoms


def prepare(medium, large, cores, agreement, output):
    import openmm
    from openmm.app.internal import customgbforces
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    implementation = out/'implementation'; implementation.mkdir()
    pins = {}
    for name in ('mace_gb.py', 'mace_hybrid.py', 'affordable_common.py'):
        shutil.copyfile(Path(__file__).with_name(name), implementation/name)
        pins[name] = record(implementation/name)
    # Pin the installed implementation without installing or altering it.
    prefix = Path(sys.prefix)
    files = [Path(openmm.__file__), Path(openmm._openmm.__file__), Path(customgbforces.__file__),
             Path(openmm.__file__).with_name('openmm.py')]
    native = Path(openmm.version.openmm_library_path)
    files += sorted(native.glob('libOpenMM.so*'))
    files += sorted((native/'plugins').glob('libOpenMM*CUDA*.so'))
    files += sorted((native/'plugins').glob('libOpenMM*Reference*.so'))
    nvidia = Path(openmm.__file__).parent.parent/'nvidia'
    files += sorted(nvidia.glob('*/lib/*.so*'))
    files = sorted(set(p.resolve() for p in files))
    write_new(out/'backend_inventory.json', {'files': [record(p) for p in files]})
    requirements = out/'requirements.txt'
    requirements.write_text('\n'.join(sorted(f'{d.metadata["Name"]}=={d.version}' for d in importlib.metadata.distributions()))+'\n')
    python_ref = record(sys.executable); python_ref['path'] = str(Path(sys.executable).absolute())
    write_new(out/'software.json', {'python': python_ref, 'requirements': record(requirements),
        'openmm_version': openmm.__version__, 'backend_source_inventory': record(out/'backend_inventory.json'),
        'available_platforms': [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())],
        'plugin_load_failures': list(openmm.Platform.getPluginLoadFailures())})
    collections = {'medium': record(medium), 'large': record(large), 'medium_core': record(cores)}
    tasks = []
    for kind, checkpoint, metal, solver, platform, variant in task_specs():
        source_id = f'1h4i_{"qm33" if kind == "core" else "full"}_{metal}'+('' if kind == 'core' else '_primary')
        sm, st, sr, atoms = source_data(verify(collections[checkpoint]), source_id)
        state = next((s for s in sm['source_states'] if read_json(verify(s))['metal'] == metal), None)
        if state is None:
            raise InvalidArtifact('source physical mapping unavailable')
        state_data = read_json(verify(state)); mapped = state_data['core_atoms' if kind == 'core' else 'physical_atoms']
        if len(atoms) != len(mapped):
            raise InvalidArtifact('source mapping length mismatch')
        for atom, physical in zip(atoms, mapped):
            element = metal if physical['element'] == 'M' else physical['element']
            if atom != (element, *physical['xyz_A']) or abs(physical['radius_A']-RADII[element]) > 1e-12:
                raise InvalidArtifact('source physical coordinate/radius mismatch')
        coords = np.array([a[1:] for a in atoms])
        center = coords[next(i for i, a in enumerate(atoms) if a[0] == metal)]
        matrix = rotation() if variant == 'rotate' else np.eye(3)
        shift = np.array([10., -7., 3.]) if variant == 'translate' else np.zeros(3)
        if variant in ('rotate', 'translate'):
            coords = (coords-center) @ matrix.T + center + shift
        name = '_'.join((checkpoint, metal, solver, variant))
        path = out/(name+'.xyz'); write_xyz(path, [(a[0], *v) for a, v in zip(atoms, coords)])
        task = {'task_id': name, 'kind': kind, 'checkpoint_label': checkpoint, 'metal': metal,
                'solver': solver, 'platform': platform, 'variant': variant,
                'charge': st['charge'], 'spin_multiplicity': st['spin_multiplicity'],
                'xyz': record(path), 'source_collection': collections[checkpoint], 'source_task_id': source_id,
                'source_state': state, 'source_density': sr['density_coefficients'],
                'source_vacuum_energy_eV': sr['energy_eV'], 'rotation_matrix': matrix.tolist(),
                'translation_A': shift.tolist(), 'solvent_dielectric': 1. if variant == 'identity' else 78.5}
        task['state'] = check_atoms(xyz(path), task['charge'])
        tasks.append(task)
    m = {'schema_version': SCHEMA, 'protocol_id': PROTOCOL, 'status': 'approved_prepared',
         'agreement': record(agreement), 'software': record(out/'software.json'), 'implementation': pins,
         'source_collections': collections, 'model': MODEL, 'tolerances': TOL, 'tasks': tasks,
         'reference': None, 'S_kcal_mol': None, 'calibrated_class': None,
         'corrected_gradient': None, 'corrected_gradient_status': 'charge_response_chain_rule_unavailable',
         'evidence_use': 'consumed_1H4I_method_development',
         'run_inventory': {'GB_energy_force_calls': 19, 'new_MACE_calls': 0, 'new_DFT_endpoints': 0},
         'compute_budget': None, 'wall_time_limit': None}
    for t in tasks:
        t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': pins})
    write_new(out/'manifest.json', m)
    return validate(out/'manifest.json')


def validate(manifest):
    m = read_json(manifest); software = read_json(verify(m['software']))
    if m['schema_version'] != SCHEMA or m['protocol_id'] != PROTOCOL or m['model'] != MODEL or m['tolerances'] != TOL:
        raise InvalidArtifact('GB protocol/model/acceptance settings changed; requires a new version')
    for ref in [m['agreement'], *m['implementation'].values(), software['python'], software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:
        verify(ref)
    expected = task_specs()
    if len(m['tasks']) != len(expected):
        raise InvalidArtifact('GB task inventory mismatch')
    source_cache = {}
    paired_inputs = {}
    for task, spec in zip(m['tasks'], expected):
        actual = tuple(task[k] for k in ('kind','checkpoint_label','metal','solver','platform','variant'))
        if actual != spec or task['task_id'] != '_'.join((spec[1],spec[2],spec[3],spec[5])):
            raise InvalidArtifact('GB task identity mismatch')
        ref = task['source_collection']; key = (ref['sha256'],task['source_task_id'])
        if ref != m['source_collections'][task['checkpoint_label']]:
            raise InvalidArtifact('source collection mismatch')
        if key not in source_cache:
            source_cache[key] = source_data(verify(ref), task['source_task_id'])
        sm, st, sr, atoms = source_cache[key]
        if task['source_state'] not in sm['source_states']:
            raise InvalidArtifact('source state mapping mismatch')
        verify(task['source_state'])
        if (task['source_density'] != sr['density_coefficients'] or task['source_vacuum_energy_eV'] != sr['energy_eV']
                or task['charge'] != st['charge'] or task['spin_multiplicity'] != st['spin_multiplicity']):
            raise InvalidArtifact('GB input changes saved electronic state')
        rows = xyz(verify(task['xyz'])); check_atoms(rows, task['charge'])
        pos = np.array([a[1:] for a in atoms]); center = pos[next(i for i,a in enumerate(atoms) if a[0] == task['metal'])]
        matrix = rotation() if task['variant'] == 'rotate' else np.eye(3)
        shift = np.array([10.,-7.,3.]) if task['variant'] == 'translate' else np.zeros(3)
        wanted = (pos-center) @ matrix.T + center + shift if task['variant'] in ('rotate','translate') else pos
        if ([r[0] for r in rows] != [r[0] for r in atoms] or
                not np.allclose(np.array([a[1:] for a in rows]),wanted,atol=1e-12,rtol=0)
                or task['rotation_matrix'] != matrix.tolist() or task['translation_A'] != shift.tolist()):
            raise InvalidArtifact('GB coordinates/rigid transform differ from declared source')
        if task['solvent_dielectric'] != (1. if task['variant'] == 'identity' else 78.5):
            raise InvalidArtifact('task dielectric mismatch')
        base = {k:v for k,v in task.items() if k != 'cache_key'}
        if task['cache_key'] != cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('GB scientific cache key mismatch')
        if task['kind'] == 'full' and task['variant'] == 'primary':
            paired_inputs[(task['checkpoint_label'],task['metal'])] = rows
    first = paired_inputs[('medium','La')]
    for rows in paired_inputs.values():
        if any(a[1:] != b[1:] or (a[0] != b[0] and {a[0],b[0]} != {'La','Ca'}) for a,b in zip(first,rows)):
            raise InvalidArtifact('full endpoint/checkpoint geometry differs')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest),'new_DFT_endpoints':0}


def build_system(atoms, charges, solver, epsilon):
    import openmm as mm
    from openmm.app.internal.customgbforces import GBSAOBC2Force
    if solver not in ('native','custom'):
        raise InvalidArtifact('unsupported GB solver')
    system = mm.System(); nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    nb.setReactionFieldDielectric(1.); nb.setForceGroup(1)
    if solver == 'native':
        gb = mm.GBSAOBCForce(); gb.setNonbondedMethod(mm.GBSAOBCForce.NoCutoff)
        gb.setSoluteDielectric(1.); gb.setSolventDielectric(epsilon); gb.setSurfaceAreaEnergy(0.)
    else:
        gb = GBSAOBC2Force(solventDielectric=epsilon, soluteDielectric=1., SA=None, cutoff=None, kappa=0.)
        gb.setNonbondedMethod(mm.CustomGBForce.NoCutoff)
    gb.setForceGroup(0)
    for atom, q in zip(atoms, charges):
        element = atom[0]
        if element not in RADII:
            raise InvalidArtifact('unsupported GB element')
        # Mass is irrelevant: no integration, minimization or velocities.
        system.addParticle(1.)
        nb.addParticle(float(q), 1., 0.)
        params = [float(q), RADII[element]/10., SCALES[element]]
        gb.addParticle(*params) if solver == 'native' else gb.addParticle(params)
    if solver == 'custom':
        gb.finalize()
    system.addForce(gb); system.addForce(nb)
    return system


def worker(manifest, task_id, output, memory_mode):
    import openmm as mm
    from openmm import unit
    m = read_json(manifest); t = next(t for t in m['tasks'] if t['task_id'] == task_id)
    out = Path(output); start = time.monotonic()
    result = {'status':'unavailable','task_id':task_id,'cache_key':t['cache_key'],
              'manifest':record(manifest),'memory_mode':memory_mode, 'platform':t['platform'],
              'energy_definition':MODEL['energy_definition'], 'force_definition':MODEL['force_definition'],
              'corrected_descriptor_gradient':None, 'solver':t['solver'],
              'peak_cuda_allocated_bytes':None,'peak_cuda_reserved_bytes':None}
    try:
        if not os.environ.get('SLURM_JOB_ID'):
            raise InvalidArtifact('GB execution requires the recorded allocation')
        rows = xyz(verify(t['xyz'])); density = np.load(verify(t['source_density'])); charges = density[:,0]
        total = float(charges.sum()); error = total-t['charge']
        if density.shape != (len(rows),4) or not np.isfinite(density).all() or abs(error) > TOL['total_charge_e']:
            raise InvalidArtifact('invalid source monopoles')
        system = build_system(rows,charges,t['solver'],t['solvent_dielectric'])
        serialized = out/'system.xml'; serialized.write_text(mm.XmlSerializer.serialize(system))
        integrator = mm.VerletIntegrator(.001)
        platform = mm.Platform.getPlatformByName(t['platform'])
        properties = {'Precision':'double','DeterministicForces':'true'} if t['platform']=='CUDA' else {}
        context_start = time.monotonic(); context = mm.Context(system,integrator,platform,properties)
        result['context_seconds'] = time.monotonic()-context_start
        result['platform_properties'] = {k:platform.getPropertyValue(context,k) for k in platform.getPropertyNames()}
        context.setPositions(np.array([a[1:] for a in rows])*.1*unit.nanometer)
        tick = time.monotonic()
        state = context.getState(getEnergy=True,getForces=True,groups={0})
        value = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        forces = np.array(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))/(EV_TO_KJ*10.)
        result['evaluation_seconds'] = time.monotonic()-tick
        if not math.isfinite(value) or not np.isfinite(forces).all() or forces.shape != (len(rows),3):
            raise InvalidArtifact('nonfinite or incomplete GB output')
        fp = out/'forces_eV_A.npy'; np.save(fp,forces)
        result.update(status='computed',energy_eV=value/EV_TO_KJ, GB_reaction_kJ_mol=value,
                      GB_reaction_kcal_mol=value/4.184, forces=record(fp), density_coefficients=t['source_density'],
                      density_total_charge_e=total,total_charge_error_e=error,charge_check=True,
                      serialized_system=record(serialized),extracted_force_groups=[0],direct_Coulomb_included=False)
        del context,integrator
    except Exception as exc:
        result.update(status='failed',reason=str(exc),exception_type=type(exc).__name__); traceback.print_exc()
    finally:
        result.update(wall_seconds=time.monotonic()-start,peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      slurm_job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                      allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),versions={'openmm':mm.__version__,'numpy':np.__version__})
        write_new(out/'result.json',result)
    print(json.dumps({k:result.get(k) for k in ('task_id','status','reason','evaluation_seconds')}),flush=True)
    return result


def rows_and_attempts(manifest):
    from mace_hybrid import accepted_attempt
    mp = Path(manifest).resolve(); m = read_json(mp); rows = {}; attempts = []
    for t in m['tasks']:
        valid = []
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a,t,mp)
            if r is not None:
                verify(r['serialized_system']); valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']] = valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    return rows,attempts


def equivalence(name, a, b, energy_tol, force_tol, matrix=None):
    if a['status'] != 'computed' or b['status'] != 'computed':
        return {'name':name,'status':'unavailable','pass':False}
    delta = abs(a['GB_reaction_kcal_mol']-b['GB_reaction_kcal_mol'])
    fa,fb = (np.load(verify(r['forces'])) for r in (a,b))
    if matrix is not None:
        fb = fb @ np.array(matrix)
    force = float(np.max(np.abs(fa-fb)))
    return {'name':name,'status':'computed','energy_difference_kcal_mol':delta,'force_difference_max_eV_A':force,
            'pass':bool(delta<=energy_tol and force<=force_tol)}


def core_checks(rows):
    return [equivalence('native_custom_'+metal, rows[f'medium_core_{metal}_native_primary'],
              rows[f'medium_core_{metal}_custom_primary'],TOL['native_custom_energy_kcal_mol'],TOL['force_max_eV_A'])
            for metal in ('La','Ca')]


def core_gate(manifest):
    checks = core_checks(rows_and_attempts(manifest)[0])
    return {'status':'pass' if all(c['pass'] for c in checks) else 'fail','checks':checks}


def collect_gb(manifest):
    m = read_json(manifest); rows, attempts = rows_and_attempts(manifest)
    checks = core_checks(rows); direct = {}
    for model in ('medium','large'):
        primaries = [rows[f'{model}_{metal}_native_primary'] for metal in ('La','Ca')]
        for metal, primary in zip(('La','Ca'),primaries):
            checks.append(equivalence(model+'_'+metal+'_rotation',primary,rows[f'{model}_{metal}_native_rotate'],
                                     TOL['energy_kcal_mol'],TOL['force_max_eV_A'],rotation()))
            identity = rows[f'{model}_{metal}_native_identity']
            if identity['status'] == 'computed':
                val = abs(identity['GB_reaction_kcal_mol']); force = float(np.max(np.abs(np.load(verify(identity['forces'])))))
                checks.append({'name':model+'_'+metal+'_identity','status':'computed','energy_kcal_mol':val,
                               'force_max_eV_A':force,'pass':val<=TOL['identity_energy_kcal_mol'] and force<=TOL['identity_force_eV_A']})
            else:
                checks.append({'name':model+'_'+metal+'_identity','status':'unavailable','pass':False})
        collection = read_json(verify(m['source_collections'][model]))
        if all(p['status']=='computed' for p in primaries):
            correction = primaries[1]['GB_reaction_kcal_mol']-primaries[0]['GB_reaction_kcal_mol']
            raw = collection['direct']['R_kcal_mol']
            direct[model] = {'MACE_vacuum_R_kcal_mol':raw,'GB_Ca_minus_La_kcal_mol':correction,
                             'descriptor_R_kcal_mol':raw+correction,'S_kcal_mol':None,'calibrated_class':None}
            rotated = [rows[f'{model}_{metal}_native_rotate'] for metal in ('La','Ca')]
            if all(r['status']=='computed' for r in rotated):
                difference = abs(rotated[1]['GB_reaction_kcal_mol']-rotated[0]['GB_reaction_kcal_mol']-correction)
                checks.append({'name':model+'_contrast_rotation','status':'computed',
                               'difference_kcal_mol':difference,'pass':difference<=TOL['energy_kcal_mol']})
        else:
            direct[model] = {'status':'unavailable','descriptor_R_kcal_mol':None,'S_kcal_mol':None,'calibrated_class':None}
    for variant in ('repeat','translate','reference'):
        checks.append(equivalence('medium_La_'+variant,rows['medium_La_native_primary'],rows[f'medium_La_native_{variant}'],
                                 TOL['energy_kcal_mol'],TOL['force_max_eV_A']))
    complete = all(r['status']=='computed' for r in rows.values())
    result = {'status':'complete' if complete else 'incomplete','protocol_id':PROTOCOL,
              'manifest':record(manifest),'collection_implementation':record(__file__),'rows':rows,'attempts':attempts,
              'checks':checks,'numerical_checks_pass':complete and all(c['pass'] for c in checks),'direct':direct,
              'reference':None,'S_kcal_mol':None,'calibrated_class':None,'corrected_gradient':None,
              'source_state_comparison':'matched_corrected_full_geometry_and_frozen_saved_densities',
              'checkpoint_disagreement_kcal_mol':None}
    if all(d.get('descriptor_R_kcal_mol') is not None for d in direct.values()):
        result['checkpoint_disagreement_kcal_mol'] = {
            key: direct['large'][key]-direct['medium'][key] for key in
            ('MACE_vacuum_R_kcal_mol','GB_Ca_minus_La_kcal_mol','descriptor_R_kcal_mol')}
    return result


def report_gb(collection, output):
    c = read_json(collection)
    lines = ['# Frozen MACE monopole / OBC-II solvent check','',f"Status: {c['status']}. Numerical checks pass: {c['numerical_checks_pass']}.",'',
             'Baseline unchanged. Consumed 1H4I development inputs; zero new MACE/DFT evaluations.',
             'G is GB reaction energy only. No duplicate direct Coulomb or CPCM contribution.',
             'Raw descriptor R = MACE vacuum R + G(Ca) - G(La); all quantities below in kcal/mol.','',
             '```json',json.dumps({'direct':c['direct'],'checkpoint_disagreement':c['checkpoint_disagreement_kcal_mol']},indent=2),'```','',
             '## Numerical checks','','```json',json.dumps(c['checks'],indent=2),'```','',
             '## Limits','',
             'Frozen monopoles omit dipolar and self-consistent solvent response. Metal cavity factors are unvalidated.',
             'GB forces hold charges fixed; they are not gradients of the combined descriptor.',
             'No compatible aquo reference, calibrated S/class, or predictive-accuracy claim.',
             f"Recorded attempts: {len(c['attempts'])}; successful tasks: {sum(r['status']=='computed' for r in c['rows'].values())}/{len(c['rows'])}.",
             'Whole-job costs are recorded separately from evaluation timings.']
    with Path(output).open('x') as f:
        f.write('\n'.join(lines)+'\n')
    return {'status':'report_written','report':record(output)}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('medium-collection','large-collection','core-collection','agreement','output'):
        p.add_argument('--'+key,required=True)
    a = p.parse_args()
    print(json.dumps(prepare(a.medium_collection,a.large_collection,a.core_collection,a.agreement,a.output),indent=2))
