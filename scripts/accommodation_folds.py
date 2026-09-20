"""Frozen PQQ fold-conditioning experiment using existing preparation/workers.

Research orchestration only. No score fitting or production change.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import json
import os
from pathlib import Path
import resource
import shutil
import time

from affordable_common import InvalidArtifact, cache_key, paired, read_json, record, verify, write_new, xyz
import mace_omol as omol

SCHEMA = 'PQQ_fold_conditioning_robustness_v1'
ADAPTER = 'explicit_Ca_source_residue_alias_v1'


def atoms(structure):
    return [(c.name, r.seqid.num, str(r.seqid.icode), a.name, a.element.name,
             str(a.altloc), a.pos.x, a.pos.y, a.pos.z)
            for c in structure[0] for r in c for a in r]


def adapt_ca(source, output):
    """Residue alias only; preserve native element, atom identity and coordinates."""
    import gemmi
    raw = verify(source['source_structure']); st = gemmi.read_structure(str(raw))
    wanted = source['raw_source_metal']; selected = []
    if len(st) != 1 or wanted['element'] != 'Ca':
        raise InvalidArtifact('adapter requires one model and explicit native Ca')
    before = atoms(st)
    for chain in st[0]:
        for res in chain:
            if (chain.name, res.seqid.num, res.seqid.icode.strip(), res.name) != (
                    wanted['chain'], wanted['resnum'], wanted['icode'], wanted['resname']):
                continue
            if len(res) != 1 or res[0].name != wanted['atom'] or res[0].element.name != 'Ca':
                raise InvalidArtifact('selected Ca residue is not the declared monatomic source')
            res.name = 'LA'; selected.append(res)
    if len(selected) != 1:
        raise InvalidArtifact('ambiguous or missing native Ca source')
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    path = out / 'source_residue_alias.cif'; st.make_mmcif_document().write_file(str(path))
    if atoms(gemmi.read_structure(str(path))) != before:
        raise InvalidArtifact('source alias changed atom identity or coordinates')
    adapted = copy.deepcopy(source)
    adapted.update(source_structure=record(path), original_source_structure=source['source_structure'],
                   source_adapter=ADAPTER)
    write_new(out / 'adapter.json', {'protocol_id': ADAPTER, 'original': source['source_structure'],
              'adapted': record(path), 'selector': wanted, 'target_residue_name': 'LA',
              'all_atom_identities_elements_coordinates_unchanged': True})
    return adapted


def prepare_one(source, config, directory):
    import pqq_fast_prepare as source_prep
    out = Path(directory); source_used = source
    try:
        if source['source_conditioning_metal'] == 'Ca':
            source_used = adapt_ca(source, out / 'adapters' / source['case_id'])
        row = source_prep.prepare_case(source_used, config, out / 'cases' / source['case_id'])
        for rep in ('core', 'context'):
            eps = row['core']['endpoints'] if rep == 'core' else row['representations'][rep]['endpoints']
            paired(verify(eps['La']['xyz']), verify(eps['Ca']['xyz']), eps['La']['charge'], eps['Ca']['charge'])
        return row
    except Exception as exc:
        row = {'case_id': source['case_id'], 'source': source, 'status': 'unsupported',
               'reason': str(exc), 'exception_type': type(exc).__name__}
        write_new(out / 'failures' / (source['case_id'] + '.json'), row)
        return row


def prepare(sources, agreement, output, workers=16):
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('parallel preparation requires allocation')
    workers = int(workers)
    if workers < 1 or workers > int(os.environ['SLURM_CPUS_ON_NODE']):
        raise InvalidArtifact('preparation process count exceeds allocation')
    m = read_json(sources); verify(record(agreement))
    ids = [s['case_id'] for s in m['cases']]
    if len(ids) != 250 or len(set(ids)) != 250:
        raise InvalidArtifact('declared all250 input inventory differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic(); rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(prepare_one, s, m['config'], str(out)) for s in m['cases']]
        for job in as_completed(jobs):
            row = job.result(); rows.append(row)
            print(json.dumps({k: row.get(k) for k in ('case_id', 'status', 'reason')}), flush=True)
    rows.sort(key=lambda r: ids.index(r['case_id']))
    result = {'schema_version': SCHEMA, 'source_manifest': record(sources), 'agreement': record(agreement),
              'config': m['config'], 'cases': rows, 'denominator': len(ids),
              'supported': sum(r['status'] == 'prepared' for r in rows),
              'wall_seconds': time.monotonic() - start, 'workers': workers,
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'source_adapter': ADAPTER}
    write_new(out / 'preparation.json', result)
    return {k: result[k] for k in ('supported', 'denominator', 'wall_seconds')}


def reconcile(preparation, output):
    """Retain cross-host floating-point replay noise, without moving any atoms."""
    p = read_json(preparation); p['original_preparation'] = record(preparation)
    fixed = []
    for row in p['cases']:
        if row['status'] != 'archive_coordinate_mismatch':
            continue
        comparisons = row['context_comparison']
        if (not row['source']['canonical_coordinate_match'] or not row['normalized_bytes_equal_archive']
                or set(comparisons) != {'Ca', 'La'}):
            continue
        if all(c['charge_identical'] and c.get('max_displacement_A') is not None
               and c['max_displacement_A'] <= 1e-12 for c in comparisons.values()):
            row['original_status'] = row['status']; row['status'] = 'prepared'
            row['coordinate_replay_status'] = 'floating_point_equal_within_1e-12_A'
            fixed.append(row['case_id'])
    p['supported'] = sum(r['status'] == 'prepared' for r in p['cases'])
    p['numeric_reconciliation'] = {'case_ids': fixed, 'tolerance_A': 1e-12,
        'reason': 'Machine last-bit cap arithmetic far below six-decimal core and three-decimal source precision; no coordinate edits or score inspection.'}
    write_new(output, p)
    return {'supported': p['supported'], 'reconciled': len(fixed)}


def prepare_mace(preparation, agreement, output, reuse_inventory=None):
    import pqq_fast_prepare as source_prep
    p = read_json(preparation); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    implementation = {}; dest = out / 'implementation'; dest.mkdir()
    for script in Path(__file__).parent.glob('*.py'):
        shutil.copyfile(script, dest / script.name); implementation[script.name] = record(dest / script.name)
    m = {'schema_version': SCHEMA, 'protocol_id': omol.PROTOCOL, 'stage': 'fold_conditioning_native',
         'preparation': record(preparation), 'agreement': record(agreement),
         'software': p['config']['software'], 'model': p['config']['model'],
         'implementation': implementation, 'tasks': [], 'reused': {}, 'cache_reuse': bool(reuse_inventory)}
    old = read_json(reuse_inventory) if reuse_inventory else None
    if old and (old['model'] != m['model'] or old['software'] != m['software']):
        raise InvalidArtifact('reuse checkpoint/software differs')
    if old:
        m['reuse_inventory'] = record(reuse_inventory)
    for row in p['cases']:
        if row['status'] != 'prepared':
            continue
        for rep in ('core', 'context'):
            eps = row['core']['endpoints'] if rep == 'core' else row['representations'][rep]['endpoints']
            for metal, ep in eps.items():
                task = {'task_id': row['case_id'] + '__' + rep + '__' + metal, 'case_id': row['case_id'],
                        'representation': rep, 'metal': metal, 'metal_index': 0, 'kind': 'core',
                        'xyz': ep['xyz'], 'charge': ep['charge'], 'spin_multiplicity': ep['multiplicity'],
                        'energy_component': omol.COMPONENT}
                task['cache_key'] = cache_key({'task': task, 'model': m['model'], 'software': m['software'],
                                              'implementation': implementation})
                m['tasks'].append(task)
                if old and row['source']['canonical_coordinate_match']:
                    previous = next(c for c in old['cases'] if c['case_id'] == row['source']['root_case_id'])
                    endpoint = previous['representations'][rep]['endpoints'][metal]
                    actual = xyz(verify(ep['xyz'])); archived = xyz(verify(endpoint['xyz']))
                    if len(actual) != len(archived) or any(a[0] != b[0] or any(abs(x-y) > 1e-12
                            for x, y in zip(a[1:], b[1:])) for a, b in zip(actual, archived)):
                        raise InvalidArtifact('canonical reuse coordinates exceed machine precision tolerance')
                    if ep['charge'] != endpoint['charge'] or ep['multiplicity'] != endpoint['multiplicity']:
                        raise InvalidArtifact('canonical reuse charge/spin differs')
                    from compact_solvation_compare import native_endpoint
                    native_endpoint(read_json(verify(endpoint['native_MACE_receipt'])), m['model'])
                    m['reused'][task['task_id']] = {'receipt': endpoint['native_MACE_receipt'],
                        'new_xyz': ep['xyz'], 'archived_xyz': endpoint['xyz'], 'coordinate_tolerance_A': 1e-12,
                        'source': record(reuse_inventory), 'model_state_and_coordinates_verified': True}
    m['verification_references'] = {}
    for task_id in list(m['reused'])[:2]:
        m['verification_references'][task_id] = m['reused'].pop(task_id)
    write_new(out / 'manifest.json', m)
    return {'tasks': len(m['tasks']) - len(m['reused']), 'reused': len(m['reused']), 'manifest': record(out / 'manifest.json')}


def run_mace(manifest, shard=0, shards=1):
    """Warm native ASE calculator; a fresh graph and cache reset for each endpoint."""
    import importlib.metadata
    import numpy as np
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():
        raise InvalidArtifact('allocated GPU required')
    m = read_json(manifest); verify(m['preparation']); verify(m['agreement'])
    for pin in m['implementation'].values(): verify(pin)
    if int(shards) < 1 or not 0 <= int(shard) < int(shards):
        raise InvalidArtifact('invalid shard')
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK'])); torch.set_default_dtype(torch.float64)
    torch.cuda.reset_peak_memory_stats(); start = time.monotonic()
    calc = mace_omol(model=str(verify(m['model']['checkpoint'])), device='cuda', default_dtype='float64')
    model = calc.models[0]
    if (type(model).__name__ != 'ScaleShiftMACE' or dict(model.embedding_specs) != m['model']['embedding_specs']
            or list(model.heads) != ['omol'] or calc.energy_units_to_eV != 1.):
        raise InvalidArtifact('native checkpoint/state embedding differs')
    for parameter in model.parameters(): parameter.requires_grad_(False)
    versions = {n: p._version for n, p in model.named_parameters()}; load_seconds = time.monotonic() - start
    tasks = sorted(m['tasks'], key=lambda t: t['task_id'] not in m['verification_references'])
    manifest_pin = record(manifest); receipts = []
    for i, task in enumerate(tasks):
        if i % int(shards) != int(shard):
            continue
        if task['task_id'] in m.get('reused', {}):
            verify(m['reused'][task['task_id']]['receipt']); continue
        out = Path(manifest).parent / 'results' / task['task_id']
        if (out / 'result.json').exists():
            old = read_json(out / 'result.json')
            if old.get('manifest') == record(manifest) and old.get('cache_key') == task['cache_key']:
                continue  # Failed receipts remain failed; retry is a separate attempt.
            raise InvalidArtifact('existing receipt does not describe this task')
        out.mkdir(parents=True, exist_ok=False)
        tick = time.monotonic()
        result = {'task_id': task['task_id'], 'cache_key': task['cache_key'], 'manifest': manifest_pin,
                  'status': 'failed', 'memory_mode': 'native', 'energy_component': omol.COMPONENT,
                  'density_coefficients': None, 'model_load_reused': True, 'energy_eV': None}
        try:
            rows = xyz(verify(task['xyz'])); omol.check_atoms(rows, task['charge'])
            a = Atoms([x[0] for x in rows], positions=[x[1:] for x in rows], pbc=False)
            a.info.update(charge=task['charge'], spin=task['spin_multiplicity'])
            calc.reset(); a.calc = calc
            result['input_state_check'] = omol.input_batch(calc, a, task['charge'], task['spin_multiplicity'], metal_index=0)
            value = float(a.get_potential_energy()); forces = np.asarray(a.get_forces(), dtype=np.float64)
            if not np.isfinite(value) or forces.shape != (len(rows), 3) or not np.isfinite(forces).all():
                raise InvalidArtifact('invalid native energy/forces')
            fp = out / 'forces_eV_A.npy'; np.save(fp, forces)
            result.update(status='computed', energy_eV=value, forces=record(fp),
                parameter_versions_unchanged=versions == {n: p._version for n, p in model.named_parameters()},
                force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy')
            if not omol.accepted_state(result, task):
                raise InvalidArtifact('native state audit failed')
            if task['task_id'] in m['verification_references']:
                old = read_json(verify(m['verification_references'][task['task_id']]['receipt']))
                delta = (value - old['energy_eV']) * omol.EV_TO_KCAL
                result['execution_verification'] = {'difference_kcal_mol': delta, 'tolerance_kcal_mol': .01,
                    'pass': abs(delta) <= .01, 'reference': m['verification_references'][task['task_id']]['receipt']}
                if abs(delta) > .01: raise InvalidArtifact('warm native replay failed existing .01 endpoint tolerance')
        except Exception as exc:
            result.update(status='failed', reason=str(exc), exception_type=type(exc).__name__)
        result.update(evaluation_seconds=time.monotonic() - tick, model_load_seconds=load_seconds,
            slurm_job_id=os.environ['SLURM_JOB_ID'], allocated_cpus=os.environ['SLURM_CPUS_PER_TASK'],
            allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'), device=torch.cuda.get_device_name(),
            versions={n: importlib.metadata.version(n) for n in ('torch', 'mace-torch', 'ase')})
        write_new(out / 'result.json', result); receipts.append(record(out / 'result.json'))
        print(json.dumps({k: result.get(k) for k in ('task_id', 'status', 'reason')}), flush=True)
        if task['task_id'] in m['verification_references'] and result['status'] != 'computed':
            raise InvalidArtifact('warm execution qualification failed; primary work stopped')
    result = {'status': 'shard_finished', 'shard': shard, 'shards': shards, 'receipts': receipts,
              'model_load_seconds': load_seconds, 'wall_seconds': time.monotonic()-start,
              'peak_host_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'peak_cuda_allocated_bytes': torch.cuda.max_memory_allocated(), 'manifest': manifest_pin,
              'slurm_job_id': os.environ['SLURM_JOB_ID']}
    write_new(Path(manifest).parent / ('execution_' + os.environ['SLURM_JOB_ID'] + '_' + str(shard) + '.json'), result)
    return {k: v for k, v in result.items() if k != 'receipts'}


def inventory(manifest, output):
    import pqq_fast_prepare as source_prep
    from compact_solvation_compare import native_endpoint, representation
    m = read_json(manifest); p = read_json(verify(m['preparation'])); cases = []; failed = []
    for row in p['cases']:
        if row['status'] != 'prepared':
            failed.append({'case_id': row['case_id'], 'reason': row.get('reason', row['status'])}); continue
        try:
            reps = {}
            for rep in ('core', 'context'):
                eps = {}
                for z in ('Ca', 'La'):
                    task_id = row['case_id'] + '__' + rep + '__' + z
                    result_path = (verify(m['reused'][task_id]['receipt']) if task_id in m.get('reused', {})
                                   else Path(manifest).parent / 'results' / task_id / 'result.json')
                    result = read_json(result_path)
                    eps[z] = native_endpoint(result, m['model'])
                prep = row['core']['parent'] if rep == 'core' else row['representations'][rep]['preparation']
                protocol = source_prep.CORE_PROTOCOL if rep == 'core' else source_prep.CONTEXT_PROTOCOL
                reps[rep] = representation(eps, prep, protocol)
            cases.append({**{k: row['source'][k] for k in ('case_id', 'biological_group', 'expected_class',
                'role', 'evidence_stratum', 'label_scope', 'root_case_id', 'source_conditioning_metal',
                'canonical_coordinate_match', 'primary_evaluation_pool')}, 'representations': reps,
                'all_evidence_consumed': True})
        except (KeyError, TypeError, OSError, InvalidArtifact) as exc:
            failed.append({'case_id': row['case_id'], 'reason': str(exc)})
    result = {'schema_version': SCHEMA, 'software': m['software'], 'model': m['model'],
              'cases': cases, 'failures': failed, 'denominator': p['denominator'],
              'mace_manifest': record(manifest), 'source_preparation': m['preparation']}
    write_new(output, result)
    return {'complete_cases': len(cases), 'failed_cases': len(failed), 'inventory': record(output)}


def split_solvent(manifest, output, shards=4):
    """Disjoint case manifests for existing ORCA executors; no changed inputs."""
    from compact_solvation import validate
    validate(manifest); m = read_json(manifest)
    if int(shards) < 1: raise InvalidArtifact('invalid shard count')
    if any(Path(t['output_path']).exists() for t in m['tasks']):
        raise InvalidArtifact('split only before any parent-manifest execution')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False); paths = []; seen = set()
    for i in range(int(shards)):
        cases = m['selection']['cases'][i::int(shards)]; part = copy.deepcopy(m)
        part['selection']['cases'] = cases
        for field in ('all_tasks', 'tasks'):
            part[field] = [t for t in m[field] if t['case_id'] in cases]
        ids = {t['task_id'] for t in part['all_tasks']}
        if seen & ids: raise InvalidArtifact('overlapping split task outputs')
        seen |= ids
        part['reused'] = {k: v for k, v in m['reused'].items() if k in ids}
        part['unexecuted_parent_manifest'] = record(manifest)
        directory = out / ('shard_' + str(i)); directory.mkdir()
        (directory / 'implementation').symlink_to(Path(manifest).resolve().parent / 'implementation', target_is_directory=True)
        for task in part['all_tasks']:
            td = directory / 'tasks' / task['task_id']; td.mkdir(parents=True)
            for field, filename in (('xyz', 'core.xyz'), ('input', 'endpoint.inp')):
                dp = td / filename; shutil.copyfile(verify(task[field]), dp); task[field] = record(dp)
            task['output_path'] = str(td / 'endpoint.out')
        part['tasks'] = [t for t in part['all_tasks'] if t['task_id'] not in part['reused']]
        path = directory / 'manifest.json'; write_new(path, part); validate(path); paths.append(record(path))
    if seen != {t['task_id'] for t in m['all_tasks']}: raise InvalidArtifact('split coverage differs')
    result = {'parent': record(manifest), 'manifests': paths, 'all_tasks': len(seen),
              'new_tasks': len(m['tasks']), 'reused': len(m['reused']), 'disjoint_outputs': True}
    write_new(out / 'inventory.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare')
    for k in ('sources', 'agreement', 'output'): q.add_argument('--' + k, required=True)
    q.add_argument('--workers', type=int, default=16)
    q = sub.add_parser('prepare-mace')
    for k in ('preparation', 'agreement', 'output'): q.add_argument('--' + k, required=True)
    q.add_argument('--reuse-inventory')
    q = sub.add_parser('reconcile')
    for k in ('preparation', 'output'): q.add_argument('--' + k, required=True)
    q = sub.add_parser('run-mace'); q.add_argument('--manifest', required=True)
    q.add_argument('--shard', type=int, default=0); q.add_argument('--shards', type=int, default=1)
    q = sub.add_parser('inventory')
    for k in ('manifest', 'output'): q.add_argument('--' + k, required=True)
    q = sub.add_parser('split-solvent')
    for k in ('manifest', 'output'): q.add_argument('--' + k, required=True)
    q.add_argument('--shards', type=int, default=4)
    args = vars(parser.parse_args()); op = args.pop('op').replace('-', '_')
    print(json.dumps(globals()[op](**args)))


if __name__ == '__main__': main()
