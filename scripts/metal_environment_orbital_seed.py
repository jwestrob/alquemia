"""Strict initial-orbital reuse for the research frozen-f embedded path.

This returns initialization provenance, never a target energy or convergence.
Target manifest remains responsible for Hamiltonian/geometry admission.
"""
import argparse
import json
import shutil
from datetime import datetime
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz, energy

PROTOCOL = 'nikasha_same_metal_frozen_f_orbital_seed_v2'
GEOMETRIC_KEYS = {'xyz_A', 'retained_xyz_A', 'omitted_xyz_A',
                  'jacobian_retained', 'jacobian_omitted'}


def identity(value):
    if isinstance(value, dict):
        return {k: identity(v) for k, v in value.items() if k not in GEOMETRIC_KEYS}
    if isinstance(value, list):
        return [identity(v) for v in value]
    return value



def boundary_identity_and_diagnostics(pin):
    """Separate only known geometric diagnostics from boundary chemistry."""
    boundary = read_json(verify(pin))
    diagnostics = {}
    for key in ('coordinate_preparation', 'original_boundary_mapping'):
        if key in boundary:
            diagnostics[key] = boundary.pop(key)
    if 'original_boundary_mapping' in diagnostics:
        verify(diagnostics['original_boundary_mapping'])
    diagnostics['ledger'] = []
    for entry in boundary['ledger']:
        geometric = {}
        for key in ('dipole_before_eA', 'retained_MM_dipole_after_eA',
                    'original_dipole_diagnostic'):
            if key in entry:
                geometric[key] = entry.pop(key)
        diagnostics['ledger'].append(geometric)
    return identity(boundary), diagnostics


def task(manifest, task_id):
    matches = [t for t in manifest['tasks'] if t['task_id'] == task_id]
    if len(matches) != 1:
        raise InvalidArtifact('source/target task must be unique')
    return matches[0]


def require(condition, message):
    if not condition:
        raise InvalidArtifact(message)


def validate_seed(source_manifest, collection, source_task_id, target_manifest, target_task_id):
    smp, tmp, cp = map(Path, (source_manifest, target_manifest, collection))
    sm, tm, col = map(read_json, (smp, tmp, cp))
    s, t = task(sm, source_task_id), task(tm, target_task_id)
    require(col['manifest'] == record(smp), 'collection does not pin source manifest')
    row = col['rows'][source_task_id]
    require(row['status'] == 'complete', 'source task not accepted')
    receipt = read_json(verify(row['receipt']))
    require(receipt['task_id'] == source_task_id and receipt['manifest'] == record(smp), 'receipt task/manifest mismatch')
    require(receipt.get('returncode') == 0 and receipt.get('normal_termination') is True
            and receipt.get('scf_converged') is True and not receipt.get('launch_error')
            and receipt.get('expected_engrad_complete') is True, 'source execution not successfully terminal')
    require(receipt['orca_executable'] == sm['orca'], 'source executable receipt mismatch')
    for pin in receipt['artifacts'].values():
        verify(pin)
    require(row['output'] == receipt['artifacts']['output'] == record(s['output_path']), 'output provenance mismatch')
    require(energy(verify(row['output'])) == row['energy_hartree'], 'accepted energy/output mismatch')
    require(row['electronic_state'] == s['electronic_state'], 'collected electronic state mismatch')
    for key in ('metal', 'charge', 'multiplicity', 'physical_multiplicity', 'electronic_state'):
        require(s[key] == t[key], f'incompatible {key}')
    require(s['electronic_state']['representation'] == 'spin_free_frozen_f_valence_model', 'only frozen-f valence path supported')
    require(sm['method'] == tm['method'] and 'NoAutostart' in tm['method'], 'method mismatch or implicit autostart')
    require(sm['orca'] == tm['orca'], 'executable mismatch')
    verify(sm['orca'])
    for kind in ('basis', 'aux'):
        a, b = sm['assets'][s['metal']][kind], tm['assets'][t['metal']][kind]
        verify(a); verify(b)
        require(a['sha256'] == b['sha256'], f'{kind}/ECP hash mismatch')
    for endpoint in (s, t):
        for key in ('input', 'xyz', 'pointcharges', 'core_mapping', 'boundary_mapping'):
            verify(endpoint[key])
    maps = [read_json(verify(x['core_mapping'])) for x in (s, t)]
    require(identity(maps[0]) == identity(maps[1]), 'ordered source/cap identity mismatch')
    sb, sd = boundary_identity_and_diagnostics(s['boundary_mapping'])
    tb, td = boundary_identity_and_diagnostics(t['boundary_mapping'])
    require(sb == tb, 'boundary identity mismatch')
    rs, rt = xyz(verify(s['xyz'])), xyz(verify(t['xyz']))
    es, et = [x[0] for x in rs], [x[0] for x in rt]
    qs, qt = [x[1:] for x in rs], [x[1:] for x in rt]
    require(es == et and len(maps[0]) == len(es), 'ordered element/map mismatch')
    require([x['qm_index'] for x in maps[0]] == list(range(len(es))), 'mapping indices not ordered')
    # GBW is not in historical execution receipts: pin it only after successful
    # immutable output/receipt checks, and reject writes after execution finished.
    gbw = Path(s['output_path']).with_name('endpoint.runtime.gbw')
    finished = datetime.fromisoformat(receipt['finished_at_utc']).timestamp()
    require(gbw.is_file() and gbw.stat().st_size > 0, 'completed source GBW unavailable')
    require(gbw.stat().st_mtime <= finished, 'GBW modified after source execution')
    pin = record(gbw)
    return dict(protocol_id=PROTOCOL, purpose='initial_guess_only', target_energy=None,
                target_convergence=None, source_manifest=record(smp), source_collection=record(cp),
                source_task_id=source_task_id, source_output=row['output'], source_receipt=row['receipt'],
                source_gbw=pin, target_manifest=record(tmp), target_task_id=target_task_id,
                electronic_state=s['electronic_state'], method=sm['method'], orca=sm['orca'],
                assets=sm['assets'][s['metal']], coordinate_changes=dict(source=s['xyz'], target=t['xyz'],
                max_displacement_A=float(np.linalg.norm(np.array(qt)-qs, axis=1).max())),
                boundary_diagnostic_changes=dict(source=s['boundary_mapping'], target=t['boundary_mapping'],
                source_diagnostics=sd, target_diagnostics=td, changed=sd != td),
                field_changes=dict(source=s['pointcharges'], target=t['pointcharges'],
                changed=s['pointcharges']['sha256'] != t['pointcharges']['sha256']),
                scf_initial_guess='Guess MORead\n MOInp "initial.gbw"',
                required_keyword='NoAutostart', automatic_fallback=False)


def stage_seed(source_manifest, collection, source_task_id, target_manifest, target_task_id, directory):
    result = validate_seed(source_manifest, collection, source_task_id, target_manifest, target_task_id)
    dest = Path(directory)
    dest.mkdir(parents=True, exist_ok=False)
    source = verify(result['source_gbw'])
    with source.open('rb') as src, (dest/'initial.gbw').open('xb') as out:
        shutil.copyfileobj(src, out)
    verify(result['source_gbw'])
    result['staged_gbw'] = record(dest/'initial.gbw')
    require(result['staged_gbw']['sha256'] == result['source_gbw']['sha256'], 'GBW copy mismatch')
    for key in ('source_manifest', 'source_collection', 'source_output', 'source_receipt', 'target_manifest'):
        verify(result[key])
    write_new(dest/'SEED.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source-manifest', 'collection', 'source-task-id', 'target-manifest', 'target-task-id'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--directory', help='new directory for initial.gbw and SEED.json; omit to validate only')
    a = vars(p.parse_args()); directory = a.pop('directory')
    print(json.dumps(stage_seed(**a, directory=directory) if directory else validate_seed(**a), indent=2))


if __name__ == '__main__':
    main()
