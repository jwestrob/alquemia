"""Initial guesses from interrupted, unconverged frozen-f jobs; never energies.

Unlike accepted-source seeding, this requires an exactly unchanged target.
The record admits a trial restart, not binary integrity or electronic stability.
"""
import argparse
import csv
import json
import shutil
from pathlib import Path
from affordable_common import read_json, record, verify, write_new
from metal_environment_orbital_seed import require, task

PROTOCOL = 'nikasha_interrupted_frozen_f_initial_guess_v1'


def validate_interrupted(source_manifest, target_manifest, task_id, evidence, inventory):
    sm, tm = read_json(source_manifest), read_json(target_manifest)
    s, t = task(sm, task_id), task(tm, task_id)
    ev, inv = read_json(evidence), read_json(inventory)
    collection_path = Path(source_manifest).parent/f"COLLECTION_{ev['job']}.json"
    collection = read_json(collection_path)
    require(collection['manifest'] == record(source_manifest), 'source manifest differs from interrupted collection')
    require(collection['rows'][task_id]['status'] != 'complete', 'source already accepted')
    accounting = list(csv.DictReader(verify(ev['accounting']).read_text().splitlines(), delimiter='|'))
    terminal = [r for r in accounting if r['JobID'] == ev['job']]
    require(len(terminal) == 1 and terminal[0]['State'].startswith('CANCELLED'), 'interrupted cancellation evidence required')
    require(inv['job'] == ev['job'], 'checkpoint job mismatch')
    require(sm['source_id'] == tm['source_id'] and sm['method'] == tm['method'] and sm['orca'] == tm['orca'], 'source/method/executable mismatch')
    require(sm.get('solver') == tm.get('solver') == 'TRAH', 'unchanged TRAH target required')
    verify(sm['orca']); verify(sm['inputs'])
    require(sm['inputs'] == tm['inputs'], 'source preparation changed')
    for k in ('metal', 'configuration', 'charge', 'multiplicity', 'physical_multiplicity', 'electronic_state'):
        require(s[k] == t[k], 'state mismatch: '+k)
    require(s['electronic_state']['representation'] == 'spin_free_frozen_f_valence_model', 'unsupported electronic representation')
    for k in ('xyz', 'pointcharges', 'core_mapping', 'boundary_mapping'):
        verify(s[k]); verify(t[k])
        require(s[k]['sha256'] == t[k]['sha256'], 'changed target '+k)
    for k in ('basis', 'aux'):
        a, b = sm['assets'][s['metal']][k], tm['assets'][t['metal']][k]
        verify(a); verify(b)
        require(a['sha256'] == b['sha256'], 'changed basis/ECP')
    output = verify(ev['rows'][s['metal']]['output'])
    require(output.resolve() == Path(s['output_path']).resolve(), 'source output mismatch')
    text = output.read_text()
    require('SCF CONVERGED AFTER' not in text and 'ORCA TERMINATED NORMALLY' not in text, 'use accepted-source path for completed sources')
    checkpoint = inv['files'][s['metal']]
    gbw = verify(checkpoint)
    require(gbw.resolve() == output.with_name('endpoint.runtime.gbw').resolve(), 'checkpoint source path mismatch')
    require(gbw.stat().st_size == checkpoint['bytes'] > 0, 'checkpoint size mismatch')
    return dict(protocol_id=PROTOCOL, purpose='unconverged_initial_guess_only',
                source_manifest=record(source_manifest), target_manifest=record(target_manifest),
                source_task_id=task_id, target_task_id=task_id,
                interrupted_collection=record(collection_path),
                evidence=record(evidence), inventory=record(inventory), source_output=record(output),
                source_gbw=record(gbw), target_energy=None, target_convergence=None,
                binary_integrity_qualified=False, automatic_fallback=False)


def stage(source_manifest, target_manifest, task_id, evidence, inventory, directory):
    result = validate_interrupted(source_manifest, target_manifest, task_id, evidence, inventory)
    dest = Path(directory).resolve(); dest.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(verify(result['source_gbw']), dest/'initial.gbw')
    result['staged_gbw'] = record(dest/'initial.gbw')
    require(result['staged_gbw']['sha256'] == result['source_gbw']['sha256'], 'checkpoint copy mismatch')
    write_new(dest/'SEED.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source-manifest', 'target-manifest', 'task-id', 'evidence', 'inventory', 'directory'):
        p.add_argument('--'+name, required=True)
    print(json.dumps(stage(**vars(p.parse_args())), indent=2))


if __name__ == '__main__':
    main()
