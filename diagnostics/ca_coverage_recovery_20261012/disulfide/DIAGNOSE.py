"""Inspect one failed disulfide preparation and two real same-protein controls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback
from collections import Counter


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(path):
    if isinstance(path, dict):
        assert pin(path['path']) == path
        path = path['path']
    return json.loads(Path(path).read_text())


def write(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def diagnose(plan):
    sys.path.insert(0, plan['scripts'])
    from second_shell_context import parent_state
    from environment_context_chemistry import cysteine_units
    from consistent_context import forced_expansion
    from openmm.app import PDBFile
    import openmm.app.topology
    old = read(plan['prior_manifest'])
    sources = {row['case_id']: row for row in read(old['source_preparation'])['cases']}
    group = next(g for g in old['groups'] if g['selection_id'] == plan['selection_id'])
    results = []
    for cid in plan['source_ids']:
        row = sources[cid]
        state = parent_state(row['core'], old['config']['topology'], require_endpoint_receipts=False)
        graph = state['graph']
        resolved = {(f['kind'], graph.locate(f['source']).key) for f in group['fragments']}
        wanted = {key for kind, key in resolved if kind == 'sidechain'
                  and graph.residues[key].canonical_resname == 'CYS'}
        topology = PDBFile(str(graph.source)).topology
        def key(atom):
            r = atom.residue
            target = graph.locate({'chain': r.chain.id, 'resnum': int(r.id),
                                  'resname': r.name, 'insertion_code': r.insertionCode.strip()})
            return graph.key(target.key, atom.name)
        bonds = [(key(a), key(b)) for a, b in topology.bonds() if a.name == b.name == 'SG']
        relevant = set(wanted)
        for a, b in bonds:
            if a[:2] in wanted or b[:2] in wanted:
                relevant.update((a[:2], b[:2]))
        residues = []
        for rkey in sorted(relevant):
            keys = graph.by_residue[rkey]
            atoms = []
            for name, atomkey in sorted(keys.items()):
                atom = graph.atoms[atomkey]
                record = {'atom': name, 'element': atom.element.name,
                          'xyz_A': [atom.pos.x, atom.pos.y, atom.pos.z]}
                if atom.element.name in ('H', 'D'):
                    record['graph_parent'] = graph.hparents.get(atomkey)
                    distances = sorted((atom.pos.dist(graph.atoms[k].pos), n)
                                       for n, k in keys.items()
                                       if graph.atoms[k].element.name not in ('H', 'D'))
                    record['nearest_heavy_atoms_A'] = [{'atom': n, 'distance_A': d} for d, n in distances[:4]]
                atoms.append(record)
            residues.append({'source': graph.residues[rkey].source_dict(), 'selected': rkey in wanted,
                             'H_parent_counts': dict(Counter(p[2] for p in graph.hparents.values() if p[:2] == rkey)),
                             'atoms': atoms})
        diagnostics = []
        for a, b in bonds:
            if a[:2] in relevant or b[:2] in relevant:
                diagnostics.append({'a': graph.meta[a], 'b': graph.meta[b],
                                    'SG_SG_distance_A': graph.atoms[a].pos.dist(graph.atoms[b].pos)})
        neighbors = []
        for rkey in sorted(wanted):
            a = graph.key(rkey, 'SG')
            for b in sorted(graph.atoms):
                if b == a or b[2] != 'SG': continue
                distance = graph.atoms[a].pos.dist(graph.atoms[b].pos)
                if distance <= 3.0:
                    neighbors.append({'a': graph.meta[a], 'b': graph.meta[b], 'distance_A': distance})
        item = {'case_id': cid, 'source': pin(graph.source), 'core_parent': row['core']['parent'],
                'raw_source': row['source'].get('original_source_structure', row['source']['source_structure']),
                'residues': residues, 'source_topology_SG_bonds': diagnostics,
                'selected_SG_neighbors_within_3A': neighbors}
        try:
            units, meta = cysteine_units(state, resolved)
            item['cysteine_gate'] = {'status': 'passed', 'closures': meta,
                                     'units': {str(k): sorted(v) for k, v in units.items()}}
        except Exception as error:
            item['cysteine_gate'] = {'status': 'failed', 'reason': str(error), 'type': type(error).__name__}
        try:
            paired, audit = forced_expansion(state, group['fragments'])
            item['full_expansion'] = {'status': 'passed', 'atoms': len(paired['La']),
                                      'added_formal_charge': audit['added_formal_charge'],
                                      'disulfide_closures': audit['disulfide_closures']}
        except Exception as error:
            item['full_expansion'] = {'status': 'failed', 'reason': str(error), 'type': type(error).__name__}
        results.append(item)
    return {'status': 'diagnostic_complete', 'sources': results, 'selection_id': plan['selection_id'],
            'source_topology_reader': pin(openmm.app.topology.__file__),
            'new_molecular_calls': 0, 'new_hydrogens': 0, 'source_coordinates_changed': False,
            'interpretation': 'Source-graph/unchanged-gate diagnostic only; no repair or score yet.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise SystemExit('Scientific parsing/preparation requires Slurm')
    plan = read(args.plan); output = Path(plan['output_directory']); started = time.time()
    receipt = {'job': os.environ['SLURM_JOB_ID'], 'plan': pin(args.plan), 'started_epoch': started}
    write(output/'START.json', receipt)
    try:
        for item in plan['pins']: assert pin(item['path']) == item, item['path']
        result = diagnose(plan); write(output/'RESULT.json', result)
        receipt.update(status='complete', result=pin(output/'RESULT.json'))
    except Exception as error:
        receipt.update(status='failed', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        receipt.update(finished_epoch=time.time(), wall_seconds=time.time()-started, new_molecular_calls=0)
        write(output/'RECEIPT.json', receipt)
