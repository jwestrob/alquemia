"""Versioned radial repair of existing unmapped CYS carbon hydrogens only."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PROTOCOL = 'Nikasha_disulfide_unmapped_source_CH_radial_repair_v1'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(path):
    if isinstance(path, dict):
        assert pin(path['path']) == path, path['path']
        path = path['path']
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def recover_source(row, diagnostic, config, out):
    import numpy as np
    from second_shell_context import parent_state
    from affordable_peptide import SourceGraph
    state = parent_state(row['core'], config['topology'], require_endpoint_receipts=False)
    graph = state['graph']; before = {k: np.array(list(a.pos)) for k, a in graph.atoms.items()}
    xml = ET.parse(config['topology']['path']).getroot()
    template = next(r for r in xml.findall('./Residues/Residue') if r.attrib['name'] == 'CYX')
    template_atoms = {a.attrib['name']: a.attrib for a in template.findall('Atom')}
    atom_types = {a.attrib['name']: a.attrib for a in xml.findall('./AtomTypes/Type')}
    audit = {'protocol_id': PROTOCOL, 'original_source': pin(graph.source), 'changes': [],
             'topology': config['topology'], 'new_protons': 0, 'new_waters': 0,
             'heavy_coordinates_changed': False, 'rule': 'unmapped existing CB-bound H in source disulfide only; same ray, ff19SB equilibrium length'}
    relevant = set()
    for bond in diagnostic['source_topology_SG_bonds']:
        assert 1.8 <= bond['SG_SG_distance_A'] <= 2.3
        for key in ('a', 'b'):
            r = graph.locate(bond[key]); relevant.add(r.key)
    changes = {}
    for rkey in sorted(relevant):
        residue = graph.residues[rkey]
        assert residue.canonical_resname == 'CYS'
        names = graph.by_residue[rkey]
        assert 'HG' not in names, 'cannot change thiol protonation'
        for name, key in sorted(names.items()):
            atom = graph.atoms[key]
            if atom.element.name != 'H' or key in graph.hparents: continue
            bonds = [b.attrib for b in template.findall('Bond') if name in (b.attrib['atomName1'], b.attrib['atomName2'])]
            partners = [b['atomName2'] if b['atomName1'] == name else b['atomName1'] for b in bonds]
            assert partners == ['CB'] and name in ('HB2', 'HB3'), 'unsupported unmapped H chemistry'
            parent_key = names['CB']; parent = graph.atoms[parent_key]
            heavy = [k for k in names.values() if graph.atoms[k].element.name not in ('H', 'D')]
            assert min(heavy, key=lambda k: atom.pos.dist(graph.atoms[k].pos)) == parent_key
            vector = before[key] - before[parent_key]; length = float(np.linalg.norm(vector))
            old_cutoff = parent.element.covalent_r + atom.element.covalent_r + .30
            assert old_cutoff < length <= 1.5, 'not a bounded overlong source C-H bond'
            classes = sorted(atom_types[template_atoms[n]['type']]['class'] for n in ('CB', name))
            parameters = [b.attrib for b in xml.findall('./HarmonicBondForce/Bond')
                          if sorted((b.attrib.get('class1'), b.attrib.get('class2'))) == classes]
            assert len(parameters) == 1, 'ambiguous ff19SB parent-H bond parameter'
            target = float(parameters[0]['length']) * 10.0
            assert 1.0 <= target <= 1.15, 'unexpected source C-H equilibrium length'
            xyz = before[parent_key] + vector * (target / length)
            changes[(residue.chain, residue.resnum, residue.insertion_code, name)] = xyz
            audit['changes'].append({'source': graph.meta[key], 'parent': graph.meta[parent_key],
                'before_A': before[key].tolist(), 'unrounded_target_A': xyz.tolist(), 'before_length_A': length,
                'original_parent_cutoff_A': old_cutoff, 'target_length_A': target,
                'forcefield_template': 'CYX', 'forcefield_classes': classes, 'forcefield_bond': parameters[0],
                'unrounded_displacement_A': float(np.linalg.norm(xyz-before[key]))})
    out.mkdir(parents=True, exist_ok=False)
    original_lines = Path(graph.source).read_text().splitlines(keepends=True); lines = []; seen = set()
    for line in original_lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[12:16].strip())
            if key in changes:
                assert key not in seen; seen.add(key)
                x, y, z = changes[key]; line = line[:30] + f'{x:8.3f}{y:8.3f}{z:8.3f}' + line[54:]
        lines.append(line)
    assert seen == set(changes)
    source = out/'protonated_repaired.pdb'; source.write_text(''.join(lines))
    new_graph = SourceGraph(source, config['topology']['path'])
    assert set(new_graph.atoms) == set(graph.atoms)
    changed_keys = []
    for key, atom in new_graph.atoms.items():
        xyz = np.array(list(atom.pos))
        if not np.array_equal(xyz, before[key]): changed_keys.append(key)
    intended = [tuple((c['source']['chain_index'], c['source']['residue_index'], c['source']['atom'])) for c in audit['changes']]
    assert sorted(changed_keys) == sorted(intended), 'unrequested coordinate changed'
    for change, key in zip(audit['changes'], intended):
        parent = (key[0], key[1], 'CB'); xyz = np.array(list(new_graph.atoms[key].pos))
        assert new_graph.hparents[key] == parent, 'unchanged graph parent gate still fails'
        length = float(np.linalg.norm(xyz-before[parent]))
        assert abs(length-change['target_length_A']) <= .001, 'PDB roundoff exceeds frozen tolerance'
        displacement = float(np.linalg.norm(xyz-before[key])); assert displacement <= .5
        change.update(after_A=xyz.tolist(), after_length_A=length, actual_displacement_A=displacement)
    audit.update(repaired_source=pin(source), all_other_coordinates_exact=True,
                 atom_inventory_exact=True, existing_H_parent_gate_passed=True)
    write(out/'REPAIR.json', audit)
    if not changes:
        assert pin(source)['sha256'] == pin(graph.source)['sha256']
        return copy.deepcopy(row), audit
    # Clone only provenance links. Every frozen electronic-core coordinate remains archived.
    parent = read(row['core']['parent']); pm = read(parent['protonation_manifest'])
    pm.update(output=pin(source), original_protonation_manifest=parent['protonation_manifest'],
              recovery_protocol=PROTOCOL, coordinate_repair=pin(out/'REPAIR.json'),
              proton_inventory_unchanged=True, heavy_coordinates_unchanged=True)
    write(out/'protonated_repair_manifest.json', pm)
    parent.update(protonation_manifest=pin(out/'protonated_repair_manifest.json'),
                  original_parent=row['core']['parent'], context_source_recovery=pin(out/'REPAIR.json'))
    write(out/'core_parent_context_repaired.json', parent)
    new = copy.deepcopy(row); new['core']['parent'] = pin(out/'core_parent_context_repaired.json')
    new['context_source_recovery'] = pin(out/'REPAIR.json')
    write(out/'SOURCE_PREPARATION.json', new)
    return new, audit


def execute(plan):
    sys.path.insert(0, plan['scripts'])
    import motion_envelope as envelope
    spec = importlib.util.spec_from_file_location('prior_ca_preparation', plan['prior_prepare_driver']['path'])
    prior = importlib.util.module_from_spec(spec); spec.loader.exec_module(prior)
    old = read(plan['prior_manifest']); sources = {r['case_id']: r for r in read(old['source_preparation'])['cases']}
    group = next(g for g in old['groups'] if g['selection_id'] == plan['selection_id'])
    diagnostic = {r['case_id']: r for r in read(plan['diagnosis'])['sources']}
    old_rows = {r['case_id']: r for r in read(plan['old_Ca_preparation'])['rows'] if r.get('selection_id') == group['selection_id']}
    anchor = next(r for r in read(old['old_preparation'])['cases']
                  if r['case_id'] == group['members'][0] and r['selection_id'] == group['selection_id'])
    assert anchor['status'] == 'prepared' and anchor['mapping_eligible']
    out = Path(plan['output_directory']); rows = []; checks = []; repairs = []
    for cid in plan['source_ids']:
        row, repair = recover_source(sources[cid], diagnostic[cid], old['config'], out/'source_repairs'/cid)
        write(out/'source_repairs'/cid/'ROW_FOR_PREPARATION.json', row)
        recovered = envelope.prepare_task(row, group, old['config'], out/'prepared')
        recovered = prior.finalized(recovered, group, anchor)
        recovered.update(recovery_protocol=PROTOCOL, repair=pin(out/'source_repairs'/cid/'REPAIR.json'),
                         original_failure=diagnostic[cid]['full_expansion'])
        if cid == plan['target_source']:
            assert len(repair['changes']) == 1, 'diagnosed one-H target no longer exact'
            assert recovered['status'] == 'prepared', recovered.get('reason')
        else:
            assert not repair['changes'], 'previously valid control would change'
            archived = anchor if cid == anchor['case_id'] else old_rows[cid]
            checks.append(prior.compare_probe(recovered, archived))
        write(out/'prepared'/'cases'/cid/'ADMISSION.json', recovered)
        rows.append(recovered); repairs.append(repair)
    result = {'protocol_id': PROTOCOL, 'status': 'prepared', 'denominator': 3, 'target_source': plan['target_source'],
              'rows': rows, 'config': old['config'], 'controls': checks, 'repairs': repairs,
              'canonical_reference': plan['canonical_reference'], 'frozen_state_anchor': anchor['state_key'],
              'new_molecular_energy_calls': 0, 'new_optimizations': 0, 'production_changed': False,
              'interpretation': 'Source-coordinate preparation recovery only; no new score or classification.'}
    write(out/'PREPARATION.json', result)
    return {'status': 'complete', 'preparation': pin(out/'PREPARATION.json'), 'prepared_sources': len(rows), 'passed_controls': len(checks)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--plan', required=True)
    args = ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise SystemExit('All scientific preparation requires Slurm')
    plan = read(args.plan); out = Path(plan['output_directory']); started = time.time()
    receipt = {'job_id': os.environ['SLURM_JOB_ID'], 'plan': pin(args.plan), 'started_epoch': started}
    write(out/'START.json', receipt)
    try:
        for p in plan['pins']: assert pin(p['path']) == p, p['path']
        receipt.update(execute(plan))
    except BaseException as error:
        receipt.update(status='failed', reason=str(error), traceback=traceback.format_exc()); raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(finished_epoch=time.time(), wall_seconds=time.time()-started,
                       process_CPU_seconds=usage.ru_utime+usage.ru_stime, peak_RSS_KiB=usage.ru_maxrss,
                       molecular_energy_calls=0)
        write(out/'RECEIPT.json', receipt)
