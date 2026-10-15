#!/usr/bin/env python3
"""Metadata-only inventory of already consumed MMOL source0/source1 artifacts.

No scientific imports, coordinate calculations, scoring or molecular calls.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
REVIEW = ROOT / 'workspaces/pqq_inconclusive_review_20261014/REVIEW_INPUT.json'
SCORES = ROOT / 'workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/RESULT.json'
NATIVE = ROOT / 'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/transfer_v1/tasks'
verified = {}


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def check(p):
    path = Path(p['path'])
    if p['path'] not in verified:
        observed = pin(path)
        if observed != p:
            raise ValueError('Changed artifact: ' + p['path'])
        verified[p['path']] = p['sha256']
    elif verified[p['path']] != p['sha256']:
        raise ValueError('Inconsistent source pins')
    return path


def read(p):
    return json.loads((check(p) if isinstance(p, dict) else Path(p)).read_text())


def pin_tree(obj):
    if isinstance(obj, dict):
        if set(obj) == {'path', 'sha256'}:
            check(obj)
        else:
            for v in obj.values():
                pin_tree(v)
    elif isinstance(obj, list):
        for v in obj:
            pin_tree(v)


def main():
    review, scores = read(REVIEW), read(SCORES)
    rows = []
    for sample, role in [(0, 'target'), (1, 'control_inventory_only')]:
        sid = f'mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-{sample}'
        entry, = [c for c in review['cases'] if c['record']['source_case_id'] == sid]
        r = entry['record']; source = r['native_source_record']
        score, = [x for x in scores['rows'] if x['record_id'] == r['record_id']]
        admission = read(source['prepared_admission']); state = read(r['state_signature'])
        preparation = read(r['source_preparation'])
        proposals = {z: read(source['proposals'][z]) for z in ('Ca', 'La')}
        proposal_manifest = proposals['Ca']['manifest']
        assert proposals['La']['manifest'] == proposal_manifest
        pm = read(proposal_manifest)
        row = {k: r[k] for k in ('case_id', 'source_case_id', 'selection_id', 'record_id',
               'biological_group', 'expected_class', 'source_conditioning_metal', 'state_signature',
               'maps', 'source_preparation', 'baseline_manifest', 'baseline_collection')}
        row.update(sample=sample, role=role, prepared_admission=source['prepared_admission'],
                   native_proposals=source['proposals'], proposal_manifest=proposal_manifest,
                   source=preparation['source'], native_source_record=source,
                   state={'atom_count': admission['atom_count'], 'cap_count': admission['cap_count'],
                     **{k: state[k] for k in ('endpoint_charges', 'multiplicities', 'water_inventory', 'core_microstates')}},
                   admission_status=admission['status'],
                   score_result=pin(SCORES), stored_score={k: score[k] for k in
                       ('geometry_choices', 'works', 'paired_candidates', 'accommodation_R_shift',
                        'accommodation_MACE_shift', 'accommodation_solvent_shift', 'selected_call')},
                   endpoints={z: r['matrix'][z]['origin']['xyz'] for z in ('Ca', 'La')},
                   proposal_metadata={}, matrix={})
        for z in ('Ca', 'La'):
            prop = proposals[z]
            task, = [t for t in pm['tasks'] if t['metal'] == z]
            row['proposal_metadata'][z] = {'task': task, 'full_q': prop['proposal']['full_q'],
                  'final_geometry': prop['final_geometry'], 'residual': prop['residual'],
                  'status': prop['status'], 'settings': pm['settings']}
            row['matrix'][z] = {}
            for candidate in ('origin', 'adaptive_Ca', 'adaptive_La'):
                cell = r['matrix'][z][candidate]
                nr = pin(NATIVE / cell['native_task_id'] / 'RESULT.json'); native = read(nr)
                mace = read(cell['MACE']); inp = native['input']
                assert native['status'] == 'available'
                assert inp['xyz'] == cell['xyz']
                assert inp['state_signature'] == row['state_signature']
                assert inp['source_preparation'] == row['source_preparation']
                assert mace['status'] == 'complete'
                row['matrix'][z][candidate] = {
                    **{k: cell[k] for k in ('xyz', 'MACE', 'MACE_eV', 'native_task_id')},
                    **{k: inp[k] for k in ('charge', 'multiplicity', 'atom_count', 'valence_electrons')},
                    'MACE_request': mace['request'], 'MACE_forces': mace['forces'],
                    'native_result': nr, 'native_input': inp,
                    'frozen_solvent_hartree': native['G_qvac_hartree'],
                    'native_water_diagonal_checks': native['checks'],
                    'native_vacuum_precision': native['vacuum_precision'],
                    'native_vacuum_components': native['vacuum']['components_Eh'],
                    'exact_native_DFT': None, 'exact_native_DFT_gradient': None,
                    'DFT_reuse_status': 'not located in the bounded declared archive search'}
        # Check directly reused metadata/raw receipts, excluding model binary duplicates.
        for key in ('prepared_admission', 'state_signature', 'maps', 'source_preparation',
                    'source', 'native_proposals', 'proposal_manifest', 'endpoints', 'matrix'):
            pin_tree(row[key])
        rows.append(row)

    paths = [
      'workspaces/accommodation_fold_DFT_20260920/run_v1/manifest.json',
      'workspaces/accommodation_torsion_20260920/prepared_v3/dft/manifest.json',
      'workspaces/exaa_donor_response_20261014/reference_v1/manifest.json',
      'workspaces/pede_solvent_response_20261013/reference_v1/manifest.json']
    hashes = {c['xyz']['sha256'] for r in rows for m in r['matrix'].values() for c in m.values()}
    search = [{'manifest': pin(ROOT / p), 'current_XYZ_hash_mentions':
                 sorted(h for h in hashes if h in (ROOT / p).read_text())} for p in paths]
    old_collection = ROOT / 'workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json'
    old_rows = [r for r in read(old_collection)['rows'] if r['case_id'] in {r['source_case_id'] for r in rows}]
    old = []
    for r in old_rows:
        for z, e in r['endpoints'].items():
            old.append({'source_case_id': r['case_id'], 'metal': z, 'endpoint': e,
                        'atom_count_from_XYZ_header': int(check(e['xyz']).read_text().splitlines()[0]),
                        'input_text': check(e['input']).read_text(),
                        'analytic_gradient_files': [str(p) for p in Path(e['output']['path']).parent.glob('*.engrad')],
                        'reuse_allowed': False,
                        'reason': '80-atom origin core/CPCM differs from current193-atom context/vacuum; no saved proposal path'})
    inventory = {'status': 'metadata_inventory_complete', 'scientific_calculations_executed': 0,
      'review_input': pin(REVIEW), 'scores': pin(SCORES), 'rows': rows,
      'available_current_cells': 12, 'current_DFT_exact_reuse_found': 0,
      'current_DFT_exact_gradient_reuse_found': 0,
      'control_rule': 'first available same-protein/same-Ca-conditioning peer by sample index after target0; sample1 retained despite known extra-Asp geometry',
      'DFT_archive_search': search, 'old_DFT_collection': pin(old_collection), 'incompatible_old_DFT': old,
      'scope_limit': 'Known completed donor, full-fold and recent PedE/ExaA reference manifests; no assertion of an exhaustive filesystem absence proof',
      'proposal': {'target_only_cell_ids': [f'sample0_{z}_{c}_vacuum' for z in ('Ca', 'La') for c in ('origin', 'adaptive_Ca', 'adaptive_La')],
        'first_scout': ['sample0_Ca_origin_vacuum', 'sample0_La_origin_vacuum'],
        'conditional_follow_on': 'four saved proposal cells after actual origin execution/admission; no optimizer or new coordinates',
        'control_status': 'inventory_only; no source1 execution implied'},
      'verified_artifact_pins': [{'path': p, 'sha256': h} for p, h in sorted(verified.items())],
      'metadata_builder': pin(__file__)}
    with (OUT / 'INVENTORY.json').open('x') as f:
        json.dump(inventory, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({'output': pin(OUT / 'INVENTORY.json'), 'rows': len(rows),
                      'cells': 12, 'verified_pins': len(verified), 'archive_exact_hash_matches': search}))


if __name__ == '__main__':
    main()
