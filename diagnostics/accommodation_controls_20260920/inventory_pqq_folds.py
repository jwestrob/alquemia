#!/usr/bin/env python3
"""Inventory existing AF3 samples for the frozen 25 reference proteins.

No confidence or energy selection, folding, preparation, or scoring. The two
proposed extra La-conditioned samples are the first two numeric sample IDs
other than the already-selected canonical geometry. All samples remain listed.
"""
from pathlib import Path
import hashlib
import json
import csv
from collections import Counter
import gemmi

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = Path('/groups/banfield/projects/multienv/corkscrew/supplementary_structures')
OUT = ROOT / 'workspaces/accommodation_controls_20260920'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def pin(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''):
            h.update(chunk)
    return {'path': str(p), 'sha256': h.hexdigest(), 'size_bytes': p.stat().st_size}


def fingerprint(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    assert len(structure) == 1
    atoms = [(c.name, r.seqid.num, r.seqid.icode.strip(), r.name, a.name,
              a.element.name, list(a.pos)) for c in structure[0] for r in c for a in r]
    sequence = ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code
                       for c in structure[0] for r in c if r.het_flag == 'A')
    assert sequence and 'X' not in sequence
    ligand_signature = [(c.name, r.seqid.num, r.name,
                         sorted((a.name, a.element.name) for a in r))
                        for c in structure[0] for r in c if r.het_flag != 'A']
    protein_b = [a.b_iso for c in structure[0] for r in c if r.het_flag == 'A' for a in r]
    protein_signature = [(a[0], a[1], a[2], a[3], a[4], a[5]) for a in atoms if a[0] == 'A']
    return {'coordinates_sha256': sha(json.dumps(atoms, separators=(',', ':')).encode()),
            'atom_identity_sha256': sha(json.dumps([a[:-1] for a in atoms], separators=(',', ':')).encode()),
            'protein_atom_identity_sha256': sha(json.dumps(protein_signature, separators=(',', ':')).encode()),
            'protein_sequence': sequence, 'sequence_sha256': sha(sequence.encode()),
            'atom_count': len(atoms), 'ligand_signature': ligand_signature,
            'software_name': [gemmi.cif.as_string(x) for x in block.find_values('_software.name')],
            'software_version': [gemmi.cif.as_string(x) for x in block.find_values('_software.version')],
            'model_group_name': [gemmi.cif.as_string(x) for x in block.find_values('_ma_model_list.model_group_name')],
            'mean_protein_atom_confidence_b_field': sum(protein_b)/len(protein_b)}


def numeric_sample(path):
    a, b = path.name.split('_')
    return int(a.split('-')[1]), int(b.split('-')[1])


def main():
    calibration_path = ROOT / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/prepared/preparation.json'
    calibration = json.loads(calibration_path.read_text())
    rows = []
    errors = []
    for ref in calibration['targets']:
        source = Path(ref['source_cif']['path'])
        assert pin(source)['sha256'] == ref['source_cif']['sha256']
        original = fingerprint(source)
        base_id = ref['panel_id'].removesuffix('-pqq-la_model')
        modes = {}
        for metal in ('La', 'Ca'):
            parent = ARCHIVE / f"{ref['class']}-binding-controls" / f'{metal}-alone'
            suffix = base_id + '-pqq-' + metal.lower()
            dirs = [p for p in parent.iterdir() if p.is_dir() and p.name.lower() in (suffix, 'xoxf-'+suffix)]
            if len(dirs) != 1:
                errors.append({'panel_id': ref['panel_id'], 'metal': metal, 'error': 'ambiguous or missing archive directory'})
                continue
            folder = dirs[0]
            data_paths = list(folder.glob('*_data.json'))
            assert len(data_paths) == 1
            data = json.loads(data_paths[0].read_text())
            proteins = [x['protein'] for x in data['sequences'] if 'protein' in x]
            ligands = [x['ligand'] for x in data['sequences'] if 'ligand' in x]
            assert len(proteins) == 1 and len(ligands) == 2
            assert proteins[0]['sequence'] == original['protein_sequence']
            expected_metal = '[La+3]' if metal == 'La' else '[Ca+2]'
            assert any(x.get('smiles') == expected_metal for x in ligands)
            pqq = [x for x in ligands if x.get('smiles') != expected_metal]
            assert len(pqq) == 1
            protein_input_hash = sha(json.dumps(proteins[0], sort_keys=True, separators=(',', ':')).encode())
            mode = {'conditioned_metal': metal, 'archive_dir': str(folder), 'input_data': pin(data_paths[0]),
                    'input_dialect': data['dialect'], 'input_format_version': data['version'],
                    'model_seeds': data['modelSeeds'], 'ligands': ligands,
                    'protein_input_including_msa_templates_sha256': protein_input_hash,
                    'ranking_file': pin(folder/'ranking_scores.csv'), 'samples': []}
            for sample in sorted(folder.glob('seed-*_sample-*'), key=numeric_sample):
                if not sample.is_dir():
                    continue
                cif = sample/'model.cif'
                fp = fingerprint(cif)
                assert fp['protein_sequence'] == original['protein_sequence']
                assert fp['protein_atom_identity_sha256'] == original['protein_atom_identity_sha256']
                if metal == 'La':
                    assert fp['ligand_signature'] == original['ligand_signature']
                summary = json.loads((sample/'summary_confidences.json').read_text())
                seed, sample_id = numeric_sample(sample)
                is_original = fp['coordinates_sha256'] == original['coordinates_sha256']
                mode['samples'].append({'seed': seed, 'sample': sample_id, 'sample_id': sample.name,
                    'cif': pin(cif), 'summary_confidences': pin(sample/'summary_confidences.json'),
                    'confidences': pin(sample/'confidences.json'), 'summary': summary,
                    'coordinate_fingerprint_sha256': fp['coordinates_sha256'],
                    'atom_identity_sha256': fp['atom_identity_sha256'],
                    'atom_count': fp['atom_count'], 'ligand_signature': fp['ligand_signature'],
                    'mean_protein_atom_confidence_b_field': fp['mean_protein_atom_confidence_b_field'],
                    'software_name': fp['software_name'], 'software_version': fp['software_version'],
                    'model_group_name': fp['model_group_name'],
                    'matches_canonical_coordinate_fingerprint': is_original,
                    'matches_canonical_file_sha256': pin(cif)['sha256'] == ref['source_cif']['sha256'],
                    'role': 'canonical_selected_geometry' if is_original else
                        ('unselected_same_La_conditioning_sample' if metal == 'La' else 'separate_Ca_conditioning_sample'),
                    'energies_inspected_here': False,
                    'prior_score_exposure': 'canonical_scores_consumed' if is_original else 'not_selected_canonical_geometry; other_historical_energy_exposure_not_globally_audited'})
            modes[metal] = mode
        if set(modes) == {'Ca', 'La'}:
            # Freeze conditioning differences rather than treating Ca as an
            # extra sample under the canonical La-conditioned protocol.
            assert [x for x in modes['La']['ligands'] if x.get('smiles') != '[La+3]'] == [x for x in modes['Ca']['ligands'] if x.get('smiles') != '[Ca+2]']
        if 'La' in modes:
            matches = [x for x in modes['La']['samples'] if x['matches_canonical_coordinate_fingerprint']]
            assert len(matches) == 1, (base_id, len(matches))
            extras = [x for x in modes['La']['samples'] if not x['matches_canonical_coordinate_fingerprint']]
            proposed = [x['sample_id'] for x in extras[:2]]
        else:
            proposed = []
        rows.append({'panel_id': ref['panel_id'], 'panel_index': ref['panel_index'],
                     'class': ref['class'], 'sequence_group': base_id,
                     'evidence_use': 'known_reference_structural_robustness; not_new_biological_validation',
                     'canonical_source': ref['source_cif'], 'canonical_fingerprint': original,
                     'canonical_core_manifest': ref['manifest'],
                     'modes': modes, 'proposed_two_extra_La_samples': proposed,
                     'same_protein_input_beyond_metal': modes.get('Ca',{}).get('protein_input_including_msa_templates_sha256') == modes.get('La',{}).get('protein_input_including_msa_templates_sha256')})
        print(f'{ref["panel_id"]}: La {len(modes.get("La",{}).get("samples",[]))}, Ca {len(modes.get("Ca",{}).get("samples",[]))}; proposed {proposed}', flush=True)
    counts = Counter(s['role'] for r in rows for m in r['modes'].values() for s in m['samples'])
    result = {'schema_version': 'alquemia.canonical_pqq_existing_fold_inventory.v1',
              'calibration_preparation': pin(calibration_path), 'implementation': pin(__file__),
              'source_archive': str(ARCHIVE), 'selection_rule': 'all available samples retained; proposed first two numeric seed/sample IDs after excluding exact canonical geometry, with no confidence or energy filtering',
              'method_family': 'same AlphaFold3 archive; samples share seed1 and input conditioning within each mode; no independent-fold-method claim',
              'rows': rows, 'counts': dict(counts), 'errors': errors,
              'new_folds': 0, 'new_energies': 0, 'scientific_energy_files_read': 0,
              'baseline_changed': False,
              'limits': ['All labels and proteins were used in calibration/development; only source-conformer robustness is tested.',
                         'Ca-conditioned samples are a separately declared perturbation, never substitutes for same-protocol La-conditioned replicates.',
                         'Model confidence is not an experimental displacement covariance or a thermal population.']}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'pqq_fold_samples_inventory.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    with (OUT/'pqq_fold_samples.tsv').open('w') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['panel_id','label','conditioning','sample_id','role','proposed_extra_La','cif','sha256','summary_confidences'])
        for r in rows:
            for metal,m in r['modes'].items():
                for s in m['samples']:
                    writer.writerow([r['panel_id'],r['class'],metal,s['sample_id'],s['role'],metal=='La' and s['sample_id'] in r['proposed_two_extra_La_samples'],s['cif']['path'],s['cif']['sha256'],s['summary_confidences']['path']])
    print(json.dumps({'counts':dict(counts),'errors':errors,'protein_inputs_equal_except_ligand_count':sum(r['same_protein_input_beyond_metal'] for r in rows)},indent=2),flush=True)


if __name__ == '__main__':
    main()
