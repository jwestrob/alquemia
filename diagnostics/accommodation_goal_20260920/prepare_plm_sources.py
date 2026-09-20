"""Map existing PLM fold records to source requests; no molecular preparation/run."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import gemmi
from affordable_common import InvalidArtifact, read_json, record, verify, write_new


def prepare(candidate_manifest, release_sources, targets, output):
    candidates = read_json(candidate_manifest)
    config = read_json(release_sources)['config']
    index = {c['target_id']: c for c in candidates['targets']}
    if not targets or len(set(targets)) != len(targets):
        raise InvalidArtifact('explicit distinct target IDs required')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    groups = []
    for tid in targets:
        c = index[tid]
        provenance = c['source_provenance']
        selection = read_json(verify(provenance['selection_review']))
        native = read_json(verify(provenance['native_AF3_input']))
        proteins = [x['protein'] for x in native['sequences'] if 'protein' in x]
        ligands = [x['ligand'] for x in native['sequences'] if 'ligand' in x]
        if (len(proteins) != 1 or proteins[0]['id'] != 'A'
                or sorted((x['id'], tuple(x['ccdCodes'])) for x in ligands) != [('B', ('LA',)), ('C', ('PQQ',))]
                or native['modelSeeds'] != [101]):
            raise InvalidArtifact('PLM source assembly/input differs: ' + tid)
        models = sorted(selection['models'], key=lambda m: m['sample'])
        if [m['sample'] for m in models] != [0, 1, 2] or any(m['seed'] != 101 for m in models):
            raise InvalidArtifact('exact three declared PLM samples required: ' + tid)
        sequence = proteins[0]['sequence']; cases = []; audits = []
        for m in models:
            raw = verify(m['source_cif']); structure = gemmi.read_structure(str(raw))
            if len(structure) != 1 or sorted(ch.name for ch in structure[0]) != ['A', 'B', 'C']:
                raise InvalidArtifact('actual source assembly differs: ' + str(raw))
            chain = structure[0]['A']
            if ([r.seqid.num for r in chain] != list(range(1, len(sequence) + 1))
                    or any(r.seqid.icode.strip() for r in chain)
                    or ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in chain) != sequence):
                raise InvalidArtifact('actual source sequence/numbering differs: ' + str(raw))
            for role, selector in c['roles'].items():
                residue = chain[selector['resnum'] - 1]
                if selector['chain'] != 'A' or selector['icode'] or residue.name != selector['resname']:
                    raise InvalidArtifact('source role identity differs: ' + role)
            metal_residues = list(structure[0]['B']); pqq_residues = list(structure[0]['C'])
            if (len(metal_residues) != 1 or metal_residues[0].name != 'LA'
                    or metal_residues[0].seqid.num != 1 or len(metal_residues[0]) != 1
                    or metal_residues[0][0].name != c['metal']['atom']
                    or metal_residues[0][0].element.name != 'La'
                    or len(pqq_residues) != 1 or pqq_residues[0].name != 'PQQ'
                    or pqq_residues[0].seqid.num != 1):
                raise InvalidArtifact('actual metal/PQQ source identity differs: ' + str(raw))
            cid = f'{tid}_AF3_sample{m["sample"]}'
            cases.append({'case_id': cid, 'root_case_id': tid, 'biological_group': tid,
                'source_structure': m['source_cif'], 'normalization': 'protenix_generic_PQQ',
                'assembly': {'kind': 'exact_source_structure', 'selection': None},
                'metal': c['metal'], 'pqq': c['pqq'], 'roles': c['roles'],
                'source_conditioning_metal': 'La', 'raw_source_metal': {**c['metal'], 'element': 'La'},
                'role': 'development_geometry_diagnostic', 'expected_class': None,
                'label_scope': 'PQQ_prediction', 'evidence_stratum': 'unresolved_PLM_prediction',
                'all_evidence_consumed': True,
                'source_provenance': {'native_AF3_input': provenance['native_AF3_input'],
                    'selection_review': provenance['selection_review'], 'summary_confidence': m['summary_confidence'],
                    'source_seed': 101, 'source_sample': m['sample']}})
            audits.append({'case_id': cid, 'raw_source': m['source_cif'], 'sequence_identical_to_input': True,
                'roles_identical_to_prior_selected_source': True, 'raw_metal_verified': 'La',
                'previous_admission_pass': m.get('admission_pass'),
                'previous_admission_used_to_select_members': False})
        request = {'schema_version': 'fixed_core_PQQ_source_to_complete_context_v1', 'config': config,
            'cases': cases, 'ensemble': {'protocol_id': 'PQQ_declared_three_La_fold_median_development_v1',
                'groups': [{'protein_id': tid, 'biological_group': tid, 'source_conditioning_metal': 'La',
                    'selection_rule': 'All saved seed101 samples0,1,2 from the existing fold selection record; no score or admission filtering.',
                    'members': [c['case_id'] for c in cases]}]}}
        path = out / (tid + '_sources.json'); write_new(path, request)
        groups.append({'target_id': tid, 'request': record(path), 'source_count': 3,
            'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'prior_selected_source': c['case_id'], 'source_audits': audits,
            'chemical_preparation_status': 'not_run', 'scoring_status': 'not_run', 'biological_label': None})
    result = {'schema': 'PLM_existing_fold_source_requests_v1', 'candidate_manifest': record(candidate_manifest),
        'released_source_configuration': record(release_sources), 'implementation': record(__file__),
        'groups': groups, 'new_fold_or_molecular_calls': 0, 'production_results_changed': False,
        'scope': 'Existing-input compatibility only; no protonation, carve, energy, label or admission-policy change.'}
    write_new(out / 'inventory.json', result)
    return {'groups': len(groups), 'raw_sources_verified': sum(g['source_count'] for g in groups),
            'inventory': record(out / 'inventory.json'), 'molecular_execution': 'not_run'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('candidate-manifest', 'release-sources', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--targets', nargs='+', required=True)
    print(json.dumps(prepare(**vars(parser.parse_args())), indent=2))
