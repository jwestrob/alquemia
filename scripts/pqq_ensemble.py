"""Opt-in, declared three-La-fold median of existing standard PQQ scores.

This is a developmental transfer of frozen single-structure bands, not a
thermal ensemble, probability, new calibration, or automatic scoring default.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import statistics

import gemmi
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
import pqq_standard as standard
from compact_solvation_scanner import decision, validate as validate_scanner
from compact_solvation_compare import mix_pair

PROTOCOL = 'PQQ_declared_three_La_fold_median_development_v1'
COMPONENTS = ('native_R_model_kcal_mol', 'GFN2_vacuum_R_kcal_mol',
              'GFN2_ALPB_R_kcal_mol', 'solvation_delta_R_kcal_mol',
              'composite_R_model_kcal_mol')


def resolve_residue(model, selector):
    matches = [r for c in model for r in c
               if (c.name, r.seqid.num, r.seqid.icode.strip(), r.name) ==
               (selector['chain'], selector['resnum'], selector.get('icode', ''), selector['resname'])]
    if len(matches) != 1:
        raise InvalidArtifact('ambiguous or missing declared source residue')
    return matches[0]


def source_identity(case):
    """Verify actual sequence/site roles and selected native source metal."""
    st = gemmi.read_structure(str(verify(case['source_structure'])))
    if len(st) != 1 or case.get('normalization') != 'protenix_generic_PQQ':
        raise InvalidArtifact('three-fold operation requires one-model supported predicted sources')
    if case.get('source_conditioning_metal') != 'La':
        raise InvalidArtifact('declared source conditioning must be La')
    raw = case.get('raw_source_metal', {})
    if raw.get('element') != 'La':
        raise InvalidArtifact('explicit native La source selector required')
    residue = resolve_residue(st[0], raw)
    if len(residue) != 1 or residue[0].name != raw['atom'] or residue[0].element.name != 'La':
        raise InvalidArtifact('actual source metal differs from declared native La')
    # Normalization changes residue aliases, not selected chain/index/atom.
    if any(case['metal'].get(k, '') != raw.get(k, '') for k in ('chain', 'resnum', 'icode', 'atom')):
        raise InvalidArtifact('raw and normalized source metal selectors disagree')
    # Exact-source assembly includes every protein chain, not only chains that
    # contribute a core role. An added/different partner is a different system.
    chains = sorted({c.name for c in st[0]
                     if any(gemmi.find_tabulated_residue(r.name).is_amino_acid() for r in c)})
    sequence = []
    for chain_name in chains:
        cc = [c for c in st[0] if c.name == chain_name]
        if len(cc) != 1:
            raise InvalidArtifact('ambiguous source protein chain')
        residues = [(r.seqid.num, r.seqid.icode.strip(), r.name) for r in cc[0]
                    if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        if not residues:
            raise InvalidArtifact('source protein sequence unavailable')
        sequence.append((chain_name, residues))
    for selector in case['roles'].values():
        resolve_residue(st[0], selector)
    return {'sequence_and_numbering_sha256': cache_key(sequence),
            'sequence_length': sum(len(v) for _, v in sequence),
            'site_roles': case['roles'], 'pqq_selector': case['pqq'],
            'assembly': case['assembly'], 'source_conditioning_metal': 'La'}


def declared_groups(request):
    if request.get('schema_version') != standard.SOURCE_PROTOCOL:
        raise InvalidArtifact('unsupported source request')
    declaration = request.get('ensemble', {})
    if declaration.get('protocol_id') != PROTOCOL:
        raise InvalidArtifact('three-fold membership must be declared in the original source request')
    cases = {c['case_id']: c for c in request['cases']}
    if len(cases) != len(request['cases']):
        raise InvalidArtifact('duplicate source case')
    seen = set(); proteins = set(); groups = []
    for group in declaration.get('groups', []):
        members = group['members']; protein = group['protein_id']
        if (not protein or protein in proteins or len(members) != 3 or len(set(members)) != 3
                or seen.intersection(members) or not set(members).issubset(cases)
                or group.get('source_conditioning_metal') != 'La' or not group.get('selection_rule')):
            raise InvalidArtifact('exactly three distinct, explicitly selected La folds per protein required')
        selected = [cases[cid] for cid in members]
        if (len({c['source_structure']['sha256'] for c in selected}) != 3
                or any(c.get('root_case_id') != protein or c.get('biological_group') != group.get('biological_group')
                       for c in selected) or not group.get('biological_group')):
            raise InvalidArtifact('source structures or declared protein identity differ')
        identities = [source_identity(c) for c in selected]
        if any(i != identities[0] for i in identities[1:]):
            raise InvalidArtifact('actual protein sequence, numbering, site roles, or assembly differ')
        groups.append({**group, 'source_identity': identities[0]})
        seen.update(members); proteins.add(protein)
    if not groups or seen != set(cases):
        raise InvalidArtifact('every source must belong to exactly one declared three-fold group')
    return groups


def prepared_state(case):
    core = read_json(verify(case['core']['parent']))
    ctx = case['representations']['context']; audit = read_json(verify(ctx['preparation']))
    fixed = core['fixed_core']
    state = {'core_protocol': core['protocol_id'], 'source_protocol': core['source_protocol'],
             'context_protocol': ctx['protocol_id'], 'chemistry_policy': case['chemistry_policy'],
             'roles': fixed['requested_roles'], 'included_roles': fixed['included_roles'],
             'water_policy': fixed['water_policy'], 'point_charge_embedding': fixed['point_charge_embedding'],
             'charge_ledger': core['charge_ledger'],
             'fragments': [{k: f.get(k) for k in ('id', 'kind', 'role', 'formal_charge', 'microstate_id')}
                           for f in core['qm_fragments']],
             'core_elements': [a[0] for a in xyz(verify(case['core']['endpoints']['La']['xyz']))],
             'multiplicities': {z: case['core']['endpoints'][z]['multiplicity'] for z in ('Ca', 'La')}}
    if not audit['water_inventory_unchanged'] or not audit['core_source_coordinates_unchanged']:
        raise InvalidArtifact('context changed frozen core or water state')
    # Context membership can legitimately vary under the same geometric rule;
    # retain those composition changes instead of silently treating it as fixed.
    context = {'added_formal_charge': audit['added_formal_charge'],
               'added_fragments': audit['added_fragments'], 'atom_count': audit['new_atom_count'],
               'charges': {z: ctx['endpoints'][z]['charge'] for z in ('Ca', 'La')},
               'preparation': ctx['preparation']}
    return state, context


def aggregate(group, rows, bands):
    members = [rows[cid] for cid in group['members']]
    components = {}
    for field in COMPONENTS:
        values = [r.get(field) for r in members]
        if any(v is not None and (not isinstance(v, (int, float)) or not math.isfinite(v)) for v in values):
            raise InvalidArtifact('nonfinite/non-numeric component score')
        missing = [r['case_id'] for r, v in zip(members, values) if v is None]
        good = not missing
        components[field] = {'status': 'available' if good else 'unavailable',
            'values': values, 'missing_members': missing,
            'median': statistics.median(values) if good else None,
            'min': min(values) if good else None, 'max': max(values) if good else None,
            'spread': max(values) - min(values) if good else None}
    median = components['composite_R_model_kcal_mol']['median']
    return {**group, 'members': members, 'required_members': 3,
            'available_members': sum(r.get('composite_R_model_kcal_mol') is not None for r in members),
            'status': 'available' if median is not None else 'unavailable', 'components': components,
            'median_R_model_kcal_mol': median,
            'developmental_frozen_band_transfer': decision(median, bands) if median is not None else 'unavailable',
            'component_medians_are_not_additive': True,
            'source_domain_status': sorted({r['source_domain_status'] for r in members})}


def collect(plan, output, collection=None):
    p = standard.checked_plan(plan)
    if p['backend'] != 'fast_pqq':
        raise InvalidArtifact('ensemble requires the named fast PQQ standard path')
    request = read_json(verify(p['request'])); groups = declared_groups(request)
    release = standard.release(verify(p['release']))
    reference = read_json(verify(release['artifacts']['calibration']))
    bands = reference['calibration']['context']['bands']
    directory = Path(plan).resolve().parent
    prep_path = directory / 'source_preparation/preparation.json'
    preparation = read_json(prep_path) if prep_path.exists() else None
    prepared = {}; contexts = {}; state = {}
    if preparation:
        if preparation['source_manifest'] != p['request'] or preparation['config'] != request['config']:
            raise InvalidArtifact('preparation differs from the declared source request')
        prepared = {c['case_id']: c for c in preparation['cases']}
        for cid, case in prepared.items():
            if case['status'] == 'prepared': state[cid], contexts[cid] = prepared_state(case)
    for group in groups:
        available_states = [state[cid] for cid in group['members'] if cid in state]
        if any(s != available_states[0] for s in available_states[1:]):
            raise InvalidArtifact('within-protein prepared core chemical state differs')
    actual = {}; component_pin = None
    if collection:
        result = read_json(collection)
        manifest = directory / 'prepared_score/scoring/manifest.json'
        if (result.get('standard_plan') != record(plan) or result.get('fast_mode') != standard.FAST_MODE
                or result['manifest'] != record(manifest) or result['published_context_bands'] != bands):
            raise InvalidArtifact('collection does not belong to the same standard plan/model/reference')
        validate_scanner(manifest)
        component_pin = result['component_collection']; original = read_json(verify(component_pin))
        if original['manifest'] != result['manifest'] or original['published_context_bands'] != bands:
            raise InvalidArtifact('component collection provenance differs')
        originals = {r['case_id']: r for r in original['rows']}
        for row in result['rows']:
            cid = row['case_id']
            if cid in actual or cid not in p['source_domain_status']:
                raise InvalidArtifact('duplicate or undeclared result member')
            if ({k: v for k, v in row.items() if k != 'source_domain_status'} != originals.get(cid)
                    or row['source_domain_status'] != p['source_domain_status'][cid] or cid not in state):
                raise InvalidArtifact('member result/preparation/domain differs from its original collection')
            if row['status'] == 'available':
                values = mix_pair({z: row['native_MACE'][z]['energy_eV'] for z in ('Ca', 'La')},
                    {z: row['low_level_endpoints'][z]['vacuum_hartree'] for z in ('Ca', 'La')},
                    {z: row['low_level_endpoints'][z]['alpb_hartree'] for z in ('Ca', 'La')})
                if any(row[f] != values[f] for f in COMPONENTS):
                    raise InvalidArtifact('component algebra differs from actual endpoint values')
            elif row['composite_R_model_kcal_mol'] is not None:
                raise InvalidArtifact('unavailable endpoint has a score')
            actual[cid] = row
    rows = {}
    for source in request['cases']:
        cid = source['case_id']; prep = prepared.get(cid, {})
        row = actual.get(cid, {'case_id': cid, 'status': 'unavailable',
            'reason': prep.get('reason', 'no matching collected score'),
            **{f: None for f in COMPONENTS}})
        rows[cid] = {**row, 'source_structure': source['source_structure'],
                     'source_domain_status': p['source_domain_status'][cid],
                     'preparation_status': prep.get('status', 'unavailable'),
                     'context_composition': contexts.get(cid)}
    summaries = [aggregate(g, rows, bands) for g in groups]
    result = {'protocol_id': PROTOCOL, 'status': 'complete' if all(g['status'] == 'available' for g in summaries) else 'partial',
        'standard_plan': record(plan), 'request': p['request'], 'implementation': record(__file__),
        'collection': record(collection) if collection else None, 'component_collection': component_pin,
        'source_preparation': record(prep_path) if preparation else None,
        'reference': release['artifacts']['calibration'], 'frozen_bands': bands,
        'groups': summaries, 'group_denominator': len(groups), 'source_denominator': len(rows),
        'band_transfer_status': 'developmental_single_structure_bands_applied_to_three_fold_median',
        'new_molecular_calls_by_this_operation': 0, 'default_changed': False, 'threshold_refitted': False,
        'thermal_ensemble': False, 'probabilities': None, 'broad_biological_validation': False}
    write_new(output, result)
    return {'status': result['status'], 'groups': len(groups), 'result': record(output)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True); parser.add_argument('--collection')
    parser.add_argument('--output', required=True)
    print(json.dumps(collect(**vars(parser.parse_args(argv)))))


if __name__ == '__main__': main()
