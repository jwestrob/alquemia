"""Source-backed preparation for the fixed-core compact PQQ scorer.

Uses existing normalization, seeded protonation, fragment and context policies.
Explicit selectors are required; no biological label is used to select chemistry.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import re
import sys
import time

import gemmi
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import check_atoms, write_xyz
from second_shell_context import parent_state
from environment_context_chemistry import complete_expansion, POLICY_ID

ROOT = Path(__file__).resolve().parents[1]
CORE_PROTOCOL = 'pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
CONTEXT_PROTOCOL = 'native_OMOL_complete_polar_context_disulfide_v2'
SOURCE_PROTOCOL = 'fixed_core_PQQ_source_to_complete_context_v1'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def helpers(config):
    cp = verify(config['calibration_implementation_pins']); pins = read_json(cp)
    fixed = load('fast_pqq_fixed_core', verify(pins['experiment_helpers']['fixed_core_carver']))
    fixed.verify_pins(cp)
    protonator = load('fast_pqq_standard_protonator', verify(pins['experiment_helpers']['protonate_standard_only']))
    crystal = load('fast_pqq_crystal_normalizer', verify(config['crystal_normalizer']))
    normalizer = load('fast_pqq_protenix_normalizer', verify(pins['canonical_helpers']['normalize_af3_cif']))
    return fixed, protonator, crystal, normalizer, cp


def selector_dict(value, fixed):
    if isinstance(value, dict):
        return value
    return fixed.parse_residue_selector(value, 'source role')


def inventory(root, output):
    root = Path(root).resolve()
    compact_path = root / 'diagnostics/compact_solvation_20260920/INVENTORY.json'
    compact = read_json(compact_path)
    context_path = root / 'workspaces/environment_pqq_20260919/prepared_v1/manifest.json'
    old = read_json(context_path)
    static = read_json(root / 'diagnostics/second_shell_20260919/CONFIG.json')
    cp = root / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/implementation_pins.json'
    config = {'calibration_implementation_pins': record(cp),
              'crystal_normalizer': record(root / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/prepare_holdouts.py'),
              'topology': static['topology'], 'software': compact['software'], 'model': compact['model']}
    fixed, _, crystal, _, _ = helpers(config)
    cases = []
    for source in compact['cases']:
        if source['label_scope'] != 'canonical_PQQ_functional_class':
            continue
        cid = source['case_id']; oldcase = next(c['source'] for c in old['cases'] if c['source']['case'] == cid)
        parent = read_json(verify(oldcase['parent'])); normal = read_json(verify(parent['normalization_manifest']))
        core = parent['fixed_core']; roles = {k: selector_dict(v, fixed) for k, v in core['requested_roles'].items()}
        if source['role'] == 'calibration':
            mode = 'protenix_generic_PQQ'; metal = parent['site_selection']['requested']
            pqq = fixed.parse_residue_selector(parent['qm_fragments'][0]['id'], 'PQQ')
            assembly = {'kind': 'exact_source_structure', 'selection': None}
        else:
            mode = 'selected_crystal_chain'
            source_metal = normal['selection']['metal_normalization']['source']
            ms = crystal.parse_atom_selector(normal['selection']['metal_normalization']['normalized'])
            metal = {'chain': ms.residue.chain, 'resname': ms.residue.resname, 'resnum': ms.residue.resnum,
                     'icode': ms.residue.icode, 'atom': ms.atom}
            pqq = fixed.parse_residue_selector(normal['source_site_audit']['pqq_selector'], 'PQQ')
            assembly = {'kind': 'explicit_selected_crystal_chain', 'selection': normal['selection'],
                        'raw_metal': source_metal, 'raw_pqq': normal['source_site_audit']['pqq_selector'],
                        'model': '1', 'chain': 'A'}
        cases.append({'case_id': cid, 'source_structure': normal['source'], 'normalization': mode,
                      'assembly': assembly, 'metal': metal, 'pqq': pqq, 'roles': roles,
                      'role': source['role'], 'biological_group': source['biological_group'],
                      'expected_class': source['expected_class'], 'evidence_stratum': source['evidence_stratum'],
                      'label_scope': source['label_scope'], 'all_evidence_consumed': True,
                      'archived_core': oldcase['parent'], 'archived_normalized': normal['output'],
                      'archived_protonated': read_json(verify(parent['protonation_manifest']))['output'],
                      'archived_context': source['representations']['context']})
    if len(cases) != 28 or sum(c['role'] == 'calibration' for c in cases) != 25:
        raise InvalidArtifact('frozen28-case scope differs')
    result = {'schema_version': SOURCE_PROTOCOL, 'config': config, 'cases': cases,
              'compact_inventory': record(compact_path), 'source_context_manifest': record(context_path),
              'new_calls': {'native_MACE': 56, 'native_GFN2': 112}, 'baseline_changed': False}
    write_new(output, result)
    return {'cases': len(cases), 'manifest': record(output)}


def coordinate_comparison(actual, old):
    a, b = xyz(actual), xyz(old)
    if len(a) != len(b) or [r[0] for r in a] != [r[0] for r in b]:
        return {'exact': False, 'reason': 'composition or ordering differs', 'max_displacement_A': None}
    diff = np.asarray([r[1:] for r in a]) - np.asarray([r[1:] for r in b])
    return {'exact': a == b, 'max_displacement_A': float(np.max(np.abs(diff))),
            'changed_atom_indices': [i for i, (x, y) in enumerate(zip(a, b)) if x != y]}


def build_core(source, protonated, protonation, normalized, output, config, fixed):
    model = gemmi.read_structure(str(protonated))[0]; base = fixed.base
    metal = source['metal']; roles = source['roles']
    if set(roles) != set(fixed.ROLE_ORDER):
        raise InvalidArtifact('exact five homologous core roles required')
    site = base.select_metal_sites(model, site_chain=metal['chain'], site_resnum=metal['resnum'],
        site_icode=metal.get('icode', ''), site_resname=metal['resname'], site_atom=metal['atom'])[0]
    pp = base._prepare_nearby_pqq(model, site.atom.pos, fixed.DEFAULT_PQQ_MICROSTATE)
    base._reject_unsupported_nearby_species(model, site, pp)
    pkey = fixed.residue_key(source['pqq'], 'PQQ')
    if set(pp) != {pkey}:
        raise InvalidArtifact('exact selected PQQ required')
    pqq, prepared = pp[pkey]
    if prepared.formula != 'C14H3N2O8' or prepared.microstate.formal_charge != -3:
        raise InvalidArtifact('unsupported PQQ state')
    keys = {role: fixed.residue_key(roles[role], role) for role in fixed.ROLE_ORDER}
    if any(k.resname not in fixed.ROLE_ALLOWED_RESNAMES[role] for role, k in keys.items()):
        raise InvalidArtifact('unsupported homologous core residue')
    cat, plus = keys['catalytic_aspartate'], keys['extra_acidic_ligand_homolog']
    if cat.chain != plus.chain or plus.resnum != cat.resnum + 2:
        raise InvalidArtifact('source homolog mapping violates fixed D+2 policy')
    residues = {role: fixed.resolve_residue(model, key, role) for role, key in keys.items()}
    partner = fixed.asp_partner_geometry(residues['catalytic_aspartate'], residues['catalytic_asp_cationic_partner'])
    acidic = plus.resname in ('ASP', 'GLU')
    included = [r for r in fixed.ROLE_ORDER if r != 'extra_acidic_ligand_homolog' or acidic]
    included_keys = {pkey, *(keys[r] for r in included)}
    contacts, _, _, _ = base._classify_contacts(model, site, pp)
    if (len(contacts) < 6 or sum(c.element.upper() == 'N' for c in contacts) > 2
            or any(c.residue not in included_keys for c in contacts)):
        raise InvalidArtifact('typed donor/core compatibility gate failed')
    fragments = [{'kind': 'fixed_core_pqq', 'role': 'pqq_cofactor', 'id': pkey.label(),
                  'formal_charge': -3, 'atom_count': len(prepared.atoms),
                  'atom_records': fixed.pqq_atom_records(pqq, prepared), 'microstate_id': fixed.DEFAULT_PQQ_MICROSTATE}]
    atoms = list(prepared.atoms)
    for role in included:
        if role == 'catalytic_asp_cationic_partner':
            rows, audit = fixed.cationic_sidechain_fragment(residues[role], keys[role], site.atom.pos)
        else:
            rows, audit = fixed.canonical_sidechain_fragment(residues[role], keys[role], role, site.atom.pos)
        atoms.extend(rows); fragments.append(audit)
    scaffold = sum(f['formal_charge'] for f in fragments); eps = {}; outputs = {}
    for z, number in (('Ca', 2), ('La', 3)):
        charge = scaffold + number; rows = [(z, site.atom.pos.x, site.atom.pos.y, site.atom.pos.z), *atoms]
        check_atoms(rows, charge); xp = output / (z + '_core.xyz')
        fixed.write_xyz(xp, stem=source['case_id'], label=z, atoms=rows, charge=charge)
        ip = output / (z + '_DFT_reference.inp'); fixed.write_orca_input(ip, xyz_name=xp.name, charge=charge, stem=source['case_id'], label=z)
        eps[z] = {'xyz': record(xp), 'input': record(ip), 'charge': charge, 'multiplicity': 1}
        outputs.update({z + '_xyz': record(xp), z + '_input': record(ip)})
    parent = {'protocol_id': CORE_PROTOCOL, 'fixed_core': {'policy_id': fixed.CORE_POLICY_ID,
                'requested_roles': roles, 'included_roles': included, 'catalytic_Asp_partner_geometry': partner,
                'water_policy': 'dry_exclude_all_source_and_synthetic_waters', 'point_charge_embedding': False},
              'source_structure': record(protonated), 'protonation_manifest': record(protonation),
              'normalized_source_structure': record(normalized), 'source_cif': source['source_structure'],
              'qm_fragments': fragments, 'charge_ledger': {'Ca_total': eps['Ca']['charge'], 'La_total': eps['La']['charge']},
              'outputs': outputs, 'assembly': source['assembly'], 'source_protocol': SOURCE_PROTOCOL}
    path = output / 'core_preparation.json'; write_new(path, parent)
    return {'case': source['case_id'], 'parent': record(path), 'endpoints': eps}


def prepare_case(source, config, output):
    start = time.monotonic(); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    fixed, protonator, crystal, normalizer, cp = helpers(config)
    raw = verify(source['source_structure']); normalized = out / 'normalized.pdb'
    if source['normalization'] == 'protenix_generic_PQQ':
        selection = normalizer.normalize(raw, normalized)
    elif source['normalization'] == 'selected_crystal_chain':
        structure = gemmi.read_structure(str(raw)); assembly = source['assembly']
        model = crystal.exact_model(structure, assembly['model']); chain = crystal.exact_chain(model, assembly['chain'])
        _, selection = crystal.write_selected_normalized_structure(structure, model, chain,
            crystal.parse_atom_selector(assembly['raw_metal']), crystal.parse_residue_selector(assembly['raw_pqq']), normalized)
    else:
        raise InvalidArtifact('unsupported source normalization/assembly policy')
    normal_time = time.monotonic() - start
    protonated = out / 'protonated.pdb'
    protonation = protonator.protonate_standard_only(normalized, protonated, pins_path=cp, ph=7.)
    pm = read_json(protonation)
    if pm['repaired_missing_atom_count'] or pm['repaired_missing_terminal_atom_count']:
        raise InvalidArtifact('heavy repair unsupported by frozen source protocol')
    htime = time.monotonic() - start - normal_time
    core = build_core(source, protonated, protonation, normalized, out, config, fixed)
    state = parent_state(core, config['topology'], require_endpoint_receipts=False)
    expanded, audit = complete_expansion(state)
    ap = out / 'context_preparation.json'; write_new(ap, audit)
    endpoints = {}; comparisons = {}
    for z in ('Ca', 'La'):
        xp = out / (z + '_context.xyz'); write_xyz(xp, expanded[z])
        charge = core['endpoints'][z]['charge'] + audit['added_formal_charge']; check_atoms(expanded[z], charge)
        endpoints[z] = {'xyz': record(xp), 'charge': charge, 'multiplicity': 1, 'metal_index': 0}
        if source.get('archived_context'):
            comparisons[z] = coordinate_comparison(xp, verify(source['archived_context']['endpoints'][z]['xyz']))
            comparisons[z]['charge_identical'] = charge == source['archived_context']['endpoints'][z]['charge']
    result = {'case_id': source['case_id'], 'role': source.get('role', 'prediction'),
              'biological_group': source.get('biological_group'), 'expected_class': source.get('expected_class'),
              'label_scope': source.get('label_scope', 'PQQ_prediction'),
              'all_evidence_consumed': source.get('all_evidence_consumed', False),
              'source': source, 'representations': {'context': {'protocol_id': CONTEXT_PROTOCOL,
                'preparation': record(ap), 'endpoints': endpoints}},
              'core': core, 'selection': selection, 'source_preparation_seconds': time.monotonic() - start,
              'normalization_seconds': normal_time, 'protonation_seconds': htime,
              'normalized_bytes_equal_archive': normalized.read_bytes() == verify(source['archived_normalized']).read_bytes() if source.get('archived_normalized') else None,
              'context_comparison': comparisons, 'chemistry_policy': POLICY_ID,
              'preparation_reused': False, 'status': 'prepared' if all(v['exact'] and v['charge_identical'] for v in comparisons.values()) else 'archive_coordinate_mismatch'}
    write_new(out / 'preparation_result.json', result)
    return result


def prepare(manifest, output):
    source = read_json(manifest)
    if source.get('schema_version') != SOURCE_PROTOCOL:
        raise InvalidArtifact('unsupported source preparation schema')
    ids = [case['case_id'] for case in source['cases']]
    if (not ids or len(ids) != len(set(ids))
            or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', cid) for cid in ids)):
        raise InvalidArtifact('source case IDs must be unique, nonempty safe path components')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    rows = []; start = time.monotonic()
    for case in source['cases']:
        try:
            row = prepare_case(case, source['config'], out / case['case_id'])
        except Exception as exc:
            row = {'case_id': case['case_id'], 'status': 'unsupported', 'reason': str(exc), 'exception_type': type(exc).__name__}
            write_new(out / case['case_id'] / 'failure.json', row)
        rows.append(row); print(json.dumps({k: row.get(k) for k in ('case_id', 'status', 'reason')}), flush=True)
    result = {'source_manifest': record(manifest), 'config': source['config'], 'cases': rows,
              'denominator': len(source['cases']), 'supported': sum(r['status'] == 'prepared' for r in rows),
              'fresh_preparation_seconds': time.monotonic() - start, 'preparation_reused': False}
    write_new(out / 'preparation.json', result)
    return {'supported': result['supported'], 'denominator': result['denominator'], 'preparation': record(out / 'preparation.json')}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('inventory'); q.add_argument('--root', required=True); q.add_argument('--output', required=True)
    q = sub.add_parser('prepare'); q.add_argument('--manifest', required=True); q.add_argument('--output', required=True)
    a = vars(p.parse_args()); op = a.pop('op'); print(json.dumps(globals()[op](**a)))


if __name__ == '__main__': main()
