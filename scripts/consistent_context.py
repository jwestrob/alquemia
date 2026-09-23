"""Opt-in fixed source-fragment union across existing folds; no default changes."""
from __future__ import annotations
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import json
from pathlib import Path
import shutil
import time

import numpy as np
from affordable_common import InvalidArtifact, cache_key, paired, read_json, record, verify, write_new, xyz
from second_shell_context import parent_state, atom_key, cap_key, sidechain
from environment_context_chemistry import cysteine_units
from mace_hybrid import check_atoms, write_xyz

PROTOCOL = 'native_OMOL_GFN2_fixed_group_source_fragment_union_v1'
PILOT_GROUPS = ('a0a3f2yly8-pqq-la_model', 'a0acd6b9f2-pqq-la_model',
                'q9z4j7-pqq-la_model', 'q88jh5-pqq-la_model')
STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def residue_id(meta):
    return (meta['chain'], meta['resnum'], meta.get('insertion_code', ''), meta['resname'])


def source_id(meta):
    return (*residue_id(meta), meta['atom'], meta['element'])


def fragment_id(fragment):
    return (fragment['kind'], *residue_id(fragment['source']))


def fragment_record(key):
    kind, chain, resnum, icode, name = key
    return {'kind': kind, 'source': {'chain': chain, 'resnum': resnum,
            'insertion_code': icode, 'resname': name}}


def forced_expansion(state, fragments):
    """Same complete source graph chemistry; only fragment selection is replaced."""
    g = state['graph']; resolved = set()
    for f in fragments:
        if f['kind'] not in ('sidechain', 'peptide'): raise InvalidArtifact('unsupported union fragment kind')
        resolved.add((f['kind'], g.locate(f['source']).key))
    units, disulfides = cysteine_units(state, resolved)
    selected = set(state['selected']); charges = dict(state['charges']); added = []
    for kind, rkey in sorted(resolved):
        if kind == 'peptide': selected |= g.peptide(rkey); charge = 0
        elif rkey in units:
            selected |= units[rkey]; charge = 0
            for key in units[rkey]: charges[key[:2]] = 0
        else:
            nodes, charge = sidechain(g, rkey); selected |= nodes; charges[rkey] = charge
        added.append({'kind': kind, 'source': g.residues[rkey].source_dict(), 'formal_charge': charge})
    atoms, mapping = g.materialize(selected)
    for a in state['opaque']:
        mapping['source_to_qm'].append(dict(a, qm_index=len(atoms)+1))
        atoms.append(state['original']['La'][a['qm_index']])
    for a in state['maps']:
        if a['kind'] == 'source' and atom_key(a['source']) in state['water_keys']:
            mapping['source_to_qm'].append(dict(a, qm_index=len(atoms)+1))
            atoms.append(state['original']['La'][a['qm_index']])
    caps = {cap_key(a): a['qm_index'] for a in state['maps'] if a['kind'] == 'sigma_link_H'}
    opaque = {(a['fragment'], a['name']): a['qm_index'] for a in state['opaque']}
    index = {0: 0}; new_sources = []
    for a in mapping['source_to_qm']:
        if a['kind'] == 'opaque_cofactor': old = opaque[(a['fragment'], a['name'])]
        elif a['kind'] == 'source':
            old = state['old_sources'].get(atom_key(a['source']))
            if old is None: new_sources.append(a['source'])
        else: old = caps.get(cap_key(a))
        if old is not None: index[old] = a['qm_index']
    if set(state['old_sources'].values())-set(index): raise InvalidArtifact('union lost original core source atoms')
    endpoints = {}
    for z, original in state['original'].items():
        rows = [original[0]]+list(atoms)
        for old, new in index.items(): rows[new] = original[old]
        pos = np.asarray([r[1:] for r in rows]); d = np.linalg.norm(pos[:, None]-pos[None, :], axis=2)
        np.fill_diagonal(d, np.inf)
        if d.min() < .45: raise InvalidArtifact('overlapping union context source atoms/caps')
        endpoints[z] = rows
    if endpoints['Ca'][1:] != endpoints['La'][1:]: raise InvalidArtifact('paired nonmetal source coordinates differ')
    for a in mapping['source_to_qm']: a['xyz_A'] = list(endpoints['La'][a['qm_index']][1:])
    charge_rows = [{'source': g.residues[key].source_dict(), 'formal_charge': value} for key, value in sorted(charges.items())]
    audit = {'protocol_id': PROTOCOL, 'source': state['source'], 'mapping': mapping, 'core_to_context': index,
             'added_fragments': added, 'source_formal_charge_ledger': charge_rows,
             'added_formal_charge': sum(charges.values())-sum(state['charges'].values()),
             'added_source_atoms': new_sources, 'disulfide_closures': disulfides,
             'removed_caps': [a for a in state['maps'] if a['kind'] == 'sigma_link_H' and a['qm_index'] not in index],
             'original_atom_count': len(state['original']['La']), 'new_atom_count': len(endpoints['La']),
             'core_source_coordinates_unchanged': True, 'water_inventory_unchanged': True,
             'paired_differing_nonmetal_indices': [], 'mapping_coordinate_endpoint': 'La',
             'selector': 'all_supported_saved_folds_fragment_identity_union', 'new_protons_or_waters': False}
    return endpoints, audit


def state_signature(state, audit, charges):
    g = state['graph']; atom_ids = []; hydrogen_parents = []; caps = []
    for a in audit['mapping']['source_to_qm']:
        if a['kind'] == 'source':
            atom_ids.append(('source', *source_id(a['source'])))
            if a['source']['element'] in ('H', 'D'):
                key = atom_key(a['source'])
                if key not in g.hparents: raise InvalidArtifact('selected H has no real covalent parent')
                hydrogen_parents.append((source_id(a['source']), source_id(g.meta[g.hparents[key]])))
        elif a['kind'] == 'opaque_cofactor':
            atom_ids.append(('cofactor', a['fragment'], a['name'], a['element']))
        else:
            caps.append((source_id(a['retained']), source_id(a['omitted']), a.get('length_A', 1.09)))
    bonds = sorted(tuple(sorted((source_id(b['a']), source_id(b['b'])))) for b in audit['mapping']['retained_bonds'])
    return {'source_atoms_and_cofactor': sorted(atom_ids), 'H_covalent_parents': sorted(hydrogen_parents),
            'caps': sorted(caps), 'retained_bonds': bonds,
            'formal_charges': sorted((residue_id(a['source']), a['formal_charge']) for a in audit['source_formal_charge_ledger']),
            'endpoint_charges': charges, 'multiplicities': {'Ca': 1, 'La': 1},
            'core_microstates': [(a['id'], a['formal_charge'], a.get('microstate_id')) for a in state['parent']['qm_fragments']],
            'water_inventory': sorted(source_id(g.meta[k]) for k in state['water_keys'])}


def protein_identity(state):
    return [residue_id(r.source_dict()) for r in state['graph'].residues.values() if r.canonical_resname in STANDARD]


def prepare_one(row, group, config, output):
    """Reuse real protonated/source core; replay old selection before union."""
    cid = row['case_id']; out = Path(output)/'cases'/cid
    out.mkdir(parents=True, exist_ok=False)
    result = {'case_id': cid, 'source': row['source'], 'original_status': row['status'],
              'protocol_id': PROTOCOL, 'status': 'unsupported', 'union': group['union']}
    try:
        state = parent_state(row['core'], config['topology'], require_endpoint_receipts=False)
        old = read_json(verify(row['representations']['context']['preparation']))
        replay, _ = forced_expansion(state, old['added_fragments'])
        replay_max = 0.
        for z in ('Ca', 'La'):
            saved = xyz(verify(row['representations']['context']['endpoints'][z]['xyz']))
            if len(saved) != len(replay[z]) or [r[0] for r in saved] != [r[0] for r in replay[z]]:
                raise InvalidArtifact('old fragment-list replay changes composition/order')
            diff = float(np.max(np.abs(np.asarray([a[1:] for a in saved])-np.asarray([a[1:] for a in replay[z]]))))
            if diff > 1e-12: raise InvalidArtifact('old fragment-list replay changes coordinates')
            replay_max = max(replay_max, diff)
        expanded, audit = forced_expansion(state, group['fragments'])
        ep = {}; charges = {}
        for z in ('Ca', 'La'):
            charge = row['core']['endpoints'][z]['charge']+audit['added_formal_charge']; charges[z] = charge
            check_atoms(expanded[z], charge); xp = out/(z+'_context.xyz'); write_xyz(xp, expanded[z])
            ep[z] = {'xyz': record(xp), 'charge': charge, 'multiplicity': 1, 'metal_index': 0}
        paired(verify(ep['La']['xyz']), verify(ep['Ca']['xyz']), ep['La']['charge'], ep['Ca']['charge'])
        write_new(out/'context_preparation.json', audit)
        signature = state_signature(state, audit, charges); identity = protein_identity(state)
        write_new(out/'state_signature.json', signature); write_new(out/'protein_identity.json', identity)
        result.update(status='prepared_awaiting_group_check', state_signature=record(out/'state_signature.json'),
                      state_key=cache_key(signature), protein_identity=record(out/'protein_identity.json'),
                      protein_key=cache_key(identity), old_context=row['representations']['context'],
                      original_core=row['core'], original_selection_replay_max_A=replay_max,
                      representations={'context': {'protocol_id': PROTOCOL,
                        'preparation': record(out/'context_preparation.json'), 'endpoints': ep}},
                      atom_count=len(expanded['La']), cap_count=len(audit['mapping']['cut_bonds_and_caps']),
                      original_atom_count=old['new_atom_count'], original_caps=len(old['mapping']['cut_bonds_and_caps']),
                      original_La_charge=row['representations']['context']['endpoints']['La']['charge'],
                      new_La_charge=charges['La'])
    except Exception as exc:
        result.update(reason=str(exc), exception_type=type(exc).__name__)
    write_new(out/'preparation_result.json', result)
    return result


def prepare(preparation, agreement, output, workers=4):
    start = time.monotonic(); old = read_json(preparation); out = Path(output).resolve()
    if old['denominator'] != 250 or len(old['cases']) != 250: raise InvalidArtifact('fixed250 source inventory required')
    source_manifest = read_json(verify(old['source_manifest']))
    if {r['case_id'] for r in old['cases']} != {r['case_id'] for r in source_manifest['cases']}:
        raise InvalidArtifact('original source denominator differs')
    out.mkdir(parents=True, exist_ok=False); grouped = defaultdict(list)
    for row in old['cases']: grouped[row['source']['root_case_id']].append(row)
    if len(grouped) != 25: raise InvalidArtifact('expected25 proteins')
    groups = []; pending = []
    for gid, rows in sorted(grouped.items()):
        if len(rows) != 10 or {r['source']['source_conditioning_metal'] for r in rows} != {'Ca', 'La'}:
            raise InvalidArtifact('expectedfive sources for each folding metal')
        canonical = [r['case_id'] for r in rows if r['source']['canonical_coordinate_match']]
        if len(canonical) != 1: raise InvalidArtifact('one original canonical match perprotein required')
        selected = defaultdict(list)
        for row in rows:
            if row['status'] != 'prepared': continue
            audit = read_json(verify(row['representations']['context']['preparation']))
            for f in audit['added_fragments']: selected[fragment_id(f)].append(row['case_id'])
        fragments = [fragment_record(f) for f in sorted(selected)]
        union = {'root_case_id': gid, 'fragments': fragments,
                 'contributors': [{'fragment': fragment_record(k), 'source_cases': sorted(set(v))} for k, v in sorted(selected.items())],
                 'source_case_ids': sorted(r['case_id'] for r in rows), 'canonical_case_id': canonical[0],
                 'rule': 'union_of_prior_supported_source_fragment_identities_before_new_energies'}
        up = out/'groups'/gid/'UNION.json'; write_new(up, union)
        group = {'root_case_id': gid, 'canonical_case_id': canonical[0], 'union': record(up), 'fragments': fragments}
        groups.append(group)
        pending.extend((row, group) for row in rows if row['status'] == 'prepared')
    # Freeze ALL group selections before preparing/scoring any new endpoint.
    write_new(out/'GROUP_SELECTION.json', {'groups': groups, 'source_preparation': record(preparation), 'agreement': record(agreement)})
    results = {}
    with ProcessPoolExecutor(max_workers=int(workers)) as pool:
        futures = {pool.submit(prepare_one, row, group, old['config'], str(out)): row['case_id'] for row, group in pending}
        for future in as_completed(futures):
            row = future.result(); results[row['case_id']] = row
            print(json.dumps({'case_id': row['case_id'], 'status': row['status'], 'reason': row.get('reason')}), flush=True)
    rows = []
    for original in old['cases']:
        row = results.get(original['case_id'])
        if row is None:
            row = {'case_id': original['case_id'], 'source': original['source'], 'status': 'prior_unsupported',
                   'reason': original.get('reason'), 'original_failure': original}
        else:
            group = next(g for g in groups if g['root_case_id'] == row['source']['root_case_id'])
            canonical = results.get(group['canonical_case_id'])
            if row['status'] == 'prepared_awaiting_group_check':
                if canonical is None or 'state_key' not in canonical:
                    row.update(status='canonical_union_state_unavailable', reason='no canonical group state; no replacement selected')
                elif row['protein_key'] != canonical['protein_key']:
                    row.update(status='protein_identity_mismatch', reason='actual source sequence/numbering differs from canonical')
                elif row['state_key'] != canonical['state_key']:
                    row.update(status='selected_state_mismatch', reason='source atoms/H parents/caps/bonds/charges/microstate differ from canonical')
                    own = read_json(verify(row['state_signature'])); ref = read_json(verify(canonical['state_signature']))
                    row['mismatching_fields'] = [k for k in own if own[k] != ref[k]]
                else: row['status'] = 'prepared'
            row['canonical_signature'] = canonical.get('state_signature') if canonical else None
        row['pilot_selected'] = bool(row['source']['canonical_coordinate_match'] or row['source']['root_case_id'] in PILOT_GROUPS)
        rows.append(row)
    result = {'protocol_id': PROTOCOL, 'source_preparation': record(preparation), 'source_manifest': old['source_manifest'],
              'agreement': record(agreement), 'config': old['config'], 'groups': groups, 'cases': rows,
              'denominator': 250, 'supported': sum(r['status'] == 'prepared' for r in rows),
              'pilot_denominator': sum(r['pilot_selected'] for r in rows),
              'pilot_supported': sum(r['pilot_selected'] and r['status'] == 'prepared' for r in rows),
              'new_molecular_calls': 0, 'wall_seconds': time.monotonic()-start,
              'source_protonation_reused': True, 'new_geometry_optimization': False, 'production_changed': False}
    write_new(out/'PREPARATION.json', result)
    return {'preparation': record(out/'PREPARATION.json'), **{k: result[k] for k in ('supported', 'denominator', 'pilot_denominator', 'pilot_supported')}}


def prepare_crystals(preparation, agreement, output):
    old = read_json(preparation); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    rows = []
    for original in old['cases']:
        if original['case_id'] not in ('1H4I', '4MAE', '1KB0'): continue
        row = copy.deepcopy(original); cid = row['case_id']
        row['source'].update(root_case_id=cid, source_conditioning_metal=None, canonical_coordinate_match=False,
                             primary_evaluation_pool=False, source_kind='crystal_singleton')
        audit = read_json(verify(row['representations']['context']['preparation']))
        fragments = [fragment_record(k) for k in sorted({fragment_id(f) for f in audit['added_fragments']})]
        up = out/'groups'/cid/'UNION.json'; write_new(up, {'fragments': fragments, 'source_case_ids': [cid], 'singleton': True})
        result = prepare_one(row, {'union': record(up), 'fragments': fragments}, old['config'], out)
        if result['status'] == 'prepared_awaiting_group_check':
            result['status'] = 'prepared'; result['canonical_signature'] = result['state_signature']
            if result['atom_count'] != result['original_atom_count'] or result['new_La_charge'] != result['original_La_charge']:
                raise InvalidArtifact('singleton union changed original chemistry')
            for z in ('Ca', 'La'):
                a = xyz(verify(result['representations']['context']['endpoints'][z]['xyz']))
                b = xyz(verify(row['representations']['context']['endpoints'][z]['xyz']))
                if not same_geometry(a, b): raise InvalidArtifact('singleton union changed coordinates')
        rows.append(result)
    if {r['case_id'] for r in rows} != {'1H4I', '4MAE', '1KB0'}: raise InvalidArtifact('fixed3 crystal controls absent')
    result = {'protocol_id': PROTOCOL, 'source_preparation': record(preparation), 'agreement': record(agreement),
              'cases': rows, 'denominator': 3, 'config': old['config'], 'new_molecular_calls': 0}
    write_new(out/'CRYSTALS.json', result)
    return {'manifest': record(out/'CRYSTALS.json'), 'supported': sum(r['status'] == 'prepared' for r in rows)}


def same_geometry(a, b):
    return (len(a) == len(b) and [r[0] for r in a] == [r[0] for r in b] and
            np.max(np.abs(np.asarray([r[1:] for r in a])-np.asarray([r[1:] for r in b]))) <= 1e-12)


def reusable_state(new, old):
    return (new['charge'] == old['charge'] and new['multiplicity'] == old['multiplicity'] and
            same_geometry(xyz(verify(new['xyz'])), xyz(verify(old['xyz']))))


def score_manifest(preparation, agreement, numerical_qualification, output, crystals,
                   stage='calibration28', reuse_inventory=None, reuse_collection=None, frozen_reference=None,
                   reuse_union_collection=None):
    """Finite pilot inputs for existing warmMACE and native ORCA executors only."""
    import mace_omol as omol
    import compact_solvation as solvent
    p = read_json(preparation); qualification = read_json(numerical_qualification); crystal = read_json(crystals)
    qm = read_json(verify(qualification['manifest']))
    if (p['protocol_id'] != PROTOCOL or p['pilot_denominator'] != 61 or
            qualification['protocol_id'] != 'adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1' or
            not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete']):
        raise InvalidArtifact('preparation or actual MaxIter500 qualification differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = {}; dest = out/'implementation'; dest.mkdir()
    for script in Path(__file__).parent.glob('*.py'):
        target = dest/script.name; shutil.copyfile(script, target); impl[script.name] = record(target)
    if crystal['protocol_id'] != PROTOCOL or crystal['denominator'] != 3 or crystal['config'] != p['config']:
        raise InvalidArtifact('singleton crystal/configuration mismatch')
    if stage not in ('calibration28', 'pilot_transfers36', 'primary_transfers225'): raise InvalidArtifact('undeclared scoring stage')
    if stage != 'calibration28':
        if frozen_reference is None:
            raise InvalidArtifact('transfer scoring requires a valid frozen canonical reference')
        ref = read_json(frozen_reference)
        if (ref.get('protocol_id') != PROTOCOL or ref.get('status') != 'available' or ref.get('available') != 25
                or not ref.get('bands') or ref.get('transfers_used') is not False):
            raise InvalidArtifact('incomplete/overlapping reference cannot authorize transfer expansion')
        cal = read_json(verify(ref['collection'])); cal_input = read_json(verify(cal['inventory']))
        if cal_input['preparation'] != record(preparation):
            raise InvalidArtifact('reference belongs to a different union preparation')
    if stage == 'calibration28':
        declared = [r for r in p['cases'] if r['source']['canonical_coordinate_match']]+crystal['cases']
    else:
        declared = [r for r in p['cases'] if r['source']['primary_evaluation_pool'] and
                    (stage == 'primary_transfers225' or r['pilot_selected'])]
    denominator = {'calibration28': 28, 'pilot_transfers36': 36, 'primary_transfers225': 225}[stage]
    if len(declared) != denominator: raise InvalidArtifact('stage denominator changed')
    archived = read_json(reuse_inventory) if reuse_inventory else None
    old_collection = read_json(reuse_collection) if reuse_collection else None
    if archived and (archived['model'] != p['config']['model'] or archived['software'] != p['config']['software']):
        raise InvalidArtifact('archive checkpoint/software differs')
    if old_collection and (archived is None or old_collection['inventory'] != record(reuse_inventory)):
        raise InvalidArtifact('archived solvent collection inventory differs')
    union_prior = read_json(reuse_union_collection) if reuse_union_collection else None
    if union_prior:
        prior_input = read_json(verify(union_prior['inventory']))
        if (union_prior['protocol_id'] != PROTOCOL or prior_input['preparation'] != record(preparation)
                or prior_input['reference'] != record(frozen_reference)):
            raise InvalidArtifact('prior union scores have another preparation/reference')
    cases = []; failures = []
    for r in declared:
        if r['status'] != 'prepared':
            failures.append({'case_id': r['case_id'], 'source': r['source'], 'status': r['status'], 'reason': r.get('reason')}); continue
        c = {k: r['source'][k] for k in ('case_id', 'root_case_id', 'biological_group', 'expected_class',
                 'label_scope', 'evidence_stratum', 'canonical_coordinate_match', 'source_conditioning_metal', 'primary_evaluation_pool')}
        c.update(role='calibration' if c['canonical_coordinate_match'] else
                 'crystal_transfer' if r['source'].get('source_kind') == 'crystal_singleton' else 'structural_transfer',
                 sequence_accession_group=c['root_case_id'], representations=r['representations'],
                 group_state_signature=r['state_signature'], old_context=r['old_context'])
        cases.append(c)
    inventory = {'protocol_id': PROTOCOL, 'stage': stage, 'cases': cases, 'failures': failures, 'denominator': denominator,
                 'model': p['config']['model'], 'software': p['config']['software'],
                 'preparation': record(preparation), 'agreement': record(agreement), 'crystals': record(crystals),
                 'reference': record(frozen_reference) if frozen_reference else None,
                 'union_reuse_collection': record(reuse_union_collection) if union_prior else None}
    ip = out/'SCORE_INPUTS.json'; write_new(ip, inventory)
    mace = {'protocol_id': omol.PROTOCOL, 'stage': PROTOCOL, 'preparation': record(preparation),
            'agreement': record(agreement), 'software': inventory['software'], 'model': inventory['model'],
            'implementation': impl, 'tasks': [], 'reused': {}, 'verification_references': {}}
    tasks = []; reused = {}; reuse_checks = {}
    for c in cases:
        rep = c['representations']['context']
        alias = c['root_case_id'] if c['canonical_coordinate_match'] else c['case_id']
        previous = next((r for r in archived['cases'] if r['case_id'] == alias), None) if archived else None
        union_row = next((r for r in union_prior['rows'] if r['case_id'] == c['case_id']), None) if union_prior else None
        for z, ep in rep['endpoints'].items():
            mt = {'task_id': c['case_id']+'__context__'+z, 'case_id': c['case_id'], 'representation': 'context',
                  'metal': z, 'metal_index': 0, 'kind': 'core', 'xyz': ep['xyz'], 'charge': ep['charge'],
                  'spin_multiplicity': ep['multiplicity'], 'energy_component': omol.COMPONENT}
            mt['cache_key'] = cache_key({'task': mt, 'model': mace['model'], 'software': mace['software'], 'implementation': impl})
            mace['tasks'].append(mt)
            old_ep = previous['representations']['context']['endpoints'][z] if previous else None
            native_source = record(reuse_inventory) if previous else None
            if union_row and union_row.get('native_endpoints', {}).get(z, {}).get('status') == 'complete':
                old_ep = union_row['native_endpoints'][z]; native_source = record(reuse_union_collection)
            if old_ep and reusable_state(ep, old_ep):
                from compact_solvation_compare import native_endpoint
                native = native_endpoint(read_json(verify(old_ep['native_MACE_receipt'])), mace['model'])
                if not reusable_state(ep, native): raise InvalidArtifact('archive scalar task geometry differs')
                mace['reused'][mt['task_id']] = {'receipt': old_ep['native_MACE_receipt'], 'source': native_source,
                    'new_xyz': ep['xyz'], 'archived_xyz': old_ep['xyz'], 'coordinate_tolerance_A': 1e-12,
                    'model_state_and_coordinates_verified': True}
            for medium in ('vacuum', 'alpb'):
                tid = mt['task_id']+'__'+medium+'__native'; directory = out/'solvent'/'tasks'/tid; directory.mkdir(parents=True)
                xp = directory/'core.xyz'; shutil.copyfile(verify(ep['xyz']), xp)
                inp = directory/'endpoint.inp'
                body = solvent.input_text(ep['charge'], ep['multiplicity'], medium, 'native').replace('%scf\n', '%scf\n MaxIter 500\n')
                inp.write_text(body)
                task = {'task_id': tid, 'case_id': c['case_id'], 'case': c['case_id'], 'representation': 'context',
                        'metal': z, 'medium': medium, 'solver': 'native', 'charge': ep['charge'], 'multiplicity': ep['multiplicity'],
                        'source_endpoint': ep, 'source_preparation': rep['preparation'], 'xyz': record(xp),
                        'input': record(inp), 'output_path': str(directory/'endpoint.out')}
                task['scientific_key'] = cache_key({'xyz_sha256': ep['xyz']['sha256'], 'charge': ep['charge'],
                     'multiplicity': ep['multiplicity'], 'medium': medium, 'solver': 'native', 'method': solvent.METHOD,
                     'input_body': body, 'orca': qm['orca']})
                tasks.append(task)
                old_end = None
                if old_collection:
                    old_row = next((r for r in old_collection['rows'] if r['case_id'] == alias and
                         r['representation'] == 'context' and r['metal'] == z and r['solver'] == 'native' and r['status'] == 'complete'), None)
                    if old_row: old_end = old_row['endpoints'][medium]
                if union_row and union_row.get('solvent_endpoints', {}).get(z, {}).get(medium, {}).get('status') == 'complete':
                    old_end = union_row['solvent_endpoints'][z][medium]
                if old_end:
                    old_mp = verify(old_end['manifest']); old_m = read_json(old_mp)
                    old_t = next(t for t in old_m['all_tasks'] if t['task_id'] == old_end['task_id'])
                    old_body = verify(old_t['input']).read_text()
                    if (old_m['method_id'] == solvent.METHOD and old_m['orca'] == qm['orca'] and reusable_state(task, old_t) and
                            old_body.replace(' MaxIter 500\n', '') == body.replace(' MaxIter 500\n', '')):
                        receipt = solvent.completed(old_mp, old_t['task_id'])
                        if receipt and receipt['energy_hartree'] == old_end['energy_hartree']:
                            reused[tid] = receipt
                            reuse_checks[tid] = {'old_task_id': old_t['task_id'], 'old_scientific_key': old_t['scientific_key'],
                                'new_scientific_key': task['scientific_key'], 'coordinate_tolerance_A': 1e-12,
                                'only_permitted_recipe_difference': 'qualified_MaxIter500_ceiling_on_already_converged_native_state',
                                'actual_old_receipt_retained': True}
    sm = {'protocol_id': solvent.PROTOCOL, 'method_id': solvent.METHOD, 'inventory': record(ip), 'agreement': record(agreement),
          'orca': qm['orca'], 'implementation': impl, 'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8},
          'execution_policy': {'task_runner': impl['run_orca_task_manifest.py'], 'runtime_renderer': impl['render_orca_runtime_input.py']},
          'selection': {'cases': [c['case_id'] for c in cases], 'representations': ['context'], 'solver': 'native'},
          'all_tasks': tasks, 'tasks': [t for t in tasks if t['task_id'] not in reused], 'reused': reused,
          'reuse_checks': reuse_checks, 'numerical_qualification': record(numerical_qualification),
          'GFN2_maxiter': 500, 'reference': None, 'baseline_changed': False}
    write_new(out/'mace'/'manifest.json', mace); write_new(out/'solvent'/'manifest.json', sm)
    dry = validate_solvent(out/'solvent'/'manifest.json'); write_new(out/'DRY_RUN.json', dry)
    ready = {'protocol_id': PROTOCOL, 'preparation': record(preparation), 'inventory': record(ip),
             'MACE_manifest': record(out/'mace'/'manifest.json'), 'GFN2_manifest': record(out/'solvent'/'manifest.json'),
             'stage': stage, 'stage_denominator': denominator, 'prepared_cases': len(cases), 'failures': failures,
             'MACE_calls': len(mace['tasks'])-len(mace['reused']), 'MACE_reused': len(mace['reused']),
             'GFN2_calls': len(sm['tasks']), 'GFN2_reused': len(reused), 'launched': False,
             'original_core_or_geometry_changed': False, 'new_membership_protocol': True}
    write_new(out/'READY.json', ready)
    return ready


def validate_solvent(manifest):
    import compact_solvation as solvent
    from affordable_workflow import dry_run
    m = read_json(manifest); inv = read_json(verify(m['inventory'])); verify(m['agreement']); verify(m['orca'])
    for pin in m['implementation'].values(): verify(pin)
    q = read_json(verify(m['numerical_qualification']))
    if not q['both_controls_pass'] or not q['formerly_failed_cell_complete'] or m['GFN2_maxiter'] != 500:
        raise InvalidArtifact('native numerical qualification unavailable')
    wanted = {(c['case_id'], z, medium) for c in inv['cases'] for z in ('Ca', 'La') for medium in ('vacuum', 'alpb')}
    actual = {(t['case_id'], t['metal'], t['medium']) for t in m['all_tasks']}
    if actual != wanted or len(actual) != len(m['all_tasks']): raise InvalidArtifact('score cell denominator mismatch')
    for t in m['all_tasks']:
        rep = next(c for c in inv['cases'] if c['case_id'] == t['case_id'])['representations']['context']; ep = rep['endpoints'][t['metal']]
        if t['source_endpoint'] != ep or t['source_preparation'] != rep['preparation'] or not reusable_state(t, ep):
            raise InvalidArtifact('declared source/state changes')
        body = solvent.input_text(t['charge'], t['multiplicity'], t['medium'], 'native').replace('%scf\n', '%scf\n MaxIter 500\n')
        if verify(t['input']).read_text() != body: raise InvalidArtifact('native recipe changed')
        key = cache_key({'xyz_sha256': ep['xyz']['sha256'], 'charge': ep['charge'], 'multiplicity': ep['multiplicity'],
             'medium': t['medium'], 'solver': 'native', 'method': solvent.METHOD, 'input_body': body, 'orca': m['orca']})
        if key != t['scientific_key']: raise InvalidArtifact('scientific key mismatch')
        if t['task_id'] in m['reused']:
            old = m['reused'][t['task_id']]; oldp = verify(old['manifest']); oldm = read_json(oldp)
            prior = next(x for x in oldm['all_tasks'] if x['task_id'] == old['task_id'])
            if (solvent.completed(oldp, old['task_id']) != old or oldm['method_id'] != solvent.METHOD or
                    oldm['orca'] != m['orca'] or not reusable_state(t, prior) or
                    verify(prior['input']).read_text().replace(' MaxIter 500\n', '') != body.replace(' MaxIter 500\n', '')):
                raise InvalidArtifact('archive reuse failed physical/numerical equivalence check')
    if m['tasks'] != [t for t in m['all_tasks'] if t['task_id'] not in m['reused']]: raise InvalidArtifact('reuse/task partition changed')
    if m['tasks']: dry_run(manifest)
    return {'status': 'validated', 'manifest': record(manifest), 'tasks': len(m['tasks']), 'reused': len(m['reused']), 'energies_run': 0}


def join_archive_collections(collections, inventory, output):
    """Combine disjoint real old collections for reuse; no relabeled receipts."""
    import compact_solvation as solvent
    rows = []; seen = set(); pins = []
    for path in collections:
        m = read_json(path)
        if m['inventory'] != record(inventory) or m['protocol_id'] != solvent.PROTOCOL:
            raise InvalidArtifact('archive collection source/protocol differs')
        pins.append(record(path))
        for row in m['rows']:
            key = (row['case_id'], row['representation'], row['metal'], row['solver'])
            if key in seen: raise InvalidArtifact('overlapping archive collection cells')
            seen.add(key); rows.append(row)
    result = {'protocol_id': solvent.PROTOCOL, 'inventory': record(inventory), 'rows': rows,
              'source_collections': pins, 'new_molecular_calls': 0}
    write_new(output, result)
    return {'rows': len(rows), 'collection': record(output)}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare')
    for k in ('preparation', 'agreement', 'output'): q.add_argument('--'+k, required=True)
    q.add_argument('--workers', type=int, default=4)
    q = sub.add_parser('prepare-crystals')
    for k in ('preparation', 'agreement', 'output'): q.add_argument('--'+k, required=True)
    q = sub.add_parser('score-manifest')
    for k in ('preparation', 'agreement', 'numerical-qualification', 'output', 'crystals'): q.add_argument('--'+k, required=True)
    q.add_argument('--stage', choices=('calibration28', 'pilot_transfers36', 'primary_transfers225'), default='calibration28')
    for k in ('reuse-inventory', 'reuse-collection', 'frozen-reference', 'reuse-union-collection'): q.add_argument('--'+k)
    q = sub.add_parser('validate-solvent'); q.add_argument('--manifest', required=True)
    q = sub.add_parser('join-archive-collections'); q.add_argument('--collections', nargs='+', required=True)
    q.add_argument('--inventory', required=True); q.add_argument('--output', required=True)
    a = vars(p.parse_args()); op = a.pop('op')
    result = globals()[op.replace('-', '_')](**a)
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
