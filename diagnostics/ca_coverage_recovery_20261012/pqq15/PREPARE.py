#!/usr/bin/env python3
"""Finite PQQ failure diagnosis and metadata-only recovery; molecular work on Slurm."""
from pathlib import Path
import argparse
import ast
import copy
import hashlib
import json
import os
import resource
import socket
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
PROTOCOL = 'Nikasha_Ca_PQQ15_existing_source_metadata_recovery_v1'
PROTEINS = ('c5axv8', 'q88jh5', 'q9z4j7')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(item):
    if pin(item['path']) != item:
        raise ValueError('pinned bytes changed: ' + item['path'])
    return Path(item['path'])


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False); f.write('\n')


def code_dependencies():
    """Static local import closure; do not import scientific packages on login."""
    pending = ['pqq_fast_prepare', 'motion_envelope', 'normalize_af3_cif', 'carve_with_pqq']
    pins = {}
    while pending:
        name = pending.pop(); path = ROOT/'scripts'/(name+'.py')
        if name in pins or not path.is_file(): continue
        pins[name] = pin(path)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import): pending.extend(n.name.split('.')[0] for n in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module: pending.append(node.module.split('.')[0])
    return pins


def freeze(workspace):
    """Metadata and byte pins only; no scientific imports or coordinate parsing."""
    out = Path(workspace).resolve(); out.mkdir(parents=True, exist_ok=False)
    inventory = ROOT/'diagnostics/orca_parameter_audit_20261012/ca_conditioned_inventory/INVENTORY.json'
    prior = ROOT/'workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json'
    old = ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json'
    sources = {r['case_id']: r for r in read(prior)['cases']}
    previous = read(old); selection_pin = previous['selection']; selection = read(verify(selection_pin))
    groups = {g['selection_id']: g for g in selection['groups']}
    previous_rows = {(r['case_id'], r['selection_id']): r for r in previous['cases']}
    requests = [r for r in read(inventory)['rows'] if r['protein'].split('-pqq-')[0] in PROTEINS]
    expected = {p+'-pqq-la_model__conditioned_Ca__seed-1_sample-'+str(i) for p in PROTEINS for i in range(5)}
    assert {r['case_id'] for r in requests} == expected and len(requests) == 15
    rows = []
    base = ROOT/'workspaces/accommodation_goal_20260920/folds_v1/source_preparation'
    for r in sorted(requests, key=lambda r: r['case_id']):
        cid = r['case_id']; oldrow = sources[cid]
        assert oldrow['status'] == 'unsupported' and oldrow['reason'] == 'exact selected PQQ required'
        group = groups[r['existing_target_selection_id']]
        anchor = previous_rows[group['members'][0], group['selection_id']]
        assert anchor['status'] == 'prepared' and anchor['mapping_eligible']
        files = {name: pin(base/'cases'/cid/name) for name in
                 ('normalized.pdb', 'protonated.pdb', 'protonated_protonation_manifest.json')}
        files['alias'] = pin(base/'adapters'/cid/'source_residue_alias.cif')
        files['alias_receipt'] = pin(base/'adapters'/cid/'adapter.json')
        files['failure'] = pin(base/'failures'/(cid+'.json'))
        verify(oldrow['source']['source_structure'])
        rows.append({'case_id': cid, 'source': oldrow['source'], 'old_failure': oldrow,
                     'files': files, 'group': group, 'anchor': anchor})
    scripts = code_dependencies()
    manifest = {'protocol_id': PROTOCOL, 'driver': pin(__file__), 'plan': pin(HERE/'PLAN.md'),
                'repository': str(ROOT), 'python': pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
                'inventory': pin(inventory), 'source_preparation': pin(prior), 'old_envelope': pin(old),
                'selection': selection_pin, 'config': selection['config'], 'rows': rows,
                'source_code': scripts, 'denominator': 15, 'numerical_threads': 1,
                'new_energy_calls': 0, 'new_optimizations': 0, 'new_protonation': False,
                'output': str(out/'execution_v1'), 'original_results_changed': False}
    write(out/'manifest.json', manifest)
    print(json.dumps({'manifest': pin(out/'manifest.json'), 'sources': 15, 'executed': False}))


def atom_rows(structure):
    rows = []
    if len(structure) != 1:
        raise ValueError('one source model required')
    for chain in structure[0]:
        for res in chain:
            for atom in res:
                rows.append({'chain': chain.name, 'resnum': res.seqid.num, 'icode': res.seqid.icode.strip(),
                             'resname': res.name, 'atom': atom.name, 'element': atom.element.name,
                             'altloc': str(atom.altloc), 'occupancy': atom.occ,
                             'xyz_A': [atom.pos.x, atom.pos.y, atom.pos.z]})
    return rows


def atom_id(a):
    return tuple(a[k] for k in ('chain', 'resnum', 'icode', 'atom', 'element', 'altloc'))


def coordinate_equality(before, after, heavy_only=False):
    if heavy_only:
        before = [a for a in before if a['element'] not in ('H', 'D')]
        after = [a for a in after if a['element'] not in ('H', 'D')]
    a = {atom_id(x): x for x in before}; b = {atom_id(x): x for x in after}
    if len(a) != len(before) or len(b) != len(after) or set(a) != set(b):
        raise ValueError('source atom identities/composition changed')
    error = max((abs(a[k]['xyz_A'][i]-b[k]['xyz_A'][i]) for k in a for i in range(3)), default=0.)
    if error != 0. or any(a[k]['occupancy'] != b[k]['occupancy'] for k in a):
        raise ValueError('source coordinates or occupancy changed')
    return {'atom_count': len(a), 'max_coordinate_difference_A': error, 'identities_and_occupancy_equal': True}


def audit_and_prepare(item, config, output, modules):
    gemmi, prep, envelope, normalizer, fixed = modules
    source = copy.deepcopy(item['source']); cid = item['case_id']; out = output/'cases'/cid
    out.mkdir(parents=True, exist_ok=False)
    result = {'case_id': cid, 'source': source, 'selection_id': item['group']['selection_id'],
              'status': 'unsupported', 'recovery_protocol': PROTOCOL,
              'prior_failure': item['files']['failure'], 'production_changed': False}
    audit = {'case_id': cid, 'source': source['source_structure'], 'files': item['files'],
             'repair_actions': [], 'energy_calls': 0, 'new_hydrogens': False, 'new_waters': False}
    try:
        for rec in item['files'].values(): verify(rec)
        raw = gemmi.read_structure(str(verify(source['source_structure'])))
        oldnorm = gemmi.read_structure(str(verify(item['files']['normalized.pdb'])))
        protonated = gemmi.read_structure(str(verify(item['files']['protonated.pdb'])))
        audit['raw_to_normalized_heavy_atoms'] = coordinate_equality(atom_rows(raw), atom_rows(oldnorm), True)
        audit['raw_to_protonated_heavy_atoms'] = coordinate_equality(atom_rows(raw), atom_rows(protonated), True)
        normalized = out/'normalized_replay.pdb'
        audit['normalizer_replay'] = normalizer.normalize(verify(item['files']['alias']), normalized)
        replay = gemmi.read_structure(str(normalized))
        audit['normalizer_replay_heavy_atoms'] = coordinate_equality(atom_rows(raw), atom_rows(replay), True)
        audit['normalizer_replay_bytes_equal_saved'] = normalized.read_bytes() == verify(item['files']['normalized.pdb']).read_bytes()
        locators = {(source[k]['chain'], source[k]['resnum'], source[k].get('icode', '')) for k in ('metal', 'pqq')}
        before = atom_rows(protonated)
        names = {(c.name, r.seqid.num, r.seqid.icode.strip()): r.name for c in replay[0] for r in c}
        for c in protonated[0]:
            for r in c:
                key = (c.name, r.seqid.num, r.seqid.icode.strip())
                if key in locators and key in names and r.name != names[key]:
                    if names[key] not in ('LA', 'PQQ'):
                        raise ValueError('normalizer did not establish exact supported ligand identity')
                    audit['repair_actions'].append({'kind': 'residue_name_only', 'locator': list(key), 'before': r.name, 'after': names[key]})
                    r.name = names[key]
        ppath = out/'protonated_reused.pdb'
        if audit['repair_actions']:
            protonated.write_pdb(str(ppath))
        else:
            ppath.write_bytes(verify(item['files']['protonated.pdb']).read_bytes())
        audit['all_protonated_atoms_preserved'] = coordinate_equality(before, atom_rows(gemmi.read_structure(str(ppath))))
        model = gemmi.read_structure(str(ppath))[0]; ms = source['metal']; base = fixed.base
        site = base.select_metal_sites(model, site_chain=ms['chain'], site_resnum=ms['resnum'],
               site_icode=ms.get('icode', ''), site_resname=ms['resname'], site_atom=ms['atom'])[0]
        audit['metal'] = {'selector': ms, 'element': site.element, 'xyz_A': list(site.atom.pos)}
        pkey = fixed.residue_key(source['pqq'], 'PQQ')
        pqq = fixed.resolve_residue(model, pkey, 'PQQ')
        prepared = base.prepare_pqq(pqq, fixed.DEFAULT_PQQ_MICROSTATE)
        audit['PQQ'] = {'selector': source['pqq'], 'atom_count': len(pqq), **prepared.metadata(),
                        'nearest_atom_A': min(site.atom.pos.dist(a.pos) for a in pqq),
                        'donor_distances_A': {a.name: site.atom.pos.dist(a.pos) for a in pqq if a.name in prepared.direct_donors}}
        if prepared.formula != 'C14H3N2O8' or prepared.microstate.formal_charge != -3:
            raise ValueError('original exact PQQ state gate failed')
        nearby = base._prepare_nearby_pqq(model, site.atom.pos, fixed.DEFAULT_PQQ_MICROSTATE)
        audit['nearby_PQQ_keys'] = [k.label() for k in nearby]
        audit['cutoffs_A'] = {'PQQ_association': base.DIAGNOSTIC_CUTOFF_A, 'typed_coordination': base.FIRST_SHELL_CUT}
        qualified, included, sulfur, untyped = base._classify_contacts(model, site, nearby)
        audit['typed_CN'] = len(qualified)
        audit['typed_contacts'] = [c.metadata() for c in qualified]
        audit['untyped_ON_within_coordination_cutoff'] = [c.metadata() for c in untyped]
        protein_on = []
        for c in model:
            for r in c:
                if r.name not in base.STANDARD_AA: continue
                for a in r:
                    if a.element.name in ('O', 'N'):
                        protein_on.append({'chain': c.name, 'resname': r.name, 'resnum': r.seqid.num,
                                           'atom': a.name, 'distance_A': site.atom.pos.dist(a.pos)})
        audit['nearest_protein_ON'] = sorted(protein_on, key=lambda r: r['distance_A'])[:10]
        audit['declared_roles'] = {}
        for role, selector in source['roles'].items():
            residue = fixed.resolve_residue(model, fixed.residue_key(selector, role), role)
            distances = [{'atom': a.name, 'element': a.element.name, 'distance_A': site.atom.pos.dist(a.pos)}
                         for a in residue if a.element.name in ('O', 'N')]
            audit['declared_roles'][role] = {'selector': selector,
                'nearest_ON': sorted(distances, key=lambda r: r['distance_A'])[:2]}
        audit['frozen_context'] = {'selection_id': item['group']['selection_id'], 'union': item['group']['union'],
                                 'state_anchor_source': item['group']['members'][0],
                                 'state_anchor_key': item['anchor']['state_key'],
                                 'membership_changed': False, 'compatibility': 'not_evaluated_until_core_admission'}
        if set(nearby) != {pkey}:
            if not nearby and audit['PQQ']['nearest_atom_A'] > base.DIAGNOSTIC_CUTOFF_A:
                result.update(status='requires_new_site_conditioned_structure',
                              reason='intact exact PQQ is beyond unchanged 4.0 A metal association gate; coordinates preserved')
            else:
                result.update(status='unsupported', reason='exact selected PQQ remains missing or ambiguous after canonical normalization replay')
            return result, audit
        # Reuse every actual protonation coordinate and retain its original receipt.
        pm = read(verify(item['files']['protonated_protonation_manifest.json']))
        if pm['repaired_missing_atom_count'] or pm['repaired_missing_terminal_atom_count']:
            raise ValueError('historical heavy atom repair violates original protocol')
        pm.update(original_protonation_manifest=item['files']['protonated_protonation_manifest.json'],
                  output=pin(ppath), source=pin(normalized), recovery_protocol=PROTOCOL,
                  protonation_reused=True, residue_name_repairs=audit['repair_actions'])
        pm_path = out/'protonated_reuse_manifest.json'; write(pm_path, pm)
        core = prep.build_core(source, ppath, pm_path, normalized, out, config, fixed)
        state = prep.parent_state(core, config['topology'], require_endpoint_receipts=False)
        expanded, context = prep.complete_expansion(state)
        cp = out/'context_preparation.json'; write(cp, context); endpoints = {}
        for metal in ('Ca', 'La'):
            xp = out/(metal+'_context.xyz'); prep.write_xyz(xp, expanded[metal])
            charge = core['endpoints'][metal]['charge'] + context['added_formal_charge']
            prep.check_atoms(expanded[metal], charge)
            endpoints[metal] = {'xyz': pin(xp), 'charge': charge, 'multiplicity': 1, 'metal_index': 0}
        row = {'case_id': cid, 'source': source, 'status': 'prepared', 'core': core,
               'representations': {'context': {'protocol_id': prep.CONTEXT_PROTOCOL,
                                   'preparation': pin(cp), 'endpoints': endpoints}}}
        write(out/'SOURCE_PREPARATION.json', row)
        result = envelope.prepare_task(row, item['group'], config, output/'expanded'/item['group']['selection_id'])
        result.update(selection_id=item['group']['selection_id'], source_preparation=pin(out/'SOURCE_PREPARATION.json'),
                      recovery_protocol=PROTOCOL, state_anchor_source=item['group']['members'][0],
                      prior_failure=item['files']['failure'], roles=['Ca_conditioned_fixed_canonical_envelope_transfer'], triple_ids=[])
        if result['status'] == 'prepared_awaiting_group_check':
            anchor = item['anchor']
            if result['state_key'] != anchor['state_key'] or result['protein_key'] != anchor['protein_key']:
                result.update(status='source_state_mismatch', reason='unchanged canonical state-anchor check failed')
            else:
                result['status'] = 'prepared'
        audit['frozen_context']['compatibility'] = result['status']
        return result, audit
    except Exception as exc:
        result.update(status='unsupported', reason=str(exc), exception_type=type(exc).__name__)
        audit['exception'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
        return result, audit


def execute(manifest_path):
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):
        raise RuntimeError('all molecular parsing/preparation requires a real Slurm allocation')
    for k in THREAD_KEYS: os.environ[k] = '1'
    m = read(manifest_path); verify(m['driver']); verify(m['plan']); verify(m['python'])
    for rec in m['source_code'].values(): verify(rec)
    for key in ('inventory', 'source_preparation', 'old_envelope', 'selection'): verify(m[key])
    assert m['denominator'] == len(m['rows']) == 15 and m['numerical_threads'] == 1
    out = Path(m['output']); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic(); started = time.time()
    receipt = {'protocol_id': PROTOCOL, 'manifest': pin(manifest_path), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
               'hostname': socket.gethostname(), 'pid': os.getpid(), 'started_epoch': started,
               'thread_environment': {k: os.environ[k] for k in THREAD_KEYS}, 'scientific_energy_calls': 0,
               'molecular_optimizations': 0, 'new_protonation_calls': 0, 'status': 'running'}
    write(out/'START.json', receipt)
    try:
        sys.path.insert(0, str(ROOT/'scripts'))
        import gemmi
        import pqq_fast_prepare as prep
        import motion_envelope as envelope
        fixed, _, _, normalizer, _ = prep.helpers(m['config'])
        rows = []; audits = []
        for item in m['rows']:
            row, audit = audit_and_prepare(item, m['config'], out, (gemmi, prep, envelope, normalizer, fixed))
            write(out/'cases'/item['case_id']/'AUDIT.json', audit)
            row['recovery_audit'] = pin(out/'cases'/item['case_id']/'AUDIT.json')
            write(out/'cases'/item['case_id']/'ADMISSION.json', row)
            rows.append(row); audits.append(audit)
            print(json.dumps({'case_id': item['case_id'], 'status': row['status'], 'reason': row.get('reason')}), flush=True)
        from collections import Counter
        summary = {'protocol_id': PROTOCOL, 'manifest': pin(manifest_path), 'terminal': True, 'rows': rows,
                   'counts': dict(Counter(r['status'] for r in rows)), 'denominator': 15,
                   'config': m['config'], 'new_energies': 0, 'production_changed': False, 'original_results_changed': False}
        write(out/'AUDIT.json', {'protocol_id': PROTOCOL, 'manifest': pin(manifest_path), 'rows': audits})
        write(out/'PREPARATION.json', summary)
        lines = ['# PQQ15 existing-source recovery', '', '| Source | Status | Nearest PQQ atom, A | Typed CN |', '|---|---|---:|---:|']
        for row, audit in zip(rows, audits):
            lines.append('| '+row['case_id']+' | '+row['status']+' | '+str(audit.get('PQQ', {}).get('nearest_atom_A', 'unavailable'))+' | '+str(audit.get('typed_CN', 'unavailable'))+' |')
        lines += ['', 'All original outcomes and bands preserved. Distances describe source placement; they do not establish biological specificity.',
                  'No atom relocation, new fold, reprotonation, scoring or gate relaxation. Off-site cases require a separately declared structure-regeneration protocol.']
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
        receipt.update(status='complete', counts=summary['counts'], preparation=pin(out/'PREPARATION.json'))
    except Exception as exc:
        receipt.update(status='failed', exception_type=type(exc).__name__, reason=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(finished_epoch=time.time(), wall_seconds=time.monotonic()-start,
                       CPU_seconds=usage.ru_utime+usage.ru_stime, peak_RSS_KiB=usage.ru_maxrss)
        write(out/'EXECUTION.json', receipt)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('freeze'); a.add_argument('--workspace', required=True)
    a = sub.add_parser('execute'); a.add_argument('--manifest', required=True)
    args = p.parse_args()
    if args.command == 'freeze': freeze(args.workspace)
    else: execute(args.manifest)
