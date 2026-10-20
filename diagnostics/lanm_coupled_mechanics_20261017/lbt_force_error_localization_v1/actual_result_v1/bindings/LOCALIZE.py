#!/usr/bin/env python3
"""Saved vacuum force-error localization; standard library only; Slurm only."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import time
import traceback

PROTOCOL = 'lbt_force_error_localization_v1'
PARENT_SHA = 'f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790'
INVENTORY_SHA = 'e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957'
IDS = tuple(s+'_'+t for s in ('7CCO', '7CCN') for t in ('before_La', 'after_La', 'before_Dy', 'after_Dy'))
PAIRS = tuple((s, g) for s in ('7CCO', '7CCN') for g in ('before', 'after'))
GROUPS = ('selected_metal', 'hydrogens', 'backbone_heavy', 'remaining_protein_heavy')
ABS, REL, TOP = 1e-10, 1e-12, 10


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as h:
        digest = hashlib.file_digest(h, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': digest}


def verify(record, checked):
    require(pin(record['path']) == record, 'File pin differs: '+record['path'])
    checked[record['path']] = record
    return Path(record['path'])


def read(record, checked):
    return json.loads(verify(record, checked).read_text())


def write(path, value):
    with path.open('x') as h:
        json.dump(value, h, indent=2, allow_nan=False)
        h.write('\n')


def vectors(value):
    require(isinstance(value, list) and len(value) == 256
            and all(isinstance(v, list) and len(v) == 3
                    and all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in v)
                    for v in value), 'Finite full256x3 vectors required')
    return value


def prepare(cfg, checked):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation'] == pin(__file__)
            and cfg['parent_result']['sha256'] == PARENT_SHA and cfg['inventory']['sha256'] == INVENTORY_SHA
            and cfg['bookkeeping_absolute_tolerance'] == ABS and cfg['bookkeeping_relative_tolerance'] == REL
            and cfg['top_n'] == TOP, 'Frozen scope or bookkeeping constants differ')
    verify(cfg['implementation'], checked); verify(cfg['specification'], checked)
    parent = read(cfg['parent_result'], checked); inv = read(cfg['inventory'], checked)
    pcfg = read(parent['config'], checked); verify(parent['implementation'], checked)
    require(parent['status'] == 'complete' and parent['saved_input_admission'] == 'all_eight_passed_before_model_load'
            and parent['inventory'] == cfg['inventory'] and pcfg['inventory'] == cfg['inventory']
            and pcfg['implementation'] == parent['implementation']
            and [r['cell_id'] for r in parent['rows']] == list(IDS)
            and [r['cell_id'] for r in inv['cells']] == list(IDS)
            and all(r['status'] == 'available' for r in parent['rows'])
            and [(r['source_id'], r['geometry']) for r in parent['matched_metal_differences']] == list(PAIRS)
            and all(r['status'] == 'available' for r in parent['matched_metal_differences']),
            'Complete eight-endpoint/four-pair parent required')
    endpoints = {}
    for slot, meta in zip(parent['rows'], inv['cells']):
        report = read(slot['force_report'], checked); atoms = read(meta['physical'], checked)
        verify(meta['positions'], checked)
        for k in ('MACE_force_receipt', 'MACE_force_task'):
            verify(report[k], checked)
        expected_receipt = slot.get('MACE', meta['MACE']['receipt'])
        expected_task = slot.get('supplement_task', meta['MACE']['task'])
        require(report['status'] == 'available' and report['cell_id'] == meta['cell_id']
                and report['source'] == meta and report['MACE_force_receipt'] == expected_receipt
                and report['MACE_force_task'] == expected_task
                and report['units'] == 'kcal/mol/Angstrom' and report['quantity'] == 'force=-gradient'
                and report['force_accuracy_qualified'] is False, 'Force-report authority differs')
        ids = [a['id'] for a in atoms]
        require(len(ids) == len(set(ids)) == 256 and ids[-1] == 'metal'
                and ids == meta['physical_ids'] == report['physical_ids']
                and [a['physical_id'] for a in report['force_table']] == ids
                and [a['index'] for a in report['force_table']] == list(range(256)),
                'Physical atom order differs')
        positions = vectors([a['xyz_A'] for a in atoms])
        errors = vectors([a['MACE_minus_PBE0']['vacuum'] for a in report['force_table']])
        for a in atoms:
            require(a['id'] == 'metal' or (a['kind'] in ('protein_source', 'terminal_completion')
                    and a['element'] in ('H', 'C', 'N', 'O')),
                    'Unexpected atom outside fixed category inventory')
        endpoints[meta['cell_id']] = dict(report=report, pin=slot['force_report'], atoms=atoms,
                                        positions=positions, ids=ids, errors=errors)
    directions = {}
    for source in ('7CCO', '7CCN'):
        d = read(inv['parents'][source]['directions'], checked)
        require(d['physical_ids'] == endpoints[source+'_before_La']['ids'], 'Direction IDs differ')
        directions[source] = {k: vectors(v) for k, v in d['directions'].items()}
    for source, geometry in PAIRS:
        a, b = (endpoints[source+'_'+geometry+'_'+m] for m in ('La', 'Dy'))
        require(a['ids'] == b['ids'] and a['positions'] == b['positions']
                and [(v['id'], v['element'], v.get('resname')) for v in a['atoms']]
                == [(v['id'], v['element'], v.get('resname')) for v in b['atoms']],
                'Matched physical metadata differs')
        saved = next(r for r in parent['matched_metal_differences']
                     if (r['source_id'], r['geometry']) == (source, geometry))
        require(saved['physical_ids'] == a['ids'], 'Matched result IDs differ')
        vectors(saved['force_vectors']['MACE_minus_PBE0']['vacuum'])
    return parent, endpoints, directions


def check(label, actual, expected, checks):
    tolerance = ABS+REL*max(abs(actual), abs(expected))
    residual = actual-expected
    entry = dict(check=label, actual=actual, expected=expected, residual=residual,
                 tolerance=tolerance, passed=math.isfinite(residual) and abs(residual) <= tolerance)
    checks.append(entry)
    require(entry['passed'], 'Bookkeeping mismatch: '+label)


def squared(v):
    return math.fsum(x*x for x in v)


def statistics(values, total=None):
    norms2 = [squared(v) for v in values]; ss = math.fsum(norms2)
    return dict(atom_count=len(values), sum_squared_norms=ss,
                RMS_atom_vector_kcal_mol_A=math.sqrt(ss/len(values)) if values else None,
                maximum_atom_kcal_mol_A=math.sqrt(max(norms2)) if values else None,
                fraction_of_full_squared_norm=ss/total if total else None)


def category(atom):
    if atom['id'] == 'metal': return 'selected_metal'
    if atom['element'] == 'H': return 'hydrogens'
    if atom['id'].rsplit('/', 1)[-1] in ('N', 'CA', 'C', 'O', 'OXT'): return 'backbone_heavy'
    return 'remaining_protein_heavy'


def localize(source, geometry, parent, endpoints, directions, checks):
    a, b = (endpoints[source+'_'+geometry+'_'+m] for m in ('La', 'Dy'))
    la, dy = a['errors'], b['errors']
    common = [[(x+y)/2 for x, y in zip(u, v)] for u, v in zip(la, dy)]
    difference = [[y-x for x, y in zip(u, v)] for u, v in zip(la, dy)]
    saved = next(r for r in parent['matched_metal_differences']
                 if (r['source_id'], r['geometry']) == (source, geometry))
    expected = saved['force_vectors']['MACE_minus_PBE0']['vacuum']
    for i, (v, old) in enumerate(zip(difference, expected)):
        for axis, (x, y) in enumerate(zip(v, old)):
            check('difference/%s/%s' % (i, axis), x, y, checks)
    stats = {k: statistics(v) for k, v in [('La', la), ('Dy', dy), ('common', common), ('difference', difference)]}
    for metal, endpoint in [('La', a), ('Dy', b)]:
        for key in ('RMS_atom_vector_kcal_mol_A', 'maximum_atom_kcal_mol_A'):
            check('endpoint/'+metal+'/'+key, stats[metal][key],
                  endpoint['report']['metrics']['MACE_minus_PBE0']['vacuum'][key], checks)
    oldmetrics = saved['metrics']['MACE_minus_PBE0']['vacuum']
    for key in ('RMS_atom_vector_kcal_mol_A', 'maximum_atom_kcal_mol_A'):
        check('matched/'+key, stats['difference'][key], oldmetrics[key], checks)
    for name, direction in directions[source].items():
        projection = math.fsum(x*y for v, d in zip(difference, direction) for x, y in zip(v, d))
        check('matched/projection/'+name, projection, oldmetrics['projections_kcal_mol_A'][name], checks)
    left = stats['La']['sum_squared_norms']+stats['Dy']['sum_squared_norms']
    right = 2*stats['common']['sum_squared_norms']+0.5*stats['difference']['sum_squared_norms']
    check('squared_norm_identity', left, right, checks)
    records = []
    for i, atom in enumerate(a['atoms']):
        records.append(dict(index=i, physical_id=atom['id'], element=atom['element'], resname=atom.get('resname'),
            group=category(atom), residue='metal' if atom['id'] == 'metal' else atom['id'].rsplit('/', 1)[0],
            e_La=la[i], e_Dy=dy[i], common=common[i], difference=difference[i],
            norms_kcal_mol_A={k: math.sqrt(squared(v[i])) for k, v in
                             [('e_La', la), ('e_Dy', dy), ('common', common), ('difference', difference)]}))
    summaries = {}
    for name, values in [('common', common), ('difference', difference)]:
        total = stats[name]['sum_squared_norms']
        summaries[name] = dict(full=statistics(values, total),
            groups={g: statistics([values[i] for i, r in enumerate(records) if r['group'] == g], total) for g in GROUPS},
            residues={r: statistics([values[i] for i, a in enumerate(records) if a['residue'] == r], total)
                      for r in sorted({a['residue'] for a in records})},
            top10=sorted(records, key=lambda r: (-r['norms_kcal_mol_A'][name], r['physical_id']))[:TOP])
    return dict(source_id=source, geometry=geometry, status='available', units='kcal/mol/Angstrom',
        squared_norm_units='(kcal/mol/Angstrom)^2', definition='e=MACE-PBE0 vacuum; common=(eLa+eDy)/2; difference=eDy-eLa',
        physical_ids=a['ids'], endpoint_force_reports={'La': a['pin'], 'Dy': b['pin']},
        original_endpoint_metrics={m: ep['report']['metrics'] for m, ep in [('La', a), ('Dy', b)]},
        original_matched_metrics=saved['metrics'], atoms=records, summaries=summaries,
        bookkeeping_checks=checks, zero_norm_fraction_policy='null when the full squared norm is zero')


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--config', required=True); p.add_argument('--output', required=True)
    args = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm allocation required before numerical arithmetic')
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    result = dict(protocol=PROTOCOL, status='unavailable', implementation=pin(__file__), config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), declared_endpoints=list(IDS),
        endpoints=[dict(cell_id=x, status='unavailable', reason='not_admitted') for x in IDS],
        pairs=[dict(source_id=s, geometry=g, status='unavailable', reason='not_attempted') for s, g in PAIRS],
        new_MACE_calls=0, new_DFT_calls=0, new_ALPB_calls=0, new_GFN_calls=0, new_MM_calls=0,
        model_loads=0, geometry_evaluations=0, optimizer_steps=0, parameter_fitting=False,
        force_accuracy_qualified=False, physical_state_qualified=False, preference_claim=False)
    checked = {}
    try:
        cfg = read(result['config'], checked); parent, endpoints, directions = prepare(cfg, checked)
        result.update(parent_result=cfg['parent_result'], inventory=cfg['inventory'], saved_input_admission='all_eight_passed',
                      inherited_caveats=parent['inherited_reference_caveats'])
        result['endpoints'] = [dict(cell_id=x, status='available', force_report=endpoints[x]['pin'],
            original_metrics=endpoints[x]['report']['metrics']) for x in IDS]
        for index, (source, geometry) in enumerate(PAIRS):
            checks = []; result['pairs'][index]['bookkeeping_checks'] = checks
            try: result['pairs'][index] = localize(source, geometry, parent, endpoints, directions, checks)
            except Exception as exc:
                result['pairs'][index].update(reason='bookkeeping_failed', failure=repr(exc))
                raise
        result['status'] = 'complete'
    except Exception as exc:
        result.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
    result.update(verified_inputs=list(checked.values()), elapsed_seconds=time.monotonic()-start,
        bookkeeping_tolerances={'absolute': ABS, 'relative': REL, 'rule': 'abs(a-b)<=absolute+relative*max(abs(a),abs(b))'},
        interpretation='Force-error localization only. Squared-force fractions are not energies or a unique interaction cause. No fitting, solvent validation, physical-state/force qualification, or preference claim.')
    write(out/'RESULT.json', result)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
