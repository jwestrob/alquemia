#!/usr/bin/env python3
"""Combine a new matched LBT pool with the unchanged Mex common5 on one Slurm CPU."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import traceback

SCOPES = ('7CCO', '7CCN', 'global')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path']) == record, 'Changed input: ' + record['path'])
    return json.loads(Path(record['path']).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def run(cfg, result):
    previous = read(cfg['previous_comparison'])
    oldcfg = read(previous['config'])
    helper_pin = previous['implementation']
    require(pin(helper_pin['path']) == helper_pin, 'Previous helper changed')
    spec = importlib.util.spec_from_file_location('previous_monomer_comparison', helper_pin['path'])
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    replay = {}
    helper.run(oldcfg, replay)
    require(previous['status'] == replay['status'] == 'complete'
            and previous['rows'] == replay['rows'], 'Previous complete15 comparison does not replay')
    pool = read(cfg['response_pool'])
    require(pool['status'] in ('complete', 'unavailable') and pool['stage'] == 'pool'
            and pool['protocol_id'] == 'lbt3_matched_fixedFF_response_pool_v1',
            'Response pool execution/admission failed')
    poolcfg = read(pool['config'])
    require(pool['implementation']['sha256'] == poolcfg['implementation']['sha256']
            and pin(pool['implementation']['path']) == pool['implementation']
            and pin(poolcfg['implementation']['path']) == poolcfg['implementation'], 'Pool implementation changed')
    require(poolcfg['plan'] == cfg['plan'] and poolcfg['matched_pool_result'] == oldcfg['lbt_pool'],
            'New pool is not the declared continuation of the same matched model')
    original_lbt = read(oldcfg['lbt_pool'])
    comparison = pool['comparison']
    require(comparison['old_comparison_unchanged'] == original_lbt['comparison'], 'Old LBT menus changed')
    new_menus = comparison['expanded_menus']
    require(len(new_menus) == 3 and {r['source_scope'] for r in new_menus} == set(SCOPES),
            'Three declared source/global scopes required')
    require((pool['status'] == 'complete') == all(r['status'] == 'available' for r in new_menus),
            'Pool/menu availability disagrees')
    old_rows = {(r['source_scope'], r['menu']): r for r in previous['rows']}
    old_lbt = {(r['source_scope'], r['menu']): r for r in original_lbt['comparison']['menus']}
    rows = []
    for scope in SCOPES:
        new = next(r for r in new_menus if r['source_scope'] == scope)
        old = old_rows[(scope, 'union')]
        row = {'source_scope': scope, 'status': 'unavailable',
               'old_R_kcal_mol': old['R_common5_Mex_kcal_mol'],
               'old_LBT_selected': old['LBT_selected'], 'new_menu': new,
               'Mex_selected_unchanged': previous['Mex_common5_selected']}
        rows.append(row)
        if new['status'] != 'available':
            row['reason'] = new.get('reason', 'Required new menu unavailable')
            continue
        require(new['missing_cells'] == [] and new['declared_search_slots'] == (4 if scope == 'global' else 2),
                'Available menu lacks its full declared denominator')
        old_parts = old_lbt[(scope, 'union')]['LBT_Dy_minus_La_components_kcal_mol']
        delta = {k: math.fsum((new['LBT_Dy_minus_La_components_kcal_mol'][k], -old_parts[k]))
                 for k in helper.NAMES}
        helper.close(math.fsum(delta.values()), new['contrast_change_from_old_union_kcal_mol'],
                     'Expanded pool contrast change')
        parts = {k: math.fsum((old['R_common5_components_kcal_mol'][k], -delta[k]))
                 for k in helper.NAMES}
        value = math.fsum(parts.values())
        change = -math.fsum(delta.values())
        row.update(status='available', new_R_kcal_mol=value, R_change_kcal_mol=change,
                   new_R_components_kcal_mol=parts, new_LBT_selected=new['selected'],
                   within_1_kcal_mol_resolution=abs(value) <= 1.,
                   absolute_change_exceeds_resolution=abs(change) > 1.)
    complete = all(r['status'] == 'available' for r in rows)
    result.update(status='complete' if complete else 'incomplete', old15_rows_unchanged=previous['rows'], rows=rows,
                  no_new_molecular_calls=True, no_new_geometries=True,
                  scope='Consumed fixed-state development; finite candidates, not populations or affinity',
                  interpretive_resolution_kcal_mol=1., resolution_is_classifier_or_uncertainty=False)
    if complete:
        old_gap = abs(rows[1]['old_R_kcal_mol'] - rows[0]['old_R_kcal_mol'])
        new_gap = abs(rows[1]['new_R_kcal_mol'] - rows[0]['new_R_kcal_mol'])
        result.update(old_source_gap_kcal_mol=old_gap, new_source_gap_kcal_mol=new_gap,
                      objective_sampling_material=any(r['absolute_change_exceeds_resolution'] for r in rows),
                      finite_objective_test_closed=True, same_search_automatic_extension=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit()
            and int(os.environ['SLURM_CPUS_PER_TASK']) == 1, 'One actual Slurm CPU required')
    cfg = json.loads(args.config.read_text())
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation changed')
    require(pin(cfg['plan']['path']) == cfg['plan'], 'Plan changed')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = {'protocol_id': 'lanm_monomer_matched_response_comparison_v1', 'status': 'running',
              'implementation': pin(__file__), 'config': pin(args.config), 'plan': cfg['plan'],
              'slurm_job_id': os.environ['SLURM_JOB_ID']}
    try:
        run(cfg, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=str(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, process_CPU_seconds=time.process_time()-cpu)
    write(out/'RESULT.json', result)
    if 'rows' in result:
        lines = ['# Matched-objective monomer response comparison', '',
                 'Positive R favors conditional Mex-La plus LBT-Dy. Units: kcal/mol.', '',
                 '| LBT scope | Previous R | New R | Change |', '|---|---:|---:|---:|']
        for row in result['rows']:
            values = [format(row[k], '.6f') if k in row else 'unavailable'
                      for k in ('old_R_kcal_mol', 'new_R_kcal_mol', 'R_change_kcal_mol')]
            lines.append('| ' + row['source_scope'] + ' | ' + ' | '.join(values) + ' |')
        lines += ['', 'All earlier15 menus are retained in RESULT.json. Same completed Mex common5 throughout.',
                  'Finite nonstationary source menus; no affinity, population or independent-protein accuracy claim.']
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
