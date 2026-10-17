#!/usr/bin/env python3
"""Saved 36 LBT + 10 Mex electronic substitution; one Slurm CPU, no model calls.
--config CONFIG --output FRESH
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lanm_monomer_fixed_pool_electronic_substitution_v1'
REPLAY_SHA = '024e2477145ed94a135c5d214e5c9ea9aa4ac79f5a30c0be1e86ef4f182a5407'
POOL_SHA = 'c9856d48a0fd32dde3df620d097eb304eb9ba81c0c4e9afd76465abb44f79452'
XTB_SHA = '481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4'
PARAM_SHA = 'f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545'
HA = 627.509474
METALS = ('La', 'Dy')
PARTS = ('GFN2_electronic', 'MM_ledger', 'fixed_FF_ALPB')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path']).resolve(strict=True)
    require(pin(p)['sha256'] == record['sha256'], 'Changed pin: '+str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def gfn_admission(record, runtime, metal, atoms, charge, xyz_pin, kind, evidence_kind):
    """Rebind accepted raw evidence; use no old solvent/composite energies."""
    r = read(record)
    t = r['task']
    require(r['status'] == 'accepted' and t['metal'] == metal and t['atom_count'] == atoms
            and t['charge'] == charge and t['multiplicity'] == 1 and t['uhf'] == 0
            and t['physical_multiplicity'] == (1 if metal == 'La' else 6)
            and t['oxidation_state'] == 3, 'GFN state differs')
    require(t['xyz']['sha256'] == xyz_pin['sha256'], 'GFN exact XYZ differs')
    verify(t['xyz']); verify(xyz_pin)
    if kind == 'vacuum':
        require(t['medium'] == 'vacuum' and t['metal_index'] == 255
                and t['representation'] == 'whole', 'LBT GFN domain differs')
    else:
        require(t['kind'] == 'embedded' and t['pointcharge_count'] == 1840,
                'Mex embedded domain differs')
    single_point = kind == 'vacuum' and evidence_kind == 'cross_cell'
    if single_point:
        require(t.get('gradient_requested') is False and r.get('gradient_requested') is False,
                'Archived cross-cell energy-only declaration differs')
    else:
        require(evidence_kind in ('origin', 'terminal_point', 'Mex')
                and t.get('gradient_requested', True) is True,
                'Archived gradient declaration differs')
    command = [runtime['executable']['path'], 'core.xyz', '--gfn', '2', '--chrg', str(charge),
               '--uhf', '0', '--etemp', '300', '--acc', '0.0001', '--iterations', '500',
               '--norestart', '--json', '--pop', '--sp' if single_point else '--grad']
    if kind == 'embedded':
        command += ['--input', 'control.inp']
    execution = read(r['execution'])
    require(t['command'] == execution['command'] == command and execution['returncode'] == 0
            and execution['owned_process_group_empty'], 'GFN execution/settings differ')
    for p in r['artifacts'].values():
        verify(p)
    task = read(r['artifacts']['TASK.json'])
    # Mex origins include their preparation-task pin; LBT adds its TASK pin after writing.
    require(task == t or task == {k: v for k, v in t.items() if k != 'task_record'},
            'GFN task receipt differs')
    obs = r['solver_observation']
    require(obs['status'] == 'available' and obs['solver_evidence']['convergence_declared']
            and r['charge_closure_e'] <= 1e-5, 'GFN convergence/charge admission differs')
    native = read(r['xtbout'])
    require(native['xtb version'].startswith('6.7.1') and native['method'] == 'GFN2-xTB'
            and native['number of electrons'] == t['valence_electrons']
            and native['number of unpaired electrons'] == 0, 'Observed GFN method/state differs')
    require(math.isfinite(r['energy_hartree'])
            and abs(native['total energy']-r['energy_hartree']) <= 1e-7,
            'GFN native energy receipt differs')
    return r


def inputs(cfg, out, result):
    prior = read(cfg['prior_comparison'])
    require(prior['status'] == 'complete' and prior['slurm_job_id'] == '12736'
            and prior['new_molecular_calls'] == 0, 'Complete12736 required')
    pcfg = read(prior['config'])
    verify(prior['implementation']); verify(pcfg['implementation'])
    require(prior['implementation']['sha256'] == pcfg['implementation']['sha256'] == REPLAY_SHA,
            'Frozen replay implementation differs')
    replay_helper = module('electronic_saved_replay', prior['implementation'])
    replay = {}
    replay_helper.run(pcfg, replay)
    for k, v in replay.items():
        require(v == prior[k], 'Exact12736 replay differs: '+k)
    old = read(pcfg['lbt_pool']); mex = read(pcfg['mex_common5'])
    mcfg = read(old['config'])
    verify(old['implementation']); verify(mcfg['implementation'])
    require(old['implementation']['sha256'] == mcfg['implementation']['sha256'] == POOL_SHA,
            'Frozen LBT admission helper differs')
    pool_helper = module('electronic_matched_pool_replay', old['implementation'])
    d = out/'old_LBT_admission'; d.mkdir()
    h, cells, _, _, _, _, basis = pool_helper.inputs(mcfg, d)
    solvent = read(cfg['GFN_origin_result']); scfg = read(solvent['config'])
    origin_cfg = read(read(mcfg['matched_origin_result'])['config'])
    require(solvent['status'] == 'complete' and solvent['slurm_job_id'] == '12510'
            and solvent['accepted'] == 8
            and solvent['preparation'] == scfg['preparation'] == origin_cfg['preparation'],
            'GFN origin parent differs')
    runtime = scfg['runtime']
    require(runtime['executable']['sha256'] == XTB_SHA and runtime['parameters']['sha256'] == PARAM_SHA,
            'GFN runtime identity differs')
    verify(runtime['executable']); verify(runtime['parameters'])
    for p in runtime['runtime_libraries']:
        verify(p)
    static = read(prior['shared_static_Mex']); stcfg = read(static['config'])
    field = read(stcfg['field_result']); fcfg = read(field['config'])
    require(fcfg['solvent_result'] == cfg['GFN_origin_result'], 'Mex/LBT GFN runtime lineage differs')
    inv = read(mcfg['inventory'])
    require(inv['candidates'] == basis['candidates'] and len(inv['cells']) == 36,
            'Frozen inventory differs')
    result.update(original_comparison_rows=prior['rows'], original_source_spreads=prior['source_local_exchange_spread'],
                  baseline_replay={'status': 'passed', 'helper': prior['implementation'],
                    'arithmetic_tolerance_kcal_mol': replay_helper.TOL,
                    'R_replay_residuals_kcal_mol': [a['R_common5_Mex_kcal_mol']-b['R_common5_Mex_kcal_mol']
                                                 for a, b in zip(replay['rows'], prior['rows'])]},
                  GFN_runtime=runtime, LBT_parent=pcfg['lbt_pool'], Mex_parent=pcfg['mex_common5'],
                  old_LBT_admission=pin(d/'PREFLIGHT.json'))
    return h, cells, old, mex, basis, inv, solvent, runtime, prior


def lbt_cells(h, cells, old, inv, solvent, runtime):
    saved = {(r['candidate_id'], r['metal']): r for r in old['rows']}
    evidence = {(r['candidate_id'], r['metal']): r for r in inv['cells']}
    rows = []
    for c, admitted in cells:
        cid, metal = admitted['candidate_id'], admitted['metal']
        row = {'candidate_id': cid, 'source_id': c['source_id'], 'metal': metal, 'status': 'unavailable'}
        rows.append(row)
        try:
            oldrow = saved[(cid, metal)]; ev = evidence[(cid, metal)]
            if ev['kind'] == 'origin':
                g = next(r for r in solvent['rows'] if r['task']['source_id'] == c['source_id']
                         and r['task']['metal'] == metal and r['task']['medium'] == 'vacuum')
                gp = pin(Path(g['task']['directory'])/'RESULT.json')
                require(read(gp) == g, 'Origin GFN row receipt differs')
            else:
                obj = read(ev['evidence_record'])
                require(obj['recipe']['GFN2_runtime'] == runtime
                        and obj['recipe']['GFN2_effective_multiplicity'] == 1, 'LBT GFN recipe differs')
                gp = obj['GFN2']['vacuum']
            g = gfn_admission(gp, runtime, metal, 256, -2, ev['xyz'], 'vacuum', ev['kind'])
            require(g['task']['source_id'] == c['source_id']
                    and g['task']['physical_source'] == ev['physical_source'], 'LBT GFN source map differs')
            sol = read(oldrow['ALPB']); h.pins_in(sol['artifacts'])
            start = read(sol['artifacts']['START.json'])
            require(h.xyz(start['xyz']) == (c['symbols'], c['x'])
                    and [float(v) for v in verify(start['charges']).read_text().split()] == c['q']
                    and start['physical_source'] == ev['physical_source'], 'Matched solvent input differs')
            require(sol['total_hartree'] == oldrow['native_energies']['ALPB_Eh']
                    and sol['components_hartree'] == oldrow['native_energies']['ALPB_parts_Eh']
                    and oldrow['native_energies']['field_Eh'] == oldrow['native_energies']['MM_kcal_mol'] == 0.,
                    'LBT unchanged solvent/partition differs')
            row.update(status='available', GFN2=gp, ALPB=oldrow['ALPB'], physical_source=ev['physical_source'],
                       xyz=ev['xyz'], coordinate_sha256=ev['coordinate_sha256'],
                       native={'electronic_Eh': g['energy_hartree'], 'MM_groups_kcal_mol': {}, 'ALPB_Eh': sol['total_hartree']})
        except Exception as exc:
            row.update(failure=repr(exc), traceback=traceback.format_exc())
    return rows


def mex_cells(mex, runtime):
    rows = []
    for geometry in mex['matrix']:
        for metal in METALS:
            cell = geometry['cells'][metal]
            row = {'aliases': geometry['aliases'], 'coordinate_sha256': geometry['coordinate_sha256'],
                   'metal': metal, 'status': 'unavailable', 'old_point': cell['point']}
            rows.append(row)
            try:
                point = read(cell['point'])
                if 'cell' in point:
                    obj = read(point['cell']); gp = obj['GFN2_embedded']
                    trial = read(point['trial']); mapping_pin = trial['mapping']; xyz_pin = trial['xyz']
                    field_pin = trial['field']; alpb_pin = obj['ALPB']
                else:
                    obj = read(point['assembled_gradient']); gp = obj['GFN2_reused']['embedded']
                    task = read(gp)['task']; mapping_pin = task['mapping']; xyz_pin = obj['source_xyz']
                    field_pin = task['field']; alpb_pin = obj['ALPB']
                g = gfn_admission(gp, runtime, metal, 209, 0, xyz_pin, 'embedded', 'Mex')
                require(g['task']['mapping'] == mapping_pin and g['task']['field'] == field_pin,
                        'Mex cap/field map differs')
                mapping = read(mapping_pin); field = read(field_pin)
                require(mapping['metal'] == metal and len(mapping['core']) == 209 and len(mapping['field']) == 1840
                        and len(mapping['physical_atoms']) == 2049 and len(mapping['caps']) == 4
                        and mapping['xyz']['sha256'] == xyz_pin['sha256']
                        and mapping['pointcharges'] == g['task']['source_pointcharges']
                        and field['parameters'] == runtime['parameters']
                        and mapping['physical_atoms'] == read(point['physical']), 'Mex physical/embedding domain differs')
                sol = read(alpb_pin)
                require(sol['total_hartree'] == cell['native_energies']['ALPB_Eh']
                        and g['energy_hartree'] == cell['native_energies']['GFN2_embedded_Eh'],
                        'Mex raw components differ')
                row.update(status='available', GFN2=gp, ALPB=alpb_pin, mapping=mapping_pin, field=field_pin,
                           physical=point['physical'], native={'electronic_Eh': g['energy_hartree'],
                           'MM_groups_kcal_mol': cell['native_energies']['MM_groups_kcal_mol'],
                           'ALPB_Eh': sol['total_hartree']})
            except Exception as exc:
                row.update(failure=repr(exc), traceback=traceback.format_exc())
    return rows


def difference(a, b):
    return {'GFN2_electronic': (a['electronic_Eh']-b['electronic_Eh'])*HA,
            'MM_ledger': math.fsum(a['MM_groups_kcal_mol'][k]-b['MM_groups_kcal_mol'][k]
                                  for k in a['MM_groups_kcal_mol']),
            'fixed_FF_ALPB': (a['ALPB_Eh']-b['ALPB_Eh'])*HA}


def choose(index, ids):
    selected, ties = {}, {}
    for metal in METALS:
        origin = index[(ids[0], metal)]['native']
        values = {i: math.fsum(difference(index[(i, metal)]['native'], origin).values()) for i in ids}
        selected[metal] = min(ids, key=values.__getitem__)
        ties[metal] = [i for i in ids if values[i] == values[selected[metal]]]
    return selected, ties


def contrast(lb, mx, ls, ms):
    l = {m: lb[(ls[m], m)]['native'] for m in METALS}
    x = {m: mx[(ms[m], m)]['native'] for m in METALS}
    ld = difference(l['Dy'], l['La']); md = difference(x['Dy'], x['La'])
    parts = {'GFN2_electronic': math.fsum((x['Dy']['electronic_Eh']-x['La']['electronic_Eh'],
                                         -(l['Dy']['electronic_Eh']-l['La']['electronic_Eh'])))*HA,
             'MM_ledger': math.fsum((md['MM_ledger'], -ld['MM_ledger'])),
             'fixed_FF_ALPB': math.fsum((x['Dy']['ALPB_Eh']-x['La']['ALPB_Eh'],
                                        -(l['Dy']['ALPB_Eh']-l['La']['ALPB_Eh'])))*HA}
    value = math.fsum(parts.values())
    return {'R_kcal_mol': value, 'components_kcal_mol': parts, 'LBT_selected': ls, 'Mex_selected': ms,
            'sign': 'positive' if value > 0 else 'negative' if value < 0 else 'zero',
            'within_1_kcal_mol': abs(value) <= 1.}


def summarize(old, mex, prior, lbrows, mxrows):
    lb = {(r['candidate_id'], r['metal']): r for r in lbrows}
    mx = {(a, r['metal']): r for r in mxrows for a in r['aliases']}
    menu = mex['declared_menu']; old_ms = {m: mex['common5']['selections'][m]['geometry_alias'] for m in METALS}
    baseline = {(r['source_scope'], r['menu']): r for r in prior['rows']}
    rows = []
    for oldrow in old['comparison']['menus']:
        ids = oldrow['candidate_ids']; original = baseline[(oldrow['source_scope'], oldrow['menu'])]
        row = {'source_scope': oldrow['source_scope'], 'menu': oldrow['menu'], 'candidate_ids': ids,
               'Mex_candidate_aliases': menu, 'status': 'unavailable',
               'original_MACE_R_kcal_mol': original['R_common5_Mex_kcal_mol'],
               'original_MACE_components_kcal_mol': original['R_common5_components_kcal_mol'],
               'original_MACE_within_1_kcal_mol': original['within_1_kcal_mol_interpretive_resolution']['common5']}
        rows.append(row)
        missing = ['LBT/'+i+'/'+m for i in ids for m in METALS if lb[(i,m)]['status'] != 'available']
        missing += ['Mex/'+a+'/'+m for a in menu for m in METALS if mx[(a,m)]['status'] != 'available']
        row['missing_cells'] = missing
        if missing:
            continue
        ls, lt = choose(lb, ids); ms, mt = choose(mx, menu)
        fixed = contrast(lb, mx, oldrow['selected'], old_ms)
        reselected = contrast(lb, mx, ls, ms)
        row.update(status='available', GFN_fixed_MACE_selections=fixed, GFN_same_pool_reselection=reselected,
                   exact_ties={'LBT': lt, 'Mex': mt},
                   fixed_electronic_substitution_change_kcal_mol=fixed['R_kcal_mol']-row['original_MACE_R_kcal_mol'],
                   reselection_change_kcal_mol=reselected['R_kcal_mol']-fixed['R_kcal_mol'])
    by = {(r['source_scope'], r['menu']): r for r in rows}; spreads = {}
    for name in ('origin', 'prior', 'control', 'candidate', 'union'):
        a, b = by[('7CCO',name)], by[('7CCN',name)]
        spread = {'status':'unavailable', 'original_MACE_gap_kcal_mol':
                  prior['source_local_exchange_spread'][name]['absolute_gap_kcal_mol']['common5']}
        if a['status'] == b['status'] == 'available':
            spread.update(status='available', **{k+'_gap_kcal_mol': abs(a[k]['R_kcal_mol']-b[k]['R_kcal_mol'])
                          for k in ('GFN_fixed_MACE_selections','GFN_same_pool_reselection')})
        spreads[name] = spread
    return rows, spreads


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('config','output'):
        p.add_argument('--'+k, type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','0')) == 1,
            'One Slurm CPU required')
    cfg = json.loads(a.config.read_text()); verify(cfg['plan'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic(); cpu = time.process_time()
    result = {'status':'unavailable','protocol_id':PROTOCOL,'implementation':pin(__file__),
              'config':pin(a.config),'plan':cfg['plan'],'slurm_job_id':os.environ['SLURM_JOB_ID'],
              'declared_LBT_cells':36,'declared_Mex_cells':10,'declared_menus':15,
              'new_molecular_calls':0,'new_geometries':0,'optimizer_steps':0,'prior_comparison':cfg['prior_comparison']}
    try:
        h,cells,old,mex,basis,inv,solvent,runtime,prior = inputs(cfg,out,result)
        lb = lbt_cells(h,cells,old,inv,solvent,runtime); mx = mex_cells(mex,runtime)
        rows, spreads = summarize(old,mex,prior,lb,mx)
        result.update(status='complete' if all(r['status']=='available' for r in rows) else 'unavailable',
                      LBT_cells=lb,Mex_cells=mx,rows=rows,source_local_spreads=spreads,
                      interpretation_resolution_kcal_mol=1.,classifier=False,
                      scope='Same archived pools; f-in-core GFN electronic substitution only. No new SCC, geometry, affinity or model promotion.')
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-tick; result['process_CPU_seconds'] = time.process_time()-cpu
    write(out/'RESULT.json',result)
    lines = ['# Same-pool electronic substitution','', 'Positive R favors Mex-La plus LBT-Dy; kcal/mol. No new molecular calls.','',
             '| Scope | Menu | Status | MACE | GFN at MACE selections | GFN reselection | GFN reselection within ±1 |',
             '|---|---|---|---:|---:|---:|---|']
    for r in result.get('rows',[]):
        f=r.get('GFN_fixed_MACE_selections',{}); s=r.get('GFN_same_pool_reselection',{})
        lines.append('| '+' | '.join(str(v) for v in (r['source_scope'],r['menu'],r['status'],
                     r['original_MACE_R_kcal_mol'],f.get('R_kcal_mol','—'),s.get('R_kcal_mol','—'),s.get('within_1_kcal_mol','—')))+' |')
    lines += ['', '| Menu | Source gap: MACE | Source gap: GFN fixed selections | Source gap: GFN reselection |',
              '|---|---:|---:|---:|']
    for name, r in result.get('source_local_spreads',{}).items():
        lines.append('| '+' | '.join(str(v) for v in (name,r['original_MACE_gap_kcal_mol'],
                     r.get('GFN_fixed_MACE_selections_gap_kcal_mol','—'),
                     r.get('GFN_same_pool_reselection_gap_kcal_mol','—')))+' |')
    lines += ['', 'All 46 cells, selections, native components, old baseline replay, and source-local gaps remain in RESULT.json.',
              'Identical mechanics/fixed-charge solvent; different electronic model. Agreement or disagreement is model sensitivity, not affinity validation.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return 0 if result['status']=='complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
