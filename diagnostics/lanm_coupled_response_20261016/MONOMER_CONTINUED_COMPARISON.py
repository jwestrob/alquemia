#!/usr/bin/env python3
"""Saved Mex common5/LBT15 comparison. --config CONFIG --output FRESH; one Slurm CPU."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import time
import traceback

NAMES=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')
SCOPES=('7CCO','7CCN','global')
MENUS=('origin','prior','control','candidate','union')
EV=23.06054783061903
HA=627.509474
TOL=1e-10


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def close(a,b,label):require(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=TOL,'Saved arithmetic replay differs: '+label)


def native_from_static(row):
    a=read(row['assembled']);ml=read(a['MACE']);sol=read(a['ALPB']);mm=read(a['classical_reused']);g={k:read(p) for k,p in a['GFN2_reused'].items()}
    require(row['status']=='complete' and all(r['status']=='accepted' for r in g.values()),'Static component unavailable')
    return {'MACE_eV':ml['energy_eV'],'GFN2_embedded_Eh':g['embedded']['energy_hartree'],
        'GFN2_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mm['components_kcal_mol'],'ALPB_Eh':sol['total_hartree']}


def work(value,origin):
    return {'MACE_vacuum':(value['MACE_eV']-origin['MACE_eV'])*EV,
        'GFN2_field_difference':math.fsum((value['GFN2_embedded_Eh']-origin['GFN2_embedded_Eh'],-(value['GFN2_isolated_Eh']-origin['GFN2_isolated_Eh'])))*HA,
        'MM_ledger':math.fsum(value['MM_groups_kcal_mol'][k]-origin['MM_groups_kcal_mol'][k] for k in origin['MM_groups_kcal_mol']),
        'fixed_FF_charge_ALPB':(value['ALPB_Eh']-origin['ALPB_Eh'])*HA}


def verify_point(cell,metal,coordhash,static_row,origin):
    p=read(cell['point']);require(p['status']=='complete' and p['metal']==cell['metal']==metal
        and p['native_energies']==cell['native_energies'] and p['assembled_gradient']==cell['assembled_gradient']
        and p['physical']==cell['physical'],'Pool point binding differs')
    ph=read(p['physical']);require(len(ph)==2049 and len({a['id'] for a in ph})==2049,'Physical inventory differs')
    xyz=[v for a in ph for v in a['xyz_A']];require(len(xyz)==6147 and all(math.isfinite(v) for v in xyz),'Physical coordinate layout differs')
    require(hashlib.sha256(struct.pack('<6147d',*xyz)).hexdigest()==coordhash,'Exact physical coordinate key differs');verify(p['assembled_gradient'])
    if 'cell' not in p:
        require(p['origin']=='admitted_static_12685' and p['native_energies']==origin and p['assembled_gradient']==static_row['assembled'],'Origin cache differs')
    else:
        raw=read(p['cell']);t=read(p['trial']);require(raw['status']=='complete' and raw['metal']==metal and raw['native_energies']==p['native_energies']
            and raw['assembled_gradient']==p['assembled_gradient'] and raw['trial']==p['trial'] and t['metal']==metal and t['physical']==p['physical'],'CELL/state binding differs')
        ml=read(raw['MACE']);sol=read(raw['ALPB']);mm=read(raw['MM_result']);mmrow=next(r for r in mm['rows'] if r['task_id']==t['task_id'])
        g={k:read(raw['GFN2_'+k]) for k in ('embedded','isolated')}
        require(mmrow['status']=='complete' and mmrow['metal']==metal and all(r['status']=='accepted' and r['task']['mapping']==t['mapping'] for r in g.values()),'Saved component unavailable or map differs')
        value={'MACE_eV':ml['energy_eV'],'GFN2_embedded_Eh':g['embedded']['energy_hartree'],
            'GFN2_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mmrow['components_kcal_mol'],'ALPB_Eh':sol['total_hartree']}
        require(value==p['native_energies'],'Native component receipts differ')
    return ph


def run(cfg,result):
    prior=read(cfg['prior_comparison']);pool=read(cfg['mex_common5']);lbt=read(cfg['lbt_pool'])
    require(prior['status']=='complete' and prior['slurm_job_id']=='12722' and prior['parents']['LBT']==cfg['lbt_pool'],'12722 comparison differs')
    require(pool['status']=='complete' and pool['stage']=='pool' and pool['slurm_job_id']=='12731'
        and pool['common5']['status']=='available' and all(s['status']=='available' for s in pool['terminal_slots']),'Common5 unavailable')
    require(lbt['status']=='complete' and lbt['slurm_job_id']=='12720' and len(lbt['rows'])==36
        and all(r['status']=='complete' for r in lbt['rows']),'Matched LBT collection differs')
    pcfg=read(prior['config']);ccfg=read(pool['config']);mex=read(prior['parents']['Mex']);lc=read(lbt['config'])
    require(pcfg['mex_response']==prior['parents']['Mex']==pool['pilot_result']==ccfg['pilot_result']
        and pcfg['lbt_pool']==cfg['lbt_pool'],'Shared pilot/matched collection differs')
    require(mex['status']=='passed' and mex['slurm_job_id']=='12717' and mex['comparison']['finite_pilot_passed'],'Pilot unavailable')
    matched=read(lc['matched_origin_result']);mcfg=read(matched['config']);require(mcfg['mex_result']==mex['static_result'],'Shared static Mex parent differs')
    static=read(mex['static_result']);require(static['status']=='complete' and static['slurm_job_id']=='12685','Shared static parent unavailable')
    static_rows={r['metal']:r for r in static['rows']};require(set(static_rows)=={'La','Dy'},'Static denominator differs')
    origins={m:native_from_static(static_rows[m]) for m in ('La','Dy')}
    require(pool['original_common3']==mex['comparison']['common3'],'Archived common3 differs')
    aliases=['origin','La_pilot','Dy_pilot','La_terminal','Dy_terminal']
    require(pool['declared_menu']==pool['common5']['declared_menu']==aliases and pool['common3_replay']['declared_menu']==aliases[:3],'Pool menu differs')
    matrix=pool['matrix'];require(len(matrix)==pool['unique_physical_geometries'] and all(r['status']=='complete' for r in matrix),'Matrix unavailable')
    allaliases=[a for r in matrix for a in r['aliases']];require(len(allaliases)==5 and set(allaliases)==set(aliases),'Complete alias denominator differs')
    byalias={a:r for r in matrix for a in r['aliases']}
    for row in matrix:
        require(set(row['cells'])=={'La','Dy'},'Unpaired matrix row');physical={}
        for metal in ('La','Dy'):
            physical[metal]=verify_point(row['cells'][metal],metal,row['coordinate_sha256'],static_rows[metal],origins[metal])
            parts=work(row['cells'][metal]['native_energies'],origins[metal]);require(set(parts)==set(row['component_work_from_original_source_kcal_mol'][metal]),'Component set differs')
            for k,v in parts.items():close(v,row['component_work_from_original_source_kcal_mol'][metal][k],'matrix '+k)
            close(math.fsum(parts.values()),row['total_work_from_original_source_kcal_mol'][metal],'matrix total')
        require([a['id'] for a in physical['La']]==[a['id'] for a in physical['Dy']]
            and [a['xyz_A'] for a in physical['La']]==[a['xyz_A'] for a in physical['Dy']],'Paired physical coordinates differ')
    deltas={}
    for key,menu in (('common3_replay',aliases[:3]),('common5',aliases)):
        selected=pool[key]['selections'];selected_rows={}
        for metal in ('La','Dy'):
            s=selected[metal];row=byalias[s['geometry_alias']];winner=min(menu,key=lambda a:byalias[a]['total_work_from_original_source_kcal_mol'][metal])
            require(s['geometry_alias']==winner and s['coordinate_sha256']==row['coordinate_sha256']
                and s['all_geometry_aliases']==row['aliases'] and s['point']==row['cells'][metal]['point'],'Selection/point/minimum differs')
            close(s['work_from_original_source_kcal_mol'],row['total_work_from_original_source_kcal_mol'][metal],'selected work');selected_rows[metal]=row
        delta={k:selected_rows['Dy']['component_work_from_original_source_kcal_mol']['Dy'][k]-selected_rows['La']['component_work_from_original_source_kcal_mol']['La'][k] for k in NAMES}
        close(math.fsum(delta.values()),pool[key]['conditional_Dy_minus_La_accommodation_change_kcal_mol'],'selected delta');deltas[key]=delta
    for k in NAMES:close(deltas['common3_replay'][k],prior['Mex_accommodation_component_change_kcal_mol'][k],'12722 component')
    oldrows={(r['source_scope'],r['menu']):r for r in prior['rows']};menus=lbt['comparison']['menus'];expected={(s,m) for s in SCOPES for m in MENUS}
    require(len(prior['rows'])==len(menus)==15 and set(oldrows)==expected and {(r['source_scope'],r['menu']) for r in menus}==expected
        and all(r['status']=='available' for r in menus),'Full 15-menu denominator differs')
    rows=[]
    for r in menus:
        old=oldrows[(r['source_scope'],r['menu'])];require(old['LBT_selected']==r['selected'] and old['R_static_Mex_kcal_mol']==r['R_kcal_mol'],'Old static column/selection differs')
        close(math.fsum(r['R_components_kcal_mol'].values()),r['R_kcal_mol'],'LBT static sum')
        c3={k:math.fsum((r['R_components_kcal_mol'][k],deltas['common3_replay'][k])) for k in NAMES}
        for k in NAMES:close(c3[k],old['R_accommodated_components_kcal_mol'][k],'Old common3 component')
        close(math.fsum(c3.values()),old['R_accommodated_Mex_kcal_mol'],'Old common3 column')
        changed={k:math.fsum((r['R_components_kcal_mol'][k],deltas['common5'][k])) for k in NAMES};value=math.fsum(changed.values())
        columns={'static':old['R_static_Mex_kcal_mol'],'common3':old['R_accommodated_Mex_kcal_mol'],'common5':value}
        rows.append({'source_scope':r['source_scope'],'menu':r['menu'],'LBT_selected':r['selected'],
            'Mex_common3_selected':old['Mex_selected'],'Mex_common5_selected':pool['common5']['selections'],
            'R_static_Mex_kcal_mol':columns['static'],'R_common3_Mex_kcal_mol':columns['common3'],'R_common5_Mex_kcal_mol':value,
            'R_common3_components_kcal_mol':old['R_accommodated_components_kcal_mol'],'R_common5_components_kcal_mol':changed,
            'within_1_kcal_mol_interpretive_resolution':{k:abs(v)<=1. for k,v in columns.items()}})
    by={(r['source_scope'],r['menu']):r for r in rows};spreads={}
    for menu in MENUS:
        scopes={s:{k:by[(s,menu)]['R_'+k+'_Mex_kcal_mol'] for k in ('static','common3','common5')} for s in ('7CCO','7CCN')}
        gaps={k:abs(scopes['7CCN'][k]-scopes['7CCO'][k]) for k in ('static','common3','common5')}
        close(gaps['static'],prior['source_local_exchange_spread_kcal_mol'][menu],'Old source gap')
        spreads[menu]={'source_local_R_kcal_mol':scopes,'absolute_gap_kcal_mol':gaps}
    result.update(status='complete',rows=rows,source_local_exchange_spread=spreads,shared_static_Mex=mex['static_result'],
        Mex_common3_component_change_kcal_mol=deltas['common3_replay'],Mex_common5_component_change_kcal_mol=deltas['common5'],
        Mex_common5_accommodation_change_kcal_mol=math.fsum(deltas['common5'].values()),
        Mex_common3_to_common5_change_kcal_mol=pool['common3_to_common5_Dy_minus_La_change_kcal_mol'],
        Mex_common3_selected=pool['common3_replay']['selections'],Mex_common5_selected=pool['common5']['selections'],
        interpretive_resolution_kcal_mol=1.,interpretive_flag_is_classifier=False,arithmetic_replay_tolerance_kcal_mol=TOL,
        scope='Consumed structures, fixed occupancy/chemistry and nonstationary bounded response; no affinity or independent-protein accuracy claim')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);a=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ['SLURM_CPUS_PER_TASK'])==1,'One Slurm CPU required')
    cfg=json.loads(a.config.read_text());require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs');verify(cfg['plan'])
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'running','slurm_job_id':os.environ['SLURM_JOB_ID'],'implementation':pin(__file__),'config':pin(a.config),'plan':cfg['plan'],
        'parents':{k:cfg[k] for k in ('prior_comparison','mex_common5','lbt_pool')},'new_molecular_calls':0,'new_geometries':0}
    try:run(cfg,result)
    except Exception as exc:result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick;result['process_CPU_seconds']=time.process_time()-cpu;write(out/'RESULT.json',result)
    if result['status']=='complete':
        lines=['# Saved monomer comparison after Mex continuation','','Positive R favors Mex-La plus LBT-Dy; kcal/mol. All 15 existing LBT menus retained. No molecular calls.','',
            '| LBT scope | Menu | Static Mex | Common3 Mex | Common5 Mex | Common5 within ±1 |','|---|---|---:|---:|---:|---|']
        lines += [f"| {r['source_scope']} | {r['menu']} | {r['R_static_Mex_kcal_mol']:.6f} | {r['R_common3_Mex_kcal_mol']:.6f} | {r['R_common5_Mex_kcal_mol']:.6f} | {r['within_1_kcal_mol_interpretive_resolution']['common5']} |" for r in result['rows']]
        lines += ['','| LBT menu | Source gap, static | Source gap, common3 | Source gap, common5 |','|---|---:|---:|---:|']
        lines += [f"| {m} | {v['absolute_gap_kcal_mol']['static']:.6f} | {v['absolute_gap_kcal_mol']['common3']:.6f} | {v['absolute_gap_kcal_mol']['common5']:.6f} |" for m,v in result['source_local_exchange_spread'].items()]
        lines += ['','The ±1 kcal/mol flag is an existing interpretive resolution, not classifier calibration or an uncertainty interval. Source dependence, fixed occupancy and nonstationary bounded endpoints remain. No affinity or independent-protein accuracy claim.']
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
