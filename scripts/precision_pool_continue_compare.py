"""Report both fixed passes and calibrate only the declared25 stable references."""
from __future__ import annotations
import argparse,copy,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new
from precision_pool_continue_run import source_rows,tid,SETTINGS
from precision_pool_continuation import PROTOCOL,CANDIDATES,data
from nikasha_pool import choose_rows
from nikasha_pool_compare import extrema_reference
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome


def compare(stage1,stage2,old_reference,output,reference_output):
    a,b=read_json(stage1),read_json(stage2);ma=data(a['manifest']);mb=data(b['manifest'])
    if a['stage']!=1 or b['stage']!=2 or mb['previous']!=record(stage1) or ma['inventory']!=mb['inventory'] or ma['settings']!=SETTINGS or mb['settings']!=SETTINGS:
        raise InvalidArtifact('exact two-pass protocol/linkage required')
    inv,source=source_rows(verify(ma['inventory']));ref=read_json(old_reference)
    if ref!=read_json(verify(read_json(verify(inv['precision34']))['new_reference'])):
        raise InvalidArtifact('old reference differs from exact precision34 authority')
    lookup=[{r['task_id']:r for r in c['rows']} for c in (a,b)];source_idx={tid(r):r for r in source}
    if any(len(c['rows'])!=384 or set(l)!=set(source_idx) for c,l in zip((a,b),lookup)):
        raise InvalidArtifact('complete384-cell denominators required')
    cells=[]
    for src in source:
        rr=[l[tid(src)] for l in lookup];ok=all(r['status']=='confirmed_restart' for r in rr)
        values=[r.get('energy_hartree') for r in rr]
        delta=(values[1]-values[0])*HA_TO_KCAL if ok else None
        cells.append({k:src[k] for k in ('case_id','candidate','metal','medium')}|{
            'status':'complete' if ok else 'unavailable','cold_hartree':src['source_energy_hartree'],
            'stage1_hartree':values[0],'stage2_hartree':values[1],'pass2_minus_pass1_kcal_mol':delta,
            'pass2_minus_cold_kcal_mol':(values[1]-src['source_energy_hartree'])*HA_TO_KCAL if ok else None,
            'energy_settling_pass':abs(delta)<=.1 if ok else None,
            'actual_outputs':[r.get('actual') for r in rr],'reused':bool(src['continued_reuse']),
            'density_checks':[{k:r.get('details',{}).get(k) for k in ('energy','max_density','rms_density','cycles')} for r in rr]})
    cell_idx={tid(r):r for r in cells};rows=[]
    for group in inv['cases']:
        cid=group['case_id'];cold=next(c for c in data(group['collection'])['cases'] if c['case_id']==cid)
        matrices=[];pools=[]
        for l in lookup:
            matrix={z:{} for z in ('Ca','La')}
            for z in matrix:
                for q in CANDIDATES:
                    original=cold['matrix'][z][q];rr={s:l['__'.join((cid,q,z,s))] for s in ('vacuum','alpb')}
                    ok=all(r['status']=='confirmed_restart' for r in rr.values())
                    matrix[z][q]={'status':'complete' if ok else 'unavailable','xyz':original['xyz'],
                        'components':{'MACE_eV':original['components']['MACE_eV'],
                            'GFN2_vacuum_hartree':rr['vacuum'].get('energy_hartree'),
                            'GFN2_ALPB_hartree':rr['alpb'].get('energy_hartree')} if ok else None,
                        'actual_low':{s:r.get('actual') for s,r in rr.items()}}
            matrices.append(matrix);pools.append(choose_rows(matrix,list(CANDIDATES)))
        oldpool=choose_rows(cold['matrix'],list(CANDIDATES));complete=all(p['status']=='available' for p in pools)
        cc=[c for c in cells if c['case_id']==cid];contrasts={};pooldiff={}
        for q in CANDIDATES:
            xs={(z,s):cell_idx['__'.join((cid,q,z,s))]['pass2_minus_pass1_kcal_mol'] for z in ('Ca','La') for s in ('vacuum','alpb')}
            diff=(xs['Ca','alpb']-xs['Ca','vacuum'])-(xs['La','alpb']-xs['La','vacuum']) if all(x is not None for x in xs.values()) else None
            contrasts[q]={'pass2_minus_pass1_R_kcal_mol':diff,'pass':abs(diff)<=.2 if diff is not None else None}
        for v in ('mathematical','operational'):
            diff=pools[1][v]['composite_R_model_kcal_mol']-pools[0][v]['composite_R_model_kcal_mol'] if complete else None
            pooldiff[v]={'pass2_minus_pass1_R_kcal_mol':diff,'pass':abs(diff)<=.2 if diff is not None else None}
        stable=complete and all(c['energy_settling_pass'] for c in cc) and all(c['pass'] for c in contrasts.values()) and all(c['pass'] for c in pooldiff.values())
        row={**group,'raw_status':'available' if complete else 'unavailable','qualified_status':'available' if stable else 'unavailable',
             'cold_pool':oldpool,'stage1_pool':pools[0],'stage2_pool':pools[1],'stage1_matrix':matrices[0],'stage2_matrix':matrices[1],
             'cells_stable':sum(c['energy_settling_pass'] is True for c in cc),'cell_denominator':12,'same_geometry_R_checks':contrasts,'pool_R_checks':pooldiff,
             'qualified_R':{v:pools[1][v]['composite_R_model_kcal_mol'] if stable else None for v in ('mathematical','operational')},
             'raw_old_reference_transfer':{v:{stage:decision(pool[v]['composite_R_model_kcal_mol'] if pool['status']=='available' else None,ref['variants'][v]['bands'])
                for stage,pool in [('cold',oldpool),('stage1',pools[0]),('stage2',pools[1])]} for v in ('mathematical','operational')}}
        rows.append(row)
    variants={}
    for v in ('mathematical','operational'):
        calibration=[{k:r[k] for k in ('case_id','expected_class','role')}|{'R_model_kcal_mol':r['qualified_R'][v]} for r in rows if r['role']=='calibration']
        variants[v]=extrema_reference(calibration,v,'precision_two_native_continuations_canonical25_v1')
        for r in rows:
            call=decision(r['qualified_R'][v],variants[v]['bands']) if variants[v]['bands'] else 'unavailable'
            r.setdefault('own_reference',{})[v]={'decision':call,'outcome':outcome(call,r['expected_class'])}
    reference={'protocol_id':PROTOCOL,'reference_id':'precision_two_native_continuations_canonical25_v1','variants':variants,
        'calibration_denominator':25,'crystals_or_noncanonical_used_for_fit':False,'reported_stage':2,'settings':SETTINGS,
        'inventory':ma['inventory'],'stage1':record(stage1),'stage2':record(stage2),'old_reference':record(old_reference),
        'frozen_UTC':datetime.now(timezone.utc).isoformat(),'implementation':record(__file__),'production_changed':False}
    write_new(reference_output,reference);counts={}
    for name,selected in [('all32',rows),('canonical25',[r for r in rows if r['role']=='calibration']),('crystals3',[r for r in rows if r['role']=='consumed_crystal_transfer']),('folds4',[r for r in rows if r['role']=='consumed_noncanonical_development'])]:
        counts[name]={'denominator':len(selected),'raw_complete':sum(r['raw_status']=='available' for r in selected),
            'qualified':sum(r['qualified_status']=='available' for r in selected),
            'own_reference':{v:dict(Counter(r['own_reference'][v]['outcome'] for r in selected)) for v in variants},
            'raw_old_reference':{v:{stage:dict(Counter(outcome(r['raw_old_reference_transfer'][v][stage],r['expected_class']) for r in selected)) for stage in ('cold','stage1','stage2')} for v in variants}}
    result={'protocol_id':PROTOCOL,'inventory':ma['inventory'],'stage1':record(stage1),'stage2':record(stage2),'reference':record(reference_output),
        'reported_stage':2,'case_denominator':32,'cell_denominator':384,'rows':rows,'cells':cells,'counts':counts,
        'stable_cells':sum(r['energy_settling_pass'] is True for r in cells),'settings':SETTINGS,'implementation':record(__file__),
        'production_changed':False,'new_molecular_calls_in_analysis':0,'biological_validation_claim':False}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','cells')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('stage1','stage2','old_reference','output','reference_output'):p.add_argument('--'+k.replace('_','-'),required=True,type=Path)
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
