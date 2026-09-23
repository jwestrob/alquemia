"""Compare both strict native initializations, preserving every384cell/32case."""
from __future__ import annotations
import argparse,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new
from strict_native_pool import PROTOCOL,BRANCHES,SETTINGS
from precision_pool_continue_run import source_rows,tid
from precision_pool_continuation import CANDIDATES,data
from nikasha_pool import choose_rows
from nikasha_pool_compare import extrema_reference
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome


def compare(fresh,cold_seed,old_reference,output,reference_output):
    collections={b:read_json(p) for b,p in [('fresh',fresh),('cold_seed',cold_seed)]}
    manifests={b:data(c['manifest']) for b,c in collections.items()}
    a,b=[manifests[k] for k in BRANCHES]
    if a['inventory']!=b['inventory'] or a['agreement']!=b['agreement']:raise InvalidArtifact('matched experiment required')
    for branch in BRANCHES:
        if collections[branch]['protocol_id']!=PROTOCOL or collections[branch]['branch']!=branch or manifests[branch]['settings']!=SETTINGS:
            raise InvalidArtifact('branch/protocol/settings differ')
    inv,source=source_rows(verify(a['inventory']));ref=read_json(old_reference)
    if ref!=data(data(inv['precision34'])['new_reference']):raise InvalidArtifact('old reference differs from precision34 authority')
    lookup={b:{tid(r):r for r in c['rows']} for b,c in collections.items()};original={tid(r):r for r in source}
    if any(len(collections[b]['rows'])!=384 or set(l)!=set(original) for b,l in lookup.items()):raise InvalidArtifact('384-cell denominator differs')
    cells=[]
    for src in source:
        rs={b:lookup[b][tid(src)] for b in BRANCHES};ok=all(r['status']=='complete' for r in rs.values())
        vals={b:r.get('energy_hartree') for b,r in rs.items()};delta=(vals['cold_seed']-vals['fresh'])*HA_TO_KCAL if ok else None
        cells.append({k:src[k] for k in ('case_id','candidate','metal','medium')}|{
            'status':'complete' if ok else 'unavailable','cold_hartree':src['source_energy_hartree'],'strict_hartree':vals,
            'cold_seed_minus_fresh_kcal_mol':delta,'initialization_agreement_pass':abs(delta)<=.1 if ok else None,
            'strict_minus_original_cold_kcal_mol':{b:(vals[b]-src['source_energy_hartree'])*HA_TO_KCAL if rs[b]['status']=='complete' else None for b in BRANCHES},
            'branches':{b:{'status':r['status'],'reason':r.get('reason'),'actual':r.get('actual'),'reused':r.get('reused',False),
                'details':r.get('details'),'observed_TolE_hartree':r.get('observed_TolE_hartree'),'initial_guess':r.get('initial_guess')} for b,r in rs.items()}})
    cell_idx={tid(r):r for r in cells};rows=[]
    for group in inv['cases']:
        cid=group['case_id'];cold=next(c for c in data(group['collection'])['cases'] if c['case_id']==cid)
        matrices={};pools={}
        for branch,l in lookup.items():
            matrix={z:{} for z in ('Ca','La')}
            for z in matrix:
                for q in CANDIDATES:
                    original_cell=cold['matrix'][z][q];rr={s:l['__'.join((cid,q,z,s))] for s in ('vacuum','alpb')};ok=all(r['status']=='complete' for r in rr.values())
                    matrix[z][q]={'status':'complete' if ok else 'unavailable','xyz':original_cell['xyz'],
                        'components':{'MACE_eV':original_cell['components']['MACE_eV'],
                            'GFN2_vacuum_hartree':rr['vacuum'].get('energy_hartree'),'GFN2_ALPB_hartree':rr['alpb'].get('energy_hartree')} if ok else None,
                        'actual_low':{s:r.get('actual') for s,r in rr.items()}}
            matrices[branch]=matrix;pools[branch]=choose_rows(matrix,list(CANDIDATES))
        oldpool=choose_rows(cold['matrix'],list(CANDIDATES));complete=all(p['status']=='available' for p in pools.values())
        cc=[c for c in cells if c['case_id']==cid];contrasts={};pooldiff={}
        for q in CANDIDATES:
            xs={(z,s):cell_idx['__'.join((cid,q,z,s))]['cold_seed_minus_fresh_kcal_mol'] for z in ('Ca','La') for s in ('vacuum','alpb')}
            delta=(xs['Ca','alpb']-xs['Ca','vacuum'])-(xs['La','alpb']-xs['La','vacuum']) if all(x is not None for x in xs.values()) else None
            contrasts[q]={'cold_seed_minus_fresh_R_kcal_mol':delta,'pass':abs(delta)<=.2 if delta is not None else None}
        for v in ('mathematical','operational'):
            delta=pools['cold_seed'][v]['composite_R_model_kcal_mol']-pools['fresh'][v]['composite_R_model_kcal_mol'] if complete else None
            pooldiff[v]={'cold_seed_minus_fresh_R_kcal_mol':delta,'pass':abs(delta)<=.2 if delta is not None else None}
        stable=complete and all(c['initialization_agreement_pass'] is True for c in cc) and all(c['pass'] is True for c in contrasts.values()) and all(c['pass'] is True for c in pooldiff.values())
        row={**group,'qualified_status':'available' if stable else 'unavailable','cold_pool':oldpool,
            'branches':{b:{'raw_status':'available' if pools[b]['status']=='available' else 'unavailable','matrix':matrices[b],'pool':pools[b],
                'qualified_R':{v:pools[b][v]['composite_R_model_kcal_mol'] if stable else None for v in ('mathematical','operational')}} for b in BRANCHES},
            'cells_agree':sum(c['initialization_agreement_pass'] is True for c in cc),'cell_denominator':12,
            'same_geometry_R_checks':contrasts,'pool_R_checks':pooldiff,
            'raw_old_reference_transfer':{v:{b:decision(p[v]['composite_R_model_kcal_mol'] if p['status']=='available' else None,ref['variants'][v]['bands'])
                for b,p in [('cold',oldpool)]+list(pools.items())} for v in ('mathematical','operational')}}
        rows.append(row)
    references={}
    for b in BRANCHES:
        variants={}
        for v in ('mathematical','operational'):
            cal=[{k:r[k] for k in ('case_id','expected_class','role')}|{'R_model_kcal_mol':r['branches'][b]['qualified_R'][v]} for r in rows if r['role']=='calibration']
            variants[v]=extrema_reference(cal,v,f'strict_native_{b}_canonical25_v1')
            for r in rows:
                call=decision(r['branches'][b]['qualified_R'][v],variants[v]['bands']) if variants[v]['bands'] else 'unavailable'
                r['branches'][b].setdefault('own_reference',{})[v]={'decision':call,'outcome':outcome(call,r['expected_class'])}
        references[b]={'protocol_id':PROTOCOL,'branch':b,'reference_id':f'strict_native_{b}_canonical25_v1','variants':variants,
            'calibration_denominator':25,'crystals_or_noncanonical_used_for_fit':False,'settings':SETTINGS,
            'inventory':a['inventory'],'collections':{b:record(p) for b,p in [('fresh',fresh),('cold_seed',cold_seed)]},
            'old_reference':record(old_reference),'frozen_UTC':datetime.now(timezone.utc).isoformat(),
            'implementation':record(__file__),'production_changed':False}
    write_new(reference_output,{'protocol_id':PROTOCOL,'branches':references});counts={}
    for name,selected in [('all32',rows),('canonical25',[r for r in rows if r['role']=='calibration']),('crystals3',[r for r in rows if r['role']=='consumed_crystal_transfer']),('folds4',[r for r in rows if r['role']=='consumed_noncanonical_development'])]:
        counts[name]={'denominator':len(selected),'qualified':sum(r['qualified_status']=='available' for r in selected),
            'raw_complete':{b:sum(r['branches'][b]['raw_status']=='available' for r in selected) for b in BRANCHES},
            'own_reference':{b:{v:dict(Counter(r['branches'][b]['own_reference'][v]['outcome'] for r in selected)) for v in ('mathematical','operational')} for b in BRANCHES},
            'raw_old_reference':{v:{b:dict(Counter(outcome(r['raw_old_reference_transfer'][v][b],r['expected_class']) for r in selected)) for b in ('cold',)+BRANCHES} for v in ('mathematical','operational')}}
    result={'protocol_id':PROTOCOL,'inventory':a['inventory'],'collections':{b:record(p) for b,p in [('fresh',fresh),('cold_seed',cold_seed)]},
        'reference':record(reference_output),'case_denominator':32,'cell_denominator':384,'rows':rows,'cells':cells,'counts':counts,
        'agreeing_cells':sum(r['initialization_agreement_pass'] is True for r in cells),'settings':SETTINGS,'implementation':record(__file__),
        'production_changed':False,'new_molecular_calls_in_analysis':0,'biological_validation_claim':False}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','cells')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fresh','cold_seed','old_reference','output','reference_output'):p.add_argument('--'+k.replace('_','-'),required=True,type=Path)
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
