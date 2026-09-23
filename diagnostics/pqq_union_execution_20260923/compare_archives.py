"""Report-only reproduction against already inspected source/score artifacts."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import HA_TO_KCAL,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL
import pqq_union_execution as execution


def geometry(a,b):
    x,y=xyz(verify(a)),xyz(verify(b));same=[v[0] for v in x]==[v[0] for v in y]
    return {'elements_and_order_identical':same,'max_coordinate_difference_A':float(np.max(np.abs(np.array([v[1:] for v in x])-np.array([v[1:] for v in y])))) if same else None,
            'bytes_identical':verify(a).read_bytes()==verify(b).read_bytes()}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--static',type=Path,required=True);ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--precision-comparison',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    st,can=read_json(a.static),read_json(a.candidate);_,req,release,ref=execution.checked(verify(can['plan']))
    arc=read_json(verify(req['archive']));previous=read_json(a.precision_comparison)
    if previous['reference']!=record(verify(read_json(verify(can['plan']))['reference'])):raise ValueError('reference mismatch')
    comparison={r['case_id']:r for r in previous['rows']}
    canonical_static={r['case_id']:r for r in read_json(verify(release['artifacts']['calibration']))['rows'] if r['representation']=='context'}
    old_sources={r['case_id']:r for p in arc['source_preparations'] for r in read_json(verify(p))['cases']}
    old_pools={};pins={p['sha256']:p for p in [*previous['collections'],*ref['collections']]}
    # Precision34 reuses may be outside the four fresh transfer collections.
    for row in previous['rows']:
        if row.get('precision_collection'):pins[row['precision_collection']['sha256']]=row['precision_collection']
    for pin in pins.values():
        for c in read_json(verify(pin))['cases']:
            if c['case_id'] in old_pools and old_pools[c['case_id']]!=c:raise ValueError('ambiguous archived pool')
            old_pools[c['case_id']]=c
    static_preparation=read_json(Path(verify(st['plan'])).parent/'source_preparation/preparation.json')
    fresh_sources={c['case_id']:c for c in static_preparation['cases']};rows=[]
    for source,x,y in zip(req['cases'],st['rows'],can['rows']):
        cid=source['case_id'];alias=source['root_case_id'] if source['canonical_coordinate_match'] else cid
        old_static=(canonical_static[alias]['composite_R_model_kcal_mol'] if source['canonical_coordinate_match'] else comparison[cid]['methods']['context_composite']['R'])
        old=old_pools[alias];oldR=old['pool']['operational']['composite_R_model_kcal_mol'] if old['pool']['status']=='available' else None
        r={'case_id':cid,'archived_pool_case_id':alias,'old_static_R':old_static,'fresh_static_R':x['R'],
           'static_R_difference':x['R']-old_static if x['R'] is not None and old_static is not None else None,
           'old_candidate_R':oldR,'fresh_candidate_R':y['R'],
           'candidate_R_difference':y['R']-oldR if y['R'] is not None and oldR is not None else None,
           'static_source_geometry':{},'candidate_cells':[],'candidate_selection_unchanged':None}
        fresh=fresh_sources[cid];prior=old_sources[cid]
        if fresh['status']=='prepared' and prior['status']=='prepared':
            r['static_source_geometry']={z:geometry(fresh['representations']['context']['endpoints'][z]['xyz'],prior['representations']['context']['endpoints'][z]['xyz']) for z in ('Ca','La')}
        if y['matrix']:
            for z in ('Ca','La'):
                for candidate in ('origin','adaptive_Ca','adaptive_La'):
                    new=y['matrix'][z].get(candidate,{});prior=old['matrix'][z].get(candidate,{})
                    cell={'metal':z,'candidate':candidate,'new_status':new.get('status','missing'),'prior_status':prior.get('status','missing')}
                    if new.get('status')=='complete' and prior.get('status')=='complete':
                        cell['geometry']=geometry(new['xyz'],prior['xyz'])
                        cell['component_differences_kcal_mol']={k:(new['components'][k]-prior['components'][k])*(EV_TO_KCAL if k=='MACE_eV' else HA_TO_KCAL) for k in ('MACE_eV','GFN2_vacuum_hartree','GFN2_ALPB_hartree')}
                    r['candidate_cells'].append(cell)
            if y['pool']['status']=='available' and old['pool']['status']=='available':
                r['candidate_selection_unchanged']=all(y['pool']['rows'][z]['operational_candidate']==old['pool']['rows'][z]['operational_candidate'] for z in ('Ca','La'))
        rows.append(r)
    result={'static_result':record(a.static),'candidate_result':record(a.candidate),'precision_comparison':record(a.precision_comparison),'archive_pools':list(pins.values()),'rows':rows,'new_molecular_calls':0,'thresholds_changed':False,
            'interpretation':'Report-only repeatability; source/candidate differences and numerical recipe differences remain explicit, not a new calibration or biological validation.'}
    write_new(a.output,result);print(json.dumps({'rows':len(rows),'max_static_R_difference':max((abs(r['static_R_difference']) for r in rows if r['static_R_difference'] is not None),default=None),'max_candidate_R_difference':max((abs(r['candidate_R_difference']) for r in rows if r['candidate_R_difference'] is not None),default=None)}))

if __name__=='__main__':main()
