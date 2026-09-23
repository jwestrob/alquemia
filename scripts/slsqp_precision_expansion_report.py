"""Report fixed34 precision comparisons and distinct canonical25 calibration."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
from nikasha_pool_compare import extrema_reference
import nikasha_pool
import slsqp_precision_expansion as expansion
import slsqp_precision as pilot


def endpoint(old_pin,new_pin):
    old=read_json(verify(old_pin));new=read_json(verify(new_pin)) if new_pin else None
    if new and (new['task_id']!=old['task_id'] or new['metal']!=old['metal']):raise InvalidArtifact('different endpoint identities')
    oldok=old['status']=='proposal_available';newok=bool(new and new['status']=='proposal_available')
    r={'old_receipt':old_pin,'new_receipt':new_pin,'old_status':old['status'],'new_status':new['status'] if new else 'not_run',
       'successful_old_comparison':oldok and newok,'native_energy_pass':None,'native_proposal_delta_kcal_mol':None,
       'maximum_coordinate_delta_A':None,'failed_old_last_iterate_diagnostic':None,
       'old_optimizer':old['optimizer'],'new_optimizer':new.get('optimizer') if new else None,
       'old_wall_seconds':old['wall_seconds'],'new_wall_seconds':new.get('wall_seconds') if new else None,
       'old_boundary':old.get('boundary_flag'),'new_boundary':new.get('boundary_flag') if new else None,
       'old_geometry':old.get('final_geometry'),'new_geometry':new.get('final_geometry') if new else None}
    if newok and old.get('proposal'):
        a=np.array([x[1:] for x in xyz(verify(old['proposal']['coordinate']))]);b=np.array([x[1:] for x in xyz(verify(new['proposal']['coordinate']))])
        d={'native_proposal_delta_kcal_mol':(new['proposal']['MACE_eV']-old['proposal']['MACE_eV'])*EV_TO_KCAL,
           'maximum_coordinate_delta_A':float(np.linalg.norm(b-a,axis=1).max()),
           'full_q_delta':(np.asarray(new['proposal']['full_q'])-old['proposal']['full_q']).tolist()}
        if oldok:r.update(d,native_energy_pass=abs(d['native_proposal_delta_kcal_mol'])<=pilot.GATES['proposal_energy_kcal_mol'])
        else:r['failed_old_last_iterate_diagnostic']=d
    return r


def compare(collection,output,reference_output):
    fresh=read_json(collection);pm=read_json(verify(fresh['manifest']));expansion.validate_pool(verify(fresh['manifest']))
    m=read_json(verify(pm['source_manifest']));oldref=read_json(verify(m['reference']));reuse=read_json(verify(m['reuse']));rpm=read_json(verify(reuse['manifest']))
    if fresh['denominator']!=30 or reuse['denominator']!=4:raise InvalidArtifact('all34 denominator differs')
    for k in ('protocol_id','settings','model','software','orca','GFN2_maxiter'):
        if pm[k]!=rpm[k]:raise InvalidArtifact('reused precise pool method differs')
    rows=[]
    for source in m['all_sources']:
        cid=source['case_id'];reused=source['reused_precision_pilot'];col=reuse if reused else fresh;parent=rpm if reused else pm
        c=next(x for x in col['cases'] if x['case_id']==cid);old=next(x for x in read_json(verify(source['old_pool']))['cases'] if x['case_id']==cid) if source['old_pool'] else None
        if c['source']!=next(x for x in read_json(verify(source['parent']))['cases'] if x['case_id']==cid):raise InvalidArtifact('source row changed')
        if c['status']=='prepared' and nikasha_pool.choose_rows(c['matrix'],[x['id'] for x in c['candidates']])!=c['pool']:raise InvalidArtifact('pool selection differs')
        endpoints={};cells=[];both=bool(old and old['pool']['status']=='available' and c['pool']['status']=='available')
        for z in ('Ca','La'):
            npth=verify(parent['source_manifest']).parent/'proposals'/(cid+'__'+z)/'result.json'
            endpoints[z]=endpoint(source['old_receipts'][z],record(npth) if npth.exists() else None)
            if both:
                for name in ('adaptive_Ca','adaptive_La'):
                    oc=old['matrix'][z][old['aliases'][name]['representative']]['components'];nc=c['matrix'][z][c['aliases'][name]['representative']]['components']
                    dv=(nc['GFN2_vacuum_hartree']-oc['GFN2_vacuum_hartree'])*HA_TO_KCAL;da=(nc['GFN2_ALPB_hartree']-oc['GFN2_ALPB_hartree'])*HA_TO_KCAL
                    cells.append({'metal':z,'candidate':name,'MACE_delta_kcal_mol':(nc['MACE_eV']-oc['MACE_eV'])*EV_TO_KCAL,
                      'vacuum_delta_kcal_mol':dv,'ALPB_delta_kcal_mol':da,'transfer_delta_kcal_mol':da-dv,
                      'vacuum_pass':abs(dv)<=pilot.GATES['GFN2_cell_work_kcal_mol'],'ALPB_pass':abs(da)<=pilot.GATES['GFN2_cell_work_kcal_mol']})
        scores={v:c['pool'][v]['composite_R_model_kcal_mol'] if c['pool']['status']=='available' else None for v in ('mathematical','operational')}
        oldscores={v:old['pool'][v]['composite_R_model_kcal_mol'] if old and old['pool']['status']=='available' else None for v in scores}
        delta={v:scores[v]-oldscores[v] if scores[v] is not None and oldscores[v] is not None else None for v in scores}
        passes=bool(both and all(e['native_energy_pass'] for e in endpoints.values()) and len(cells)==4 and all(x['vacuum_pass'] and x['ALPB_pass'] for x in cells) and abs(delta['operational'])<=pilot.GATES['pooled_R_kcal_mol']) if old is not None else None
        calls={v:{'old':decision(oldscores[v],oldref['variants'][v]['bands']),'new':decision(scores[v],oldref['variants'][v]['bands'])} for v in scores}
        rows.append({**source,'expected_class':source['known_class'],'status':c['pool']['status'],'scores':scores,'old_scores':oldscores,'delta_R':delta,
          'endpoints':endpoints,'components':cells,'numerical_pass':passes,'old_reference_transfer':calls,
          'new_pool':c['pool'],'old_pool_result':old['pool'] if old else None,
          'selected_candidates':{v:{z:c['pool']['rows'][z][v+'_candidate'] for z in ('Ca','La')} for v in scores} if c['pool']['status']=='available' else None,
          'collection':m['reuse'] if reused else record(collection)})
    if len(rows)!=34 or len({r['case_id'] for r in rows})!=34:raise InvalidArtifact('fixed denominator differs')
    variants={}
    for v in ('mathematical','operational'):
        chosen=[{**r,'R_model_kcal_mol':r['scores'][v]} for r in rows if r['role']=='calibration']
        variants[v]=extrema_reference(chosen,v,'union_ftol1e8_canonical25_v1')
        for r in rows:
            call=decision(r['scores'][v],variants[v]['bands']) if variants[v]['bands'] else 'unavailable'
            r.setdefault('own_reference',{})[v]={'decision':call,'outcome':outcome(call,r['expected_class'])}
    ref={'protocol_id':pilot.POOL_PROTOCOL,'reference_id':'union_ftol1e8_canonical25_v1','variants':variants,'calibration_denominator':25,
      'crystals_or_noncanonical_used_for_fit':False,'frozen_UTC':datetime.now(timezone.utc).isoformat(),'source_manifest':pm['source_manifest'],
      'collections':[record(collection),m['reuse']],'old_reference':m['reference'],'model':m['model'],'settings':pm['settings'],
      'optimizer_settings':m['settings'],'implementation':record(__file__),'production_changed':False}
    write_new(reference_output,ref)
    counts={}
    for group,selected in [('all34',rows),('canonical25',[r for r in rows if r['role']=='calibration']),('crystals3',[r for r in rows if r['case_id'] in ('1H4I','4MAE','1KB0')]),('noncanonical4',[r for r in rows if r['case_id'] in expansion.NONCANONICAL]),('old_failures2',[r for r in rows if r['case_id'] in expansion.FAILED])]:
        counts[group]={'denominator':len(selected),'available':sum(r['status']=='available' for r in selected),
          'new_own_reference':{v:dict(Counter(r['own_reference'][v]['outcome'] for r in selected)) for v in variants},
          'new_old_reference':{v:dict(Counter(outcome(r['old_reference_transfer'][v]['new'],r['expected_class']) for r in selected)) for v in variants}}
    result={'protocol_id':pilot.PROTOCOL,'collection':record(collection),'manifest':pm['source_manifest'],'reused':m['reuse'],'old_reference':m['reference'],
      'new_reference':record(reference_output),'gates':pilot.GATES,'denominator':34,'rows':rows,'counts':counts,
      'numerical_comparison_denominator':sum(r['numerical_pass'] is not None for r in rows),'numerical_passes':sum(r['numerical_pass'] is True for r in rows),
      'original225_unchanged':True,'new_calls_in_analysis':0,'production_changed':False,'implementation':record(__file__)}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('collection','output','reference_output'):p.add_argument('--'+n.replace('_','-'),required=True,type=Path)
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
