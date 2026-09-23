"""Compare actual finite-pilot pools; retain old-band transfer as development only."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
from accommodation_nonlinear import relative_components
from nikasha_pool import choose_rows, pinned
from nikasha_finite_candidates import PROTOCOL, sources


def result_row(case, prior, expected, reference):
    names=[q['id'] for q in case['candidates']]
    result=case['pool'];before=prior['pool']
    if case['status']=='prepared' and choose_rows(case['matrix'],names)!=result:
        raise InvalidArtifact('actual candidate matrix does not reproduce collected result')
    if case['prior_pool']!=before:raise InvalidArtifact('candidate lost its actual prior result')
    value=None if result['status']!='available' else result['operational']['composite_R_model_kcal_mol']
    old=before['operational']['composite_R_model_kcal_mol']
    static=prior['old_result']['R0']['composite_R_model_kcal_mol']
    bands={'released':reference['old_frozen_bands'],'adaptive':reference['variants']['operational']['bands']}
    row={'case_id':case['case_id'],'expected_class':expected,'status':result['status'],
         'reason':case.get('reason') or result.get('reason'), 'static_R':static,'adaptive_R':old,
         'candidate_R':value,'delta_R_from_adaptive':None if value is None else value-old,
         'decisions':{},'endpoint_work_from_adaptive':None,'selected':None,
         'calibration_status':'unavailable_new_protocol; old_band_transfer_only'}
    for label,b in bands.items():
        row['decisions'][label]={}
        for name,v in (('static',static),('adaptive',old),('candidate',value)):
            call=decision(v,b)
            row['decisions'][label][name]={'decision':call,'outcome':outcome(call,expected)}
    if value is not None:
        works={};selected={}
        for z in ('Ca','La'):
            selected[z]=result['rows'][z]['operational_candidate']
            oldname=before['rows'][z]['operational_candidate']
            works[z]=relative_components(case['matrix'][z][selected[z]]['components'],prior['matrix'][z][oldname]['components'])
        delta=works['Ca']['composite_kcal_mol']-works['La']['composite_kcal_mol']
        if abs(delta-(value-old))>1e-6:raise InvalidArtifact('paired work sign or unit mismatch')
        row.update(endpoint_work_from_adaptive=works,selected=selected,
                   mathematical_R=result['mathematical']['composite_R_model_kcal_mol'])
    return row


def compare(inputs,collections,output):
    config=read_json(inputs);reference=pinned(config['adaptive_reference']);branches={}
    allowed={r['case_id']:r for r in config['cases']}
    for path in collections:
        data=read_json(path);m=pinned(data['manifest']);branch=m['branch']
        if data['protocol_id']!=PROTOCOL or m['protocol_id']!=PROTOCOL or m['inputs']!=record(inputs):
            raise InvalidArtifact('not a matching finite pilot collection')
        if branch in branches:raise InvalidArtifact('duplicate branch collection')
        rows=[]
        for case in data['cases']:
            if case['case_id'] not in allowed:raise InvalidArtifact('undeclared source')
            source,prior,_,_=sources(config,case['case_id'])
            rows.append(result_row(case,prior,source['known_class'],reference))
        expected=[r['case_id'] for r in config['cases'] if branch!='local_basin_breadth' or r['basin_four']]
        if [r['case_id'] for r in rows]!=expected or data['denominator']!=len(rows):
            raise InvalidArtifact('changed pilot denominator')
        complete=[r for r in rows if r['candidate_R'] is not None]
        counts={band:{method:dict(Counter(r['decisions'][band][method]['outcome'] for r in rows))
                      for method in ('static','adaptive','candidate')} for band in ('released','adaptive')}
        matched={band:{method:dict(Counter(r['decisions'][band][method]['outcome'] for r in complete))
                       for method in ('static','adaptive','candidate')} for band in ('released','adaptive')}
        branches[branch]={'collection':record(path),'rows':rows,'denominator':len(rows),'available':len(complete),
                          'full_denominator_old_band_counts':counts,'matched_old_band_counts':matched,
                          'MACE_calls_started':data['MACE_calls_started'],'MACE_calls_complete':data['MACE_calls_complete'],
                          'required_new_GFN2_calls':data['required_new_GFN2_calls'],'GFN2_complete':data['GFN2_complete'],
                          'local_basin_note':'row minima only; conditional integration reported separately' if branch=='local_basin_breadth' else None}
    result={'inputs':record(inputs),'implementation':record(__file__),'branches':branches,
            'new_molecular_calls_in_comparison':0,'new_calibration':None,'production_changed':False,
            'interpretation':'selected consumed development pilot; old-band transfer, no new calibrated accuracy estimate'}
    write_new(output,result)
    lines=['# Nikasha parallel pilot comparison','',result['interpretation'], '']
    for branch,b in branches.items():
        lines += ['## '+branch,'',f"Available: {b['available']}/{b['denominator']}. Both band sets are transferred, not recalibrated.",'',
                  '| Source | ΔR vs adaptive | Released-band adaptive → candidate | Adaptive-band adaptive → candidate |',
                  '|---|---:|---|---|']
        for r in b['rows']:
            delta='unavailable' if r['delta_R_from_adaptive'] is None else f"{r['delta_R_from_adaptive']:+.6f}"
            calls=[r['decisions'][band]['adaptive']['decision']+' → '+r['decisions'][band]['candidate']['decision'] for band in ('released','adaptive')]
            lines.append(f"| {r['case_id']} | {delta} | {calls[0]} | {calls[1]} |")
        lines.append('')
    path=Path(output).with_suffix('.md')
    with path.open('x') as stream:stream.write('\n'.join(lines)+'\n')
    return {'output':record(output),'branches':{k:{'available':v['available'],'denominator':v['denominator']} for k,v in branches.items()}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',required=True)
    p.add_argument('--collections',nargs='+',required=True);p.add_argument('--output',required=True)
    print(json.dumps(compare(**vars(p.parse_args()))))
