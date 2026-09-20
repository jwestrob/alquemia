"""Preserved DFT on the exact existing100 three-of-four La-source subsets."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
from itertools import combinations
import json
from pathlib import Path
import statistics

from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new
from accommodation_fold_dft import checked_endpoint,decision

METHODS=('DFT_baseline','MACE_native_core','MACE_native_context','MACE_composite_core','MACE_composite_context')
OUTCOMES=('correct','wrong','inconclusive','unavailable')


def method_outcome(call,expected):
    if call==expected+'-supported':return 'correct'
    return 'wrong' if call in ('Ca-supported','La-supported') else call


def protein_flags(rows,method):
    cc=Counter(r['methods'][method]['outcome'] for r in rows)
    return {'counts':{k:cc[k] for k in OUTCOMES},'denominator':len(rows),
            'all_four_correct':len(rows)==4 and cc['correct']==4,
            'any_wrong':bool(cc['wrong']),'any_inconclusive':bool(cc['inconclusive']),
            'any_unavailable':bool(cc['unavailable']),
            'worst_outcome':next(k for k in ('wrong','unavailable','inconclusive','correct') if cc[k]),
            'worst_scored_outcome':next((k for k in ('wrong','inconclusive','correct') if cc[k]),'unavailable')}


def compare(dft_collection,mace_triples,plan,output):
    d=read_json(dft_collection);m=read_json(mace_triples)
    if d['schema_version']!='accommodation_fold_DFT_collection_v1' or m['schema_version']!='PQQ_all_three_of_four_La_source_development_v1':
        raise InvalidArtifact('incompatible saved comparison schemas')
    original=read_json(verify(m['comparison']));reference=read_json(verify(d['reference']))
    master=read_json(verify(d['manifest']))
    if any(master[k]!=d[k] for k in ('preparation','sources','reference')):
        raise InvalidArtifact('DFT collection differs from its frozen run manifest')
    tasks={(t['case_id'],t['metal']):t for t in master['new_tasks']}
    if d['preparation']!=original['preparation'] or d['sources']!=original['sources']:
        raise InvalidArtifact('DFT and MACE source preparation differs')
    expected_bands={'Ca_max':reference['calibration']['U_max_Ca_kcal_mol'],'La_min':reference['calibration']['L_min_La_kcal_mol']}
    if d['bands']!=expected_bands:raise InvalidArtifact('DFT bands differ from released original reference')
    source=read_json(verify(d['sources']));declared={r['case_id']:r for r in source['cases']}
    index={r['case_id']:r for r in d['rows']}
    if len(index)!=250 or set(index)!=set(declared):raise InvalidArtifact('full250-source DFT denominator changed')
    mace_index={};groups=defaultdict(list)
    for row in m['rows']:
        selected=tuple(sorted(row['members']));key=(row['root_case_id'],selected)
        rep=row['representation']
        if len(selected)!=3 or len(set(selected))!=3 or rep not in ('core','context'):raise InvalidArtifact('invalid existing triple')
        if (key,rep) in mace_index:raise InvalidArtifact('duplicate existing MACE triple')
        mace_index[key,rep]=row
        if rep=='core':groups[row['root_case_id']].append(key)
    if len(mace_index)!=200 or len(groups)!=25:raise InvalidArtifact('original triple denominator changed')
    validated={};endpoint_checks=[]
    for root,keys in groups.items():
        la4=sorted(cid for cid,s in declared.items() if s['root_case_id']==root and
                   s['source_conditioning_metal']=='La' and not s['canonical_coordinate_match'])
        if len(la4)!=4 or {key[1] for key in keys}!=set(combinations(la4,3)):
            raise InvalidArtifact('not all four exact declared triples')
        for cid in la4:
            row=index[cid];s=declared[cid]
            if any(row[k]!=s[k] for k in ('root_case_id','biological_group','expected_class','source_conditioning_metal','canonical_coordinate_match')):
                raise InvalidArtifact('DFT source metadata differs')
            if row['R_kcal_mol'] is None:
                validated[cid]=None;continue
            eps=row['endpoints'];energies={}
            for z in ('Ca','La'):
                ep=eps[z]
                if ep['status']!='complete':raise InvalidArtifact('available DFT score lacks complete endpoint')
                if any(ep[k]!=tasks[cid,z][k] for k in ('input','xyz','manifest')):
                    raise InvalidArtifact('DFT endpoint is not this exact declared source task')
                actual=checked_endpoint(ep['input'],ep['xyz'],ep['output'],ep['receipt'],ep['manifest'],ep['task_id'])
                if actual['energy_hartree']!=ep['energy_hartree']:raise InvalidArtifact('saved DFT energy differs from output')
                energies[z]=actual['energy_hartree'];endpoint_checks.append({'case_id':cid,'metal':z,'receipt':ep['receipt']})
            r=(energies['Ca']-energies['La'])*HA_TO_KCAL
            if r!=row['R_kcal_mol'] or row['decision']!=decision(r,expected_bands):raise InvalidArtifact('DFT sign/units/decision differ')
            validated[cid]=r
    rows=[];proteins=[]
    for root,keys in sorted(groups.items()):
        current=[]
        for key in sorted(keys):
            members=list(key[1]);core=mace_index[key,'core'];context=mace_index.get((key,'context'))
            if context is None or core['expected_class']!=context['expected_class'] or core['omitted_member']!=context['omitted_member']:
                raise InvalidArtifact('MACE representations disagree on triple identity')
            if any(declared[c]['expected_class']!=core['expected_class'] for c in members):raise InvalidArtifact('MACE label differs')
            missing=[cid for cid in members if validated[cid] is None]
            value=None if missing else statistics.median(validated[c] for c in members)
            call=decision(value,expected_bands)
            rr={'root_case_id':root,'biological_group':declared[members[0]]['biological_group'],
                'expected_class':core['expected_class'],'members':members,'omitted_member':core['omitted_member'],
                'methods':{'DFT_baseline':{'median_R_kcal_mol':value,'decision':call,
                   'outcome':method_outcome(call,core['expected_class']),'missing_members':missing,
                   'required_members':3,'available_members':3-len(missing),
                   'member_R_kcal_mol':{cid:validated[cid] for cid in members},
                   'missing_statuses':{cid:{'preparation_status':index[cid]['preparation_status'],
                      'preparation_reason':index[cid]['preparation_reason'],
                      'endpoint_statuses':{z:e['status'] for z,e in index[cid]['endpoints'].items()}} for cid in missing}}}}
            for rep,old in [('core',core),('context',context)]:
                for method in ('native','composite'):
                    rr['methods'][f'MACE_{method}_{rep}']=old['methods'][method]
            current.append(rr);rows.append(rr)
        proteins.append({'root_case_id':root,'biological_group':current[0]['biological_group'],
                         'expected_class':current[0]['expected_class'],
                         'methods':{method:protein_flags(current,method) for method in METHODS}})
    if len(rows)!=100 or len(proteins)!=25:raise InvalidArtifact('DFT triple denominator changed')
    totals={}
    for method in METHODS:
        cc=Counter(r['methods'][method]['outcome'] for r in rows);pp=[p['methods'][method] for p in proteins]
        totals[method]={'triple_denominator':100,'protein_denominator':25,'triples':{k:cc[k] for k in OUTCOMES},
                       'proteins_all_four_correct':sum(p['all_four_correct'] for p in pp),
                       'proteins_any_wrong':sum(p['any_wrong'] for p in pp),
                       'proteins_any_inconclusive':sum(p['any_inconclusive'] for p in pp),
                       'proteins_any_unavailable':sum(p['any_unavailable'] for p in pp),
                       'protein_worst_outcomes':dict(Counter(p['worst_outcome'] for p in pp))}
    for rep in ('core','context'):
        for method in ('native','composite'):
            if totals[f'MACE_{method}_{rep}']!=m['counts'][rep][method]:raise InvalidArtifact('copied MACE summary changed')
    result={'schema_version':'preserved_DFT_vs_MACE_same100_triples_v1','plan':record(plan),
            'DFT_collection':record(dft_collection),'MACE_triples':record(mace_triples),
            'preparation':d['preparation'],'sources':d['sources'],'DFT_reference':d['reference'],
            'DFT_bands':expected_bands,'MACE_frozen_bands':m['frozen_bands'],
            'DFT_new_endpoint_statuses':d['new_endpoint_statuses'],'verified_DFT_endpoint_receipts':endpoint_checks,
            'rows':rows,'proteins':proteins,'counts':totals,'new_molecular_calls':0,
            'threshold_refitted':False,'production_changed':False,'all_evidence_consumed':True,
            'implementation':record(__file__),'interpretation':'identical100 correlated development triples; partial DFT execution does not establish method superiority'}
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    report(result,out);return result


def report(result,out):
    lines=['# Preserved DFT and MACE on the same100 three-fold descriptors','',
           f"DFT endpoint snapshot: {result['DFT_new_endpoint_statuses']}. Zero new molecular calls for this comparison.",'',
           '| Method | Correct | Wrong | Inconclusive | Unavailable | Proteins with all four correct |',
           '|---|---:|---:|---:|---:|---:|']
    for method,c in result['counts'].items():
        t=c['triples'];lines.append(f"| {method} | {t['correct']} | {t['wrong']} | {t['inconclusive']} | {t['unavailable']} | {c['proteins_all_four_correct']}/25 |")
    lines+=['','Every method has100 declared triples and25 proteins. Each median requires all three sources; no favorable subset or successful-only median is used. DFT uses its original released bands; all four MACE arms are copied unchanged from the completed addendum.',
            '', 'Partial DFT availability depends on execution order, so these partial counts cannot establish overall superiority. Triples overlap within consumed calibration proteins and are not independent biological observations, affinity measurements or validation of actual PLM folds.',
            '', 'All members, medians, decisions, missing statuses and per-protein worst-scored/missing-coverage flags are retained in result.json and all_triples.csv. Original fold and MACE triple reports remain unchanged.', '']
    (out/'REPORT.md').write_text('\n'.join(lines))
    with (out/'all_triples.csv').open('w',newline='') as handle:
        w=csv.writer(handle);w.writerow(['root_case_id','biological_group','expected_class','members','omitted_member','method','median_R','score_units','decision','outcome','missing_members'])
        for r in result['rows']:
            for method,e in r['methods'].items():w.writerow([r['root_case_id'],r['biological_group'],r['expected_class'],';'.join(r['members']),r['omitted_member'],method,e.get('median_R_kcal_mol',e.get('median_R_model_kcal_mol')),'kcal/mol' if method=='DFT_baseline' else 'model kcal/mol',e['decision'],e['outcome'],';'.join(e['missing_members'])])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('dft-collection','mace-triples','plan','output'):p.add_argument('--'+key,required=True)
    result=compare(**vars(p.parse_args()));print(json.dumps(result['counts']))


if __name__=='__main__':main()
