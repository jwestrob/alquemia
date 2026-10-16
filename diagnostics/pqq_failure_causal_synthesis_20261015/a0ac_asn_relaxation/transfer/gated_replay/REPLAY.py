#!/usr/bin/env python3
"""Retrospective label-free operational routing over saved results; Slurm only."""
import argparse, collections, hashlib, importlib.util, json, os, subprocess, sys, time, traceback
from pathlib import Path
POLICIES=('operational','mathematical'); POOLS=('common3','common4','common6')
SUMMARY_SHA='8f40e005bb5e59f1e219455af6dde4e92310b57960d9b01c4aef67089634e4a5'

def require(v,s):
    if not v:raise ValueError(s)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    require(pin(p['path'])==p,'Changed input: '+p['path']);return Path(p['path'])
def read(p):return json.loads(check(p).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(p):
    spec=importlib.util.spec_from_file_location('_saved_transfer_decision',check(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def route(available,ca_candidate,la_candidate):
    """Only screen availability and candidate IDs enter the routing decision."""
    if not available:return {'status':'screen_unavailable','trigger':None,'required_pool':None}
    require(ca_candidate is not None and la_candidate is not None,'Available screen lacks both selections')
    trigger='amide_flip' in (ca_candidate,la_candidate)
    return {'status':'triggered' if trigger else 'not_triggered','trigger':trigger,'required_pool':'common6' if trigger else 'common3'}

def endpoint(r):
    return {'available':r['R_kcal_mol'] is not None,'R_kcal_mol':r['R_kcal_mol'],
        'selected':{z:r[z+'_candidate'] for z in ('Ca','La')},'decision':r['decision'],'outcome':r['outcome'],'saved_pool_status':r['pool_status']}

def replay(c):
    h=load(c['summary_helper']);require(c['summary_helper']['sha256']==SUMMARY_SHA,'Wrong saved-summary helper')
    d=read(c['decision_result']);di=read(c['decision_inputs']);de=read(c['decision_execution'])
    require(de['status']=='complete' and de['implementation']==di['implementation']==c['summary_helper'] and de['result']==c['decision_result'],'Completed decision-summary lineage differs')
    require(read(di['declaration'])['files']==di['files'],'Decision input declaration differs')
    require(h.summarize(di)==d,'Exact completed decision-summary replay differs')
    require(pin(sys.executable)==di['python']==c['python'],'Wrong stdlib interpreter')
    unique={};aliases=collections.defaultdict(set);cohorts=collections.defaultdict(set)
    for r in h.table(di['files']['MEMBERS.tsv']):
        if r['policy'] not in POLICIES:continue
        key=(r['context_id'],r['pool'],r['policy']);payload=h.without_alias(r)
        require(key not in unique or unique[key]==payload,'Aliased scientific rows differ')
        unique[key]=payload;aliases[r['context_id']].add(r['record_id']);cohorts[r['context_id']].update(r['cohorts'])
    require(len(aliases)==257 and sum(len(v) for v in aliases.values())==260 and len(unique)==257*6,'Current context/membership denominator differs')
    rows=[]
    for cid in sorted(aliases):
        screen=unique[cid,'common4','operational']
        routing=route(screen['R_kcal_mol'] is not None,screen['Ca_candidate'],screen['La_candidate'])
        # Known class, source identity and decisions are consulted only after routing.
        row={k:screen[k] for k in h.IDENTITY};row.update(record_ids=sorted(aliases[cid]),cohorts=sorted(cohorts[cid]),
            routing=routing,screen_selected={z:screen[z+'_candidate'] for z in ('Ca','La')},comparisons={})
        for policy in POLICIES:
            full=endpoint(unique[cid,'common6',policy]);pool=routing['required_pool'];reason=None
            if pool is None:
                gated={'available':False,'R_kcal_mol':None,'selected':{'Ca':None,'La':None},'decision':None,'outcome':'unavailable','saved_pool_status':None}
                reason='required_common4_operational_screen_unavailable'
            else:
                gated=endpoint(unique[cid,pool,policy])
                if not gated['available']:reason='required_'+pool+'_'+policy+'_result_unavailable'
            available=gated['available'] and full['available']
            comparison={'required_pool':pool,'gated':gated,'full_common6':full,'unavailable_reason':reason,
                'gated_minus_full_R_kcal_mol':gated['R_kcal_mol']-full['R_kcal_mol'] if available else None,
                'score_differs':gated['R_kcal_mol']!=full['R_kcal_mol'],
                'selection_differs':gated['selected']!=full['selected'],
                'decision_or_outcome_differs':(gated['decision'],gated['outcome'])!=(full['decision'],full['outcome']),
                'coverage_lost':full['available'] and not gated['available'],
                'lost_correct':full['outcome']=='correct' and gated['outcome']!='correct',
                'new_wrong':gated['outcome']=='wrong' and full['outcome']!='wrong'}
            comparison['any_difference']=any(comparison[k] for k in ('score_differs','selection_differs','decision_or_outcome_differs'))
            row['comparisons'][policy]=comparison
        rows.append(row)
    routes=collections.Counter(r['routing']['status'] for r in rows);summaries={}
    for policy in POLICIES:
        by_method={};membership_counts={}
        for method in ('gated','full_common6'):
            by_method[method]=h.counts([r['comparisons'][policy][method] for r in rows])
            membership_counts[method]=h.counts([r['comparisons'][policy][method] for r in rows for unused in r['record_ids']])
        require(by_method['full_common6']==d['unique_context_counts']['common6'][policy],'Full-common6 count replay differs')
        summaries[policy]={'unique_context_outcomes':by_method,'membership_outcomes':membership_counts,
            'difference_context_ids':[r['context_id'] for r in rows if r['comparisons'][policy]['any_difference']],
            'untriggered_difference_context_ids':[r['context_id'] for r in rows if r['routing']['trigger'] is False and r['comparisons'][policy]['any_difference']],
            'regression_context_ids':[r['context_id'] for r in rows if any(r['comparisons'][policy][k] for k in ('lost_correct','new_wrong','coverage_lost'))],
            'gated_unavailable_context_ids':[r['context_id'] for r in rows if not r['comparisons'][policy]['gated']['available']]}
        summaries[policy]['counts']={k:len(v) for k,v in summaries[policy].items() if k.endswith('_context_ids')}
    failures=[r for r in h.table(c['candidate_table']) if r['status']!='available']
    require(all(r['context_id'] in aliases for r in failures),'Candidate failure context differs')
    return {'schema':'retrospective_operational_amide_gate.v1','post_inspection':True,'production_changed':False,
        'decision_summary':c['decision_result'],'report_inputs':di['files'],'declared_contexts':257,'declared_memberships':260,
        'rule':'common4 operational either-metal amide_flip -> common6; otherwise common3; missing screen or required result -> unavailable',
        'policy_screen':'operational for both reporting policies; no class, conditioning or desired sign in routing',
        'routing_context_counts':{k:routes[k] for k in ('triggered','not_triggered','screen_unavailable')},
        'trigger_context_ids':[r['context_id'] for r in rows if r['routing']['trigger'] is True],
        'unknown_screen_context_ids':[r['context_id'] for r in rows if r['routing']['trigger'] is None],
        'replay_search_slot_counts':{'eligible_searches':2*routes['triggered'],'avoidable_searches_in_replay':2*routes['not_triggered'],
            'unknown_screen_searches':2*routes['screen_unavailable'],'basis':'two own-metal searches per context; all cross-scoring required on triggered branch; counts only, not measured savings'},
        'summaries':summaries,'rows':rows,'alias_contexts':d['alias_contexts'],
        'fixed_canonical_bands':d['fixed_canonical_bands'],'bands_equal_old_and_all_pools':d['bands_equal_old_and_all_pools'],
        'superseded_primary_records':d['superseded_primary_records'],'saved_candidate_failures':failures,'saved_failed_new_cells':d['failed_new_cells'],
        'new_molecular_calls':0,'scores_recomputed':False,'bands_refitted':False,'measured_cost_savings':None,
        'limits':'Retrospective consumed and correlated evidence. An initially unfavorable fixed flip may relax usefully: this replay cannot establish unseen-basin recall. No automatic promotion or runtime-saving claim. Missing required data never triggers a baseline fallback.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--inputs-sha',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Saved-result aggregation requires Slurm')
    placement=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True);fields=dict(v.split('=',1) for v in placement.split() if '=' in v)
    require(fields.get('JobState')=='RUNNING' and fields.get('NumCPUs')=='1','One running Slurm CPU required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'status':'failed','job_id':os.environ['SLURM_JOB_ID'],'placement':placement,'implementation':pin(__file__),'new_molecular_calls':0}
    try:
        ip={'path':str(a.inputs.resolve()),'sha256':a.inputs_sha};c=read(ip)
        require(pin(__file__)==c['implementation'],'Replay source changed');check(c['plan']);check(c['wrapper'])
        write(out/'INPUTS.json',{'declaration':ip,**c});result=replay(c);write(out/'RESULT.json',result);receipt.update(status='complete',result=pin(out/'RESULT.json'))
    except BaseException as e:receipt.update(error=str(e),traceback=traceback.format_exc())
    finally:receipt['elapsed_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':raise SystemExit(main())
