#!/usr/bin/env python3
"""Summarize pinned saved transfer exports; no scoring, calibration or molecular calls."""
import argparse, collections, csv, hashlib, json, math, os, subprocess, sys, time, traceback
from pathlib import Path

POOLS=('common3','common4','common6'); POLICIES=('operational','mathematical')
PAIRS=(('common3','common4'),('common4','common6'),('common3','common6'))
IDENTITY=('context_id','source_case_id','case_id','protein','biological_group','source_conditioning_metal','expected_class')
OUTCOMES=('correct','wrong','inconclusive','unavailable','calibration_unavailable')
FLAGS=('lost_correct','new_wrong','lost_score_coverage')
REPORTER_SHA='5c526e22a182ae24ef6c0f6cfefdd506963765b4eafbcefab8fc990df0004f99'

def require(v,s):
    if not v: raise ValueError(s)
def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    require(pin(p['path'])==p,'Changed input: '+p['path']); return Path(p['path'])
def read(p): return json.loads(check(p).read_text())
def write(p,v):
    with Path(p).open('x') as f: json.dump(v,f,indent=2,allow_nan=False); f.write('\n')
def parsed(k,v):
    if v=='': return None
    if v in ('True','False'): return v=='True'
    if v.startswith(('[','{')) or v=='null': return json.loads(v)
    if k.endswith('_kcal_mol'):
        x=float(v); require(math.isfinite(x),'Nonfinite saved scalar'); return x
    if k in ('declared_memberships','available_scores'): return int(v)
    return v
def table(p):
    with check(p).open(newline='') as f:
        return [{k:parsed(k,v) for k,v in r.items()} for r in csv.DictReader(f,delimiter='\t')]
def without_alias(r): return {k:v for k,v in r.items() if k not in ('record_id','cohorts')}
def counts(rows):
    c=collections.Counter(r['outcome'] for r in rows)
    require(set(c)<=set(OUTCOMES),'Unknown outcome'); return {k:c[k] for k in OUTCOMES}

def summarize(inputs):
    p=inputs['files']; result=read(p['RESULT.json']); execution=read(p['EXECUTION.json']); source=read(p['INPUTS.json'])
    require(execution['status']=='complete' and execution['implementation']==source['implementation'] and source['implementation']['sha256']==REPORTER_SHA,'Expected completed saved reporter')
    check(source['implementation'])
    require(result['status'] in ('complete','partial') and result['declared_contexts']==257 and result['declared_memberships']==260,'Wrong declared denominator')
    require(result['source_collection']==source['collection']==execution['source_collection'] and result['source_execution']==source['execution']==execution['source_execution'],'Reporter source join differs')
    require(result['new_calibration'] is False and result['new_class_decisions'] is False and result['new_molecular_calls']==0 and result['consumed_correlated_development_evidence'],'Unexpected scientific scope')
    members=[r for r in table(p['MEMBERS.tsv']) if r['policy'] in POLICIES]
    records={}; unique={}; aliases=collections.defaultdict(set); cohorts=collections.defaultdict(set)
    for r in members:
        rid,cid=r['record_id'],r['context_id']; key=(rid,r['pool'],r['policy']); uk=(cid,r['pool'],r['policy'])
        require(r['pool'] in POOLS and key not in records,'Duplicate/unknown membership pool')
        require(r['outcome'] in OUTCOMES and ((r['R_kcal_mol'] is None)==(r['outcome']=='unavailable')),'Score availability differs from saved outcome')
        records[key]=r; payload=without_alias(r)
        require(uk not in unique or unique[uk]==payload,'Context aliases disagree; refusing deduplication')
        unique[uk]=payload; aliases[cid].add(rid); cohorts[cid].update(r['cohorts'])
    record_ids={k[0] for k in records}; context_ids=set(aliases)
    require(len(records)==260*6 and len(record_ids)==260 and len(unique)==257*6 and len(context_ids)==257,'Full257/260 denominator required')
    expected={(cid,pool,policy) for cid in context_ids for pool in POOLS for policy in POLICIES}
    require(set(unique)==expected and sum(len(v) for v in aliases.values())==260,'Missing pool or record alias collision')
    identities={}
    for cid in sorted(context_ids):
        r=unique[cid,'common3','operational']; ident={k:r[k] for k in IDENTITY}
        require(all({k:unique[cid,pool,policy][k] for k in IDENTITY}==ident for pool in POOLS for policy in POLICIES),'Context identity changed across pools')
        identities[cid]={**ident,'record_ids':sorted(aliases[cid]),'cohorts':sorted(cohorts[cid])}
    unique_counts={pool:{policy:counts([unique[cid,pool,policy] for cid in context_ids]) for policy in POLICIES} for pool in POOLS}
    for pool in POOLS:
        for policy in POLICIES:
            require(sum(unique[cid,pool,policy]['R_kcal_mol'] is not None for cid in context_ids)==result['unique_context_coverage'][pool][policy],'Unique coverage differs from reporter')
    transitions={}; seen=set()
    for r in table(p['TRANSITIONS.tsv']):
        if r['policy'] not in POLICIES: continue
        rid,cid=r['record_id'],r['context_id']; before,after=r['before_pool'],r['after_pool']; policy=r['policy']
        key=(rid,before,after,policy); uk=(cid,before,after,policy)
        require((before,after) in PAIRS and key not in seen,'Duplicate/unknown transition'); seen.add(key)
        b,a=records[rid,before,policy],records[rid,after,policy]
        require(all(r[k]==b[k]==a[k] for k in IDENTITY),'Transition identity differs')
        for side,v in (('before',b),('after',a)):
            require(r[side+'_R_kcal_mol']==v['R_kcal_mol'] and r[side+'_decision']==v['decision'] and r[side+'_outcome']==v['outcome'],'Saved transition/endpoint join differs')
            if r[side+'_selected'] is not None: require(r[side+'_selected']=={z:v[z+'_candidate'] for z in ('Ca','La')},'Selected candidates differ')
        expected_flags=(b['outcome']=='correct' and a['outcome']!='correct',a['outcome']=='wrong' and b['outcome']!='wrong',b['R_kcal_mol'] is not None and a['R_kcal_mol'] is None)
        require(tuple(r[k] for k in FLAGS)==expected_flags,'Recorded regression flag differs')
        payload=without_alias(r); require(uk not in transitions or transitions[uk]==payload,'Transition aliases disagree')
        transitions[uk]=payload
    require(len(seen)==260*6 and len(transitions)==257*6,'Missing transition denominator')
    changes={}; regression_rows=[]
    for before,after in PAIRS:
        label=before+'_to_'+after; changes[label]={}
        for policy in POLICIES:
            rows=[transitions[cid,before,after,policy] for cid in sorted(context_ids)]
            def selected(pred): return [{**identities[r['context_id']],**r} for r in rows if pred(r)]
            changes[label][policy]={
                'decision_changes':selected(lambda r:r['before_decision']!=r['after_decision'] or r['before_outcome']!=r['after_outcome']),
                'score_changes':selected(lambda r:r['before_R_kcal_mol']!=r['after_R_kcal_mol']),
                'selection_changes':selected(lambda r:r['before_selected']!=r['after_selected']),
                'newly_correct_context_ids':[r['context_id'] for r in rows if r['after_outcome']=='correct' and r['before_outcome']!='correct'],
                'regressions':selected(lambda r:any(r[k] for k in FLAGS))}
            changes[label][policy]['counts']={k:len(changes[label][policy][k]) for k in ('decision_changes','score_changes','selection_changes','newly_correct_context_ids','regressions')}
            regression_rows.extend(changes[label][policy]['regressions'])
    cal=result['recorded_calibrations']; bands={pool:{policy:cal['calibrations'][pool][policy]['bands'] for policy in POLICIES} for pool in POOLS}
    unchanged={policy:all(bands[pool][policy]==bands['common3'][policy]==cal['old_calibrations'][policy]['bands'] for pool in POOLS) for policy in POLICIES}
    require(all(unchanged.values()),'Canonical bands unexpectedly changed')
    spreads={}
    for r in table(p['SOURCE_SPREADS.tsv']):
        if r['policy'] not in POLICIES: continue
        key=(r['cohort'],r['biological_group'],r['conditioning'],r['policy']); g=spreads.setdefault(key,{})
        require(r['pool'] in POOLS and r['pool'] not in g,'Duplicate spread/pool')
        require(set(r['available_record_ids']+r['unavailable_record_ids'])<=record_ids,'Unknown spread member')
        g[r['pool']]=r
    spread_rows=[]
    for key,group in spreads.items():
        require(set(group)==set(POOLS),'Missing saved spread pool')
        row={'cohort':key[0],'biological_group':key[1],'conditioning':key[2],'policy':key[3],'saved_pools':group,'changes':{}}
        for before,after in PAIRS:
            b,a=group[before],group[after]; same=set(b['available_record_ids'])==set(a['available_record_ids'])
            row['changes'][before+'_to_'+after]={'same_available_memberships':same,'before_range_kcal_mol':b['range_kcal_mol'],'after_range_kcal_mol':a['range_kcal_mol'],
                'range_change_kcal_mol':a['range_kcal_mol']-b['range_kcal_mol'] if same and a['range_kcal_mol'] is not None and b['range_kcal_mol'] is not None else None,
                'coverage_changed':not same}
        spread_rows.append(row)
    triples={}; triple_counts={}; triple_changes={}
    for r in table(p['TRIPLES.tsv']):
        if r['policy'] not in POLICIES: continue
        key=(r['triple_id'],r['pool'],r['policy']); require(key not in triples,'Duplicate triple'); triples[key]=r
    tids={k[0] for k in triples}; require(len(tids)==100 and len(triples)==600,'Existing100 triples required per pool/policy')
    for tid in tids:
        b=triples[tid,'common3','operational']
        require(all(all(triples[tid,pool,policy][k]==b[k] for k in ('protein','members','pool_case_ids','prior_preparation_status')) for pool in POOLS for policy in POLICIES),'Triple source definition changed')
    for policy in POLICIES:
        triple_counts[policy]={pool:counts([triples[tid,pool,policy] for tid in tids]) for pool in POOLS}
        triple_changes[policy]={b+'_to_'+a:[{'triple_id':tid,'before':triples[tid,b,policy],'after':triples[tid,a,policy]} for tid in sorted(tids) if any(triples[tid,b,policy][k]!=triples[tid,a,policy][k] for k in ('available','decision','outcome','median_R_kcal_mol'))] for b,a in PAIRS}
    excluded=table(p['SUPERSEDED_PRIMARY.tsv']); require(len(excluded)==len({r['record_id'] for r in excluded})==18 and not ({r['record_id'] for r in excluded}&record_ids),'Superseded rows differ or inflate current denominator')
    final_status={policy:{status:[{**identities[cid],'saved':unique[cid,'common6',policy]} for cid in sorted(context_ids) if unique[cid,'common6',policy]['outcome']==status] for status in ('inconclusive','wrong','unavailable','calibration_unavailable')} for policy in POLICIES}
    return {'schema':'uniform_amide_saved_decision_summary.v1','scientific_status':result['status'],
        'denominators':{'current_unique_contexts':257,'current_memberships':260,'superseded_historical_records':18,'original_ledger_records':278,'existing_triples_per_pool':100},
        'alias_contexts':[identities[cid] for cid in sorted(context_ids) if len(aliases[cid])>1],
        'unique_context_counts':unique_counts,'recorded_membership_cohort_counts':result['cohorts'],'changes':changes,
        'marginal_relaxation':{'definition':'common4_to_common6; common3_to_common4 is the fixed-flip increment','policies':changes['common4_to_common6'],
            'full_pool_and_fixed_flip_newly_correct_contexts_equal':{policy:changes['common3_to_common4'][policy]['newly_correct_context_ids']==changes['common3_to_common6'][policy]['newly_correct_context_ids'] for policy in POLICIES}},
        'remaining_common6_statuses':final_status,'any_regression_or_coverage_loss':bool(regression_rows),'regressions':regression_rows,
        'fixed_canonical_bands':bands,'bands_equal_old_and_all_pools':unchanged,'calibration_evidence':cal,
        'saved_source_spreads':spread_rows,'existing_triple_counts':triple_counts,'existing_triple_changes':triple_changes,
        'existing_triples_unavailable':{policy:{pool:[triples[tid,pool,policy] for tid in sorted(tids) if not triples[tid,pool,policy]['available']] for pool in POOLS} for policy in POLICIES},
        'superseded_primary_records':excluded,'copied_actual_costs':result['costs'],'copied_COSTS_rows':table(p['COSTS.tsv']),
        'resource_receipt_links':result['resource_receipt_links'],'candidate_availability':result['candidate_availability'],'failed_new_cells':result['failed_new_cells'],
        'saved_algebra_checks':result['recorded_algebra_checks'],'new_molecular_calls':0,'scores_recomputed':False,'bands_refitted':False,
        'interpretation':'Consumed, correlated development evidence. Fixed-flip and subsequent relaxation are separate increments; score/range changes are not additional correct calls. Saved range differences are comparable only at identical available membership. No biological independence, equilibrium or production claim.'}

def main():
    a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path,required=True);a.add_argument('--inputs-sha',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'All saved-data numerical aggregation requires Slurm')
    placement=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True);fields=dict(v.split('=',1) for v in placement.split() if '=' in v)
    require(fields.get('JobState')=='RUNNING' and fields.get('NumCPUs')=='1','One running Slurm CPU required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'status':'failed','job_id':os.environ['SLURM_JOB_ID'],'placement':placement,'implementation':pin(__file__),'new_molecular_calls':0}
    try:
        ip={'path':str(args.inputs.resolve()),'sha256':args.inputs_sha};inputs=read(ip);require(pin(sys.executable)==inputs['python'],'Pinned stdlib interpreter differs')
        for value in inputs['files'].values():check(value)
        write(out/'INPUTS.json',{'declaration':ip,**inputs,'implementation':pin(__file__)})
        result=summarize(inputs);write(out/'RESULT.json',result);receipt.update(status='complete',result=pin(out/'RESULT.json'))
    except BaseException as e:receipt.update(error=str(e),traceback=traceback.format_exc())
    finally:receipt['elapsed_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':raise SystemExit(main())
