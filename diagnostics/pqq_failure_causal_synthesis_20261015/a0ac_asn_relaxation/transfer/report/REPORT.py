#!/usr/bin/env python3
"""Export recorded transfer results; descriptive aggregation only, on one Slurm CPU."""
import argparse,collections,csv,hashlib,json,math,os,subprocess,time,traceback
from pathlib import Path
POOLS=('common3','common4','common6');POLICIES=('operational','mathematical','origin');METALS=('Ca','La')
OUTCOMES=('correct','wrong','inconclusive','unavailable','calibration_unavailable')
TERMS=('MACE','frozen_solvent','total','Gelec','Ghb','Gsasa','Gshift')
TRANSITIONS=(('common3','common4'),('common4','common6'),('common3','common6'))
def require(v,s):
    if not v:raise ValueError(s)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
    if isinstance(p,dict):require(pin(p['path'])==p,'Changed artifact: '+p['path']);p=p['path']
    return json.loads(Path(p).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def table(p,rows,fields=None):
    fields=fields or list(rows[0])
    with Path(p).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()})
def score(r,pool,policy):
    v=r['pools'][pool][policy];return v['R_model_kcal_mol'] if v is not None else None
def outcome(r,pool,policy):return r['decisions_by_pool'][pool][policy]['outcome']
def identity(r):return {k:r[k] for k in ('record_id','context_id','source_case_id','case_id','protein','biological_group','source_conditioning_metal','expected_class')}
def components(value,prefix):return {prefix+k+'_kcal_mol':None if value is None else value[k] for k in TERMS}

def validate(d,e,c):
    require(d['protocol_id']=='uniform_terminal_amide_transfer_solvent_v1' and d['terminal'],'Expected terminal transfer solvent collection')
    require(e['config']==d['config'] and e['status'] in ('complete','partial','failed') and e['action']=='execute','Source execution/config mismatch')
    require(d['declared_contexts']==257 and d['declared_memberships']==260 and len(d['rows'])==260 and len(d['unique_context_rows'])==257,'Incomplete declared population')
    require(d['new_calibration'] is False and d['production_changed'] is False and d['development_evidence'] is True,'Unexpected scientific scope')
    require(d['pool_orders']==c['pool_orders'] and tuple(d['pool_orders'])==POOLS,'Pool order changed')
    ledger=read(d['historical_ledger']);cohorts={'expanded135':ledger['cohorts']['expanded135'],'current_ca125':[r['record_id'] for r in ledger['cohorts']['current_ca125_annotated_overlay']],'two_ca15':ledger['cohorts']['two_ca15']}
    lookup={r['record_id']:r for r in d['rows']};unique={r['context_id']:r for r in d['unique_context_rows']}
    require(len(lookup)==260 and len(unique)==257 and set(lookup)==set(cohorts['expanded135'])|set(cohorts['current_ca125']),'Membership identities differ')
    require({k:len(v) for k,v in cohorts.items()}=={'expanded135':135,'current_ca125':125,'two_ca15':15} and set(cohorts['two_ca15'])<=set(cohorts['current_ca125']),'Cohort identities differ')
    require(len(ledger['records'])==278 and len(d['superseded_primary_records'])==18 and {r['record_id'] for r in d['superseded_primary_records']}=={r['record_id'] for r in ledger['records']}-set(lookup),'Historical exclusions differ')
    require(d['calibration_reference']==c['canonical_gate'],'Recorded canonical calibration differs')
    declared=read(c['inventory'])['cohorts']['existing_100_triples'];require(len(declared)==100,'Exact existing100 triples required')
    for pool in POOLS:
        triples=d['existing_triples'][pool];require(len(triples)==100 and [t['triple_id'] for t in triples]==[t['triple_id'] for t in declared],'Triple membership/order differs')
        for t,old in zip(triples,declared):require(t['members']==old['members'] and t['pool_case_ids']==old['pool_case_ids'] and t['prior_preparation_status']==old['status'],'Triple source definition changed')
    for r in d['rows']:
        u=unique[r['context_id']];require(r['record_id'] in u['memberships'],'Context alias missing')
        for k in ('pools','changes','all_origin_relative_works_kcal_mol','evidence'):require(r[k]==u[k],'Membership scientific record differs from its context')
        for pool in POOLS:
            for policy in POLICIES:
                v=score(r,pool,policy);require(v is None or math.isfinite(v),'Nonfinite recorded score')
                for field in ('decisions_by_pool','old_band_transfer'):
                    dec=r[field][pool][policy];require(dec['outcome'] in OUTCOMES,'Unknown recorded outcome')
                    require((v is None)==(dec['outcome']=='unavailable'),'Missing score/decision inconsistency')
        require(all(set(r['evidence'][z])==set(d['pool_orders']['common6']) for z in METALS),'Candidate evidence omitted')
    return cohorts,lookup

def known_class_ordering(cohort,rows,cal):
    ordering=[]
    for policy in POLICIES:
        matched=[r for r in rows if all(score(r,p,policy) is not None for p in POOLS)]
        for pool in POOLS:
            for coverage,members in (('available',[r for r in rows if score(r,pool,policy) is not None]),('matched_common3_common4_common6',matched)):
                ca=[r for r in members if r['expected_class']=='Ca'];la=[r for r in members if r['expected_class']=='La'];cmax=max((score(r,pool,policy) for r in ca),default=None);lmin=min((score(r,pool,policy) for r in la),default=None)
                ordering.append({'cohort':cohort,'pool':pool,'policy':policy,'coverage':coverage,'declared_memberships':len(rows),'included_memberships':len(members),'included_unique_contexts':len({r['context_id'] for r in members}),'Ca_count':len(ca),'La_count':len(la),'other_class_count':len(members)-len(ca)-len(la),
                  'maximum_known_Ca_R_kcal_mol':cmax,'minimum_known_La_R_kcal_mol':lmin,'observed_ordering_gap_kcal_mol':lmin-cmax if ca and la else None,
                  'maximum_Ca_record_ids':[r['record_id'] for r in ca if score(r,pool,policy)==cmax],'minimum_La_record_ids':[r['record_id'] for r in la if score(r,pool,policy)==lmin],
                  'fixed_canonical_calibration_status':cal[pool][policy]['status'],'fixed_canonical_gap_kcal_mol':cal[pool][policy]['gap_kcal_mol'],'fixed_canonical_bands':cal[pool][policy]['bands'],
                  'diagnostic_only_not_calibration':True})
    return ordering

def aggregates(d,cohorts,lookup):
    counts=[];spreads=[];ordering=[];cal=d['calibration_reference']['calibrations']
    for cohort,ids in cohorts.items():
        rows=[lookup[i] for i in ids]
        for policy in POLICIES:
            matched=[r for r in rows if all(score(r,p,policy) is not None for p in POOLS)]
            for pool in POOLS:
                observed=collections.Counter(outcome(r,pool,policy) for r in rows)
                require(dict(observed)==d['cohorts'][cohort][pool][policy],'Recorded cohort count differs')
                common=collections.Counter(outcome(r,pool,policy) for r in matched)
                counts.append({'cohort':cohort,'pool':pool,'policy':policy,'declared_memberships':len(rows),'available_scores':sum(score(r,pool,policy) is not None for r in rows),**{k:observed[k] for k in OUTCOMES},'matched_common3_common4_common6_scores':len(matched),**{'matched_'+k:common[k] for k in OUTCOMES}})
                groups=collections.defaultdict(list)
                for r in rows:groups[(r['biological_group'],r['source_conditioning_metal'])].append(r)
                for (protein,conditioning),members in groups.items():
                    available=[r for r in members if score(r,pool,policy) is not None];values=[score(r,pool,policy) for r in available]
                    spreads.append({'cohort':cohort,'biological_group':protein,'protein_names':sorted({r['protein'] for r in members}),'conditioning':conditioning,'pool':pool,'policy':policy,
                      'declared_memberships':len(members),'available_scores':len(values),'minimum_R_kcal_mol':min(values) if values else None,'maximum_R_kcal_mol':max(values) if values else None,
                      'range_kcal_mol':max(values)-min(values) if values else None,'available_record_ids':[r['record_id'] for r in available],'unavailable_record_ids':[r['record_id'] for r in members if score(r,pool,policy) is None]})
        ordering.extend(known_class_ordering(cohort,rows,cal))
    combined=set(cohorts['expanded135'])|set(cohorts['current_ca125'])
    ordering.extend(known_class_ordering('combined_current260',[lookup[rid] for rid in sorted(combined)],cal))
    return counts,spreads,ordering

def exports(d,cohorts):
    memberships=[];transitions=[];cells=[];triples=[]
    for r in d['rows']:
        ident={**identity(r),'cohorts':[k for k,ids in cohorts.items() if r['record_id'] in ids]}
        for pool in POOLS:
            for policy in POLICIES:
                v=r['pools'][pool][policy];dec=r['decisions_by_pool'][pool][policy];old=r['old_band_transfer'][pool][policy]
                x={**ident,'pool':pool,'policy':policy,'pool_status':r['pools'][pool]['status'],'R_kcal_mol':score(r,pool,policy),'decision':dec['decision'],'outcome':dec['outcome'],'old_band_decision':old['decision'],'old_band_outcome':old['outcome']}
                for z in METALS:
                    q=v[z]['candidate'] if v else None;x[z+'_candidate']=q;x.update(components(r['all_origin_relative_works_kcal_mol'][z][q] if q else None,z+'_origin_relative_'))
                memberships.append(x)
        for before,after in TRANSITIONS:
            for policy in POLICIES:
                b=r['decisions_by_pool'][before][policy];a=r['decisions_by_pool'][after][policy];delta=r['changes'][before+'_to_'+after][policy]
                x={**ident,'before_pool':before,'after_pool':after,'policy':policy,'before_decision':b['decision'],'after_decision':a['decision'],'before_outcome':b['outcome'],'after_outcome':a['outcome'],
                  'old_band_before':r['old_band_transfer'][before][policy],'old_band_after':r['old_band_transfer'][after][policy],
                  'lost_correct':b['outcome']=='correct' and a['outcome']!='correct','new_wrong':a['outcome']=='wrong' and b['outcome']!='wrong',
                  'lost_score_coverage':score(r,before,policy) is not None and score(r,after,policy) is None,
                  'before_R_kcal_mol':score(r,before,policy),'after_R_kcal_mol':score(r,after,policy),'recorded_R_change_kcal_mol':delta['R_model_change_kcal_mol'] if delta else None,
                  'before_selected':delta['before'] if delta else None,'after_selected':delta['after'] if delta else None}
                x.update(components(delta['differential_change_kcal_mol'] if delta else None,'paired_change_'))
                for z in METALS:x.update(components(delta['per_metal_change_kcal_mol'][z] if delta else None,z+'_change_'))
                transitions.append(x)
    for r in d['unique_context_rows']:
        for z in METALS:
            for q in d['pool_orders']['common6']:
                v=r['evidence'][z][q];w=r['all_origin_relative_works_kcal_mol'][z][q]
                cells.append({'context_id':r['context_id'],'memberships':r['memberships'],'source_case_id':r['source_case_id'],'metal':z,'candidate':q,'status':v['status'],'reason':v.get('reason'),
                  'MACE_eV':v.get('MACE_eV'),'G_qvac_hartree':v.get('G_qvac_hartree'),'native_components_hartree':v.get('components_hartree'),
                  'component_v2_accepted':v.get('component_v2',{}).get('accepted'),'historical_absolute_point_accepted':v.get('historical_absolute_point_accepted'),
                  'native_result':v.get('native_result'),'MACE_receipt':v.get('MACE'),'xyz':v.get('xyz'),'reuse_route':v.get('exact_geometry_reuse'),**components(w,'origin_relative_')})
    for pool in POOLS:
        for t in d['existing_triples'][pool]:
            for policy in POLICIES:triples.append({'pool':pool,'policy':policy,**{k:t[k] for k in ('triple_id','protein','members','pool_case_ids','prior_preparation_status')},
              'available':t['available_by_policy'][policy],'median_R_kcal_mol':t['median_R_kcal_mol'][policy],**t['decisions'][policy]})
    return memberships,transitions,cells,triples

def costs(c,e,ep):
    phases=[('transfer_solvent',ep,e)];scoreconfig=read(c['score_config']);phases.append(('transfer_MACE',c['score_execution'],read(c['score_execution'])))
    if scoreconfig.get('scout_execution'):phases.append(('transfer_MACE_scout',scoreconfig['scout_execution'],read(scoreconfig['scout_execution'])))
    rows=[]
    for label,p,v in phases:
        fields=('status','action','elapsed_seconds','wall_seconds','parent_CPU_seconds','allocated_CPU_seconds','allocated_GPU_seconds','actual_model_forwards','native_model_requests','native_failed_requests','new_MACE_calls','new_SCC_calls','new_DFT_calls','new_optimizations','SCC_started','SCC_terminal','native_component_calls','reason','collection_error')
        rows.append({'phase':label,'receipt':p,'job_id':v.get('job_id',v.get('placement',{}).get('job_id')),**{k:v.get(k) for k in fields}})
    return rows,{'measured_child_receipts':e.get('measured_resource_receipts'),'aggregation':'No recursive child reanalysis or cumulative campaign sum. Allocation CPU-seconds are not measured child CPU time; START counts are intents; terminal/native receipts include failed calls. Reused canonical/A0AC costs and earlier preparation/preflight costs are not charged again here.'}

def run(collection,execution,output):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before saved-result aggregation')
    placement=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True);fields=dict(v.split('=',1) for v in placement.split() if '=' in v)
    require(fields.get('JobState')=='RUNNING' and fields.get('NumCPUs')=='1','One running Slurm CPU required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'status':'failed','job_id':os.environ['SLURM_JOB_ID'],'implementation':pin(__file__),'new_molecular_calls':0,'placement':placement}
    try:
        cp=pin(collection);ep=pin(execution);d=read(cp);e=read(ep);c=read(d['config']);write(out/'INPUTS.json',{'collection':cp,'execution':ep,'config':d['config'],'implementation':pin(__file__)})
        cohorts,lookup=validate(d,e,c);counts,spreads,ordering=aggregates(d,cohorts,lookup);members,transitions,cells,triples=exports(d,cohorts);cost,links=costs(c,e,ep)
        regressions=[r for r in transitions if r['lost_correct'] or r['new_wrong'] or r['lost_score_coverage']]
        for name,rows in [('COHORTS',counts),('MEMBERS',members),('TRANSITIONS',transitions),('CANDIDATES',cells),('SOURCE_SPREADS',spreads),('ORDERING_DIAGNOSTIC',ordering),('TRIPLES',triples),('COSTS',cost)]:table(out/(name+'.tsv'),rows)
        table(out/'REGRESSIONS.tsv',regressions,list(transitions[0]));exclusions=[{k:r[k] for k in ('record_id','source_case_id','family','preparation_status','preparation_reason','input_status','origin_input_status','missing_cells','membership')} for r in d['superseded_primary_records']];table(out/'SUPERSEDED_PRIMARY.tsv',exclusions)
        coverage={pool:{policy:sum(score(r,pool,policy) is not None for r in d['unique_context_rows']) for policy in POLICIES} for pool in POOLS}
        result={'status':'complete' if e['status']=='complete' and coverage['common6']['operational']==257 else 'partial','source_collection':cp,'source_execution':ep,'source_execution_status':e['status'],
          'declared_contexts':257,'declared_memberships':260,'unique_context_coverage':coverage,'cohorts':counts,'recorded_calibrations':d['calibration_reference'],'ordering_diagnostic':ordering,
          'regressions':regressions,'candidate_availability':dict(collections.Counter(v['status'] for v in cells)),'failed_new_cells':d['failed_new_cells'],
          'declared_new_unique_cells':d['declared_new_unique_cells'],'available_new_unique_cells':d['available_new_unique_cells'],'recorded_algebra_checks':d['algebra_checks'],
          'costs':cost,'resource_receipt_links':links,'new_calibration':False,'new_class_decisions':False,'new_molecular_calls':0,'consumed_correlated_development_evidence':True}
        write(out/'RESULT.json',result)
        md=['# Uniform terminal-amide transfer: saved result report','',f"Report coverage: **{result['status']}**. Source execution: **{e['status']}**. All 260 memberships / 257 contexts remain; common6 operational scores are available for {coverage['common6']['operational']}/257 contexts. Earlier complete pools remain visible when additions fail.",'',
          '| Cohort | Pool | Available scores | Correct | Wrong | Inconclusive | Unavailable | Calibration unavailable |','|---|---|---:|---:|---:|---:|---:|---:|']
        for v in counts:
            if v['policy']=='operational':md.append('| '+ ' | '.join(str(v[k]) for k in ('cohort','pool','available_scores','correct','wrong','inconclusive','unavailable','calibration_unavailable'))+' |')
        md+=['','Mathematical and origin policies, old-band decisions and all membership transitions are retained in the TSVs. REGRESSIONS.tsv includes loss of a correct call, newly wrong calls or lost score coverage; no failure is removed from its denominator. The two-Ca15 cohort is a subset of currentCa125. The two main cohorts share three contexts and include the designated canonical25; their counts must not be added as independent observations.','',
          '## Descriptive class ordering','', 'These consumed-cohort extrema are diagnostics, not fitted bands or new decisions. Positive raw gap means available known classes remain ordered; a negative gap means overlap. Matched coverage uses the same members across common3/4/6. ORDERING_DIAGNOSTIC.tsv retains all tied extremal identities, coverage and fixed canonical gaps.','',
          '| Cohort | Pool | Coverage | Known Ca / La | Observed gap, kcal/mol | Fixed canonical gap, kcal/mol |','|---|---|---|---:|---:|---:|']
        for v in ordering:
            if v['policy']=='operational':md.append(f"| {v['cohort']} | {v['pool']} | {v['coverage']} | {v['Ca_count']} / {v['La_count']} | {v['observed_ordering_gap_kcal_mol']} | {v['fixed_canonical_gap_kcal_mol']} |")
        md+=['','## Scope and costs','',f"Native stage: {d['available_new_unique_cells']}/{d['declared_new_unique_cells']} declared new unique endpoints available; all earlier MACE/preparation gaps remain in CANDIDATES.tsv separately from failed_new_cells. Reported native admission and historical absolute-v1 failures are copied separately, not reclassified. All 100 pre-existing triples per pool, including prior exclusions, remain in TRIPLES.tsv.",'',
          'COSTS.tsv copies actual execution receipts for the solvent stage and its MACE/scout producers. RESULT.json links measured child receipts. Allocation time is not measured child CPU time; cached earlier campaigns and concurrent child wall times are not summed into a new total. Missing measured forward counts remain unavailable.','',
          'Source spreads show available-only ranges with explicit missing identities. Range contraction with different coverage is not evidence of reduced fold sensitivity. These bounded four-mode proposals are not unconstrained minima or equilibrium populations. No class-band fit, electronic calculation or production change occurs here. All sources were consumed during development; correlated source folds and overlapping triples are not independent biological validation.','',
          f"Source collection: `{cp['path']}` (SHA256 `{cp['sha256']}`).",f"Source execution: `{ep['path']}` (SHA256 `{ep['sha256']}`)."]
        (out/'REPORT.md').write_text('\n'.join(md)+'\n');receipt.update(status='complete',reported_scientific_status=result['status'],source_collection=cp,source_execution=ep)
    except BaseException as error:
        receipt.update(reason=str(error),traceback=traceback.format_exc());(out/'FAILURE.md').write_text('# Reporting failed\n\n'+str(error)+'\n\nNo missing scientific values were synthesized. See EXECUTION.json.\n')
    finally:receipt['elapsed_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--collection',required=True);p.add_argument('--execution',required=True);p.add_argument('--output',required=True);a=p.parse_args();raise SystemExit(run(a.collection,a.execution,a.output))
