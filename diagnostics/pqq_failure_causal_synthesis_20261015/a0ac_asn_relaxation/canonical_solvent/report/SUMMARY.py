#!/usr/bin/env python3
"""Saved scalar reporting only; all aggregation and real-fixture checks on Slurm."""
import argparse, collections, csv, hashlib, importlib.util, json, math, os, shutil, subprocess, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_solvent/report'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
POOLS=('common3','common4','common6');POLICIES=('operational','mathematical','origin');METALS=('Ca','La')
TERMS={'Gelec','Ghb','Gsasa','Gshift'}
def require(v,why):
    if not v:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def freeze(source_config,output):
    sp=pin(source_config);sc=read(sp);base=Path(source_config).resolve().parent
    ep=pin(base/'execution/EXECUTION.json');ex=read(ep)
    require(ex['config']==sp and ex['status'] in ('complete','partial','failed'),'Actual terminal source execution required')
    final=base/'execution/FINAL_COLLECTION.json';fp=pin(final) if final.exists() else None
    if fp:require(read(fp)['config']==sp,'Terminal collection/config mismatch')
    summary_base=WORK/'asn_relaxation_solvent_summary_v1';summary=pin(summary_base/'execution/RESULT.json')
    require(summary['sha256']=='2d5ab6278f90b0705e96727399a5fc60ff440d2c6077ef34213f6492c7d31649','Wrong actual5030 summary')
    oldc=read(WORK/'asn_relaxation_solvent_v1/CONFIG.json');fixture=WORK/'asn_canonical_solvent_v1/preflight_5135'
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir();impl={}
    for name in ('SUMMARY.py','RUN.sbatch','PLAN.md'):
        shutil.copyfile(HERE/name,out/'implementation'/name);impl[name]=pin(out/'implementation'/name)
    c={'implementation':impl,'source_config':sp,'source_execution':ep,'source_collection':fp,'source_collection_path':str(final),
       'old_canonical':sc['canonical_collection'],'canonical_search':sc['score_collection'],'canonical_search_execution':sc['score_execution'],
       'scorer':sc['scorer'],'A0AC_summary':summary,'A0AC_summary_execution':pin(summary_base/'execution/EXECUTION.json'),
       'A0AC_collection':pin(WORK/'asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json'),
       'A0AC_execution':pin(WORK/'asn_relaxation_solvent_v1/execution/EXECUTION.json'),
       'A0AC_search':oldc['score_collection'],'A0AC_search_execution':oldc['score_execution'],
       'contact':pin(WORK/'asn_relaxation_contact_review_v1/execution/RESULT.json'),
       'failed_fixture':{'config':pin(fixture.parent/'CONFIG.json'),'execution':pin(fixture/'EXECUTION.json'),'collection_path':str(fixture/'FINAL_COLLECTION.json')},
       'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11'),
       'accounting_job_ids':['5135','5157','5158','5089','5057','4971','5030','4831','4844'],
       'output':str(out/'execution'),'new_molecular_calls':0,'fit_calibration':False}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def calibration(collection,pool,policy,execution_complete):
    if not collection or not execution_complete:return None
    v=collection['calibrations'][pool][policy]
    if v['status']!='available':require(v['bands'] is None,'Unavailable calibration has bands');return None
    require(v['available']==25 and v['bands'] is not None,'Available calibration lacks exact25')
    return v['bands']

def canonical_rows(sc,collection,old,search,selector,execution_complete):
    sources=sc['rows'];require(len(sources)==25 and len({r['record_id'] for r in sources})==25,'Exact25 declaration required')
    originals={r['record_id']:r for r in old['rows']};actual={r['record_id']:r for r in collection['rows']} if collection else {}
    require([r['record_id'] for r in sources]==[r['record_id'] for r in old['rows']],'Designated canonical order differs')
    if collection:
        require(len(actual)==len(collection['rows']) and set(actual)<=set(originals),'Unexpected/duplicate canonical rows')
        require(collection['old_calibrations']==old['calibrations'] and collection['calibrations']['common3']==old['calibrations'],'Released common3 calibration differs')
    sr={r['source_case_id']:r for r in search['rows']};rows=[]
    for source in sources:
        rid=source['record_id'];saved=actual.get(rid);previous=originals[rid]
        if saved:require(all(saved[k]==source[k] for k in source),'Saved canonical identity differs')
        item={**source,'collection_row_available':saved is not None,'pools':{},'changes':saved['changes'] if saved else None,
              'all_origin_relative_works_kcal_mol':saved['all_origin_relative_works_kcal_mol'] if saved else None,
              'search_boundary':boundaries(sr.get(source['source_case_id']))}
        for pool in POOLS:
            item['pools'][pool]={}
            for policy in POLICIES:
                v=saved['pools'][pool][policy] if saved else previous['pool'][policy] if pool=='common3' else None
                score=v['R_model_kcal_mol'] if v else None
                require(score is None or math.isfinite(score),'Nonfinite recorded score')
                oldcall=selector.decision(score,old['calibrations'][policy]['bands'],source['expected_class'])
                newcall=selector.decision(score,calibration(collection,pool,policy,execution_complete),source['expected_class'])
                if saved:
                    require(oldcall==saved['old_band_transfer'][pool][policy],'Saved old-band decision differs')
                    if execution_complete:require(newcall==saved['decisions_by_pool'][pool][policy],'Saved designated25 decision differs')
                item['pools'][pool][policy]={'R_kcal_mol':score,'selected':{z:v[z]['candidate'] for z in METALS} if v else None,
                    'old_band_sensitivity':oldcall,'designated25_calibration':newcall,'recorded_pool':v}
        rows.append(item)
    return rows

def boundaries(row):
    if row is None:return {'status':'unavailable','reason':'No saved search row'}
    result={'source_status':row.get('status'),'source_reason':row.get('reason'),'proposals':{}}
    for z in METALS:
        p=row.get('proposals',{}).get(z)
        if not p:result['proposals'][z]={'status':'unavailable','reason':'Missing actual proposal'};continue
        v=read(p);g=v.get('final_geometry',{})
        result['proposals'][z]={'receipt':p,'status':v.get('status'),'boundary_flag':v.get('boundary_flag'),
            'angular_boundary_flag':g.get('angular_boundary_flag'),'boundary_heavy_atoms':g.get('boundary_heavy_atoms'),
            'physical_feasible':g.get('physical_feasible'),'optimizer_success':v.get('optimizer',{}).get('success'),
            'unconstrained_minimum_claimed':v.get('unconstrained_minimum_claimed')}
    return result

def admission(collection):
    rows=[]
    for source in collection.get('rows',[]) if collection else []:
        for z,evidence in source['evidence'].items():
            for candidate,v in evidence.items():
                item={'source_case_id':source['source_case_id'],'metal':z,'candidate':candidate,'status':v['status'],'reason':v.get('reason')}
                if v['status']=='available':
                    n=read(v['native_result']);require(n['status']=='available' and n['G_qvac_hartree']==v['G_qvac_hartree'],'Native saved scalar mismatch')
                    require(v['component_v2']['accepted'],'Accepted cell lacks component-v2 gate')
                    require({a['component'] for a in n['checks']}==TERMS and all(a['passed'] and a['informative'] for a in n['checks']),'Native water diagonal failed')
                    require({a['component'] for a in n['vacuum_precision']}==TERMS and all(a['informative'] for a in n['vacuum_precision']),'Vacuum precision failed')
                    require(all(n['water'][k]['sha256']==n['vacuum'][k]['sha256'] for k in ('cavity','born_matrix_with_HB')),'Native cavity changed')
                    item.update(native_result=v['native_result'],component_v2_accepted=True,diagonal_passed=True,precision_informative=True,
                        historical_absolute_point_accepted=v['historical_absolute_point_accepted'])
                rows.append(item)
    return {'declared_cells':300,'available_cells':sum(x['status']=='available' for x in rows),'recorded_cells':len(rows),'cells':rows}

def overlap(collection,a0ac):
    sid='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3'
    found=[r for r in collection['rows'] if r['source_case_id']==sid] if collection else []
    result={'source_case_id':sid,'status':'unavailable','reason':'Canonical overlap row missing or incomplete'}
    if not found:return result
    require(len(found)==1,'Duplicate canonical overlap');left=found[0];right,=[r for r in a0ac['rows'] if r['source_case_id']==sid]
    order=collection['pool_orders']['common6'];require(order==a0ac['pool_orders']['common6'],'Overlap pool order differs')
    if any(left['evidence'][z].get(q,{}).get('status')!='available' or right['evidence'][z].get(q,{}).get('status')!='available' for z in METALS for q in order):return result
    cells=[]
    for z in METALS:
        for q in order:
            a,b=left['evidence'][z][q],right['evidence'][z][q];na,nb=read(a['native_result']),read(b['native_result'])
            require(a['xyz']['sha256']==b['xyz']['sha256'] and a['native_result']==b['native_result'],'Overlap exact geometry/native reuse differs')
            require(na['manifest']==nb['manifest'] and all(na['input'][k]==nb['input'][k] for k in ('maps','state_signature','source_preparation','charge','multiplicity','atom_count','valence_electrons')),'Overlap physical state/maps differ')
            cells.append({'metal':z,'candidate':q,'native_result':a['native_result'],'canonical_MACE':a['MACE'],'A0AC_MACE':b['MACE'],'xyz_sha256':a['xyz']['sha256']})
    tol=collection['algebra_checks']['tolerance_kcal_mol'];require(tol==a0ac['saved_algebra_checks']['tolerance_kcal_mol'],'Overlap frozen scalar tolerance differs')
    residuals=[]
    for pool in POOLS:
        for p in POLICIES:
            a,b=left['pools'][pool][p],right['pools'][pool][p]
            if a is None or b is None:return result
            delta=a['R_model_kcal_mol']-b['R_model_kcal_mol'];require(abs(delta)<=tol and all(a[z]['candidate']==b[z]['candidate'] for z in METALS),'Overlap pool score/selection differs')
            residuals.append({'pool':pool,'policy':p,'residual_kcal_mol':delta})
    return {**result,'status':'passed','reason':None,'cells':cells,'score_residuals':residuals,'tolerance_kcal_mol':tol,'exact_state_maps_source_and_native_results':True}

def table(path,rows):
    with Path(path).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def outcome_counts(values):
    found=collections.Counter(values)
    return {k:found[k] for k in ('correct','wrong','inconclusive','unavailable','calibration_unavailable')}
def counts(rows,field):return outcome_counts(v[field]['outcome'] for v in rows)
def receipt(p):
    v=read(p)
    keys=('status','job_id','action','elapsed_seconds','wall_seconds','allocated_CPU_seconds','SCC_started','SCC_terminal','native_component_calls','native_failed_requests','new_MACE_calls','new_DFT_calls','new_optimizations','reason','collection_error')
    return {'receipt':p,**{k:v[k] for k in keys if k in v}}
def costs(c):
    cmd=['sacct','-n','-P','-j',','.join(c['accounting_job_ids']),'--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,CPUTimeRAW,AllocTRES']
    try:
        raw=subprocess.check_output(cmd,text=True,stderr=subprocess.STDOUT,timeout=30)
        records=[dict(zip(('job_id','state','exit_code','elapsed_seconds','allocated_CPUs','allocated_CPU_seconds','allocated_TRES'),line.rstrip('|').split('|'))) for line in raw.splitlines() if line.split('|')[0] in c['accounting_job_ids']]
        return {'status':'available' if len(records)==len(c['accounting_job_ids']) else 'partial','queried_epoch':time.time(),'command':cmd,'rows':records,'raw':raw,'aggregation':'Separate job receipts; no sum across reused campaigns'}
    except Exception as e:return {'status':'unavailable','reason':str(e),'command':cmd}

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Real Slurm required before aggregation')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['implementation']['SUMMARY.py']==pin(__file__) and c['python']==pin(sys.executable),'Implementation/interpreter differs')
    placement=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True)
    fields=dict(x.split('=',1) for x in placement.split() if '=' in x);require((fields.get('JobState'),fields.get('NumNodes'),fields.get('NumCPUs'))==('RUNNING','1','1'),'One running CPU/node required')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();ex={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'placement':placement,'new_molecular_calls':0}
    try:
        for p in c['implementation'].values():check(p)
        sc=read(c['source_config']);se=read(c['source_execution']);old=read(c['old_canonical']);search=read(c['canonical_search'])
        require(se['config']==c['source_config'],'Source execution identity differs')
        final=read(c['source_collection']) if c['source_collection'] else None
        if final:
            require(final['config']==c['source_config'] and final['declared_sources']==25 and final['pool_orders']==sc['pool_orders'],'Source collection identity/scope differs')
        else:require(not Path(c['source_collection_path']).exists(),'Missing source collection appeared after freeze')
        spec=importlib.util.spec_from_file_location('_canonical_report_selector',check(c['scorer']));selector=importlib.util.module_from_spec(spec);spec.loader.exec_module(selector)
        require(sc['scorer']==old['scorer']==c['scorer'],'Saved decision policy differs')
        completed=se['status']=='complete' and final is not None
        rows=canonical_rows(sc,final,old,search,selector,completed)
        fx=c['failed_fixture'];fe=read(fx['execution']);fc=read(fx['config'])
        require(fe['config']==fx['config'] and fe['status']=='failed' and all(fe[k]==0 for k in ('SCC_started','SCC_terminal','native_component_calls')),'Real5135 fixture differs')
        require(not Path(fx['collection_path']).exists(),'Failed5135 fixture unexpectedly has collection')
        missing=canonical_rows(fc,None,old,search,selector,False)
        require(len(missing)==25 and all(r['pools']['common6'][p]['R_kcal_mol'] is None and r['pools']['common3'][p]['R_kcal_mol'] is not None for r in missing for p in POLICIES),'Actual failed-fixture denominator handling failed')
        prior=read(c['A0AC_summary']);ac=read(c['A0AC_collection']);ars=read(c['A0AC_search']);contact=read(c['contact'])
        require(prior['status']=='complete' and len(prior['rows'])==len(ac['rows'])==10 and read(c['A0AC_summary_execution'])['status']=='complete','Actual5030 complete fixture unavailable')
        require(prior['old_bands']==old['calibrations'],'A0AC old calibration differs')
        a_by={r['source_case_id']:r for r in ac['rows']};as_by={r['source_case_id']:r for r in ars['rows']};a_rows=[];flat=[]
        for r in prior['rows']:
            sid=r['source_case_id'];actual=a_by[sid];v={**r,'search_boundary':boundaries(as_by.get(sid)),'fixed_common6_calibration':{}}
            for pool in POOLS:
                for p in POLICIES:
                    s=r['pools'][pool][p];require(s['R_kcal_mol']==actual['pools'][pool][p]['R_model_kcal_mol'] and s['selected']=={z:actual['pools'][pool][p][z]['candidate'] for z in METALS},'Actual5030 score/selection replay differs')
                    require(selector.decision(s['R_kcal_mol'],old['calibrations'][p]['bands'],r['expected_class'])==s['old_band_sensitivity'],'Actual5030 old decision differs')
                    new=selector.decision(s['R_kcal_mol'],calibration(final,'common6',p,completed),r['expected_class']) if pool=='common6' else None
                    if new:v['fixed_common6_calibration'][p]=new
                    flat.append({'source':sid,'conditioning':r['conditioning'],'sample':r['sample'],'pool':pool,'policy':p,'R_kcal_mol':s['R_kcal_mol'],
                        'Ca_candidate':s['selected']['Ca'],'La_candidate':s['selected']['La'],'old_outcome':s['old_band_sensitivity']['outcome'],'new_common6_outcome':new['outcome'] if new else None})
            a_rows.append(v)
        # Copy ranges qualified by the actual5030 helper; do not fit on these10.
        for group in ('Ca','La','all10'):
            for pool in POOLS:
                for p in POLICIES:
                    values=[r['pools'][pool][p]['R_kcal_mol'] for r in a_rows if group=='all10' or r['conditioning']==group]
                    stored,=[x for x in prior['ranges'] if (x['conditioning_group'],x['pool'],x['policy'])==(group,pool,p)]
                    require((len(values),min(values),max(values),max(values)-min(values))==(stored['count'],stored['minimum_R_kcal_mol'],stored['maximum_R_kcal_mol'],stored['range_kcal_mol']),'5030 range regression differs')
        summaries={pool:{p:{'coverage':sum(r['pools'][pool][p]['R_kcal_mol'] is not None for r in rows),
            'old_band':counts([r['pools'][pool][p] for r in rows],'old_band_sensitivity'),
            'designated25':counts([r['pools'][pool][p] for r in rows],'designated25_calibration')} for p in POLICIES} for pool in POOLS}
        a_counts={p:{'old_common3':outcome_counts(r['pools']['common3'][p]['old_band_sensitivity']['outcome'] for r in a_rows),
            'old_common6':outcome_counts(r['pools']['common6'][p]['old_band_sensitivity']['outcome'] for r in a_rows),
            'fixed_new_common6':outcome_counts(r['fixed_common6_calibration'][p]['outcome'] for r in a_rows)} for p in POLICIES}
        canonical_flat=[{'source':r['source_case_id'],'protein':r['protein'],'expected_class':r['expected_class'],'pool':pool,'policy':p,
            'R_kcal_mol':r['pools'][pool][p]['R_kcal_mol'],'Ca_candidate':(r['pools'][pool][p]['selected'] or {}).get('Ca'),'La_candidate':(r['pools'][pool][p]['selected'] or {}).get('La'),
            'old_outcome':r['pools'][pool][p]['old_band_sensitivity']['outcome'],'designated25_outcome':r['pools'][pool][p]['designated25_calibration']['outcome']} for r in rows for pool in POOLS for p in POLICIES]
        changed={transition:{p:[{'record_id':r['record_id'],'source_case_id':r['source_case_id'],**r['changes'][transition][p]}
            for r in rows if r['changes'] and r['changes'][transition][p] is not None and r['changes'][transition][p]['before']!=r['changes'][transition][p]['after']]
            for p in POLICIES} for transition in ('common3_to_common4','common4_to_common6','common3_to_common6')}
        report={'config':cp,'status':'complete' if completed and final['all25_common6_available'] else 'partial','source_execution_status':se['status'],
            'canonical_rows':rows,'canonical_counts':summaries,'recorded_calibrations':final['calibrations'] if final else None,'old_calibrations':old['calibrations'],
            'changed_canonical_selections':changed,
            'common6_calibration_application_enabled':{p:calibration(final,'common6',p,completed) is not None for p in POLICIES},
            'A0AC_rows':a_rows,'A0AC_counts':a_counts,'A0AC_conditioning_ranges':prior['ranges'],'canonical_A0AC_overlap':overlap(final,ac),
            'saved_native_admission':admission(final),'A0AC_saved_admission':prior['saved_admission'],'contact_all_geometries_available':contact['all_candidate_geometries_available'],
            'actual_execution_receipts':{k:receipt(c[k]) for k in ('source_execution','canonical_search_execution','A0AC_execution','A0AC_search_execution')},'accounting':costs(c),
            'real_artifact_checks':{'actual5135_missing_collection_preserves25_and_old_baselines':True,'actual5030_all10_scores_selections_decisions_and_ranges_replay':True},
            'new_molecular_calls':0,'calibration_refit':False,'production_changed':False,'interpretation':'Consumed development sources; fixed designated25 bands only. Boundary proposals are constrained-domain results, not unconstrained minima.'}
        table(out/'CANONICAL.tsv',canonical_flat);table(out/'A0AC.tsv',flat);table(out/'A0AC_RANGES.tsv',prior['ranges']);write(out/'RESULT.json',report)
        ex.update(status=report['status'],result=pin(out/'RESULT.json'))
    except BaseException as e:ex.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:ex['wall_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',ex)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('freeze');a.add_argument('--source-config',required=True);a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('config');a.add_argument('sha');a=p.parse_args()
    if a.action=='freeze':freeze(a.source_config,a.output)
    else:execute(a.config,a.sha)
