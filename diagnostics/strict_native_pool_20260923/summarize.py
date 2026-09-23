"""Actual strict32 scalar agreement, raw discrimination and measured costs."""
import json,re,subprocess,sys
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
W=ROOT/'workspaces/strict_native_pool_20260923/run_v1';comparison=read_json(W/'COMPARISON.json')

def runtime(text):
    m=re.search(r'TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec',text)
    if not m:return None
    d,h,minute,s,ms=map(int,m.groups());return d*86400+h*3600+minute*60+s+ms/1000
branches={}
for b in ('fresh','cold_seed'):
    c=read_json(W/b/'COLLECTION.json');execution=read_json(W/b/'EXECUTION.json');times=[];source_times=[];density={k:0 for k in ('energy','max_density','rms_density')};m=read_json(W/b/'manifest.json');ti={t['task_id']:t for t in m['tasks']}
    for r in c['rows']:
        if r.get('details'):
            for k in density:density[k]+=bool(r['details'].get(k,{}).get('within_printed_tolerance'))
        if r.get('reused'):continue
        path=verify(r['actual']['output']) if r.get('actual') else Path(ti[r['task_id']]['output_path'])
        t=runtime(path.read_text()) if path.exists() else None
        times.append({'task_id':r['task_id'],'ORCA_seconds':t,'status':r['status']})
        src=ti[r['task_id']]['source'];source_times.append(runtime(verify(src['output']).read_text()))
    branches[b]={'new_calls':len(ti),'exact_reused':len(m['reused']),'logical_denominator':384,'complete':c['complete'],
        'status_counts':dict(Counter(r['status'] for r in c['rows'])),'printed_density_checks':density,
        'executor_wall_seconds':execution['wall_seconds'],'executor_allocated_core_seconds':execution['wall_seconds']*execution['allocated_cpus'],
        'summed_new_ORCA_seconds':sum(t['ORCA_seconds'] for t in times if t['ORCA_seconds'] is not None),
        'new_ORCA_times_available':sum(t['ORCA_seconds'] is not None for t in times),
        'historical_loose_cold_ORCA_seconds_for_new_tasks':sum(t for t in source_times if t is not None),
        'production_recipe':('one fresh strict call' if b=='fresh' else 'loose cold plus one strict continuation'),
        'production_call_count_per_cell':1 if b=='fresh' else 2,'historical_seed_cost_in_new_job':False,'tasks':times}
raw=subprocess.check_output(['sacct','-j','1211270','--parsable2','--noheader','--format=JobID,State,ElapsedRaw,AllocCPUS,ReqMem,AllocTRES,MaxRSS'],text=True)
r=next(x.split('|') for x in raw.splitlines() if x.split('|')[0]=='1211270')
if r[1] not in ('COMPLETED','FAILED'):raise RuntimeError('job not terminal')
cost={'job_id':'1211270','state':r[1],'wall_seconds':int(r[2]),'allocated_cpus':int(r[3]),'allocated_core_seconds':int(r[2])*int(r[3]),'requested_memory':r[4],'GPU_seconds':0,
    'new_scalar_calls':sum(b['new_calls'] for b in branches.values()),'exact_reused_calls':10,'branches':branches,'sacct':raw,'local_preparation_tests_report_in_allocation':False,
    'new_MACE_DFT_optimization_calls':0,'branches_are_development_validation_not_both_routine':True}
write_new(W/'COSTS.json',cost)
reference=read_json(W/'REFERENCES.json');failedcells=[r for r in comparison['cells'] if r['initialization_agreement_pass'] is not True]
failedcases=[{'case_id':r['case_id'],'cells_agree':r['cells_agree'],'same_geometry_R_checks':r['same_geometry_R_checks'],'pool_R_checks':r['pool_R_checks']} for r in comparison['rows'] if r['qualified_status']!='available']
summary={'comparison':record(W/'COMPARISON.json'),'reference':record(W/'REFERENCES.json'),'cost':record(W/'COSTS.json'),
    'counts':comparison['counts'],'agreeing_cells':comparison['agreeing_cells'],'cell_denominator':384,
    'max_cell_difference_kcal_mol':max((abs(r['cold_seed_minus_fresh_kcal_mol']) for r in comparison['cells'] if r['cold_seed_minus_fresh_kcal_mol'] is not None),default=None),
    'max_same_geometry_R_difference_kcal_mol':max((abs(c['cold_seed_minus_fresh_R_kcal_mol']) for r in comparison['rows'] for c in r['same_geometry_R_checks'].values() if c['cold_seed_minus_fresh_R_kcal_mol'] is not None),default=None),
    'max_pool_R_difference_kcal_mol':max((abs(c['cold_seed_minus_fresh_R_kcal_mol']) for r in comparison['rows'] for c in r['pool_R_checks'].values() if c['cold_seed_minus_fresh_R_kcal_mol'] is not None),default=None),
    'failed_cells':failedcells,'failed_cases':failedcases,
    'references':{b:reference['branches'][b]['variants'] for b in ('fresh','cold_seed')},'production_changed':False}

summary['candidate_choice_changes']=[]
for r in comparison['rows']:
    for b in ('fresh','cold_seed'):
        new=r['branches'][b]['pool'];old=r['cold_pool']
        if new['status']!='available' or old['status']!='available':continue
        for z in ('Ca','La'):
            for v in ('mathematical','operational'):
                field=v+'_candidate'
                if new['rows'][z][field]!=old['rows'][z][field]:
                    summary['candidate_choice_changes'].append({'case_id':r['case_id'],'branch':b,'metal':z,'variant':v,'cold':old['rows'][z][field],'strict':new['rows'][z][field]})
summary['branch_physical_observation']={b:{k:v for k,v in costs.items() if k not in ('tasks',)} for b,costs in branches.items()}

write_new(ROOT/'diagnostics/strict_native_pool_20260923/RESULT.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k not in ('failed_cells','failed_cases','references')},indent=2))
print('failed cells',[(r['case_id'],r['candidate'],r['metal'],r['medium'],r['cold_seed_minus_fresh_kcal_mol']) for r in failedcells])
