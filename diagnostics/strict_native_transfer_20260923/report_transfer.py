"""Compact report and exact accounting for the already executed strict225 scope."""
import argparse,datetime,json,re,statistics,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new


def seconds(s):
    days,clock=(s.split('-',1) if '-' in s else ('0',s));v=list(map(float,clock.split(':')))
    return int(days)*86400+sum(x*60**i for i,x in enumerate(reversed(v)))


def build(comparison,submission,output):
    d=read_json(comparison);sub=read_json(submission);root=Path(comparison).parent
    ids=[r['job_id'] for r in sub['jobs']];fields=['JobIDRaw','State','ExitCode','AllocCPUS','ElapsedRaw','CPUTimeRAW','TotalCPU','MaxRSS','NodeList']
    cmd=['sacct','-n','-P','-j',','.join(ids),'--format='+','.join(fields)]
    proc=subprocess.run(cmd,capture_output=True,text=True,check=True)
    records=[dict(zip(fields,line.split('|'))) for line in proc.stdout.splitlines() if line.strip()]
    top=[r for r in records if r['JobIDRaw'] in ids]
    if len(top)!=4 or any(r['State'] in ('RUNNING','PENDING','COMPLETING') for r in top):raise RuntimeError('all four terminal accounting rows required')
    cpu_steps=[r for r in records if '.' in r['JobIDRaw']]
    new_cells=[];batches=[]
    for i in range(4):
        b=root/f'shard_{i}';c=read_json(b/'COLLECTION.json');e=read_json(b/'EXECUTION.json');batches.append({'shard':i,'execution':record(b/'EXECUTION.json'),'wall_seconds':e['wall_seconds'],'error':e['error']})
        for r in c['rows']:
            actual=r.get('actual');item={k:r[k] for k in ('task_id','case_id','candidate','metal','medium','status')};item['reason']=r.get('reason')
            if actual:
                rr=read_json(verify(actual['receipt']));start=datetime.datetime.fromisoformat(rr['started_at_utc']);end=datetime.datetime.fromisoformat(rr['finished_at_utc']);item['receipt_wall_seconds']=(end-start).total_seconds()
                text=verify(actual['output']).read_text();match=re.search(r'TOTAL RUN TIME:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+(\d+) seconds\s+(\d+) msec',text)
                item['ORCA_seconds']=sum(int(x)*factor for x,factor in zip(match.groups(),(86400,3600,60,1,.001))) if match else None
            new_cells.append(item)
    cost={'submission':record(submission),'job_accounting_command':cmd,'sacct':records,'jobs':top,
        'allocated_core_seconds':sum(int(r['CPUTimeRAW']) for r in top),'requested_GPU_seconds':0,
        'actual_CPU_seconds_sum_steps':sum(seconds(r['TotalCPU']) for r in cpu_steps if r['TotalCPU']),
        'actual_CPU_policy':'Sum .batch/.extern steps, exclude top-level duplicate aggregate; no srun computation steps in these jobs.',
        'new_scalar_task_denominator':2400,'new_completed':sum(r['status']=='complete' for r in new_cells),
        'new_status_counts':dict(__import__('collections').Counter(r['status'] for r in new_cells)),
        'sum_ORCA_seconds':sum(r.get('ORCA_seconds') or 0 for r in new_cells),'batches':batches,'cells':new_cells,
        'reused_scalar_cells':96,'reused_pools':8,'reuse_cost_in_new_totals':False,
        'local_preparation_wall_seconds':read_json(root/'PREFLIGHT.json')['preparation_wall_seconds'],
        'local_implementation_tests_reporting_additional_unmetered':True,'new_MACE_DFT_search_calls':0}
    write_new(root/'COSTS.json',cost)
    methods=('context_composite','union_precision','union_precision_strict');names={'context_composite':'Released static','union_precision':'Prior precision tenfold','union_precision_strict':'Strict fresh tenfold'}
    summary={'comparison':record(comparison),'costs':record(root/'COSTS.json'),'all225':{m:d['counts']['all225'][m] for m in methods},'conditioning':{k:{m:d['counts'][k][m] for m in methods} for k in ('La100','Ca125')},
        'aggregates':{k:{m:d['counts'][k][m] for m in methods} for k in ('La4','Ca5','balanced','La100_triples')},
        'matched':{k:{m:d['matched'][k][m] for m in methods[:-1]} for k in d['matched']},
        'changed_from_precision':[{'case_id':r['case_id'],'expected_class':r['expected_class'],'old':r['methods']['union_precision'],'new':r['methods']['union_precision_strict'],'delta_R':r['strict_delta_R']['operational'],'old_band_new_call':r['strict_old_reference_transfer']['operational']} for r in d['rows'] if r['methods']['union_precision']['decision']!=r['methods']['union_precision_strict']['decision']],
        'wrong_or_inconclusive':[{'case_id':r['case_id'],'expected_class':r['expected_class'],'result':r['methods']['union_precision_strict']} for r in d['rows'] if r['methods']['union_precision_strict']['outcome'] in ('wrong','inconclusive')]}
    comp=[abs(x[k]) for r in d['rows'] for x in r['strict_components_delta'] for k in ('vacuum_delta_kcal_mol','ALPB_delta_kcal_mol') if x[k] is not None]
    shifts=[abs(r['strict_delta_R']['operational']) for r in d['rows'] if r['strict_delta_R']['operational'] is not None]
    summary['old_loose_comparison']={'component_denominator':len(comp),'components_within_0p1':sum(x<=.1 for x in comp),'max_component_abs_kcal_mol':max(comp) if comp else None,
        'pool_denominator':len(shifts),'pools_within_0p2':sum(x<=.2 for x in shifts),'max_pool_abs_kcal_mol':max(shifts) if shifts else None}
    spread=[]
    for r in d['pools']:
        if r['pool'] not in ('La4','Ca5'):continue
        a=r['methods']['union_precision']['within_protein_range'];b=r['methods']['union_precision_strict']['within_protein_range']
        if a is not None and b is not None:spread.append({'root_case_id':r['root_case_id'],'conditioning':r['pool'],'old_range':a,'strict_range':b,'delta':b-a})
    summary['structural_spread']={'rows':spread,'count':len(spread),'shrink':sum(r['delta']<0 for r in spread),'grow':sum(r['delta']>0 for r in spread),'unchanged':sum(r['delta']==0 for r in spread),
        'median_old_range':statistics.median(r['old_range'] for r in spread) if spread else None,'median_strict_range':statistics.median(r['strict_range'] for r in spread) if spread else None}
    write_new(root/'SUMMARY.json',summary)
    def count(x):return f"{x['correct']} / {x['wrong']} / {x['inconclusive']} / {x['unavailable']}"
    a=d['counts']['all225']['union_precision_strict'];lines=['# Full225 strict-native scalar transfer', '',f"**{a['correct']} correct, {a['wrong']} wrong, {a['inconclusive']} inconclusive and {a['unavailable']} unavailable.** The method uses the same completed tenfold-precision geometry pools and MACE energies, with fresh strict native solvent calculations. Production is unchanged.", '', '| Method | Correct / wrong / inconclusive / unavailable |','|---|---|']
    lines += [f"| {names[m]} | {count(d['counts']['all225'][m])} |" for m in methods]
    lines += ['', '## Matched coverage and conditioning', '']
    for m in methods[:-1]:
        c=d['matched']['all225'][m];lines.append(f"On the same {c['common']} available sources, {names[m]}: {count(c['counts'][m])}; strict: {count(c['counts']['union_precision_strict'])}.")
    for k in ('La100','Ca125'):lines.append(f"{k}: strict {count(d['counts'][k]['union_precision_strict'])}.")
    lines += ['', '| Strict aggregate | Released | Prior precision | Strict fresh |','|---|---|---|---|']
    lines += [f"| {k} | "+' | '.join(count(d['counts'][k][m]) for m in methods)+' |' for k in ('La4','Ca5','balanced','La100_triples')]
    lines += ['', 'Each aggregate requires every declared member. The100 triples and225 folds are correlated samples of25 consumed proteins, not independent biological validations. All17 original preparation exclusions remain; no failed cell is filled with an old energy.', '', '## Decision changes versus prior precision', '']
    for r in summary['changed_from_precision']:lines.append(f"- `{r['case_id']}`: {r['old']['decision']} → {r['new']['decision']}; ΔR={r['delta_R']:.9f} model kcal/mol; new score under old bands: {r['old_band_new_call']}.")
    if not summary['changed_from_precision']:lines.append('No decision changes.')
    lines += ['', 'Remaining wrong/inconclusive calls:', '']+[f"- `{r['case_id']}`: {r['result']['decision']} ({r['expected_class']}-class reference), R={r['result']['R']:.9f}." for r in summary['wrong_or_inconclusive']]
    if not summary['wrong_or_inconclusive']:lines.append('None.')
    n=summary['old_loose_comparison'];s=summary['structural_spread'];lines += ['', '## Numerical interpretation', '',f"All geometry and MACE cells were reused unchanged. {n['components_within_0p1']}/{n['component_denominator']} scalar cells remain within0.1 kcal/mol of their original loose result; {n['pools_within_0p2']}/{n['pool_denominator']} pools remain within0.2. Maximum absolute component/pool shifts: {n['max_component_abs_kcal_mol']:.9f}/{n['max_pool_abs_kcal_mol']:.9f} kcal/mol. These comparisons describe the change; reproducing the loose solution is not an acceptance criterion.", '',f"Among{s['count']} complete conditioning groups, spread shrinks in{s['shrink']}, grows in{s['grow']} and is unchanged in{s['unchanged']}. Median range: {s['median_old_range']:.6f}→{s['median_strict_range']:.6f} model kcal/mol. Whole-group medians, raw contrasts, selected candidates and each component remain in the full JSON.", '', 'The frozen fresh strict32 reference uses only the original25 calibration sources. No threshold was refitted here. Strict32 fresh/seeded agreement is a separate numerical qualification; these transfer sources receive one fresh start each. Scalar consistency does not qualify forces, prove a unique electronic solution or turn electronic contrasts into binding free energies.', '', '## Actual cost and reproducibility', '',f"Four CPU-only32-core/64GiB jobs: {', '.join(ids)}. New calls: {cost['new_completed']}/2400 complete; eight whole pools/96 cells reused separately. No new MACE, searches, DFT or preparation chemistry. Allocated {cost['allocated_core_seconds']} core-seconds; actual accounting steps {cost['actual_CPU_seconds_sum_steps']:.3f} CPU-seconds; zero requested GPU-seconds. Summed native ORCA time {cost['sum_ORCA_seconds']:.3f}s. Full job/step receipts, failures and per-cell times are in COSTS.json.", '',f"Local source/manifest audit took{cost['local_preparation_wall_seconds']:.3f}s; implementation, tests and reporting add unmetered local work. Reused qualification/comparator allocations and historical MACE work are excluded from these new-job totals, not declared free. This is not a matched hardware speed comparison.", '', 'Full source/aggregate results: `workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json`; compact summary, costs and pinned receipts beside it. See COMMANDS.md for the exact no-molecule comparison operation.']
    Path(output).write_text('\n'.join(lines)+'\n');return {'report':record(output),'summary':record(root/'SUMMARY.json'),'costs':record(root/'COSTS.json')}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('comparison','submission','output'):p.add_argument('--'+k,required=True,type=Path)
    print(json.dumps(build(**vars(p.parse_args())),indent=2))
