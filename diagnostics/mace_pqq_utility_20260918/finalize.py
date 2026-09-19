#!/usr/bin/env python3
"""Join the original and cached-runtime comparisons without double-counting DFT."""
import argparse
import csv
import json
from pathlib import Path
from benchmark import read, record, write
from report import generate, scheduler, cpu_seconds


def finalize(workspace, output):
    work=Path(workspace).resolve();out=Path(output).resolve()
    runs={'DFT':work/'timing_v2','MACE_original':work/'timing_v2','MACE_cached':work/'timing_cached_v1'}
    jobs={name:read(path/('DFT_submission.json' if name=='DFT' else 'MACE_submission.json'))['stdout'].strip().split(';')[0]
          for name,path in runs.items()}
    jobs.update(failed_DFT_startup=read(work/'timing_v1/DFT_submission.json')['stdout'].strip().split(';')[0],
                report_only_check=read(work/'report_cache_matched_v1/submission.json')['stdout'].strip().split(';')[0])
    account=scheduler(list(jobs.values()))
    # Fail before creating outputs if either complete comparison is unavailable.
    for label,directory in runs.items():
        method='DFT' if label=='DFT' else 'MACE'
        receipts=list((directory/method).glob('*/timing.json'))
        if len(receipts)!=25 or any(read(p)['status']!='complete' for p in receipts):
            raise ValueError('incomplete campaign: '+label)
    out.mkdir(parents=True,exist_ok=False)
    generate(work,out/'original')
    generate(work,out/'cached',mace_run=work/'timing_cached_v1')
    old=read(out/'original/result.json');new=read(out/'cached/result.json')
    if old['method_manifests']['DFT']!=new['method_manifests']['DFT']:
        raise ValueError('DFT results were not reused identically')
    old_rows={r['case_id']:r for r in old['rows']};rows=[]
    for r in new['rows']:
        a=old_rows.pop(r['case_id'])
        if (a['expected_class'],a['DFT_score'],a['DFT_seconds'])!=(r['expected_class'],r['DFT_score'],r['DFT_seconds']):
            raise ValueError('matched baseline case differs')
        rows.append({'case_id':r['case_id'],'expected_class':r['expected_class'],
                     'DFT_seconds':r['DFT_seconds'],'original_MACE_seconds':a['MACE_seconds'],
                     'cached_MACE_seconds':r['MACE_seconds'],
                     'DFT_S_kcal_mol':r['DFT_score'],'DFT_class':r['DFT_class'],
                     'original_MACE_R_model_kcal':a['MACE_score'],'original_MACE_class':a['MACE_class'],
                     'cached_MACE_R_model_kcal':r['MACE_score'],'cached_MACE_class':r['MACE_class'],
                     'cache_score_change_model_kcal':r['MACE_score']-a['MACE_score']})
    if old_rows:raise ValueError('case inventories differ')
    accuracy_by_class={}
    for expected in ('Ca','La'):
        selected=[r for r in rows if r['expected_class']==expected]
        accuracy_by_class[expected]={'total':len(selected)}
        for method in ('DFT','original_MACE','cached_MACE'):
            accuracy_by_class[expected][method]={
                'correct':sum(r[method+'_class']==expected for r in selected),
                'inconclusive':sum(r[method+'_class']=='inconclusive' for r in selected),
                'opposite_class':sum(r[method+'_class'] in ('Ca','La') and r[method+'_class']!=expected for r in selected)}
    costs={}
    for name,job in jobs.items():
        s=account['rows'][job];elapsed=int(s['elapsed_seconds']);cpus=int(s['allocated_cpus'])
        costs[name]={'job_id':job,'state':s['state'],'allocation_seconds':elapsed,'Slurm_CPUs':cpus,
                     'allocated_CPU_seconds':elapsed*cpus,'reported_CPU_seconds':cpu_seconds(s['CPU_time']),
                     'GPU_allocation_seconds':elapsed if name.startswith('MACE') else 0}
    result={'status':'complete','original_comparison':record(out/'original/result.json'),
            'cached_comparison':record(out/'cached/result.json'),
            'reference_fidelity':new['reference_fidelity'],'recommendation':new['recommendation'],
            'scientific_protocol_changed':False,'production_default_changed':False,
            'rows':rows,'accuracy_by_class':accuracy_by_class,'unique_cluster_costs':costs,
            'development_cluster_totals':{key:sum(v[key] for v in costs.values()) for key in
                                        ('allocated_CPU_seconds','reported_CPU_seconds','GPU_allocation_seconds')},
            'local_development_receipts':[record(work/p) for p in
                ('preflight_v2_receipt.json','report_profile_v1/receipt.json','cached_report_v1/receipt.json')],
            'local_development_CPU_total_status':'partial receipts; not presented as a complete total',
            'original_speed_gate_pass':old['speed_gate_pass'],'cached_speed_gate_pass':new['speed_gate_pass'],
            'original_fidelity_gate_pass':old['fidelity_gate_pass'],'cached_fidelity_gate_pass':new['fidelity_gate_pass'],
            'max_cache_score_change_model_kcal':max(abs(r['cache_score_change_model_kcal']) for r in rows),
            'CPU_unit_note':'Slurm allocated CPU units; the host has two hardware threads per physical core.',
            'new_model_evaluations':100,'successful_new_DFT_endpoints':50,'failed_DFT_startup_attempts':4,
            'preparation_scope':'already prepared labelled structures; original folding/protonation not newly timed'}
    write(out/'accounting.json',account);write(out/'result.json',result)
    with (out/'all_scores_and_timings.tsv').open('x') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    baseline=new['methods']['DFT'];original=old['methods']['MACE'];cached=new['methods']['MACE']
    lines=['# PQQ reference fidelity and measured MACE utility','',
           f"Recommendation: **{new['recommendation']}**. The production DFT default is unchanged.",'',
           '## Prediction fidelity','',
           '| Reference class | Cases | DFT correct | Original MACE correct | Cached MACE correct |',
           '|---|---:|---:|---:|---:|']
    for expected,a in accuracy_by_class.items():
        lines.append(f"| {expected} | {a['total']} | {a['DFT']['correct']} | {a['original_MACE']['correct']} | {a['cached_MACE']['correct']} |")
    lines+=['', 'Archived crystal transfers: DFT 3/3; MACE 2/2 scored, with 1KB0 unsupported (2/3 of the complete transfer set). '
            'Reference replay preserves performance; it does not demonstrate improved accuracy on unseen proteins.',
            '', '## Full workflow timing','',
           '| Workflow | Median seconds / protein | Total case seconds | Fresh correct / 25 | Scores reproduced / 25 |',
           '|---|---:|---:|---:|---:|']
    for label,s in [('DFT',baseline),('Original masked MACE',original),('Parser-cache masked MACE',cached)]:
        t=s['latency_seconds'];lines.append(f"| {label} | {t['median']:.2f} | {t['sum']:.2f} | {s['literal_correct']} | {s['within_declared_reproduction_tolerance']} |")
    lines+=['',f"Ratio of median times versus DFT: original **{old['ratio_of_median_latencies']:.3f}x**, "
            f"cached **{new['ratio_of_median_latencies']:.3f}x**. The predeclared engineering target was 1.5x plus lower total time.",
            '', '## What this establishes','',
            '- The frozen masked-MACE descriptor retains the labelled PQQ reference task; the cache changes only repeated coordinate parsing.',
            '- Archived calibration: 25/25 for both methods. Consumed crystals: DFT 3/3; MACE 2/3, with 1KB0 unsupported.',
            '- The benchmark uses each method\'s own frozen calibration. These are consumed references, not new independent biological validation.',
            '- PQQ functional association is the target. Broad affinity discrimination, physiological occupancy and within-lanthanide preference are not established.',
            '', '## Numerical boundary and coverage','',
            'Literal decisions and score-reproduction checks are reported separately. Small DFT rerun changes can cross a calibration extremum; '
            'the full table retains every inconclusive call. No threshold was moved. 1KB0 remains unscorable by whole-chain MACE, with no baseline substitution.',
            '', '## Measured resources','',
            '| Job purpose | Allocation seconds | Allocated CPU-seconds | GPU allocation seconds |',
            '|---|---:|---:|---:|']
    for label,c in costs.items():lines.append(f"| {label} | {c['allocation_seconds']} | {c['allocated_CPU_seconds']} | {c['GPU_allocation_seconds']} |")
    lines+=['','DFT used 32 Slurm CPU units; each MACE run used one A5000 and 16 CPU units on the same host. '
            'Each main job requested 64474 MiB. CPU units are scheduler units, not an asserted count of physical cores. '
            'A GPU-second is not equated with a CPU-second or a monetary price.',
            '', 'Timing includes validation, loading, execution and reporting from existing prepared inputs. '
            'It excludes queue waiting and original folding/protonation. Each workflow has one fresh run per case; '
            'shared-host variability was not estimated through repetitions.',
            '', 'All 150 successful fresh endpoint receipts, four failed MPI startups, exact scheduler costs, '
            'local development receipts, and original/cached comparisons remain linked in result.json. '
            'No new scientific model, reference, threshold or production default was introduced.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'result':record(out/'result.json'),'report':record(out/'REPORT.md'),
            'recommendation':result['recommendation'],'cached_median_speedup':new['ratio_of_median_latencies']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',required=True);p.add_argument('--output',required=True)
    print(json.dumps(finalize(**vars(p.parse_args())),indent=2))
