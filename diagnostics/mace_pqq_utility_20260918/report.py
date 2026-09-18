#!/usr/bin/env python3
"""Final utility report from completed real jobs; refuses partial success claims."""
import argparse
import csv
import json
from pathlib import Path
import statistics
import subprocess
from benchmark import read, record, verify, write, collect


def cpu_seconds(value):
    days, _, clock = value.rpartition('-')
    parts = [float(x) for x in clock.split(':')]
    if len(parts) > 3:raise ValueError('invalid scheduler CPU time')
    return (int(days)*86400 if days else 0)+sum(v*60**i for i,v in enumerate(reversed(parts)))


def scheduler(jobs):
    command=['sacct','-j',','.join(jobs),'-n','-P',
             '--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,AllocTRES']
    result=subprocess.run(command,check=True,capture_output=True,text=True)
    rows={}
    for line in result.stdout.splitlines():
        values=line.split('|')
        if len(values)!=7:raise ValueError('unexpected scheduler columns')
        rows[values[0]]=dict(zip(('job_id','state','elapsed_seconds','allocated_cpus',
                                 'CPU_time','max_RSS','allocated_TRES'),values))
    for job in jobs:
        if job not in rows or rows[job]['state'] not in ('COMPLETED','FAILED','CANCELLED','OUT_OF_MEMORY','TIMEOUT'):
            raise ValueError('job not terminal; final report unavailable: '+job)
    return {'command':command,'raw_output':result.stdout,'rows':rows}


def generate(workspace, output):
    work=Path(workspace).resolve();run=work/'timing_v2';out=Path(output).resolve()
    m=read(run/'manifest.json');accuracy=read(work/'accuracy_v1/result.json')
    for ref in accuracy['source_records']:verify(ref)
    if accuracy['counts']!={'calibration':{'total':25,'DFT_correct':25,'MACE_correct':25,'MACE_unavailable':0},
                           'transfer':{'total':3,'DFT_correct':3,'MACE_correct':2,'MACE_unavailable':1}}:
        raise ValueError('reference inventory differs; investigate before reporting')
    submissions={k:read(run/(k+'_submission.json')) for k in ('DFT','MACE')}
    failed_submission=read(work/'timing_v1/DFT_submission.json')
    jobs={k:r['stdout'].strip().split(';')[0] for k,r in submissions.items()}
    failed_job=failed_submission['stdout'].strip().split(';')[0]
    accounting=scheduler([*jobs.values(),failed_job])
    for method in ('DFT','MACE'):
        paths=list((run/method).glob('*/timing.json'))
        if len(paths)!=25 or any(read(p)['status']!='complete' for p in paths):
            raise ValueError('all 25 successful cases per method required for final comparison')
    out.mkdir(parents=True,exist_ok=False)
    collect(run/'manifest.json',out/'paired.json');paired=read(out/'paired.json')
    if paired['complete_pairs']!=25:raise ValueError('incomplete paired inventory')
    rows=[];endpoints={'DFT':[],'MACE':[]};mace_results=[]
    for c in paired['rows']:
        row={k:c[k] for k in ('case_id','expected_class')}
        for method in ('DFT','MACE'):
            r=c[method];result=r['result'];case=run/method/c['case_id']
            if method=='DFT':
                for metal in ('Ca','La'):
                    p=case/(metal+'.out.execution.json');e=read(p)
                    verify(e['artifacts']['output'])
                    if not (e['returncode']==0 and e['scf_converged'] and e['normal_termination']):
                        raise ValueError('unsuccessful DFT execution')
                    endpoints[method].append(record(p))
            else:
                score=read(verify(result['scientific_report']))
                if not score['numerical_gate_pass']:raise ValueError('failed MACE numerical checks')
                for p in (case/'score/execution').glob('*/attempt_*/receipt.json'):
                    e=read(p);computed=read(verify(e['result']));verify(e['resource_usage'])
                    if e['returncode']!=0 or computed['status']!='computed':raise ValueError('failed MACE execution')
                    mace_results.append(computed)
                    endpoints[method].append(record(p))
            row.update({method+'_seconds':r['wall_seconds'],method+'_score':result['score'],
                        method+'_class':result['class'],method+'_difference':r['score_difference'],
                        method+'_reproduced':r['reproduced'],
                        method+'_literal_correct':r['literal_class_correct']})
        row['speedup']=row['DFT_seconds']/row['MACE_seconds'];rows.append(row)
    if any(len(v)!=50 for v in endpoints.values()):raise ValueError('fresh endpoint count differs')
    stats={}
    for method,job in jobs.items():
        s=accounting['rows'][job];times=[r[method+'_seconds'] for r in rows]
        resources=run/(method+'_job_'+job+'.resources.txt')
        if not resources.is_file():raise ValueError('final process resource receipt missing')
        stats[method]={'job_id':job,'job_state':s['state'],
                       'latency_seconds':{'min':min(times),'median':statistics.median(times),'max':max(times),'sum':sum(times)},
                       'allocation_seconds':int(s['elapsed_seconds']),
                       'allocated_cpus':int(s['allocated_cpus']),
                       'allocated_core_seconds':int(s['elapsed_seconds'])*int(s['allocated_cpus']),
                       'scheduler_CPU_seconds':cpu_seconds(s['CPU_time']),
                       'GPU_allocation_seconds':int(s['elapsed_seconds']) if method=='MACE' else 0,
                       'process_resources':record(resources),
                       'per_case_phase_totals_seconds':{
                           phase:sum(r[method]['result'][phase] for r in paired['rows'])
                           for phase in ('preparation_seconds','execution_seconds','report_seconds')},
                       'allocated_host_memory_MiB':sorted({int(r[method]['allocated_memory_MiB']) for r in paired['rows']}),
                       'literal_correct':sum(r[method+'_literal_correct'] for r in rows),
                       'within_declared_reproduction_tolerance':sum(r[method+'_reproduced'] for r in rows),
                       'max_score_difference':max(abs(r[method+'_difference']) for r in rows)}
    stats['MACE'].update(
        summed_model_evaluation_seconds=sum(r['evaluation_seconds'] for r in mace_results),
        peak_cuda_allocated_bytes=max(r['peak_cuda_allocated_bytes'] for r in mace_results),
        peak_cuda_reserved_bytes=max(r['peak_cuda_reserved_bytes'] for r in mace_results),
        maximum_worker_process_RSS_KiB=max(r['peak_host_RSS_KiB'] for r in mace_results))
    speed=statistics.median(r['speedup'] for r in rows)
    median_latency_ratio=stats['DFT']['latency_seconds']['median']/stats['MACE']['latency_seconds']['median']
    total_ratio=stats['DFT']['latency_seconds']['sum']/stats['MACE']['latency_seconds']['sum']
    raw_job_ratio=stats['DFT']['allocation_seconds']/stats['MACE']['allocation_seconds']
    speed_pass=median_latency_ratio>=m['speed_target'] and total_ratio>1
    reproduction_pass=all(s['within_declared_reproduction_tolerance']==25 for s in stats.values())
    fidelity_pass=stats['MACE']['literal_correct']==25 and reproduction_pass
    recommendation=('use_frozen_masked_MACE_as_opt_in_PQQ_screen' if speed_pass and fidelity_pass else
                    'do_not_adopt_for_claimed_PQQ_speed_utility_on_this_evidence')
    result={'status':'complete','manifest':record(run/'manifest.json'),
            'reference_fidelity':record(work/'accuracy_v1/result.json'),
            'rows':rows,'methods':stats,'endpoint_receipts':endpoints,
            'median_pair_speedup':speed,'total_case_speedup':total_ratio,
            'ratio_of_median_latencies':median_latency_ratio,
            'whole_allocation_speedup':raw_job_ratio,
            'speed_target':m['speed_target'],'speed_gate_pass':speed_pass,'fidelity_gate_pass':fidelity_pass,
            'recommendation':recommendation,'baseline_changed':False,'default_promotion':False,
            'failed_startup_job':accounting['rows'][failed_job],
            'failed_startup_receipts':[record(p) for p in (work/'timing_v1/DFT').glob('*/*.out.execution.json')],
            'development_preflight':record(work/'preflight_v2_receipt.json'),
            'raw_source_preparation_newly_timed':False,
            'limits':['All labels/results consumed; canonical calibration is not prospective accuracy.',
                      'One fresh timing per case; shared-host variation has not been quantified by repetitions.',
                      'Different representations and resources: whole-chain MACE on GPU, fixed-core DFT on CPU.',
                      'No claim that GPU-seconds equal CPU-seconds or that latency ratio equals monetary savings.',
                      'Prepared-input-to-score timing excludes original folding and protonation.',
                      '1KB0 unsupported by MACE; no silent DFT fallback or universal replacement.',
                      'Literal threshold crossings within numerical tolerance remain visible.']}
    write(out/'accounting.json',accounting);write(out/'result.json',result)
    with (out/'scores_and_timings.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    lines=['# Frozen masked-MACE PQQ utility comparison','',f'Recommendation: **{recommendation}**.',
           '',f'Ratio of median scoring times **{median_latency_ratio:.3f}×**; '
           f'median paired speedup **{speed:.3f}×**; total case speedup **{total_ratio:.3f}×**; '
           f'whole-allocation speedup **{raw_job_ratio:.3f}×**.',
           '', '| Method | Literal correct / 25 reruns | Reproduced within 0.01 / 25 | Median seconds | Total seconds |',
           '|---|---:|---:|---:|---:|']
    for name,s in stats.items():
        t=s['latency_seconds'];lines.append(f"| {name} | {s['literal_correct']} | {s['within_declared_reproduction_tolerance']} | {t['median']:.2f} | {t['sum']:.2f} |")
    lines+=['','Archived references:25/25canonical both; crystalsDFT3/3,MACE2/3 with1KB0unsupported.',
            'Each method retains its own calibration; the scales are not interchangeable.',
            '', '## Limits','',*['- '+x for x in result['limits']],
            '', 'Exact costs, startup failures, endpoint receipts and all28archived references are linked in result.json.',
            'All25fresh paired rows are in scores_and_timings.tsv. Production default unchanged.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'result':record(out/'result.json'),'recommendation':recommendation,
            'ratio_of_median_latencies':median_latency_ratio,'median_pair_speedup':speed}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',required=True);p.add_argument('--output',required=True)
    print(json.dumps(generate(**vars(p.parse_args())),indent=2))
