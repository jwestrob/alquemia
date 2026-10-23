#!/usr/bin/env python3
"""One saved review of the terminal eighteen-control collection; zero native calls."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'runtime_recovery_v1/REFERENCE.py'
SOURCE_SHA='8ee0d091373913fe7e95e3731503a38c9da2c1435560adaee310b95a92b736c3'
s=importlib.util.spec_from_file_location('qualified_local_reference',SOURCE)
r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
b=r.b


def main(path,out):
    tick,cpu=time.monotonic(),time.process_time()
    b.check(b.pin(SOURCE)['sha256']==SOURCE_SHA,'Frozen collector source changed')
    actual=json.loads(Path(path).read_text());b.check(actual['implementation']==b.pin(SOURCE) and actual['stage']=='controls','Unexpected collection producer/stage')
    p=r.load(b.verify(actual['panel']));inv=r.read(actual['inventory']);cells=[c for c in p['cells'] if c['stage']=='controls']
    b.check(len(cells)==actual['declared']==18 and set(inv)==set(actual['rows'])=={c['cell_id'] for c in cells},'Eighteen-cell denominator differs')
    out.mkdir(parents=True,exist_ok=False);admitted={};details={};recovered=0;native_cpu=0.
    for c in cells:
        cid=c['cell_id'];row=actual['rows'][cid]
        if row['status']!='complete':
            details[cid]=dict(status='unavailable',collector_record=row);continue
        b.check(row['result']==b.pin(inv[cid]),'Collection result authority differs')
        native=r.read(row['native_result']);start=r.read(native['start']);raw=b.verify(native['output']).read_text()
        b.check(native['cell']==c and native['panel']==actual['panel'] and native['returncode']==0,'Native physical inventory differs')
        b.check(start['cell']==c and start['input']==native['input'] and start['panel']==actual['panel'] and
                start['ranks']==4 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Native input/resources differ')
        b.check(b.verify(native['input']).read_text()==r.deck(p,c,4,start['maxcore_MiB']),'Exact executed deck differs')
        b.check(native['native_program_invocations_attempted']==1,'Native invocation denominator differs')
        producer=p['adapter'] if cid in r.old.CONTROL_SCOUTS else b.pin(SOURCE)
        b.check(native['adapter']==row['native_producer']==producer,'Original/versioned native producer differs')
        q=r.parse(raw,c,0);serialized=json.loads(json.dumps(q))
        b.check(q['status']=='complete' and all(row[k]==v for k,v in serialized.items()),'Raw energy/SCF/local/CC/triples replay differs: '+cid)
        if cid in r.old.CONTROL_SCOUTS:
            saved=r.read(row['result']);b.check(saved['status']=='saved_native_recovery_complete' and
                saved['original_result']==row['native_result'] and saved['implementation']==b.pin(r.rec.__file__),'Recovered original authority differs')
            recovered+=1
        elapsed=native['elapsed_seconds'];cost=native['allocated_cpu_seconds'];b.check(abs(cost-4*elapsed)<1e-6,'Native allocation cost differs')
        native_cpu+=cost;admitted[cid]=row
        details[cid]=dict(status='raw_replay_passed',result=row['result'],native_result=row['native_result'],native_producer=producer,
            output=native['output'],input=native['input'],start=native['start'],job_id=native['job_id'],
            energy_hartree=q['energy_hartree'],CC_terminal_residuals=q['CC_terminal_residuals'],triples_terminal_deltas=q['triples_terminal_deltas'],
            T1_diagnostic=q['T1_diagnostic'],native_scf_stop_signals=q['native_scf_stop_signals'],printed_scf_criteria=q['printed_scf_criteria'],
            warning_and_core_lines=[x for x in q['local_diagnostics'] if re.search('warning|core',x['text'],re.I)],
            elapsed_seconds=elapsed,elapsed_times_four_MPI_CPU_seconds=cost)
    b.check(len(admitted)==actual['complete'],'Complete count differs')
    score=None;failed={}
    if len(admitted)==18:
        b.check(recovered==4,'Original four producer denominator differs')
        score=r.control_score(p,admitted);b.check(len(score['checks'])==35,'Frozen comparison denominator differs')
        b.check(all(actual[k]==v for k,v in score.items()),'Frozen thirty-five-check score replay differs')
        failed={k:v for k,v in score['checks'].items() if not v}
    b.check(abs(native_cpu-actual['native_allocated_CPU_seconds'])<1e-6,'Collector native cost sum differs')
    result=dict(status='passed_saved_control_collection_review' if len(admitted)==18 else 'incomplete_control_collection_review',
        implementation=b.pin(__file__),collector_result=b.pin(path),collector_source=b.pin(SOURCE),parser=b.pin(r.rec.__file__),
        original_native_input_adapter=p['adapter'],panel=actual['panel'],inventory=actual['inventory'],rows=details,
        declared_controls=18,admitted_controls=len(admitted),recovered_original_scouts=recovered,new_remaining_native=len(admitted)-recovered,
        declared_comparisons=35,comparison_score=score,failed_comparison_names=list(failed),
        all_method_controls_passed=bool(score and score['all_controls_passed']),new_native_calls=0,
        native_elapsed_times_allocated_rank_CPU_seconds=native_cpu,
        collector_saved_cost={k:actual[k] for k in ('wall_seconds','parent_cpu_seconds')},
        review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,local_threads=1),
        cluster_release_authorized_by_this_review=False,physical_accuracy_qualified=False,cluster_CP_or_B_qualified=False)
    rp=b.write(out/'REVIEW.json',result)
    lines=['# Eighteen correlated-reference controls','',f'Raw admission: {len(admitted)}/18. Preserved original recovered scouts: {recovered}/4.',
           f'Method comparisons: {35-len(failed)}/35 pass.' if score else 'Method comparison denominator unavailable; no partial qualification.',
           '',f'Native elapsed × allocated MPI ranks: {native_cpu:.3f} CPU-s ({native_cpu/3600:.3f} CPU-h). This is a runtime allocation proxy, separate from scheduler accounting.',
           '', 'All native SCF, correlation, core/semicore and T1 warnings remain in the structured review. Small-system qualification is a combined RI/local-correlation check; it does not certify cluster locality error or physical reference accuracy. Root retains the largest-cluster scout decision.']
    if score:
        lines.extend(['','| Comparison | Local minus canonical, kcal/mol |','| --- | ---: |'])
        lines.extend(f'| CP {k} | {v:+.9f} |' for k,v in score['CP_errors_kcal_mol'].items())
        for label,key in [('Paired CP','paired_CP_local_canonical'),('Paired donor–water','paired_exchange_local_canonical')]:
            lines.extend(f'| {label} {k} | {v[0]-v[1]:+.9f} |' for k,v in score[key].items())
        lines.extend(['','Failed frozen comparisons: '+(', '.join(failed) if failed else 'none')+'.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(review=rp,status=result['status'],all_controls_passed=result['all_method_controls_passed'])))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--result',required=True);a.add_argument('--out',type=Path,required=True);x=a.parse_args();main(x.result,x.out)
