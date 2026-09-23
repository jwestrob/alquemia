"""Compare the actual fresh two-crystal execution to pinned prior artifacts."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,xyz
import nikasha_pool as pool
from accommodation_nonlinear import relative_components


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    run=args.run;result=read_json(run/'result.json');prep=read_json(run/'source_preparation/preparation.json')
    original_path=ROOT/'workspaces/adaptive_completion_20260922/original30_pool_v1/final_collection.json'
    original=read_json(original_path);fresh=read_json(run/'pool/collection.json');pm=read_json(run/'proposals/manifest.json')
    rows=[]
    for c in fresh['cases']:
        old=next(x for x in original['cases'] if x['case_id']==c['case_id']);src=next(x for x in prep['cases'] if x['case_id']==c['case_id'])
        endpoints={}
        for z in ('Ca','La'):
            t=next(t for t in pm['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
            a=c['matrix'][z]['origin'];b=old['matrix'][z]['origin'];x=xyz(verify(a['xyz']));y=xyz(verify(b['xyz']))
            proposal=read_json(run/'proposals/proposals'/t['task_id']/'result.json')
            actual=xyz(verify(proposal['proposal']['coordinate']));old_candidate=next(q for q in old['candidates'] if q['id']=='adaptive_'+z)
            archive=xyz(verify(old_candidate['xyz']))
            endpoints[z]={'source_atom_order_and_coordinates_exact':x==y,'charge':t['charge'],'multiplicity':t['multiplicity'],
                'fresh_origin_components':a['components'],'archived_origin_components':b['components'],
                'fresh_minus_archived_origin_kcal_mol':relative_components(a['components'],b['components']),
                'active_mode_ids':t['active_mode_ids'],'proposal_status':proposal['status'],
                'proposal_max_coordinate_difference_A':float(np.max(np.abs(np.array([a[1:] for a in actual])-np.array([a[1:] for a in archive])))),
                'optimizer':proposal['optimizer'],'boundary_flag':proposal['boundary_flag']}
        old_min=pool.choose_rows(old['matrix'],['origin','adaptive_Ca','adaptive_La']);new=c['pool']
        rows.append({'case_id':c['case_id'],'preparation_seconds':src['source_preparation_seconds'],'endpoints':endpoints,
            'old_minimal_pool':old_min,'fresh_pool':new,
            'fresh_minus_archived_operational_R_model_kcal_mol':new['operational']['composite_R_model_kcal_mol']-old_min['operational']['composite_R_model_kcal_mol'],
            'reported':next(r for r in result['rows'] if r['case_id']==c['case_id'])})
    summaries=[read_json(p) for p in run.glob('*/gpu_summary_*.json')]
    low=[]
    for stage in ('origins','pool'):
        lp=run/stage/'solvent/shard_0/manifest.json';m=read_json(lp)
        for t in m['tasks']:
            rp=Path(t['output_path']+'.execution.json');receipt=read_json(rp)
            low.append({'task_id':t['task_id'],'receipt':record(rp),'returncode':receipt['returncode'],'normal_termination':receipt['normal_termination'],'scf_converged':receipt['scf_converged']})
    report={'fresh_result':record(run/'result.json'),'fresh_preparation':record(run/'source_preparation/preparation.json'),'old_pool':record(original_path),
        'rows':rows,'denominator':2,'available':result['available'],
        'MACE_calls':sum(x['new_MACE_calls'] for x in summaries),'MACE_failed_requests':sum(x['failed_requests'] for x in summaries),
        'GPU_summaries':[record(p) for p in run.glob('*/gpu_summary_*.json')],
        'GFN2_attempts':len(low),'GFN2_complete':sum(r['returncode']==0 and r['normal_termination'] and r['scf_converged'] for r in low),
        'low_receipts':low,'DFT_calls':0,'fresh_preparation_seconds':sum(x['source_preparation_seconds'] for x in prep['cases']),
        'max_abs_origin_composite_difference_kcal_mol':max(abs(e['fresh_minus_archived_origin_kcal_mol']['composite_kcal_mol']) for r in rows for e in r['endpoints'].values()),
        'max_abs_pooled_R_difference_kcal_mol':max(abs(r['fresh_minus_archived_operational_R_model_kcal_mol']) for r in rows),
        'production_changed':False,'independent_biological_validation':False}
    write_new(args.output,report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','low_receipts','GPU_summaries')},indent=2))

if __name__=='__main__':main()
