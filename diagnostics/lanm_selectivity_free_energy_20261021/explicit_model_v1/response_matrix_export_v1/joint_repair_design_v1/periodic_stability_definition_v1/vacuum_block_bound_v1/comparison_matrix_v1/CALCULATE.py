#!/usr/bin/env python3
"""Fixed block comparison spectra and one conditional Cw=1 solve each."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
PARTITIONS=('site_3','intact_water_9_plus_metal_3')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,'Changed pin '+p['path']
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time();assert all(os.environ.get(k)=='1' for k in THREADS),'One-thread context'
    inputs=json.loads((HERE/'INPUTS.json').read_text());old=read(inputs['unweighted_result']);origin=read(inputs['unweighted_inputs'])
    assert old['status']=='completed_exploratory_vacuum_block_bounds' and old['inputs']==inputs['unweighted_inputs'] and old['unique_snapshots']==6
    for p in (inputs['unweighted_declaration'],inputs['unweighted_source'],origin['qualification_result'],origin['qualification_review'],origin['qualification_manifest']):verify(p)
    assert set(old['rows'])==set(s['id'] for s in origin['snapshots'])
    for s in origin['snapshots']:
        verify(s['export']);verify(s['xyz']);assert old['rows'][s['id']]['source']==s['export'] and old['rows'][s['id']]['geometry']==s['xyz']
    rows={};spectra=0;solves=0
    for cid,row in old['rows'].items():
        out={'source':row['source'],'geometry':row['geometry'],'waters':row['waters'],'metal':row['metal'],'saved_dense_min':row['saved_dense_min'],'saved_dense_max':row['saved_dense_max'],'saved_dense_ratio':row['saved_dense_ratio'],'partitions':{}}
        assert set(row['partitions'])==set(PARTITIONS)
        for name in PARTITIONS:
            original=row['partitions'][name];blocks=original['blocks'];norms=np.array(original['off_block_spectral_norms']);n=len(blocks)
            assert norms.shape==(n,n) and np.array_equal(norms,norms.T) and np.array_equal(np.diag(norms),np.zeros(n)) and np.all(norms>=0)
            lo=np.array([b['diagonal_lambda_min'] for b in blocks]);hi=np.array([b['diagonal_lambda_max'] for b in blocks]);C=np.diag(lo)-norms
            spectrum=np.linalg.eigvalsh(C);spectra+=1;minimum=float(spectrum[0]);maximum=float(spectrum[-1])
            arithmetic=128*n*np.finfo(float).eps*max(1.,float(np.max(abs(C))),abs(row['saved_dense_min']))
            assert minimum<=row['saved_dense_min']+arithmetic,'Comparison spectral lower bound enclosure'
            rec={'blocks_atoms':[b['atoms'] for b in blocks],'comparison_matrix':C.tolist(),'comparison_spectrum':spectrum.tolist(),'comparison_lambda_min':minimum,'comparison_lambda_max':maximum,'original_unweighted_L':original['L'],'original_unweighted_U':original['U'],'original_unweighted_passed':original['criterion_passed'],'positive_weight_witness_possible_in_exact_comparison_matrix':bool(minimum>0),'weighted_criterion_passed':False,'weights':None,'L_weighted':None,'U_weighted':None,'no_weight_solve_reason':None,'inherited_arithmetic_enclosure_bound':arithmetic}
            if minimum>0:
                w=np.linalg.solve(C,np.ones(n));solves+=1;cw=np.array([math.fsum(float(C[i,j])*float(w[j]) for j in range(n)) for i in range(n)])
                finite_positive=bool(np.isfinite(w).all() and np.all(w>0));rec.update(weights=w.tolist(),all_weights_finite_positive=finite_positive,linear_solve_residual_max=float(np.max(abs(cw-1.))))
                if finite_positive:
                    lower_rows=[float(cw[i]/w[i]) for i in range(n)];upper_rows=[float(hi[i])+math.fsum(float(norms[i,j])*float(w[j])/float(w[i]) for j in range(n) if i!=j) for i in range(n)]
                    L,U=min(lower_rows),max(upper_rows)
                    rec.update(lower_per_block=lower_rows,upper_per_block=upper_rows,L_weighted=L,U_weighted=U,ratio_L_over_U=L/U if U>0 else None,criterion_margin=L-1e-6*U,weighted_criterion_passed=bool(U>0 and L>1e-6*U),floating_enclosure_passed=bool(L<=row['saved_dense_min']+arithmetic and U>=row['saved_dense_max']-arithmetic))
                    assert rec['floating_enclosure_passed'],'Weighted bound enclosure'
                else:rec['no_weight_solve_reason']='Single prescribed solve did not produce finite positive weights; no recovery.'
            else:rec['no_weight_solve_reason']='Comparison minimum is nonpositive; no positive block-weight Gershgorin witness can certify positive definiteness from this comparison matrix. No solve attempted.'
            out['partitions'][name]=rec
        rows[cid]=out
    assert spectra==12 and solves<=12
    summary={name:{'comparison_positive':sum(r['partitions'][name]['comparison_lambda_min']>0 for r in rows.values()),'weighted_criterion_passed':sum(r['partitions'][name]['weighted_criterion_passed'] for r in rows.values()),'total':6,'failed_ids':[cid for cid,r in rows.items() if not r['partitions'][name]['weighted_criterion_passed']]} for name in PARTITIONS}
    result={'status':'completed_exploratory_comparison_matrix_bounds','implementation':pin(__file__),'declaration':pin(HERE/'DECLARATION.md'),'inputs':pin(HERE/'INPUTS.json'),'rows':rows,'summary':summary,'small_comparison_spectra':spectra,'conditional_linear_solves':solves,'new_native_calls':0,'optimizer_calls':0,'rounded_certificate':False,'native_guard_qualified':False,'periodic_qualification':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    (args.out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Fixed comparison-matrix result','','Same six snapshots and two partitions; first unweighted attempt preserved. Ordinary floating-point exploration only.','','| System | Partition | min eig(C) | Weighted L | Weighted U | Frozen criterion |','|---|---|---:|---:|---:|:---:|']
    for cid,r in rows.items():
        for name,p in r['partitions'].items():
            L='not solved' if p['L_weighted'] is None else f"{p['L_weighted']:.9g}";U='not solved' if p['U_weighted'] is None else f"{p['U_weighted']:.9g}"
            lines.append(f"| {cid.removeprefix('matrix__')} | {name} | {p['comparison_lambda_min']:.9g} | {L} | {U} | {p['weighted_criterion_passed']} |")
    lines+=['',f'{spectra} comparison spectra; {solves} prescribed linear solves; zero optimization or native calls.','', 'Nonpositive comparison eigenvalues rule out positive-weight Gershgorin certificates for that compressed comparison matrix in exact arithmetic. They do not imply instability of the original operator. All numerical signs and solves here remain uncertified floating-point results.', '',f"Cost: {result['cost']['wall_seconds']:.3f} wall seconds / {result['cost']['CPU_seconds']:.3f} CPU seconds, one thread."]
    (args.out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':pin(args.out/'RESULT.json'),'summary':summary,'spectra':spectra,'solves':solves,'cost':result['cost']}));print('\n'.join(lines))

if __name__=='__main__':main()
