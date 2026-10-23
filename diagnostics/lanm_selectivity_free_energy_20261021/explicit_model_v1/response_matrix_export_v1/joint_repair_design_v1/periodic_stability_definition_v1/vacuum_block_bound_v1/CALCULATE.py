#!/usr/bin/env python3
"""Two fixed block-Gershgorin evaluations on six admitted vacuum matrices."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    assert pin(p['path'])==p,'Changed pin '+p['path']
    return Path(p['path'])

def read(p):return json.loads(verify(p).read_text())

def partitions(parsed,xyz):
    atoms={}
    for line in xyz.read_text().splitlines()[1:]:
        t=line.split()
        if t:atoms[int(t[0])]={'type':int(t[5]),'bonds':list(map(int,t[6:]))}
    sites=parsed['sites'];assert len(atoms)==len(sites)
    position={s['atom']:i for i,s in enumerate(sites)}
    assert set(position)==set(atoms) and all(atoms[s['atom']]['type']==s['type'] for s in sites)
    metal=[s['atom'] for s in sites if s['type'] in (500,501)]
    assert len(metal)==1 and not atoms[metal[0]]['bonds']
    waters=[]
    for s in sites:
        if s['type']==36:
            oxygen=s['atom'];hydrogens=sorted(atoms[oxygen]['bonds'])
            assert len(hydrogens)==2 and all(atoms[h]['type']==37 and atoms[h]['bonds']==[oxygen] for h in hydrogens)
            waters.append([oxygen]+hydrogens)
    covered=metal+[a for w in waters for a in w]
    assert len(covered)==len(set(covered))==len(atoms) and set(covered)==set(atoms)
    def expand(group):return [3*position[a]+c for a in group for c in range(3)]
    return {
        'site_3': [{'atoms':[s['atom']],'dofs':expand([s['atom']])} for s in sites],
        'intact_water_9_plus_metal_3': [{'atoms':metal,'dofs':expand(metal)}]+[{'atoms':w,'dofs':expand(w)} for w in waters]
    },len(waters)

def block_bound(S,partition):
    count=len(partition);norms=np.zeros((count,count));diagonal=[]
    for i,b in enumerate(partition):
        ev=np.linalg.eigvalsh(S[np.ix_(b['dofs'],b['dofs'])]);diagonal.append([float(ev[0]),float(ev[-1])])
        for j in range(i):
            sigma=np.linalg.svd(S[np.ix_(b['dofs'],partition[j]['dofs'])],compute_uv=False)
            norms[i,j]=norms[j,i]=float(sigma[0])
    rows=[]
    for i,b in enumerate(partition):
        radius=math.fsum(float(x) for x in norms[i]);lo,hi=diagonal[i]
        rows.append({'atoms':b['atoms'],'dofs_zero_based':b['dofs'],'diagonal_lambda_min':lo,'diagonal_lambda_max':hi,'off_block_radius':radius,'lower':lo-radius,'upper':hi+radius})
    lower=min(r['lower'] for r in rows);upper=max(r['upper'] for r in rows)
    return {'blocks':rows,'off_block_spectral_norms':norms.tolist(),'L':lower,'U':upper,'ratio_L_over_U':lower/upper if upper>0 else None,'criterion_margin_L_minus_1e_minus6_U':lower-1e-6*upper,'criterion_passed':bool(upper>0 and lower>1e-6*upper)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time();assert all(os.environ.get(k)=='1' for k in THREADS),'One-thread context'
    sources=json.loads((HERE/'INPUTS.json').read_text());result=read(sources['qualification_result']);review=read(sources['qualification_review']);manifest=read(sources['qualification_manifest'])
    assert review['status']=='passed_local_response_native_qualification' and review['result']==sources['qualification_result'] and result['numerical_checks_passed'] and result['native_program_invocations']==57
    assert result['manifest']==sources['qualification_manifest']
    spec=importlib.util.spec_from_file_location('frozen_block_bound_parser',verify(sources['matrix_parser']));helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    assert len(sources['snapshots'])==6 and len({x['export']['path'] for x in sources['snapshots']})==6
    rows={};numeric_checks=[]
    for item in sources['snapshots']:
        cid=item['id'];rec=result['rows'][cid];row=manifest['rows'][rec['job']['cell']]
        assert cid.startswith('matrix__') and item['export']==rec['exports'][0] and Path(item['export']['path']).name=='response_matrix_001.dat'
        assert item['xyz']==rec['runtime_input']['xyz'] and rec['status']=='complete'
        parsed=helper.parse_export(item['export'],row,manifest['units']);saved=rec['matrix_qualification']['exports'][0]
        assert parsed['sequence']==1 and all(parsed[k]==saved[k] for k in ('A','vectors','sites','checks','restored','peek'))
        A=np.array(parsed['A']);alpha=np.repeat([s['alpha_A3'] for s in parsed['sites']],3);sqrt=np.sqrt(alpha);B=sqrt[:,None]*A*sqrt[None,:];S=.5*(B+B.T)
        actual=np.array(saved['spectrum']['normalized_A']);dense=np.linalg.eigvalsh(S);roundoff=helper.arithmetic(len(B),B,actual,dense)
        difference=float(np.max(abs(actual-dense)));assert difference<=roundoff,'Saved dense replay outside inherited arithmetic'
        groups,waters=partitions(parsed,verify(item['xyz']));bounds={}
        for name,group in groups.items():
            value=block_bound(S,group)
            value.update(lower_gap_to_saved_min=float(actual[0]-value['L']),upper_gap_above_saved_max=float(value['U']-actual[-1]),floating_enclosure_passed=bool(value['L']<=actual[0]+roundoff and value['U']>=actual[-1]-roundoff))
            assert value['floating_enclosure_passed'],'Floating bound does not enclose saved spectrum'
            bounds[name]=value
        rows[cid]={'source':item['export'],'geometry':item['xyz'],'profile':rec['job']['profile'],'metal':row['metal'],'waters':waters,'sites':len(parsed['sites']),'dimension':len(B),'normalization':'sqrt(alpha_A3) A_A^-3 sqrt(alpha_A3)','antisymmetric_max':float(np.max(abs(B-B.T))),'symmetrization_operator_norm':float(np.linalg.norm(.5*(B-B.T),2)),'saved_dense_min':float(actual[0]),'saved_dense_max':float(actual[-1]),'saved_dense_ratio':float(actual[0]/actual[-1]),'saved_dense_criterion_passed':bool(actual[0]>1e-6*actual[-1]>0),'symmetrized_dense_min':float(dense[0]),'symmetrized_dense_max':float(dense[-1]),'saved_dense_max_replay_difference':difference,'inherited_arithmetic_bound':roundoff,'partitions':bounds}
    summary={name:{'passed':sum(r['partitions'][name]['criterion_passed'] for r in rows.values()),'total':6,'failed_ids':[cid for cid,r in rows.items() if not r['partitions'][name]['criterion_passed']]} for name in ('site_3','intact_water_9_plus_metal_3')}
    output={'status':'completed_exploratory_vacuum_block_bounds','implementation':pin(__file__),'declaration':pin(HERE/'DECLARATION.md'),'inputs':pin(HERE/'INPUTS.json'),'rows':rows,'summary':summary,'unique_snapshots':6,'partitions_per_snapshot':2,'criterion':'L > 1e-6 U with U > 0; unchanged','new_native_calls':0,'new_fitted_parameters':0,'rigorous_rounded_certificate':False,'periodic_qualification':False,'limitations':'Ordinary floating-point block bounds on six fixed vacuum matrices. A failed sufficient lower bound is not instability; successful bounds are not rigorous rounded certificates or periodic qualification.','cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    (args.out/'RESULT.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    lines=['# Fixed vacuum block bounds','', 'Ordinary floating-point exploration; no rigorous rounded certificate or periodic result. Both fixed partitions and the 1e-6 criterion were declared before calculation.','', '| System | Waters | Saved dense min/max | Site L / U | Site pass | Water+metal L / U | Water+metal pass |','|---|---:|---:|---:|:---:|---:|:---:|']
    for cid,r in rows.items():
        p,q=r['partitions']['site_3'],r['partitions']['intact_water_9_plus_metal_3']
        lines.append(f"| {cid.removeprefix('matrix__')} | {r['waters']} | {r['saved_dense_min']:.8g} / {r['saved_dense_max']:.8g} | {p['L']:.8g} / {p['U']:.8g} | {p['criterion_passed']} | {q['L']:.8g} / {q['U']:.8g} | {q['criterion_passed']} |")
    lines+=['',f"Site partition: {summary['site_3']['passed']}/6 pass. Intact-water partition: {summary['intact_water_9_plus_metal_3']['passed']}/6 pass.",'','All six saved dense operators retain their original positive-definiteness criterion. Negative block lower bounds are failures of this sufficient estimate, not detected instabilities. No partition or weight search was performed.', '',f"Cost: {output['cost']['wall_seconds']:.3f} wall seconds; {output['cost']['CPU_seconds']:.3f} CPU seconds; one thread; zero native calls."]
    (args.out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':pin(args.out/'RESULT.json'),'summary':summary,'cost':output['cost']}));print('\n'.join(lines))

if __name__=='__main__':main()
