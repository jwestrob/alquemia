#!/usr/bin/env python3
"""One terminal scout review: raw readback and cost only; no native invocation."""
import argparse
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

N=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(path):return json.loads(Path(path).read_text())
def verified(p):
    if pin(p['path'])!=p:raise ValueError('Changed authority: '+p['path'])
    return Path(p['path'])
def number(x):return float(x.replace('D','E').replace('d','e'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    tick=time.monotonic();cpu=time.process_time();cfg=read(args.config);checks=[]
    def check(name,passed,detail=None):
        checks.append({'name':name,'passed':bool(passed),'detail':detail})
        if not passed:raise ValueError(name)
    for key in ('implementation','panel','ready','preflight','release','submission','recovery_inputs','recovery_preparation'):verified(cfg[key])
    check('frozen reviewer',cfg['reviewer']==pin(__file__))
    path=Path(cfg['result_path']);check('terminal native RESULT exists',path.is_file())
    result_pin=pin(path);saved=read(path)
    check('terminal producer record',saved.get('finished_utc') and saved.get('returncode') is not None)
    spec=importlib.util.spec_from_file_location('frozen_rotation_adapter',cfg['implementation']['path']);adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
    p=adapter.load(cfg['panel']['path']);cell=next(c for c in p['cells'] if c['cell_id']==adapter.SCOUT)
    out={'status':'failed_dry_acetate_correlated_rotation_scout_review','implementation':cfg['implementation'],'panel':cfg['panel'],'result':result_pin,'reviewer':pin(__file__),'configuration':pin(args.config),'checks':checks,'cost_reviewed':False,'new_native_calls':0,'remaining_seven_released':False,'error':None}
    try:
        check('exact scout and release',saved['cell_id']==adapter.SCOUT and saved['release']==cfg['release'])
        admitted=adapter.admit(p,path)
        check('unchanged full raw native admission',admitted==saved)
        raw=verified(saved['output']).read_text();native=raw.split('****END OF INPUT****')[-1];start=read(verified(saved['start']))
        check('normal native termination',saved['returncode']==0 and native.count('ORCA TERMINATED NORMALLY')==1)
        check('native RHF and NoFrozenCore',bool(re.search(r'Reference Wavefunction\s*\.{2,}\s*RHF\b',native)) and bool(re.search(r'Frozen core treatment\s*\.{2,}\s*NO frozen core\b',native)))
        for name,pattern,expected in [('electrons',r'Number of [Ee]lectrons\s+(?:NEL\s+)?\.{2,}\s*(\d+)',40),('correlated electrons',r'Number of correlated electrons\s*\.{2,}\s*(\d+)',40),('orbital basis',r'Number of basis functions\s*\.{2,}\s*(\d+)',189),('charge',r'Total Charge\s+(?:Charge\s+)?\.{2,}\s*([-+]?\d+)',2),('multiplicity',r'Multiplicity\s+(?:Mult\s+)?\.{2,}\s*(\d+)',1)]:
            values=[int(x) for x in re.findall(pattern,native,re.I)];check('native '+name,bool(values) and all(x==expected for x in values),values)
        check('La ECP46 assigned only to real metal',re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+\s+ECP group\s*=>',native)==['0'] and re.findall(r'replacing\s+(\d+)\s+core electrons',native)==['46'])
        basis=native.split('BASIS SET INFORMATION',1)[1].split('ECP PARAMETER INFORMATION',1)[0]
        groups={int(i):sum(int(n)*(2*'spdfghi'.index(l)+1) for n,l in re.findall(r'(\d+)([spdfghi])',s)) for i,s in re.findall(r'Group\s+(\d+)\s+Type\s+[A-Za-z]+:?\s*:\s*[^\n]*?contracted to\s+([0-9spdfghi]+)',basis,re.I)}
        assignment={int(i):int(g) for i,g in re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+:?\s+basis set group\s*=>\s*(\d+)',basis)}
        counts=[groups.get(assignment.get(i)) for i in range(8)];check('all eight orbital centers',counts==[47,31,31,31,31,6,6,6],counts)
        au=native.split('CARTESIAN COORDINATES (A.U.)',1)[1].split('INTERNAL COORDINATES',1)[0]
        atoms=re.findall(r'^\s*(\d+)\s+([A-Za-z]+)\s+('+N+r')\*?\s+\d+\s+'+N+r'\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*$',au,re.M)
        check('all eight real nuclei',len(atoms)==8 and [int(a[0]) for a in atoms]==list(range(8)) and [a[1] for a in atoms]==cell['symbols'] and [number(a[2]) for a in atoms]==[11,6,6,8,8,1,1,1])
        error=max(abs(number(a[j+3])*adapter.b.BOHR-cell['positions_A'][i][j]) for i,a in enumerate(atoms) for j in range(3));check('literal rotation coordinates',error<1e-6,error)
        header=re.findall(r'^Iter\s+E\(tot\)\s+E\(Corr\)\s+Delta-E\s+Residual\s+Time\s*$',native,re.M);check('recognized canonical residual table',len(header)==1)
        trace=native.split(header[0],1)[1].split('COUPLED CLUSTER ENERGY',1)[0]
        records=[list(map(number,x)) for x in re.findall(r'^\s*('+N+r')\s+('+N+r')\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*$',trace,re.M)]
        tol=[number(x) for x in re.findall(r'Convergence tolerance \(max\. residuum\)\s*\.{2,}\s*('+N+r')',native)]
        check('unchanged residual threshold and actual convergence',tol==[1e-8] and bool(records) and abs(records[-1][3])<=1e-8 and bool(re.search(r'---\s*The\s+(?:Coupled[ -]Cluster|CCSD)\s+iterations\s+have\s+converged\s*---',native,re.I)),{'threshold':tol,'last_row':records[-1] if records else None})
        tokens={}
        for key,pattern in [('final',r'FINAL SINGLE POINT ENERGY\s+'),('CCSD_T',r'^E\(CCSD\(T\)\)\s*\.{2,}\s*'),('CCSD',r'^E\(CCSD\)\s*\.{2,}\s*'),('triples',r'^Triples Correction \(T\)\s*\.{2,}\s*')]:
            values=re.findall(pattern+'('+N+')',native,re.M);check('unique finite '+key,len(values)==1 and math.isfinite(number(values[0])));tokens[key]=values[0]
        E={k:Decimal(v) for k,v in tokens.items()};closures={'final_minus_CCSD_T':str(E['final']-E['CCSD_T']),'final_minus_CCSD_minus_triples':str(E['final']-E['CCSD']-E['triples'])}
        check('independent final and component closure',all(abs(Decimal(v))<Decimal('1e-8') for v in closures.values()),closures)
        check('no D4 or local-correlation input',all(x not in verified(saved['input']).read_text() for x in ['D4','DLPNO','RIJCOSX','EnGrad']))
        check('one call, four ranks, one node and allocation-derived MaxCore',saved['new_native_calls']==1 and start['ranks']==4 and start['nodes']==1 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4))
        check('finite measured allocation cost',math.isfinite(saved['elapsed_seconds']) and saved['elapsed_seconds']>0 and saved['allocated_cpu_seconds']==4*saved['elapsed_seconds'])
        check('empty native stderr',verified(saved['stderr']).stat().st_size==0)
        old=next(x for x in p['cost_controls'] if x['cell_id']=='cc__acetate_bi__La__complex')
        out.update(status='passed_dry_acetate_correlated_rotation_scout_review',cost_reviewed=True,raw_output=saved['output'],raw_input=saved['input'],start=saved['start'],job_id=saved['job_id'],energies_hartree=tokens,energy_closures_hartree=closures,final_CC_residual=records[-1][3],CC_table_rows=len(records),T1_diagnostic=saved['T1_diagnostic'],largest_printed_amplitudes=saved['largest_printed_amplitudes'],SCF_stop_signals=saved['native_scf_stop_signals'],SCF_check_mode=saved['native_scf_mode'],printed_SCF_criteria=saved['printed_scf_criteria'],resource_cost={'MPI':4,'nodes':1,'memory_MiB':start['memory_MiB'],'MaxCore_MiB':start['maxcore_MiB'],'native_elapsed_seconds':saved['elapsed_seconds'],'native_allocation_CPU_seconds':saved['allocated_cpu_seconds'],'prior_La_bidentate_elapsed_seconds':old['elapsed_seconds'],'scout_to_prior_elapsed_ratio':saved['elapsed_seconds']/old['elapsed_seconds'],'scheduler_allocation_CPU_seconds':None,'scheduler_evidence':'not yet joined; native elapsed×ranks is not CPU utilization'},interpretation={'complete_native_cells':1,'declared_total_cells':8,'derivative_checks_available':0,'directional_or_method_conclusion':False,'physical_state_qualified':False,'method_accuracy_qualified':False,'limitation':'One exact native complex energy only. No two-step derivative, CP force, physical-state, diffuse/basis or method-accuracy claim. Root retains cost decision and remaining-seven release.'},review_role='Original Singer-authored raw scientific checks executed by Planck through explicit V2 JSON-representation recovery; no native rerun. Root owns independent source review and continuation release.')
    except Exception as exc:out['error']=repr(exc)
    out['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    args.out.mkdir(parents=True,exist_ok=False);rp=args.out/'REVIEW.json'
    with rp.open('x') as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    (args.out/'REPORT.md').write_text('# Canonical rotation scout raw review\n\n'+out['status']+'.\n\n'+('One of eight native endpoints admitted. The three full-panel derivative checks remain unavailable. Native cost and all SCF/CC warnings are retained in REVIEW.json. No remaining endpoint was released by this review.\n' if out['error'] is None else 'Failure preserved: '+out['error']+'\n'))
    print(json.dumps({'review':pin(rp),'status':out['status'],'cost':out['review_cost'],'new_native_calls':0},indent=2))

if __name__=='__main__':main()
