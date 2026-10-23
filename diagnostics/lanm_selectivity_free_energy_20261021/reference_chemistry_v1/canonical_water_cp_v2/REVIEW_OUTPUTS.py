#!/usr/bin/env python3
"""Saved native evidence review. Does not modify admission or make native calls."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import sys
import time

F = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'

def pin(p):
    p=Path(p); return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    tick=time.monotonic(); ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); ap.add_argument('results',nargs='+'); a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False); rows=[]
    for filename in a.results:
        r=json.loads(Path(filename).read_text()); output=Path(r['output']['path']); assert pin(output)==r['output']; text=output.read_text(); native=text.split('****END OF INPUT****')[-1]
        residuals=[{'name':name.strip(),'value':float(v),'tolerance':float(t),'passes_printed_tolerance':abs(float(v))<=float(t)} for name,v,t in re.findall(r'Last\s+(.+?)\s+\.{3}\s+('+F+r')\s+Tolerance\s*:\s*('+F+r')',native)]
        row={'cell_id':r['cell_id'],'result':pin(filename),'output':pin(output),'native_status':r['status'],'all_native_checks_pass':all(r['checks'].values()),'printed_scf_criteria':residuals,'all_six_printed_criteria_pass':len(residuals)==6 and all(x['passes_printed_tolerance'] for x in residuals),'physical_state_qualified':False,'method_accuracy_qualified':False,'elapsed_seconds':r['elapsed_seconds'],'allocated_cpu_seconds':r['allocated_cpu_seconds']}
        if 'COUPLED CLUSTER ENERGY' in native:
            trace=native.split('Iter       E(tot)',1)[1].split('COUPLED CLUSTER ENERGY',1)[0]
            records=re.findall(r'^\s*('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s*$',trace,re.M)
            row['cc_trace']=[dict(zip(('total_Eh','correlation_Eh','delta_Eh','residual','seconds'),map(float,v))) for v in records]
            row['CC_final_residual']=row['cc_trace'][-1]['residual']; row['CC_residual_passes_declared_1e_8']=row['CC_final_residual']<=1e-8
            row['T1_diagnostic']=[float(v) for v in re.findall(r'T1 diagnostic\s*\.{3}\s*('+F+r')',native)]
            row['triples_Eh']=[float(v) for v in re.findall(r'^Triples Correction \(T\)\s*\.{3}\s*('+F+r')',native,re.M)]
            amps=native.split('LARGEST AMPLITUDES',1)[1].split('!!!',1)[0]
            row['largest_printed_amplitudes']=[float(v) for v in re.findall(r'^\s*\d+->\s*\d+\s+\d+->\s*\d+\s+('+F+r')\s*$',amps,re.M)]
            row['correlated_electrons']=r['correlated_electrons']; row['effective_nuclear_charges']=[x['effective_Z'] for x in r['atomic_table']]
        rows.append(row)
    result={'status':'saved_native_review','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rows':rows,'molecular_calls':0,'local_threads':1,'command':sys.argv,'source':pin(__file__),'elapsed_seconds':time.monotonic()-tick,'interpretation':'Native admission is preserved. Printed SCF residual compliance, coupled-cluster convergence, physical state, and method accuracy are distinct.'}
    (a.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'native_admitted':sum(x['all_native_checks_pass'] for x in rows),'all_scf_criteria_pass':sum(x['all_six_printed_criteria_pass'] for x in rows),'result':pin(a.out/'RESULT.json')}))

if __name__=='__main__': main()
