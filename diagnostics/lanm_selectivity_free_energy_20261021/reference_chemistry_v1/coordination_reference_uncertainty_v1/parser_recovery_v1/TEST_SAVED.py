#!/usr/bin/env python3
"""Focused mutations of an actual saved scout; never invokes molecular code."""
import argparse,importlib.util,json,time
from pathlib import Path
s=importlib.util.spec_from_file_location('recollection',Path(__file__).with_name('RECOLLECT.py'))
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--result',required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();tick,cpu=time.monotonic(),time.process_time();p=m.load(a.panel);r=json.loads(Path(a.result).read_text());raw=m.b.verify(r['output']).read_text();checks=[]
    def check(value,name):m.b.check(value,name);checks.append(name)
    old=m.old.parse(raw,r['cell'],0);q=m.parse(raw,r['cell'],0)
    check(q['status']=='complete','actual native trace admitted')
    check({k for k,v in old['checks'].items() if not v}=={'local_TCutMKN','local_TCutDO','local_MP2_converged'},'exact original three grammar failures retained')
    check(all(q['checks'][k]==v for k,v in old['checks'].items() if k not in m.CHANGED),'all unrelated gates unchanged')
    check(q['energy_hartree']==r['energy_hartree'] and q['CC_terminal_residuals']==r['CC_terminal_residuals'] and q['triples_terminal_deltas']==r['triples_terminal_deltas'],'energy and CC/triples unchanged')
    e=q['stage_specific_local_evidence'];check([x['kind'] for x in e['stages']]==['semi-local MP2','full local MP2'],'two distinct MP2 stages retained')
    check({'strong Triples','weak Triples'}<={x['stage'] for x in e['map_threshold_reports']},'separate triples thresholds retained')
    cases=[
        ('TCutMKN:           1.000e-03','TCutMKN:           2.000e-03','local_TCutMKN','changed primary threshold'),
        ('    TCutDO             =        5.000e-03','    TCutDO             =        6.000e-03','local_TCutDO','changed full-local threshold'),
        ('   TCutMKN   ...     1.0e-03','   TCutMKN   ...     2.0e-03','local_TCutMKN','changed fine threshold'),
        ('TScaleMKNMP2PreScr:1.000e+01','TScaleMKNMP2PreScr:2.000e+01','local_TCutMKN','changed prescreen scale'),
        ('    TCutMKN            =        1.000e-02','    TCutMKN            =        2.000e-02','local_TCutMKN','changed coarse threshold'),
        ('.... full local MP2','.... unfamiliar local MP2','local_MP2_converged','unknown stage'),
        (' => CONVERGED',' => FINISHED','local_MP2_converged','missing positive convergence'),
        ('RMAX=   2.3e-07 => CONVERGED','RMAX=   9.0e-07 => CONVERGED','local_MP2_converged','excess terminal residual'),
        ('DE=   4.3e-09 RMAX','DE=   9.0e-07 RMAX','local_MP2_converged','excess terminal energy delta'),
        ('RMAX=   2.3e-07 => CONVERGED','RMAX=   ******** => CONVERGED','local_MP2_converged','starred terminal residual'),
        ('RMAX=   2.3e-07 => CONVERGED','RMAX=   NaN => CONVERGED','local_MP2_converged','nonfinite terminal residual'),
        ('RMAX=   7.4e-06','RMAX=   7.4e-06 => CONVERGED','local_MP2_converged','multiple convergence markers'),
        ('LMP2-Iter=   6:','LMP2-Iter=   9:','local_MP2_converged','nonsequential trace'),
        ('Residual convergence tolerance ...    5.0e-07','Residual convergence tolerance ...    5.0e-06','local_MP2_converged','relaxed printed residual tolerance'),
        ('Energy convergence tolerance   ...    1.0e-07','Energy convergence tolerance   ...    1.0e-06','local_MP2_converged','relaxed printed energy tolerance'),
    ]
    for before,after,key,name in cases:
        check(before in raw,'mutation exists: '+name)
        check(not m.stage_evidence(raw.replace(before,after),r['cell']['tier'])['checks'][key],'reject '+name)
    check(m.parse(raw,r['cell'],1)['status']=='failed_admission','nonzero exit remains rejected')
    admitted=m.admit_original(a.result,r['cell'],p,a.panel)
    check(admitted['original_status']=='failed_admission' and admitted['new_native_calls']==0,'exact original producer and failed result admitted without overwrite')
    a.out.mkdir(parents=True,exist_ok=False)
    result=dict(status='passed_saved_parser_recovery_controls',implementation=m.b.pin(m.__file__),test=m.b.pin(__file__),declaration=m.b.pin(m.HERE/'DECLARATION.md'),
        panel=m.b.pin(a.panel),original_result=m.b.pin(a.result),raw_output=r['output'],checks=checks,count=len(checks),new_native_calls=0,local_threads=1,
        wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    pin=m.b.write(a.out/'RESULT.json',result);print(json.dumps(dict(result=pin,status=result['status'],checks=len(checks))))
if __name__=='__main__':main()
