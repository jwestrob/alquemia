"""Read actual strict20 traces, charge agreement and allocation receipts."""
import json,re,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
W=ROOT/'workspaces/strict_native_stopping_20260923/run_v1';c=read_json(W/'COLLECTION.json');traces=[];total=0.
for r in c['rows']:
    text=verify(r['actual']['output']).read_text();block=text.split('Iteration    Energy (Eh)')[-1].split('**** Energy Check signals convergence ****')[0];rows=[]
    for line in block.splitlines():
        fields=line.split()
        if len(fields)==7 and fields[0].isdigit():
            try:rows.append({'cycle':int(fields[0]),'energy_hartree':float(fields[1]),'delta_E_printed':float(fields[2]),'RMSDP':float(fields[3]),'MaxDP':float(fields[4]),'Omega':float(fields[5]),'time_seconds':float(fields[6])})
            except ValueError:pass
    if not rows or rows[-1]['cycle']!=r['details']['cycles']:raise RuntimeError('trace incomplete')
    match=re.search(r'TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec',text)
    d,h,m,s,ms=map(int,match.groups());elapsed=d*86400+h*3600+m*60+s+ms/1000;total+=elapsed
    traces.append({k:r[k] for k in ('task_id','case_id','candidate','metal','medium','seed_kind')}|{'actual':r['actual'],'observed_TolE_hartree':r['observed_TolE_hartree'],'convergence':{k:r['details'][k] for k in ('energy','max_density','rms_density','cycles')},'Mulliken_charges_e':r['audit']['mulliken_charges_e'],'trace':rows,'ORCA_seconds':elapsed})
changes=[]
for cell in c['cells']:
    rr={r['seed_kind']:r for r in c['rows'] if all(r[k]==cell[k] for k in ('case_id','candidate','metal','medium'))}
    delta=np.array(rr['pass2']['audit']['mulliken_charges_e'])-np.array(rr['cold']['audit']['mulliken_charges_e'])
    changes.append({k:cell[k] for k in ('case_id','candidate','metal','medium')}|{'maximum_atomic_charge_difference_e':float(np.max(np.abs(delta))),'metal_charge_difference_e':float(delta[0]),'charge_L2_difference_e':float(np.linalg.norm(delta))})
write_new(W/'TRACE_AUDIT_v1.json',{'collection':record(W/'COLLECTION.json'),'rows':traces,'charge_agreement':changes,'new_molecular_calls':0,'implementation':record(__file__)})
raw=subprocess.check_output(['sacct','-j','1211229','--parsable2','--noheader','--format=JobID,State,ElapsedRaw,AllocCPUS,ReqMem,AllocTRES,MaxRSS'],text=True)
r=next(x.split('|') for x in raw.splitlines() if x.split('|')[0]=='1211229')
if r[1]!='COMPLETED':raise RuntimeError('job not complete')
cost={'job_id':'1211229','state':r[1],'wall_seconds':int(r[2]),'allocated_cpus':int(r[3]),'allocated_core_seconds':int(r[2])*int(r[3]),'requested_memory':r[4],'GPU_seconds':0,'actual_new_scalar_calls':20,'summed_ORCA_seconds':total,'sacct':raw,'new_MACE_DFT_optimization_calls':0,'historical_seed_cost_included':False,'local_test_report_metered':False}
write_new(W/'COSTS.json',cost)
summary={'collection':record(W/'COLLECTION.json'),'trace_audit':record(W/'TRACE_AUDIT_v1.json'),'cost':record(W/'COSTS.json'),'complete':c['complete'],'denominator':20,'two_start_agreement':c['seed_agreement_passes'],'two_start_denominator':10,'maximum_seed_energy_difference_kcal_mol':max(abs(r['seed_difference_kcal_mol']) for r in c['cells']),'maximum_seed_atomic_charge_difference_e':max(r['maximum_atomic_charge_difference_e'] for r in changes),'density_checks':{k:sum(r['details'][k]['within_printed_tolerance'] for r in c['rows']) for k in ('energy','max_density','rms_density')},'production_changed':False,'classifier_or_reference':None}
write_new(ROOT/'diagnostics/strict_native_stopping_20260923/RESULT.json',summary);print(json.dumps(summary,indent=2));print(json.dumps(cost,indent=2))
