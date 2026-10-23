#!/usr/bin/env python3
"""One saved admission and independent Decimal raw-component arithmetic."""
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent; BASE=HERE.parent
ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
WORK=ROOT/'workspaces'/BASE.relative_to(ROOT/'diagnostics'); OUT=WORK/'independent_actual_review_v1'
RESULT=WORK/'run_v1/execution/RESULT.json'
checks={}; tick=time.monotonic(); cpu=time.process_time()
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(label,value):
    checks[label]=bool(value);assert value,label
def verify(p):
    check('pin/'+p['path'],pin(p['path'])==p);return Path(p['path'])
def dec(x):return Decimal(x.replace('D','E').replace('d','e'))
def hu(x):return Decimal(5).scaleb(x.as_tuple().exponent-1)
ready=json.loads((BASE/'READY.json').read_text()); source=verify(ready['implementation'])
sp=importlib.util.spec_from_file_location('admitted_component_producer',source);a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
r=a.admit_result(verify(ready['manifest']),RESULT)
check('existing_saved_admission',r['error'] is None and r['complete_cells']==r['new_native_calls']==16)
check('actual_job',r['job_id']=='18733')
m=json.loads(verify(r['manifest']).read_text())
source_review=json.loads((WORK/'independent_review_v1/REVIEW.json').read_text())
check('source_review_joins',source_review['implementation']==r['implementation'] and source_review['manifest']==r['manifest'])
oldcheck=json.loads(verify(source_review['independent_checks']).read_text())
raw={}; receipts=[]; max_energy_closure=Decimal(0)
number=r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?'
labels=('Bond Stretching','Angle Bending','Stretch-Bend','Out-of-Plane Bend','Torsional Angle','Van der Waals','Atomic Multipoles','Polarization','Total')
for cid,rr in r['rows'].items():
    q=rr['receipts'][0];receipts.append(q);text=verify(q['output']).read_text()
    matches=re.findall(r'^\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+number+r')\s+(\d+)\s*$',text.split('Energy Component Breakdown :')[1],re.M)
    tokens={name.strip():value for name,value,count in matches}; counts={name.strip():int(count) for name,value,count in matches}
    totals=re.findall(r'Total Potential Energy\s*:\s*('+number+')',text)
    check(cid+'/raw_denominator',len(totals)==1 and len(matches)==8 and set(tokens)==set(labels[:-1]));tokens['Total']=totals[0]
    check(cid+'/exact_tokens',tokens==rr['parsed']['tokens'])
    check(cid+'/interaction_counts',counts==m['central'][m['rows'][cid]['metal']]['energy']['interactions'])
    values={k:dec(v) for k,v in tokens.items()}; bounds={k:hu(v) for k,v in values.items()}; raw[cid]=(values,bounds)
    remainder=values['Total']-sum((values[k] for k in labels[:-1]),Decimal(0)); bound=sum(bounds.values(),Decimal(0))
    check(cid+'/raw_closure',abs(remainder)<=bound and remainder==dec(rr['parsed']['component_sum_remainder_kcal_mol']) and bound==dec(rr['parsed']['component_sum_rounding_bound_kcal_mol']))
    max_energy_closure=max(max_energy_closure,abs(remainder))
    st=json.loads((Path(q['cwd'])/'START.json').read_text())
    check(cid+'/start_lineage',st['input']==rr['input'] and st['manifest']==r['manifest'] and st['release']==r['release'] and st['analyze']==r['analyze'] and st['cell_id']==cid)
summary={}; physical_checks={}; errors=[]; changes=[]; max_fd_closure=Decimal(0)
with localcontext() as ctx:
    ctx.prec=50
    for direction in ('radial','rotation'):
        summary[direction]={}
        for kind in ('La','Dy','paired'):
            coefficients={'La':1} if kind=='La' else {'Dy':1} if kind=='Dy' else {'La':-1,'Dy':1}
            saved=r['comparison']['directions'][direction][kind]
            analytic=sum(coefficients[k]*oldcheck['central_raw_observations'][k]['derivatives_kcal_mol_A'][direction] for k in coefficients)
            check(direction+'/'+kind+'/independent_central_projection',abs(analytic-saved['central_analytic_kcal_mol_A'])<=1e-12)
            steps={}
            for step in ('0.01','0.005'):
                h=Decimal(step); ds={k:Decimal(0) for k in labels}; bs={k:Decimal(0) for k in labels}
                for metal,c in coefficients.items():
                    pair={s:next(cid for cid,row in m['rows'].items() if row['metal']==metal and row['direction']==direction and row['q_A']==s*float(h)) for s in (-1,1)}
                    for k in labels:
                        ds[k]+=c*(raw[pair[1]][0][k]-raw[pair[-1]][0][k])/(2*h)
                        bs[k]+=abs(c)*(raw[pair[1]][1][k]+raw[pair[-1]][1][k])/(2*h)
                sv=saved['steps'][step]; prefix=direction+'/'+kind+'/'+step
                for k in labels:
                    check(prefix+'/'+k+'/derivative',ds[k]==dec(sv['derivative_decimal_tokens'][k]) and float(ds[k])==sv['derivatives_kcal_mol_A'][k])
                    check(prefix+'/'+k+'/force_sign',-float(ds[k])==sv['generalized_forces_kcal_mol_A'][k])
                    check(prefix+'/'+k+'/print_bound',float(bs[k])==sv['print_bounds_kcal_mol_A'][k])
                closure=ds['Total']-sum((ds[k] for k in labels[:-1]),Decimal(0)); bound=sum(bs.values(),Decimal(0));max_fd_closure=max(max_fd_closure,abs(closure))
                check(prefix+'/component_closure',abs(closure)<=bound and float(closure)==sv['component_sum_remainder_kcal_mol_A'] and float(bound)==sv['component_sum_bound_kcal_mol_A'])
                error=float(ds['Total'])-analytic; errors.append(abs(error)); flag=abs(error)<=.1;label=f'{direction}/{kind}/h={float(h)}/analytic';physical_checks[label]=flag
                check(prefix+'/resolution_flag',flag==sv['analytic_resolution_passed'] and flag==r['comparison']['numerical_checks'][label])
                check(prefix+'/central_force_closure',abs(error-sv['total_minus_analytic_kcal_mol_A'])<=1e-12)
                steps[step]={'derivatives_kcal_mol_A':{k:str(v) for k,v in ds.items()},'absolute_print_bounds_kcal_mol_A':{k:str(v) for k,v in bs.items()},'total_minus_analytic_kcal_mol_A':error}
            diffs={k:dec(steps['0.005']['derivatives_kcal_mol_A'][k])-dec(steps['0.01']['derivatives_kcal_mol_A'][k]) for k in labels}
            for k,value in diffs.items():check(direction+'/'+kind+'/'+k+'/two_step',abs(float(value)-saved['smaller_minus_larger_step_kcal_mol_A'][k])<=1e-12)
            change=float(diffs['Total']);changes.append(abs(change));flag=abs(change)<=.1;label=f'{direction}/{kind}/two_step';physical_checks[label]=flag
            check(label+'/resolution_flag',flag==saved['two_step_total_resolution_passed'] and flag==r['comparison']['numerical_checks'][label])
            summary[direction][kind]={'central_analytic_kcal_mol_A':analytic,'steps':steps,'total_two_step_difference_kcal_mol_A':change}
check('all18_flags_replayed',physical_checks==r['comparison']['numerical_checks'] and len(physical_checks)==18 and all(physical_checks.values())==r['comparison']['numerical_checks_passed'])
check('native_elapsed_cost',abs(sum(q['elapsed_seconds'] for q in receipts)-r['cost']['native_elapsed_seconds'])<=1e-12)
for field,cost in [('child_user_cpu_seconds','child_user_CPU_seconds'),('child_system_cpu_seconds','child_system_CPU_seconds')]:
    check(cost,abs(sum(q[field] for q in receipts)-r['cost'][cost])<=1e-10)
check('serial_application_allocation',r['cost']['allocated_CPUs']==1 and r['cost']['allocated_CPU_seconds']==r['cost']['wall_seconds'])
event_text=(WORK/'run_v1/completion_event.json').read_text(); event=json.loads(event_text); scheduler=event['records']['18733']
check('scheduler_complete',scheduler['state']=='COMPLETED' and scheduler['exit_code']=='0:0' and scheduler['allocated_cpus']=='1')
OUT.mkdir(parents=True,exist_ok=True);(OUT/'COMPLETION_SNAPSHOT.json').write_text(event_text)
result={'status':'passed_dry_acetate_component_actual_review','implementation':r['implementation'],'manifest':r['manifest'],'result':pin(RESULT),'release':r['release'],'analyze':r['analyze'],'source_review':pin(WORK/'independent_review_v1/REVIEW.json'),'reviewer_source':pin(__file__),'complete_cells':16,'native_programs_in_producer':16,'new_native_calls_by_review':0,'existing_admission_calls':1,'check_count':len(checks),'checks':checks,'numerical_checks':physical_checks,'numerical_checks_passed':all(physical_checks.values()),'summary':summary,'max_abs_FD_minus_analytic_kcal_mol_A':max(errors),'max_abs_two_step_total_difference_kcal_mol_A':max(changes),'max_abs_raw_component_closure_kcal_mol':str(max_energy_closure),'max_abs_FD_component_closure_kcal_mol_A':str(max_fd_closure),'producer_cost':r['cost'],'scheduler':scheduler,'scheduler_allocated_CPU_seconds':int(scheduler['allocated_cpus'])*int(scheduler['elapsed_seconds']),'completion_snapshot':pin(OUT/'COMPLETION_SNAPSHOT.json'),'printed_induction_residual_available':False,'physical_accuracy_qualified':False,'component_causation_qualified':False,'review_cost':{'wall_seconds':time.monotonic()-tick,'process_CPU_seconds':time.process_time()-cpu}}
(OUT/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'review':pin(OUT/'REVIEW.json'),'checks':len(checks),'max_analytic_error':max(errors),'max_step_change':max(changes),'small_step':{d:{k:v['steps']['0.005']['derivatives_kcal_mol_A'] for k,v in kinds.items()} for d,kinds in summary.items()}},indent=2))
