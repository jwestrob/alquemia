#!/usr/bin/env python3
"""Saved-only review of the fixed amide8 collector; no jobs or native calls."""
import argparse
from decimal import Decimal, getcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[4]
PROJECT='lanm_selectivity_free_energy_20261021'
D=ROOT/'diagnostics'/PROJECT/'reference_chemistry_v1'
W=ROOT/'workspaces'/PROJECT/'reference_chemistry_v1'
SOURCE=D/'hydrated_amide_force_qualification_v1/REFERENCE.py'
BASE=W/'hydrated_amide_force_qualification_v1'
PANEL=BASE/'prepared_v2/PANEL.json'
COLLECTION=BASE/'collection_v1/RESULT.json'
SUBMISSIONS=BASE/'run_v1/SUBMISSIONS.json'
EXPECTED_SOURCE='7de9c17aa9198ba07c66c3f934390fc2a67d94a77860eca8ba4c0378b25ad677'
EXPECTED_PANEL='38b42c1020c47b12512076d59d66e299d4c8eb7f14754f8cfe5f71bec502e4cc'
NUMBER=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
TERMS={'electronic':'electronic_energy_hartree','D4':'D4_energy_hartree','PBE0_D4':'energy_hartree'}
HA=Decimal('627.5094740631')
checks=[]
def check(ok,label):
    checks.append({'check':label,'passed':bool(ok)})
    if not ok:raise ValueError(label)
def pin(path):
    path=Path(path).resolve();return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def decimal(token):return Decimal(token.replace('D','E').replace('d','e'))
def forbid(*args,**kwargs):raise ValueError('Native or scheduler execution forbidden in review')
def close(actual,expected,label):check(abs(float(actual)-float(expected))<=1e-7,label+' (floating/decimal replay only)')

def execute(out):
    tick=time.monotonic();out.mkdir(parents=True,exist_ok=False);getcontext().prec=40
    result={'status':'started','implementation':pin(__file__),'checks':checks,'new_native_calls':0,'new_submissions':0,'local_threads':1,'command':sys.argv}
    try:
        check(pin(SOURCE)['sha256']==EXPECTED_SOURCE and pin(PANEL)['sha256']==EXPECTED_PANEL,'Frozen adapter/panel')
        sp=importlib.util.spec_from_file_location('amide8_actual',SOURCE);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);m.subprocess.run=forbid
        global HA
        HA=Decimal(str(m.b.HA))
        p=m.load(PANEL);submitted=json.loads(SUBMISSIONS.read_text());saved=json.loads(COLLECTION.read_text())
        check(len(submitted)==8 and {x['cell_id'] for x in submitted}=={c['cell_id'] for c in p['cells']},'Exact submitted eight-slot inventory')
        paths=[x['result'] for x in submitted];replay=m.collect(PANEL,paths)
        check(replay==saved,'Original collector exact saved replay')
        result.update(collector=pin(COLLECTION),panel=pin(PANEL),submissions=pin(SUBMISSIONS),collector_status=saved['status'],slots=saved['slots'])
        if saved['status']!='complete':
            result.update(status='unavailable_native_denominator',reason='No derivative conclusion: all eight native slots must be admitted.')
            return result
        values={};native=[];costs=[]
        for entry in submitted:
            filename=Path(entry['result']);r=json.loads(filename.read_text());cid=entry['cell_id'];c=next(x for x in p['cells'] if x['cell_id']==cid)
            start=json.loads((filename.parent/'START.json').read_text())
            check(str(r['job_id'])==str(entry['job_id'])==str(start['job_id']) and start['ranks']==4,'Job/cell/allocation '+cid)
            check(r['new_native_calls']==1 and r['allocated_cpu_seconds']==r['elapsed_seconds']*4,'One call/native allocation accounting '+cid)
            raw=m.b.verify(r['output']).read_text();text=raw.split('****END OF INPUT****')[-1]
            token_lists={
                'electronic':re.findall(r'^Total Energy\s*:\s*('+NUMBER+r')\s+Eh',text,re.M),
                'D4':re.findall(r'^Dispersion correction\s+('+NUMBER+r')',text,re.M),
                'PBE0_D4':re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+r')',text)}
            check(all(token_lists.values()),'Raw native energy tokens '+cid)
            values[cid]={k:decimal(v[-1]) for k,v in token_lists.items()}
            for term,field in TERMS.items():check(float(values[cid][term])==r[field],'Raw/stored '+term+' '+cid)
            check(abs(values[cid]['electronic']+values[cid]['D4']-values[cid]['PBE0_D4'])<Decimal('1e-8'),'Unchanged electronic/D4 closure '+cid)
            strict=m.parse(raw,c,r['returncode'])
            check(strict['status']=='complete' and all(strict['checks'].values()),'Actual raw admission '+cid)
            check(strict['native_energy_stages']==r['native_energy_stages'],'Retained stage evidence '+cid)
            m.b.verify(r['stderr'])
            native.append({'result':pin(filename),'input':r['input'],'output':r['output'],'stderr':r['stderr'],'raw_energy_tokens':token_lists,'checks':strict['checks'],'printed_scf_criteria':strict['printed_scf_criteria'],'SCF_mode':strict['native_scf_mode'],'SCF_stop_signals':strict['native_scf_stop_signals'],'energy_stages':strict['native_energy_stages']})
            costs.append({'cell_id':cid,'job_id':entry['job_id'],'native_elapsed_seconds':r['elapsed_seconds'],'allocated_cpus':4,'native_elapsed_times_allocated_cpus':r['allocated_cpu_seconds']})
        derivatives={}
        for metal in ('La','Dy'):
            derivatives[metal]={}
            for h,suffix in ((Decimal('.010'),'010'),(Decimal('.005'),'005')):
                a=values[f'pbe0__{metal}__plus_{suffix}__complex'];z=values[f'pbe0__{metal}__minus_{suffix}__complex']
                derivatives[metal][str(float(h))]={term:(a[term]-z[term])*HA/(2*h) for term in TERMS}
        derivatives['Dy_minus_La']={h:{term:derivatives['Dy'][h][term]-derivatives['La'][h][term] for term in TERMS} for h in ('0.01','0.005')}
        summaries={}
        for name,steps in derivatives.items():
            original=saved['Dy_minus_La'] if name=='Dy_minus_La' else saved['per_metal'][name]
            if name=='Dy_minus_La':analytic=saved['Dy_minus_La']['central_analytic_derivative_kcal_mol_A']
            else:
                center=p['central'][name];g=center['original_cell']['geometry'];s=json.loads(m.b.verify(center['saved_result']).read_text())
                analytic=-sum(s['forces_kcal_mol_A'][i][j]*g['radial_direction'][j] for i in g['ligand_indices'] for j in range(3))
                close(analytic,original['central_analytic']['coordinate_derivative_kcal_mol_A'],name+' actual central ligand-force projection')
            errors={h:float(v['PBE0_D4'])-analytic for h,v in steps.items()}
            changes={term:float(steps['0.005'][term]-steps['0.01'][term]) for term in TERMS}
            for h,vals in steps.items():
                for term,v in vals.items():close(v,original['derivatives_kcal_mol_A'][h][term],name+'/'+h+'/'+term)
                close(errors[h],original['PBE0_D4_FD_minus_central_analytic'][h],name+'/'+h+' FD-minus-analytic')
            for term,change in changes.items():
                close(change,original['smaller_minus_larger_step'][term],name+'/'+term+' two-step difference')
                check((abs(change)<=.1)==original['two_step_resolved'][term],name+'/'+term+' declared two-step gate')
            check(all(abs(x)<=.1 for x in errors.values())==original['sampled_analytic_consistency_at_both_steps'],name+' declared analytic gate')
            summaries[name]={'central_analytic_kcal_mol_A':analytic,'derivatives_kcal_mol_A':{h:{k:float(v) for k,v in vals.items()} for h,vals in steps.items()},'smaller_minus_larger_step':changes,'two_step_resolved':{k:abs(v)<=.1 for k,v in changes.items()},'FD_minus_analytic':errors,'both_analytic_checks_pass':all(abs(x)<=.1 for x in errors.values())}
        changed=[{'metal':metal,'shift':shift,'comparison':v} for metal in ('La','Dy') for shift,v in saved['per_metal'][metal]['energy_stage_comparisons'].items() if v['changed_energy_stage_signature']]
        scheduler={};receipts=[]
        for receipt in [BASE/'run_v1/completion_event.json',BASE/'collection_v1/completion_event.json']+[Path(x['result']).parents[1]/'completion_event.json' for x in submitted]:
            if receipt.is_file():
                event=json.loads(receipt.read_text());receipts.append(pin(receipt))
                for jid,row in event.get('records',{}).items():scheduler[jid]=row
        jobids={str(x['job_id']) for x in submitted};scheduled={k:v for k,v in scheduler.items() if k in jobids}
        complete_cost=all(k in scheduled and scheduled[k].get('state')=='COMPLETED' for k in jobids)
        result.update(status='passed_saved_result_review',native_reviews=native,independent_derivatives=summaries,changed_energy_stage_signatures=changed,coarse_screen_preserved=saved['coarse_screen_preserved'],native_cost_rows=costs,native_elapsed_times_allocated_cpus=sum(x['native_elapsed_times_allocated_cpus'] for x in costs),scheduler_receipts=receipts,scheduler_records=scheduled,scheduler_complete_denominator=complete_cost,scheduler_allocated_cpu_seconds=sum(int(v['elapsed_seconds'])*int(v['allocated_cpus']) for v in scheduled.values()) if complete_cost else None,cost_measurement_limit='Allocation wall time products are not measured CPU utilization.',physical_state_qualified=False,physical_accuracy_qualified=False,derivative_resolution_kcal_mol_A=.1)
        return result
    except Exception as exc:
        result.update(status='failed_saved_review',error=repr(exc));raise
    finally:
        result.update(passed_checks=sum(x['passed'] for x in checks),elapsed_seconds=time.monotonic()-tick)
        (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();r=execute(a.out)
    print(json.dumps({'status':r['status'],'checks':r['passed_checks'],'result':pin(a.out/'RESULT.json')}))
