#!/usr/bin/env python3
"""Saved-only four-field scout replay; no native execution or continuation release."""
import decimal
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

ROOT=Path(__file__).resolve().parents[5]
BRANCH='lanm_selectivity_free_energy_20261021/reference_chemistry_v1/conditional_susceptibility_v1'
D=decimal.Decimal
decimal.getcontext().prec=50
NUMBER=r'[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?'

def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    assert pin(p['path'])==p,'Changed pin '+p['path']
    return Path(p['path'])

def read(p):return json.loads(verify(p).read_text())

def property_tokens(path):
    text=Path(path).read_text()
    e=re.findall(r'&FinalEnergy\s+\[&Type "Double"\]\s+('+NUMBER+r')',text)
    assert len(e)==1,'Property energy denominator'
    tail=text.split('&dipoleTotal ',1)[1].split('\n   &',1)[0]
    a=re.findall(r'^\s*(\d+)\s+('+NUMBER+r')\s*$',tail,re.M)
    assert [int(i) for i,_ in a]==[0,1,2],'Property dipole denominator'
    return D(e[0]),[D(v) for _,v in a]

def main():
    tick=time.monotonic();cpu=time.process_time()
    collection_path=ROOT/'workspaces'/BRANCH/'scout_collection_v1/RESULT.json'
    collection=json.loads(collection_path.read_text())
    out=ROOT/'workspaces'/BRANCH/'independent_scout_actual_review_v1'
    out.mkdir(exist_ok=False)
    report={'status':'failed_conditional_susceptibility_scout_review','collection':pin(collection_path),'review_source':pin(__file__),'new_native_calls':0,'error':None,'checks':[]}
    def check(ok,label):
        report['checks'].append({'name':label,'passed':bool(ok)})
        if not ok:raise ValueError(label)
    try:
        source=collection['implementation'];panel=read(collection['panel'])
        report.update(implementation=source,panel=collection['panel'])
        check(source['sha256']=='4741a45c6cbb7d30e5c543f87b6716db8f57d51ff829329e6d92ac4e0c632a19','frozen reference source')
        spec=importlib.util.spec_from_file_location('independent_saved_field_scout',verify(source))
        adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
        check(collection['stage']=='scout' and collection['declared']==4 and collection['cell_ids']==panel['scout_cell_ids'],'exact four La-dimer fields')
        check(set(collection['rows'])==set(panel['scout_cell_ids']),'full four denominator')
        check(collection['complete']==4,'all four actual native cells admitted')
        report['raw_rows']={};rows={};native_cost=0.
        for cid in panel['scout_cell_ids']:
            saved=collection['rows'][cid];cell=next(c for c in panel['cells'] if c['cell_id']==cid)
            actual=adapter.admit_result(verify(saved['result']),cell,panel,collection['panel']['path'])
            check(actual==saved,cid+' exact saved input/native/property replay')
            e,mu=property_tokens(verify(saved['properties']['property_file']))
            check(all(v.is_finite() for v in [e,*mu]),cid+' finite independent raw property tokens')
            check(float(e)==saved['properties']['energy_hartree'] and list(map(float,mu))==saved['properties']['dipole_au'],cid+' independent property-token match')
            rows[cell['field_tag']]={'E':e,'mu':mu}
            parsed=saved['parsed'];native_cost+=saved['allocated_elapsed_CPU_seconds_proxy']
            report['raw_rows'][cid]={'result':saved['result'],'input':saved['input'],'output':saved['output'],'property_file':saved['properties']['property_file'],'job_id':saved['job_id'],'field_au':cell['field_au'],'energy_hartree':float(e),'dipole_au':list(map(float,mu)),'printed_scf_criteria':parsed.get('printed_scf_criteria'),'native_scf_stop_signals':parsed.get('native_scf_stop_signals'),'energy_stage_metadata':parsed.get('energy_stage_metadata'),'final_COSX_correction_hartree':parsed.get('final_COSX_correction_hartree'),'total_energy_after_final_integration_hartree':parsed.get('total_energy_after_final_integration_hartree'),'native_program_invocations':saved['native_program_invocations'],'native_elapsed_seconds':saved['elapsed_seconds'],'allocated_CPU_seconds_proxy':saved['allocated_elapsed_CPU_seconds_proxy']}
        zero_id=panel['cells'][0]['zero_id'];zero=panel['zero_replays'][zero_id]
        zero_raw=verify(panel['zero_field_records'][zero_id]['output']).read_text()
        tokens=re.findall(r'^Total Energy\s*:\s*('+NUMBER+r')\s+Eh',zero_raw,re.M)
        check(len(tokens)==1 and float(tokens[0])==zero['energy_hartree'],'same-stage saved zero raw energy token')
        e0=D(tokens[0]);mu0=D(str(zero['dipole_au'][2]));estimates={};gates={}
        for step,h,minus,plus in [('h',D('.001'),'minus_h','plus_h'),('half',D('.0005'),'minus_half','plus_half')]:
            p,n=rows[plus],rows[minus]
            ce=-(p['E']+n['E']-2*e0)/(h*h)
            cm=[(a-b)/(2*h) for a,b in zip(p['mu'],n['mu'])]
            me=-(p['E']-n['E'])/(2*h)
            estimates[step]={'chi_energy_a0_cubed':float(ce),'chi_moment_column_a0_cubed':list(map(float,cm)),'odd_energy_moment_ea0':float(me),'zero_moment_ea0':float(mu0),'odd_minus_zero_ea0':float(me-mu0),'even_moment_numerator_ea0':[float(a+b-2*D(str(z))) for a,b,z in zip(p['mu'],n['mu'],zero['dipole_au'])]}
            stored=collection['response']['per_state'][zero_id]['estimates'][step]
            check(max(abs(float(ce)-stored['chi_energy']),abs(float(cm[2])-stored['chi_moment']),abs(float(me)-stored['odd_energy_moment']))<=1e-6,step+' independent decimal response replay')
            gates['energy_vs_moment_'+step]=abs(float(ce-cm[2]))<=.5
            gates['odd_slope_'+step]=abs(float(me-mu0))<=1e-4
        gates['chi_energy_step']=abs(estimates['h']['chi_energy_a0_cubed']-estimates['half']['chi_energy_a0_cubed'])<=.5
        gates['chi_moment_step']=abs(estimates['h']['chi_moment_column_a0_cubed'][2]-estimates['half']['chi_moment_column_a0_cubed'][2])<=.5
        report.update(estimates=estimates,numerical_checks=gates,numerical_checks_passed=all(gates.values()),zero_reference=panel['zero_field_records'][zero_id]['result'],admitted_native_programs=4,allocated_native_CPU_seconds_proxy=native_cost,collector_cost=collection['cost'])
        check(gates==collection['response']['per_state'][zero_id]['checks'],'all six gate classifications replay')
        check(all(gates.values())==collection['numerical_checks_passed'],'collector numerical classification agrees')
        report['status']='passed_conditional_susceptibility_scout_review' if all(gates.values()) else 'numerically_unqualified_conditional_susceptibility_scout_review'
        report['limitations']='Four-field La-dimer native/numerical scout only. Preserves failed printed SCF criteria and final-COSX stages; passing finite-step consistency is not reference accuracy or physical-state qualification. An unqualified response leaves nonlinearity, energy/property consistency and numerical noise unresolved. No remaining12 release or native work performed.'
    except Exception as exc:
        report['error']=repr(exc)
    report['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1,'new_native_calls':0}
    p=out/'REVIEW.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'review':pin(p),'status':report['status'],'error':report['error'],'estimates':report.get('estimates')}))
    return 0 if report['error'] is None else 1

if __name__=='__main__':raise SystemExit(main())
