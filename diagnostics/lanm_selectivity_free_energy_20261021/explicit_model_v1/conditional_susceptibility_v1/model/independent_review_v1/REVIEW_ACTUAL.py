#!/usr/bin/env python3
"""One saved replay of19 model outputs and independent decimal response arithmetic."""
import decimal
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

ROOT=Path(__file__).resolve().parents[6]
D=decimal.Decimal
decimal.getcontext().prec=50
HA=D('627.509474063'); BOHR=D('0.529177210544'); DEBYE=D('4.803206802')
BRANCH='lanm_selectivity_free_energy_20261021/explicit_model_v1/conditional_susceptibility_v1/model'

def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(ok,label):
        checks.append({'name':label,'passed':bool(ok)})
        if not ok:raise ValueError(label)
    result_path=ROOT/'workspaces'/BRANCH/'run_v1/execution/RESULT.json'
    r=json.loads(result_path.read_text());cfg=r['config'];manifest=json.loads(Path(r['manifest']['path']).read_text())
    check(pin(cfg['implementation']['path'])==cfg['implementation'],'frozen implementation')
    check(pin(r['manifest']['path'])==r['manifest'],'frozen preparation')
    spec=importlib.util.spec_from_file_location('independent_actual_field_adapter',cfg['implementation']['path'])
    adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
    admitted=adapter.admit_result(cfg,manifest,r)
    check(admitted['status']=='saved_conditional_model_admitted' and admitted['new_native_calls']==0,'all19 raw admissions')
    check(r['job_id']=='18473' and r['native_program_invocations']==len(r['rows'])==19,'actual finite native ledger')
    raw={};max_ind=0.;max_closure=0.;native_cost=0.
    for jid,rec in r['rows'].items():
        text=Path(rec['moments']['analyze_output']['path']).read_text()
        tokens=re.findall(r'Total Potential Energy\s*:\s*([-+0-9.]+)',text)
        check(len(tokens)==1,jid+' raw energy token')
        lines=Path(rec['moments']['induced_output']['path']).read_text().splitlines()[1:]
        vectors=[[D(t) for t in line.split()[2:]] for line in lines if line.strip()]
        check(len(vectors)==manifest['rows'][rec['job']['state']]['atoms'] and all(len(v)==3 for v in vectors),jid+' raw site denominator')
        raw[jid]={'E':D(tokens[0]),'mu_ind':[sum(v[j] for v in vectors) for j in range(3)]}
        max_ind=max(max_ind,rec['moments']['induction']['maximum_printed_RMS_residual_Debye'])
        max_closure=max(max_closure,abs(rec['moments']['component_closure_kcal_mol']))
        check(rec['moments']['induction']['solutions']==2,jid+' both induction solutions')
        check(len(rec['receipts'])==1,jid+' one native program')
        native_cost+=rec['receipts'][0]['child_user_cpu_seconds']+rec['receipts'][0]['child_system_cpu_seconds']
    controls={}
    z=r['rows']['zero_Dy_trimer']['moments'];old=manifest['rows']['Dy_trimer']['baseline']
    controls['zero_identity']={'energy_exact':z['energy_kcal_mol']==old['energy_kcal_mol'],'induced_exact':z['induced_Debye']==old['induced_Debye']}
    check(all(controls['zero_identity'].values()),'zero-field printed identity exact')
    for sign,label in [(1,'ion_La_plus_h'),(-1,'ion_La_minus_h')]:
        field=D(str(manifest['jobs'][label]['field_MV_cm'][2]))/D('1439.96455')
        expected_e=-D('.5')*D('332.0637133')*D('1.6597')*field**2
        expected_u=D('4.803206802')*D('1.6597')*field
        errors={'E_kcal_mol':float(raw[label]['E']-expected_e),'u_Debye':float(raw[label]['mu_ind'][2]-expected_u)}
        check(abs(errors['E_kcal_mol'])<=1e-7 and abs(errors['u_Debye'])<=1e-6 and raw[label]['mu_ind'][:2]==[0,0],label+' independent analytic sign/scale')
        check(sign*raw[label]['mu_ind'][2]>0 and raw[label]['E']<0,label+' nonzero field applied')
        controls[label]=errors
    states={};max_replay=0.;max_em=0.;max_odd=0.;gates=[]
    for label in ('La_dimer','La_trimer','Dy_dimer','Dy_trimer'):
        zero=manifest['rows'][label]['baseline'];e0=D(str(zero['energy_kcal_mol']))
        mu0=D(str(zero['total_precise_lab_au'][2]));values={}
        for step,h in [('h',D('.001')),('half',D('.0005'))]:
            plus=raw[label+'__plus_'+step];minus=raw[label+'__minus_'+step]
            ce=-(plus['E']+minus['E']-2*e0)/(HA*h*h)
            cm=[(a-b)/(2*h*DEBYE*BOHR) for a,b in zip(plus['mu_ind'],minus['mu_ind'])]
            me=-(plus['E']-minus['E'])/(2*h*HA)
            values[step]={'chi_E':float(ce),'chi_mu_column':list(map(float,cm)),'mu_E':float(me),'mu_E_minus_mu0':float(me-mu0)}
            saved=r['summary']['states'][label]['steps'][step]
            differences=[abs(float(ce)-saved['chi_E_a0_cubed']),abs(float(me)-saved['mu_E_au'])]+[abs(float(x)-y) for x,y in zip(cm,saved['chi_mu_column_a0_cubed'])]
            max_replay=max(max_replay,*differences);check(max(differences)<=1e-8,label+step+' independent raw-token estimates')
            em=abs(float(ce-cm[2]));odd=abs(float(me-mu0));max_em=max(max_em,em);max_odd=max(max_odd,odd)
            gates.extend([em<=.5,odd<=1e-4])
        gates.extend([abs(values['h']['chi_E']-values['half']['chi_E'])<=.5,abs(values['h']['chi_mu_column'][2]-values['half']['chi_mu_column'][2])<=.5])
        states[label]=values
    increments={}
    for name,terms in {'La':{'La_trimer':1,'La_dimer':-1},'Dy':{'Dy_trimer':1,'Dy_dimer':-1},'Dy_minus_La':{'Dy_trimer':1,'Dy_dimer':-1,'La_trimer':-1,'La_dimer':1}}.items():
        values={}
        for step in ('h','half'):
            values[step]={'chi_E':math.fsum(c*states[k][step]['chi_E'] for k,c in terms.items()),'chi_mu_column':[math.fsum(c*states[k][step]['chi_mu_column'][j] for k,c in terms.items()) for j in range(3)],'mu_E':math.fsum(c*states[k][step]['mu_E'] for k,c in terms.items()),'mu_E_minus_mu0':math.fsum(c*states[k][step]['mu_E_minus_mu0'] for k,c in terms.items())}
            saved=r['summary']['increments'][name]['steps'][step]
            check(abs(values[step]['chi_E']-saved['chi_E_a0_cubed'])<=1e-8 and max(abs(a-b) for a,b in zip(values[step]['chi_mu_column'],saved['chi_mu_column_a0_cubed']))<=1e-8,name+step+' independent signed contrast')
            gates.append(abs(values[step]['chi_E']-values[step]['chi_mu_column'][2])<=1.)
        gates.extend([abs(values['h']['chi_E']-values['half']['chi_E'])<=1.,abs(values['h']['chi_mu_column'][2]-values['half']['chi_mu_column'][2])<=1.])
        increments[name]={'terms':terms,'steps':values}
    check(len(gates)==36 and all(gates),'all36 independent response gates')
    check(math.isfinite(max_ind) and 0<=max_ind<=1e-8,'all induction residuals within frozen gate')
    reported_native=r['cost']['child_user_CPU_seconds']+r['cost']['child_system_CPU_seconds']
    check(abs(native_cost-reported_native)<=1e-8,'sequential native CPU sum agrees with process-wide delta')
    output={'status':'passed_independent_actual_model_susceptibility_review','implementation':r['implementation'],'manifest':r['manifest'],'input_result':pin(result_path),'release':r['release'],'review_source':pin(__file__),'job_id':r['job_id'],'full_native_denominator':19,'initial_controls':controls,'response_checks_passed':36,'response_checks_declared':36,'checks':checks,'states':states,'increments':increments,'max_source_arithmetic_difference':max_replay,'max_state_energy_vs_moment_chi_difference_a0_cubed':max_em,'max_state_odd_slope_difference_ea0':max_odd,'maximum_printed_induction_residual_Debye':max_ind,'maximum_energy_component_closure_kcal_mol':max_closure,'actual_native_cost':r['cost'],'independent_sum_native_CPU_seconds':native_cost,'interpretation':'Numerical consistency of original baseline at four fixed systems and two z-field steps only. Quantum scout results are not inspected; no model/reference accuracy or alpha3 claim. Print bounds exclude solver and physical-model error. Prototype external-field guard remains unchanged.','review_cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1,'new_native_calls':0,'synthetic_control_reruns':0}}
    out=ROOT/'workspaces'/BRANCH/'independent_actual_review_v1';out.mkdir(exist_ok=False)
    path=out/'REVIEW.json';path.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'review':pin(path),'max_energy_moment':max_em,'max_odd':max_odd,'maximum_induction':max_ind,'increments':increments,'cost':output['review_cost']}))

if __name__=='__main__':main()
