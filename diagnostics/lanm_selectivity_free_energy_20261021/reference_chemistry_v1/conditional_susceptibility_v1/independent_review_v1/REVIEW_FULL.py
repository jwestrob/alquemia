#!/usr/bin/env python3
"""Extend the frozen scout's saved replay to all16 fields and all36 response gates."""
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('prior_scout_review',HERE/'REVIEW_SCOUT.py')
scout=importlib.util.module_from_spec(s);s.loader.exec_module(s)
ROOT=Path(__file__).resolve().parents[5]
D=scout.D;pin=scout.pin;read=scout.read;verify=scout.verify
BRANCH=scout.BRANCH

def main():
    tick=time.monotonic();cpu=time.process_time()
    path=ROOT/'workspaces'/BRANCH/'full_collection_v1/RESULT.json'
    collection=json.loads(path.read_text())
    out=ROOT/'workspaces'/BRANCH/'independent_full_actual_review_v1';out.mkdir(exist_ok=False)
    report={'status':'failed_conditional_susceptibility_full_review','collection':pin(path),'review_source':pin(__file__),'prior_reviewer':pin(HERE/'REVIEW_SCOUT.py'),'new_native_calls':0,'error':None,'checks':[]}
    def check(ok,label):
        report['checks'].append({'name':label,'passed':bool(ok)})
        if not ok:raise ValueError(label)
    try:
        source=collection['implementation'];panel=read(collection['panel']);report.update(implementation=source,panel=collection['panel'])
        check(source['sha256']=='4741a45c6cbb7d30e5c543f87b6716db8f57d51ff829329e6d92ac4e0c632a19','frozen reference source')
        spec=importlib.util.spec_from_file_location('independent_full_field_adapter',verify(source))
        adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
        ids=[c['cell_id'] for c in panel['cells']]
        check(collection['stage']=='full' and collection['declared']==16 and collection['cell_ids']==ids,'exact16 declared fields')
        check(collection['complete']==len(collection['rows'])==16 and set(collection['rows'])==set(ids),'complete16 denominator')
        rows={};report['raw_rows']={};native_cost=0.
        for c in panel['cells']:
            cid=c['cell_id'];saved=collection['rows'][cid]
            actual=adapter.admit_result(verify(saved['result']),c,panel,collection['panel']['path'])
            check(actual==saved,cid+' exact input/native/property saved replay')
            energy,mu=scout.property_tokens(verify(saved['properties']['property_file']))
            check(all(x.is_finite() for x in [energy,*mu]),cid+' finite raw property tokens')
            check(float(energy)==saved['properties']['energy_hartree'] and list(map(float,mu))==saved['properties']['dipole_au'],cid+' independent token match')
            rows.setdefault(c['zero_id'],{})[c['field_tag']]={'E':energy,'mu':mu}
            q=saved['parsed'];native_cost+=saved['allocated_elapsed_CPU_seconds_proxy']
            report['raw_rows'][cid]={'result':saved['result'],'input':saved['input'],'output':saved['output'],'property_file':saved['properties']['property_file'],'job_id':saved['job_id'],'field_au':c['field_au'],'energy_hartree':float(energy),'dipole_au':list(map(float,mu)),'printed_scf_criteria':q.get('printed_scf_criteria'),'native_scf_stop_signals':q.get('native_scf_stop_signals'),'energy_stage_metadata':q.get('energy_stage_metadata'),'final_COSX_correction_hartree':q.get('final_COSX_correction_hartree'),'total_energy_after_final_integration_hartree':q.get('total_energy_after_final_integration_hartree'),'native_elapsed_seconds':saved['elapsed_seconds'],'allocated_CPU_seconds_proxy':saved['allocated_elapsed_CPU_seconds_proxy']}
        check(set(rows)==set(adapter.BASES) and all(len(x)==4 for x in rows.values()),'four complete state quartets')
        states={};allgates=[];maxdiff=0.
        for gid in adapter.BASES:
            zero=panel['zero_replays'][gid];text=verify(panel['zero_field_records'][gid]['output']).read_text()
            tokens=re.findall(r'^Total Energy\s*:\s*('+scout.NUMBER+r')\s+Eh',text,re.M)
            check(len(tokens)==1 and float(tokens[0])==zero['energy_hartree'],gid+' exact same-stage zero energy')
            e0=D(tokens[0]);mu0=D(str(zero['dipole_au'][2]));values={};gates={}
            for step,h,minus,plus in [('h',D('.001'),'minus_h','plus_h'),('half',D('.0005'),'minus_half','plus_half')]:
                p,n=rows[gid][plus],rows[gid][minus]
                ce=-(p['E']+n['E']-2*e0)/(h*h);cm=[(a-b)/(2*h) for a,b in zip(p['mu'],n['mu'])];me=-(p['E']-n['E'])/(2*h)
                values[step]={'chi_energy':float(ce),'chi_moment':float(cm[2]),'chi_moment_column':list(map(float,cm)),'odd_energy_moment':float(me),'zero_moment':float(mu0),'odd_minus_zero':float(me-mu0),'even_moment_numerator':[float(a+b-2*D(str(z))) for a,b,z in zip(p['mu'],n['mu'],zero['dipole_au'])]}
                stored=collection['response']['per_state'][gid]['estimates'][step]
                differences=[abs(values[step][k]-stored[k]) for k in ['chi_energy','chi_moment','odd_energy_moment']]
                maxdiff=max(maxdiff,*differences);check(max(differences)<=1e-6,gid+step+' independent decimal estimates')
                gates['energy_vs_moment_'+step]=abs(float(ce-cm[2]))<=.5;gates['odd_slope_'+step]=abs(float(me-mu0))<=1e-4
            gates.update(chi_energy_step=abs(values['h']['chi_energy']-values['half']['chi_energy'])<=.5,chi_moment_step=abs(values['h']['chi_moment']-values['half']['chi_moment'])<=.5)
            check(gates==collection['response']['per_state'][gid]['checks'],gid+' all six gate classifications')
            allgates.extend(gates.values());states[gid]={'estimates':values,'checks':gates,'passed':all(gates.values())}
        contrasts={}
        ld,lt,dd,dt=adapter.BASES
        for name,terms in [('La',{lt:1,ld:-1}),('Dy',{dt:1,dd:-1}),('Dy_minus_La',{dt:1,dd:-1,lt:-1,ld:1})]:
            values={};gates={};stored=collection['response']['paired'] if name=='Dy_minus_La' else collection['response']['increments'][name]
            for step in ['h','half']:
                values[step]={key:math.fsum(c*states[g]['estimates'][step][key] for g,c in terms.items()) for key in ['chi_energy','chi_moment','odd_energy_moment','zero_moment','odd_minus_zero']}
                values[step]['chi_moment_column']=[math.fsum(c*states[g]['estimates'][step]['chi_moment_column'][j] for g,c in terms.items()) for j in range(3)]
                differences=[abs(values[step][k]-stored['estimates'][step][k]) for k in ['chi_energy','chi_moment','odd_energy_moment','zero_moment']]
                maxdiff=max(maxdiff,*differences);check(max(differences)<=2e-6,name+step+' independent signed contrast')
                gates['energy_vs_moment_'+step]=abs(values[step]['chi_energy']-values[step]['chi_moment'])<=1.
            gates.update(chi_energy_step=abs(values['h']['chi_energy']-values['half']['chi_energy'])<=1.,chi_moment_step=abs(values['h']['chi_moment']-values['half']['chi_moment'])<=1.)
            check(gates==stored['checks'],name+' all four gate classifications');allgates.extend(gates.values());contrasts[name]={'terms':terms,'estimates':values,'checks':gates,'passed':all(gates.values())}
        check(len(allgates)==36 and all(allgates)==collection['numerical_checks_passed'],'all36 numerical classifications agree')
        report.update(status='passed_conditional_susceptibility_full_review' if all(allgates) else 'numerically_unqualified_conditional_susceptibility_full_review',states=states,contrasts=contrasts,numerical_checks_passed=all(allgates),numerical_checks_count=36,max_independent_arithmetic_difference=maxdiff,admitted_native_programs=16,allocated_native_CPU_seconds_proxy=native_cost,collector_cost=collection['cost'])
        report['limitations']='Fixed geometries, one lab-z response column and two field sizes only; conditional increment is not alpha3 or full-shellB. Failed printed SCF criteria and field-dependent final-COSX stages remain visible. Passing numerical consistency is not an absolute reference-error bound, physical-state qualification or unique-mechanism evidence. No model/reference table or native calls generated here.'
    except Exception as exc:report['error']=repr(exc)
    report['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1,'new_native_calls':0}
    dest=out/'REVIEW.json';dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'review':pin(dest),'status':report['status'],'error':report['error'],'contrasts':report.get('contrasts')}))
    return 0 if report['error'] is None else 1

if __name__=='__main__':raise SystemExit(main())
