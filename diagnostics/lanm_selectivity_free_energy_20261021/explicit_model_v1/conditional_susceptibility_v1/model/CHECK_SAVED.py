#!/usr/bin/env python3
"""Targeted no-native controls for the finite original-baseline field adapter."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('conditional_model',HERE/'RUN.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
OUT=a.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/conditional_susceptibility_v1/model/saved_controls_v1'
PREP=OUT.parent/'preparation_v1'
checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def rejects(name,fn):
    try:fn()
    except (ValueError,KeyError):checks.append(name);return
    raise AssertionError('Did not reject: '+name)

def main():
    tick=time.monotonic();cpu=time.process_time();OUT.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((PREP/'CONFIG_RUN_TEMPLATE.json').read_text());m,h=a.context(cfg)
    check('exact19_and_three_first',len(m['jobs'])==19 and list(m['jobs'])[:3]==list(a.CONTROL_IDS))
    for jid,j in m['jobs'].items():
        text=a.verify(j['input']['key']).read_text();base=a.verify(m['rows'][j['state']]['input']['key']).read_text()
        check(jid+'_field_only_edit',text==a.render(base,j['field_au'][2]) and text.count('EXTERNAL-FIELD')==1)
        for kind in ('xyz','parameters'):check(jid+'_'+kind+'_unchanged',a.verify(j['input'][kind]).read_bytes()==a.verify(m['rows'][j['state']]['input'][kind]).read_bytes())
    for label,row in m['rows'].items():
        check(label+'_native_M_agreement',max(abs(x-y) for x,y in zip(row['baseline']['total_precise_lab_Debye'],row['baseline']['total_lab_Debye']))<=.0005+row['atoms']*.5e-6+2e-9)
        shifted=copy.deepcopy(row);translation=[.1,-.2,.3]
        for sig in shifted['frame_signatures']:
            sig['position_A']=[x+y for x,y in zip(sig['position_A'],translation)]
            sig['axis_positions_A']=[[x+y for x,y in zip(p,translation)] for p in sig['axis_positions_A']]
        change=[x-y for x,y in zip(a.permanent_moment(shifted),row['permanent_lab_Debye'])]
        check(label+'_charge_origin_translation',max(abs(x-a.DEBYE*row['charge']*t) for x,t in zip(change,translation))<1e-12)
    for field in (float('nan'),.002):rejects('undeclared_field_'+str(field),lambda f=field:a.render(a.verify(m['rows']['La_dimer']['input']['key']).read_text(),f))
    rejects('existing_field',lambda:a.render('EXTERNAL-FIELD 0 0 0\npolar-print\nsave-uinduce\n',.001))
    rejects('prototype_keyword',lambda:a.render('POLAR-PAIR-RESPONSE 500 36 .05 0\npolar-print\nsave-uinduce\n',.001))
    bad=copy.deepcopy(cfg);bad['maximum_native_programs']=20;rejects('20_program_config',lambda:a.context(bad))
    synthetic={}
    polar={'La_dimer':11.,'La_trimer':24.,'Dy_dimer':12.,'Dy_trimer':26.}
    for jid,job in m['jobs'].items():
        row=m['rows'][job['state']];value=copy.deepcopy(row['baseline']);f=job['field_au'][2]
        if job['state']=='ion_La':
            internal=f*a.FIELD_FACTOR/a.ELEFIELD
            value['energy_kcal_mol']=-.5*a.ELECTRIC*1.6597*internal**2
            value['induced_sum_Debye']=[0.,0.,a.DEBYE*1.6597*internal]
        elif f:
            chi=polar[job['state']];mu0=value['total_precise_lab_au'][2]
            value['energy_kcal_mol']+=a.HARTREE*(-mu0*f-.5*chi*f*f)
            value['total_precise_lab_au'][2]+=chi*f
        synthetic[jid]={'status':'complete','job':job,'moments':value}
        if jid in a.CONTROL_IDS:check(jid+'_analytic_control',a.control(job,row,value)['passed'])
    result=a.summarize(m,synthetic);check('linear_fixture_all_gates',result['numerical_checks_passed'])
    for state,expected in polar.items():
        for step,v in result['states'][state]['steps'].items():check(state+'_'+step+'_curvature',abs(v['chi_E_a0_cubed']-expected)<1e-6 and abs(v['chi_mu_column_a0_cubed'][2]-expected)<1e-9)
    for label,expected in [('La',13.),('Dy',14.),('Dy_minus_La',1.)]:
        for step,v in result['increments'][label]['steps'].items():check(label+'_'+step+'_balanced_increment',abs(v['chi_E_a0_cubed']-expected)<1e-6)
    missing=copy.deepcopy(synthetic);missing.pop(next(iter(missing)));rejects('missing_slot',lambda:a.summarize(m,missing))
    nonlinear=copy.deepcopy(synthetic);nonlinear['La_trimer__plus_half']['moments']['energy_kcal_mol']+=.001
    check('inconsistent_field_not_qualified',not a.summarize(m,nonlinear)['numerical_checks_passed'])
    ionjob=m['jobs']['ion_La_plus_h'];ion=m['rows']['ion_La'];good=synthetic['ion_La_plus_h']['moments']
    for name,change in [('field_noop',lambda v:v.update(energy_kcal_mol=0.,induced_sum_Debye=[0.,0.,0.])),('wrong_energy_factor',lambda v:v.update(energy_kcal_mol=2*v['energy_kcal_mol'])),('wrong_field_sign',lambda v:v.update(induced_sum_Debye=[-x for x in v['induced_sum_Debye']]))]:
        bad=copy.deepcopy(good);change(bad);rejects(name,lambda v=bad:a.control(ionjob,ion,v))
    # Exercise the actual serial runner with an explicitly fake invoke/parser.
    old_context,old_release,old_parse=a.context,a.release_gate,a.parse
    class FakeNative:
        def __init__(self):self.calls=[]
        def invoke(self,argv,cwd,output,receipts):
            self.calls.append(cwd.name);output.write_text('SYNTHETIC: no native process was started.\n')
            receipts.append({'argv':argv,'cwd':str(cwd),'returncode':0,'threads':{k:'1' for k in a.THREADS},'output':a.pin(output)})
    try:
        os.environ.update({k:'1' for k in a.THREADS});os.environ.update(SLURM_JOB_ID='SYNTHETIC_NO_NATIVE',SLURM_CPUS_PER_TASK='1',SLURM_NTASKS='1')
        for name,fail in [('complete',None),('control_failure','ion_La_plus_h')]:
            folder=OUT/name;folder.mkdir();native=FakeNative()
            a.context=lambda cfg:(m,(None,None,native));a.release_gate=lambda cfg,manifest:{'execution_root':str(folder)}
            def parsed(output,row,field,helpers):
                jid=output.parent.name;v=copy.deepcopy(synthetic[jid]['moments'])
                if jid==fail:v.update(energy_kcal_mol=0.,induced_sum_Debye=[0.,0.,0.])
                return v
            a.parse=parsed
            value=a.run(cfg,folder)
            if fail:
                check('control_failure_stops_before_fields',native.calls==list(a.CONTROL_IDS[:2]) and value['native_program_invocations']==2 and value['global_stop'] and all(value['rows'][k]['status']=='not_attempted' for k in list(m['jobs'])[2:]))
            else:check('serial_complete19_control_order',value['status']=='conditional_model_susceptibility_complete' and native.calls==list(m['jobs']) and value['native_program_invocations']==19 and value['controls_passed'])
    finally:a.context,a.release_gate,a.parse=old_context,old_release,old_parse
    receipt={'status':'passed_saved_controls','implementation':a.pin(HERE/'RUN.py'),'control_source':a.pin(__file__),'manifest':a.pin(PREP/'MANIFEST.json'),'checks':checks,'passed':len(checks),'new_native_calls':0,'synthetic_invocations_are_not_processes':True,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    a.save(OUT/'RESULT.json',receipt);print(json.dumps(a.pin(OUT/'RESULT.json')))

if __name__=='__main__':main()
