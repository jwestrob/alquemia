"""Narrow binder controls and mocked serial orchestration; zero molecular calls."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import time
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('regression_controls',HERE/'RUN.py');R=importlib.util.module_from_spec(s);s.loader.exec_module(R)
checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append({'name':name,'passed':True})
def reject(name,fn):
    try:fn()
    except (ValueError,KeyError):check(name,True);return
    raise AssertionError('Expected rejection '+name)
def fn_ast(path,name):
    return ast.dump(next(n for n in ast.parse(Path(path).read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name),include_attributes=False)

def main(manifest_path,out):
    wall,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in R.THREADS)
    m=json.loads(Path(manifest_path).read_text());inputs=json.loads((HERE/'SOURCE_INPUTS.json').read_text())
    check('literal_cycles_AST_unchanged',fn_ast(__file__.replace('TEST.py','RUN.py'),'cycles')==fn_ast(R.verify(inputs['previous_adapter']),'cycles'))
    old=R.read(R.read(inputs['previous_regression'])['manifest'])
    check('six_prepared_inputs_and_baselines_unchanged',m['rows']==old['rows'] and m['baseline_cells']==old['baseline_cells'])
    check('all_eight_exact_cycles_unchanged',m['cycles']==old['cycles'] and len(m['cycles'])==8)
    template=json.loads((Path(manifest_path).parent/'CONFIG_RUN_TEMPLATE.json').read_text())
    check('future_candidate_and_release_unbound',all(template[k] is None for k in ('evaluator_config','fit_result','fit_review','source_review','hostname','output_directory','root_release')))
    pin=R.pin(HERE/'RUN.py');params=[.1]*12;key=R.hashlib.sha256(json.dumps(params,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    cfg={'fit_result':pin,'evaluator_config':pin};mm={'inputs':{'evaluator':pin},'campaign_declaration':pin}
    fc={k:pin for k in ('implementation','adapter','core','declaration','execution_context','parent_evaluator','parent_configuration')}
    selected={'status':'complete','parameters':params}
    fr={'status':'bounded_component_fit_complete','mode':'run','candidate_frozen':True,'error':None,'implementation':pin,'configuration':pin,'selected':pin,'selected_key':key,'current_key':key,'vectors':{key:{'artifact':pin,'status':'complete','parameters':params}},'selected_metrics':{'key':key},'reserved_outcomes_opened':False,'automatic_extension':False}
    review={'status':'passed_component_fit_actual_review','result':pin,'configuration':pin,'implementation':pin,'selected':pin,'selected_native_admission_passed':True}
    R.candidate_contract(cfg,mm,fr,fc,selected,review);check('faithful_terminal_component_schema',True)
    for label,field,value in [('failed','status','failed_component_fit'),('preparation','mode','prepare'),('unfrozen','candidate_frozen',False),('terminal_error','error','failure')]:
        changed=copy.deepcopy(fr);changed[field]=value
        reject(label+'_candidate_rejected',lambda:R.candidate_contract(cfg,mm,changed,fc,selected,review))
    wrong=copy.deepcopy(fr);wrong['status']='bounded_joint_response_fit_complete'
    reject('old_fit_status_not_coerced',lambda:R.candidate_contract(cfg,mm,wrong,fc,selected,review))
    wrong=copy.deepcopy(review);wrong['selected_native_admission_passed']=False
    reject('unadmitted_candidate_rejected',lambda:R.candidate_contract(cfg,mm,fr,fc,selected,wrong))
    wrong=copy.deepcopy(selected);wrong['parameters'][0]=.2
    reject('selected_parameter_mismatch_rejected',lambda:R.candidate_contract(cfg,mm,fr,fc,wrong,review))
    with tempfile.TemporaryDirectory(prefix='component_regression_controls_') as td:
        td=Path(td);host=socket.gethostname()
        cfg={'hostname':host,'output_directory':str(td)};release={'execution_mode':'login','aggregate_login_limit':16,'workers':1}
        with patch.dict(os.environ,{'SLURM_JOB_ID':'','SLURM_STEP_ID':''}):R.login_contract(cfg,release,td)
        check('truthful_serial_login_contract',True)
        with patch.dict(os.environ,{'SLURM_JOB_ID':'123'}):reject('fabricated_scheduler_identity_rejected',lambda:R.login_contract(cfg,release,td))
        # Actual run loop, synthetic science records, native helper replaced completely.
        fake=copy.deepcopy(m);fake['initial_vector']=[0.,0.]
        for cid in R.IDS:fake['rows'][cid]['metal']='La' if '_La__' in cid else 'Dy'
        allterms={k for t in m['cycles'].values() for k in t['model_terms']}
        static={k:{'energy_kcal_mol':0.,'mu_Debye':[0.,0.,0.]} for k in allterms-set(R.IDS)}
        previous={'parameters':[2.,2.],'cycles':{},'rows':{cid:{'status':'complete','inherited_native':True,'energy':0.,'mu_Debye':[0.,0.,0.]} for cid in R.IDS}}
        em={'static':static};calls=[]
        def parse(proxy,record,row):return {'energy':0.,'mu_Debye':[0.,0.,0.]}
        def invoke(argv,folder,out,receipts):
            calls.append(argv);receipts.append({'argv':argv,'cwd':str(folder),'returncode':0})
        ev=SimpleNamespace(effective=lambda v,metal,domain:v[0 if metal=='La' else 1],render=lambda *a,**k:('synthetic key\n',{}),Evaluator=SimpleNamespace(_parse=parse))
        proxy=SimpleNamespace(h={'resolver':None,'native_helper':SimpleNamespace(invoke=invoke)},previous_regression=previous)
        ec={'executables':{'analyze':{'path':'MOCK_ONLY_NO_EXECUTABLE'}}};frr={'selected':pin,'selected_metrics':{}}
        sel={'score':{'original_EF_checks_passed':False,'old_EF':{}},'cells':{}}
        def run_case(name,values,fail=False):
            folder=td/name;folder.mkdir();cc={'implementation':pin,'manifest':R.pin(manifest_path),'fit_result':pin,'fit_review':pin,'evaluator_config':pin,'hostname':host,'output_directory':str(folder)}
            sr={'status':'passed_component_postselection_regression_source_review','implementation':pin,'manifest':cc['manifest']};R.save(td/(name+'_source.json'),sr);cc['source_review']=R.pin(td/(name+'_source.json'))
            rel=dict(cc,status='released_component_postselection_regression',selected=pin,maximum_native_calls=6,workers=1,reserved_outcomes_authorized=False,execution_mode='login',aggregate_login_limit=16);R.save(td/(name+'_release.json'),rel);cc['root_release']=R.pin(td/(name+'_release.json'))
            returned=(fake,ev,ec,em,frr,sel,values,{},proxy)
            def failcall(*args):calls.append(['synthetic failure']);raise RuntimeError('synthetic native failure')
            calls.clear()
            with patch.object(R,'context',return_value=returned),patch.dict(os.environ,{'SLURM_JOB_ID':'','SLURM_STEP_ID':''}),patch.object(proxy.h['native_helper'],'invoke',side_effect=failcall if fail else invoke):
                r=R.run(cc,folder)
            return r,len(calls)
        r,n=run_case('six',[1.,1.]);check('six_serial_calls_all_eight_cycles',n==r['native_program_invocations']==6 and r['status']=='consumed_trimer_regression_complete' and len(r['cycles'])==8)
        check('truthful_result_metadata_and_failed_EF_reporting',r['job_id'] is None and r['executor']['mode']=='login' and r['cost']['scheduler_allocated_CPU_seconds']==0 and r['original_EF_checks_passed'] is False)
        r,n=run_case('none',[0.,0.]);check('zero_exact_baseline_reuse',n==0 and r['native_program_invocations']==0 and r['status']=='consumed_trimer_regression_complete')
        r,n=run_case('three',[0.,1.]);check('three_one_metal_reuse',n==3 and r['native_program_invocations']==3)
        r,n=run_case('prior',[2.,2.]);check('prior_admitted_effective_parameter_reuse',n==0 and all(x.get('reused_native') for x in r['rows'].values()))
        r,n=run_case('failed',[1.,1.],True);check('first_failure_stops_five_untouched',n==r['native_program_invocations']==1 and r['status']=='consumed_trimer_regression_failed' and sum(x['status']=='not_attempted' for x in r['rows'].values())==5 and len(r['cycles'])==8)
    result={'status':'passed_component_postselection_saved_mock_controls','implementation':R.pin(HERE/'RUN.py'),'test_source':R.pin(__file__),'manifest':R.pin(manifest_path),'checks':checks,'passed':len(checks),'new_native_calls':0,'new_reference_calls':0,'actual_candidate_bound':False,'cost':{'wall_seconds':time.monotonic()-wall,'process_CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    with Path(out).open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'status':result['status'],'passed':len(checks),'cost':result['cost']}))

if __name__=='__main__':main(sys.argv[1],sys.argv[2])
