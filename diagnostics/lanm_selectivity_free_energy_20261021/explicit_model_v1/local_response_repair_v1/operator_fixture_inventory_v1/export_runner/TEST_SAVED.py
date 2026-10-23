#!/usr/bin/env python3
"""Narrow orchestration checks. Every executor below is a fake Python function."""
import argparse, copy, importlib.util, json, time
from pathlib import Path
from types import SimpleNamespace

ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
p=Path(__file__).with_name('RUN.py');sp=importlib.util.spec_from_file_location('remaining_export_controls',p);run=importlib.util.module_from_spec(sp);sp.loader.exec_module(run)
out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
cfg=run.read(run.pin(a.config));m,h=run.context(cfg);checks={}
checks['actual_prepared_54_rows']=len(m['rows'])==54
checks['actual_remaining53_excludes_scout']=len(m['remaining_cell_ids'])==53 and run.SCOUT not in m['remaining_cell_ids']
checks['reporting_flags_only']=all(run.verify(r['input']['key']).read_text()==run.render(run.verify(r['source_input']['key']).read_text())[0] for r in m['rows'].values())
checks['ten_original_keys_add_both_old_reporting_flags']=sum(len(r['reporting_flags_added'])==3 for r in m['rows'].values())==10
checks['48_saved_u_baselines']=sum(r['baseline_u_available'] for r in m['rows'].values())==48
checks['two_ion_controls']=sum(r['kind']=='ion' for r in m['rows'].values())==2
altered=copy.deepcopy(cfg);altered['resolver']=cfg['native_helper']
try:run.context(altered)
except ValueError as e:checks['different_valid_helper_pin_rejected']='helper implementation differs' in str(e)
else:checks['different_valid_helper_pin_rejected']=False
try:run.release_gate(cfg,m,out,{})
except ValueError as e:checks['unbound_gate_before_executor']='remain unbound' in str(e)
else:checks['unbound_gate_before_executor']=False
try:run.render('parameters amoeba09.prm\nEXTERNAL-FIELD 0 0 1\n')
except ValueError:checks['preexisting_field_rejected']=True
else:checks['preexisting_field_rejected']=False

def case(name,mode):
    folder=out/name;folder.mkdir();counts={'fake_invokes':0,'fake_qualifications':0}
    def invoke(argv,cwd,output,receipts):
        counts['fake_invokes']+=1
        output.write_text('MOCK ONLY; no native program\nResponse Matrix Written To : response_matrix_001.dat\nResponse Matrix Written To : response_matrix_002.dat\n')
        receipts.append({'synthetic':True,'output':run.pin(output),'argv':argv})
        for i in (1,2):(cwd/f'response_matrix_{i:03d}.dat').write_text('MOCK ONLY; no scientific matrix\n')
        if mode=='native_failure':raise ValueError('synthetic native failure')
    def parse(output,induced,row,induction):return {'synthetic':True,'energy_kcal_mol':row['baseline']['energy_kcal_mol']}
    def qualify(*args):counts['fake_qualifications']+=1;return {'synthetic':True,'passed':True}
    fake={'native_helper':SimpleNamespace(invoke=invoke),
          'moment_adapter':SimpleNamespace(induced_file=lambda p:p,parse_moments=parse),
          'induction_helper':SimpleNamespace(induction=None),
          'scout_helper':SimpleNamespace(parity=lambda *x,**y:{'synthetic':True,'passed':mode!='parity_failure'},qualify_exports=qualify)}
    cc=copy.deepcopy(cfg);cc['analyze']={'path':'SYNTHETIC_NO_EXECUTABLE','sha256':'not_a_native_pin'}
    result=run.execute_remaining(cc,m,folder,fake,{'status':'complete','synthetic':True,'inherited':True})
    return result,counts

r,c=case('synthetic_full','success')
checks['synthetic_full_exact53_invocations']=c['fake_invokes']==r['native_program_invocations']==53
checks['synthetic_full_54_108_denominator']=r['status']=='baseline_matrix_inventory_complete' and r['complete_fixtures']==54 and r['admitted_snapshots']==108
checks['synthetic_no_scout_reexecution']=not (out/'synthetic_full/cells'/run.SCOUT).exists()
for name,mode in [('synthetic_native_failure','native_failure'),('synthetic_parity_failure','parity_failure')]:
    r,c=case(name,mode)
    checks[name+'_stops_after_one']=c['fake_invokes']==r['native_program_invocations']==1 and c['fake_qualifications']==0
    checks[name+'_all54_slots_retained']=len(r['rows'])==54 and sum(v['status']=='not_attempted' for v in r['rows'].values())==52
    checks[name+'_unavailable']=r['status']=='baseline_matrix_inventory_incomplete' and r['global_stop'] and r['error'] is not None
    checks[name+'_raw_failure_files_preserved']=len(r['rows'][m['remaining_cell_ids'][0]]['exports_present'])==2
result={'status':'passed' if all(checks.values()) else 'failed','implementation':run.pin(p),'test_implementation':run.pin(__file__),
        'configuration':run.pin(a.config),'manifest':cfg['manifest'],'checks':checks,'passed':sum(checks.values()),'total':len(checks),
        'new_native_calls':0,'mock_executor_only':True,'cost':{'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'local_threads':1}}
run.save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'result':run.pin(out/'RESULT.json'),'cost':result['cost']}))
if not all(checks.values()):raise SystemExit(1)
