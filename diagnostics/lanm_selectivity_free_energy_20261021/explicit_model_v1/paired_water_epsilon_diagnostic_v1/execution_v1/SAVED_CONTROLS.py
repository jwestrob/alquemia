"""Saved-input and synthetic zero-response integration controls; no engines."""
import importlib.util
import json
from pathlib import Path
import sys
import time

start=time.monotonic();cpu=time.process_time()
here=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('epsilon_controls',here/'DIAGNOSTIC.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
cfg=json.loads((here/'CONFIG_PREPARE.json').read_text());prep=m.read(m.pin(sys.argv[1]));manifest=m.read(prep['manifest'])
adapter=m.module(cfg['authority_helper'],'control_authority');data=adapter.loaded(m.verify(cfg['original32']))
training=m.read(cfg['training']);body=m.module(cfg['body_helper'],'control_body');objective=m.module(cfg['objective'],'control_objective')
checks={}
def check(name,ok):
    checks[name]=bool(ok);m.require(ok,name)
check('final_source_binding',prep['implementation']==cfg['implementation']==m.pin(here/'DIAGNOSTIC.py'))
check('zero_call_preparation',prep['status']=='prepared'and prep['new_native_calls']==0)
check('four_fixed_probes',[p['id']for p in manifest['probes']]==['La_minus','La_plus','Dy_minus','Dy_plus'])
check('sixteen_changed_cells',sum(len(p['inputs'])for p in manifest['probes'])==16)
check('exact_opposite_metal_reuse',all(v==manifest['baseline']['model'][gid]for p in manifest['probes']for gid,v in p['reused_cells'].items()))
check('complete_derivative_authority_denominator',len(manifest['derivative_keys'])==8 and len(manifest['admitted_native_authorities'])==40)
records={p['id']:{'epsilon':p['epsilon'],'eta':p['eta'],'model':manifest['baseline']['model'],'score':manifest['baseline']['score']}for p in manifest['probes']}
analysis=m.analyze(manifest,records,data,training,body,objective)
check('synthetic_zero_response_has_no_new_rank',analysis['retained_increment_rank']==0 and all(x==0 for row in analysis['epsilon_columns_per_eta']for x in row))
check('synthetic_zero_response_passes_midpoint',analysis['parameter_linearity_passed'])
check('nested_span_unchanged_without_new_response',analysis['weighted_residual_losses']['after_original_four_span']==analysis['weighted_residual_losses']['after_nested_six_span'])
check('all241_residual_entries_reported',all(len(x)==241 for x in analysis['residual_components'].values()))
check('no_parameter_solution_or_selection',not analysis['parameter_solution_computed']and not analysis['new_candidate_selected'])
for p in manifest['protected_parents']:m.verify(p)
out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False)
m.save(out/'RESULT.json',{'status':'passed','implementation':m.pin(here/'DIAGNOSTIC.py'),'control_source':m.pin(__file__),'preparation':m.pin(sys.argv[1]),'checks':checks,'new_native_calls':0,'synthetic_analysis_is_not_probe_evidence':True,'elapsed_seconds':time.monotonic()-start,'parent_cpu_seconds':time.process_time()-cpu})
print(json.dumps({'status':'passed','checks':len(checks),'new_native_calls':0,'result':m.pin(out/'RESULT.json')}))
