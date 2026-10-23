#!/usr/bin/env python3
"""Focused collector joins/unavailability controls; no native/outcome access."""
import copy
import importlib.util
import json
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('collector',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time();out=m.WORK/'collection/saved_controls_v1';out.mkdir(parents=True,exist_ok=False)
cfg=m.load(HERE/'CONFIG_SOURCE.json');geometry=m.read(cfg['geometry_manifest']);checks=[]
def check(name,ok):
    checks.append({'name':name,'passed':bool(ok)});assert ok,name
def reject(name,fn):
    try:fn()
    except (ValueError,KeyError):check(name,True)
    else:check(name,False)
check('source_configuration_unbound',all(cfg[k] is None for k in m.FUTURE))
reject('unbound_run_rejected',lambda:m.collect(HERE/'CONFIG_SOURCE.json',out/'forbidden_unbound_run'))
check('unbound_run_creates_no_output',not (out/'forbidden_unbound_run').exists())
synthetic={'path':'synthetic-only','sha256':'not-a-real-authority'}
rows={g:{'geometry':r['geometry'],'input':r['input'],'bodies':r['bodies'],'files':{'xyz':r['input']['xyz']}} for g,r in geometry['rows'].items()}
zero={'energy':0.,'forces':[[0.,0.,0.] for _ in range(28)]}
anchors={metal:copy.deepcopy(zero) for metal in ('La','Dy')}
bound={'rows':rows,'candidate_result':synthetic,'candidate_review':synthetic,'parameters':[1]*6,'geometry_manifest':cfg['geometry_manifest'],'implementation':cfg['model_source'],'anchor_reference':anchors,'anchor_model':anchors}
mcfg={'candidate_result':synthetic,'candidate_review':synthetic,'geometry_manifest':cfg['geometry_manifest'],'implementation':cfg['model_source'],'body_helper':cfg['body_helper']}
release={'candidate_result':synthetic,'selected_review':synthetic,'parameters':[1]*6,'prepared_manifest':cfg['reference_manifest'],'configuration':cfg['reference_configuration'],'implementation':cfg['reference_source']}
m.check_binding(cfg,bound,mcfg,release,geometry);check('exact_synthetic_binding_join',True)
for key in ('candidate_result','selected_review','parameters'):
    bad=copy.deepcopy(release);bad[key]=None
    reject('reject_changed_release_'+key,lambda bad=bad:m.check_binding(cfg,bound,mcfg,bad,geometry))
bad=copy.deepcopy(bound);bad['rows'][m.GIDS[0]]['geometry']['positions_A'][0][0]+=1
reject('reject_changed_geometry',lambda:m.check_binding(cfg,bad,mcfg,release,geometry))
bad=copy.deepcopy(bound);bad['rows'][m.GIDS[0]]['input']['types'][0]=999
reject('reject_changed_typing',lambda:m.check_binding(cfg,bad,mcfg,release,geometry))
refs={'declared':4,'rows':{g:{**copy.deepcopy(r),'status':'complete','reference':dict(copy.deepcopy(zero),changed_central_energy_stage_signature=False)} for g,r in rows.items()},'anchor_reference':anchors,'allocated_native_CPU_seconds':0.}
admitted={'rows':{g:{'status':'complete','model':copy.deepcopy(zero)} for g in m.GIDS},'anchor_model':anchors}
score=m.module(cfg['score_source'],'test_frozen_score');body=m.module(cfg['body_helper'],'test_frozen_body')
complete=m.compare(bound,refs,admitted,score,body)
check('delegated_twelve_synthetic_self_checks',complete['declared_checks']==12 and complete['available_checks']==12 and complete['all_validation_checks_passed'])
partial=copy.deepcopy(refs);partial['rows'][m.GIDS[0]].update(status='unavailable',reference=None)
result=m.compare(bound,partial,admitted,score,body)
check('missing_reference_retains_twelve_and_cannot_pass',result['declared_checks']==12 and result['available_checks']<12 and not result['all_validation_checks_passed'])
bad=copy.deepcopy(refs);bad['anchor_reference']['La']['energy']=1
reject('reject_changed_central_reference_anchor',lambda:m.compare(bound,bad,admitted,score,body))
counter={'bound':0,'admit':0}
def bound_inputs(*args):counter['bound']+=1
def admit_result(*args):counter['admit']+=1;return copy.deepcopy(admitted)
fake=SimpleNamespace(bound_inputs=bound_inputs,admit_result=admit_result)
testcfg={**cfg,'model_result_path':str(out/'missing.json'),'model_bound':synthetic}
a,c,p=m.admit_model(testcfg,bound,mcfg,fake)
check('missing_model_keeps_all_four_without_raw_admission',a['complete']==0 and set(a['rows'])==set(m.GIDS) and counter['admit']==0)
path=out/'malformed_synthetic.json';path.write_text('{"truncated":');testcfg['model_result_path']=str(path)
a,c,p=m.admit_model(testcfg,bound,mcfg,fake)
check('malformed_model_keeps_raw_pin_and_four_slots',a['complete']==0 and a['raw_result']==m.pin(path) and 'JSONDecodeError' in a['reason'] and counter['admit']==0)
result=m.compare(bound,refs,a,score,body)
check('malformed_model_twelve_unavailable_not_pass',result['available_checks']==0 and result['declared_checks']==12 and result['all_validation_checks_passed'] is False)
raw={'bound':synthetic,'candidate_result':synthetic,'candidate_review':synthetic,'native_program_invocations':8,'elapsed_seconds':1.}
path=out/'synthetic_model.json';m.save(path,raw);testcfg['model_result_path']=str(path)
a,c,p=m.admit_model(testcfg,bound,mcfg,fake)
check('parseable_model_uses_existing_admission_once',counter['admit']==1 and a['rows']==admitted['rows'])
check('saved_cost_fields_preserved',c['native_program_invocations']==8 and c['elapsed_seconds']==1.)
def failed(*args):raise ValueError('synthetic native admission failure')
fake.admit_result=failed;a,c,p=m.admit_model(testcfg,bound,mcfg,fake)
check('native_admission_failure_keeps_four_unavailable',a['complete']==0 and set(a['rows'])==set(m.GIDS) and 'synthetic native admission failure' in a['reason'])
text=m.report({'status':'validation_unavailable','score':result,'reference_admission':refs})
check('report_keeps_twelve_rows_and_limits',text.count('| unavailable |')==12 and 'New native calls in this collection: 0' in text and 'provisional frozen-f' in text)
receipt={'status':'saved_collector_controls_passed','implementation':m.pin(HERE/'RUN.py'),'configuration':m.pin(HERE/'CONFIG_SOURCE.json'),'test_source':m.pin(__file__),'checks':checks,'check_count':len(checks),'new_native_calls':0,'fit_or_fresh_outcomes_read':False,'controls_use':'frozen geometry metadata plus explicitly synthetic energies/forces/authorities only','local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
print(json.dumps(m.save(out/'RESULT.json',receipt)))
