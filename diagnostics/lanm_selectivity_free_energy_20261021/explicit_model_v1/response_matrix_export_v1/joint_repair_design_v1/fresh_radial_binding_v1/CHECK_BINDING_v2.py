#!/usr/bin/env python3
"""Narrow new-interface controls; no actual candidate or fresh output is opened."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('binding_controls',HERE/'BIND.py');B=importlib.util.module_from_spec(s);s.loader.exec_module(B)

def run(preparation,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((Path(preparation)/'CONFIG_BIND_TEMPLATE.json').read_text());m,inp,ec,ev=B.load(cfg);checks=[]
    def check(ok,name):B.require(ok,name);checks.append(dict(name=name,passed=True))
    def rejected(call,needle,name):
        try:call()
        except ValueError as exc:check(needle in str(exc),name);return
        raise AssertionError(name)
    rejected(lambda:B.candidate(cfg),'intentionally unbound','All actual candidate gates remain closed')
    for missing in B.FUTURE:
        c=copy.deepcopy(cfg);c.update({k:dict(path='synthetic_unused',sha256='0'*64) for k in B.FUTURE});c[missing]=None
        rejected(lambda:B.candidate(c),'intentionally unbound','Missing '+missing+' closes binder before reading outcomes')
    wrong=copy.deepcopy(cfg);wrong['implementation']['sha256']='0'*64
    rejected(lambda:B.load(wrong),'configuration','Changed binder source rejected')
    # Exercise the new collector's missing-result path. Only the reference
    # collector is replaced by explicit unavailable synthetic slots; the old
    # pure comparison and scorer execute unchanged on the actual geometry graph.
    C=B.module(m['sources']['COLLECT'],'binding_missing_collection');real_module=C.B.module
    bound=dict(rows={g:dict(r,files=dict(xyz=r['input']['xyz'])) for g,r in m['rows'].items()},anchor_reference=m['anchor_reference'],
               anchor_model={x:dict(energy=0.,forces=[[0.,0.,0.]]*28) for x in ('La','Dy')},candidate_result=None,candidate_review=None,selected=None,parameters=[0.]*12,root_acceptance=None)
    synthetic_bound=dict(path=str(out/'SYNTHETIC_BOUND_NOT_READ.json'),sha256='0'*64)
    reference=dict(declared=4,available=0,rows={g:dict(geometry=r['geometry'],input=r['input'],bodies=r['bodies'],reference=None,status='unavailable') for g,r in m['rows'].items()},anchor_reference=m['anchor_reference'],allocated_native_CPU_seconds=0.,binding=dict(bound=synthetic_bound,selected=None,parameters=[0.]*12,root_acceptance=None))
    body=B.module(ec['body_helper'],'binding_control_body');ctx=(m,inp,None,None,None,None,None,None,None,SimpleNamespace(h={'body_helper':body}))
    C.B.bound_context=lambda unused:(bound,ctx)
    def module(p,name):
        if p==m['sources']['REFERENCE']:
            return SimpleNamespace(collect=lambda *args:B.save(args[-1],reference))
        if p==m['sources']['MODEL']:return SimpleNamespace()
        return real_module(p,name)
    C.B.module=module
    release=B.save(out/'SYNTHETIC_UNAVAILABLE_RELEASE.json',dict(synthetic=True));inventory=B.save(out/'SYNTHETIC_UNAVAILABLE_INVENTORY.json',dict(synthetic=True))
    c=dict(implementation=B.pin(HERE/'COLLECT.py'),bound=synthetic_bound,reference_release=release,reference_result_inventory=inventory,model_result_path=str(out/'intentionally_missing_model.json'))
    r=C.collect(c,out/'missing_result')
    check(r['status']=='validation_unavailable','Missing model remains unavailable')
    check(set(r['model_admission']['rows'])==set(B.IDS) and r['model_admission']['complete']==0,'Four model slots retained')
    check(r['score']['declared_checks']==12 and r['score']['available_checks']==0 and r['score']['passed_checks']==0 and r['score']['all_validation_checks_passed'] is False,'All twelve missing checks retained, none passes')
    check(r['new_native_calls']==0,'Saved collector makes zero native calls')
    result=dict(status='passed_narrow_binding_controls',implementation=B.pin(__file__),sources=m['sources'],manifest=cfg['manifest'],checks=checks,
                synthetic_missing_result_test=True,actual_candidate_read=False,fresh_outcomes_read=False,new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    return B.save(out/'RESULT.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preparation',required=True);p.add_argument('--out',required=True);a=p.parse_args();print(json.dumps(run(a.preparation,a.out)))
