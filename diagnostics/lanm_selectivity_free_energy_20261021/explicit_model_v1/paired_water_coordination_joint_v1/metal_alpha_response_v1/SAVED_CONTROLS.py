#!/usr/bin/env python3
"""Prepared inputs and bounded negative controls; never invokes native programs."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import re
import time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--preparation',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    start=time.monotonic();cpu=time.process_time();out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
    p=Path(__file__).with_name('PROBE.py');s=importlib.util.spec_from_file_location('alpha_controls',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    cfg=json.loads(Path(a.config).read_text());prep=m.read(m.pin(a.preparation));man=m.read(prep['manifest']);data=m.read(man['data']);checks=[]
    def check(ok,label):m.require(ok,label);checks.append(label)
    def rejects(fn,label):
        try:fn()
        except (ValueError,KeyError,IndexError,TypeError):checks.append(label);return
        raise AssertionError('Expected rejection: '+label)
    resolver=m.module(cfg['resolver'],'alpha_test_resolver');check(cfg['implementation']==m.pin(p) and man['implementation']==m.pin(p),'final source pins')
    check(len(man['probes'])==8 and sum(len(q['inputs']) for q in man['probes'])==40,'eight vectors and forty new cells')
    check(len(data['six_column_Jacobian'])==413 and all(len(v)==6 for v in data['six_column_Jacobian']),'413 by six saved Jacobian')
    check(data['baseline_score']['declared_energy_controls']==10 and data['baseline_score']['declared_force_controls']==15 and data['baseline_score']['development_checks_passed'] is False,'unchanged failed scientific controls')
    for q in man['probes']:
        check(len(q['inputs'])==len(q['reused_cells'])==5 and set(q['inputs'])|set(q['reused_cells'])==set(m.CELLS),'complete denominator '+q['id'])
        check(all(not g.endswith('_'+q['metal']) for g in q['reused_cells']),'opposite-metal exact reuse '+q['id'])
        for gid,inp in q['inputs'].items():
            base=m.verify(inp['baseline_input']['key']).read_text();prm=m.verify(inp['baseline_input']['parameters']).read_text();key,render=m.alpha_key(base,prm,q['metal'],q['alpha'],resolver)
            check(key==m.verify(inp['files']['key']).read_text() and render==inp['render'],'exact rendered key '+q['id']+'/'+gid)
            b=re.findall(r'\S+',base);n=re.findall(r'\S+',key);diff=[j for j,(x,y) in enumerate(zip(b,n)) if x!=y]
            check(len(b)==len(n) and len(diff)==1 and float(n[diff[0]])==q['alpha'] and float(b[diff[0]])==m.ALPHA[q['metal']],'exactly one numeric alpha token '+q['id']+'/'+gid)
            check(all(inp['files'][k]['sha256']==inp['baseline_input'][k]['sha256'] for k in ('xyz','parameters')),'unchanged geometry and parameters '+q['id']+'/'+gid)
    example=man['probes'][0];inp=next(iter(example['inputs'].values()));base=m.verify(inp['baseline_input']['key']).read_text();prm=m.verify(inp['baseline_input']['parameters']).read_text()
    for val in (0,-1,float('nan'),float('inf'),True):rejects(lambda val=val:m.alpha_key(base,prm,'La',val,resolver),'invalid alpha '+repr(val))
    rejects(lambda:m.alpha_key(base+'\npolarize 500 1.6597 0.25\n',prm,'La',1.7,resolver),'duplicate metal record')
    rejects(lambda:m.alpha_key(base.replace('1.6597','1.6600'),prm,'La',1.7,resolver),'changed alpha origin')
    rejects(lambda:m.alpha_key(base,prm,'X',1.7,resolver),'unknown metal')
    c=m.Counter()
    for i in range(80):check(c.reserve()==i+1,'finite program reservation '+str(i+1))
    rejects(c.reserve,'eighty-first program blocked');check(c.count==80,'failed reservation does not increase counter')
    c=m.Counter();c.stop.set();rejects(c.reserve,'global stop blocks first program');check(c.count==0,'stopped count unchanged')
    review=m.read(cfg['selected_review']);r={'status':'bounded_joint_coordination_water_complete','job_id':'17676','error':None,'stopping_reason':'no_actual_improvement_at_complete_trial_batch',
        'selected':{'parameters':data['baseline_parameters'],'status':'complete','model':data['baseline_model'],'score':data['baseline_score']},'development_checks_passed':False}
    m.admitted_terminal(r,review,cfg['candidate_result'],cfg['reviewer']);check(True,'actual raw-pass scientific-fail distinction admitted')
    for key,val in (('all_selected_raw_replays_exact',False),('exact_selected_score_replayed',False),('all_development_controls_passed',True),('complete',9),('new_native_calls',1)):
        bad=copy.deepcopy(review);bad[key]=val;rejects(lambda bad=bad:m.admitted_terminal(r,bad,cfg['candidate_result'],cfg['reviewer']),'invalid review '+key)
    bad=copy.deepcopy(r);bad['development_checks_passed']=True;rejects(lambda:m.admitted_terminal(bad,review,cfg['candidate_result'],cfg['reviewer']),'passing-science substitute rejected')
    result={'status':'saved_controls_passed','implementation':m.pin(__file__),'probe_source':m.pin(p),'config':m.pin(a.config),'preparation':m.pin(a.preparation),
        'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    m.save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'output':str(out)}))

if __name__=='__main__':main()
