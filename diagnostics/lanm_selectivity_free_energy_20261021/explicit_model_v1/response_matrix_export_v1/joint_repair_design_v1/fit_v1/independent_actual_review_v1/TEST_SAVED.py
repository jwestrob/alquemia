#!/usr/bin/env python3
"""Small saved-initial and selection controls; never opens the live fit."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import time
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    q=json.loads(a.config.read_text());s=importlib.util.spec_from_file_location('review',q['implementation']['path']);r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
    assert all(os.environ.get(k)=='1' for k in r.THREADS)
    start,cpu=time.monotonic(),time.process_time();checks={}
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError,AssertionError):checks[name]=True;return
        raise AssertionError('Accepted '+name)
    for value in ({'stage':'run','status':'started','cost':{}},{'stage':'run','status':'bounded_joint_response_fit_complete'}, {'stage':'prepare','status':'prepared_joint_response_fit','cost':{}}):
        reject('terminal_guard_'+str(len(checks)),lambda v=value:r.terminal_guard(v))
    r.terminal_guard({'stage':'run','status':'bounded_joint_response_fit_complete','cost':{}});checks['terminal_accepts_completed']=True
    scores={'current':{'objective':3.},'derivative':{'objective':-10.},'tie':{'objective':3.},'trial':{'objective':2.},'same':{'objective':2.}}
    assert r.choose('current',['tie'],scores)=='current';checks['current_wins_tie']=True
    assert r.choose('current',['tie','trial','same','trial'],scores)=='trial';checks['ordered_strict_trials_only']=True
    assert r.choose('current',['current','tie'],scores)!='derivative';checks['derivative_only_point_excluded']=True
    lower=np.zeros(12);upper=np.ones(12);center=np.array([0.,1.]+[.5]*10);points=r.probes(center,lower,upper)
    assert len(points)==24 and points[0]==center.tolist() and points[3]==center.tolist();checks['24_clipped_derivatives']=True
    assert all(points[2*j+1][j]>points[2*j][j] for j in range(12));checks['positive_one_sided_spans']=True
    dc=r.read(q['driver_config']);ec=r.read(dc['evaluator_config']);prep=r.read(dc['preparation_result']);v=r.read(prep['vectors'][prep['initial_key']]['artifact'])
    targets=r.read(ec['targets']);data=r.read(ec['development_data']);h={k:r.module(ec[k],'control_'+k) for k in ('objective','old_objective','body_helper')}
    r.score_audit(v,targets,data,h);checks['actual_saved_initial_531_and_25_exact_replay']=True
    bad=copy.deepcopy(v);bad['score']['block_objectives']['chi']+=1
    reject('changed_response_block',lambda:r.score_audit(bad,targets,data,h))
    bad=copy.deepcopy(v);bad['score']['weighted_residuals'][413]+=1
    reject('changed_pair_weight_or_order',lambda:r.score_audit(bad,targets,data,h))
    bad=copy.deepcopy(v);bad['response']['chi_states_a0_cubed'].pop('Dy_trimer')
    reject('missing_state_denominator',lambda:r.score_audit(bad,targets,data,h))
    bad=copy.deepcopy(v);bad['score']['original_EF_checks_passed']=not bad['score']['original_EF_checks_passed']
    reject('hidden_EF_failure',lambda:r.score_audit(bad,targets,data,h))
    assert q['fit_result'] is None;checks['terminal_authority_unbound']=True
    a.out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed_saved_terminal_reviewer_controls','implementation':r.pin(__file__),'reviewer':q['implementation'],'configuration':r.pin(a.config),'initial_artifact':prep['vectors'][prep['initial_key']]['artifact'],'checks':checks,'passed':sum(checks.values()),'initial_objective':v['score']['objective'],'new_native_calls':0,'new_parameter_evaluations':0,'live_fit_read':False,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (a.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'result':r.pin(a.out/'RESULT.json'),'passed':result['passed']}))


if __name__=='__main__':main()
