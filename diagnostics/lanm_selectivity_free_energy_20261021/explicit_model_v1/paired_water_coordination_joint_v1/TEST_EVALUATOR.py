#!/usr/bin/env python3
"""Saved-only focused ten-cell evaluator checks. Native invoke is never allowed."""
import copy
import json
import os
from pathlib import Path
import sys
import time
import importlib.util
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('joint_evaluator_controls',H/'EVALUATE.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)

def main(out):
    start=time.monotonic();cpu=time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    e.require(all(os.environ.get(k)=='1' for k in e.THREADS),'One thread required')
    cfg=json.loads((H/'EVALUATOR_CONFIG.json').read_text());ev=e.Evaluator(cfg,out/'cache');checks=[]
    def check(v,label):
        if not v:raise AssertionError(label)
        checks.append(label)
    def rejects(fn,label):
        try:fn()
        except (ValueError,RuntimeError):checks.append(label)
        else:raise AssertionError(label)
    native_attempts=[]
    def forbidden(*args,**kwargs):native_attempts.append('unexpected');raise AssertionError('No actual native allowed')
    ev.native.invoke=forbidden
    base=ev.evaluate(e.START);again=ev.evaluate(e.START)
    for gid in e.CELLS:
        check(all(base[gid][k]==ev.initial_model[gid][k]==again[gid][k] for k in ('energy','forces','native')),gid+' exact initial replay')
        check(base[gid]['provenance']['inherited_native'],gid+' prior producer retained')
        check(base[gid]['provenance']['native_implementation']!=cfg['implementation'],gid+' producer not relabeled')
    check(ev.snapshot()['native_program_invocations']==0 and ev.snapshot()['distinct_parameter_vectors']==1,'Repeated initial vector zero cost')
    baseline={g:ev.cache_identity(e.START,g)[0] for g in e.CELLS};derivative_inputs=set()
    for j,p in enumerate(ev.family['parameters']):
        width=p['upper']-p['lower'];active='La' if j in (0,1,4) else 'Dy'
        for sign in (-1,1):
            v=e.START[:];v[j]+=sign*.005*width
            changed=[]
            for gid in e.CELLS:
                key=ev.cache_identity(v,gid)[0]
                if key!=baseline[gid]:changed.append(gid);derivative_inputs.add((gid,key))
                else:
                    same=ev._cell(v,gid)
                    check(same['energy']==base[gid]['energy'] and same['forces']==base[gid]['forces'] and same['provenance']['requested_vector']==v,'Opposite-metal reuse '+str(j)+str(sign)+gid)
            check(len(changed)==5 and all(g.endswith('_'+active) for g in changed),'Exactly five active cells '+str(j)+str(sign))
        for vbad in (float('nan'),float('inf'),True,p['lower']-1e-8,p['upper']+1e-8):
            v=e.START[:];v[j]=vbad;rejects(lambda:ev.checked_values(v),'Reject nonfinite/type/bound '+str(j)+' '+repr(vbad))
    check(len(derivative_inputs)==60 and len(derivative_inputs)*2==120,'Full six-column batch at most120 programs')
    proposal_inputs=set()
    for fraction in (-.02,-.01,.01,.02):
        v=[a+fraction*(p['upper']-p['lower']) for a,p in zip(e.START,ev.family['parameters'])]
        for gid in e.CELLS:proposal_inputs.add((gid,ev.cache_identity(v,gid)[0]))
    check(len(proposal_inputs)==40 and 2*len(proposal_inputs)==80,'Four full proposals at most80 programs')
    # A real historical17605 point distinct from the selection is reused by
    # exact input identity, including the actual old producer metadata.
    prior=e.read(cfg['prior_result']);historical=next(v for v in prior['vectors'].values() if v['parameters']!=e.START)
    for gid in e.OLD_CELLS:
        value=ev._cell(historical['parameters'],gid)
        check(all(value[k]==historical['model'][gid][k] for k in ('energy','forces','native')),'Historical real-cell reuse '+gid)
    check(not native_attempts and ev.snapshot()['native_program_invocations']==0,'No native for cache controls')
    untouched=copy.deepcopy(ev._ledger);v=e.START[:];v[0]+=.000123
    ev._ledger['vector_keys']=['synthetic_'+str(i) for i in range(385)]
    rejects(lambda:ev.evaluate(v),'385 total-vector hard ceiling before evaluation');ev._ledger=copy.deepcopy(untouched)
    # Synthetic call-cap and failure controls exercise the real _cell path;
    # numerical programs are replaced, never run, and this cache is test-only.
    ev._gate=lambda:e.require(not ev.stop_event.is_set(),'Synthetic stopped gate')
    os.environ['SLURM_JOB_ID']='saved_control_only'
    ev._ledger['native_program_invocations']=5120
    rejects(lambda:ev._cell(v,e.CELLS[0]),'5120 program ceiling prevents invoke')
    check(not native_attempts and ev.stop_event.is_set(),'Ceiling fail closes before invoke')
    ev._ledger=copy.deepcopy(untouched);ev.stop_event.clear()
    simulated=[]
    def failing(*args,**kwargs):simulated.append(args[0]);raise RuntimeError('Declared synthetic native failure')
    ev.native.invoke=failing;v[0]+=.000111
    rejects(lambda:ev.evaluate(v),'Synthetic native failure propagates')
    check(len(simulated)==1 and ev.stop_event.is_set() and ev._ledger['stopped'],'First failure stops native starts')
    request=ev.snapshot()['requests'][-1];check(request['status']=='failed' and set(request['cells'])==set(e.CELLS),'Failure retains all ten slots')
    v[0]+=.000111;rejects(lambda:ev._cell(v,e.CELLS[0]),'Later program rejected after failure');check(len(simulated)==1,'No second invoke after failure')
    result={'status':'passed','implementation':e.pin(__file__),'evaluator':e.pin(H/'EVALUATE.py'),'config':e.pin(H/'EVALUATOR_CONFIG.json'),
        'checks':checks,'checks_passed':len(checks),'actual_native_program_invocations':0,'synthetic_invoke_attempts':len(simulated),
        'training_result':cfg['training_result'],'prior_files_unchanged':True,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}
    for p in ev.data['sources'].values():e.verify(p)
    e.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':e.pin(out/'RESULT.json')}))
if __name__=='__main__':main(sys.argv[1])
