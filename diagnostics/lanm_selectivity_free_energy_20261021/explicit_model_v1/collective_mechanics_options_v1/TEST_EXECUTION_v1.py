#!/usr/bin/env python3
"""Focused execution guards; all invocations below are explicit Python fakes."""
import copy
import importlib.util
import json
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('trimer_executor_controls',HERE/'RUN.py')
r=importlib.util.module_from_spec(s);s.loader.exec_module(r)


def main():
    tick=time.monotonic();cpu=time.process_time();checks=[]
    out=r.p.W/'explicit_model_v1/collective_mechanics_options_v1/execution_controls_v1';out.mkdir(exist_ok=False)
    cfg=json.loads((HERE/'CONFIG_EXECUTE.json').read_text());m,common,parent,model,ind=r.context(cfg)
    def check(v,label):r.require(v,label);checks.append(label)
    def rejects(fn,label):
        try:fn()
        except (ValueError,FileExistsError,KeyError):checks.append(label);return
        raise AssertionError(label)
    fake_review=out/'FAKE_REVIEW.json';r.save(fake_review,{'status':'passed_TEST_ONLY_not_a_native_release'})
    release={'status':'approved_consumed_trimer_energy_moment_diagnostic','implementation':cfg['implementation'],
        'configuration':r.pin(HERE/'CONFIG_EXECUTE.json'),'manifest':cfg['manifest'],'preparation':cfg['preparation'],'declaration':cfg['declaration'],
        'reference_cell_ids':[c['cell_id'] for c in m['reference_cells']],'model_cell_ids':list(m['model_rows']),
        'maximum_reference_programs':16,'maximum_model_programs':16,'reference_ranks':4,'model_workers':4,
        'fit_authorized':False,'force_or_FD_authorized':False,'finite_field_authorized':False,
        'preparation_review':r.pin(r.p.W/'explicit_model_v1/collective_mechanics_options_v1/independent_review_v1/REVIEW.json'),
        'execution_review':r.pin(fake_review),'reference_execution_root':str(out/'fake_reference'),'model_execution_root':str(out/'fake_model_success')}
    r.release(cfg,m,release);checks.append('exact finite release admitted without execution')
    for key,value in [('maximum_reference_programs',17),('maximum_model_programs',17),('reference_ranks',8),('model_workers',8),('fit_authorized',True),('force_or_FD_authorized',True),('finite_field_authorized',True),('reference_cell_ids',release['reference_cell_ids'][:-1]),('model_cell_ids',release['model_cell_ids'][::-1])]:
        bad=dict(release);bad[key]=value;rejects(lambda:r.release(cfg,m,bad),'changed release rejects '+key)
    counter=r.Counter()
    for _ in range(16):counter.reserve()
    check(counter.count==16,'exact16 reservations')
    rejects(counter.reserve,'seventeenth reservation rejected');counter=r.Counter();counter.stop.set();rejects(counter.reserve,'stop prevents reservation')
    rr=out/'claimed_reference';r.reference_claim(rr,rr/'a','a');checks.append('one reference slot claim')
    rejects(lambda:r.reference_claim(rr,rr/'a','a'),'duplicate reference attempt rejected')
    rejects(lambda:r.reference_claim(rr,out/'elsewhere','b'),'undeclared reference attempt path rejected')
    r.reference_stop(rr,{'error':'TEST'});rejects(lambda:r.reference_claim(rr,rr/'b','b'),'recorded failure blocks next reference claim')
    r.reference_stop(rr,{'error':'SECOND'});check(json.loads((rr/'STOP.json').read_text())=={'error':'TEST'},'first reference failure preserved')
    check(len(m['reference_cells'])==len(m['model_rows'])==16,'actual prepared complete input inventory')
    original_context=r.context;original_runtime=r.runtime;original_module=r.module
    calls=[];fail=False
    def invoke(argv,cwd,output,receipts):
        calls.append(str(cwd));output.write_text('TEST ONLY no native execution\n');(cwd/'system.uind').write_text('TEST ONLY\n')
        receipts.append({'argv':argv,'cwd':str(cwd),'returncode':0,'output':r.pin(output),'threads':{'OMP_NUM_THREADS':'1'}})
    def parse(output,induced,row,induction):
        if fail:raise ValueError('TEST parser failure after fake native completion')
        return {'TEST_ONLY':True,'energy_kcal_mol':0.0,'complete_sites':row['atoms']}
    fake_model=SimpleNamespace(parse_moments=parse,induced_file=lambda output:output.parent/'system.uind')
    r.context=lambda config:(m,common,parent,fake_model,ind)
    r.runtime=lambda kind:None
    r.os.environ['SLURM_JOB_ID']='TEST_ONLY_NO_NATIVE'
    r.module=lambda pin,name:SimpleNamespace(invoke=invoke) if name=='trimer_execution_native' else original_module(pin,name)
    try:
        rp=out/'FAKE_RELEASE_SUCCESS.json';r.save(rp,release)
        success=r.run_model(cfg,rp,Path(release['model_execution_root']))
        check(success['status']=='complete' and success['native_program_invocations']==16 and len(calls)==16,'actual driver executes exactly16 fake rows')
        check(all(v['status']=='complete' for v in success['rows'].values()) and not success['global_stop'],'actual driver complete denominator')
        rejects(lambda:r.run_model(cfg,rp,Path(release['model_execution_root'])),'model whole-attempt rerun rejected')
        fail=True;calls.clear();release['model_execution_root']=str(out/'fake_model_failure');rp=out/'FAKE_RELEASE_FAILURE.json';r.save(rp,release)
        failed=r.run_model(cfg,rp,Path(release['model_execution_root']))
        check(failed['status']=='failed' and failed['global_stop'] and len(failed['rows'])==16,'actual parser failure stops with all16 slots')
        check(1<=len(calls)<=4 and failed['native_program_invocations']==len(calls),'failure permits only already admitted worker calls')
        check(sum(v['status']=='not_attempted' for v in failed['rows'].values())==16-len(calls),'untouched denominator explicit')
        check(all(v['status'] in ('failed','not_attempted') for v in failed['rows'].values()),'failed calls not relabeled/retried')
    finally:r.context=original_context;r.runtime=original_runtime;r.module=original_module
    result={'status':'passed_focused_execution_controls','implementation':r.pin(__file__),'executor':cfg['implementation'],
        'configuration':r.pin(HERE/'CONFIG_EXECUTE.json'),'manifest':cfg['manifest'],'checks_passed':len(checks),'checks':checks,
        'fake_success_rows':16,'fake_failure_started':len(calls),'new_native_calls':0,
        'fake_artifacts_are_not_native_authority':True,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    r.save(out/'RESULT.json',result);print(json.dumps({'result':r.pin(out/'RESULT.json'),'checks':len(checks),'native_calls':0}))


if __name__=='__main__':main()
