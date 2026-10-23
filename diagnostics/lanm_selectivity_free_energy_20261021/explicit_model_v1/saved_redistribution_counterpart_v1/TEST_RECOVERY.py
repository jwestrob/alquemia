#!/usr/bin/env python3
"""Narrow saved4/readback/staging controls; zero molecular calls."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('moment_recovery_test',HERE/'RUN_v2.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def run(config,preparation,out):
    start=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config).read_text());prep=json.loads(Path(preparation).read_text());manifest=m.read(prep['manifest']);checks=[]
    def ok(name,value):
        if not value:raise AssertionError(name)
        checks.append(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,KeyError,IndexError):checks.append(name);return
        raise AssertionError('Accepted '+name)
    m.bound_inputs(cfg,manifest);helper=m.module(cfg['induction_helper'],'recovery_test_induction')
    ok('four_reused_74_pending',len(manifest['reused_rows'])==4 and len(manifest['pending_ids'])==74 and not set(manifest['reused_rows'])&set(manifest['pending_ids']))
    ok('old78_input_slots',len(manifest['rows'])==78)
    ok('ions_first_then72',manifest['stages']['bare_ions']==['ion_La','ion_Dy'] and len(manifest['stages']['remaining'])==72)
    for cid,r in manifest['reused_rows'].items():
        ok('old_producer_'+cid,r['inherited_native_authority']['result']==cfg['previous_result'] and r['inherited_native_authority']['producer_implementation']==cfg['previous_adapter'])
        ok('actual_uind_'+cid,Path(r['moments']['induced_output']['path']).name=='system.uind')
        ok('unchanged_moment_admission_'+cid,r['moments']['induction']['solutions']==2 and r['moments']['induction']['maximum_printed_RMS_residual_Debye']<=1e-8)
    first=next(iter(manifest['reused_rows']));r=manifest['reused_rows'][first];row=manifest['rows'][first];folder=out/'explicit_filename_controls';folder.mkdir();a=folder/'analyze.out';u=folder/'system.uind';text=m.verify(r['moments']['analyze_output']).read_text();utext=m.verify(r['moments']['induced_output']).read_text();a.write_text(text);u.write_text(utext)
    ok('archive_announcement',m.induced_file(a)==u)
    u.rename(folder/'system.001u');a.write_text(text.replace('system.uind','system.001u'));ok('numbered_announcement',m.induced_file(a)==folder/'system.001u')
    u.write_text(utext);rejects('two_files_rejected',lambda:m.induced_file(a));(folder/'system.001u').unlink();a.write_text(text)
    badrow=copy.deepcopy(row);badrow['native_names']=badrow['names'];rejects('XYZ_La_not_native_La3',lambda:m.parse_moments(a,u,badrow,helper.induction))
    u.write_text(utext+utext);rejects('second_frame_rejected',lambda:m.parse_moments(a,u,row,helper.induction));u.write_text(utext)
    a.write_text(text.replace('system.uind','../system.uind'));rejects('external_filename_rejected',lambda:m.induced_file(a));a.write_text(text)
    counter=m.Counter()
    for i in range(74):counter.reserve()
    rejects('new_call75_rejected',counter.reserve);counter=m.Counter();counter.stop.set();rejects('global_stop',counter.reserve)
    rejects('unreleased_recovery',lambda:m.release_gate(cfg,prep,False))
    class SerialExecutor:
        def __init__(self,**kw):pass
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def map(self,fn,items):return map(fn,items)
    original={n:getattr(m,n) for n in ('checked_config','bound_inputs','release_gate','read','module')};pool=m.concurrent.futures.ThreadPoolExecutor;job=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='explicit_synthetic_recovery'
    testcfg=dict(cfg,preparation_result={'synthetic':'preparation'})
    try:
        m.checked_config=lambda *a,**k:None;m.bound_inputs=lambda *a:manifest;m.release_gate=lambda *a,**k:None;m.read=lambda p:prep if p==testcfg['preparation_result'] else manifest;m.concurrent.futures.ThreadPoolExecutor=SerialExecutor
        for scenario in ('ion_failure','two_ions_then_failure'):
            admitted=[]
            def invoke(argv,directory,output,receipts):
                cid=directory.name;admitted.append(cid)
                if scenario=='ion_failure' or not cid.startswith('ion_'):raise ValueError('EXPLICIT SYNTHETIC failure, no process started')
                rr=manifest['rows'][cid]
                output.write_text('EXPLICIT SYNTHETIC\nTinker Version 26.2\nInduced Dipoles : Iterations 1 RMS Residual 0.0000000000\nTotal Potential Energy : 0.00000000 Kcal/mole\nEnergy Component Breakdown :\nInduced Dipoles : Iterations 1 RMS Residual 0.0000000000\nTotal Electric Charge : 3.00000 Electrons\nDipole X,Y,Z-Components : 0.000 0.000 0.000\nInduced Dipoles Written To : system.uind\n')
                (directory/'system.uind').write_text('1 EXPLICIT SYNTHETIC\n1 '+rr['native_names'][0]+' 0.000000 0.000000 0.000000\n');receipts.append({'explicit_synthetic':True,'process_started':False})
            m.module=lambda p,n:SimpleNamespace(invoke=invoke) if p==cfg['native_helper'] else helper
            dest=out/scenario;dest.mkdir();result=m.run(testcfg,dest)
            ok(scenario+'_full78',len(result['rows'])==78 and result['status']=='failed_or_unavailable' and result['global_stop'])
            ok(scenario+'_preserved4',all(result['rows'][cid]==value for cid,value in manifest['reused_rows'].items()))
            ok(scenario+'_no_duplicate',not set(admitted)&set(manifest['reused_rows']))
            ok(scenario+'_ledger',result['native_program_invocations']==len(admitted) and result['inherited_native_program_invocations']==4 and result['total_native_program_invocations']==4+len(admitted))
            ok(scenario+'_stage_order',admitted==(['ion_La'] if scenario=='ion_failure' else ['ion_La','ion_Dy',manifest['stages']['remaining'][0]]))
    finally:
        for name,fn in original.items():setattr(m,name,fn)
        m.concurrent.futures.ThreadPoolExecutor=pool
        if job is None:os.environ.pop('SLURM_JOB_ID',None)
        else:os.environ['SLURM_JOB_ID']=job
    result={'status':'passed','implementation':m.pin(HERE/'RUN_v2.py'),'test_implementation':m.pin(__file__),'specification':m.pin(HERE/'SPEC_v2.md'),'config':m.pin(config),'preparation':m.pin(preparation),'manifest':prep['manifest'],'checks_passed':len(checks),'checks':checks,'new_native_calls':0,'synthetic_processes_started':0,'max_local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    m.save(out/'RESULT.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=run(a.config,a.preparation,out)
    except Exception as exc:m.save(out/'FAILURE.json',{'status':'failed','error':repr(exc),'new_native_calls':0});raise
    print(json.dumps({'status':r['status'],'checks_passed':r['checks_passed'],'new_native_calls':0}))
