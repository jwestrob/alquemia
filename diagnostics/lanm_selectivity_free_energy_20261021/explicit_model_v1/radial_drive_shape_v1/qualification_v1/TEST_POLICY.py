#!/usr/bin/env python3
"""Finite metadata and serial mock controls, never a native invocation."""
import argparse,copy,hashlib,json,math,os,subprocess,tempfile,time,traceback
from pathlib import Path
import numpy as np
import PREPARE as P
import RUN as R
import CHECKS as C

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();checks=[]
    def ck(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    result={'status':'started','implementation':P.pin(__file__),'native_calls':0,'solver_calls':0,'checks':checks,'mock_results_are_not_scientific_evidence':True}
    try:
        cfg=json.loads(Path(args.config).read_text());m,h,ev=R.context(cfg,False)
        ck('49_total_unique',len(m['programs'])==len({j['id'] for j in m['programs']})==49)
        ck('132_FD_energies',sum(6*m['rows'][j['cell']]['atoms'] for j in m['programs'] if j['fd'])==132)
        ck('21_EM_42_exports',sum(j['program']=='analyze' and j['stage']!='directions' for j in m['programs'])==21)
        ck('exact_stage_counts',{s:sum(j['stage']==s for j in m['programs']) for s in ('absent','zero','active','directions','fields','bypass')}==dict(absent=4,zero=10,active=12,directions=16,fields=4,bypass=3))
        for kind,cells in (('metal_x',P.CENTERS),('water_rotation_x',P.DIMERS)):
            for g in cells:
                row=m['rows'][g];center=np.asarray(row['positions_A']);p,t=P.path_geometry(row,kind,1e-6);n,_=P.path_geometry(row,kind,-1e-6);err=np.max(abs((np.asarray(p)-n)/2e-6-t));ck(g+'_'+kind+'_analytic_tangent',err<1e-8)
                if kind=='water_rotation_x':
                    shifted=np.asarray(P.path_geometry(row,kind,1e-4)[0]);ck(g+'_rigid_water_distances',max(abs(np.linalg.norm(shifted[i]-shifted[j])-np.linalg.norm(center[i]-center[j])) for i in range(1,4) for j in range(i+1,4))<1e-14)
        bad=copy.deepcopy(cfg);bad['limits']['native_programs']=50
        try:R.context(bad,False)
        except ValueError:ck('reject_changed_program_cap',True)
        else:ck('reject_changed_program_cap',False)
        inv=R.invocation(h);ck('only_timeout_changes_inherited_invocation',120 in inv.__code__.co_consts and 300 not in inv.__code__.co_consts)
        # Runtime scheduling is exercised serially, using no worker pool or child.
        class SerialPool:
            def __init__(self,**kw):ck('mock_requested_workers12',kw=={'max_workers':12})
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def map(self,fn,seq):return [fn(x) for x in seq]
        def forbidden(*a,**kw):raise AssertionError('No subprocess allowed in source controls')
        subprocess.run=forbidden;R.context=lambda c,runtime=False:(m,h,ev);R.concurrent.futures.ThreadPoolExecutor=SerialPool
        C.operator_checks=lambda *a:{};C.field_checks=lambda *a:{}
        def parse(rec,j,row,*a):
            rec['value']={'energy':0.,'gradient':[[0.,0.,0.] for _ in range(row['atoms'])]};rec['moments']={'energy_kcal_mol':0.};rec['response_readback']={'native_solve_blocks':1+6*row['atoms'] if j['fd'] else 1}
        C.parse=parse;os.environ['SLURM_JOB_ID']='MOCK_ONLY'
        c=copy.deepcopy(cfg);c.update(executables={k:{'path':'MOCK_NOT_EXECUTABLE'} for k in ('analyze','testgrad')},mock_only=True)
        def invoke_ok(argv,folder,output,receipts):output.write_text('MOCK ONLY\n');receipts.append({'mock_only':True,'returncode':0})
        R.invocation=lambda _:invoke_ok;p=out/'mock_success';p.mkdir();success=R.run(c,p)
        ck('mock_complete49_and132',success['native_program_invocation_attempts']==success['native_programs_returned_zero']==success['native_programs_completed']==49 and success['internal_FD_energies_observed']==success['internal_FD_energies_completed']==132)
        def invoke_fail(argv,folder,output,receipts):output.write_text('MOCK FAILURE\n');receipts.append({'mock_only':True,'returncode':1});raise RuntimeError('declared synthetic failure')
        R.invocation=lambda _:invoke_fail;p=out/'mock_failure';p.mkdir();failure=R.run(c,p)
        ck('mock_stop_after_first_failure',failure['native_program_invocation_attempts']==1 and failure['native_programs_returned_zero']==0 and failure['global_stop'])
        ck('mock_48_untouched_slots',sum(v['status']=='not_attempted' for v in failure['rows'].values())==48)
        ck('failure_preserved',failure['status']=='radial_drive_numerical_qualification_failed' and failure['error']=="RuntimeError('declared synthetic failure')")
        ck('cache_bridge_not_auto_admitted',success['cache_bridge_admitted'] is False and success['physical_accuracy_qualified'] is False)
        result['status']='passed_radial_qualification_policy_controls'
    except Exception as e:result.update(status='failed',error=repr(e),traceback=traceback.format_exc())
    finally:
        result['check_count']=len(checks);result['cost']={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1,'mock_workers_executed_serially':True};P.save(out/'RESULT.json',result);print(json.dumps({k:result[k] for k in ('status','check_count','cost')}))
    if result['status']=='failed':raise SystemExit(1)
if __name__=='__main__':main()
