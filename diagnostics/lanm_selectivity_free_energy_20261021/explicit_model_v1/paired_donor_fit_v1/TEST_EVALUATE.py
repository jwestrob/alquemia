#!/usr/bin/env python3
"""Bounded saved-only baseline/render/cache/ceiling checks. Never native execution."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import time

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    path=Path(__file__).with_name('EVALUATE.py');s=importlib.util.spec_from_file_location('donor_eval_controls',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    result={'status':'running','implementation':m.pin(path),'test_source':m.pin(__file__),'checks':{},'native_invocations':0,'local_threads':1,'error':None}
    (out/'EVALUATE.py').write_bytes(path.read_bytes())
    def check(label,yes):
        result['checks'][label]=bool(yes);m.require(yes,label)
    def rejects(label,fn,needle):
        try:fn()
        except Exception as exc:check(label,needle in str(exc));return
        check(label,False)
    try:
        check('one_thread_controls',all(os.environ.get(k)=='1' for k in THREADS))
        cfg=m.default_config();m.save(out/'CONFIG.json',cfg);ev=m.Evaluator(cfg,out/'cache')
        def forbidden(*a,**k):
            result['native_invocations']+=1;raise AssertionError('Forbidden native invocation')
        ev.native.invoke=forbidden
        check('original28_data_unchanged',ev.data==m.read(cfg['training_result']))
        transfer=m.read(cfg['transfer_result']);initial=ev.evaluate(ev.initial_vector)
        check('complete28',set(initial)==set(m.CELLS))
        for gid in m.CELLS:
            check('baseline_'+gid,all(initial[gid][k]==transfer['rows'][gid]['model'][k] for k in ('energy','forces','native')))
            check('baseline_key_'+gid,ev.cache_identity(ev.initial_vector,gid)[2]==m.verify(cfg['baseline_inputs'][gid]['key']).read_text())
        check('all8_water_admitted',ev.water_invariants['passed'] and set(ev.water_invariants['records'])==set(m.WATER_CELLS) and all(r['unchanged'] for r in ev.water_invariants['records'].values()))
        check('water_center_baselines',set(ev.water_baseline)==set(m.WATER_CENTERS) and all(all(ev.water_baseline[g][k]==initial[g][k] for k in ('energy','forces','native')) for g in m.WATER_CENTERS))
        check('exact_six_central_FD_metadata',sum(initial[g]['native']['finite_difference'] is not None for g in m.CELLS)==6)
        changed_inventory={}
        for j in range(8):
            point=ev.initial_vector[:];point[j]+=.005*(ev.family['parameters'][j]['upper']-ev.family['parameters'][j]['lower'])
            changed=[]
            for gid in m.CELLS:
                token,identity,key,render=ev.cache_identity(point,gid);old=ev.cache_identity(ev.initial_vector,gid)
                family,index,cl,typ=ev.recipient(gid);affected=index is not None and index<=j<index+2
                check('active_identity_'+str(j)+'_'+gid,(token!=old[0])==affected)
                check('active_key_'+str(j)+'_'+gid,(key!=old[2])==affected)
                check('other_pairs_'+str(j)+'_'+gid,render['fixed_water_and_other_parameters_unchanged'])
                if affected:
                    changed.append(gid);check('two_effective_values_'+str(j)+'_'+gid,len(identity['effective_values'])==2)
                else:
                    value=ev._cell(point,gid)
                    check('unaffected_native_reuse_'+str(j)+'_'+gid,all(value[k]==initial[gid][k] for k in ('energy','forces','native')) and value['provenance']['reused_baseline'])
            changed_inventory[str(j)]=changed
            check('affected_denominator_'+str(j),len(changed)==(8 if j in (0,1,4,5) else 4))
        result['changed_cell_inventory']=changed_inventory
        check('zero_initial_and_unaffected_calls',ev.snapshot()['native_program_invocations']==0 and result['native_invocations']==0)
        check('one_initial_vector',ev.snapshot()['distinct_parameter_vectors']==1)
        for j in range(8):
            bad=ev.initial_vector[:];bad[j]=ev.family['parameters'][j]['upper']+1e-6
            rejects('bound_'+str(j),lambda bad=bad:ev.checked_values(bad),'outside frozen bounds')
        rejects('wrong_dimension',lambda:ev.checked_values([1]*7),'Eight parameters')
        bad=ev.initial_vector[:];bad[0]=float('nan');rejects('nonfinite',lambda:ev.checked_values(bad),'outside frozen bounds')
        bad[0]=True;rejects('boolean',lambda:ev.checked_values(bad),'outside frozen bounds')
        point=ev.initial_vector[:];point[0]+=.001
        saved=copy.deepcopy(ev._ledger);ev._ledger['vector_keys']=[str(i) for i in range(481)]
        rejects('481_total_vector_ceiling',lambda:ev.evaluate(point),'Parameter-vector ceiling reached');ev._ledger=saved
        release={'status':'approved_paired_donor_development','evaluator':cfg['implementation'],'evaluator_config_sha256':ev.fingerprint,
                 'max_parameter_vectors':481,'max_native_program_invocations':9216,'validation_scoring_authorized':False,'physical_accuracy_claimed':False}
        m.save(out/'SYNTHETIC_GUARD_RELEASE.json',release);ev.root_release=m.pin(out/'SYNTHETIC_GUARD_RELEASE.json')
        previous={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_NTASKS','SLURM_CPUS_PER_TASK')}
        try:
            os.environ.update(SLURM_JOB_ID='0',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='8')
            ev._ledger['native_program_invocations']=9216
            rejects('9216_native_ceiling_before_invoke',lambda:ev._cell(point,'hydrated_acetate_mono_minus_La'),'Native call ceiling reached')
            check('cap_failure_global_stop',ev.stop_event.is_set() and ev._ledger['stopped'])
            check('cap_failure_no_native_calls',result['native_invocations']==0)
        finally:
            for k,v in previous.items():
                if v is None:os.environ.pop(k,None)
                else:os.environ[k]=v
        result['guard_ledger_is_synthetic_not_execution']=True
        result['status']='passed'
    except Exception as exc:result.update(status='failed',error=repr(exc))
    finally:
        result.update(elapsed_seconds=time.monotonic()-start,parent_CPU_seconds=time.process_time()-cpu,check_count=len(result['checks']))
        m.save(out/'RESULT.json',result)
    print(json.dumps({k:result[k] for k in ('status','error','check_count','native_invocations','elapsed_seconds','parent_CPU_seconds')}))
    return result['status']=='passed'

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',required=True);args=a.parse_args();raise SystemExit(0 if main(args.output) else 1)
