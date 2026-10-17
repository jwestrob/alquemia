#!/usr/bin/env python3
"""Saved-only actual failing-frame diagnosis; Slurm, no molecular calls."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def require(ok,message):
    if not ok:raise ValueError(message)


def read(record):
    require(pin(record['path'])['sha256']==record['sha256'],'Changed pin: '+record['path'])
    return json.loads(Path(record['path']).read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();cfg=json.loads(a.config.read_text())
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Diagnostic source differs')
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'Python executable differs')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    result={'status':'execution_failure','config':pin(a.config),'implementation':pin(__file__),
        'failed_preflight':cfg['failed_preflight'],'slurm_job_id':os.environ['SLURM_JOB_ID'],
        'new_molecular_calls':0,'model_initialized':False,'rows':[]}
    try:
        failed=read(cfg['failed_preflight']);parent_cfg=read(failed['config'])
        require(failed['slurm_job_id']=='12954' and failed['status']=='execution_failure'
                and failed['failure']=="ValueError('Exact ten-mode terminal lineage differs')"
                and failed['new_composite_evaluations']==0,'Exact precall failure required')
        require(pin(failed['implementation']['path'])['sha256']==failed['implementation']['sha256']
                ==parent_cfg['implementation']['sha256'],'Failed adapter source differs')
        spec=importlib.util.spec_from_file_location('unchanged_cartesian_diagnostic',failed['implementation']['path'])
        adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
        try:
            adapter.load_inputs(parent_cfg,out)
            raise ValueError('Original failure did not reproduce')
        except ValueError as exc:
            require(str(exc)=='Exact ten-mode terminal lineage differs','Different failure: '+str(exc))
            frames=[];tb=exc.__traceback__
            while tb:
                frame=tb.tb_frame
                if frame.f_code.co_name=='load_inputs' and 'old_starts' in frame.f_locals:
                    frames.append(frame)
                tb=tb.tb_next
            require(len(frames)==1,'Unique original ten-mode frame absent')
            f=frames[0];v=f.f_locals;np=v['np'] if 'np' in v else f.f_globals['np']
            result['original_exception']=repr(exc);result['original_traceback']=traceback.format_exc()
            result['failed_frame_source']=pin(f.f_code.co_filename)
            for sid,metal in f.f_globals['STATES']:
                old_start=v['old_starts'][(sid,metal,'candidate')];cid=f.f_globals['START_IDS'][(sid,metal)]
                slot=next(r for r in v['old']['searches'] if (r['source_id'],r['metal'],r['arm'])==(sid,metal,'candidate'))
                search=read(slot['search']);row=v['by'][(cid,metal)]
                checks={'slot_status':slot['status']=='admitted','candidate_id':slot['candidate_id']==cid,
                    'point_pin':slot['point']==search['candidate'],
                    'MACE_evidence_pin':row['MACE_evidence']['evidence_record']==slot['point'],
                    'candidate_admitted':search['candidate_admitted'],'final_geometry':search['final_geometry']['pass'],
                    'source_rebased_false':search['source_rebased'] is False,
                    'active_indices':search['active_indices']==old_start['active'],
                    'coordinate_scaling_T_exact':search['coordinate_scaling_T']==old_start['metric']['coordinate_scaling_T']}
                old=np.asarray(search['coordinate_scaling_T']);new=np.asarray(old_start['metric']['coordinate_scaling_T'])
                require(old.shape==new.shape==(10,10) and np.isfinite(old).all() and np.isfinite(new).all(),'Invalid compared T')
                _,_,jac,_=old_start['source']['kin'].evaluate(np.zeros(64))
                matrix=jac[old_start['active']].reshape(10,-1).T
                residuals={name:float(np.max(np.abs((matrix@T).T@(matrix@T)-np.eye(10)))) for name,T in [('archived',old),('recomputed',new)]}
                result['rows'].append({'source_id':sid,'metal':metal,'search':slot['search'],'checks':checks,
                    'archived_T':old.tolist(),'recomputed_T':new.tolist(),
                    'T_maximum_absolute_difference':float(np.max(np.abs(old-new))),
                    'T_maximum_absolute_value':float(max(np.max(np.abs(old)),np.max(np.abs(new)))),
                    'T_differing_entries':int(np.count_nonzero(old!=new)),
                    'original_whitening_maximum_residual':residuals,
                    'existing_whitening_tolerance':1e-10,
                    'recomputed_metric':old_start['metric']})
            result['status']='complete'
            result['runtime']={'python':sys.version,'numpy':np.__version__,
                'thread_environment':{k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')}}
    except Exception as exc:
        result.update(failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
