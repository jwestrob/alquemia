#!/usr/bin/env python3
"""One foreground AF3 invocation per declared target; Q9 scout precedes other targets."""
from pathlib import Path
import argparse,json,os,resource,subprocess,time,traceback
from PREPARE import read,pin,verify,write,require_slurm,command_prefix,launcher_identity

def run(manifest,target):
    allocated=int(os.environ.get('SLURM_CPUS_PER_TASK','0'))
    affinity=sorted(os.sched_getaffinity(0))
    if allocated!=12 or len(affinity)<allocated:raise RuntimeError('twelve allocated CPU slots required')
    if len(affinity)>allocated:
        os.sched_setaffinity(0,affinity[:allocated]);affinity=sorted(os.sched_getaffinity(0))
    require_slurm(allocated);m=read(manifest)
    for rec in m['implementation'].values():verify(rec)
    for key in ('container','weights','python'):verify(m[key])
    launcher=launcher_identity(m)
    if not os.environ.get('CUDA_VISIBLE_DEVICES') and not os.environ.get('SLURM_JOB_GPUS'):
        raise RuntimeError('declared GPU allocation required')
    assert m['model_seeds']==[1] and m['num_diffusion_samples']==5 and m['num_recycles']==10
    prep=read(Path(m['prepared_output'])/'PREPARED.json')
    assert prep['status']=='prepared' and prep['manifest']==pin(manifest)
    chosen=next(t for t in prep['targets'] if t['accession']==target)
    if target!=m['scout_accession']:
        scout_path=Path(m['fold_output'])/m['scout_accession']/'EXECUTION.json';scout=read(scout_path)
        if scout['status']!='execution_complete_not_structurally_validated' or len(scout['models'])!=5:
            raise RuntimeError('complete real five-model Q9 executable scout required')
        if scout['manifest']!=pin(manifest):raise RuntimeError('scout manifest differs')
        for model in scout['models']:verify(model['model'])
    verify(chosen['input']);out=Path(m['fold_output'])/target;out.mkdir(parents=True,exist_ok=False)
    command=command_prefix(m,True)+['python','/app/alphafold/run_alphafold.py',
             '--json_path='+chosen['input']['path'],'--output_dir='+str(out/'outputs'),
             '--model_dir='+str(Path(m['weights']['path']).parent),'--run_data_pipeline=false',
             '--run_inference=true','--num_diffusion_samples=5','--num_recycles=10']
    receipt={'manifest':pin(manifest),'prepared':pin(Path(m['prepared_output'])/'PREPARED.json'),
             'container_launcher':launcher,'target':target,'input':chosen['input'],'command':command,'job_id':os.environ['SLURM_JOB_ID'],
             'status':'running','models':[],'new_scores':0,'structural_admission':'not_yet_evaluated',
             'CPU_affinity':affinity,'allocated_CPUs':allocated,'numerical_thread_limit':allocated,
             'SLURM_resources':{k:v for k,v in os.environ.items() if k in ('SLURM_JOB_NODELIST','SLURM_JOB_CPUS_PER_NODE','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE','SLURM_MEM_PER_CPU','SLURM_JOB_GPUS','CUDA_VISIBLE_DEVICES')}}
    write(out/'START.json',receipt);start=time.monotonic()
    try:
        with (out/'af3.log').open('x') as log:
            process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        receipt['returncode']=process.returncode
        for i in range(5):
            matches=list((out/'outputs').glob('*/seed-1_sample-'+str(i)+'/*model.cif'))
            if len(matches)==1:
                receipt['models'].append({'sample':i,'model':pin(matches[0]),'new_source_id':chosen['new_source_ids'][i],
                                           'original_source_id':chosen['original_source_ids'][i],'state_id':chosen['new_state_id']})
            else:
                receipt.setdefault('unavailable_samples',[]).append({'sample':i,'model_matches':len(matches)})
        receipt['status']='execution_complete_not_structurally_validated' if process.returncode==0 and len(receipt['models'])==5 else 'failed_or_partial'
    except Exception as exc:receipt.update(status='failed_or_partial',reason=str(exc),traceback=traceback.format_exc())
    finally:
        child=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(wall_seconds=time.monotonic()-start,child_CPU_seconds=child.ru_utime+child.ru_stime,
                       child_peak_RSS_KiB=child.ru_maxrss,log=pin(out/'af3.log') if (out/'af3.log').exists() else None)
        write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='execution_complete_not_structurally_validated' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--target',choices=('q9z4j7','c5axv8','q88jh5'),required=True)
    a=p.parse_args();raise SystemExit(run(a.manifest,a.target))
