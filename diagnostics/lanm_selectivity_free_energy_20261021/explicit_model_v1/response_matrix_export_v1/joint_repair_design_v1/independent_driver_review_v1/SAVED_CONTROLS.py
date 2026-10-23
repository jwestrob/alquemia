#!/usr/bin/env python3
"""Three synthetic driver paths; no real evaluator, chemistry or fitted data."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
s=importlib.util.spec_from_file_location('driver',ROOT/'REFINE.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks=[]
def check(name,test):
    if not test:raise AssertionError(name)
    checks.append(name)
fake=out/'SYNTHETIC_EVALUATOR.py'
fake.write_text('''import threading
class Evaluator:
    def __init__(self,cfg,cache,root_release=None):
        self.cfg=cfg;self.initial_vector=[.5]*12;self.family={'parameters':cfg['parameters']};self.stop_event=threading.Event();self.count=0;self.lock=threading.Lock()
    def snapshot(self):return {'native_program_invocations':0,'synthetic_evaluations':self.count}
    def evaluate(self,p):
        with self.lock:self.count+=1
        if self.stop_event.is_set():raise RuntimeError('synthetic stopped')
        if self.cfg['mode']=='failure' and p[0]>.5:raise RuntimeError('intentional synthetic vector failure')
        v=(p[0]-.505) if self.cfg['mode']=='improve' else 1.
        r=[v]+[0.]*530;j=v*v
        return {'score':{'weighted_residuals':r,'objective':j,'block_objectives':{'synthetic':j},'original_EF_checks_passed':True},'synthetic':True}
''')
parameters=[dict(name='synthetic_'+str(i),lower=0.,upper=1.,initial=.5)for i in range(12)]
domain=out/'SYNTHETIC_DOMAIN.json';d.save(domain,{'parameters':parameters})
os.environ.update({k:'1'for k in d.THREADS})
# Exercise only the allocation syntax guard; the fake evaluator cannot run programs.
os.environ.update(SLURM_JOB_ID='999999999',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='12')
summaries={}
for mode in ('flat','improve','failure'):
    folder=out/mode;folder.mkdir();ec=folder/'EVALUATOR.json';d.save(ec,{'parameters':parameters,'mode':mode})
    cfg={'implementation':d.pin(ROOT/'REFINE.py'),'protocol':d.pin(ROOT/'FIT_PROTOCOL.md'),'root_declaration':d.pin(ROOT/'ROOT_DECLARATION.md'),'parameter_domain':d.pin(domain),'directions':d.pin(ROOT/'DIRECTIONS.py'),'evaluator':d.pin(fake),'evaluator_config':d.pin(ec),'objective_ready':d.pin(domain),'limits':d.LIMITS,'root_release':None,'preparation_result':None}
    cp=folder/'CONFIG_PREPARE.json';d.save(cp,cfg)
    with contextlib.redirect_stdout(io.StringIO()):rc=d.run('prepare',str(cp),str(folder/'prepare'))
    prep=json.loads((folder/'prepare/RESULT.json').read_text());check(mode+': saved initial replay',rc==0 and prep['status']=='prepared_joint_response_fit' and prep['evaluator']['native_program_invocations']==0 and prep['new_parameter_vectors']==0)
    release=folder/'SYNTHETIC_RELEASE.json';d.save(release,{'status':'released_bounded_joint_response_fit','limits':d.LIMITS,'scientific_fingerprint':d.fingerprint(cfg),'synthetic':True});cfg.update(root_release=d.pin(release),preparation_result=d.pin(folder/'prepare/RESULT.json'));cp=folder/'CONFIG_RUN.json';d.save(cp,cfg)
    with contextlib.redirect_stdout(io.StringIO()):rc=d.run('run',str(cp),str(folder/'run'))
    result=json.loads((folder/'run/RESULT.json').read_text());summaries[mode]={'result':d.pin(folder/'run/RESULT.json'),'status':result['status'],'new_parameter_vectors':result['new_parameter_vectors'],'iterations':len(result['iterations'])}
    check(mode+': all exact full vector files retained',all(d.read(v['artifact'])['parameters']==v['parameters'] for v in result['vectors'].values()))
    check(mode+': no native programs',result['evaluator']['native_program_invocations']==0)
    if mode=='failure':
        check('failure: fail closed without selected candidate',rc==1 and result['status']=='failed_or_unavailable' and result['candidate_frozen'] is False and 'selected'not in result)
        check('failure: unavailable vectors retained',any(v['status']=='unavailable' for v in result['vectors'].values()))
        check('failure: no proposal batch after derivative failure',len(result['iterations'])==1 and 'trial_points'not in result['iterations'][0])
    else:
        check(mode+': complete declared proposal batches',all(len(x['trial_points'])==len(x['trial_keys'])==4 and x['status']=='complete'for x in result['iterations']))
        check(mode+': no EF-only early stop',len(result['iterations'])>=1 and result['original_EF_checks_passed'] is True)
        check(mode+': deterministic no-improvement termination',result['stopping_reason']=='no_actual_improvement_at_complete_trial_batch')
        if mode=='flat':check('flat: current wins exact ties',result['selected_key']==result['initial_key'] and len(result['iterations'])==1)
        else:
            key=result['iterations'][0]['selected_key'];check('improve: reused derivative must be declared as trial',result['vectors'][key]['role']=='derivative' and key in result['iterations'][0]['trial_keys'] and result['iterations'][0]['actual_improvement'])
            check('improve: strictly lower joint objective',result['selected_objective']<result['vectors'][result['initial_key']]['objective'] and len(result['iterations'])==2)
direction_control=Path(str(ROOT).replace('/diagnostics/','/workspaces/'))/'direction_controls_v1/RESULT.json'
previous=json.loads(direction_control.read_text());check('prior direction controls passed',previous['status']=='passed_direction_source_controls')
receipt={'status':'passed','implementation':d.pin(__file__),'driver':d.pin(ROOT/'REFINE.py'),'directions':d.pin(ROOT/'DIRECTIONS.py'),'protocol':d.pin(ROOT/'FIT_PROTOCOL.md'),'prior_direction_controls':d.pin(direction_control),'checks':checks,'check_count':len(checks),'synthetic_runs':summaries,'synthetic_only':True,'new_native_calls':0,'new_scientific_evaluations':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
d.save(out/'RESULT.json',receipt);print(json.dumps({'checks':len(checks),'result':d.pin(out/'RESULT.json')}))
