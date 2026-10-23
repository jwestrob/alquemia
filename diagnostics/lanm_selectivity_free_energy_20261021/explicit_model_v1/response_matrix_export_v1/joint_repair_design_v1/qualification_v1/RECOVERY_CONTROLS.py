#!/usr/bin/env python3
"""Original functions unchanged; mock only the finite recovery orchestration."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import threading
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('recovery',HERE/'RUN_v2.py');q=importlib.util.module_from_spec(s);s.loader.exec_module(q)
out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks=[]
def check(n,v):
    if not v:raise AssertionError(n)
    checks.append(n)
def funcs(p):return {n.name:ast.dump(n,include_attributes=False)for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
a,b=funcs(HERE/'RUN.py'),funcs(HERE/'RUN_v2.py')
for name in ('inventory','render','response_readback','regularized','matrix_checks','field_checks','gate','prepare'):
    check('unchanged '+name,a[name]==b[name])
check('only context/run changed among old functions',[n for n in a if a[n]!=b[n]]==['context','run'])
# Reproduce the exact iterator invalidation, then prove detachment removes it.
local={'status':'started'};published=local;it=iter(published);next(it);local['value']=1
try:next(it);raise AssertionError('Old alias failed to reproduce')
except RuntimeError:checks.append('old shared record iterator failure reproduced')
local={'status':'started'};published=copy.deepcopy(local);it=iter(published);next(it);local['value']=1
try:next(it)
except StopIteration:checks.append('detached published record remains stable')
W=Path(str(HERE).replace('/diagnostics/','/workspaces/'));m=q.read(q.pin(W/'recovery_preparation_v2/MANIFEST.json'));cfg=q.read(q.pin(W/'recovery_preparation_v2/CONFIG_RUN_TEMPLATE.json'));old=q.read(m['original_manifest'])
check('all 57 exact original inputs retained',m['programs']==old['programs'])
check('27 exact reused prefix',len(m['reused_rows'])==27 and all(v['status']=='complete'for v in m['reused_rows'].values()))
check('30 untouched future slots',len([j for j in m['programs']if j['id']not in m['reused_rows']])==30)
check('remaining FD168',sum(6*m['rows'][j['cell']]['atoms']for j in m['programs']if j['fd'] and j['id']not in m['reused_rows'])==168)
jobs={j['id']:j for j in m['programs']};h=q.helpers(m['sources']);counter=[];guard=threading.Lock();q.context=lambda cfg,runtime:(m,h)
def invoke(argv,cell,output,receipts):
    with guard:counter.append(cell.name)
    time.sleep(.002)
    output.write_text('SYNTHETIC ORCHESTRATION ONLY; NO NATIVE PROGRAM\n')
    (cell/'system.uind').write_text('synthetic\n')
    for i in (1,2):(cell/f'response_matrix_{i:03d}.dat').write_text('synthetic\n')
    receipts.append({'synthetic':True,'new_native_calls':0})
h['native'].invoke=invoke;h['induction'].induction=lambda output:{'synthetic':True}
h['gradient'].gradient=lambda output,n,fd,native:{'energy':0.,'fd_passed':True,'synthetic':True}
h['gradient'].program_admission=lambda *args:None
def energy(output):
    j=jobs[output.parent.name];g=j['cell'];center=m['reused_rows']['combined__'+g]['value'];e=center['energy']
    if 'q'in j:
        _,t=h['old_harness'].path_geometry(m['rows'][g]['positions_A'],j['kind'],0.)
        e+=j['q']*sum(x*y for a,b in zip(center['gradient'],t)for x,y in zip(a,b))
    return {'energy_kcal_mol':e,'synthetic':True}
h['moments'].parse_energy=energy;h['moments'].induced_file=lambda output:output.parent/'system.uind';h['moments'].parse_moments=lambda output,*args:energy(output)
q.response_readback=lambda *args:{'synthetic':True};q.matrix_checks=lambda *args:None;q.field_checks=lambda *args:{'synthetic':True,'passed':True}
os.environ['SLURM_JOB_ID']='999999999';r=q.run(cfg,out/'mock_execution') if False else None
(out/'mock_execution').mkdir();r=q.run(cfg,out/'mock_execution')
check('actual orchestration completed mocked remainder',r['status']=='local_response_numerical_qualification_complete')
check('exact30 mock starts, zero prefix repeats',len(counter)==len(set(counter))==30 and set(counter)==set(jobs)-set(m['reused_rows']))
check('new30 plus reused27 equals57',r['native_program_invocations']==30 and r['inherited_native_program_invocations']==27 and r['total_native_program_invocations']==57)
check('new168 plus inherited468 equals636',r['new_internal_FD_energies_completed']==168 and r['inherited_internal_FD_energies_completed']==468 and r['internal_FD_energies_completed']==636)
check('full57 published rows and no global stop',len(r['rows'])==57 and all(v['status']=='complete'for v in r['rows'].values()) and r['global_stop']is False)
check('published inherited records unchanged',all(r['rows'][k]==v for k,v in m['reused_rows'].items()))
q.save(out/'RESULT.json',{'status':'passed','implementation':q.pin(__file__),'recovery_source':q.pin(HERE/'RUN_v2.py'),'preparation':q.pin(W/'recovery_preparation_v2/MANIFEST.json'),'checks':checks,'check_count':len(checks),'mock_execution':q.pin(out/'mock_execution/RESULT.json'),'mock_native_starts':30,'real_native_calls':0,'new_scientific_points':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}});print(json.dumps(q.pin(out/'RESULT.json')))
