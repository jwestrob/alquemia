#!/usr/bin/env python3
"""Bounded metadata, independent scalar/tensor arithmetic and failure controls; no native calls."""
import argparse
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import types

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('exch_control_subject',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True);a=p.parse_args()
out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def rejects(name,fn):
    try:fn()
    except (ValueError,KeyError,AssertionError):checks.append(name);return
    raise AssertionError(name+' admitted unexpectedly')
cfg=json.loads(Path(a.config).read_text());prep=m.read(m.pin(a.preparation));manifest=m.read(prep['manifest']);data=m.read(cfg['data'])
try:
    check('actual zero-call preparation complete',prep['status']=='prepared' and prep['new_native_calls']==0)
    check('final source/config joined',prep['implementation']==m.pin(HERE/'RUN.py') and prep['config']==m.pin(a.config))
    check('actual92/116/636 inventory',len(manifest['rows'])==92 and sum(r['native_programs'] for r in manifest['rows'].values())==116 and sum(r['internal_FD_energies'] for r in manifest['rows'].values())==636)
    check('stage counts4/4/84',list(map(lambda s:len(manifest['stages'][s]),('controls','derivatives','matrix')))==[4,4,84])
    check('FD only q05 consumedclusters',set(manifest['stages']['derivatives'])=={'q05__'+g for g in m.CLUSTERS})
    for q in m.Q:
        rows=[r for r in manifest['rows'].values() if r['kind'] in ('fixture','pair') and r['q']==q]
        check('same amplitude and44 literal cells '+str(q),len(rows)==44 and sum(r['kind']=='pair' for r in rows)==34 and sum(r['kind']=='fixture' for r in rows)==10)
    resolver=m.module(cfg['resolver'],'exch_test_resolver');pair=m.module(cfg['pair_adapter'],'exch_test_pair')
    # Independent native overlap and its derivative, not a force-field evaluation.
    for r in (2.0,2.5,3.0,4.0):
        t=2*r;s=(1+t+t*t/3)*math.exp(-t);expected=s*s
        check('native overlap at '+str(r),m.overlap(r)==expected)
        ds=-4/3*s*(t+t*t)*math.exp(-t);h=1e-5
        check('overlap derivative finite difference '+str(r),abs((m.overlap(r+h)-m.overlap(r-h))/(2*h)-ds)<1e-10)
    for q in m.Q:
        check('normalized spring '+str(q),abs(manifest['overlap']['O_springs'][m.label(q)]*m.overlap(2.5)-q)<1e-16)
    for cid,row in manifest['rows'].items():
        text=m.verify(row['input']['key']).read_text();base=m.verify(row['source_input']['key']).read_text()
        check('unchanged original prefix '+cid,text.startswith(base+'\n'))
        check('positive full tensors '+cid,row['tensor_control']['positive_definite'] and all(min(s['eigenvalues'])>=1 for s in row['tensor_control']['sites']))
        if row['kind']=='shell':check('metal absent identity '+cid,all(s['matrix']==[[float(i==j) for j in range(3)] for i in range(3)] for s in row['tensor_control']['sites']))
    fixture=manifest['rows']['q05__aquo9_central_La'];sites=fixture['tensor_control']['sites']
    check('nonzero offdiagonal tested',max(s['maximum_absolute_off_diagonal'] for s in sites)>0.01)
    oxygen=[s for s in sites if s['type']==36]
    check('radial eigenvalue1.05',all(abs(s['eigenvalues'][2]-1.05)<1e-12 for s in oxygen))
    base=m.verify(fixture['source_input']['key']).read_text();prm=m.verify(fixture['source_input']['parameters']).read_text();xyz=m.verify(fixture['source_input']['xyz']).read_text()
    rejects('negative amplitude',lambda:m.render(base,-.01));rejects('undeclared amplitude',lambda:m.render(base,.02));rejects('nonfinite amplitude',lambda:m.render(base,float('nan')))
    rejects('existing term collision',lambda:m.render(base+'\nexchpol 34 1 1 4 0\n',.05))
    rejects('zero alpha current site',lambda:m.tensors(xyz,.05,pair,resolver,prm,base+'\npolarize 37 0 0.39 36\n'))
    rejects('incorrect class map',lambda:m.tensors(xyz,.05,pair,resolver,prm,base+'\natom 36 73 O "wrong class" 8 15.999 2\n'))
    # Synthetic numerical receipt values test tolerances before molecular execution.
    rows={}
    for metal in ('La','Dy'):rows['zero_aquo9_'+metal]={'status':'complete','model':copy.deepcopy(data['baseline_model']['aquo9_central_'+metal])}
    shell={'energy':1.,'forces':[[0.,0.,0.]]*27}
    rows.update(shell_zero={'status':'complete','model':copy.deepcopy(shell)},shell_q05={'status':'complete','model':copy.deepcopy(shell)})
    check('exact descriptive flags',all(c['printed_energy_exact'] and c['printed_forces_exact'] for c in m.prechecks(rows,data)['comparisons'].values()))
    near=copy.deepcopy(rows);near['shell_q05']['model']['energy']+=5e-5;near['shell_q05']['model']['forces'][0][0]=5e-4
    check('prospective inherited tolerances',m.prechecks(near,data)['passed'])
    far=copy.deepcopy(rows);far['shell_q05']['model']['energy']+=2e-4;rejects('energy replay failure',lambda:m.prechecks(far,data))
    far=copy.deepcopy(rows);far['shell_q05']['model']['forces'][0][0]=.002;rejects('force replay failure',lambda:m.prechecks(far,data))
    missing=copy.deepcopy(rows);missing['shell_q05']['status']='not_attempted';rejects('missing precheck not pass',lambda:m.prechecks(missing,data))
    nan=copy.deepcopy(rows);nan['shell_q05']['model']['energy']=float('nan');rejects('nonfinite replay not pass',lambda:m.prechecks(nan,data))
    f=out/'synthetic_induction.out';f.write_text(' Induced Dipoles : Iterations 7 RMS Residual 0.0000000001\n Induced Dipoles : Iterations 8 RMS Residual 0.0000000100\n')
    check('all induction summaries admitted',m.induction(f)['solutions']==2)
    f.write_text(' Induced Dipoles : Iterations 8 RMS Residual 0.0000000200\n');rejects('residual above tolerance',lambda:m.induction(f))
    f.write_text(' Induced Dipoles : Iterations 8 RMS Residual -0.0000000001\n');rejects('negative residual',lambda:m.induction(f))
    f.write_text(' Induced Dipoles : Iterations 8 RMS Residual NaN\n');rejects('nonfinite residual',lambda:m.induction(f))
    f.write_text(' Induced Dipoles : Iterations 7 RMS Residual 0.0000000001\n Induced Dipoles : Iterations 8 RMS Residual ********\n');rejects('mixed valid/starred residual',lambda:m.induction(f))
    f.write_text(' Induced Dipoles : Iterations 7 RMS Residual 0.0000000001\n Induced Dipoles : Iterations 8 RMS Residual NaN\n');rejects('mixed valid/nonfinite residual',lambda:m.induction(f))
    c=m.Counter()
    for _ in range(116):c.reserve()
    rejects('hard116 program ceiling',c.reserve)
    c=m.Counter();c.reserve(636);rejects('hard636 internal energy ceiling',lambda:c.reserve(1))
    c=m.Counter();c.stop.set();rejects('global stop blocks next program',c.reserve)
    rejects('unreleased execution gate',lambda:m.release_gate(cfg,prep,manifest,runtime=False))
    # Exercise the actual run worker, with no executable and an inline single-thread pool.
    # The first analyze exits zero and prints a valid induction residual but has fatal energy output.
    invoked=[]
    def fake_invoke(argv,cwd,output,receipts):
        invoked.append(argv);Path(output).write_text('Tinker Version 26.2\n Induced Dipoles : Iterations 7 RMS Residual 0.0000000001\n FATAL synthetic analyze energy failure\n');receipts.append({'synthetic':True,'argv':argv,'returncode':0})
    class InlinePool:
        def __init__(self,**kw):pass
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def map(self,fn,items):return map(fn,items)
    original_module=m.module;original_gate=m.release_gate;original_pool=m.concurrent.futures.ThreadPoolExecutor;job=os.environ.get('SLURM_JOB_ID')
    try:
        m.module=lambda pin,name:types.SimpleNamespace(invoke=fake_invoke) if pin==cfg['native_helper'] else original_module(pin,name)
        m.release_gate=lambda *args,**kwargs:None;m.concurrent.futures.ThreadPoolExecutor=InlinePool;os.environ['SLURM_JOB_ID']='0'
        rcfg=copy.deepcopy(cfg);rcfg['preparation_result']=m.pin(a.preparation);fakeout=out/'synthetic_analyze_failure';fakeout.mkdir();failed=m.run(rcfg,fakeout)
    finally:
        m.module=original_module;m.release_gate=original_gate;m.concurrent.futures.ThreadPoolExecutor=original_pool
        if job is None:os.environ.pop('SLURM_JOB_ID',None)
        else:os.environ['SLURM_JOB_ID']=job
    check('analyze energy failure blocks testgrad',len(invoked)==1 and invoked[0][0]==cfg['executables']['analyze']['path'])
    check('first failure stops later stages',failed['status']=='failed_or_unavailable' and failed['global_stop'] and failed['summaries'] is None and len(failed['rows'])==92 and sum(r['status']=='not_attempted' for r in failed['rows'].values())==91)
    report={'status':'passed','implementation':m.pin(HERE/'TEST.py'),'subject':m.pin(HERE/'RUN.py'),'config':m.pin(a.config),'preparation':m.pin(a.preparation),'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'molecular_calculations':False,'threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
except Exception as exc:
    report={'status':'failed','error':repr(exc),'checks_passed':len(checks),'checks':checks,'new_native_calls':0,'subject':m.pin(HERE/'RUN.py')}
m.save(out/'RESULT.json',report);print(json.dumps({k:report[k] for k in ('status','checks_passed','new_native_calls')}));raise SystemExit(0 if report['status']=='passed' else 1)
