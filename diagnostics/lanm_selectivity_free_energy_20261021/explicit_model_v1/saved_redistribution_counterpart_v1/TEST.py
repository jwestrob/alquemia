#!/usr/bin/env python3
"""One-thread saved/explicit-synthetic controls; never invokes native chemistry."""
import argparse
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('moments_test_adapter',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def run(config,preparation,out):
    start=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config).read_text());prep=json.loads(Path(preparation).read_text());manifest=m.read(prep['manifest']);checks=[]
    def ok(name,value):
        if not value:raise AssertionError(name)
        checks.append(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,KeyError,TypeError,IndexError):checks.append(name);return
        raise AssertionError('Accepted negative control: '+name)
    m.bound_inputs(cfg,manifest);helper=m.module(cfg['induction_helper'],'saved_moment_induction')
    ok('exact78',len(manifest['rows'])==78);ok('old38_energy_authorities',sum(r['energy_replay_authority'] is not None for r in manifest['rows'].values())==38)
    ok('all12_cancellation_proofs',all(p['permanent_cancellation_exact'] and set(p['permanent_atom_coefficients'])=={0} for g in manifest['permanent_cancellation'].values() for p in g.values()))
    ok('kind_inventory',{k:sum(r['kind']==k for r in manifest['rows'].values()) for k in ('cluster','dimer','water','shell','ion')}=={'cluster':4,'dimer':34,'water':34,'shell':4,'ion':2})
    counter=m.Counter()
    for i in range(78):counter.reserve()
    rejects('program79_rejected',counter.reserve);counter=m.Counter();counter.stop.set();rejects('global_stop_rejected',counter.reserve)
    rejects('unreleased_execution_rejected',lambda:m.release_gate(cfg,prep,False))
    bad=copy.deepcopy(manifest['rows']);bad[m.GIDS[0]+'__water_01']['frame_signatures'][0]['axis_parent_indices'].reverse();rejects('axis_order_change_rejected',lambda:m.cancellation(bad))
    bad=copy.deepcopy(manifest['rows']);bad[m.GIDS[0]+'__dimer_01']['frame_signatures'][1]['permanent_record']['components'][0]+=.1;rejects('permanent_charge_change_rejected',lambda:m.cancellation(bad))
    rejects('cut_H_bond_rejected',lambda:m.extract([['1','O','0','0','0','36','2'],['2','H','1','0','0','37','1']],[0],'bad'))
    for token in ('external-field 0 0 0','exchpol 34 1 1 4 1','active 1','mass 36 1'):
        rejects('forbidden_'+token.split()[0],lambda t=token:m.render('parameters amoeba09.prm\n'+t,'polarization MUTUAL'))
    ok('prose_apostrophe_guard',m.render('parameters amoeba09.prm',"A water molecule's parameters\npolarization MUTUAL").endswith('save-uinduce\n'))
    results={};fixture_paths={}
    for cid,row in manifest['rows'].items():
        folder=out/'explicit_synthetic_fixtures'/cid;folder.mkdir(parents=True)
        v=[]
        for a in row['parent_atom_indices']:
            v.append([0.,0.,0.] if row['kind'] in ('water','ion') else [round(.00001*(a+1)*(j+1)*(1 if row['kind']=='cluster' else -1),6) for j in range(3)])
        induced=[math.fsum(x[j] for x in v) for j in range(3)]
        # Arbitrary fixed per-parent-atom permanent vectors cancel symbolically.
        perm=[math.fsum((a*.002*(j+1)) for a in row['parent_atom_indices']) for j in range(3)]
        mu=[perm[j]+induced[j]-m.DEBYE*row['charge']*row['COM_A'][j] for j in range(3)]
        energy=row['energy_replay_kcal_mol'] if row['energy_replay_kcal_mol'] is not None else .12345678
        text='EXPLICIT SYNTHETIC CONTROL, NOT A MOLECULAR RESULT\nTinker Version 26.2\n'
        text+='Induced Dipoles : Iterations 2 RMS Residual 0.0000000010\n'
        text+=f'Total Potential Energy : {energy:.8f} Kcal/mole\nEnergy Component Breakdown :\nPolarization {energy:.8f} 1\n'
        text+='Induced Dipoles : Iterations 2 RMS Residual 0.0000000000\n'
        text+=f'Total Electric Charge : {row["charge"]:.5f} Electrons\nDipole X,Y,Z-Components : '+' '.join(f'{x:.3f}' for x in mu)+'\n'
        a=folder/'analyze.out';u=folder/'system.001u';a.write_text(text);u.write_text(str(row['atoms'])+' EXPLICIT SYNTHETIC\n'+'\n'.join(f'{i+1} {name} '+' '.join(f'{x:.6f}' for x in vec) for i,(name,vec) in enumerate(zip(row['names'],v)))+'\n')
        parsed=m.parse_moments(a,u,row,helper.induction);ok('synthetic_complete_'+cid,parsed['induced_Debye']==v and parsed['charge_e']==row['charge']);results[cid]={'status':'complete','moments':parsed};fixture_paths[cid]=(a,u)
    cyc=m.cycles(manifest,results);ok('all_M_origin_precision_checks',cyc['all_M_consistency_checks_passed']);ok('four_parent_cycles',len(cyc['parents'])==4)
    gid=m.GIDS[0];a,u=fixture_paths[gid];original=a.read_text();uoriginal=u.read_text();row=manifest['rows'][gid]
    for name,badtext in [('NaN_residual',original.replace('0.0000000010','NaN')),('starred_residual',original.replace('0.0000000010','********')),('negative_residual',original.replace('0.0000000010','-0.0000000010')),('large_residual',original.replace('0.0000000010','0.0000000101')),('missing_solution',original.replace('Induced Dipoles : Iterations 2 RMS Residual 0.0000000000\n','')),('wrong_charge',original.replace('3.00000 Electrons','2.00000 Electrons'))]:
        a.write_text(badtext);rejects(name,lambda:m.parse_moments(a,u,row,helper.induction))
    a.write_text(original)
    for name,badtext in [('missing_site','\n'.join(uoriginal.splitlines()[:-1])+'\n'),('duplicate_site',uoriginal+'1 La 0 0 0\n'),('nonfinite_site',uoriginal.replace('0.000010','NaN',1)),('site_name',uoriginal.replace(row['names'][0],'XX',1))]:
        u.write_text(badtext);rejects(name,lambda:m.parse_moments(a,u,row,helper.induction))
    u.write_text(uoriginal)
    badrow=copy.deepcopy(row);badrow['energy_replay_kcal_mol']+=.01;rejects('old_energy_mismatch',lambda:m.parse_moments(a,u,badrow,helper.induction))
    missing=copy.deepcopy(results);del missing[gid];rejects('missing_cycle_row',lambda:m.cycles(manifest,missing))
    # Execute the actual driver failure path with a serial executor and a stub
    # that raises before starting any process. Gates are bypassed ONLY here.
    class SerialExecutor:
        def __init__(self,**kw):pass
        def __enter__(self):return self
        def __exit__(self,*a):return False
        def map(self,fn,items):return map(fn,items)
    admissions=[]
    def fail(*args):admissions.append(args[0]);raise ValueError('EXPLICIT SYNTHETIC pre-process native failure')
    oldfuncs={n:getattr(m,n) for n in ('checked_config','bound_inputs','release_gate','read','module')};oldpool=m.concurrent.futures.ThreadPoolExecutor;oldjob=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='synthetic_no_native'
    testcfg=copy.deepcopy(cfg);testcfg['preparation_result']={'synthetic':'preparation'}
    try:
        m.checked_config=lambda *a,**k:None;m.bound_inputs=lambda *a:manifest;m.release_gate=lambda *a,**k:None;m.read=lambda p:prep if p==testcfg['preparation_result'] else manifest
        m.module=lambda p,n:SimpleNamespace(invoke=fail) if p==cfg['native_helper'] else helper;m.concurrent.futures.ThreadPoolExecutor=SerialExecutor
        dest=out/'explicit_mock_failure';dest.mkdir();r=m.run(testcfg,dest)
        ok('actual_driver_first_failure_stop',len(admissions)==1 and r['native_program_invocations']==1 and r['global_stop'] and r['status']=='failed_or_unavailable')
        ok('actual_driver_full78_denominator',len(r['rows'])==78 and sum(x['status']=='not_attempted' for x in r['rows'].values())==77)
    finally:
        for name,fn in oldfuncs.items():setattr(m,name,fn)
        m.concurrent.futures.ThreadPoolExecutor=oldpool
        if oldjob is None:os.environ.pop('SLURM_JOB_ID',None)
        else:os.environ['SLURM_JOB_ID']=oldjob
    result={'status':'passed','implementation':m.pin(HERE/'RUN.py'),'test_implementation':m.pin(__file__),'specification':m.pin(HERE/'SPEC.md'),'config':m.pin(config),'preparation':m.pin(preparation),'manifest':prep['manifest'],
        'checks_passed':len(checks),'checks':checks,'new_native_calls':0,'synthetic_invocations':len(admissions),'synthetic_processes_started':0,'max_local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    m.save(out/'RESULT.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=run(a.config,a.preparation,out)
    except Exception as exc:m.save(out/'FAILURE.json',{'status':'failed','implementation':m.pin(__file__),'error':repr(exc),'new_native_calls':0});raise
    print(json.dumps({'status':r['status'],'checks_passed':r['checks_passed'],'new_native_calls':0}))
