#!/usr/bin/env python3
"""Narrow saved/serialization controls; all native invocation is mocked."""
import argparse
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('donor_control_runner',HERE/'RUN.py');q=importlib.util.module_from_spec(s);s.loader.exec_module(q)

def run(manifest_path,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    m=json.loads(Path(manifest_path).read_text());h=q.helpers(m);checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    def rejects(name,call):
        try:call()
        except (ValueError,AssertionError):check(name,True)
        else:check(name,False)
    result={'status':'running','implementation':q.pin(__file__),'runner':q.pin(HERE/'RUN.py'),'matrix_adapter':q.pin(HERE/'MATRIX.py'),
            'manifest':q.pin(manifest_path),'checks':checks,'new_native_calls':0,'native_invocation_mocked':True}
    try:
        check('121 programs /576 requested FD /20 exports',len(m['programs'])==121 and sum(j['internal_FD_energies'] for j in m['programs'])==576 and sum(j['exports'] for j in m['programs'])==20)
        check('matrix donor qualifier has no equal-dual gates','equal_dual_RHS' not in (HERE/'MATRIX.py').read_text() and 'equal_dual_solution' not in (HERE/'MATRIX.py').read_text())
        check('actual native induced filename API',callable(h['moments'].induced_file))
        for cid in ('La_central','Dy_central','acetate7','ion_La','ion_Dy','aquo9_central_La__near__two_water','aquo9_central_La__dimer_01','aquo9_central_Dy__dimer_01'):
            job=next(j for j in m['programs'] if j['cell']==cid)
            r=q.metadata(m['rows'][cid],job)
            check(cid+' metadata site denominator',len(r['native_names'])==r['atoms'] and len(r['alpha_A3'])==r['atoms'])
            if cid=='acetate7':check('literal acetate charge -1',abs(r['charge']+1)<2e-12)
            if cid.endswith('central'):
                b=r['baseline'];parsed=h['gradient'].gradient(q.verify(b['outputs']['gradient_output']),r['atoms'],False,h['native'])
                check(cid+' old raw E/F replay',parsed['energy']==b['energy'] and parsed['forces']==b['forces'])
            if r.get('water_profile_active'):
                old=q.read(m['sources']['water_result'])['rows']['combined__'+cid]
                rb=h['water'].response_readback(q.verify(old['actual_output']),{'profile':'combined'},r)
                check(cid+' old nonzero water readback',rb['passed'] and not rb['bypass'])
        job=next(j for j in m['programs'] if j['cell']=='La_central' and j['profile']=='combined' and j['fd_total_separation_A']==.0001)
        row=q.metadata(m['rows']['La_central'],job);ex=q.regularized(row,'combined');n=row['atoms']
        def fmt(label,v):return label+' '+' '.join(format(float(x),'.17e') for x in v)+'\n'
        def block(value):
            text=''
            for mt in (500,501):text+=fmt('DONOR_RESPONSE_PARAMETERS',[mt,403,.1,.1])
            for pair in ex['pairs']:text+=fmt('DONOR_RESPONSE_PAIR',pair)
            for i in (3,4):text+=fmt('DONOR_RESPONSE_FIELD',[i+1,*ex['h'][i]])
            text+=fmt('DONOR_RESPONSE_PERMANENT_LAB_EA',[.3,-.2,.7])+fmt('DONOR_RESPONSE_CHARGE',[2])
            text+=fmt('LOCAL_RESPONSE_SPECTRUM',[.7,1.5,0,0,0])
            for i in range(n):text+=fmt('LOCAL_RESPONSE_INDUCED',[i+1,value,0,0,value+.01,0,0])
            text+=fmt('LOCAL_RESPONSE_DUAL',[.02,.01]);return text
        text=block(.123)+fmt('DONOR_RESPONSE_TOTAL_ENERGY',[-300])
        for i in range(n):text+=fmt('DONOR_RESPONSE_TOTAL_GRADIENT',[i+1,0,0,0])
        for i in range(3):text+=fmt('DONOR_RESPONSE_TOTAL_VIRIAL',[i+1,0,0,0])
        text+=block(.456)*48
        p=out/'synthetic_FD_readback.out';p.write_text(text);rb=q.readback(p,job,row,h)
        check('49 solve blocks, central u/v precede48 displaced solves',rb['all_solve_count']==49 and rb['u_eA'][0][0]==.123 and rb['v_eA'][0][0]==.133)
        check('actual token rounding retained',rb['u_print_bounds_eA'][0][0]>0 and rb['gradient_print_bounds'][0][0]>0)
        pbad=out/'synthetic_missing_solve.out';pbad.write_text(text[:-len(block(.456))]);rejects('missing displaced solve rejects',lambda:q.readback(pbad,job,row,h))
        pbad=out/'synthetic_bad_kappa.out';pbad.write_text(text.replace(fmt('DONOR_RESPONSE_PAIR',ex['pairs'][0]),fmt('DONOR_RESPONSE_PAIR',[*ex['pairs'][0][:3],ex['pairs'][0][3]+.01,*ex['pairs'][0][4:]]),1));rejects('changed fixed kappa rejects',lambda:q.readback(pbad,job,row,h))
        for metal in ('La','Dy'):
            j=next(j for j in m['programs'] if j['cell']==metal+'_central' and j['profile']=='combined' and j['stage']=='C' and j['program']=='analyze');r=q.metadata(m['rows'][metal+'_central'],j)
            a=q.regularized(r,'drive');c=q.regularized(r,'stiffness');b=q.regularized(r,'combined')
            check(metal+' independent drive/stiffness support',np.count_nonzero(a['beta'])==0 and np.count_nonzero(c['h'])==0 and np.array_equal(a['h'],b['h']) and np.array_equal(c['beta'],b['beta']) and np.count_nonzero(b['beta'])==2)
        Ha=m['units']['hartree_kcal_mol'];muq=-.73;step=.0005;coord=.005
        Fq=lambda f:4+Ha*f*muq
        mu=lambda qq:.2+muq*qq
        check('Maxwell positive sign /Hartree units',abs((Fq(step)-Fq(-step))/(2*step*Ha)-(mu(coord)-mu(-coord))/(2*coord))<1e-12)
        # Exercise the real orchestration and detached publication under mocks.
        originals={k:getattr(q,k) for k in ('context','parse_cell','baseline_checks','operator_checks','field_checks','mechanics_checks')}
        fake=out/'fake_executable';fake.write_text('not executable; mocked invocation only\n');fp=q.pin(fake)
        mockm=copy.deepcopy(m);mockm['sources']['original_testgrad']=fp;calls=[]
        cfg={'workers':12,'manifest':q.pin(manifest_path),'analyze':fp,'testgrad':fp}
        def invoke(argv,cwd,output,receipts):
            calls.append(str(cwd));Path(output).write_text('mock only\n');receipts.append({'argv':argv,'returncode':0,'output':q.pin(output)})
        q.context=lambda cfg,runtime:(mockm,{'native':SimpleNamespace(invoke=invoke)})
        q.parse_cell=lambda *args:{'numerical_checks_passed':True,'value':{'energy':0.}}
        for k in ('baseline_checks','operator_checks','field_checks','mechanics_checks'):setattr(q,k,lambda *args:{'passed':True})
        oldenv=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='source_mock_no_allocation'
        try:
            full=q.run(cfg,out/'mock_full')
            check('mock121 unique starts and576 completed FD',len(calls)==len(set(calls))==121 and full['native_program_invocations']==121 and full['internal_FD_energies_completed']==576 and full['numerical_checks_passed'])
            check('mock all published complete, no shared mutable rows',all(r['status']=='complete' for r in full['rows'].values()))
            calls.clear();cfg['workers']=1
            def fail(argv,cwd,output,receipts):
                calls.append(str(cwd));Path(output).write_text('synthetic native failure\n');receipts.append({'argv':argv,'returncode':1,'output':q.pin(output)});raise ValueError('synthetic first native failure')
            q.context=lambda cfg,runtime:(mockm,{'native':SimpleNamespace(invoke=fail)})
            failed=q.run(cfg,out/'mock_native_failure')
            check('first failure leaves120 untouched',len(calls)==1 and failed['native_program_invocations']==1 and sum(r['status']=='not_attempted' for r in failed['rows'].values())==120 and failed['global_stop'])
            broken=copy.deepcopy(mockm);broken['programs'][0]['input']['xyz']['sha256']='0'*64
            q.context=lambda cfg,runtime:(broken,{'native':SimpleNamespace(invoke=invoke)})
            calls.clear();pre=q.run(cfg,out/'mock_precall_failure')
            check('precall failure reserved1 actual0',not calls and pre['requested_slots']==1 and pre['native_program_invocations']==0 and sum(r['status']=='not_attempted' for r in pre['rows'].values())==120)
        finally:
            for k,v in originals.items():setattr(q,k,v)
            if oldenv is None:os.environ.pop('SLURM_JOB_ID',None)
            else:os.environ['SLURM_JOB_ID']=oldenv
        result['status']='passed_donor_source_controls'
    except BaseException as exc:result['status']='failed_donor_source_controls';result['error']=repr(exc)
    finally:
        result['checks_passed']=sum(c['passed'] for c in checks);result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1};q.save(out/'RESULT.json',result)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();r=run(a.manifest,a.out);print(json.dumps({'status':r['status'],'checks':r['checks_passed'],'error':r.get('error')}));raise SystemExit(0 if r['status'].startswith('passed') else 1)
