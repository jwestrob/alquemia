#!/usr/bin/env python3
"""One immutable MaxDIIS12 recovery; original canonical chemistry and gates."""
import argparse, difflib, importlib.util, json, math, os, re, socket, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('canonical_rotation_v2',HERE.parent/'dry_acetate_correlated_rotation_v2/REFERENCE.py')
old=importlib.util.module_from_spec(sp);sp.loader.exec_module(old)
b=old.b;read=old.read
CID='cc__Dy__rotation__plus_005__complex'
WORK=b.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'/HERE.name
OUTPUT=WORK/'run_v1/execution'
EXPECTED=dict(MaxIter=100,STol=1e-8,LShift=.2,MaxDIIS=12)
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
F=old.old.donor.F

def write(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
    return b.pin(path)

def scope(cell_id,out):
    b.check(cell_id==CID,'Only frozen failed Dy +005 cell is authorized')
    b.check(Path(out).resolve()==OUTPUT,'Only new fixed recovery output is authorized')

def deck(p,c,ranks,maxcore):
    b.check(c['cell_id']==CID and ranks==4,'Recovery cell/ranks differ')
    raw=old.deck(p,c,ranks,maxcore)
    b.check(raw.count('%mdci\n')==1 and 'MaxDIIS' not in raw,'Unexpected original MDCI block')
    return raw.replace('%mdci\n','%mdci\n MaxDIIS 12\n',1)

def solver_readback(text):
    native=text.split('****END OF INPUT****')[-1]
    names={'MaxIter':'Maximum number of iterations','STol':r'Convergence tolerance \(max\. residuum\)',
           'LShift':'Level shift for amplitude update','MaxDIIS':'Maximum number of DIIS vectors'}
    result={}
    for name,pattern in names.items():
        values=re.findall(r'^\s*'+pattern+r'\s*\.{2,}\s*('+F+r')\s*$',native,re.M)
        b.check(len(values)==1,'Missing/ambiguous native solver setting '+name)
        value=old.old.donor.cp.number(values[0]);b.check(math.isfinite(value),'Nonfinite solver setting')
        result[name]=value
    return result

def parse(text,c,returncode):
    # Preserve the established raw parser, including failed SCF criteria as flags.
    try:result=old.parse(text,c,returncode)
    except Exception as exc:
        result=old.old.donor.parse(text,c,returncode)
        result.update(status='failed_admission',error=repr(exc))
    try:settings=solver_readback(text);valid=settings==EXPECTED
    except Exception as exc:settings=None;valid=False;result['solver_readback_error']=repr(exc)
    result['solver_readback']=settings;result['checks']['recovery_solver_settings']=valid
    if result.get('error') or not all(result['checks'].values()):result['status']='failed_admission'
    return result

def context():
    inputs=read(HERE/'INPUTS.json');old.old.verify_tree(inputs)
    b.check(inputs['original_implementation']==b.pin(old.__file__),'Original v2 source changed')
    p=old.load(b.verify(inputs['original_panel']));c=next(x for x in p['cells'] if x['cell_id']==CID)
    failed=read(b.verify(inputs['failed_result']));start=read(b.verify(inputs['failed_start']))
    b.check(failed['status']=='failed_admission' and failed['job_id']=='18866' and failed['cell']==c and failed['panel']==inputs['original_panel'],'Original failed cell differs')
    b.check(failed['adapter']==inputs['original_implementation'] and failed['new_native_calls']==1,'Original producer/call differs')
    for key in ('input','output','stderr','start'):b.check(failed[key]==inputs['failed_'+key],'Original raw '+key+' differs')
    b.check(failed['release']==inputs['original_release'],'Original release differs')
    release=old.release_gate(p,inputs['original_release'],CID)
    b.check(release['task_output_directories'][CID]==str(Path(inputs['failed_result']['path']).parent),'Original native path differs')
    for key in ('cell','panel','adapter','release','input','job_id'):b.check(start[key]==failed[key],'Original START '+key+' differs')
    b.check(start['ranks']==4 and start['nodes']==1 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Original resources differ')
    b.check(start['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Original executable/command differs')
    raw=b.verify(inputs['failed_input']).read_text();b.check(raw==old.deck(p,c,4,start['maxcore_MiB']),'Original scientific input differs')
    text=b.verify(inputs['failed_output']).read_text();settings=solver_readback(text)
    b.check(settings==dict(EXPECTED,MaxDIIS=7),'Original native solver settings differ')
    parsed=old.old.donor.parse(text,c,failed['returncode'])
    b.check(len(parsed['CC_residual_trace'])==100 and parsed['CC_residual_trace'][-1][3]==1.38e-6,'Original 100-step failure trace differs')
    b.check(not parsed['checks']['cc_converged'] and not parsed['checks']['normal_termination'] and not parsed['checks']['declared_CC_residual'],'Original failure no longer explicit')
    try:old.parse(text,c,failed['returncode'])
    except ValueError as exc:b.check(repr(exc)==failed['error'],'Original parser failure differs')
    else:raise ValueError('Original failed output unexpectedly admitted')
    return inputs,p,c,failed,settings

def load(path):
    inputs,p,c,failed,settings=context();m=read(path)
    b.check(m['status']=='prepared_one_cell_CC_DIIS_recovery' and m['implementation']==b.pin(__file__) and m['inputs']==b.pin(HERE/'INPUTS.json') and m['declaration']==b.pin(HERE/'DECLARATION.md'),'Recovery source/manifest differs')
    b.check(m['original_panel']==inputs['original_panel'] and m['failed_result']==inputs['failed_result'] and m['cell']==c and m['maximum_native_calls']==1 and m['solver_settings']==EXPECTED,'Recovery cell/scope differs')
    scope(m['cell_id'],m['output_directory'])
    b.check(b.verify(m['template_input']).read_text()==deck(p,c,4,3000),'Prepared recovery input differs')
    b.check(b.verify(m['deck_diff']).read_text()==''.join(difflib.unified_diff(old.deck(p,c,4,3000).splitlines(True),deck(p,c,4,3000).splitlines(True),fromfile='original.inp',tofile='MaxDIIS12.inp')),'Prepared deck delta differs')
    return m,p,c

def prepare(out):
    tick,cpu=time.monotonic(),time.process_time()
    b.check(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')
    b.check(out.resolve().parent==WORK and out.name.startswith('preparation_'),'Preparation must use owned versioned path')
    inputs,p,c,failed,settings=context();before=old.deck(p,c,4,3000);after=deck(p,c,4,3000)
    b.check(after.replace(' MaxDIIS 12\n','',1)==before and after.count('MaxDIIS')==1,'Only MaxDIIS insertion allowed')
    checks=['original_failed_raw_pins','unchanged_original_panel_and_cell','original_failure_readmission','one_line_deck_delta','native_call_cap_1']
    for cid,path in [('cc__Dy__rotation__plus_010__complex',OUTPUT),(CID,Path(inputs['failed_result']['path']).parent)]:
        try:scope(cid,path)
        except ValueError:checks.append('rejected_wrong_cell' if cid!=CID else 'rejected_old_output_path')
        else:raise ValueError('Unauthorized scope accepted')
    text=b.verify(inputs['failed_output']).read_text();synthetic=re.sub(r'(Maximum number of DIIS vectors\s*\.{2,}\s*)7\b',r'\g<1>12',text)
    b.check(solver_readback(synthetic)==EXPECTED,'Synthetic DIIS12 readback failed');checks.append('synthetic_DIIS12_readback')
    changed=re.sub(r'(Level shift for amplitude update\s*\.{2,}\s*)2.000e-01',r'\g<1>3.000e-01',synthetic)
    b.check(solver_readback(changed)!=EXPECTED,'Changed level shift accepted');checks.append('changed_level_shift_rejected')
    missing=re.sub(r'^.*Maximum number of DIIS vectors.*\n','',synthetic,flags=re.M)
    try:solver_readback(missing)
    except ValueError:checks.append('missing_solver_setting_rejected')
    else:raise ValueError('Missing solver readback accepted')
    b.check(parse(text,c,failed['returncode'])['status']=='failed_admission','Old failure accepted as recovery');checks.append('old_failure_rejected')
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(after)
    (out/'DECK.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original.inp',tofile='MaxDIIS12.inp')))
    manifest=dict(status='prepared_one_cell_CC_DIIS_recovery',implementation=b.pin(__file__),inputs=b.pin(HERE/'INPUTS.json'),declaration=b.pin(HERE/'DECLARATION.md'),original_panel=inputs['original_panel'],failed_result=inputs['failed_result'],cell_id=CID,cell=c,maximum_native_calls=1,solver_settings=EXPECTED,original_solver_settings=settings,ranks=4,nodes=1,output_directory=str(OUTPUT),template_input=b.pin(out/'endpoint.inp'),deck_diff=b.pin(out/'DECK.diff'),original_native_cost_proxy_CPU_seconds=failed['allocated_cpu_seconds'])
    mp=write(out/'MANIFEST.json',manifest);load(mp['path'])
    q=dict(status='passed_one_cell_CC_DIIS_preflight',implementation=b.pin(__file__),manifest=mp,checks=checks,synthetic_readback_is_not_native_evidence=True,new_native_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    pp=write(out/'PREFLIGHT.json',q);print(json.dumps(dict(manifest=mp,preflight=pp,status=q['status'])))

def release_gate(manifest_path,release_path,cell_id,out):
    scope(cell_id,out);m,p,c=load(manifest_path);r=read(release_path)
    b.check(r['status']=='released_one_cell_CC_DIIS_recovery' and r['implementation']==b.pin(__file__) and r['manifest']==b.pin(manifest_path) and r['inputs']==b.pin(HERE/'INPUTS.json'),'Root provenance release differs')
    b.check(r['cell_ids']==[CID] and r['maximum_native_calls']==1 and r['ranks']==4 and r['nodes']==1 and r['output_directory']==str(OUTPUT),'Root release scope differs')
    q=read(b.verify(r['preflight']));b.check(q['status']=='passed_one_cell_CC_DIIS_preflight' and q['manifest']==r['manifest'] and q['implementation']==r['implementation'] and q['new_native_calls']==0,'Released preflight differs')
    return m,p,c,r

def run(manifest_path,cell_id,release_path,out):
    m,p,c,release=release_gate(manifest_path,release_path,cell_id,out)
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));nodes=int(os.environ.get('SLURM_JOB_NUM_NODES',os.environ.get('SLURM_NNODES','0')))
    b.check(ranks==4 and nodes==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks on one node required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant missing');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for key in (*THREADS,'BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];rp=b.pin(release_path);mp=b.pin(manifest_path)
    start=dict(adapter=b.pin(__file__),manifest=mp,panel=m['original_panel'],failed_result=m['failed_result'],cell=c,release=rp,input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=ranks,nodes=nodes,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,utc=b.now())
    sp=write(out/'START.json',start);tick=time.monotonic();returncode=None;calls=0
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            proc=subprocess.Popen(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr);calls=1;returncode=proc.wait()
        result=parse((out/'endpoint.out').read_text(errors='replace'),c,returncode)
    except Exception as exc:result=dict(status='failed_admission',error=repr(exc))
    elapsed=time.monotonic()-tick
    result.update(cell_id=cell_id,cell=c,panel=m['original_panel'],manifest=mp,failed_result=m['failed_result'],adapter=b.pin(__file__),release=rp,start=sp,input=start['input'],output=b.pin(out/'endpoint.out') if (out/'endpoint.out').exists() else None,stderr=b.pin(out/'endpoint.err') if (out/'endpoint.err').exists() else None,returncode=returncode,job_id=start['job_id'],elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,new_native_calls=calls,finished_utc=b.now())
    write(out/'RESULT.json',result);print(json.dumps(dict(status=result['status'],cell_id=cell_id,seconds=elapsed)));return 0 if result['status']=='complete' else 1

def admit(manifest_path,result_path):
    result_path=Path(result_path);r=read(result_path);out=result_path.resolve().parent
    m,p,c,release=release_gate(manifest_path,b.verify(r['release']),r['cell_id'],out)
    b.check(r['status']=='complete' and r['adapter']==b.pin(__file__) and r['manifest']==b.pin(manifest_path) and r['cell']==c and r['panel']==m['original_panel'] and r['failed_result']==m['failed_result'] and r['new_native_calls']==1,'Recovered result identity differs')
    for key,name in [('input','endpoint.inp'),('output','endpoint.out'),('stderr','endpoint.err'),('start','START.json')]:b.check(r[key]==b.pin(out/name),'Native artifact differs '+key)
    s=read(b.verify(r['start']))
    for key in ('adapter','manifest','panel','failed_result','cell','release','input','job_id'):b.check(s[key]==r[key],'START differs '+key)
    b.check(s['ranks']==4 and s['nodes']==1 and s['maxcore_MiB']==math.floor(.75*s['memory_MiB']/4) and s['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native resources differ')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,s['maxcore_MiB']),'Executed recovery input differs')
    replay=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(replay['status']=='complete','Raw recovery admission failed');old.raw_equal(r,replay)
    b.check(math.isfinite(r['elapsed_seconds']) and r['elapsed_seconds']>=0 and r['allocated_cpu_seconds']==4*r['elapsed_seconds'],'Native cost differs');return r

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--out',type=Path,required=True)
    q=sub.add_parser('run');q.add_argument('--manifest',required=True);q.add_argument('--cell',required=True);q.add_argument('--release',required=True);q.add_argument('--out',type=Path,required=True)
    q=sub.add_parser('admit');q.add_argument('--manifest',required=True);q.add_argument('--result',required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    elif a.action=='run':return run(a.manifest,a.cell,a.release,a.out)
    else:print(json.dumps(dict(status=admit(a.manifest,a.result)['status'],new_native_calls=0)))
if __name__=='__main__':sys.exit(main())
