#!/usr/bin/env python3
"""One-node executor and saved three-role collection; unchanged native helpers."""
import argparse, importlib.util, json, os, socket, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_tighter_dy_adapter',HERE.parent/'REFERENCE.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
b=old.b
CID='tighter__cc__acetate_bi__Dy__ligand_ghost_ion'
OLD_SHA='b925ec776e87e750956b47af43ca5d03c2098795a192ee6007f1af2ff8fa8ce9'

def load(manifest):
    m=json.loads(Path(manifest).read_text())
    b.check(m['implementation']==b.pin(__file__) and m['original_adapter']==b.pin(old.__file__)
        and m['original_adapter']['sha256']==OLD_SHA,'Executor/source changed')
    for key in ('configuration','declaration','run_wrapper','collect_wrapper','original_panel','progress_review','startup_comparison','root_decision','original_start','original_input','input'):
        b.verify(m[key])
    cfg=old.read(m['configuration'])
    b.check(all(m[k]==v for k,v in cfg.items()),'Preparation configuration differs')
    p=old.load(b.verify(m['original_panel']))
    c=next(c for c in p['cells'] if c['cell_id']==CID)
    st=old.read(m['original_start']);decision=old.read(m['root_decision'])
    b.check(st['job_id']=='18518' and st['cell']==c and st['panel']==m['original_panel'] and st['input']==m['original_input'],'Original endpoint lineage differs')
    b.check(decision['status']=='declared_owned_startup_stall_recovery' and decision['job_id']==18518
        and decision['scientific_settings_changed'] is False,'Root decision differs')
    b.check(m['cell']==c and m['maximum_new_native_calls']==1 and m['nodes']==1 and m['ranks']==8,'Recovery scope differs')
    b.check(b.verify(m['input']).read_bytes()==b.verify(m['original_input']).read_bytes()
        and b.verify(m['input']).read_text()==old.deck(p,c,8,3000),'Input bytes changed')
    b.check(set(m['original_results'])=={x['cell_id'] for x in p['cells'] if x['cell_id']!=CID},'Two retained roles differ')
    for authority in m['original_results'].values():b.verify(authority)
    return m,p,c

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__) and cfg['original_adapter']==b.pin(old.__file__),'Source differs')
    p=old.load(b.verify(cfg['original_panel']));c=next(c for c in p['cells'] if c['cell_id']==CID)
    out.mkdir(parents=True,exist_ok=False)
    target=out/'endpoint.inp';target.write_bytes(b.verify(cfg['original_input']).read_bytes())
    b.check(target.read_text()==old.deck(p,c,8,3000),'Original byte identity differs')
    m=dict(cfg,configuration=b.pin(config),cell=c,input=b.pin(target),maximum_new_native_calls=1,nodes=1,ranks=8)
    mp=b.write(out/'MANIFEST.json',m);load(mp['path'])
    for cid,authority in m['original_results'].items():
        r=old.read(authority)
        b.check(r['status']=='complete' and r['cell_id']==cid and r['panel']==m['original_panel'] and all(r['checks'].values()),'Retained result metadata differs')
    q=dict(status='passed_single_node_recovery_preparation',implementation=b.pin(__file__),manifest=mp,
        exact_original_input_bytes=True,retained_completed_roles=2,new_native_calls=0,threads=1,
        wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    result=b.write(out/'RESULT.json',q);print(json.dumps(dict(result=result,manifest=mp)));return q

def release_contract(m,manifest,r):
    b.check(r['status']=='approved_single_node_tighter_attribution' and r['implementation']==m['implementation']
        and r['manifest']==b.pin(manifest) and r['root_decision']==m['root_decision'],'Release authority differs')
    b.check(r['cell_id']==CID and r['maximum_new_native_calls']==1 and r['nodes']==1 and r['ranks']==8
        and r['execution_out']==m['execution_out'],'Release scope differs')
    prep=old.read(r['preparation']);review=old.read(r['source_review'])
    b.check(prep['status']=='passed_single_node_recovery_preparation' and prep['manifest']==r['manifest'] and prep['implementation']==m['implementation'] and prep['new_native_calls']==0,'Preparation missing')
    b.check(review['status']=='passed_single_node_recovery_source_review' and review['implementation']==m['implementation']
        and review['manifest']==r['manifest'] and review['preparation']==r['preparation'],'Reviewed source missing')
    terminal=old.read(r['original_terminal_receipt'])['records']['18518']
    b.check(terminal['state'].startswith('CANCELLED') and str(terminal['allocated_cpus'])=='8','Canceled original terminal receipt missing')
    original=old.read(r['original_collection'])
    b.check(original['status']=='unavailable' and original['declared']==3 and original['complete']==2
        and original['panel']==m['original_panel'] and original['implementation']==m['original_adapter'],'Preserved original incomplete collection missing')
    b.check(set(original['rows'])==set(m['original_results'])|{CID} and original['rows'][CID]['status']=='unavailable','Original missing role differs')
    for cid,pin in m['original_results'].items():
        b.check(original['rows'][cid]['status']=='complete' and original['rows'][cid]['result']==pin,'Original retained-result identity differs')
    b.check(r['full_eighteen_qualified'] is False and r['cluster_release_authorized'] is False,'Attribution cannot release clusters')
    return r

def resources(env):
    ranks,mem,maxcore=old.resources(env)
    nodes=int(env.get('SLURM_JOB_NUM_NODES',env.get('SLURM_NNODES','0')))
    b.check(nodes==1 and maxcore==3000 and mem==32000,'Require one node and unchanged grant-derived MaxCore')
    return ranks,mem,maxcore,nodes

def run(manifest,release):
    m,p,c=load(manifest);r=release_contract(m,manifest,json.loads(Path(release).read_text()))
    ranks,mem,maxcore,nodes=resources(os.environ);out=Path(m['execution_out'])
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_bytes(b.verify(m['input']).read_bytes())
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    provenance={k:p[k] for k in ('adapter','native_input_adapter','parser')}
    provenance.update(runtime_adapter=b.pin(__file__),manifest=b.pin(manifest),supersedes_original_job='18518')
    st=b.write(out/'START.json',dict(cell=c,panel=m['original_panel'],release=b.pin(release),input=b.pin(out/'endpoint.inp'),
        job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=ranks,nodes=nodes,memory_MiB=mem,maxcore_MiB=maxcore,
        command=cmd,maximum_native_calls=1,utc=b.now(),**provenance))
    tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=old.parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick
    q.update(cell_id=CID,cell=c,panel=m['original_panel'],**provenance,release=b.pin(release),start=st,returncode=code,
        native_program_invocations_attempted=1,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),
        input=b.pin(out/'endpoint.inp'),output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'))
    b.write(out/'RESULT.json',q);print(json.dumps(dict(status=q['status'],cell_id=CID,seconds=elapsed)));return 0 if q['status']=='complete' else 1

def admit_replacement(path,m,p,c,manifest,release):
    b.check(Path(path).resolve()==Path(m['execution_out']).resolve()/'RESULT.json','Replacement path differs')
    r=json.loads(Path(path).read_text());st=old.read(r['start'])
    b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Replacement unavailable')
    b.check(r['cell']==c and r['cell_id']==CID and r['panel']==m['original_panel'] and r['manifest']==b.pin(manifest)
        and r['runtime_adapter']==m['implementation'] and r['release']==b.pin(release),'Replacement provenance differs')
    b.check(all(r[k]==p[k] for k in ('adapter','native_input_adapter','parser')) and r['supersedes_original_job']=='18518','Native helper identity differs')
    b.check(r['start']==b.pin(Path(path).parent/'START.json') and all(st[k]==r[k] for k in ('cell','panel','manifest','runtime_adapter','release','input','adapter','native_input_adapter','parser','job_id','supersedes_original_job')),'START join differs')
    b.check(st['ranks']==8 and st['nodes']==1 and st['memory_MiB']==32000 and st['maxcore_MiB']==3000 and st['maximum_native_calls']==r['native_program_invocations_attempted']==1,'Resource/call receipt differs')
    b.check(st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native invocation differs')
    b.check(b.verify(r['input']).read_bytes()==b.verify(m['input']).read_bytes(),'Runtime input bytes differ');b.verify(r['stderr'])
    q=json.loads(json.dumps(old.parse(b.verify(r['output']).read_text(),c,0)))
    b.check(q['status']=='complete' and all(r[k]==v for k,v in q.items()),'Raw replacement replay differs')
    return dict(q,result=b.pin(path),output=r['output'],input=r['input'],start=r['start'],release=r['release'],native_producer=r['adapter'],runtime_adapter=r['runtime_adapter'],
        elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'])

def collect(manifest,release,inventory,out):
    tick,cpu=time.monotonic(),time.process_time();m,p,c=load(manifest)
    r=release_contract(m,manifest,json.loads(Path(release).read_text()));inv=json.loads(Path(inventory).read_text())
    expected={cid:a['path'] for cid,a in m['original_results'].items()};expected[CID]=str(Path(m['execution_out'])/'RESULT.json')
    b.check(inv==expected,'Exact three-role inventory differs');out.mkdir(parents=True,exist_ok=False)
    rows={};good={}
    for cell in p['cells']:
        cid=cell['cell_id']
        try:
            if cid==CID:q=admit_replacement(inv[cid],m,p,c,manifest,release)
            else:q=old.admit_result(inv[cid],cell,p,m['original_panel']['path'])
            rows[cid]=dict(q);good[cid]=q
        except Exception as e:rows[cid]=dict(status='unavailable',path=inv[cid],error=repr(e))
    q=dict(status='unavailable',implementation=b.pin(__file__),original_adapter=m['original_adapter'],manifest=b.pin(manifest),panel=m['original_panel'],
        release=b.pin(release),inventory=b.pin(inventory),original_incomplete_collection=r['original_collection'],
        original_canceled_receipt=r['original_terminal_receipt'],declared=3,complete=len(good),rows=rows,
        full_eighteen_qualified=False,cluster_release_authorized=False,physical_accuracy_qualified=False,new_native_calls=0)
    if len(good)==3:q.update(status='dy_acetate_tighter_attribution_complete',**old.comparison(p,good))
    q.update(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,
        selected_native_allocated_CPU_seconds=sum(x['allocated_cpu_seconds'] for x in good.values()),
        abandoned_attempt_scheduler_CPU_seconds=int(old.read(r['original_terminal_receipt'])['records']['18518']['elapsed_seconds'])*8)
    b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','declared','complete','new_native_calls')}));return q

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',required=True,type=Path)
    for name in ('run','collect'):
        a=sub.add_parser(name);a.add_argument('--manifest',required=True);a.add_argument('--release',required=True)
        if name=='collect':a.add_argument('--inventory',required=True);a.add_argument('--out',required=True,type=Path)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='run':return run(a.manifest,a.release)
    else:collect(a.manifest,a.release,a.inventory,a.out)
if __name__=='__main__':sys.exit(main())
