#!/usr/bin/env python3
"""Reuse tighter3; run15 remaining canonical controls; unchanged full35 score."""
import argparse,copy,importlib.util,json,math,os,socket,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('qualified_tighter_recovery',HERE.parent/'tighter_dy_acetate_v1/one_node_recovery_v1/RECOVER.py')
recovery=importlib.util.module_from_spec(sp);sp.loader.exec_module(recovery)
tight=recovery.old;old=tight.old;b=tight.b;read=tight.read;deck=tight.deck;parse=tight.parse
REUSED=tuple('control__cc__acetate_bi__Dy__'+role for role in tight.ROLES)
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def sources(cfg):
    b.check(cfg['implementation']==b.pin(__file__) and cfg['declaration']==b.pin(HERE/'DECLARATION.md') and cfg['run_wrapper']==b.pin(HERE/'RUN.sbatch'),'Source/declaration changed')
    for key in ('tighter_adapter','recovery_adapter','parent_panel','scout_collection','scout_actual_review','failed_tight_controls'):b.verify(cfg[key])
    b.check(cfg['tighter_adapter']==b.pin(tight.__file__) and cfg['recovery_adapter']==b.pin(recovery.__file__),'Frozen helper changed')
    parent=tight.rec.load(b.verify(cfg['parent_panel']));scout=read(cfg['scout_collection']);review=read(cfg['scout_actual_review'])
    b.check(scout['status']=='dy_acetate_tighter_attribution_complete' and scout['complete']==scout['declared']==3 and scout['all_scout_diagnostics_passed'],'Complete passing tighter3 required')
    b.check(review['status']=='passed_tighter_dy_acetate_recovery_actual_review' and review['collection']==cfg['scout_collection'] and review['implementation']==cfg['recovery_adapter'] and review['diagnostic_limit_passed'],'Actual tighter3 review differs')
    failed=read(cfg['failed_tight_controls']);b.check(failed['status']=='local_controls_failed' and failed['complete']==18 and not failed['all_controls_passed'] and failed['panel']==cfg['parent_panel'],'Original failedTight18 must remain preserved')
    m,sp,sc=recovery.load(b.verify(scout['manifest']));recovery.release_contract(m,scout['manifest']['path'],read(scout['release']))
    b.check(sp['parent_panel']==cfg['parent_panel'] and sp['failed_controls']==cfg['failed_tight_controls'],'Scout canonical/development lineage differs')
    b.check(parent['assets']==sp['assets'] and parent['orca']==sp['orca'] and old.TIERS['tighter']==sp['tier'],'Same native recipe/assets required')
    return parent,scout,m,sp,sc

def definitions(parent):
    cells=[]
    for c in parent['cells']:
        if c['stage']!='controls':continue
        c=copy.deepcopy(c);c.pop('template_input');c['tier']='tighter';cells.append(c)
    b.check(len(cells)==18 and {c['cell_id'] for c in cells}>=set(REUSED),'Original18 controls required');return cells

def reused_rows(cfg,p=None):
    parent,scout,m,sp,sc=sources(cfg);rows={}
    for cell in sp['cells']:
        cid=cell['cell_id'];expected='control__'+cid.split('tighter__',1)[1];stored=scout['rows'][cid]
        if cid==recovery.CID:q=recovery.admit_replacement(stored['result']['path'],m,sp,sc,scout['manifest']['path'],scout['release']['path'])
        else:q=tight.admit_result(stored['result']['path'],cell,sp,scout['panel']['path'])
        q=json.loads(json.dumps(q));b.check(q==stored,'Exact admitted tighter3 raw replay differs')
        if p is not None:
            same=next(c for c in p['cells'] if c['cell_id']==expected);b.check(deck(p,same,8,3000)==deck(sp,cell,8,3000),'Reused recipe/geometry differs')
        rows[expected]=dict(q,reused_from_tighter_scout=True,original_scout_cell_id=cid)
    b.check(set(rows)==set(REUSED),'Three reused roles differ');return rows

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text());parent,scout,m,sp,sc=sources(cfg);cells=definitions(parent);out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()
    p=dict(adapter=b.pin(__file__),configuration=b.pin(config),declaration=cfg['declaration'],run_wrapper=cfg['run_wrapper'],tighter_adapter=cfg['tighter_adapter'],recovery_adapter=cfg['recovery_adapter'],parser=b.pin(tight.rec.__file__),native_input_adapter=parent['adapter'],parent_panel=cfg['parent_panel'],scout_collection=cfg['scout_collection'],scout_actual_review=cfg['scout_actual_review'],failed_tight_controls=cfg['failed_tight_controls'],assets=parent['assets'],orca=parent['orca'],controls=parent['controls'],cells=cells,tier=old.TIERS['tighter'],reused_cell_ids=list(REUSED),new_cell_ids=[c['cell_id'] for c in cells if c['cell_id'] not in REUSED],maximum_new_native_calls=15,declared_controls=18,qualification_check_count=35,ranks=8,nodes=1,execution_root=cfg['execution_root'],created_utc=b.now())
    deltas={}
    for c in cells:
        previous=next(x for x in parent['cells'] if x['cell_id']==c['cell_id']);raw=b.verify(previous['template_input']).read_text();expected=raw.replace('%pal nprocs 4 end','%pal nprocs 8 end').replace(' TCutPairs 1e-05\n',' TCutPairs 1e-06\n').replace(' TCutPNO 1e-07\n',' TCutPNO 1e-08\n').replace(' TCutMKN 0.001\n',' TCutMKN 0.0001\n');actual=deck(p,c,8,3000);b.check(actual==expected,'Unexpected deck delta')
        dst=out/'inputs'/(c['cell_id']+'.inp');dst.write_text(actual);c['template_input']=b.pin(dst);deltas[c['cell_id']]={'original_input':previous['template_input'],'tighter_input':c['template_input'],'only_three_cutoffs_and_MPI_changed':True,'new_native_required':c['cell_id'] not in REUSED}
    p['reused']=reused_rows(cfg,p)
    # Reuse the original exact canonical parser; only these18 archived control roles.
    canonical={k:old.saved_canonical(v['result']) for k,v in p['controls'].items()};b.check(canonical==p['controls'] and len(canonical)==18,'Original canonical18 raw evidence differs')
    panel=b.write(out/'PANEL.json',p);load(panel['path']);q=dict(status='passed_tighter_controls_preparation',implementation=b.pin(__file__),panel=panel,configuration=b.pin(config),new_native_calls=0,new_declared=15,reused=3,full_declared=18,original_qualification_checks=35,input_deltas=deltas,canonical_raw_replays=18,tighter_raw_replays=3,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1);rp=b.write(out/'RESULT.json',q);print(json.dumps({'result':rp,'panel':panel,'cost':{k:q[k] for k in ('wall_seconds','CPU_seconds')}}))

def load(path):
    p=json.loads(Path(path).read_text());cfg=read(p['configuration']);parent,scout,m,sp,sc=sources(cfg)
    b.check(p['adapter']==b.pin(__file__) and p['parser']==b.pin(tight.rec.__file__) and p['native_input_adapter']==parent['adapter'],'Prepared source differs')
    for k in ('declaration','run_wrapper','parent_panel','tighter_adapter','recovery_adapter','scout_collection','scout_actual_review','failed_tight_controls','execution_root'):b.check(p[k]==cfg[k],'Prepared configuration differs: '+k)
    b.check(p['assets']==parent['assets'] and p['orca']==parent['orca'] and p['controls']==parent['controls'] and p['tier']==old.TIERS['tighter'],'Original chemistry/canonical evidence differs')
    b.check(p['reused_cell_ids']==list(REUSED) and p['new_cell_ids']==[c['cell_id'] for c in p['cells'] if c['cell_id'] not in REUSED] and len(p['new_cell_ids'])==p['maximum_new_native_calls']==15 and p['declared_controls']==18 and p['qualification_check_count']==35 and p['ranks']==8 and p['nodes']==1,'Finite15+3 scope differs')
    for c,e in zip(p['cells'],definitions(parent)):
        b.check({k:v for k,v in c.items() if k!='template_input'}==e,'Exact control role differs');b.check(b.verify(c['template_input']).read_text()==deck(p,c,8,3000),'Prepared deck differs')
    b.check(len(p['cells'])==18 and set(p['reused'])==set(REUSED),'Full18/reuse denominator differs')
    for cid,q in p['reused'].items():
        oldid=q['original_scout_cell_id'];stored=scout['rows'][oldid];b.check({k:v for k,v in q.items() if k not in ('reused_from_tighter_scout','original_scout_cell_id')}==stored and q['reused_from_tighter_scout'],'Stored reuse authority differs');b.verify(q['result'])
        cell=next(c for c in p['cells'] if c['cell_id']==cid);oldcell=next(c for c in sp['cells'] if c['cell_id']==oldid);b.check(deck(p,cell,8,3000)==deck(sp,oldcell,8,3000),'Reused deck byte identity differs')
    return p

def release_gate(p,path,r):
    b.check(r['status']=='released_tighter_controls15' and r['implementation']==p['adapter'] and r['panel']==b.pin(path) and r['parser']==p['parser'] and r['declaration']==p['declaration'],'Root source release differs')
    b.check(r['cell_ids']==p['new_cell_ids'] and r['maximum_new_native_calls']==15 and r['ranks']==8 and r['nodes']==1 and r['execution_root']==p['execution_root'] and r['cluster_release_authorized'] is False,'Release finite scope differs')
    b.check(r['scout_actual_review']==p['scout_actual_review'] and r['scout_collection']==p['scout_collection'],'Actual three-role admission must bind before release')
    prep=read(r['preparation']);review=read(r['source_review']);b.check(prep['status']=='passed_tighter_controls_preparation' and prep['implementation']==p['adapter'] and prep['panel']==b.pin(path) and prep['new_native_calls']==0,'Preparation missing')
    b.check(review['status']=='passed_tighter_controls_source_review' and review['implementation']==p['adapter'] and review['panel']==b.pin(path) and review['preparation']==r['preparation'],'Independent source review differs');return r

def resources(env):
    ranks,mem,maxcore=tight.resources(env);nodes=int(env.get('SLURM_JOB_NUM_NODES',env.get('SLURM_NNODES','0')));b.check(nodes==1,'New15 require one node');return ranks,mem,maxcore,nodes

def run(path,release,cell_id,out):
    p=load(path);r=release_gate(p,path,json.loads(Path(release).read_text()));b.check(cell_id in p['new_cell_ids'],'Only new15 may start');c=next(c for c in p['cells'] if c['cell_id']==cell_id);b.check(out.resolve()==Path(p['execution_root']).resolve()/'tasks'/cell_id/'execution' and not out.exists(),'One-attempt path differs or exists')
    ranks,mem,maxcore,nodes=resources(os.environ);out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for k in THREADS:env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];provenance={k:p[k] for k in ('adapter','parser','native_input_adapter','tighter_adapter')}
    st=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=ranks,nodes=nodes,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now(),**provenance));tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick;q.update(cell_id=cell_id,cell=c,panel=b.pin(path),release=b.pin(release),start=st,returncode=code,native_program_invocations_attempted=1,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),input=b.pin(out/'endpoint.inp'),output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),**provenance);b.write(out/'RESULT.json',q);print(json.dumps({'status':q['status'],'cell_id':cell_id,'elapsed_seconds':elapsed}));return 0 if q['status']=='complete' else 1

def admit_result(path,c,p,panel):
    path=Path(path);r=json.loads(path.read_text());b.check(c['cell_id'] in p['new_cell_ids'] and path.resolve()==Path(p['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','Native result role/path differs')
    b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()) and r['cell']==c and r['cell_id']==c['cell_id'] and r['panel']==b.pin(panel),'Complete native cell required')
    b.check(all(r[k]==p[k] for k in ('adapter','parser','native_input_adapter','tighter_adapter')),'Native source lineage differs');release_gate(p,panel,read(r['release']));st=read(r['start'])
    b.check(r['start']==b.pin(path.parent/'START.json') and all(st[k]==r[k] for k in ('cell','panel','release','input','job_id','adapter','parser','native_input_adapter','tighter_adapter')),'Native START differs')
    b.check(st['ranks']==8 and st['nodes']==1 and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/8) and st['maximum_native_calls']==r['native_program_invocations_attempted']==1 and st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native resource/call scope differs')
    b.check(b.verify(r['input']).read_text()==deck(p,c,8,st['maxcore_MiB']),'Executed scientific deck differs');b.verify(r['stderr']);q=json.loads(json.dumps(parse(b.verify(r['output']).read_text(),c,r['returncode'])))
    b.check(q['status']=='complete' and all(r[k]==v for k,v in q.items()),'Full raw native replay differs');return dict(q,result=b.pin(path),input=r['input'],output=r['output'],start=r['start'],release=r['release'],native_producer=r['adapter'],elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'],reused_from_tighter_scout=False)

def score(p,rows):
    b.check(set(rows)=={c['cell_id'] for c in p['cells']} and len(rows)==18 and all(r['status']=='complete' for r in rows.values()),'All18 admissions required; no partial qualification');q=old.control_score(p,rows);b.check(len(q['checks'])==35,'Original35 qualification denominator differs');return q

def collect(path,inventory,out):
    tick,cpu=time.monotonic(),time.process_time();p=load(path);inv=json.loads(Path(inventory).read_text());b.check(set(inv)==set(p['new_cell_ids']),'Exact15 new result slots required');out.mkdir(parents=True,exist_ok=False);rows={c['cell_id']:dict(status='unavailable',error=None) for c in p['cells']};good={}
    try:reuse=reused_rows(read(p['configuration']),p)
    except Exception as e:reuse={};reuse_error=repr(e)
    for c in p['cells']:
        cid=c['cell_id']
        try:
            if cid in REUSED:b.check(cid in reuse,'Reused tighter3 unavailable: '+str(locals().get('reuse_error')));q=reuse[cid]
            else:b.check(inv[cid] is not None,'Missing new native result');q=admit_result(inv[cid],c,p,path)
            rows[cid]=dict(q);good[cid]=q
        except Exception as e:rows[cid]=dict(status='unavailable',path=inv.get(cid),error=repr(e))
    q=dict(status='unavailable',implementation=p['adapter'],parser=p['parser'],panel=b.pin(path),inventory=b.pin(inventory),declared=18,new_declared=15,reuse_declared=3,complete=len(good),new_complete=sum(cid in good for cid in p['new_cell_ids']),reused_complete=sum(cid in good for cid in REUSED),rows=rows,all_controls_passed=False,cluster_release_authorized=False,physical_accuracy_qualified=False,new_native_calls=0)
    if len(good)==18:q.update(score(p,good));q['status']='tighter_local_controls_passed' if q['all_controls_passed'] else 'tighter_local_controls_failed'
    scout=read(p['scout_collection']);q.update(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,new_native_allocated_CPU_seconds=sum(v['allocated_cpu_seconds'] for k,v in good.items() if k not in REUSED),reused_native_allocated_CPU_seconds=sum(v['allocated_cpu_seconds'] for k,v in good.items() if k in REUSED),historical_abandoned_scheduler_CPU_seconds=scout['abandoned_attempt_scheduler_CPU_seconds']);b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','declared','complete','new_complete','reused_complete','new_native_calls')}));return q

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True);a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',required=True,type=Path);a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',required=True,type=Path);a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--out',required=True,type=Path);args=ap.parse_args()
    if args.mode=='prepare':prepare(args.config,args.out)
    elif args.mode=='run':sys.exit(run(args.panel,args.release,args.cell,args.out))
    else:collect(args.panel,args.inventory,args.out)
