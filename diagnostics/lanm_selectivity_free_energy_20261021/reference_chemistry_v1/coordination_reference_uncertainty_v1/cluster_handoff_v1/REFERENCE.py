#!/usr/bin/env python3
"""Eight existing cluster decks; tighter W9 Dy first; root releases 1 then7."""
import argparse,copy,importlib.util,json,math,os,socket,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('admitted_tighter_controls',HERE.parent/'tighter_controls_v1/REFERENCE.py')
controls=importlib.util.module_from_spec(spec);spec.loader.exec_module(controls)
old=controls.old;rec=controls.tight.rec;b=old.b;read=old.read;deck=old.deck;parse=rec.parse
SCOUT='cluster__tighter__aquo9_central_Dy';RANKS=24;NODES=1
THREADS=controls.THREADS
AUTH=('parent_panel','qualified_controls','control_actual_review','control_root_review','failed_tight_controls','cost_accounting')

def sources(cfg):
    b.check(cfg['implementation']==b.pin(__file__) and cfg['declaration']==b.pin(HERE/'DECLARATION.md') and cfg['run_wrapper']==b.pin(HERE/'RUN.sbatch'),'Handoff source changed')
    b.check(cfg['native_input_adapter']==b.pin(old.__file__) and cfg['parser']==b.pin(rec.__file__) and cfg['controls_adapter']==b.pin(controls.__file__),'Frozen helpers changed')
    parent=rec.load(b.verify(cfg['parent_panel']))
    q=read(cfg['qualified_controls']);rv=read(cfg['control_actual_review']);root=read(cfg['control_root_review']);failed=read(cfg['failed_tight_controls']);b.verify(cfg['cost_accounting'])
    b.check(q['status']=='tighter_local_controls_passed' and q['implementation']==cfg['controls_adapter'] and q['parser']==cfg['parser'] and q['complete']==q['declared']==18 and q['new_complete']==15 and q['reused_complete']==3 and q['all_controls_passed'] and len(q['checks'])==35 and all(q['checks'].values()),'Actual full tighter18/35 required')
    qp=read(q['panel']);b.check(qp['parent_panel']==cfg['parent_panel'] and qp['tier']==old.TIERS['tighter'] and qp['failed_tight_controls']==cfg['failed_tight_controls'],'Tighter controls have different parent/recipe')
    b.check(rv['status']=='passed_tighter_controls_actual_review' and rv['result']==cfg['qualified_controls'] and rv['native_adapter']==cfg['controls_adapter'] and rv['panel']==q['panel'] and rv['admitted_cells']==18 and rv['new_admitted']==15 and rv['reused_admitted']==3 and rv['all35_passed'] and rv['qualification_checks']==q['checks'],'Independent tighter18 admission differs')
    b.check(root['status']=='admitted_fixed18_tighter_local_reference_controls' and root['result']==cfg['qualified_controls'] and root['review']==cfg['control_actual_review'] and root['checks']==35,'Root tighter18 admission differs')
    b.check(failed['status']=='local_controls_failed' and failed['panel']==cfg['parent_panel'] and failed['complete']==18 and len(failed['checks'])==35 and sum(failed['checks'].values())==32 and not failed['all_controls_passed'],'Failed original Tight18 must remain failed')
    b.check(set(q['rows'])=={c['cell_id'] for c in qp['cells']} and len(q['rows'])==18,'Tighter raw denominator differs')
    for row in q['rows'].values():
        b.check(row['status']=='complete','Unadmitted tighter control')
        for key in ('result','input','output','start','release'):b.verify(row[key])
    return parent,q

def definitions(parent):
    return [{k:copy.deepcopy(v) for k,v in c.items() if k!='template_input'} for c in parent['cells'] if c['stage']=='cluster']

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text());parent,q=sources(cfg)
    b.check(cfg['ranks']==RANKS and cfg['nodes']==NODES and cfg['maximum_new_native_calls']==8 and cfg['scout_cell_id']==SCOUT,'Declared handoff scope differs')
    cells=definitions(parent);b.check(len(cells)==8 and SCOUT in {c['cell_id'] for c in cells},'Original eight cluster identities required')
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()
    p=dict(adapter=b.pin(__file__),configuration=b.pin(config),declaration=cfg['declaration'],run_wrapper=cfg['run_wrapper'],native_input_adapter=cfg['native_input_adapter'],parser=cfg['parser'],controls_adapter=cfg['controls_adapter'],**{k:cfg[k] for k in AUTH},assets=parent['assets'],orca=parent['orca'],native_modules=parent['native_modules'],saved_complexes=parent['saved_complexes'],tiers=parent['tiers'],cells=cells,original_inputs={},ranks=RANKS,nodes=NODES,maximum_new_native_calls=8,scout_cell_id=SCOUT,original_scout_cell_id=parent['scout_cell_id'],execution_root=cfg['execution_root'],created_utc=b.now())
    deltas={}
    for c in cells:
        previous=next(x for x in parent['cells'] if x['cell_id']==c['cell_id']);raw=b.verify(previous['template_input']).read_text();actual=deck(p,c,RANKS,3000)
        b.check(actual==raw.replace('%pal nprocs 4 end','%pal nprocs 24 end'),'Only MPI count may change in prepared deck')
        dst=out/'inputs'/(c['cell_id']+'.inp');dst.write_text(actual);c['template_input']=b.pin(dst);p['original_inputs'][c['cell_id']]=previous['template_input'];deltas[c['cell_id']]=dict(original=previous['template_input'],prepared=c['template_input'],only_MPI_count_changed=True)
    pp=b.write(out/'PANEL.json',p);load(pp['path'])
    q=dict(status='passed_cluster_handoff_preparation',implementation=b.pin(__file__),panel=pp,configuration=b.pin(config),declared=8,scout_cell_id=SCOUT,original_scout_cell_id=parent['scout_cell_id'],new_native_calls=0,deck_deltas=deltas,qualified_control_count=18,qualified_control_checks=35,original_Tight_pass_count=32,ranks=RANKS,nodes=NODES,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1)
    rp=b.write(out/'RESULT.json',q);print(json.dumps(dict(result=rp,panel=pp,status=q['status'])));return q

def load(path):
    p=json.loads(Path(path).read_text());cfg=read(p['configuration']);parent,q=sources(cfg)
    b.check(p['adapter']==b.pin(__file__),'Prepared source differs')
    for key in ('declaration','run_wrapper','native_input_adapter','parser','controls_adapter',*AUTH,'execution_root','ranks','nodes','maximum_new_native_calls','scout_cell_id'):b.check(p[key]==cfg[key],'Prepared configuration differs: '+key)
    b.check(p['ranks']==RANKS and p['nodes']==NODES and p['maximum_new_native_calls']==len(p['cells'])==8 and p['scout_cell_id']==SCOUT and p['original_scout_cell_id']==parent['scout_cell_id'],'Prepared finite scope differs')
    for key in ('assets','orca','native_modules','saved_complexes','tiers'):b.check(p[key]==parent[key],'Scientific authority differs: '+key)
    b.check(set(p['original_inputs'])=={c['cell_id'] for c in p['cells']},'Original input inventory differs')
    for c,e in zip(p['cells'],definitions(parent)):
        b.check({k:v for k,v in c.items() if k!='template_input'}==e,'Original cluster definition changed')
        previous=next(x for x in parent['cells'] if x['cell_id']==c['cell_id']);b.check(p['original_inputs'][c['cell_id']]==previous['template_input'],'Original input pin differs')
        b.check(b.verify(c['template_input']).read_text()==deck(p,c,RANKS,3000)==b.verify(previous['template_input']).read_text().replace('%pal nprocs 4 end','%pal nprocs 24 end'),'Prepared deck delta differs')
    return p

def subsets(p):
    ids=[c['cell_id'] for c in p['cells']]
    return dict(scout=[SCOUT],remaining=[x for x in ids if x!=SCOUT])

def release_gate(p,path,r):
    b.check(r['status']=='released_consumed_cluster_subset' and r['implementation']==p['adapter'] and r['panel']==b.pin(path) and r['parser']==p['parser'] and r['declaration']==p['declaration'],'Root source release differs')
    b.check(r['maximum_new_native_calls']==8 and r['ranks']==RANKS and r['nodes']==NODES and r['execution_root']==p['execution_root'],'Released budget/resources differ')
    for k in AUTH:b.check(r[k]==p[k],'Released authority differs: '+k)
    sets=subsets(p);b.check(r['stage'] in sets and r['cell_ids']==sets[r['stage']],'Only exact1/7 subset may execute')
    prep=read(r['preparation']);review=read(r['source_review'])
    b.check(prep['status']=='passed_cluster_handoff_preparation' and prep['implementation']==p['adapter'] and prep['panel']==b.pin(path) and prep['new_native_calls']==0 and prep['declared']==8,'Prepared handoff differs')
    b.check(review['status']=='passed_cluster_handoff_source_review' and review['implementation']==p['adapter'] and review['panel']==b.pin(path) and review['preparation']==r['preparation'],'Independent source review differs')
    if r['stage']=='remaining':
        rv=read(r['scout_actual_review']);root=read(r['scout_root_review']);col=read(r['scout_collection'])
        b.check(col['status']=='cluster_handoff_scout_complete' and col['implementation']==p['adapter'] and col['panel']==b.pin(path) and col['complete']==col['declared']==1 and col['cell_ids']==[SCOUT],'Complete scout collection required')
        rp=col['rows'][SCOUT]['result'];b.check(rv['status']=='passed_cluster_handoff_scout_review' and rv['implementation']==p['adapter'] and rv['panel']==b.pin(path) and rv['collection']==r['scout_collection'] and rv['result']==rp and rv['native_checks_passed'] is True and rv['cost_reviewed'] is True,'Independent actual scout native/cost review required')
        b.check(root['status']=='accepted_cluster_handoff_scout_feasibility' and root['collection']==r['scout_collection'] and root['result']==rp and root['independent_review']==r['scout_actual_review'] and root['cost_and_memory_accepted'] is True and root['scratch_reviewed'] is True,'Root actual cost/memory/scratch acceptance required')
        actual=admit_result(b.verify(rp),next(c for c in p['cells'] if c['cell_id']==SCOUT),p,path)
        b.check(actual==col['rows'][SCOUT],'Scout raw replay differs')
    return r

def resources(env):
    b.check(env.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(env.get('SLURM_NTASKS','0'));nodes=int(env.get('SLURM_JOB_NUM_NODES',env.get('SLURM_NNODES','0')))
    b.check(ranks==RANKS and nodes==NODES and int(env.get('SLURM_CPUS_PER_TASK','1'))==1,'Exactly24 single-thread ranks on one node required')
    mem=int(env.get('SLURM_MEM_PER_NODE','0')) or ranks*int(env.get('SLURM_MEM_PER_CPU','0'));maxcore=math.floor(.75*mem/ranks)
    b.check(mem>0 and maxcore>0,'Positive actual default memory grant required');return ranks,nodes,mem,maxcore

def run(path,release,cell_id,out):
    p=load(path);r=release_gate(p,path,json.loads(Path(release).read_text()));b.check(cell_id in r['cell_ids'],'Cell not released')
    c=next(c for c in p['cells'] if c['cell_id']==cell_id);b.check(out.resolve()==Path(p['execution_root']).resolve()/'tasks'/cell_id/'execution' and not out.exists(),'Unique one-attempt path required')
    ranks,nodes,mem,maxcore=resources(os.environ);out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for k in THREADS:env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    provenance={k:p[k] for k in ('adapter','native_input_adapter','parser')};st=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=ranks,nodes=nodes,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now(),**provenance))
    tick=time.monotonic();code=None;invoked=0
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            invoked=1;code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick;q.update(cell_id=cell_id,cell=c,panel=b.pin(path),release=b.pin(release),start=st,returncode=code,native_program_invocations_attempted=invoked,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),input=b.pin(out/'endpoint.inp'),output=b.pin(out/'endpoint.out') if (out/'endpoint.out').exists() else None,stderr=b.pin(out/'endpoint.err') if (out/'endpoint.err').exists() else None,artifact_sizes={x.name:x.stat().st_size for x in out.iterdir() if x.is_file()},**provenance)
    b.write(out/'RESULT.json',q);print(json.dumps(dict(status=q['status'],cell_id=cell_id,elapsed_seconds=elapsed)));return 0 if q['status']=='complete' else 1

def admit_result(path,c,p,panel):
    path=Path(path);r=json.loads(path.read_text());b.check(path.resolve()==Path(p['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','Native result path differs')
    b.check(r['status']=='complete' and r['returncode']==0 and r['cell_id']==c['cell_id'] and r['cell']==c and r['panel']==b.pin(panel) and all(r['checks'].values()),'Complete native role required')
    b.check(all(r[k]==p[k] for k in ('adapter','native_input_adapter','parser')),'Native source differs');rel=release_gate(p,panel,read(r['release']));b.check(c['cell_id'] in rel['cell_ids'],'Native cell not released')
    st=read(r['start']);b.check(r['start']==b.pin(path.parent/'START.json') and all(st[k]==r[k] for k in ('cell','panel','release','input','job_id','adapter','native_input_adapter','parser')),'Native START differs')
    b.check(st['ranks']==RANKS and st['nodes']==NODES and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/RANKS) and st['maximum_native_calls']==r['native_program_invocations_attempted']==1 and st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native resources/invocation differ')
    b.check(b.verify(r['input']).read_text()==deck(p,c,RANKS,st['maxcore_MiB']),'Executed scientific input differs');b.verify(r['stderr']);q=json.loads(json.dumps(parse(b.verify(r['output']).read_text(),c,r['returncode'])))
    b.check(q['status']=='complete' and all(r[k]==v for k,v in q.items()),'Full raw native replay differs')
    return dict(q,result=b.pin(path),input=r['input'],output=r['output'],start=r['start'],release=r['release'],native_producer=r['adapter'],elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'],artifact_sizes=r['artifact_sizes'])

def contrast(values):
    return values['aquo9_central_Dy']-values['aquo9_central_La']-values['aquo8_central_Dy']+values['aquo8_central_La']

def score(p,rows):
    b.check(set(rows)=={c['cell_id'] for c in p['cells']} and len(rows)==8 and all(x['status']=='complete' for x in rows.values()),'All eight native rows required')
    vals={t:{g:rows[f'cluster__{t}__{g}']['energy_hartree'] for g in old.h.GIDS} for t in old.TIERS}
    cc={t:contrast(v)*b.HA for t,v in vals.items()};saved=p['saved_complexes'];pe=contrast({g:saved[g]['electronic_energy_hartree'] for g in old.h.GIDS})*b.HA;d4=contrast({g:saved[g]['D4_energy_hartree'] for g in old.h.GIDS})*b.HA
    diff=cc['tighter']-cc['tight']
    return dict(raw_CC_C_kcal_mol=cc,PBE0_electronic_C_kcal_mol=pe,original_D4_C_kcal_mol=d4,PBE0_D4_C_kcal_mol=pe+d4,original_PBE0_final_C_kcal_mol=contrast({g:saved[g]['final_energy_hartree'] for g in old.h.GIDS})*b.HA,tighter_minus_tight_C_kcal_mol=diff,target_threshold_sensitivity_passed=abs(diff)<=.25,endpoint_tighter_minus_tight_kcal_mol={g:(vals['tighter'][g]-vals['tight'][g])*b.HA for g in old.h.GIDS},CC_minus_PBE0_electronic_C_kcal_mol={t:v-pe for t,v in cc.items()},CC_minus_PBE0_D4_C_kcal_mol={t:v-pe-d4 for t,v in cc.items()},cutoff_sensitivity_limit_kcal_mol=.25)

def collect(path,inventory,out,stage):
    tick,cpu=time.monotonic(),time.process_time();p=load(path);inv=json.loads(Path(inventory).read_text());cells=[c for c in p['cells'] if stage=='full' or c['cell_id']==SCOUT]
    b.check(set(inv)=={c['cell_id'] for c in cells},'Exact scout/full inventory required');out.mkdir(parents=True,exist_ok=False);rows={};good={}
    for c in cells:
        cid=c['cell_id']
        try:b.check(inv[cid] is not None,'Native result missing');q=admit_result(inv[cid],c,p,path);rows[cid]=q;good[cid]=q
        except Exception as e:rows[cid]=dict(status='unavailable',path=inv.get(cid),error=repr(e))
    q=dict(status='cluster_handoff_unavailable',implementation=p['adapter'],panel=b.pin(path),inventory=b.pin(inventory),stage=stage,cell_ids=[c['cell_id'] for c in cells],declared=len(cells),complete=len(good),rows=rows,qualified_controls=p['qualified_controls'],control_actual_review=p['control_actual_review'],control_root_review=p['control_root_review'],failed_tight_controls=p['failed_tight_controls'],original_Tight_controls_passed=False,tighter_controls_passed=True,physical_accuracy_qualified=False,cluster_CP_or_B_qualified=False,new_native_calls=0)
    if len(good)==len(cells):
        if stage=='scout':q['status']='cluster_handoff_scout_complete'
        else:q.update(score(p,good));q['status']='cluster_handoff_complete' if q['target_threshold_sensitivity_passed'] else 'cluster_handoff_cutoff_unqualified'
    q.update(native_allocated_CPU_seconds_admitted_rows=sum(v['allocated_cpu_seconds'] for v in good.values()),cost_scope='Application costs of admitted rows only; root accounts every scheduler allocation including failures and collectors, without step duplication.',wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','declared','complete','new_native_calls')}));return q

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',required=True,type=Path)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',required=True,type=Path)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--stage',choices=['scout','full'],required=True);a.add_argument('--out',required=True,type=Path)
    args=ap.parse_args()
    if args.mode=='prepare':prepare(args.config,args.out)
    elif args.mode=='run':sys.exit(run(args.panel,args.release,args.cell,args.out))
    else:q=collect(args.panel,args.inventory,args.out,args.stage);sys.exit(0 if q['status'] in ('cluster_handoff_scout_complete','cluster_handoff_complete') else 1)
