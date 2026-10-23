#!/usr/bin/env python3
"""Versioned 14/1/7 continuation; unchanged prepared native decks and scoring."""
import argparse,importlib.util,json,math,os,socket,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('reviewed_saved_recovery',HERE.parent/'parser_recovery_v1/RECOLLECT.py')
rec=importlib.util.module_from_spec(s);s.loader.exec_module(rec)
old=rec.old;b=old.b;read=old.read;deck=old.deck;parse=rec.parse;control_score=old.control_score
PARSER_SHA='a25d8a7dffacb7a7a2a663fecf8928ee844ed5c181db433e8a87610abaf59fbd'

def load(path):
    b.check(b.pin(rec.__file__)['sha256']==PARSER_SHA,'Reviewed parser differs')
    return rec.load(path)

def subsets(p):
    ids=[c['cell_id'] for c in p['cells']]
    return dict(controls_remaining=[x for x in ids if x.startswith('control__') and x not in old.CONTROL_SCOUTS],
        cluster_scout=[old.SCOUT],cluster_remaining=[x for x in ids if x.startswith('cluster__') and x!=old.SCOUT])

def admit_recovered(path,c,p,panel):
    r=json.loads(Path(path).read_text())
    b.check(r['status']=='saved_native_recovery_complete' and r['implementation']==b.pin(rec.__file__)
        and r['declaration']==b.pin(rec.HERE/'DECLARATION.md') and r['new_native_calls']==0,'Saved recovery source/status differs')
    q=rec.admit_original(b.verify(r['original_result']),c,p,panel)
    canonical=json.loads(json.dumps(q))
    b.check(all(r[k]==v for k,v in canonical.items() if k!='status'),'Complete saved recovery replay differs')
    return dict(q,status='complete',result=b.pin(path),native_result=r['original_result'],native_producer=r['original_adapter'])

def release_contract(p,path,r):
    b.check(r['status']=='approved_coordination_local_cc_subset' and r['implementation']==b.pin(__file__)
        and r['native_input_adapter']==p['adapter'] and r['parser']==b.pin(rec.__file__) and r['panel']==b.pin(path),'Versioned release source differs')
    ids=[c['cell_id'] for c in p['cells']]
    b.check(r['cell_ids']==ids and r['maximum_native_calls']==26 and r['ranks']==4,'Released finite denominator differs')
    pre=read(r['preparation']);b.check(pre['status']=='passed_static_preflight' and pre['panel']==b.pin(path),'Original prepared decks unavailable')
    prep=read(r['runtime_preparation']);b.check(prep['status']=='passed_runtime_handoff_preparation' and prep['implementation']==b.pin(__file__)
        and prep['parser']==b.pin(rec.__file__) and prep['panel']==b.pin(path) and prep['all_26_original_decks_identical'] and prep['new_native_calls']==0,'Versioned preparation differs')
    review=read(r['independent_review']);b.check(review['status']=='passed_runtime_handoff_review' and review['implementation']==b.pin(__file__)
        and review['parser']==b.pin(rec.__file__) and review['panel']==b.pin(path),'Versioned independent review differs')
    recovery_review=read(r['parser_recovery_review']);b.check(recovery_review['status']=='passed_saved_parser_recovery_review'
        and recovery_review['implementation']==b.pin(rec.__file__) and recovery_review['panel']==b.pin(path),'Saved parser review differs')
    stage=r['stage'];sets=subsets(p)
    b.check(stage in sets and r['released_cell_ids']==sets[stage],'Only the exact uncomputed 14/1/7 subsets may execute')
    if stage=='controls_remaining':
        sr=read(r['control_scout_review'])
        b.check(sr['status']=='passed_control_scout_review' and sr['panel']==b.pin(path) and sr['parser']==b.pin(rec.__file__)
            and sr['parser_recovery_review']==r['parser_recovery_review'] and set(sr['results'])==set(sr['recoveries'])==set(old.CONTROL_SCOUTS),'Four recovered original scouts required')
        for cid,pin in sr['recoveries'].items():
            q=admit_recovered(b.verify(pin),next(c for c in p['cells'] if c['cell_id']==cid),p,path)
            b.check(q['original_result']==sr['results'][cid],'Scout recovery/producer mismatch')
    if stage.startswith('cluster_'):
        qual=read(r['control_qualification'])
        b.check(qual['status']=='local_controls_passed' and qual['panel']==b.pin(path) and qual['complete']==qual['declared']==18
            and qual['all_controls_passed'] and len(qual['checks'])==35 and all(qual['checks'].values())
            and qual['implementation']==b.pin(__file__),'All eighteen unchanged method controls required')
    if stage=='cluster_remaining':
        scout=read(r['cluster_scout_review'])
        b.check(scout['status']=='passed_cluster_scout_review' and scout['panel']==b.pin(path) and scout['cost_and_memory_accepted'] is True,'Actual largest scout cost/native admission missing')
        c=next(c for c in p['cells'] if c['cell_id']==old.SCOUT);admit_result(b.verify(scout['result']),c,p,path)
    return r

def admit_result(path,c,p,panel):
    path=Path(path);r=json.loads(path.read_text())
    if r.get('status')=='saved_native_recovery_complete':return admit_recovered(path,c,p,panel)
    b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Native result not admitted')
    b.check(r['cell_id']==c['cell_id'] and r['cell']==c and r['panel']==b.pin(panel) and r['adapter']==b.pin(__file__)
        and r['native_input_adapter']==p['adapter'] and r['parser']==b.pin(rec.__file__),'Versioned cell/source differs')
    st=read(r['start']);rel=read(r['release']);release_contract(p,panel,rel)
    b.check(r['start']==b.pin(path.parent/'START.json') and st['cell']==c and st['panel']==r['panel'] and st['release']==r['release']
        and st['input']==r['input'] and st['adapter']==r['adapter'] and st['native_input_adapter']==r['native_input_adapter'] and st['parser']==r['parser'],'START lineage differs')
    b.check(c['cell_id'] in rel['released_cell_ids'] and path.resolve()==Path(rel['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','One-attempt path/subset differs')
    b.check(st['ranks']==4 and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/4),'Resources differ')
    b.check(st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and r['native_program_invocations_attempted']==1,'Invocation differs')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,st['maxcore_MiB']),'Executed input differs');b.verify(r['stderr'])
    q=parse(b.verify(r['output']).read_text(),c,r['returncode'])
    b.check(q['status']=='complete' and all(r[k]==v for k,v in json.loads(json.dumps(q)).items()),'Complete serialized raw replay differs')
    return dict(q,result=b.pin(path),native_result=b.pin(path),native_producer=r['adapter'],output=r['output'],input=r['input'],start=r['start'],release=r['release'],elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'])

def prepare(path,out):
    tick,cpu=time.monotonic(),time.process_time();p=load(path)
    for c in p['cells']:b.check(b.verify(c['template_input']).read_text()==deck(p,c,4,3000),'Original prepared deck changed')
    out.mkdir(parents=True,exist_ok=False)
    q=dict(status='passed_runtime_handoff_preparation',implementation=b.pin(__file__),parser=b.pin(rec.__file__),native_input_adapter=p['adapter'],panel=b.pin(path),
        declaration=b.pin(HERE/'DECLARATION.md'),subsets=subsets(p),all_26_original_decks_identical=True,original_decks={c['cell_id']:c['template_input'] for c in p['cells']},
        new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    pin=b.write(out/'RESULT.json',q);print(json.dumps(dict(result=pin,status=q['status'])));return q

def run(path,release,cell_id,out):
    p=load(path);r=release_contract(p,path,json.loads(Path(release).read_text()));b.check(cell_id in r['released_cell_ids'],'Cell not released')
    b.check(out.resolve()==Path(r['execution_root']).resolve()/'tasks'/cell_id/'execution' and not out.exists(),'One-attempt path differs/existing')
    c=next(c for c in p['cells'] if c['cell_id']==cell_id);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,4,maxcore));env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    provenance=dict(adapter=b.pin(__file__),native_input_adapter=p['adapter'],parser=b.pin(rec.__file__))
    st=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=4,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now(),**provenance))
    tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick;q.update(cell_id=cell_id,cell=c,panel=b.pin(path),**provenance,release=b.pin(release),start=st,returncode=code,
        native_program_invocations_attempted=1,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*4,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),
        output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),artifact_sizes={x.name:x.stat().st_size for x in out.iterdir() if x.is_file()})
    b.write(out/'RESULT.json',q);print(json.dumps(dict(status=q['status'],cell_id=cell_id,seconds=elapsed)));return 0 if q['status']=='complete' else 1

def collect(path,inventory,out,stage):
    # Same arithmetic and complete-denominator logic as the frozen producer;
    # admission alone dispatches explicit recovered versus new producer records.
    tick,cpu=time.monotonic(),time.process_time();p=load(path);inv=json.loads(Path(inventory).read_text());cells=[c for c in p['cells'] if stage=='full' or c['stage']=='controls']
    b.check(set(inv)=={c['cell_id'] for c in cells},'Exact result inventory including unavailable slots required')
    out.mkdir(parents=True,exist_ok=False);rows={c['cell_id']:dict(status='unavailable',error=None) for c in cells};good={}
    for c in cells:
        try:
            b.check(inv[c['cell_id']] is not None,'Missing native result');admitted=admit_result(inv[c['cell_id']],c,p,path);rows[c['cell_id']]=dict(admitted,status='complete');good[c['cell_id']]=admitted
        except Exception as e:rows[c['cell_id']]=dict(status='unavailable',error=repr(e),path=inv[c['cell_id']])
    q=dict(status='unavailable',implementation=b.pin(__file__),parser=b.pin(rec.__file__),native_input_adapter=p['adapter'],panel=b.pin(path),inventory=b.pin(inventory),stage=stage,declared=len(cells),complete=len(good),rows=rows,all_controls_passed=False,new_native_calls=0)
    if all(c['cell_id'] in good for c in p['cells'] if c['stage']=='controls'):
        sc=control_score(p,good);q.update(sc);q['status']='local_controls_passed' if stage=='controls' and sc['all_controls_passed'] else 'local_controls_failed'
    if stage=='full' and len(good)==26:
        def contrast(e):return e['aquo9_central_Dy']-e['aquo9_central_La']-e['aquo8_central_Dy']+e['aquo8_central_La']
        values={t:{g:good[f'cluster__{t}__{g}']['energy_hartree'] for g in old.h.GIDS} for t in old.TIERS}
        cc={t:contrast(v)*b.HA for t,v in values.items()};saved=p['saved_complexes'];pe=contrast({g:saved[g]['electronic_energy_hartree'] for g in old.h.GIDS})*b.HA;d4=contrast({g:saved[g]['D4_energy_hartree'] for g in old.h.GIDS})*b.HA
        diff=cc['tighter']-cc['tight'];q.update(status='local_coordination_method_sensitivity_complete',raw_CC_C_kcal_mol=cc,PBE0_electronic_C_kcal_mol=pe,original_D4_C_kcal_mol=d4,PBE0_D4_C_kcal_mol=pe+d4,
            tighter_minus_tight_C_kcal_mol=diff,target_threshold_sensitivity_passed=abs(diff)<=.25,
            endpoint_tighter_minus_tight_kcal_mol={g:(values['tighter'][g]-values['tight'][g])*b.HA for g in old.h.GIDS},
            CC_minus_PBE0_D4_C_kcal_mol={t:v-pe-d4 for t,v in cc.items()},method_sensitivity_controls_passed=q['all_controls_passed'] and abs(diff)<=.25)
    q.update(physical_accuracy_qualified=False,cluster_CP_or_B_qualified=False,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    q['native_allocated_CPU_seconds']=sum(v['allocated_cpu_seconds'] for v in good.values())
    b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','complete','declared','new_native_calls')}));return q

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--panel',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--stage',choices=['controls','full'],required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.panel,a.out)
    elif a.action=='run':return run(a.panel,a.release,a.cell,a.out)
    else:collect(a.panel,a.inventory,a.out,a.stage)
if __name__=='__main__':sys.exit(main())
