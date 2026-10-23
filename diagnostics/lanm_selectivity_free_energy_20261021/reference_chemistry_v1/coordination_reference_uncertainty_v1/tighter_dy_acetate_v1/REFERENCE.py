#!/usr/bin/env python3
"""Three tighter Dy-acetate roles; frozen deck/parser helpers, root releases."""
import argparse, copy, importlib.util, json, math, os, re, socket, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('qualified_local_cc_recovery',HERE.parent/'parser_recovery_v1/RECOLLECT.py')
rec=importlib.util.module_from_spec(s);s.loader.exec_module(rec)
old=rec.old;b=old.b;read=old.read;deck=old.deck
PARSER_SHA='a25d8a7dffacb7a7a2a663fecf8928ee844ed5c181db433e8a87610abaf59fbd'
ROLES=['complex','ion_ghost_ligand','ligand_ghost_ion']
TRIPLES_KEYS=['TCutTNO','TCutMP2Pairs','TCutDOStrong','TCutMKNStrong','TCutDOWeak','TCutMKNWeak',
              'TCutFock','TCutWeak','TCutWeakGuess','TCutIter','TNOScaleS','TNOScaleW','ETol']

def definitions(parent):
    result=[]
    for role in ROLES:
        cid='control__cc__acetate_bi__Dy__'+role
        c=copy.deepcopy(next(c for c in parent['cells'] if c['cell_id']==cid))
        c.pop('template_input');c.update(cell_id='tighter__cc__acetate_bi__Dy__'+role,
            source_tight_cell_id=cid,tier='tighter',stage='attribution_scout')
        result.append(c)
    return result

def triples_settings(text):
    native=text.split('****END OF INPUT****')[-1]
    sections=native.split('DLPNO BASED TRIPLES CORRECTION')
    part=sections[1].split('Fragment selection',1)[0] if len(sections)==2 else ''
    values={k:rec.numbers(part,k,r'\.{2,}') for k in TRIPLES_KEYS}
    passed=all(len(v)==1 and v[0] is not None and math.isfinite(v[0]) and v[0]>0 for v in values.values())
    return dict(values=values,raw_section=part,complete_and_finite=passed)

def parse(text,c,returncode):
    q=rec.parse(text,c,returncode);settings=triples_settings(text)
    q['checks']['complete_triples_settings_readback']=settings['complete_and_finite']
    q['actual_triples_settings']=settings
    q['status']='complete' if all(q['checks'].values()) else 'failed_admission'
    return q

def load(path):
    b.check(b.pin(rec.__file__)['sha256']==PARSER_SHA,'Frozen parser differs')
    p=json.loads(Path(path).read_text())
    b.check(p['adapter']==b.pin(__file__) and p['parser']==b.pin(rec.__file__)
        and p['declaration']==b.pin(HERE/'DECLARATION.md') and p['run_wrapper']==b.pin(HERE/'RUN.sbatch'),'Source binding differs')
    for key in ('configuration','declaration','parent_panel','failed_controls','saved_review','components','interpretation','run_wrapper'):b.verify(p[key])
    cfg=read(p['configuration'])
    b.check(cfg['implementation']==p['adapter'] and cfg['parser']==p['parser'] and all(cfg[k]==p[k]
        for k in ('parent_panel','failed_controls','saved_review','components','interpretation','execution_root')),'Configuration join differs')
    parent=rec.load(b.verify(p['parent_panel']))
    b.check(p['native_input_adapter']==parent['adapter'] and p['assets']==parent['assets'] and p['orca']==parent['orca'],'Original native assets differ')
    expected=definitions(parent)
    b.check(len(p['cells'])==p['maximum_native_calls']==3 and p['ranks']==8 and p['tier']==old.TIERS['tighter'],'Finite scope differs')
    for c,e in zip(p['cells'],expected):
        b.check({k:v for k,v in c.items() if k!='template_input'}==e,'Literal cell differs')
        b.check(b.verify(c['template_input']).read_text()==deck(p,c,8,3000),'Prepared deck differs')
    failed=read(p['failed_controls']);review=read(p['saved_review']);comp=read(p['components'])
    b.check(failed['status']=='local_controls_failed' and failed['complete']==18 and not failed['all_controls_passed'],'Original failed controls must remain failed')
    b.check(review['status']=='passed_saved_control_collection_review' and review['collector_result']==p['failed_controls']
        and review['panel']==p['parent_panel'] and not review['all_method_controls_passed'],'Saved raw review differs')
    b.check(comp['status']=='saved_component_accounting_complete' and comp['review']==p['saved_review']
        and comp['collector_result']==p['failed_controls'] and comp['panel']==p['parent_panel'],'Component authority differs')
    for c in p['cells']:
        v=comp['rows'][c['source_tight_cell_id']];q=failed['rows'][c['source_tight_cell_id']]
        b.check(v['local']['output']==q['output'] and v['canonical']['output']==parent['controls'][c['canonical_id']]['output'],'Component native join differs')
    b.check(p['baseline_components']=={c['source_tight_cell_id']:comp['rows'][c['source_tight_cell_id']] for c in p['cells']},'Baseline components changed')
    return p

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__) and cfg['parser']==b.pin(rec.__file__),'Preparation source differs')
    parent=rec.load(b.verify(cfg['parent_panel']));comp=read(cfg['components']);cells=definitions(parent)
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()
    p=dict(adapter=b.pin(__file__),parser=b.pin(rec.__file__),native_input_adapter=parent['adapter'],configuration=b.pin(config),
        declaration=b.pin(HERE/'DECLARATION.md'),run_wrapper=b.pin(HERE/'RUN.sbatch'),
        **{k:cfg[k] for k in ('parent_panel','failed_controls','saved_review','components','interpretation','execution_root')},
        assets=parent['assets'],orca=parent['orca'],cells=cells,maximum_native_calls=3,ranks=8,tier=old.TIERS['tighter'],
        baseline_components={c['source_tight_cell_id']:comp['rows'][c['source_tight_cell_id']] for c in cells},created_utc=b.now())
    deltas={}
    for c in cells:
        oldc=next(x for x in parent['cells'] if x['cell_id']==c['source_tight_cell_id'])
        before=b.verify(oldc['template_input']).read_text();after=deck(p,c,8,3000)
        expected=before.replace('%pal nprocs 4 end','%pal nprocs 8 end').replace(' TCutPairs 1e-05\n',' TCutPairs 1e-06\n').replace(' TCutPNO 1e-07\n',' TCutPNO 1e-08\n').replace(' TCutMKN 0.001\n',' TCutMKN 0.0001\n')
        b.check(after==expected,'Unexpected scientific/native input delta')
        dest=out/'inputs'/(c['cell_id']+'.inp');dest.write_text(after);c['template_input']=b.pin(dest)
        deltas[c['cell_id']]=dict(original_input=oldc['template_input'],new_input=c['template_input'],only_declared_three_cutoffs_and_MPI_changed=True)
    pp=b.write(out/'PANEL.json',p);load(pp['path'])
    q=dict(status='passed_tighter_dy_acetate_preparation',implementation=b.pin(__file__),parser=b.pin(rec.__file__),panel=pp,
        configuration=b.pin(config),declared=3,input_deltas=deltas,template_resources_only=dict(ranks=8,memory_MiB=32000,maxcore_MiB=3000),
        new_native_calls=0,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,threads=1)
    pin=b.write(out/'RESULT.json',q);print(json.dumps(dict(result=pin,panel=pp)));return q

def release_contract(p,path,r):
    b.check(r['status']=='approved_tighter_dy_acetate_attribution' and r['implementation']==b.pin(__file__)
        and r['parser']==p['parser'] and r['panel']==b.pin(path) and r['declaration']==p['declaration'],'Release source differs')
    b.check(r['cell_ids']==[c['cell_id'] for c in p['cells']] and r['maximum_native_calls']==3 and r['ranks']==8
        and r['execution_root']==p['execution_root'],'Release finite inventory/resources differ')
    prep=read(r['preparation']);review=read(r['source_review'])
    b.check(prep['status']=='passed_tighter_dy_acetate_preparation' and prep['implementation']==p['adapter'] and prep['panel']==b.pin(path)
        and prep['declared']==3 and prep['new_native_calls']==0,'Actual preparation missing')
    b.check(review['status']=='passed_tighter_dy_acetate_source_review' and review['implementation']==p['adapter']
        and review['panel']==b.pin(path) and review['preparation']==r['preparation'],'Root source/input review missing')
    b.check(r['full_eighteen_qualified'] is False and r['cluster_release_authorized'] is False,'Attribution scope cannot release clusters')
    return r

def resources(env):
    ranks=int(env.get('SLURM_NTASKS','0'))
    b.check(env.get('SLURM_JOB_ID','').isdigit() and ranks==8 and int(env.get('SLURM_CPUS_PER_TASK','1'))==1,'Eight single-thread Slurm ranks required')
    mem=int(env.get('SLURM_MEM_PER_NODE','0')) or ranks*int(env.get('SLURM_MEM_PER_CPU','0'))
    b.check(mem>0,'Actual positive memory grant required');maxcore=math.floor(.75*mem/ranks)
    b.check(maxcore>0,'Actual grant too small');return ranks,mem,maxcore

def run(path,release,cell_id,out):
    p=load(path);r=release_contract(p,path,json.loads(Path(release).read_text()))
    b.check(cell_id in r['cell_ids'],'Unknown role');c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    b.check(out.resolve()==Path(p['execution_root']).resolve()/'tasks'/cell_id/'execution' and not out.exists(),'One-attempt path differs or exists')
    ranks,mem,maxcore=resources(os.environ);out.mkdir(parents=True,exist_ok=False)
    (out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];provenance=dict(adapter=p['adapter'],native_input_adapter=p['native_input_adapter'],parser=p['parser'])
    st=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],
        host=socket.gethostname(),ranks=ranks,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now(),**provenance))
    tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick;q.update(cell_id=cell_id,cell=c,panel=b.pin(path),**provenance,release=b.pin(release),start=st,returncode=code,
        native_program_invocations_attempted=1,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),
        output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),artifact_sizes={x.name:x.stat().st_size for x in out.iterdir() if x.is_file()})
    b.write(out/'RESULT.json',q);print(json.dumps(dict(status=q['status'],cell_id=cell_id,seconds=elapsed)));return 0 if q['status']=='complete' else 1

def admit_result(path,c,p,panel):
    path=Path(path);r=json.loads(path.read_text())
    b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Native result unavailable')
    b.check(r['cell']==c and r['cell_id']==c['cell_id'] and r['adapter']==p['adapter'] and r['native_input_adapter']==p['native_input_adapter']
        and r['parser']==p['parser'] and r['panel']==b.pin(panel),'Native producer/cell differs')
    st=read(r['start']);rel=release_contract(p,panel,read(r['release']))
    b.check(r['start']==b.pin(path.parent/'START.json') and st['cell']==c and st['panel']==r['panel'] and st['release']==r['release']
        and st['input']==r['input'] and all(st[k]==r[k] for k in ('adapter','native_input_adapter','parser')),'START lineage differs')
    b.check(path.resolve()==Path(p['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','One-attempt result path differs')
    b.check(st['ranks']==8 and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/8)
        and st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and r['native_program_invocations_attempted']==1,'Resources/invocation differ')
    b.check(b.verify(r['input']).read_text()==deck(p,c,8,st['maxcore_MiB']),'Executed deck differs');b.verify(r['stderr'])
    q=parse(b.verify(r['output']).read_text(),c,r['returncode'])
    b.check(q['status']=='complete' and all(r[k]==v for k,v in json.loads(json.dumps(q)).items()),'Complete raw replay differs')
    return dict(q,result=b.pin(path),output=r['output'],input=r['input'],start=r['start'],release=r['release'],native_producer=r['adapter'],
        elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'])

def comparison(p,rows):
    b.check(set(rows)=={c['cell_id'] for c in p['cells']},'All three CP roles required')
    sets={name:{k:0. for k in ('CCSD','triples','total')} for name in ('canonical','tight','tighter')};hf={};aux=[]
    for c in p['cells']:
        q=rows[c['cell_id']];previous=p['baseline_components'][c['source_tight_cell_id']];sign=1 if c['role']=='complex' else -1
        new=dict(CCSD=q['CCSD_energies_hartree'][0],triples=q['triples_hartree'][0],total=q['energy_hartree'])
        for name,values in [('tighter',new),('tight',previous['local']['values']),('canonical',previous['canonical']['values'])]:
            for k in sets[name]:sets[name][k]+=sign*float(values[k])*b.HA
        hf[c['cell_id']]=q['RHF_energies_hartree'][0]-float(previous['canonical']['values']['RHF']);aux.append(q['aux_C'])
    errors={name:{k:sets[name][k]-sets['canonical'][k] for k in sets[name]} for name in ('tight','tighter')}
    aux_ok=all(x['contractions']==aux[0]['contractions'] and x['counts']==aux[0]['counts'] and x['totals']==aux[0]['totals'] for x in aux)
    checks={**{'RHF__'+k:abs(v)<=1e-6 for k,v in hf.items()},'same_parent_AuxC':aux_ok,'original_0p25_CP_diagnostic':abs(errors['tighter']['total'])<=.25}
    return dict(CP_components_kcal_mol=sets,CP_errors_vs_canonical_kcal_mol=errors,
        tighter_minus_tight_CP_components_kcal_mol={k:sets['tighter'][k]-sets['tight'][k] for k in sets['tight']},
        RHF_errors_hartree=hf,diagnostic_checks=checks,diagnostic_limit_passed=checks['original_0p25_CP_diagnostic'],all_scout_diagnostics_passed=all(checks.values()))

def collect(path,inventory,out):
    tick,cpu=time.monotonic(),time.process_time();p=load(path);inv=json.loads(Path(inventory).read_text())
    b.check(set(inv)=={c['cell_id'] for c in p['cells']},'Exact three result slots required');out.mkdir(parents=True,exist_ok=False)
    rows={};good={}
    for c in p['cells']:
        try:
            b.check(inv[c['cell_id']] is not None,'Missing native result');q=admit_result(inv[c['cell_id']],c,p,path)
            rows[c['cell_id']]=dict(q,status='complete');good[c['cell_id']]=q
        except Exception as e:rows[c['cell_id']]=dict(status='unavailable',path=inv[c['cell_id']],error=repr(e))
    q=dict(status='unavailable',implementation=p['adapter'],parser=p['parser'],panel=b.pin(path),inventory=b.pin(inventory),declared=3,complete=len(good),rows=rows,
        full_eighteen_qualified=False,cluster_release_authorized=False,physical_accuracy_qualified=False,new_native_calls=0)
    if len(good)==3:q.update(status='dy_acetate_tighter_attribution_complete',**comparison(p,good))
    q.update(wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,native_allocated_CPU_seconds=sum(x['allocated_cpu_seconds'] for x in good.values()))
    b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','declared','complete','new_native_calls')}));return q

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',required=True,type=Path)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',required=True,type=Path)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--out',required=True,type=Path)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='run':return run(a.panel,a.release,a.cell,a.out)
    else:collect(a.panel,a.inventory,a.out)
if __name__=='__main__':sys.exit(main())
