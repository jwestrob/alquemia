#!/usr/bin/env python3
"""Exact legacy-scout JSON repair; seven original native continuations only."""
import argparse,importlib.util,json,math,os,socket,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('original_canonical_rotation',HERE.parent/'dry_acetate_correlated_rotation_v1/REFERENCE.py')
old=importlib.util.module_from_spec(sp);sp.loader.exec_module(old)
b=old.b;read=old.read;deck=old.deck;parse=old.parse;score=old.score
ALL=old.ALL;SCOUT=old.SCOUT;REMAINING=old.REMAINING;METALS=old.METALS
ORIGINAL_SHA='21679bf648b546117b459d91bc728b497f3b56941a078c77c3d336a846f82634'

def canonical(value):return json.loads(json.dumps(value,allow_nan=False))
def raw_equal(saved,replay):
    for k,v in replay.items():b.check(canonical(saved[k])==canonical(v),'Stored/raw parse differs: '+k)

def load(path):
    b.check(b.pin(old.__file__)['sha256']==ORIGINAL_SHA,'Original source changed');p=old.load(path);inputs=read(HERE/'INPUTS.json');old.verify_tree(inputs)
    b.check(inputs['legacy_implementation']==b.pin(old.__file__) and inputs['panel']==b.pin(path),'Legacy source/panel differs')
    legacy=read(b.verify(inputs['legacy_scout_result']))
    for field,key in [('adapter','legacy_implementation'),('release','legacy_scout_release'),('output','legacy_scout_output'),('input','legacy_scout_input'),('start','legacy_scout_start')]:b.check(legacy[field]==inputs[key],'Legacy '+field+' changed')
    b.check(legacy['cell_id']==SCOUT and legacy['panel']==inputs['panel'],'Exact legacy scout differs')
    p['_recovery_inputs']=inputs;return p

def release_gate(p,release_pin,cell_id):
    r=read(b.verify(release_pin));pp=b.pin(p['_path'])
    b.check(r['status']=='released_dry_acetate_correlated_rotation_v2_remaining' and r['stage']=='remaining' and r['implementation']==b.pin(__file__) and r['panel']==pp and r['recovery_inputs']==b.pin(HERE/'INPUTS.json') and r['declaration']==b.pin(HERE/'DECLARATION.md'),'V2 source release differs')
    b.check(cell_id in REMAINING and r['cell_ids']==list(REMAINING) and r['maximum_total_native_calls']==8 and r['maximum_stage_calls']==7 and r['ranks']==4 and r['nodes']==1,'Only exact remaining seven may run')
    b.check(r['scout_result']==p['_recovery_inputs']['legacy_scout_result'],'Wrong legacy scout result')
    scout=admit(p,b.verify(r['scout_result']));sr=read(b.verify(r['scout_review']))
    b.check(sr['status']=='passed_dry_acetate_correlated_rotation_scout_review' and sr['result']==r['scout_result'] and sr['implementation']==b.pin(__file__) and sr['panel']==pp and sr['cost_reviewed'] is True,'Actual corrected scout/cost review missing')
    first=read(b.verify(scout['release']));b.check(first['execution_root']==r['execution_root'],'Continuation moved original execution root')
    root=Path(r['execution_root']).resolve();b.check(r['task_output_directories']=={cid:str(root/'tasks'/cid/'execution') for cid in REMAINING},'Immutable seven output directories differ')
    review=read(b.verify(r['source_review']));b.check(review['status']=='passed_dry_acetate_correlated_rotation_v2_source_review' and review['implementation']==b.pin(__file__) and review['panel']==pp and review['recovery_inputs']==r['recovery_inputs'],'V2 source review differs')
    return r

def admit(p,filename):
    r=read(filename);cid=r['cell_id'];c=next(c for c in p['cells'] if c['cell_id']==cid);pp=b.pin(p['_path'])
    if cid==SCOUT:
        inputs=p['_recovery_inputs'];b.check(b.pin(filename)==inputs['legacy_scout_result'],'Only exact original scout may be reused')
        expected=inputs['legacy_implementation'];release=old.release_gate(p,r['release'],cid)
        b.check(r['output']==inputs['legacy_scout_output'] and r['input']==inputs['legacy_scout_input'] and r['release']==inputs['legacy_scout_release'],'Legacy raw pins differ')
    else:
        b.check(cid in REMAINING,'Unknown continuation role');expected=b.pin(__file__);release=release_gate(p,r['release'],cid)
    b.check(r['status']=='complete' and r['adapter']==expected and r['panel']==pp and r['cell']==c and r['new_native_calls']==1,'Native result identity/status differs')
    directory=Path(filename).resolve().parent;b.check(str(directory)==release['task_output_directories'][cid],'Native directory differs')
    for key,name in [('input','endpoint.inp'),('output','endpoint.out'),('stderr','endpoint.err'),('start','START.json')]:b.check(r[key]==b.pin(directory/name),'Native artifact join differs: '+key)
    s=read(b.verify(r['start']));b.check(s['adapter']==r['adapter'] and s['cell']==c and s['panel']==pp and s['release']==r['release'] and s['input']==r['input'] and s['job_id']==r['job_id'],'START identity differs')
    b.check(s['ranks']==4 and s['nodes']==1 and s['maxcore_MiB']==math.floor(.75*s['memory_MiB']/4) and s['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native resources/command differ')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,s['maxcore_MiB']),'Executed input differs')
    replay=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(replay['status']=='complete','Raw canonical admission failed');raw_equal(r,replay)
    b.check(math.isfinite(r['elapsed_seconds']) and r['elapsed_seconds']>=0 and r['allocated_cpu_seconds']==4*r['elapsed_seconds'],'Native cost differs');return r

def prepare(path,out):
    tick,cpu=time.monotonic(),time.process_time();p=load(path);v=old.validate(path);r=admit(p,b.verify(p['_recovery_inputs']['legacy_scout_result']))
    out.mkdir(parents=True,exist_ok=False);q=dict(status='passed_canonical_rotation_v2_preparation',implementation=b.pin(__file__),panel=b.pin(path),recovery_inputs=b.pin(HERE/'INPUTS.json'),declaration=b.pin(HERE/'DECLARATION.md'),all_original8_templates_unchanged=True,original_validation=v,legacy_scout_result=p['_recovery_inputs']['legacy_scout_result'],legacy_raw_status=r['status'],new_native_calls=0,remaining_ids=list(REMAINING),cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1));b.write(out/'RESULT.json',q);print(json.dumps(dict(result=b.pin(out/'RESULT.json'),status=q['status'])))

def run(path,cell_id,release_path,out):
    p=load(path);rp=b.pin(release_path);r=release_gate(p,rp,cell_id);c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    b.check(str(out.resolve())==r['task_output_directories'][cell_id],'Output differs from release')
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));nodes=int(os.environ.get('SLURM_JOB_NUM_NODES',os.environ.get('SLURM_NNODES','0')))
    b.check(ranks==4 and nodes==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks on one node required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant missing');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];pp=b.pin(path)
    start={'adapter':b.pin(__file__),'cell':c,'panel':pp,'release':rp,'input':b.pin(out/'endpoint.inp'),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'nodes':nodes,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()}
    sp=b.write(out/'START.json',start);tick=time.monotonic();returncode=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
        returncode=proc.returncode;result=parse((out/'endpoint.out').read_text(errors='replace'),c,returncode)
    except Exception as exc:result={'status':'failed_admission','error':repr(exc)}
    elapsed=time.monotonic()-tick
    result.update(cell_id=cell_id,cell=c,panel=pp,adapter=b.pin(__file__),release=rp,start=sp,input=start['input'],output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),returncode=returncode,job_id=start['job_id'],elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,new_native_calls=1,finished_utc=b.now())
    b.write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'cell_id':cell_id,'seconds':elapsed}));return 0 if result['status']=='complete' else 1


def collect(path,inventory_path,out):
    tick=time.monotonic();cpu=time.process_time();p=load(path);inventory=read(inventory_path);b.check(set(inventory)==set(ALL),'Exact eight-slot inventory required')
    rows={};slots={cid:{'status':'unavailable'} for cid in ALL};costs=[]
    for cid,filename in inventory.items():
        if filename is None or not Path(filename).is_file():continue
        try:
            raw=read(filename)
            if raw.get('new_native_calls')==1:costs.append({'cell_id':cid,'job_id':raw.get('job_id'),'elapsed_seconds':raw.get('elapsed_seconds'),'allocated_cpu_seconds':raw.get('allocated_cpu_seconds'),'reported_status':raw.get('status')})
            admitted=admit(p,filename);b.check(admitted['cell_id']==cid,'Result inventory mismatch');rows[cid]=admitted;slots[cid]={'status':'complete','result':b.pin(filename)}
        except Exception as exc:slots[cid]={'status':'unavailable','result':b.pin(filename),'error':repr(exc)}
    scored=score(rows);qm=read(b.verify(p['comparison_authorities']['reference_derivative_result']));model=read(b.verify(p['comparison_authorities']['model_derivative_result']))
    b.check(qm['admitted_cells']==16 and qm['stage_numerical_checks_passed'] is True and model['complete_cells']==16 and model['error'] is None,'Saved comparison unavailable')
    comparison={}
    for m in (*METALS,'paired'):
        comparison[m]={}
        for q in ('0.01','0.005'):
            cc=None if scored['groups'][m] is None else scored['groups'][m]['steps'][q]['CCSD_T']['derivative_kcal_mol_A']
            ref=qm['scoring']['directions']['rotation'][m]['derivatives_kcal_mol_A'][q];mm=model['comparison']['directions']['rotation'][m]['steps'][q]['derivatives_kcal_mol_A']
            comparison[m][q]={'canonical_CCSD_T':cc,'PBE0_electronic':ref['electronic'],'PBE0_D4':ref['PBE0_D4'],'D4_only':ref['D4'],'fixed_model':mm['Total'],'fixed_model_components':mm,'model_minus_CC':None if cc is None else mm['Total']-cc,'PBE0_electronic_minus_CC':None if cc is None else ref['electronic']-cc,'PBE0_D4_minus_CC':None if cc is None else ref['PBE0_D4']-cc,'units':'kcal/mol/Angstrom','cross_method_accuracy_gate':None}
    status='canonical_rotation_complete' if len(rows)==8 else 'canonical_rotation_incomplete'
    result={'status':status,'implementation':b.pin(__file__),'panel':b.pin(path),'inventory':b.pin(inventory_path),'slots':slots,'admitted_cells':len(rows),'declared_cells':8,'scoring':scored,'comparison':comparison,'comparison_authorities':p['comparison_authorities'],'native_admissions':rows,'native_costs':costs,'native_allocated_CPU_seconds':sum(x['allocated_cpu_seconds'] for x in costs if isinstance(x['allocated_cpu_seconds'],(int,float)) and math.isfinite(x['allocated_cpu_seconds'])),'new_native_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False,'CP_force_qualified':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    out.mkdir(parents=True,exist_ok=False);rp=b.write(out/'RESULT.json',result)
    lines=['# Consumed dry-acetate canonical rotation','',f'{len(rows)}/8 admitted; all three two-step numerical checks pass: {scored["all_checks_passed"]}. No new native calls in collection.','','| Comparison | q (Å) | CCSD(T) | PBE0 electronic | PBE0-D4 | Fixed model |','|---|---:|---:|---:|---:|---:|']
    for m,steps in comparison.items():
        for q,v in steps.items():lines.append('| '+m+' | '+q+' | '+('unavailable' if v['canonical_CCSD_T'] is None else f'{v["canonical_CCSD_T"]:.9f}')+f' | {v["PBE0_electronic"]:.9f} | {v["PBE0_D4"]:.9f} | {v["fixed_model"]:.9f} |')
    lines+=['','Slopes are kcal/mol/Å, q=1Å×rotation angle. The0.1 resolution applies only to CC two-step consistency. Cross-method differences are descriptive. D4 is separate and never added to CC. Frozen-f/state/basis/diffuse/single-reference limitations and raw SCF/CC signals remain in RESULT; this is not a CP-corrected force or physical-truth claim.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':rp,'status':status,'numerical_checks_passed':scored['all_checks_passed']},indent=2))

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--panel',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--release',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.panel,a.out)
    elif a.action=='run':return run(a.panel,a.cell,a.release,a.out)
    else:collect(a.panel,a.inventory,a.out)
if __name__=='__main__':sys.exit(main())
