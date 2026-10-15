#!/usr/bin/env python3
"""Twelve fixed MACE points using the qualified warm worker; no pooled choice."""
import argparse,copy,fcntl,hashlib,importlib.util,json,math,os,shutil,signal,sys,time,traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/donor_H_sensitivity/scoring'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
METALS=('Ca','La');BASES=('origin','amide_flip','relaxed_flip_La');EV=23.06054783061903
SOURCES=tuple('a0acd6b9f2-pqq-la_model__conditioned_'+z+'__seed-1_sample-4' for z in METALS)
PROTOCOL='A0AC_two_H_radial_fixed_MACE_same_path_response_v1'
def require(v,why):
    if not v:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def tid(s,z,g):return s+'__'+z+'__'+g+'__two_H_radial'

def bind(preparation,output):
    p=read(preparation);pc=read(p['config']);oldpin=pin(WORK/'asn_rotamer_scoring_v1/CONFIG.json');old=read(oldpin);template=read(old['worker_manifest'])
    require(p['status']=='complete' and p['declared_endpoints']==12 and [r['source_case_id'] for r in p['rows']]==list(SOURCES),'Actual complete two-source preparation required')
    cells=[]
    for r in p['rows']:
        require(r['status']=='prepared','Source preparation unavailable')
        for z in METALS:
            for g in BASES:
                v=r['matrix'][z][g];require(v['status']=='prepared' and read(v['geometry_admission'])['pass'],'Endpoint geometry unavailable')
                cells.append({'task_id':tid(r['source_case_id'],z,g),'case_id':r['source_case_id'],'source_case_id':r['source_case_id'],'metal':z,'candidate':g,
                    **{k:v[k] for k in ('xyz','charge','multiplicity','atom_count','geometry_admission')},'state_signature':r['state_signature'],'maps':r['maps'],'source_preparation':r['source_preparation']})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for n in ('SCORE.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/n,out/n)
    worker=copy.deepcopy(template);worker.update(protocol_id=PROTOCOL,tasks=cells,cases=[],fixtures=[],origin_task_count=12,preparation=pin(preparation),
        finite_scope={'MACE_calls_maximum':12,'qualification_calls':0,'optimizer_starts':0,'SCC_calls':0,'DFT_calls':0})
    write(out/'worker_manifest.json',worker)
    c={'protocol':PROTOCOL,'implementation':pin(out/'SCORE.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(out/'PLAN.md'),'preparation':pin(preparation),
       'preparation_config':p['config'],'preparation_helper':pc['implementation']['PREPARE.py'],'rows':p['rows'],'cells':cells,'worker_manifest':pin(out/'worker_manifest.json'),
       'old_config':oldpin,'old_scorer':old['implementation'],'template':old['worker_manifest'],'python':old['python'],'placement_helper':old['prior_adapter'],
       'declared_new_cells':12,'reused_old_cells':12,'source_ids':list(SOURCES),'candidate_order':list(BASES),'resolution_kcal_mol':1.0,'kcal_per_eV':EV,'CPUs':7,'output':str(out/'execution')}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,b,scorer,np):
    require(c['protocol']==PROTOCOL and c['source_ids']==list(SOURCES) and c['candidate_order']==list(BASES),'Finite protocol differs')
    require(c['declared_new_cells']==c['reused_old_cells']==12 and c['resolution_kcal_mol']==1.0 and c['kcal_per_eV']==EV,'Units/resolution/scope differs')
    for k in ('plan','wrapper','preparation_helper'):check(c[k])
    require(pin(sys.executable)==c['python'],'Qualified parent interpreter differs')
    p=read(c['preparation']);pc=read(c['preparation_config']);m=read(c['worker_manifest']);template=read(c['template']);old=read(c['old_config'])
    require(p['status']=='complete' and p['rows']==c['rows'] and p['config']==c['preparation_config'] and pc['implementation']['PREPARE.py']==c['preparation_helper'],'Actual preparation join differs')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        require(m[k]==template[k],'Qualified worker changed: '+k)
    require(m['tasks']==c['cells'] and len(c['cells'])==12 and not m['fixtures'] and not m['cases'] and m['threads_per_model']==7 and m['lanes']==1,'Worker call/resource scope differs')
    pf=scorer.preflight(old,b);baseline={};checks=[]
    for row in c['rows']:
        sid=row['source_case_id'];saved,=[s for s in pc['sources'] if s['source_case_id']==sid];orig=saved['prepared']['original_row']
        require(row['status']=='prepared' and all(row[k]==orig[k] for k in ('maps','state_signature','source_preparation','protonated_source')),'Prepared source provenance changed')
        state=read(row['state_signature']);cp=read(row['source_preparation']);entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']}
        named={x['source']['atom']:i for i,x in entries.items() if x.get('source') and (x['source']['chain'],x['source']['resnum'],x['source']['resname'])==('A',290,'ASN')}
        moved=[named[h] for h in ('HD21','HD22')];require(row['moved_indices_zero_based']==moved,'Wrong changed physical H')
        targets=read(row['coordinate_derivation']['targets']);require(targets['forcefield']==pc['forcefield'] and row['coordinate_derivation']['kernel']==pc['radial_kernel']['mace_hydrogen.py'],'Projection parameter/kernel differs')
        baseline[sid]={z:{} for z in METALS}
        for g in BASES:
            pair={}
            for z in METALS:
                t,=[t for t in c['cells'] if t['task_id']==tid(sid,z,g)];v=row['matrix'][z][g];o=saved['baseline'][z][g]
                require(v['original']==o and v['original_xyz']==o['xyz'] and v['original_MACE']==o['receipt'] and t['xyz']==v['xyz'] and t['geometry_admission']==v['geometry_admission'],'Original/new point join differs')
                gate=read(v['geometry_admission']);require(gate['pass'] and gate['status']=='prepared' and gate['candidate']==g and gate['source_case_id']==sid,'Actual geometry admission missing')
                require((t['charge'],t['multiplicity'],t['atom_count'])==(state['endpoint_charges'][z],state['multiplicities'][z],204),'Endpoint chemical state differs')
                result=read(o['receipt']);request=read(result['request']) if result.get('request') else next(x for x in read(result['manifest'])['tasks'] if x['task_id']==result['task_id'])
                model=result['model'] if result.get('request') else read(result['manifest'])['model'];ot={**t,'xyz':o['xyz']}
                scorer.physical_result(result,request,ot,m['model'],model)
                if result['device']=='cpu':scorer.admitted_result(result,request,ot,m['model'],pf['versions'])
                else:require(result['device']=='NVIDIA H200','Unqualified archived device')
                require(o['status']=='complete' and result['energy_eV']==o['energy_eV'] and result['forces']==o['forces'],'Archived energy/forces changed')
                force=np.load(check(o['forces']),allow_pickle=False);require(force.shape==(204,3) and np.isfinite(force).all(),'Archived force array invalid')
                before=b.xyz(check(o['xyz']));after=b.xyz(check(t['xyz']));before_lines=check(o['xyz']).read_text().splitlines(keepends=True);after_lines=check(t['xyz']).read_text().splitlines(keepends=True)
                require(len(after)==len(before)==204 and [a[0] for a in after]==[a[0] for a in before] and after[0][0]==z,'Endpoint atom identity differs')
                require(all(before_lines[i+2]==after_lines[i+2] for i in range(204) if i not in moved),'Nonselected coordinate line changed')
                for h,term in zip(('HD21','HD22'),targets['bonds']):
                    require(term['atom_names']==['ND2',h],'Named target order differs');i=named[h];j=named['ND2'];anchor=np.asarray(before[j][1:]);vector=np.asarray(before[i][1:])-anchor
                    expected=anchor+vector*(term['forcefield_equilibrium_A']/np.linalg.norm(vector))
                    require(np.max(np.abs(np.asarray(after[i][1:])-expected))<=pc['settings']['replay_tolerance_A'],'Actual projected H does not replay')
                require(all(math.isfinite(x) for atom in after for x in atom[1:]),'Nonfinite projected coordinates');pair[z]=after
                baseline[sid][z][g]=o
            require(pair['Ca'][1:]==pair['La'][1:] and pair['Ca'][0][1:]==pair['La'][0][1:],'Paired metal geometry differs')
        checks.append({'source_case_id':sid,'six_original_receipts_verified':True,'six_new_geometries_verified':True,'moved_indices_zero_based':moved})
    return {'status':'passed','baseline':baseline,'versions':pf['versions'],'checks':checks,'model':m['model'],'prior_archive_receipts_rechecked':pf['validated_real204_atom_receipts'],
        'reused_experiment_cells':12,'new_model_calls':0,'precision':'Actual unchanged float64 model; raw serialized energies retained'}

def collect(c,pf,results):
    rows=[];complete=0
    for row in c['rows']:
        sid=row['source_case_id'];old={z:{g:row['matrix'][z][g]['original'] for g in BASES} for z in METALS};matrix={z:{} for z in METALS};paths={}
        for z in METALS:
            for g in BASES:
                new=results.get(tid(sid,z,g),{'status':'unavailable','energy_eV':None,'reason':'Not completed','xyz':row['matrix'][z][g]['xyz']})
                ok=new['status']=='complete';complete+=int(ok)
                matrix[z][g]={'original':old[z][g],'projected_H':new,'H_correction_kcal_mol':(new['energy_eV']-old[z][g]['energy_eV'])*EV if ok else None}
        for g in BASES[1:]:
            values={}
            for label in ('original','projected_H'):
                cells=[matrix[z][q][label] for z in METALS for q in ('origin',g)]
                if all(v['status']=='complete' for v in cells):
                    w={z:(matrix[z][g][label]['energy_eV']-matrix[z]['origin'][label]['energy_eV'])*EV for z in METALS}
                    values[label]={'per_metal_work_kcal_mol':w,'Ca_minus_La_differential_work_kcal_mol':w['Ca']-w['La']}
                else:values[label]=None
            change=values['projected_H']['Ca_minus_La_differential_work_kcal_mol']-values['original']['Ca_minus_La_differential_work_kcal_mol'] if all(values.values()) else None
            paths[g]={**values,'differential_work_change_kcal_mol':change,'absolute_change_at_least_resolution':abs(change)>=c['resolution_kcal_mol'] if change is not None else None}
        rows.append({'source_case_id':sid,'conditioning':row['conditioning'],'sample':4,'status':'complete' if all(matrix[z][g]['projected_H']['status']=='complete' for z in METALS for g in BASES) else 'unavailable',
                     'matrix':matrix,'same_path_responses':paths})
    return {'protocol':PROTOCOL,'preparation':c['preparation'],'declared_sources':2,'declared_new_cells':12,'complete_new_cells':complete,'reused_original_cells':12,'rows':rows,
        'status':'complete' if complete==12 else 'partial','kcal_per_eV':EV,'resolution_kcal_mol':c['resolution_kcal_mol'],'model':read(c['worker_manifest'])['model'],
        'baseline_runtime_admitted':pf is not None,'new_optimizer_calls':0,'new_SCC_calls':0,'new_DFT_calls':0,
        'interpretation':'Same-path vacuum response to two N-H lengths; no pool, classifier, affinity or force-field energy claim'}

def run(path,sha,action):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(pin(__file__)==c['implementation'],'Wrong frozen scorer')
    placement=load('_H_placement',c['placement_helper']).placement();scorer=load('_H_old_scorer',c['old_scorer']);m=read(c['worker_manifest']);b=load('_H_worker',m['implementation'])
    out=Path(c['output'] if action=='run' else str(Path(c['output']).parent/('preflight_'+os.environ['SLURM_JOB_ID'])));start=time.monotonic();cpu=time.process_time();transport=None;handlers={};results={};pf=None
    receipt={'status':'failed','config':cp,'action':action,'placement':placement,'new_optimizer_calls':0,'new_SCC_calls':0,'new_DFT_calls':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            b.parent_env();b.modules(m)
            import numpy as np
            pf=preflight(c,b,scorer,np);write(out/'PREFLIGHT.json',pf)
            if action=='run':
                attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers();transport=b.WarmCPU(check(c['worker_manifest']),m,attempt,0)
                require(transport.ready['threads']==7 and transport.ready['versions']==pf['versions'],'Actual warm runtime differs')
                for t in c['cells']:
                    req={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
                    v,rp=b.cell(attempt/'cells'/t['task_id'],req,transport,m['model']);r={'status':'unavailable','receipt':rp,'xyz':t['xyz'],'energy_eV':None,'reason':v.get('reason')}
                    try:
                        scorer.admitted_result(v,read(v['request']),t,m['model'],pf['versions']);r.update(status='complete',energy_eV=v['energy_eV'],forces=v['forces'],request=v['request'])
                    except Exception as exc:r['reason']=str(exc)
                    results[t['task_id']]=r
                transport.close();transport=None
            receipt['status']='complete' if action=='preflight' or len(results)==12 and all(v['status']=='complete' for v in results.values()) else 'partial'
        except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc());raise
        finally:
            try:write(out/'PARTIAL_RESULTS.json',{'config':cp,'declared_cells':12,'results':results})
            except BaseException as exc:receipt.update(status='failed',partial_reporting_error=str(exc))
            try:
                if transport:transport.abort()
                for sig,handler in handlers.items():signal.signal(sig,handler)
            except BaseException as exc:receipt.update(status='failed',cleanup_error=str(exc))
            try:
                if action=='run':write(out/'FINAL_COLLECTION.json',{'config':cp,**collect(c,pf,results)})
            except BaseException as exc:receipt.update(status='failed',collection_error=str(exc),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,parent_CPU_seconds=time.process_time()-cpu)
                try:receipt.update(actual_model_forwards=sum(read(p).get('actual_model_forwards',0) for p in out.glob('mace/cells/*/mace_result.json')),worker_receipts=[pin(p) for p in out.glob('mace/lane*/CPU_WORKER_SUMMARY.json')])
                except BaseException as exc:receipt.update(status='failed',cost_reporting_error=str(exc))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);b=s.add_parser('bind');b.add_argument('--preparation',required=True);b.add_argument('--output',required=True)
    for op in ('preflight','run'):
        q=s.add_parser(op);q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='bind':bind(a.preparation,a.output)
    else:sys.exit(run(a.config,a.config_sha,a.action))
