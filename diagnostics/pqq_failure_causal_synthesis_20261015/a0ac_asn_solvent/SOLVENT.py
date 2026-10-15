#!/usr/bin/env python3
"""Conditional Asn flip: existing fixed-charge native ALPB, full common3→4."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_solvent'
TEMPLATE=ROOT/'workspaces/exaa_metal_midpoint_20261014/scoring_v1/CONFIG.json'
BASES=('origin','adaptive_Ca','adaptive_La');NEW='amide_flip';ORDER=BASES+(NEW,);METALS=('Ca','La');POLICIES=('operational','mathematical','origin')
PROTOCOL='A0ACD6B9F2_Asn290_fixed_180degree_flip_frozen_ALPB_common4_v1'
KEYS=('matrix','runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model','corrected_ORCA_parameter_provenance','ORCA_correction_policy')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def require(v,reason):
    if not v:raise ValueError(reason)

def gate(sc,sr,ex,sp):
    require(sr['config']==sp==ex['config'] and sr['terminal'],'Exact terminal MACE artifact required')
    require(len(sc['cells'])==sc['admitted_new_cells']==sr['complete_new_cells']==20 and len(sr['rows'])==10,'Actual complete ten-source/twenty-cell MACE matrix required')
    require(ex['status']=='complete' and ex['actual_model_forwards']==20,'MACE execution incomplete')
    require(sc['new_order']==list(ORDER) and sr['stage_gate']['decision']=='qualifies_for_root_solvent_review' and sr['stage_gate']['threshold_kcal_mol']==1.0,'Recorded declared solvent gate not passed')

def prepare(score_config,score_collection,output):
    """Bind actual completed artifacts only; scientific gate replay stays in Slurm."""
    sc=read(score_config);sr=read(score_collection);old=read(TEMPLATE);inv=read(sc['inventory'])
    execution=Path(score_collection).with_name('EXECUTION.json');ex=read(execution);gate(sc,sr,ex,pin(score_config))
    cells=[];baselines=[]
    for source in inv['rows']:
        sid=source['source_case_id'];row,=[r for r in sr['rows'] if r['source_case_id']==sid]
        for z in METALS:
            for q in BASES:
                v=source['matrix'][z][q];native=read(v['native_frozen_ALPB']);t=native['input']
                original={**v,**{k:t[k] for k in ('charge','multiplicity','atom_count','valence_electrons')},
                    'native_frozen_ALPB_result':v['native_frozen_ALPB'],'native_frozen_ALPB_hartree':native['G_qvac_hartree']}
                baselines.append({'source_case_id':sid,'metal':z,'candidate':q,'original':original,
                    **{k:source[k] for k in ('maps','state_signature','source_preparation')}})
            task,=[t for t in sc['cells'] if t['source_case_id']==sid and t['metal']==z]
            value=row['matrix'][z][NEW];require(value['status']=='complete' and value['xyz']==task['xyz'],'Completed flipped MACE cell missing')
            origin=read(source['matrix'][z]['origin']['native_frozen_ALPB'])['input']
            verify(value['receipt']);verify(task['xyz'])
            cells.append({**task,'valence_electrons':origin['valence_electrons'],'MACE':value['receipt'],'MACE_eV':value['energy_eV']})
    review=read(inv['review_input']);historical=read(review['completed_scores'])
    require(all(len([r for r in historical['rows'] if r['record_id']==s['record_id'] and r==s['current_score']])==1 for s in inv['rows']),'Current original pools differ from authoritative collection')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('SOLVENT.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    c={k:old[k] for k in KEYS}
    c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),wrapper=pin(out/'RUN.sbatch'),
        plan=pin(HERE.parent/'ASN_ROTAMER_PLAN.md'),adapter_plan=pin(HERE/'PLAN.md'),template=pin(TEMPLATE),
        score_config=pin(score_config),score_collection=pin(score_collection),score_execution=pin(execution),score_adapter=sc['implementation'],
        preparation=sc['preparation'],preparation_manifest=sc['preparation_manifest'],inventory=sc['inventory'],rows=sc['rows'],cells=cells,baselines=baselines,
        historical_collection=review['completed_scores'],solvent_adapter=old['solvent_adapter'],pool_adapter=old['pool_adapter'],scorer=old['scorer'],
        numerical_python=sc['python'],CPUs=7,protein='A0ACD6B9F2',declared_sources=10,declared_new_cells=20,reused_original_cells=60,
        original_order=list(BASES),new_order=list(ORDER),maximum_SCC_calls=80,maximum_native_component_calls=40,
        new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0,classification=None,
        parameter_role='Published standalone GFN2 parameter file and qualified internal native ALPB; corrected full ORCA JSON is audited provenance only, never default ORCA GFN2')
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,f,r,n,a,h,selector,score):
    for k in ('plan','adapter_plan','wrapper','corrected_ORCA_parameter_provenance','ORCA_correction_policy'):verify(c[k])
    require(pin(sys.executable)==c['numerical_python'],'Qualified numerical interpreter differs')
    old=read(c['template'])
    for k in KEYS:require(c[k]==old[k],'Qualified solvent kernel/settings changed: '+k)
    require(read(c['ORCA_correction_policy'])['candidate_sha256']==c['corrected_ORCA_parameter_provenance']['sha256'],'Corrected full La parameter lineage differs')
    require((c['protocol_id'],c['protein'],c['declared_sources'],c['declared_new_cells'],c['reused_original_cells'])==(PROTOCOL,'A0ACD6B9F2',10,20,60),'Finite scope differs')
    require(c['original_order']==list(BASES) and c['new_order']==list(ORDER),'Pool order differs')
    sc=read(c['score_config']);sr=read(c['score_collection']);ex=read(c['score_execution']);prep=read(c['preparation']);gate(sc,sr,ex,c['score_config'])
    require(sr['preparation']==c['preparation']==sc['preparation'] and prep['rows']==c['rows']==sc['rows'],'Preparation/MACE joins differ')
    require(sc['implementation']==c['score_adapter'] and sc['plan']==c['plan'] and sc['inventory']==c['inventory'] and sc['preparation_manifest']==c['preparation_manifest'],'MACE protocol provenance differs')
    worker=read(sc['worker_manifest']);b=load('_asn_saved_worker_metadata',worker['implementation']);pf=score.preflight(sc,b)
    require(worker['model']==c['MACE_model'],'Native/MACE model lineage differs')
    mr={};saved_results={};require(len(c['cells'])==20 and len(c['baselines'])==60,'Finite endpoint count differs')
    for cell in c['cells']:
        t,=[t for t in sc['cells'] if t['task_id']==cell['task_id']]
        require(all(cell[k]==t[k] for k in t),'Completed MACE input changed')
        row,=[row for row in sr['rows'] if row['source_case_id']==t['source_case_id']];v=row['matrix'][t['metal']][NEW]
        require(v['status']=='complete' and v['receipt']==cell['MACE'] and v['energy_eV']==cell['MACE_eV'],'Completed MACE receipt differs')
        result=read(cell['MACE']);score.admitted_result(result,read(result['request']),t,c['MACE_model'],pf['versions'])
        saved_results[t['task_id']]=v;mr[t['task_id']]={'status':'complete','receipt':cell['MACE'],'energy_eV':cell['MACE_eV']}
        symbols,_=f.coordinates(cell['xyz']);require(len(symbols)==204 and sum(f.VALENCE[z] for z in symbols)-cell['charge']==cell['valence_electrons'],'Native204-atom electron inventory differs')
    replay=score.collect(sc,pf,saved_results)
    require(replay['rows']==sr['rows'] and replay['stage_gate']==sr['stage_gate'],'Actual MACE common4/gate replay differs')
    require(all(row['common3_to4']['mathematical'] is not None for row in replay['rows']) and any(abs(row['common3_to4']['mathematical']['paired_contrast_change_kcal_mol'])>=1.0 for row in replay['rows']),'Declared complete sign-blind stage gate failed')
    runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for k in ('xtb','xtb_parameter_file','tblite_python'):verify(runtime[k])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:verify(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:verify(p)
    for p in list(nm['implementation'].values())+nm['build']['pins']+build['linked_libraries']:verify(p)
    verify(build['binary']);require(build['manifest']==c['native_manifest'] and build['binary']==c['native_binary'] and nm['MACE_model']==c['MACE_model'],'Qualified native build/model differs')
    inv=read(c['inventory']);historical=read(c['historical_collection']);old_rows={};count=0;residuals=[]
    for source in inv['rows']:
        sid=source['source_case_id'];stored,=[row for row in historical['rows'] if row['record_id']==source['record_id']]
        require(stored==source['current_score'],'Authoritative baseline pool differs');evidence={z:{} for z in METALS}
        for cell in (cell for cell in c['baselines'] if cell['source_case_id']==sid):
            z,q=cell['metal'],cell['candidate'];original=source['matrix'][z][q]
            require(all(cell['original'][k]==v for k,v in original.items()),'Original cell metadata differs')
            for k in ('maps','state_signature','source_preparation'):require(cell[k]==source[k],'Original source state differs')
            v=a.reusable(c,cell,n,r)
            require(v['MACE']==pf['baseline'][sid][z][q]['receipt'] and v['MACE_eV']==pf['baseline'][sid][z][q]['energy_eV'],'Original MACE/native pairing differs')
            evidence[z][q]=v;count+=1
        require(all(set(evidence[z])==set(BASES) for z in METALS),'Incomplete original common3')
        matrix={z:{q:{'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} for q,v in evidence[z].items()} for z in METALS}
        h.reproduce(h.choose(selector,matrix,BASES),stored['pool'],residuals)
        old_rows[sid]={'evidence':evidence,'historical3':stored['pool']}
    require(count==60,'Exactly sixty existing solvent cells must be reused')
    return {'status':'passed','stage_gate':replay['stage_gate'],'MACE_collection':c['score_collection'],'reused_MACE':mr,'old_rows':old_rows,
        'reused_original_cells':60,'historical_pool_replay':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
        'new_molecular_calls':0,'qualification':'Unchanged component-v2, native water diagonal, cavity equality and vacuum charge precision; historical absolute v1 failures remain visible'}

def solvent(c,base,f,r,n,mr):
    """MATRIX.solvent dispatch copy: only protein metadata generalized to c['protein']."""
    directory=base/'solvent';directory.mkdir();tasks=[]
    for cell in c['cells']:
        for medium in ('vacuum','water'):
            tasks.append({'id':cell['task_id']+'__'+medium,'protein':c['protein'],'context':cell['case_id'],
                'geometry':cell['candidate'],'medium':medium,'uhf':0,'reference':None,'scout':False,'reuse':None,
                **{k:cell[k] for k in ('xyz','charge','multiplicity','valence_electrons','atom_count','metal')}})
    m={'protocol_id':c['protocol_id']+'_GFN2','implementation':c['factorial'],'runtime':c['runtime'],'settings':c['settings'],'tasks':tasks,'config':pin(base.parent/'CONFIG.json')}
    write(directory/'manifest.json',m);runtime=read(c['runtime'])
    with ThreadPoolExecutor(max_workers=7) as pool:points=list(pool.map(lambda t:f.execute_point(t,m,directory,runtime),tasks))
    native=read(c['native_manifest']);out=directory/'replay';out.mkdir();rows=[];env=f.clean_env(runtime['xtb_parameter_file']['path'])
    for cell in c['cells']:
        tid=cell['task_id'];record={'task_id':tid,'status':'unavailable','reason':None,'G_qvac_hartree':None}
        try:
            v,w=[next(t for t in tasks if t['id']==tid+'__'+medium) for medium in ('vacuum','water')]
            component=r.component(directory,v,w);record['component_parity']=component
            require(component['accepted'],component.get('reason','component qualification failed'));require(mr[tid]['status']=='complete','MACE unavailable')
            t={k:cell[k] for k in ('xyz','atom_count','metal','charge','multiplicity','valence_electrons','state_signature','source_preparation','maps')}
            t.update(task_id=tid,scientific_input_key={'runtime':c['runtime']},uses=[{'MACE':mr[tid]['receipt'],'MACE_eV':mr[tid]['energy_eV']}],endpoints={})
            for medium in ('vacuum','water'):
                p=directory/'points'/(tid+'__'+medium);point=read(p/'RESULT.json');ep=p/'xtb';execution=read(ep/'EXECUTION.json')
                t['endpoints'][medium]={'point':pin(p/'RESULT.json'),'point_manifest':point['manifest'],'task':point['task'],
                    'execution':pin(ep/'EXECUTION.json'),'stdout':execution['stdout'],'charges':pin(ep/'charges'),'executed_xyz':pin(ep/'core.xyz')}
            record.update(n.run_task(native,verify(c['native_manifest']),t,c['native_binary'],out,env,None))
        except Exception as e:record['reason']=type(e).__name__+': '+str(e)
        rows.append(record)
    write(directory/'COLLECTION.json',{'rows':rows,'raw_inputs':len(points),'maximum_backend_calls':2*len(tasks)});return rows

def summarize(c,pf,sr,a,h,n,r,selector,score):
    native={v['task_id']:v for v in sr};rows=[];residuals=[];available=0
    for source in c['rows']:
        sid=source['source_case_id'];old=pf['old_rows'][sid];evidence=copy.deepcopy(old['evidence']);fail=[]
        for cell in (cell for cell in c['cells'] if cell['source_case_id']==sid):
            z=cell['metal'];v=native.get(cell['task_id']);evidence[z][NEW]={'status':'unavailable','reason':'not executed'}
            try:
                require(v is not None and v['status']=='available','Native/component endpoint unavailable: '+str(v.get('reason') if v else 'missing'))
                original={**cell,'native_frozen_ALPB_result':v['result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
                evidence[z][NEW]=a.reusable(c,{'original':original,**{k:cell[k] for k in ('state_signature','maps','source_preparation')}},n,r);available+=1
            except Exception as e:evidence[z][NEW]={'status':'unavailable','reason':str(e)};fail.append({'task_id':cell['task_id'],'metal':z,'native_record':v,'reason':str(e)})
        matrix={z:{q:({'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} if v['status']=='available' else None) for q,v in evidence[z].items()} for z in METALS}
        pools={name:h.choose(selector,matrix,order) for name,order in (('original3',BASES),('expanded4',ORDER))};h.reproduce(pools['original3'],old['historical3'],residuals)
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for policy in POLICIES:
            before,after=pools['original3'][policy],pools['expanded4'][policy]
            if before is None or after is None:changes[policy]=None;continue
            metals={z:a.difference(works[z][after[z]['candidate']],works[z][before[z]['candidate']]) for z in METALS}
            delta=a.difference(metals['Ca'],metals['La']);raw=after['R_model_kcal_mol']-before['R_model_kcal_mol'];h.close(raw,delta['total'],residuals)
            changes[policy]={'per_metal_change_kcal_mol':metals,'differential_change_kcal_mol':delta,'R_model_change_kcal_mol':raw,
                'before':{z:before[z]['candidate'] for z in METALS},'after':{z:after[z]['candidate'] for z in METALS}}
        flip={z:works[z][NEW] for z in METALS}
        rows.append({'source_case_id':sid,'sample':source['sample'],'conditioning':source['conditioning'],'pools':pools,'common3_to4':changes,
            'all_origin_relative_works_kcal_mol':works,'flip_work_kcal_mol':flip,'Ca_minus_La_flip_work_kcal_mol':a.difference(flip['Ca'],flip['La']),
            'evidence':evidence,'failed_new_cells':fail,'classification':None})
    grouped={}
    for condition in METALS:
        group=[row for row in rows if row['conditioning']==condition];matched=[row for row in group if row['pools']['expanded4']['status']=='available']
        grouped[condition]={'declared_sources':5,'complete_expanded_sources':len(matched),'matched_source_ids':[row['source_case_id'] for row in matched],
            'by_policy':{p:{'original3_all5':score.span([row['pools']['original3'][p]['R_model_kcal_mol'] for row in group]),
                'original3_matched':score.span([row['pools']['original3'][p]['R_model_kcal_mol'] for row in matched]),
                'expanded4_matched':score.span([row['pools']['expanded4'][p]['R_model_kcal_mol'] for row in matched])} for p in POLICIES}}
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_sources':10,'declared_new_cells':20,'available_new_cells':available,'reused_original_cells':60,
        'rows':rows,'conditioning_spreads':grouped,'MACE_stage_gate':pf['stage_gate'],'pool_orders':{'original3':list(BASES),'expanded4':list(ORDER)},
        'saved_algebra_checks':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
        'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'classification':None,'new_calibration':False,'production_changed':False,
        'interpretation':'Finite fixed flip with represented short H contacts and no relaxation; not all Asn conformations, physical equilibrium, affinity or classifier validation'}

def execute(path,action):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required before scientific imports')
    c=read(path);require(c['implementation']==pin(__file__),'Wrong adapter snapshot');matrix=load('_asn_solvent_placement',c['matrix']);placement=matrix.placement()
    out=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID']);start=time.monotonic();pf=None;sr=[]
    receipt={'status':'failed','config':pin(path),'placement':placement,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            f=load('_asn_factorial',c['factorial']);r=load('_asn_component',c['component']);n=load('_asn_native',c['native_adapter'])
            a=load('_asn_prior_solvent_checks',c['solvent_adapter']);h=load('_asn_pool_checks',c['pool_adapter']);selector=load('_asn_selector',c['scorer']);score=load('_asn_vacuum_gate',c['score_adapter'])
            pf=preflight(c,f,r,n,a,h,selector,score);write(out/'PREFLIGHT.json',pf)
            if action=='execute':sr=solvent(c,out,f,r,n,pf['reused_MACE'])
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            try:
                if pf and action=='execute':
                    if not sr:
                        for cell in c['cells']:
                            p=out/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                            if p.exists():
                                v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    write(out/'PARTIAL_RESULTS.json',{'native_rows':sr,'config':pin(path),'declared_new_cells':20})
                    result={'config':pin(path),**summarize(c,pf,sr,a,h,n,r,selector,score)};write(out/'FINAL_COLLECTION.json',result)
                    if result['available_new_cells']!=20:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
                try:
                    receipt.update(SCC_started=sum(1 for _ in out.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in out.glob('solvent/points/*/*/EXECUTION.json')),
                        native_component_calls=sum(1 for _ in out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
                        measured_resource_receipts={'SCC':[pin(p) for p in sorted(out.glob('solvent/points/*/*/EXECUTION.json'))],
                            'native_ALPB':[pin(p) for p in sorted(out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]})
                except BaseException as e:receipt.update(status='failed',resource_reporting_failure={'reason':str(e),'traceback':traceback.format_exc()})
                write(out/'EXECUTION.json',receipt)
    if receipt['status']=='failed':raise RuntimeError('Execution or reporting failed; inspect EXECUTION.json')

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--score-config',required=True);q.add_argument('--score-collection',required=True);q.add_argument('--output',required=True)
    for name in ('preflight','execute'):
        q=sub.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.score_config,a.score_collection,a.output)
    else:execute(verify({'path':str(Path(a.config).resolve()),'sha256':a.config_sha}),a.action)
