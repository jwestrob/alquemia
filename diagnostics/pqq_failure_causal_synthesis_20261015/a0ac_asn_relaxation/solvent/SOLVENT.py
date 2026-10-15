#!/usr/bin/env python3
"""Conditional common6 adapter; unchanged executed Asn solvent kernels."""
import argparse, copy, fcntl, hashlib, importlib.util, json, os, shutil, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/solvent'
OLD=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1'
METALS=('Ca','La');ORDER3=('origin','adaptive_Ca','adaptive_La');ORDER4=ORDER3+('amide_flip',);NEW=('relaxed_flip_Ca','relaxed_flip_La');ORDER6=ORDER4+NEW
PROTOCOL='A0AC_relaxed_amide_frozen_distribution_common6_v1'
def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def require(v,why):
    if not v:raise ValueError(why)
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def gate(rc,sr,ex,cp):
    require(sr['config']==ex['config']==cp and rc['stage']==sr['stage']=='remaining','Actual remaining-stage provenance required')
    require(ex['status']=='complete' and sr['stage_complete'] and sr['all20_searches_and_common6_complete'],'All twenty credible searches/common6 required')
    require(sr['declared_sources']==10 and sr['declared_new_searches']==20 and sr['reused_endpoint_cells']==80,'Declared population differs')
    require(len(sr['rows'])==10 and [r['source_case_id'] for r in sr['rows']]==rc['all_source_ids'],'All ten ordered sources required')
    require(sr['pool_orders']==rc['pool_orders']=={'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)},'Common-pool identity differs')
    g=sr['solvent_review'];require(g['eligible'] and g['threshold_kcal_mol']==rc['solvent_review_threshold_kcal_mol']==1.0 and g['source_ids_above_sign_blind_gate'],'Declared sign-blind review gate failed')
    for row in sr['rows']:
        require(row['status']=='complete' and set(row['proposals'])==set(METALS),'Incomplete proposal pair')
        for z in METALS:
            p=read(row['proposals'][z]);require(p['status']=='proposal_available' and p['optimizer']['success'] and p['final_geometry']['physical_feasible'],'Noncredible proposal')
            for q in ORDER6:require(row['matrix'][z][q]['status']=='complete','Incomplete common6 MACE matrix')

def bind(score_config,score_collection,output):
    """Metadata only. Gate numerics and physical replay occur later in Slurm."""
    rc=read(score_config);sr=read(score_collection);expath=Path(score_collection).with_name('EXECUTION.json');gate(rc,sr,read(expath),pin(score_config))
    prior=read(OLD/'CONFIG.json');saved=read(OLD/'execution/FINAL_COLLECTION.json');prep=read(rc['preparation']);inv=read(prior['inventory'])
    scout=read(rc['scout_collection']);sexpath=Path(rc['scout_collection']['path']).with_name('EXECUTION.json');sex=read(sexpath)
    require(sex['config']==scout['config'] and sex['status']=='complete' and scout['stage_complete'],'Actual scout execution unavailable')
    require(saved['config']==pin(OLD/'CONFIG.json') and saved['available_new_cells']==20 and read(OLD/'execution/EXECUTION.json')['status']=='complete','Actual complete common4 solvent required')
    require(rc['old_config']==prior['score_config'] and rc['old_collection']==prior['score_collection'],'Relaxation and old solvent MACE lineage differ')
    require([r['source_case_id'] for r in prep['rows']]==[r['source_case_id'] for r in inv['rows']]==rc['all_source_ids'],'Prepared/old source order differs')
    cells=[]
    for source,row,p in zip(inv['rows'],sr['rows'],prep['rows']):
        require(p['status']=='prepared' and row['prepared_mapping']==p['mapping'] and row['admission']==p['admission'],'Prepared mapping/admission differs')
        for z in METALS:
            origin=read(source['matrix'][z]['origin']['native_frozen_ALPB'])['input']
            for q in NEW:
                v=row['matrix'][z][q];require(v['proposal_receipt']==row['proposals'][q.removeprefix('relaxed_flip_')],'Proposal identity differs')
                for k in ('xyz','receipt','full_source_admission'):check(v[k])
                require(read(v['full_source_admission'])['pass'],'Final source geometry unadmitted')
                cells.append({'task_id':row['source_case_id']+'__'+q+'__'+z,'source_case_id':row['source_case_id'],'case_id':p['original_row']['case_id'],
                    'candidate':q,'metal':z,'sample':row['sample'],'conditioning':row['conditioning'],
                    'xyz':v['xyz'],'MACE':v['receipt'],'MACE_eV':v['energy_eV'],
                    **{k:origin[k] for k in ('charge','multiplicity','atom_count','valence_electrons')},
                    'maps':{metal:p['mapping'] for metal in METALS},'state_signature':source['state_signature'],'source_preparation':source['source_preparation'],
                    'proposal_receipt':v['proposal_receipt'],'full_source_admission':v['full_source_admission'],'prepared_admission':p['admission']})
    require(len(cells)==40 and len({t['task_id'] for t in cells})==40,'Exact forty-cell scope required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('SOLVENT.py','RUN.sbatch'):shutil.copy2(HERE/name,out/name)
    deps=('matrix','runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model',
          'corrected_ORCA_parameter_provenance','ORCA_correction_policy','solvent_adapter','pool_adapter','scorer','numerical_python')
    c={k:prior[k] for k in deps};c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),wrapper=pin(out/'RUN.sbatch'),plan=pin(HERE/'PLAN.md'),
        relaxation_plan=rc['plan'],prior_config=pin(OLD/'CONFIG.json'),prior_collection=pin(OLD/'execution/FINAL_COLLECTION.json'),
        prior_execution=pin(OLD/'execution/EXECUTION.json'),prior_solvent=pin(OLD/'execution/solvent/COLLECTION.json'),prior_adapter=prior['implementation'],
        score_config=pin(score_config),score_collection=pin(score_collection),score_execution=pin(expath),score_adapter=rc['implementation'],
        scout_collection=rc['scout_collection'],scout_execution=pin(Path(rc['scout_collection']['path']).with_name('EXECUTION.json')),
        preparation=rc['preparation'],inventory=prior['inventory'],rows=[{k:r[k] for k in ('source_case_id','sample','conditioning','record_id','case_key')} for r in sr['rows']],cells=cells,
        protein='A0ACD6B9F2',CPUs=7,declared_sources=10,declared_new_cells=40,reused_original_cells=80,maximum_SCC_calls=160,maximum_native_component_calls=80,
        pool_orders=rc['pool_orders'],stage_gate=sr['solvent_review'],automatic_submission=False,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,old,f,r,n,a,h,selector,score,relax):
    for k in ('plan','relaxation_plan','wrapper','prior_execution','scout_execution'):check(c[k])
    prior=read(c['prior_config']);rc=read(c['score_config']);sr=read(c['score_collection']);gate(rc,sr,read(c['score_execution']),c['score_config'])
    require(c['protocol_id']==PROTOCOL and c['declared_sources']==10 and c['declared_new_cells']==len(c['cells'])==40 and c['reused_original_cells']==80,'Finite scope differs')
    require(c['preparation']==rc['preparation'] and c['inventory']==prior['inventory'] and c['score_adapter']==rc['implementation'] and c['prior_adapter']==prior['implementation'],'Dependency joins differ')
    require(c['stage_gate']==sr['solvent_review'] and c['pool_orders']==rc['pool_orders'],'Recorded gate/order differs')
    for k in old.KEYS+('solvent_adapter','pool_adapter','scorer','numerical_python'):require(c[k]==prior[k],'Qualified dependency differs: '+k)
    priorpf=old.preflight(prior,f,r,n,a,h,selector,score)
    replay=old.summarize(prior,priorpf,read(c['prior_solvent'])['rows'],a,h,n,r,selector,score);saved=read(c['prior_collection'])
    require(replay['rows']==saved['rows'] and replay['available_new_cells']==20,'Actual eighty-cell common4 replay differs')
    worker=read(rc['worker_manifest']);b=load('_relaxed_solvent_cpu_metadata',worker['implementation']);b.parent_env();b.modules(worker)
    import numpy as np
    from mace_site_kinematics import Kinematics
    _,completion=relax.preflight(rc,b,score,np)
    require(relax.collect(rc,Path(c['score_collection']['path']).parent,score)==sr,'Actual common6/gate numerical replay differs')
    scout=read(c['scout_collection']);sex=read(c['scout_execution']);require(sex['config']==scout['config'] and sex['status']=='complete' and scout['stage_complete'],'Actual scout execution unavailable')
    helper=load('_relaxed_solvent_fullsource',rc['source_preparation_helper']);prep=read(c['preparation']);mr={};checks=[]
    for row,pr,original in zip(sr['rows'],prep['rows'],replay['rows']):
        sid=row['source_case_id'];require(sid==pr['source_case_id']==original['source_case_id'],'Source order differs')
        require(row['prepared_mapping']==pr['mapping'] and row['admission']==pr['admission'] and read(pr['admission'])['pass'],'Actual preparation binding differs')
        for z in METALS:
            for q in ORDER4:
                v=row['matrix'][z][q];e=original['evidence'][z][q]
                require((v['receipt'],v['xyz'],v['energy_eV'])==(e['MACE'],e['xyz'],e['MACE_eV']),'Old common4/MACE cell differs')
        kin=Kinematics(read(pr['mapping'])['context']);contact=pr['contact'];ids=read(contact['identities'])['centers']
        with np.load(check(contact['origin']),allow_pickle=False) as data:union=data['positions'].copy();symbols=data['symbols'].tolist()
        for q in NEW:
            metal=q.removeprefix('relaxed_flip_');proposal=read(row['proposals'][metal]);point=proposal['proposal'];task,=[t for t in pr['tasks'] if t['metal']==metal]
            require(proposal['task_id']==task['task_id'] and proposal['case_id']==task['case_id'],'Proposal task/state differs')
            require(read(proposal['manifest'])['tasks']==pr['tasks'],'Proposal manifest source tasks differ')
            active=np.asarray(point['active_q_radian']);full=completion.angular.full_q(task,active)
            require(np.array_equal(full,np.asarray(point['full_q'])),'Final q mapping differs')
            coords=kin.evaluate(full)[1];geometry=completion.angular.final_geometry(kin,task,active,[x[0] for x in b.xyz(check(point['coordinate']))])
            require(geometry['physical_feasible'],'Final physical geometry failed replay')
            changed=union.copy();changed[:204]=coords;contact_result=helper.contact_check(np,union,changed,symbols,contact['moving_indices'],ids,set(contact['caps']))
            require(contact_result['pass'],'Final full-source contact replay failed')
            pair={}
            for z in METALS:
                cell,=[t for t in c['cells'] if t['source_case_id']==sid and t['candidate']==q and t['metal']==z];v=row['matrix'][z][q]
                require(cell['proposal_receipt']==row['proposals'][metal] and cell['full_source_admission']==v['full_source_admission'],'Final proposal/admission pin differs')
                sg=read(cell['full_source_admission']);require(sg['pass'] and sg['task_id']==task['task_id'] and sg['active_q_radian']==point['active_q_radian'] and sg['contact_inputs']==contact,'Final contact provenance differs')
                require(cell['maps']=={zz:pr['mapping'] for zz in METALS} and cell['prepared_admission']==pr['admission'],'Actual cap/source map differs')
                for k in ('state_signature','source_preparation'):require(cell[k]==pr['original_row'][k],'Original physical state differs')
                require(cell['MACE']==v['receipt'] and cell['xyz']==v['xyz'] and cell['MACE_eV']==v['energy_eV'],'Final cell pin differs')
                if z==metal:require(cell['MACE']==point['MACE'] and cell['xyz']==point['coordinate'],'Own-metal final differs')
                result=read(cell['MACE']);score.admitted_result(result,read(result['request']),cell,c['MACE_model'],prep['versions'])
                pair[z]=b.xyz(check(cell['xyz']));require(len(pair[z])==cell['atom_count']==204 and pair[z][0][0]==z,'Atom/metal inventory differs')
                require(np.max(np.abs(np.asarray([x[1:] for x in pair[z]])-coords))<=1e-12,'Actual scored frame differs from source/cap replay')
                state=read(cell['state_signature']);require(cell['charge']==task['charge']+(0 if z==metal else (1 if z=='La' else -1)) and cell['charge']==state['endpoint_charges'][z] and cell['multiplicity']==1,'Endpoint state differs')
                require(sum(f.VALENCE[x[0]] for x in pair[z])-cell['charge']==cell['valence_electrons'],'Valence inventory differs')
                mr[cell['task_id']]={'status':'complete','receipt':cell['MACE'],'energy_eV':cell['MACE_eV']}
            require(pair['Ca'][1:]==pair['La'][1:] and pair['Ca'][0][1:]==pair['La'][0][1:],'Paired physical coordinates differ')
            checks.append({'source_case_id':sid,'candidate':q,'proposal':row['proposals'][metal],'physical_pass':True,'full_source_pass':True})
    require(len(checks)==20 and len(mr)==40,'Exact forty final frames required')
    return {'status':'passed','MACE_stage_gate':sr['solvent_review'],'geometry_checks':checks,'reused_MACE':mr,'old_rows':replay['rows'],'reused_original_cells':80,'new_molecular_calls':0}

def summarize(c,pf,sr,a,h,n,r,selector):
    native={v['task_id']:v for v in sr};rows=[];residuals=[];available=0
    for source,old in zip(c['rows'],pf['old_rows']):
        sid=source['source_case_id'];require(sid==old['source_case_id'],'Summary source order differs');evidence=copy.deepcopy(old['evidence']);fail=[]
        for cell in (t for t in c['cells'] if t['source_case_id']==sid):
            z,q=cell['metal'],cell['candidate'];v=native.get(cell['task_id'])
            try:
                require(v is not None and v['status']=='available','Native/component endpoint unavailable: '+str(v.get('reason') if v else 'not executed'))
                original={**cell,'native_frozen_ALPB_result':v['result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
                evidence[z][q]=a.reusable(c,{'original':original,**{k:cell[k] for k in ('state_signature','maps','source_preparation')}},n,r);available+=1
            except Exception as exc:evidence[z][q]={'status':'unavailable','reason':str(exc)};fail.append({'task_id':cell['task_id'],'native_record':v,'reason':str(exc)})
        matrix={z:{q:({'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} if v['status']=='available' else None) for q,v in evidence[z].items()} for z in METALS}
        pools={label:h.choose(selector,matrix,order) for label,order in (('common3',ORDER3),('common4',ORDER4),('common6',ORDER6))}
        require(pools['common3']==old['pools']['original3'] and pools['common4']==old['pools']['expanded4'],'Old composite pools changed')
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for label,before,after in (('common3_to4','common3','common4'),('common4_to6','common4','common6'),('common3_to6','common3','common6')):
            changes[label]={}
            for policy in ('operational','mathematical','origin'):
                b,e=pools[before][policy],pools[after][policy];changes[label][policy]=None
                if b is None or e is None:continue
                metals={z:a.difference(works[z][e[z]['candidate']],works[z][b[z]['candidate']]) for z in METALS};delta=a.difference(metals['Ca'],metals['La']);raw=e['R_model_kcal_mol']-b['R_model_kcal_mol'];h.close(raw,delta['total'],residuals)
                changes[label][policy]={'per_metal_change_kcal_mol':metals,'differential_change_kcal_mol':delta,'R_model_change_kcal_mol':raw,'before':{z:b[z]['candidate'] for z in METALS},'after':{z:e[z]['candidate'] for z in METALS}}
        rows.append({**source,'pools':pools,**changes,'all_origin_relative_works_kcal_mol':works,'evidence':evidence,'failed_new_cells':fail,'classification':None})
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_sources':10,'declared_new_cells':40,'available_new_cells':available,'reused_original_cells':80,'rows':rows,
        'pool_orders':c['pool_orders'],'MACE_stage_gate':pf['MACE_stage_gate'],'saved_algebra_checks':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
        'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'classification':None,'new_calibration':False,'production_changed':False,
        'interpretation':'Finite alternative-basin composite response; no equilibrium, affinity, class or composite-stationary-point claim'}

def execute(path,sha,action):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports');cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['implementation']==pin(__file__),'Wrong frozen adapter');matrix=load('_relaxed_solvent_placement',c['matrix']);placement=matrix.placement()
    out=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID']);start=time.monotonic();pf=None;sr=[];receipt={'status':'failed','config':cp,'placement':placement,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'));old=load('_relaxed_solvent_prior',c['prior_adapter']);prior=read(c['prior_config'])
            f=load('_relaxed_solvent_factorial',c['factorial']);r=load('_relaxed_solvent_component',c['component']);n=load('_relaxed_solvent_native',c['native_adapter'])
            a=load('_relaxed_solvent_reuse',c['solvent_adapter']);h=load('_relaxed_solvent_pool',c['pool_adapter']);selector=load('_relaxed_solvent_selector',c['scorer'])
            score=load('_relaxed_solvent_old_score',prior['score_adapter']);relax=load('_relaxed_solvent_relaxation',c['score_adapter'])
            pf=preflight(c,old,f,r,n,a,h,selector,score,relax);write(out/'PREFLIGHT.json',pf)
            if action=='execute':sr=old.solvent(c,out,f,r,n,pf['reused_MACE'])
            receipt['status']='complete'
        except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc())
        finally:
            try:
                if pf and action=='execute':
                    if not sr:
                        for cell in c['cells']:
                            p=out/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                            if p.exists():
                                v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    write(out/'PARTIAL_RESULTS.json',{'config':cp,'native_rows':sr,'declared_new_cells':40})
                    result={'config':cp,**summarize(c,pf,sr,a,h,n,r,selector)};write(out/'FINAL_COLLECTION.json',result)
                    if result['available_new_cells']!=40 and receipt['status']=='complete':receipt['status']='partial'
            except BaseException as exc:receipt.update(status='failed',collection_failure={'reason':str(exc),'traceback':traceback.format_exc()})
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
                try:receipt.update(SCC_started=sum(1 for _ in out.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in out.glob('solvent/points/*/*/EXECUTION.json')),
                    native_component_calls=sum(1 for _ in out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
                    measured_resource_receipts={'SCC':[pin(p) for p in sorted(out.glob('solvent/points/*/*/EXECUTION.json'))],'native_ALPB':[pin(p) for p in sorted(out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]})
                except BaseException as exc:receipt.update(status='failed',resource_reporting_error=str(exc))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);a=s.add_parser('bind');a.add_argument('--score-config',required=True);a.add_argument('--score-collection',required=True);a.add_argument('--output',required=True)
    for action in ('preflight','execute'):
        a=s.add_parser(action);a.add_argument('--config',required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='bind':bind(a.score_config,a.score_collection,a.output)
    else:sys.exit(execute(a.config,a.config_sha,a.action))
