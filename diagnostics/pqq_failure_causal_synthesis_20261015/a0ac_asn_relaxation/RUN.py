#!/usr/bin/env python3
"""Finite flipped-origin searches through unchanged completion and WarmCPU kernels."""
import argparse, copy, fcntl, hashlib, importlib.util, json, os, shutil, signal, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation'
METALS=('Ca','La');ORDER3=('origin','adaptive_Ca','adaptive_La');ORDER4=ORDER3+('amide_flip',);ORDER6=ORDER4+('relaxed_flip_Ca','relaxed_flip_La')
PROTOCOL='A0AC_admitted_amide_basin_frozen_four_mode_relaxation_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(v,why):
    if not v:raise ValueError(why)
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def bind(preparation,stage,output,scout=None):
    p=read(preparation);pc=read(p['config']);old=read(pc['old_config']);template=read(old['worker_manifest'])
    require(p['protocol']==PROTOCOL and len(p['rows'])==10,'Wrong preparation scope');ids=pc['source_ids']
    require([r['source_case_id'] for r in p['rows']]==ids and ids[0]==pc['scout_source_id'],'Ordered source identity changed')
    selected=ids[:1] if stage=='scout' else ids[1:];rows=[r for r in p['rows'] if r['source_case_id'] in selected]
    require(all(r['status']=='prepared' for r in rows),'Selected source preparation unavailable')
    scout_pin=pin(scout) if scout else None
    if stage=='remaining':
        require(scout_pin is not None,'Actual first-source collection required');s=read(scout_pin)
        require(s['stage']=='scout' and s['stage_complete'] and s['selected_source_ids']==ids[:1],'Initial pair did not complete')
        sr=next(r for r in s['rows'] if r['source_case_id']==ids[0]);pr=p['rows'][0]
        require(sr['status']=='complete' and sr['prepared_mapping']['sha256']==pr['mapping']['sha256'],'Scout map differs or common pool unavailable')
    else:require(scout is None,'Scout cannot reuse another execution')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('RUN.py','RUN.sbatch'):shutil.copy2(HERE/name,out/name)
    tasks=[t for r in rows for t in r['tasks']];worker=copy.deepcopy(template)
    worker.update(protocol_id=PROTOCOL,tasks=tasks,cases=[],fixtures=[],origin_task_count=len(tasks),preparation=pin(preparation),
        finite_scope={'declared_all_sources':10,'stage_sources':len(rows),'stage_optimizer_starts':len(tasks),'stage_maximum_cross_calls':len(tasks),'qualification_calls':0,'SCC_calls':0,'DFT_calls':0})
    write(out/'worker_manifest.json',worker)
    c={'protocol':PROTOCOL,'implementation':pin(out/'RUN.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pc['plan'],'preparation':pin(preparation),'preparation_config':p['config'],
       'old_config':pc['old_config'],'old_collection':pc['old_collection'],'old_scorer':pc['old_scorer'],'source_preparation_helper':pc['source_preparation_helper'],
       'worker_manifest':pin(out/'worker_manifest.json'),'template':old['worker_manifest'],'python':pc['python'],'stage':stage,'selected_source_ids':selected,
       'all_source_ids':ids,'scout_collection':scout_pin,'rows':rows,'output':str(out/'execution'),'declared_all_searches':20,'reused_searches':20,'reused_endpoint_cells':80,
       'pool_orders':{'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)},'solvent_review_threshold_kcal_mol':1.0,'new_SCC_calls':0,'new_DFT_calls':0}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,b,scorer,np):
    p=read(c['preparation']);pc=read(c['preparation_config']);m=read(c['worker_manifest']);old=read(c['old_config']);template=read(c['template'])
    require(p['config']==c['preparation_config'] and c['protocol']==p['protocol']==PROTOCOL,'Preparation/config mismatch')
    require(c['plan']==pc['plan'] and pin(sys.executable)==c['python'],'Plan/interpreter changed');check(c['plan']);check(c['wrapper'])
    require([r['source_case_id'] for r in p['rows']]==c['all_source_ids']==pc['source_ids'],'Finite denominator changed')
    require(c['selected_source_ids']==(c['all_source_ids'][:1] if c['stage']=='scout' else c['all_source_ids'][1:]),'Stage sources changed')
    require(c['rows']==[r for r in p['rows'] if r['source_case_id'] in c['selected_source_ids']],'Prepared stage rows changed')
    require(c['pool_orders']=={'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)} and c['solvent_review_threshold_kcal_mol']==1.0,'Pool or gate changed')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        require(m[k]==template[k],'Frozen worker setting changed: '+k)
    require(m['tasks']==[t for r in c['rows'] for t in r['tasks']] and m['origin_task_count']==len(m['tasks']) and not m['cases'] and not m['fixtures'],'Worker task scope differs')
    oldpf=scorer.preflight(old,b);cached=read(c['old_collection']);require(cached['complete_new_cells']==20 and cached['config']==c['old_config'],'Old flip matrix differs')
    import slsqp_precision
    import scipy
    import scipy.optimize._slsqp_py as scipy_wrapper
    import scipy.optimize._slsqplib as scipy_kernel
    from mace_site_kinematics import Kinematics
    from adaptive_force_diagnostic import project
    completion=slsqp_precision.engine()[0].completion;require(completion.SETTINGS==m['settings'],'Completion numerical settings differ')
    require(scipy.__version__==m['optimizer_software']['version'] and pin(scipy_wrapper.__file__)==m['optimizer_software']['wrapper'] and pin(scipy_kernel.__file__)==m['optimizer_software']['kernel'],'Actual pinned SLSQP runtime differs')
    checked_origins=[]
    for r in p['rows']:
        baseline=next(x['matrix'] for x in cached['rows'] if x['source_case_id']==r['source_case_id'])
        for z in METALS:
            cell=baseline[z]['amide_flip'];value=read(cell['receipt']);req=read(value['request']);task={'xyz':cell['xyz'],'charge':-4 if z=='Ca' else -3,'multiplicity':1,'atom_count':204}
            scorer.admitted_result(value,req,task,m['model'],oldpf['versions']);require(value['energy_eV']==cell['energy_eV'] and value['forces']==cell['forces'],'Cached flip changed')
        if r['source_case_id'] not in c['selected_source_ids']:continue
        require(r['status']=='prepared' and read(r['admission'])['pass'],'Stage source admission unavailable');read(r['contact']['identities']);check(r['contact']['origin'])
        for z,t in zip(METALS,r['tasks']):
            require(t['metal']==z and t['mapping']==r['mapping'] and t['contact']==r['contact'],'Prepared task mismatch');g=read(t['mapping']);kin=Kinematics(g['context']);zero=np.zeros(len(kin.modes))
            oldtask=r['old_searches'][z]['task'];oldm=read(r['old_searches'][z]['manifest']);oldr=read(r['old_searches'][z]['result'])
            require(oldm['settings']==m['settings'] and oldr['status']=='proposal_available','Old search/settings changed')
            require(t['active_indices']==oldtask['active_indices'] and t['active_mode_ids']==oldtask['active_mode_ids'],'Original selected modes changed')
            point=t['origin_reuse']['point'];cell=baseline[z]['amide_flip'];require(point['MACE']==cell['receipt'] and point['forces']==cell['forces'] and point['coordinate']==cell['xyz']==t['xyz'],'Exact flipped q0 pin changed')
            force=np.load(check(point['forces']),allow_pickle=False);coords,_,raw,_,_=project(kin.data,zero,force)
            require(np.max(np.abs(coords-np.asarray([a[1:] for a in b.xyz(check(t['xyz']))])))<=1e-12,'Actual q0 replay differs')
            require(np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad']) and point['MACE_eV']==cell['energy_eV'],'Cached q0 gradient/energy changed')
            require(not any(point['full_q']) and kin.check(zero)['pass'],'Rebased q0/Jacobian unavailable');completion.angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in b.xyz(check(t['xyz']))])
            checked_origins.append(t['task_id'])
    if c['stage']=='remaining':
        scout=read(c['scout_collection']);require(scout['stage_complete'] and scout['selected_source_ids']==c['all_source_ids'][:1],'Scout not complete')
        sr=next(r for r in scout['rows'] if r['source_case_id']==c['all_source_ids'][0]);require(sr['prepared_mapping']['sha256']==p['rows'][0]['mapping']['sha256'] and sr['status']=='complete','Scout geometry differs')
    return {'status':'passed','versions':oldpf['versions'],'cached_q0_tasks':checked_origins,'all80_saved_cells_verified':True,'new_model_calls':0},completion

def guarded_class(completion,helper,np):
    original=completion.CompletionProposal
    class Guarded(original):
        def __init__(self,*args):
            super().__init__(*args);c=self.task['contact'];self.contact_spec=c
            with np.load(check(c['origin']),allow_pickle=False) as n:self.union=n['positions'].copy();self.symbols=n['symbols'].tolist()
            self.ids=read(c['identities'])['centers'];self.moving=c['moving_indices'];self.caps=set(c['caps'])
        def evaluate(self,q,purpose):
            q=np.asarray(q,dtype=float);require(q.shape==(4,) and np.isfinite(q).all() and np.max(np.abs(q))<=.8+1e-12,'Outside frozen angular domain')
            full=completion.angular.full_q(self.task,q);_,coords,_,_=self.kin.evaluate(full)
            require(np.max(np.abs(coords[list(self.caps)]-self.kin.evaluate(np.zeros(self.task['mode_count']))[1][list(self.caps)]))<=1e-12,'Cap anchors moved')
            changed=self.union.copy();changed[:204]=coords
            contacts=helper.contact_check(np,self.union,changed,self.symbols,self.moving,self.ids,self.caps)
            name=completion.cache_key({'manifest':self.manifest,'task':self.task['task_id'],'q':q.tolist()})
            gp=self.directory/'full_source_checks'/(name+'.json')
            gate={'task_id':self.task['task_id'],'active_q_radian':q.tolist(),'contact_inputs':self.contact_spec,'contacts':contacts,'pass':contacts['pass']}
            if gp.exists():require(read(gp)==gate,'Cached source-contact gate differs')
            else:write(gp,gate)
            require(contacts['pass'],'Full-source new severe overlap; exact failure retained')
            return super().evaluate(q,purpose)
    return original,Guarded

def source_run(c,row,worker,m,b,scorer,transport,completion,out):
    sid=row['source_case_id'];d=out/'sources'/sid;d.mkdir(parents=True);matrix=copy.deepcopy(row['baseline']);proposals={};reason=None
    result={'source_case_id':sid,'conditioning':row['conditioning'],'sample':row['sample'],'case_key':row['case_key'],'record_id':row['record_id'],
       'status':'unavailable','matrix':matrix,'proposals':proposals,'prepared_mapping':row['mapping'],'admission':row['admission']}
    for z in METALS:
        for name in ORDER6[4:]:matrix[z][name]={'status':'unavailable','energy_eV':None,'reason':'proposal not complete'}
    pm={'kind':'proposals','batch_manifest':worker,'model':m['model'],'settings':m['settings'],'tasks':row['tasks'],'prepared_admission':row['admission'],
        'original_searches':row['old_searches'],'frozen_selector_policy':'exact original four modes; no reselection'};pp=d/'proposal_manifest.json';write(pp,pm)
    try:
        for t in row['tasks']:
            rp=b.proposal(pp,t,transport,completion);r=read(rp);z=t['metal'];name='relaxed_flip_'+z;proposals[z]=rp
            if r['status']!='proposal_available':continue
            require(r['optimizer']['success'] and r['final_geometry']['physical_feasible'],'Noncredible optimizer proposal')
            point=r['proposal'];v=read(point['MACE']);task={**t,'xyz':point['coordinate']};scorer.admitted_result(v,read(v['request']),task,m['model'],read(c['preparation'])['versions'])
            qp=completion.cache_key({'manifest':pin(pp),'task':t['task_id'],'q':point['active_q_radian']});gate=Path(pp).parent/'proposals'/t['task_id']/'full_source_checks'/(qp+'.json')
            require(read(gate)['pass'],'Final full-source gate missing')
            matrix[z][name]={'status':'complete','receipt':point['MACE'],'xyz':point['coordinate'],'energy_eV':v['energy_eV'],'forces':v['forces'],'charge':t['charge'],'multiplicity':1,
                'atom_count':204,'proposal_receipt':rp,'full_source_admission':pin(gate),'boundary_flag':r['boundary_flag']}
            other='La' if z=='Ca' else 'Ca';ot=next(t0 for t0 in row['tasks'] if t0['metal']==other);cd=d/'cross_pool'/name;cd.mkdir(parents=True)
            lines=check(point['coordinate']).read_text().splitlines(keepends=True);lines[2]=other+' '+lines[2].split(maxsplit=1)[1];xp=cd/'context.xyz';xp.write_text(''.join(lines))
            ct={**ot,'task_id':sid+'__'+name+'__'+other,'xyz':pin(xp)};cm={'kind':'cross_pool','batch_manifest':worker,'model':m['model'],'tasks':[ct],'candidate':point['coordinate']};cp=cd/'manifest.json';write(cp,cm)
            req={'manifest':pin(cp),'task_id':ct['task_id'],**{k:ct[k] for k in ('xyz','charge','multiplicity')}}
            v,vr=b.cell(cd,req,transport,m['model']);matrix[other][name].update(receipt=vr,xyz=ct['xyz'],reason=v.get('reason'))
            if v['status']=='complete':
                scorer.admitted_result(v,read(v['request']),ct,m['model'],read(c['preparation'])['versions'])
                require(b.xyz(xp)[1:]==b.xyz(check(point['coordinate']))[1:] and b.xyz(xp)[0][1:]==b.xyz(check(point['coordinate']))[0][1:],'Cross physical geometry changed')
                matrix[other][name].update(status='complete',energy_eV=v['energy_eV'],forces=v['forces'],charge=ct['charge'],multiplicity=1,atom_count=204,proposal_receipt=rp,full_source_admission=pin(gate))
        result['status']='complete' if all(matrix[z][n]['status']=='complete' for z in METALS for n in ORDER6) and len(proposals)==2 else 'unavailable'
    except BaseException as e:
        result.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:write(d/'SOURCE_RESULT.json',result)
    return result

def collect(c,out,scorer):
    baseline=read(c['old_collection']);prepared=read(c['preparation']);prior=read(c['scout_collection']) if c['scout_collection'] else None;rows=[]
    for old in baseline['rows']:
        sid=old['source_case_id'];path=out/'sources'/sid/'SOURCE_RESULT.json';pr=next(r for r in prepared['rows'] if r['source_case_id']==sid)
        r=read(path) if path.exists() else (copy.deepcopy(next(x for x in prior['rows'] if x['source_case_id']==sid)) if prior and sid==c['all_source_ids'][0] else
            {'source_case_id':sid,'conditioning':old['conditioning'],'sample':old['sample'],'status':'not_executed','matrix':copy.deepcopy(old['matrix']),'preparation_status':pr['status'],'reason':pr.get('reason'),'proposals':{}})
        for z in METALS:
            for name in ORDER6[4:]:r['matrix'][z].setdefault(name,{'status':'unavailable','energy_eV':None,'reason':'not executed'})
        energies={z:{n:(v.get('energy_eV') if v['status']=='complete' else None) for n,v in r['matrix'][z].items()} for z in METALS}
        r['pools']={label:scorer.choose(energies,order) for label,order in [('common3',ORDER3),('common4',ORDER4),('common6',ORDER6)]};r['common4_to6']={}
        for policy in ('mathematical','operational'):
            a,e=r['pools']['common4'],r['pools']['common6'];r['common4_to6'][policy]=None
            if a and e:
                change={z:e[policy]['per_metal_work_kcal_mol'][z]-a[policy]['per_metal_work_kcal_mol'][z] for z in METALS}
                r['common4_to6'][policy]={'per_metal_change_kcal_mol':change,'paired_change_kcal_mol':change['Ca']-change['La']}
        rows.append(r)
    complete=all(r['status']=='complete' for r in rows);stage_complete=all(r['status']=='complete' for r in rows if r['source_case_id'] in c['selected_source_ids'])
    flag=[r['source_case_id'] for r in rows if r['common4_to6']['mathematical'] is not None and abs(r['common4_to6']['mathematical']['paired_change_kcal_mol'])>=1.0]
    return {'protocol':PROTOCOL,'config':pin(Path(out).parent/'CONFIG.json'),'preparation':c['preparation'],'stage':c['stage'],'selected_source_ids':c['selected_source_ids'],
       'stage_complete':stage_complete,'all20_searches_and_common6_complete':complete,'declared_sources':10,'declared_new_searches':20,'reused_old_searches':20,'reused_endpoint_cells':80,
       'pool_orders':c['pool_orders'],'rows':rows,'solvent_review':{'eligible':complete and bool(flag),'source_ids_above_sign_blind_gate':flag,'threshold_kcal_mol':1.0,
           'criterion':'complete credible common6 and absolute mathematical common4-to6 paired change','automatic_submission':False},'new_SCC_calls':0,'new_DFT_calls':0,'classification':None}

def run(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before numerical imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(pin(__file__)==c['implementation'],'Wrong frozen runner')
    old=read(c['old_config']);placement=load('_asn_placement',old['prior_adapter']).placement();scorer=load('_asn_scoring',c['old_scorer'])
    out=Path(c['output']);m=read(c['worker_manifest']);b=load('_asn_warm_worker',m['implementation']);start=time.monotonic();cpu=time.process_time();transport=None;handlers={};completion=None;original_class=None
    receipt={'status':'failed','config':cp,'placement':placement,'stage':c['stage'],'new_SCC_calls':0,'new_DFT_calls':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            b.parent_env();b.modules(m)
            import numpy as np
            pf,completion=preflight(c,b,scorer,np);write(out/'PREFLIGHT.json',pf)
            helper=load('_asn_full_source_guard',c['source_preparation_helper']);original_class,guarded=guarded_class(completion,helper,np);completion.CompletionProposal=guarded
            attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers();transport=b.WarmCPU(check(c['worker_manifest']),m,attempt,0)
            require(transport.ready['threads']==7 and transport.ready['versions']==pf['versions'],'Warm runtime differs')
            for row in c['rows']:source_run(c,row,c['worker_manifest'],m,b,scorer,transport,completion,out)
            transport.close();transport=None;receipt['status']='complete'
        except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc());raise
        finally:
            try:
                if transport:transport.abort()
            except BaseException as e:receipt.update(status='failed',cleanup_error=str(e))
            try:
                if completion is not None and original_class is not None:completion.CompletionProposal=original_class
                for sig,handler in handlers.items():signal.signal(sig,handler)
                collection=collect(c,out,scorer);write(out/'FINAL_COLLECTION.json',collection)
                if not collection['stage_complete']:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,parent_CPU_seconds=time.process_time()-cpu)
                try:receipt.update(actual_model_forwards=sum(read(p).get('actual_model_forwards',0) for p in out.rglob('mace_result.json')),worker_receipts=[pin(p) for p in out.glob('mace/lane0/CPU_WORKER_SUMMARY.json')])
                except BaseException as e:receipt.update(status='failed',cost_reporting_error=str(e))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);b=s.add_parser('bind');b.add_argument('--preparation',required=True);b.add_argument('--stage',choices=['scout','remaining'],required=True);b.add_argument('--output',required=True);b.add_argument('--scout-collection')
    r=s.add_parser('run');r.add_argument('--config',required=True);r.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='bind':bind(a.preparation,a.stage,a.output,a.scout_collection)
    else:sys.exit(run(a.config,a.config_sha))
