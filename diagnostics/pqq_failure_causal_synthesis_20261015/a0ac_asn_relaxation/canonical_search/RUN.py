#!/usr/bin/env python3
"""Canonical25 adapter for unchanged warm MACE and bounded CompletionProposal."""
import argparse, copy, fcntl, hashlib, importlib.util, json, math, os, shutil, signal, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_search'
PROTOCOL='canonical25_terminal_amide_common6_v1'
METALS=('Ca','La');ORDER3=('origin','adaptive_Ca','adaptive_La');ORDER4=ORDER3+('amide_flip',);ORDER6=ORDER4+('relaxed_flip_Ca','relaxed_flip_La')
ORDERS={'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)}
REUSE='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
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

def task_template(row,t):
    t=copy.deepcopy(t);t.update(task_id=row['source_case_id']+'__canonical_amide_flip__'+t['metal'],case_id=row['case_id']+'__amide_flip')
    return t

def bind(preparation,stage,output,scout=None):
    """Only metadata/hashes/copies; no scientific imports or coordinate parsing."""
    p=read(preparation);gm=read(p['manifest']);prior=read(gm['prior_config']);old=read(prior['old_config']);template=read(gm['worker_manifest'])
    require(p['protocol']=='canonical25_fixed_flip_rebased_geometry_only_v1' and p['status']=='complete' and p['prepared_sources']==25,'Actual all25 geometry required')
    ids=[r['source_case_id'] for r in p['rows']];require(len(ids)==len(set(ids))==25 and ids==[s['row']['source_case_id'] for s in gm['sources']],'Canonical geometry denominator differs')
    require(gm['reuse']['source_case_id']==REUSE and read(gm['completed_search_execution'])['status']=='complete','Exact executed A0AC reuse unavailable')
    reuse=next(r for r in p['rows'] if r['source_case_id']==REUSE);require(read(reuse['regression'])['status']=='passed' and reuse['executed_result_alias']==gm['reuse']['executed_row'],'Actual A0AC regression/alias required')
    cc=read(gm['canonical_collection']);inv=read(cc['inventory']);designated=inv['calibration_rule']['designated_record_ids']
    require([r['record_id'] for r in p['rows']]==designated==[r['record_id'] for r in cc['rows']],'Exact designated canonical membership/order changed')
    newids=[sid for sid in ids if sid!=REUSE];selected=newids[:1] if stage=='scout' else newids[1:]
    scoutpin=pin(scout) if scout else None;scoutexec=None
    if stage=='remaining':
        require(scoutpin is not None,'Actual scout required');s=read(scoutpin);scoutexec=pin(Path(scout).parent/'EXECUTION.json');se=read(scoutexec);sc=read(s['config'])
        require(s['protocol']==PROTOCOL and s['stage']=='scout' and s['stage_complete'] and s['selected_source_ids']==newids[:1],'Incomplete scout cannot release remaining sources')
        require(se['status']=='complete' and se['config']==s['config'] and sc['preparation']==pin(preparation) and sc['plan']['sha256']==pin(HERE/'PLAN.md')['sha256'],'Scout execution/protocol differs')
    else:require(scout is None,'Scout cannot reuse another stage')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('RUN.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/name,out/name)
    rows=[r for r in p['rows'] if r['source_case_id'] in selected];tasks=[task_template(r,t) for r in rows for t in r['task_templates']]
    worker=copy.deepcopy(template);worker.update(protocol_id=PROTOCOL,tasks=tasks,cases=[],fixtures=[],origin_task_count=len(tasks),preparation=pin(preparation),
        finite_scope={'all_sources':25,'new_sources':24,'stage_sources':len(rows),'fixed_flip_calls':len(tasks),'optimizer_starts':len(tasks),'maximum_cross_calls':len(tasks),'maxiter':200,'new_SCC_calls':0,'new_DFT_calls':0})
    write(out/'worker_manifest.json',worker)
    c={'protocol':PROTOCOL,'implementation':pin(out/'RUN.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(out/'PLAN.md'),'preparation':pin(preparation),
       'preparation_manifest':p['manifest'],'geometry_execution':pin(Path(preparation).parent/'EXECUTION.json'),'old_config':prior['old_config'],'old_scorer':prior['old_scorer'],
       'generic_contact_helper':gm['generic_contact_helper'],'worker_manifest':pin(out/'worker_manifest.json'),'template':gm['worker_manifest'],'python':gm['python'],
       'canonical_collection':gm['canonical_collection'],'calibration_inventory':cc['inventory'],'calibration_scorer':cc['scorer'],
       'calibration_plan':pin(ROOT/'diagnostics/pede_solvent_response_20261013/frozen_distribution_candidate/PLAN.md'),
       'stage':stage,'selected_source_ids':selected,'all_source_ids':ids,'new_source_ids':newids,'reuse_source_id':REUSE,'rows':rows,
       'scout_collection':scoutpin,'scout_execution':scoutexec,'pool_orders':ORDERS,'output':str(out/'execution'),
       'new_fixed_flip_calls_total':48,'new_optimizer_starts_total':48,'maximum_new_cross_calls_total':48,'reused_original_cells':150,'reused_additional_A0AC_cells':6,
       'new_SCC_calls':0,'new_DFT_calls':0,'classification':None,'automatic_follow_on':False}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def physical_result(result,request,task,model,observed_model):
    """Original physical predicate, generalized only to the actual source count."""
    require(result['status'] in ('computed','complete') and observed_model==model,'Result/model unavailable')
    require((request['charge'],request.get('multiplicity',request.get('spin_multiplicity')))==(task['charge'],task['multiplicity']),'Request state differs')
    require(request['xyz']['sha256']==task['xyz']['sha256'],'Exact request coordinates differ');check(request['xyz']);check(result['forces'])
    state=result.get('state',result.get('input_state_check'));n=task['atom_count']
    require((state['charge'],state['spin_multiplicity'],state['atoms'],state['graph_count'],state['head'],state['selected_metal_index'])==(task['charge'],1,n,1,'omol',0),'Actual source-count OMOL state differs')
    require(result['force_definition'] in ('negative Cartesian total vacuum OMOL energy gradient, eV/angstrom','negative Cartesian energy gradient, eV/angstrom','negative_Cartesian_gradient_of_total_vacuum_OMOL_energy') and math.isfinite(result['energy_eV']),'Force convention/nonfinite energy')

def archive(cell,task,m,scorer,versions,np):
    rp=cell.get('MACE',cell.get('receipt'));r=read(rp)
    if r.get('request'):request=read(r['request']);observed=r['model']
    else:
        manifest=read(r['manifest']);request,=[t for t in manifest['tasks'] if t['task_id']==r['task_id']];observed=manifest['model']
    physical_result(r,request,task,m['model'],observed)
    require(r['energy_eV']==cell.get('MACE_eV',cell.get('energy_eV')),'Archived energy changed')
    if cell.get('forces'):require(cell['forces']==r['forces'],'Archived force pin changed')
    force=np.load(check(r['forces']),allow_pickle=False);require(force.shape==(task['atom_count'],3) and np.isfinite(force).all(),'Actual archived force shape/values invalid')
    if r['device']=='cpu':scorer.admitted_result(r,request,task,m['model'],versions)
    else:require(r['device']=='NVIDIA H200','Unexpected archived execution device')
    return {**cell,'status':'complete','receipt':rp,'xyz':task['xyz'],'energy_eV':r['energy_eV'],'forces':r['forces'],'request':r.get('request'),
       'historical_manifest':r.get('manifest'),'historical_device':r['device'],'charge':task['charge'],'multiplicity':1,'atom_count':task['atom_count'],'reused':True}

def preflight(c,b,scorer,np):
    p=read(c['preparation']);gm=read(c['preparation_manifest']);m=read(c['worker_manifest']);old=read(c['old_config']);template=read(c['template'])
    require(pin(sys.executable)==c['python'] and p['manifest']==c['preparation_manifest'],'Interpreter/preparation differs')
    for k in ('plan','wrapper','calibration_plan','calibration_scorer','generic_contact_helper'):check(c[k])
    require(read(c['geometry_execution'])['status']=='complete' and p['prepared_sources']==25 and c['pool_orders']==ORDERS,'Geometry/pool gate differs')
    require([r['source_case_id'] for r in p['rows']]==c['all_source_ids'] and c['new_source_ids']==[sid for sid in c['all_source_ids'] if sid!=REUSE],'All25 denominator changed')
    require(c['selected_source_ids']==(c['new_source_ids'][:1] if c['stage']=='scout' else c['new_source_ids'][1:]),'Stage scope differs')
    require(c['rows']==[r for r in p['rows'] if r['source_case_id'] in c['selected_source_ids']],'Stage geometry changed')
    cc=read(c['canonical_collection']);inv=read(c['calibration_inventory'])
    require(cc['inventory']==c['calibration_inventory'] and [r['record_id'] for r in p['rows']]==inv['calibration_rule']['designated_record_ids']==[r['record_id'] for r in cc['rows']],'Canonical identity changed')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        require(m[k]==template[k],'Frozen worker setting differs: '+k)
    for q in m['source_modules'].values():check(q)
    for k in ('implementation','calculator_source','model_executable','optimizer_executable'):check(m[k])
    check(m['model']['checkpoint'])
    require(m['tasks']==[task_template(r,t) for r in c['rows'] for t in r['task_templates']] and m['origin_task_count']==len(m['tasks']) and m['threads_per_model']==7 and m['lanes']==1 and not m['cases'] and not m['fixtures'],'Worker dispatch differs')
    qualification=read(old['qualification']);summary=read(old['qualified_worker_summary'])
    require(qualification['pass'] and summary['threads']==7 and summary['manifest']==old['template'],'Actual CPU qualification differs')
    for q in qualification['checks']:
        r=read(q['receipt']);check(r['forces']);read(r['request']);require(q['pass'] and r['model']==m['model'] and r['versions']==summary['versions'],'Actual qualification fixture differs')
    versions=summary['versions']
    import slsqp_precision, scipy
    import scipy.optimize._slsqp_py as sw
    import scipy.optimize._slsqplib as sk
    from mace_site_kinematics import Kinematics
    from adaptive_force_diagnostic import project
    completion=slsqp_precision.engine()[0].completion
    require(completion.SETTINGS==m['settings'] and scipy.__version__==m['optimizer_software']['version'] and pin(sw.__file__)==m['optimizer_software']['wrapper'] and pin(sk.__file__)==m['optimizer_software']['kernel'],'Actual optimizer differs')
    for module in list(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename and Path(filename).name in m['source_modules']:require(pin(filename)==m['source_modules'][Path(filename).name],'Actually imported scientific source differs')
    baseline={};checks=[];nchecked=0
    for r in p['rows']:
        sid=r['source_case_id'];s,=[s['scope'] for s in gm['sources'] if s['row']['source_case_id']==sid];n=r['atom_count'];st=read(s['state_signature'])
        original_search=read(s['search_manifest']);require(original_search['settings']==m['settings'] and [t for t in original_search['tasks'] if t['case_id']==s['case_id']]==s['original_tasks'],'Original executed task/settings differ')
        require(st['endpoint_charges']==r['endpoint_charges'] and st['multiplicities']=={'Ca':1,'La':1} and st['water_inventory']==r['water_inventory'],'Source state differs')
        require(r['original_row']['state_signature']==s['state_signature'] and r['original_row']['maps']==s['maps'],'Original paired source mapping differs')
        for z in METALS:check(s['maps'][z])
        check(s['source_preparation']);check(s['protonated_source']);require(read(r['admission'])['pass'] and r['status']=='prepared','Prepared source unavailable')
        baseline[sid]={z:{} for z in METALS};coordinates={}
        for z in METALS:
            coordinates[z]={}
            for q in ORDER3:
                cell=s['six_saved_cells'][z][q];t={'xyz':cell['xyz'],'charge':r['endpoint_charges'][z],'multiplicity':1,'atom_count':n}
                coords=b.xyz(check(t['xyz']));require(len(coords)==n and coords[0][0]==z,'Archived count/metal differs');coordinates[z][q]=coords
                baseline[sid][z][q]=archive(cell,t,m,scorer,versions,np);nchecked+=1
        for q in ORDER3:
            require(coordinates['Ca'][q][1:]==coordinates['La'][q][1:] and coordinates['Ca'][q][0][1:]==coordinates['La'][q][0][1:],'Original common-pool paired coordinates differ')
        kin=Kinematics(read(r['mapping'])['context']);zero=np.zeros(len(kin.modes));require(kin.check(zero)['pass'],'Rebased geometry/Jacobian failed')
        check(r['contact']['origin']);read(r['contact']['identities']);check(r['contact']['source'])
        for t in r['task_templates']:
            z=t['metal'];ot,=[x for x in s['original_tasks'] if x['metal']==z]
            require(t['mapping']==r['mapping'] and t['contact']==r['contact'] and t['settings']==m['settings'] and t['atom_count']==n and t['charge']==r['endpoint_charges'][z] and t['multiplicity']==1,'Task state/geometry/settings differs')
            for k in ('active_indices','active_mode_ids','active_roles','selector','mode_count','selection_id'):require(t[k]==ot[k],'Original mode selection changed: '+k)
            a=b.xyz(check(t['xyz']));require(len(a)==n and a[0][0]==z and np.max(np.abs(kin.evaluate(zero)[1]-np.asarray([v[1:] for v in a])))<=1e-12,'Rebased q0 geometry differs')
            completion.angular.final_geometry(kin,t,np.zeros(4),[v[0] for v in a])
        checks.append({'source_case_id':sid,'atom_count':n,'six_archived_cells_pass':True,'rebased_geometry_pass':True,'terminal_mode_selected_by_metal':r['terminal_mode_selected_by_metal']})
    require(nchecked==150,'Original cell denominator differs')
    reused=next(r for r in p['rows'] if r['source_case_id']==REUSE);alias=reused['executed_result_alias'];final=read(gm['completed_search_collection'])
    require(read(gm['completed_search_execution'])['status']=='complete' and alias==next(r for r in final['rows'] if r['source_case_id']==REUSE),'Executed A0AC alias differs')
    require(read(reused['regression'])['status']=='passed' and alias['status']=='complete' and alias['prepared_mapping']==reused['mapping'] and alias['admission']==reused['prepared_state_alias'],'A0AC actual geometry/receipt differs')
    for z in METALS:
        t,=[t for t in reused['reused_prepared_tasks'] if t['metal']==z]
        require(t['mapping']==reused['mapping'] and t['contact']==reused['contact'],'A0AC prepared task alias differs')
        point=t['origin_reuse']['point'];flip=alias['matrix'][z]['amide_flip'];kin=Kinematics(read(t['mapping'])['context']);zero=np.zeros(len(kin.modes))
        require(point['MACE']==flip['receipt'] and point['coordinate']==flip['xyz']==t['xyz'] and point['forces']==flip['forces'] and point['MACE_eV']==flip['energy_eV'],'A0AC cached q0 pins differ')
        _,_,raw,_,_=project(kin.data,zero,np.load(check(point['forces']),allow_pickle=False));require(np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad']) and not any(point['full_q']),'A0AC actual cached gradient differs')
        for q in ORDER6:
            cell=alias['matrix'][z][q]
            if q in ORDER3:require(cell['receipt']==baseline[REUSE][z][q]['receipt'] and cell['xyz']==baseline[REUSE][z][q]['xyz'],'A0AC baseline alias differs')
            else:archive(cell,{**t,'xyz':cell['xyz']},m,scorer,versions,np)
        proposal=read(alias['proposals'][z]);require(proposal['status']=='proposal_available' and proposal['optimizer']['success'] and proposal['final_geometry']['physical_feasible'],'Reused A0AC optimizer unavailable')
    if c['stage']=='remaining':
        scout=read(c['scout_collection']);se=read(c['scout_execution']);sc=read(scout['config'])
        require(scout['stage_complete'] and se['status']=='complete' and se['config']==scout['config'] and sc['preparation']==c['preparation'] and sc['implementation']['sha256']==c['implementation']['sha256'] and sc['plan']['sha256']==c['plan']['sha256'],'Scout execution/source differs')
        sr=next(r for r in scout['rows'] if r['source_case_id']==c['new_source_ids'][0]);pr=p['rows'][0]
        require(sr['status']=='complete' and sr['prepared_mapping']==pr['mapping'],'Scout final map differs')
        for z in METALS:
            t=next(t for t in pr['task_templates'] if t['metal']==z)
            for q in ORDER6:archive(sr['matrix'][z][q],{**t,'xyz':sr['matrix'][z][q]['xyz']},m,scorer,versions,np)
    return {'status':'passed','versions':versions,'all150_saved_cells_verified':True,'A0AC_six_additional_cells_verified':True,'source_checks':checks,'new_model_calls':0,'baseline':baseline},completion

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
            changed=self.union.copy();changed[:self.task['atom_count']]=coords
            contacts=helper.contact_check(np,self.union,changed,self.symbols,self.moving,self.ids,self.caps,self.task['atom_count'])
            name=completion.cache_key({'manifest':self.manifest,'task':self.task['task_id'],'q':q.tolist()});gp=self.directory/'full_source_checks'/(name+'.json')
            gate={'task_id':self.task['task_id'],'active_q_radian':q.tolist(),'contact_inputs':self.contact_spec,'contacts':contacts,'pass':contacts['pass']}
            if gp.exists():require(read(gp)==gate,'Cached source-contact gate differs')
            else:write(gp,gate)
            require(contacts['pass'],'Full-source new severe overlap; exact failure retained')
            return super().evaluate(q,purpose)
    return original,Guarded

def source_run(c,row,m,b,scorer,transport,completion,np,pf,out):
    from mace_site_kinematics import Kinematics
    from adaptive_force_diagnostic import project
    sid=row['source_case_id'];d=out/'sources'/sid;d.mkdir(parents=True);matrix=copy.deepcopy(pf['baseline'][sid]);proposals={};tasks=[]
    result={**{k:row[k] for k in ('source_case_id','case_id','record_id','selection_id','protein','atom_count','endpoint_charges','terminal_mode_selected_by_metal','original_four_modes')},
        'status':'unavailable','matrix':matrix,'proposals':proposals,'prepared_mapping':row['mapping'],'admission':row['admission'],'reason':'not completed'}
    for z in METALS:
        for name in ORDER6[3:]:matrix[z][name]={'status':'unavailable','energy_eV':None,'reason':'not completed'}
    try:
        kin=Kinematics(read(row['mapping'])['context']);zero=np.zeros(len(kin.modes));_,core,_,cj=kin.evaluate(zero)
        for tt in row['task_templates']:
            t=task_template(row,tt);z=t['metal'];cd=d/'fixed_flip'/z
            req={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
            v,vr=b.cell(cd,req,transport,m['model']);matrix[z]['amide_flip'].update(receipt=vr,xyz=t['xyz'],reason=v.get('reason'))
            require(v['status']=='complete','Fixed flip unavailable: '+str(v.get('reason')));scorer.admitted_result(v,read(v['request']),t,m['model'],pf['versions'])
            matrix[z]['amide_flip'].update(status='complete',energy_eV=v['energy_eV'],forces=v['forces'],charge=t['charge'],multiplicity=1,atom_count=t['atom_count'])
            force=np.load(check(v['forces']),allow_pickle=False);require(force.shape==(t['atom_count'],3) and np.isfinite(force).all(),'Fresh force array invalid')
            coords,_,raw,_,_=project(kin.data,zero,force);active=t['active_indices'];gradient=-np.einsum('mij,ij->m',cj[active],force)*scorer.EV
            require(np.max(np.abs(coords-np.asarray([a[1:] for a in b.xyz(check(t['xyz']))])))<=1e-12 and np.array_equal(gradient,raw[active]),'Fresh q0 geometry/force contraction failed')
            point={'status':'complete','coordinate':t['xyz'],'MACE':vr,'MACE_eV':v['energy_eV'],'forces':v['forces'],'full_q':zero.tolist(),'active_q_radian':[0.]*4,'gradient_kcal_mol_rad':gradient.tolist(),'reused_scientific_origin':True}
            t.update(q0_status='native_only_composite_unavailable',q0={'components':{'MACE_eV':v['energy_eV']},'native_MACE_receipt':vr,'low':None},
                origin_reuse={'proposal_receipt':vr,'point':point,'receipt_kind':'actual_completed_flip_with_rebased_projection'},selector_policy='retain_exact_original_four_IDs_no_reselection')
            tasks.append(t)
        write(d/'Q0_QUALIFICATION.json',{'status':'passed','source_case_id':sid,'tasks':tasks,'new_fixed_calls':2,'q0_repeated_model_calls':0,'same_contraction':True})
        pp=d/'proposal_manifest.json';write(pp,{'kind':'proposals','batch_manifest':c['worker_manifest'],'model':m['model'],'settings':m['settings'],'tasks':tasks,'prepared_admission':row['admission'],'original_source':row['original_row'],'frozen_selector_policy':'exact original four modes; no reselection'})
        for t in tasks:
            rp=b.proposal(pp,t,transport,completion);r=read(rp);z=t['metal'];name='relaxed_flip_'+z;proposals[z]=rp
            if r['status']!='proposal_available':matrix[z][name]['reason']=r.get('reason',r['status']);continue
            require(r['optimizer']['success'] and r['final_geometry']['physical_feasible'],'Noncredible optimizer proposal')
            point=r['proposal'];v=read(point['MACE']);scorer.admitted_result(v,read(v['request']),{**t,'xyz':point['coordinate']},m['model'],pf['versions'])
            qp=completion.cache_key({'manifest':pin(pp),'task':t['task_id'],'q':point['active_q_radian']});gate=d/'proposals'/t['task_id']/'full_source_checks'/(qp+'.json');require(read(gate)['pass'],'Final full-source gate missing')
            matrix[z][name]={'status':'complete','receipt':point['MACE'],'xyz':point['coordinate'],'energy_eV':v['energy_eV'],'forces':v['forces'],'charge':t['charge'],'multiplicity':1,
                'atom_count':t['atom_count'],'proposal_receipt':rp,'full_source_admission':pin(gate),'boundary_flag':r['boundary_flag']}
            other='La' if z=='Ca' else 'Ca';ot=next(t0 for t0 in tasks if t0['metal']==other);cd=d/'cross_pool'/name;cd.mkdir(parents=True)
            lines=check(point['coordinate']).read_text().splitlines(keepends=True);lines[2]=other+' '+lines[2].split(maxsplit=1)[1];xp=cd/'context.xyz';xp.write_text(''.join(lines))
            ct={**ot,'task_id':sid+'__'+name+'__'+other,'xyz':pin(xp)};cp=cd/'manifest.json';write(cp,{'kind':'cross_pool','batch_manifest':c['worker_manifest'],'model':m['model'],'tasks':[ct],'candidate':point['coordinate']})
            req={'manifest':pin(cp),'task_id':ct['task_id'],**{k:ct[k] for k in ('xyz','charge','multiplicity')}}
            v,vr=b.cell(cd,req,transport,m['model']);matrix[other][name].update(receipt=vr,xyz=ct['xyz'],reason=v.get('reason'))
            if v['status']=='complete':
                scorer.admitted_result(v,read(v['request']),ct,m['model'],pf['versions'])
                require(b.xyz(xp)[1:]==b.xyz(check(point['coordinate']))[1:] and b.xyz(xp)[0][1:]==b.xyz(check(point['coordinate']))[0][1:],'Cross physical geometry changed')
                matrix[other][name].update(status='complete',energy_eV=v['energy_eV'],forces=v['forces'],charge=ct['charge'],multiplicity=1,atom_count=ct['atom_count'],proposal_receipt=rp,full_source_admission=pin(gate))
        result['status']='complete' if all(matrix[z][n]['status']=='complete' for z in METALS for n in ORDER6) and len(proposals)==2 else 'unavailable'
        result['reason']=None if result['status']=='complete' else 'Required fixed/proposal/cross cells unavailable'
    except Exception as e:result.update(reason=str(e),traceback=traceback.format_exc())
    finally:write(d/'SOURCE_RESULT.json',result)
    return result

def collect(c,out,scorer,pf):
    prep=read(c['preparation']);gm=read(c['preparation_manifest']);prior=read(c['scout_collection']) if c['scout_collection'] else None;rows=[]
    canonical={r['record_id']:r for r in read(c['canonical_collection'])['rows']}
    for row in prep['rows']:
        sid=row['source_case_id'];path=out/'sources'/sid/'SOURCE_RESULT.json'
        if path.exists():r=read(path)
        elif sid==REUSE:r=copy.deepcopy(row['executed_result_alias']);r.update(reused_executed_collection=gm['completed_search_collection'],geometry_regression=row['regression'],canonical_geometry_admission=row['admission'])
        elif prior and sid==c['new_source_ids'][0]:r=copy.deepcopy(next(x for x in prior['rows'] if x['source_case_id']==sid));r['reused_scout_collection']=c['scout_collection']
        else:
            s=next(s['scope'] for s in gm['sources'] if s['row']['source_case_id']==sid)
            baseline=copy.deepcopy(pf['baseline'][sid]) if pf else {z:{q:{'status':'unavailable','reason':'numerical preflight unavailable','saved_input':s['six_saved_cells'][z][q]} for q in ORDER3} for z in METALS}
            r={'source_case_id':sid,'status':'not_executed','matrix':baseline,'reason':'Stage not executed','proposals':{}}
        for k in ('case_id','record_id','selection_id','protein','atom_count','endpoint_charges','terminal_mode_selected_by_metal','original_four_modes'):r[k]=row[k]
        for k in ('expected_class','biological_group','source_conditioning_metal'):r[k]=canonical[row['record_id']][k]
        for z in METALS:
            for name in ORDER6:r['matrix'][z].setdefault(name,{'status':'unavailable','energy_eV':None,'reason':'not executed'})
        energies={z:{n:v.get('energy_eV') if v['status']=='complete' else None for n,v in r['matrix'][z].items()} for z in METALS}
        r['pools']={label:scorer.choose(energies,order) for label,order in ORDERS.items()};r['changes']={}
        for a,e in [('common3','common4'),('common4','common6'),('common3','common6')]:
            r['changes'][a+'_to_'+e]={}
            for policy in ('mathematical','operational'):
                before,after=r['pools'][a],r['pools'][e];change=None
                if before and after:
                    dz={z:after[policy]['per_metal_work_kcal_mol'][z]-before[policy]['per_metal_work_kcal_mol'][z] for z in METALS}
                    change={'per_metal_change_kcal_mol':dz,'paired_change_kcal_mol':dz['Ca']-dz['La']}
                r['changes'][a+'_to_'+e][policy]=change
        r.pop('common4_to6',None);rows.append(r)
    stage_complete=bool(pf) and all(r['status']=='complete' for r in rows if r['source_case_id'] in c['selected_source_ids'])
    return {'protocol':PROTOCOL,'config':pin(Path(out).parent/'CONFIG.json'),'preparation':c['preparation'],'stage':c['stage'],'stage_complete':stage_complete,
        'selected_source_ids':c['selected_source_ids'],'declared_sources':25,'all25_common6_complete':bool(pf) and all(r['status']=='complete' for r in rows),
        'new_fixed_flip_calls_total':48,'new_optimizer_starts_total':48,'maximum_new_cross_calls_total':48,'reused_original_cells':150,'reused_additional_A0AC_cells':6,
        'pool_orders':ORDERS,'rows':rows,'classification':None,'new_SCC_calls':0,'new_DFT_calls':0,'automatic_follow_on':False,'stage_cost_separate_from_reused_receipts':True}

def run(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(c['protocol']==PROTOCOL and pin(__file__)==c['implementation'],'Wrong frozen runner')
    old=read(c['old_config']);placement=load('_canonical_placement',old['prior_adapter']).placement();scorer=load('_canonical_scorer',c['old_scorer']);scorer.physical_result=physical_result
    out=Path(c['output']);m=read(c['worker_manifest']);b=load('_canonical_warm',m['implementation']);start=time.monotonic();cpu=time.process_time();transport=None;handlers={};completion=None;original_class=None;pf=None
    receipt={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'placement':placement,'stage':c['stage'],'new_SCC_calls':0,'new_DFT_calls':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            b.parent_env();b.modules(m)
            import numpy as np
            pf,completion=preflight(c,b,scorer,np);write(out/'PREFLIGHT.json',pf)
            helper=load('_canonical_contacts',c['generic_contact_helper']);original_class,guarded=guarded_class(completion,helper,np);completion.CompletionProposal=guarded
            attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers();transport=b.WarmCPU(check(c['worker_manifest']),m,attempt,0)
            require(transport.ready['threads']==7 and transport.ready['versions']==pf['versions'],'Actual warm runtime differs')
            partial=[]
            for row in c['rows']:
                result=source_run(c,row,m,b,scorer,transport,completion,np,pf,out);partial.append({'source_case_id':row['source_case_id'],'result':pin(out/'sources'/row['source_case_id']/'SOURCE_RESULT.json'),'status':result['status']})
                write(out/'progress'/f'{len(partial):02d}.json',{'completed_source_attempts':partial.copy()})
            transport.close();transport=None;receipt['status']='complete'
        except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc())
        finally:
            try:
                if transport:transport.abort()
            except BaseException as e:receipt.update(status='failed',cleanup_error=str(e))
            try:
                if completion is not None and original_class is not None:completion.CompletionProposal=original_class
                for sig,handler in handlers.items():signal.signal(sig,handler)
                collection=collect(c,out,scorer,pf);write(out/'FINAL_COLLECTION.json',collection)
                if not collection['stage_complete']:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                receipt.update(elapsed_seconds=time.monotonic()-start,parent_CPU_seconds=time.process_time()-cpu)
                receipt['allocated_CPU_seconds']=7*receipt['elapsed_seconds']
                try:receipt.update(actual_model_forwards=sum(read(p).get('actual_model_forwards',0) for p in out.rglob('mace_result.json')),worker_receipts=[pin(p) for p in out.glob('mace/lane0/CPU_WORKER_SUMMARY.json')])
                except BaseException as e:receipt.update(status='failed',cost_reporting_error=str(e))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);b=s.add_parser('bind');b.add_argument('--preparation',required=True);b.add_argument('--stage',choices=['scout','remaining'],required=True);b.add_argument('--output',required=True);b.add_argument('--scout-collection')
    r=s.add_parser('run');r.add_argument('--config',required=True);r.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='bind':bind(a.preparation,a.stage,a.output,a.scout_collection)
    else:sys.exit(run(a.config,a.config_sha))
