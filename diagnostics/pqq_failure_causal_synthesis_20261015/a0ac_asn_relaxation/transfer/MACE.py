#!/usr/bin/env python3
"""Finite population adapter for the unchanged canonical search and qualified GPU bridge."""
import argparse, fcntl, os, shutil, signal, time, traceback
from COMMON import *

PROTOCOL='uniform_terminal_amide_transfer_MACE_v1'
GPU_TEMPLATE=WORK/'asn_canonical_search_GPU_v1/CONFIG.json'

def dispatch_row(row):
    """Only request/directory identifiers change; physical source identities stay recorded."""
    r=copy.deepcopy(row)
    r.update(source_case_id=row['context_id'],case_id=row['context_id'],
             original_source_case_id=row['source_case_id'],original_case_id=row['case_id'])
    return r

def selected_ids(prep,stage):
    new=[r['context_id'] for r in prep['rows'] if r['status']=='prepared' and not r.get('reuse_only')]
    return new, new[:1] if stage=='scout' else new[1:] if stage=='remaining' else new

def bind(preparation,output,stage='scout',scout_collection=None):
    """Metadata only; actual geometry and canonical release must already exist."""
    pp=pin(preparation);p=read(pp);gm=read(p['manifest']);gate=canonical_gate(gm['canonical_gate']['collection'])
    require(gate==gm['canonical_gate'] and p['status']=='complete' and read(p['execution'])['status']=='complete','Actual canonical/preparation gates required')
    require(len(p['rows'])==257 and sum(len(r['memberships']) for r in p['rows'])==260,'Finite denominator differs')
    ids=[r['context_id'] for r in p['rows']];require(len(set(ids))==257,'Duplicate context dispatch IDs')
    cached=[r['context_id'] for r in p['rows'] if r.get('reuse_only')]
    require(len(cached)==34 and all(r['status']=='prepared' for r in p['rows'] if r.get('reuse_only')),'Exact34 cached contexts required')
    new,selected=selected_ids(p,stage);sp=pin(scout_collection) if scout_collection else None;se=None
    if stage=='remaining':
        require(sp is not None,'Actual first-context scout required');s=read(sp);sc=read(s['config']);se=pin(Path(sp['path']).with_name('EXECUTION.json'));e=read(se)
        require(s['protocol']==PROTOCOL and s['stage']=='scout' and s['selected_context_ids']==new[:1] and s['all_selected_common6_complete'],'Actual credible first-context pair required')
        require(e['status']=='complete' and e['config']==s['config'] and sc['preparation']==pp and sc['implementation']['sha256']==pin(__file__)['sha256'],'Scout execution/source differs')
    else:require(sp is None,'Only remaining stage reuses a scout')
    bp=pin(GPU_TEMPLATE);base=read(bp);template=read(base['worker_manifest'])
    require(gm['kernels']['search']==base['implementation'] and gm['kernels']['gpu_bridge']==base['execution_adapter'],'Qualified scientific/transport kernels differ')
    canonical=load('_transfer_bind_search',gm['kernels']['search'])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('MACE.py','MACE.sbatch','COMMON.py'):shutil.copyfile(HERE/name,out/name)
    rows=[r for r in p['rows'] if r['context_id'] in selected]
    tasks=[canonical.task_template(dispatch_row(r),t) for r in rows for t in r['task_templates']]
    worker=copy.deepcopy(template);worker.update(protocol_id=PROTOCOL,tasks=tasks,cases=[],fixtures=[],origin_task_count=len(tasks),preparation=pp,
       finite_scope={'declared_contexts':257,'declared_memberships':260,'cached_contexts':34,'new_admitted_contexts':len(new),'stage_contexts':len(rows),
                     'fixed_flip_calls':len(tasks),'optimizer_starts':len(tasks),'maximum_cross_calls':len(tasks),'maxiter':200,'new_SCC_calls':0,'new_DFT_calls':0})
    write(out/'worker_manifest.json',worker)
    gpu=read(base['GPU_manifest']);gpu.update(tasks=copy.deepcopy(tasks),source_worker_manifest=pin(out/'worker_manifest.json'),protocol_id=base['GPU_transport'])
    write(out/'gpu_manifest.json',gpu)
    c={k:base[k] for k in ('old_config','old_scorer','python','GPU_transport','GPU_qualification','GPU_qualification_execution','GPU_qualification_config','GPU_launcher','GPU_server','GPU_qualified_runtime')}
    c.update(protocol=PROTOCOL,implementation=pin(out/'MACE.py'),common=pin(out/'COMMON.py'),wrapper=pin(out/'MACE.sbatch'),plan=gm['implementation']['PLAN.md'],
       preparation=pp,preparation_manifest=p['manifest'],geometry_execution=p['execution'],inventory=gm['inventory'],ledger=gm['ledger'],cohorts=p['cohorts'],
       canonical_gate=gate,search_kernel=gm['kernels']['search'],gpu_adapter=gm['kernels']['gpu_bridge'],generic_contact_helper=gm['generic_contact_helper'],
       template_GPU_config=bp,template=base['worker_manifest'],worker_manifest=pin(out/'worker_manifest.json'),GPU_manifest=pin(out/'gpu_manifest.json'),GPU_plan=gm['implementation']['PLAN.md'],
       stage=stage,rows=rows,all_context_ids=ids,new_context_ids=new,cached_context_ids=cached,selected_context_ids=selected,
       scout_collection=sp,scout_execution=se,pool_orders=ORDERS,output=str(out/'execution'),declared_contexts=257,declared_memberships=260,
       new_SCC_calls=0,new_DFT_calls=0,classification=None,automatic_follow_on=False)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def install_admission(c,b,scorer,np):
    """Historical GPU cells retain their actual producing config; old CPU admission stays strict."""
    canonical=load('_transfer_admission_search',c['search_kernel']);gpu=load('_transfer_admission_GPU',c['gpu_adapter'])
    scorer.physical_result=canonical.physical_result
    old=getattr(scorer,'_transfer_CPU_admission',scorer.admitted_result);scorer._transfer_CPU_admission=old
    def admit(result,request,task,model,versions):
        if result.get('transport_kind')==gpu.KIND:return gpu.gpu_admission(result,request,task,model,versions,np,canonical,b)
        return old(result,request,task,model,versions)
    scorer.admitted_result=admit
    return old

def archive(cell,task,m,scorer,versions,np,canonical):
    rp=cell.get('MACE',cell.get('receipt'));r=read(rp)
    if r.get('request'):request=read(r['request']);observed=r['model']
    else:
        manifest=read(r['manifest']);request,=[t for t in manifest['tasks'] if t['task_id']==r['task_id']];observed=manifest['model']
    canonical.physical_result(r,request,task,m['model'],observed)
    require(r['energy_eV']==cell.get('MACE_eV',cell.get('energy_eV')),'Archived energy changed')
    for k in ('forces','saved_forces'):
        if cell.get(k):require(cell[k]==r['forces'],'Archived force pin differs')
    forces=np.load(check(r['forces']),allow_pickle=False)
    require(forces.shape==(task['atom_count'],3) and np.isfinite(forces).all(),'Archived force shape/values invalid')
    if r.get('transport_kind') or r['device']=='cpu':scorer.admitted_result(r,request,task,m['model'],versions)
    else:require(r['device']=='NVIDIA H200','Unexpected legacy native device')
    return {**cell,'status':'complete','receipt':rp,'xyz':task['xyz'],'energy_eV':r['energy_eV'],'forces':r['forces'],
       'request':r.get('request'),'historical_manifest':r.get('manifest'),'historical_device':r['device'],
       'charge':task['charge'],'multiplicity':1,'atom_count':task['atom_count'],'reused':True}

def validate_cached(row,c,m,b,scorer,versions,np,canonical,baseline):
    alias=row['executed_result_alias'];coll=read(row['producing_collection']);e=read(row['producing_execution']);pc=read(row['producing_config'])
    require(e['status']=='complete' and e['config']==coll['config']==row['producing_config'],'Cached producing execution differs')
    require(pc['preparation']==row['producing_preparation'],'Cached producing preparation differs')
    actual,=[r for r in coll['rows'] if r['source_case_id']==row['source_case_id']]
    require(actual==alias and alias['status']=='complete' and alias['prepared_mapping']==row['mapping'] and alias['admission']==row['admission'],'Cached exact source/result/admission differs')
    require(row.get('repeat_search_permitted') is False,'Cached context must not run again')
    for z in METALS:
        t,=[t for t in row['task_templates'] if t['metal']==z]
        for q in ORDER6:
            cell=alias['matrix'][z][q]
            require(cell['status']=='complete','Cached common6 incomplete')
            if q in ORDER3:require(cell['receipt']==baseline[z][q]['receipt'] and cell['xyz']==baseline[z][q]['xyz'],'Cached common3 aliases differ')
            else:archive(cell,{**t,'xyz':cell['xyz'],'atom_count':row['atom_count']},m,scorer,versions,np,canonical)
        require(alias['matrix'][z]['amide_flip']['xyz']['sha256']==t['xyz']['sha256'],'Cached fixed flip frame differs')
        proposal=read(alias['proposals'][z]);require(proposal['status']=='proposal_available' and proposal['optimizer']['success'] and proposal['final_geometry']['physical_feasible'],'Cached proposal not credible')
        for metal in METALS:
            cell=alias['matrix'][metal]['relaxed_flip_'+z]
            require(cell['proposal_receipt']==alias['proposals'][z] and read(cell['full_source_admission'])['pass'],'Cached full-source proposal admission differs')
    for q in ORDER6:
        a,both=[b.xyz(check(alias['matrix'][z][q]['xyz'])) for z in METALS]
        require(a[1:]==both[1:] and a[0][1:]==both[0][1:],'Cached paired geometry differs')
    return copy.deepcopy(alias)

def preflight(c,b,scorer,np):
    p=read(c['preparation']);gm=read(c['preparation_manifest']);m=read(c['worker_manifest']);template=read(c['template']);base=read(c['template_GPU_config'])
    require(pin(sys.executable)==c['python'] and p['manifest']==c['preparation_manifest'] and p['execution']==c['geometry_execution'],'Interpreter/preparation binding differs')
    require(p['status']=='complete' and read(c['geometry_execution'])['status']=='complete','Actual terminal geometry inventory required')
    require(canonical_gate(c['canonical_gate']['collection'])==c['canonical_gate']==gm['canonical_gate'],'Canonical release evidence changed')
    for k in ('implementation','common','wrapper','plan','search_kernel','gpu_adapter','generic_contact_helper','inventory','ledger'):check(c[k])
    require(c['search_kernel']==gm['kernels']['search']==base['implementation'] and c['gpu_adapter']==gm['kernels']['gpu_bridge']==base['execution_adapter'],'Qualified kernels differ')
    require(c['pool_orders']==ORDERS and c['cohorts']==p['cohorts']==gm['cohorts'],'Frozen pool/cohort identity changed')
    ids=[r['context_id'] for r in p['rows']];require(ids==c['all_context_ids'] and len(ids)==len(set(ids))==257,'Context denominator/order differs')
    memberships=[v for r in p['rows'] for v in r['memberships']]
    cohorts=cohort_ids(read(c['ledger']))
    require(len(memberships)==len(set(memberships))==260 and set(memberships)==set(cohorts['expanded135'])|set(cohorts['current_ca125']),'Membership denominator differs')
    require([(r['context_id'],r['memberships'],r['scope']) for r in p['rows']]==[(s['context_id'],s['memberships'],s['scope']) for s in gm['sources']],'Exact source/context binding changed')
    new,selected=selected_ids(p,c['stage']);cached=[r['context_id'] for r in p['rows'] if r.get('reuse_only')]
    require(new==c['new_context_ids'] and selected==c['selected_context_ids'] and cached==c['cached_context_ids'] and len(cached)==34,'Finite dispatch/reuse selection differs')
    require(c['rows']==[r for r in p['rows'] if r['context_id'] in selected],'Selected preparation row changed')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        require(m[k]==template[k],'Frozen engine setting differs: '+k)
    require(gm['settings']==m['settings'] and gm['source_modules']==m['source_modules'] and gm['source_snapshot']==m['source_snapshot'],'Preparation/search kernels differ')
    for q in m['source_modules'].values():check(q)
    for k in ('implementation','calculator_source','model_executable','optimizer_executable'):check(m[k])
    check(m['model']['checkpoint']);canonical=load('_transfer_search_kernel',c['search_kernel'])
    require(m['tasks']==[canonical.task_template(dispatch_row(r),t) for r in c['rows'] for t in r['task_templates']] and m['origin_task_count']==len(m['tasks']) and m['threads_per_model']==7 and m['lanes']==1 and not m['cases'] and not m['fixtures'],'Worker dispatch differs')
    require(read(c['GPU_manifest'])['tasks']==m['tasks'] and c['GPU_transport']==base['GPU_transport'],'GPU authorization differs')
    old=read(c['old_config']);qualification=read(old['qualification']);summary=read(old['qualified_worker_summary'])
    require(qualification['pass'] and summary['threads']==7 and summary['manifest']==old['template'],'Actual CPU reference qualification differs')
    for q in qualification['checks']:
        r=read(q['receipt']);check(r['forces']);read(r['request']);require(q['pass'] and r['model']==m['model'] and r['versions']==summary['versions'],'CPU reference fixture differs')
    versions=summary['versions'];install_admission(c,b,scorer,np)
    import slsqp_precision, scipy
    import scipy.optimize._slsqp_py as sw
    import scipy.optimize._slsqplib as sk
    from mace_site_kinematics import Kinematics
    completion=slsqp_precision.engine()[0].completion
    require(completion.SETTINGS==m['settings'] and scipy.__version__==m['optimizer_software']['version'] and pin(sw.__file__)==m['optimizer_software']['wrapper'] and pin(sk.__file__)==m['optimizer_software']['kernel'],'Actual optimizer differs')
    for module in list(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename and Path(filename).name in m['source_modules']:require(pin(filename)==m['source_modules'][Path(filename).name],'Actually imported scientific source differs')
    baseline={};reused={};checks=[]
    for row in p['rows']:
        cid=row['context_id'];s=row['scope'];n=row['atom_count'];st=read(s['state_signature']);search=read(s['search_manifest'])
        require(search['settings']==m['settings'] and [t for t in search['tasks'] if t['case_id']==s['case_id']]==s['original_tasks'],'Original executed task/settings differ')
        require(st['endpoint_charges']==row['endpoint_charges'] and st['multiplicities']=={'Ca':1,'La':1} and st['water_inventory']==row['water_inventory'],'Source state differs')
        require(all(row[k]==s[k] for k in ('source_case_id','case_id','record_id','selection_id','protein')),'Original context identity changed')
        for z in METALS:check(s['maps'][z])
        check(s['source_preparation']);check(s['protonated_source']);baseline[cid]={z:{} for z in METALS}
        for z in METALS:
            for q in ORDER3:
                cell=s['six_saved_cells'][z][q];t={'xyz':cell['xyz'],'charge':row['endpoint_charges'][z],'multiplicity':1,'atom_count':n}
                atoms=b.xyz(check(t['xyz']));require(len(atoms)==n and atoms[0][0]==z,'Archived count/metal differs')
                baseline[cid][z][q]=archive(cell,t,m,scorer,versions,np,canonical)
        for q in ORDER3:
            a,other=[b.xyz(check(baseline[cid][z][q]['xyz'])) for z in METALS]
            require(a[1:]==other[1:] and a[0][1:]==other[0][1:],'Original paired common pool differs')
        if row['status']!='prepared':
            require(not row.get('reuse_only') and bool(row.get('reason')),'Unavailable preparation lacks reason');checks.append({'context_id':cid,'baseline_pass':True,'geometry':'unavailable','reason':row['reason']});continue
        admission=read(row['admission']);require(admission['pass'],'Rebased admission failed')
        if not row.get('reuse_only') or 'q0_geometry_checks' in admission:
            require(all(v['physical_feasible'] for v in admission['q0_geometry_checks'].values()),'Rebased q0 admission failed')
        require(all(row['original_row'][k]==s[k] for k in ('maps','state_signature','source_preparation','protonated_source')),'Prepared original source/state differs')
        kin=Kinematics(read(row['mapping'])['context']);zero=np.zeros(len(kin.modes));require(kin.check(zero)['pass'],'Existing rebased geometry/Jacobian gate failed')
        check(row['contact']['origin']);read(row['contact']['identities']);require(row['contact']['source']==s['protonated_source'],'Full-source contact identity differs')
        templates={t['metal']:t for t in row['task_templates']};require(set(templates)==set(METALS) and len(row['task_templates'])==2,'Prepared task pair differs')
        for z,t in templates.items():
            ot,=[x for x in s['original_tasks'] if x['metal']==z]
            require(t['mapping']==row['mapping'] and t['contact']==row['contact'] and (t.get('settings',m['settings']) if row.get('reuse_only') else t['settings'])==m['settings'] and t['atom_count']==n and t['charge']==row['endpoint_charges'][z] and t['multiplicity']==1,'Task state/geometry/settings differs')
            for k in ('active_indices','active_mode_ids','active_roles','selector','mode_count','selection_id'):require(t[k]==ot[k],'Original selector/mode changed: '+k)
            a=b.xyz(check(t['xyz']));require(len(a)==n and a[0][0]==z and np.max(np.abs(kin.evaluate(zero)[1]-np.asarray([v[1:] for v in a])))<=1e-12,'Rebased q0 differs')
            require(completion.angular.final_geometry(kin,t,np.zeros(4),[v[0] for v in a])['physical_feasible'],'Original q0 physical check failed')
        a,other=[b.xyz(check(templates[z]['xyz'])) for z in METALS];require(a[1:]==other[1:] and a[0][1:]==other[0][1:],'Prepared fixed-flip pair differs')
        if row.get('reuse_only'):reused[cid]=validate_cached(row,c,m,b,scorer,versions,np,canonical,baseline[cid])
        checks.append({'context_id':cid,'baseline_pass':True,'geometry':'passed','reuse_only':bool(row.get('reuse_only'))})
    require(len(baseline)==257 and len(reused)==34,'All original/cached context denominator differs')
    prior={}
    if c['stage']=='remaining':
        scoll=read(c['scout_collection']);se=read(c['scout_execution']);sc=read(scoll['config'])
        require(se['status']=='complete' and se['config']==scoll['config'] and scoll['all_selected_common6_complete'] and scoll['selected_context_ids']==new[:1],'Actual first-context scout failed')
        for k in ('preparation','preparation_manifest','search_kernel','gpu_adapter','canonical_gate','plan'):require(sc[k]==c[k],'Scout protocol/input differs: '+k)
        require(all(sc[k]['sha256']==c[k]['sha256'] for k in ('implementation','common','wrapper')),'Scout wrapper changed')
        r,=[r for r in scoll['rows'] if r['context_id']==new[0]];pr=next(r for r in p['rows'] if r['context_id']==new[0])
        require(r['status']=='complete' and r['prepared_mapping']==pr['mapping'] and r['admission']==pr['admission'],'Scout final context differs')
        for z in METALS:
            t=next(t for t in pr['task_templates'] if t['metal']==z)
            for q in ORDER6:archive(r['matrix'][z][q],{**t,'xyz':r['matrix'][z][q]['xyz']},m,scorer,versions,np,canonical)
            prop=read(r['proposals'][z]);require(prop['status']=='proposal_available' and prop['optimizer']['success'] and prop['final_geometry']['physical_feasible'],'Scout proposal not credible')
        prior[new[0]]=r
    return {'status':'passed','versions':versions,'baseline':baseline,'cached':reused,'scout':prior,'source_checks':checks,
            'verified_original_cells':257*6,'verified_cached_contexts':34,'new_model_calls':0},completion

def normalized_result(raw,row):
    r=copy.deepcopy(raw)
    for k in ('source_case_id','case_id','record_id','selection_id','protein','atom_count','endpoint_charges','terminal_mode_selected_by_metal','original_four_modes','context_id','memberships'):
        r[k]=row[k]
    r['runtime_context_id']=row['context_id'];r['preparation_status']=row['status'];r['preparation_reason']=row.get('reason')
    return r

def collect(c,out,scorer,pf):
    p=read(c['preparation']);rows=[];attempted=[]
    for row in p['rows']:
        cid=row['context_id'];path=Path(out)/'sources'/cid/'SOURCE_RESULT.json';origin='not_executed'
        if path.exists():
            raw=read(path);require(raw['source_case_id']==cid and raw['case_id']==cid,'Runtime dispatch identity differs');origin='new_stage';attempted.append(cid)
        elif pf and cid in pf['cached']:raw=pf['cached'][cid];origin='exact_executed_reuse'
        elif pf and cid in pf['scout']:raw=pf['scout'][cid];origin='completed_scout_reuse'
        else:
            matrix=copy.deepcopy(pf['baseline'][cid]) if pf else {z:{q:{'status':'unavailable','energy_eV':None,'reason':'Numerical preflight unavailable','saved_input':row['scope']['six_saved_cells'][z][q]} for q in ORDER3} for z in METALS}
            raw={'status':'unavailable' if row['status']!='prepared' else 'not_executed','reason':row.get('reason') or 'Declared stage not executed','matrix':matrix,'proposals':{}}
        r=normalized_result(raw,row);r['result_origin']=origin
        if path.exists():r['source_result']=pin(path)
        if row.get('reuse_only'):r['reuse_provenance']={k:row[k] for k in ('producing_config','producing_collection','producing_execution','producing_preparation')}
        for z in METALS:
            for q in ORDER6:r['matrix'][z].setdefault(q,{'status':'unavailable','energy_eV':None,'reason':r.get('reason') or 'Not executed'})
        energies={z:{q:v.get('energy_eV') if v['status']=='complete' else None for q,v in r['matrix'][z].items()} for z in METALS}
        r['pools']={label:scorer.choose(energies,order) for label,order in ORDERS.items()};r['changes']={}
        for a,e in (('common3','common4'),('common4','common6'),('common3','common6')):
            r['changes'][a+'_to_'+e]={}
            for policy in ('mathematical','operational'):
                before,after=r['pools'][a],r['pools'][e];change=None
                if before and after:
                    dz={z:after[policy]['per_metal_work_kcal_mol'][z]-before[policy]['per_metal_work_kcal_mol'][z] for z in METALS}
                    change={'per_metal_change_kcal_mol':dz,'paired_change_kcal_mol':dz['Ca']-dz['La']}
                r['changes'][a+'_to_'+e][policy]=change
        r.pop('common4_to6',None);rows.append(r)
    stage_complete=bool(pf) and all(cid in attempted for cid in c['selected_context_ids'])
    terminal=bool(pf) and all(r['status']!='not_executed' for r in rows)
    return {'protocol':PROTOCOL,'config':pin(Path(c['output']).parent/'CONFIG.json'),'preparation':c['preparation'],'preparation_manifest':c['preparation_manifest'],
        'stage':c['stage'],'stage_complete':stage_complete,'all_contexts_terminal':terminal,'selected_context_ids':c['selected_context_ids'],
        'all_selected_common6_complete':bool(pf) and all(r['status']=='complete' for r in rows if r['context_id'] in c['selected_context_ids']),
        'all_admitted_common6_complete':bool(pf) and all(r['status']=='complete' for r in rows if r['preparation_status']=='prepared'),
        'declared_contexts':257,'declared_memberships':260,'cached_contexts':34,'all_context_ids':c['all_context_ids'],'pool_orders':ORDERS,'cohorts':c['cohorts'],
        'canonical_gate':c['canonical_gate'],'GPU_transport':c['GPU_transport'],'rows':rows,'classification':None,
        'new_SCC_calls':0,'new_DFT_calls':0,'automatic_follow_on':False,'production_changed':False,'actual_model_forwards':None}

def run(path,sha,action='run'):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['protocol']==PROTOCOL and c['implementation']==pin(__file__) and c['common']==pin(Path(__file__).with_name('COMMON.py')),'Frozen adapter differs')
    out=Path(c['output']) if action=='run' else Path(c['output']).parent/('preflight_'+os.environ['SLURM_JOB_ID'])
    start=time.monotonic();cpu=time.process_time();transport=None;completion=None;original=None;handlers={};pf=None;scorer=None;old_admission=None
    allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_JOB_PARTITION','SLURM_JOB_GPUS','SLURM_GPUS_ON_NODE','SLURM_CPUS_ON_NODE','CUDA_VISIBLE_DEVICES')}
    receipt={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'allocation':allocation,'stage':c['stage'],'action':action,'actual_model_forwards':None,'new_SCC_calls':0,'new_DFT_calls':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            m=read(c['worker_manifest']);b=load('_transfer_batch',m['implementation']);canonical=load('_transfer_search',c['search_kernel']);gpu=load('_transfer_GPU',c['gpu_adapter'])
            scorer=load('_transfer_scorer',c['old_scorer']);old_admission=scorer.admitted_result;b.parent_env();b.modules(m)
            import numpy as np
            pf,completion=preflight(c,b,scorer,np);write(out/'PREFLIGHT.json',pf)
            if action=='run' and c['rows']:
                require(os.environ.get('SLURM_JOB_NUM_NODES')=='1' and os.environ.get('SLURM_GPUS_ON_NODE')=='1' and bool(os.environ.get('SLURM_JOB_GPUS')),'Actual one-node/one-GPU allocation required')
                require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles') and len(os.environ['CUDA_VISIBLE_DEVICES'].split(','))==1 and int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])>=1,'Actual GPU visibility/CPU grant required')
                launcher=load('_transfer_GPU_launcher',c['GPU_launcher']);runtime=gpu.runtime(c,cp,launcher,out)
                helper=load('_transfer_contact',c['generic_contact_helper']);original,guarded=canonical.guarded_class(completion,helper,np);completion.CompletionProposal=guarded
                native=load('_transfer_native_GPU',c['GPU_server']);handlers=b.handlers();transport=gpu.Bridge(c,cp,m,b,native,launcher,runtime)
                receipt.update(GPU_command=transport.gpu.command,worker_cgroup=transport.gpu.worker_cgroup)
                progress=[]
                for row in c['rows']:
                    raw=canonical.source_run(c,dispatch_row(row),m,b,scorer,transport,completion,np,pf,out)
                    result=normalized_result(raw,row);write(out/'sources'/row['context_id']/'CONTEXT_RESULT.json',result)
                    progress.append({'context_id':row['context_id'],'source_case_id':row['source_case_id'],'memberships':row['memberships'],'status':raw['status'],
                                     'result':pin(out/'sources'/row['context_id']/'CONTEXT_RESULT.json')})
                    write(out/'progress'/f'{len(progress):03d}.json',{'completed_context_attempts':progress.copy()})
                calls=len(transport.calls);transport.close();transport=None
                sp=Path(c['GPU_manifest']['path']).parent/('gpu_summary_'+os.environ['SLURM_JOB_ID']+'.json');summary=read(sp);receipt['native_GPU_summary']=pin(sp)
                require(summary['manifest']==c['GPU_manifest'] and summary['new_MACE_calls']==calls,'Native request accounting differs')
                receipt.update(native_model_requests=summary['new_MACE_calls'],native_failed_requests=summary['failed_requests'])
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
        finally:
            try:
                if transport:transport.abort()
            except BaseException as e:receipt.update(status='failed',cleanup_error=str(e))
            try:
                if completion is not None and original is not None:completion.CompletionProposal=original
                if scorer is not None and old_admission is not None:scorer.admitted_result=old_admission
                for sig,h in handlers.items():signal.signal(sig,h)
                if action=='run' and scorer is not None:
                    result=collect(c,out,scorer,pf);write(out/'FINAL_COLLECTION.json',result)
                    if not result['stage_complete']:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,parent_CPU_seconds=time.process_time()-cpu,
                    allocated_CPU_seconds=int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])*elapsed,
                    allocated_GPU_seconds=elapsed if os.environ.get('SLURM_GPUS_ON_NODE')=='1' else 0)
                if action=='preflight':receipt['new_model_calls']=0
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('bind');a.add_argument('--preparation',required=True);a.add_argument('--output',required=True);a.add_argument('--stage',choices=['scout','remaining','all'],default='scout');a.add_argument('--scout-collection')
    for op in ('run','preflight'):
        a=s.add_parser(op);a.add_argument('--config',required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.op=='bind':bind(a.preparation,a.output,a.stage,a.scout_collection)
    else:sys.exit(run(a.config,a.config_sha,a.op))
