#!/usr/bin/env python3
"""Conditional canonical25 join around the executed fixed-charge solvent kernel."""
import argparse, copy, fcntl, hashlib, importlib.util, json, os, shutil, sys, time, traceback
from pathlib import Path
from types import FunctionType
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_solvent'
OLD=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_v1'
FAILED=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_v1'
LEGACY=ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2'
REUSE='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3'
PROTOCOL='canonical25_terminal_amide_frozen_ALPB_common6_v1'
GPU_KIND='qualified_GPU_request_bridge_v1'
METALS=('Ca','La');ORDER3=('origin','adaptive_Ca','adaptive_La');ORDER4=ORDER3+('amide_flip',);ORDER6=ORDER4+('relaxed_flip_Ca','relaxed_flip_La')
ORDERS={'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)}
KEYS=('matrix','runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model','corrected_ORCA_parameter_provenance','ORCA_correction_policy','solvent_adapter','pool_adapter','scorer','numerical_python')
def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
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

def gate(sc,sr,ex,cp):
    require(sc['protocol']==sr['protocol']=='canonical25_terminal_amide_common6_v1' and sc['stage']==sr['stage']=='remaining','Completed full canonical stage required')
    require(sr['config']==ex['config']==cp and ex['status']=='complete' and sr['stage_complete'] and sr['all25_common6_complete'],'Actual complete execution and all25 common6 required')
    require(sr['declared_sources']==25 and [r['source_case_id'] for r in sr['rows']]==sc['all_source_ids'] and sr['pool_orders']==sc['pool_orders']==ORDERS,'Canonical scope/order differs')
    if sc.get('GPU_transport'):
        require(sc['GPU_transport']==sr.get('GPU_transport')==GPU_KIND,'Unqualified GPU collection marker')
        check(sc['execution_adapter']);check(sc['GPU_qualification']);check(sc['GPU_qualification_execution']);check(ex['native_GPU_summary'])
        require(read(sc['GPU_qualification'])['all_four_pass'] and read(sc['GPU_qualification_execution'])['status']=='complete','Actual GPU parity qualification missing')
        require(ex['native_failed_requests']==0,'GPU native request failure')
    else:require('GPU_transport' not in sr,'GPU collection cannot use CPU-only configuration')
    for row in sr['rows']:
        require(row['status']=='complete' and set(row['proposals'])==set(METALS),'Incomplete canonical proposal pair')
        for z in METALS:
            p=read(row['proposals'][z]);require(p['status']=='proposal_available' and p['optimizer']['success'] and p['final_geometry']['physical_feasible'],'Noncredible final proposal')
            for q in ORDER6:require(row['matrix'][z][q]['status']=='complete','Incomplete canonical MACE cell')

def descriptor(sid,z,q,cell,native_pin):
    nr=read(native_pin);t=nr['input']
    require(nr['status']=='available','Cached native endpoint unavailable')
    original={**cell,**{k:t[k] for k in ('charge','multiplicity','atom_count','valence_electrons')},'native_frozen_ALPB_result':native_pin,'native_frozen_ALPB_hartree':nr['G_qvac_hartree']}
    require(t['xyz']==original['xyz'] and any(u['MACE']==original['MACE'] and u['MACE_eV']==original['MACE_eV'] for u in t['uses']),'Cached native/MACE provenance differs')
    return {'source_case_id':sid,'metal':z,'candidate':q,'original':original,**{k:t[k] for k in ('maps','state_signature','source_preparation')}}

def identity(t):
    v=t.get('original',t)
    return (t['source_case_id'],t['metal'],v['xyz']['sha256'],v['charge'],v['multiplicity'],v['atom_count'],t['state_signature']['sha256'],t['source_preparation']['sha256'])

class LegacyComponentFacade:
    """Only the pinned pre-reuse task schema; no shared module or task mutation."""
    def __init__(self,c,r):
        self.r=r;self.pins=c['legacy_component'];self.manifest=read(self.pins['manifest']);self.saved=read(self.pins['collection']);self.base=check(self.pins['manifest']).parent
        require(self.manifest['driver']==self.pins['driver'] and self.saved['manifest']==self.pins['manifest'],'Legacy driver/collection lineage differs')
        require(self.manifest['runtime']==c['runtime'] and self.manifest['implementation']==c['factorial'],'Legacy runtime/implementation differs')
        require(all(self.manifest['settings'][k]==c['settings'][k] for k in r.SCIENTIFIC_SETTINGS),'Legacy scientific settings differ')
        driver=load('_canonical_legacy_component_reader',self.pins['driver']);require(driver.BASE==self.base,'Legacy reader base differs')
        self.raw=driver.raw_results(self.manifest);self.evidence=[]
        # Reuse the exact component function with a private raw resolver. Its
        # module globals, original raw tasks and archive files remain unchanged.
        self.kernel=FunctionType(r.component.__code__,dict(r.component.__globals__,raw=self.legacy_raw),r.component.__name__,r.component.__defaults__,r.component.__closure__)

    def legacy_raw(self,base,task):
        require(Path(base).resolve()==self.base and 'reuse' not in task,'Unrecognized legacy routing request')
        actual,result,p=self.raw.get(task['id'],(None,None,None))
        require(actual==task and result is not None and p is not None,'Legacy task is not the exact archived task')
        require(read(p)==result and result['task']==task,'Legacy RESULT/task identity differs')
        require(result['manifest']==self.pins['manifest'],'Legacy endpoint manifest differs')
        require(p==pin(self.base/'points'/task['id']/'RESULT.json'),'Legacy endpoint path differs')
        self.r.verify_raw(result,task,self.manifest|{'_base':str(self.base)})
        return result,p

    def component(self,base,v,w):
        if 'reuse' in v and 'reuse' in w:return self.r.component(base,v,w)
        require('reuse' not in v and 'reuse' not in w and Path(base).resolve()==self.base,'Unsupported mixed/unknown component schema')
        out=self.kernel(base,v,w);vr,vp=self.legacy_raw(base,v);wr,wp=self.legacy_raw(base,w)
        # Exact additional criterion from the pinned COMPONENT_V2.collect.
        ref=max(abs(point.get('reference_energy_difference_hartree',0.0)) for point in (vr,wr))
        out['accepted']=out['accepted'] and ref<=1e-6
        out['maximum_available_archived_primary_error_hartree']=ref
        if not out['accepted']:out['primary_solvation_hartree']=None;out['reason']='component_or_population_or_archived_parity_failed'
        matches=[(i,row) for i,row in enumerate(self.saved['rows']) if row['inputs']=={'vacuum':vp,'water':wp}]
        require(len(matches)==1,'Legacy component row join differs');index,stored=matches[0]
        for key in ('accepted','primary_vacuum_energy_hartree','primary_water_energy_hartree','primary_solvation_hartree','component_error_hartree','absolute_energy_errors_hartree','maximum_available_archived_primary_error_hartree'):
            require(out[key]==stored[key],'Legacy component replay differs: '+key)
        require(out['maximum_population_error_e']==stored['max_population_error_e'],'Legacy population criterion differs')
        out['legacy_component_receipt']={'manifest':self.pins['manifest'],'driver':self.pins['driver'],'collection':self.pins['collection'],'row_index':index}
        self.evidence.append({'vacuum_task':v['id'],'water_task':w['id'],'collection_row_index':index,'accepted':out['accepted'],'exact_archived_component_replay':True})
        return out

    def negative_metadata_checks(self):
        row=next(row for row in self.saved['rows'] if row['inputs']['vacuum'] is not None and read(row['inputs']['vacuum'])['manifest']==self.pins['manifest'])
        v=read(row['inputs']['vacuum'])['task'];w=read(row['inputs']['water'])['task'];checks=[]
        for label,alter in [('unknown_task_id',{'id':v['id']+'__invalid_metadata'}),('changed_xyz_pin',{'xyz':dict(v['xyz'],sha256='0'*64)})]:
            bad=copy.deepcopy(v);bad.update(alter)
            try:self.component(self.base,bad,w)
            except (AssertionError,ValueError,KeyError) as e:checks.append({'case':label,'rejected':True,'reason':str(e),'real_task':v['id']})
            else:raise ValueError('Corrupted legacy metadata was accepted: '+label)
        return checks

def bind(score_config,score_collection,output):
    """Metadata only; called by root after its actual full-MACE terminal review."""
    sc=read(score_config);sr=read(score_collection);expath=Path(score_collection).with_name('EXECUTION.json');gate(sc,sr,read(expath),pin(score_config))
    prep=read(sc['preparation']);gm=read(sc['preparation_manifest']);cc=read(sc['canonical_collection']);inv=read(cc['inventory']);nc=read(cc['canonical_native_collection'])
    require([r['record_id'] for r in sr['rows']]==inv['calibration_rule']['designated_record_ids']==[r['record_id'] for r in cc['rows']],'Exact designated25 order required')
    prior=read(OLD/'CONFIG.json');saved=read(OLD/'execution/FINAL_COLLECTION.json')
    require(saved['config']==pin(OLD/'CONFIG.json') and saved['available_new_cells']==40 and read(OLD/'execution/EXECUTION.json')['status']=='complete','Actual A0AC solvent reuse unavailable')
    cached=[];records={r['record_id']:r for r in inv['records']};origins={}
    for row in sr['rows']:
        source=records[row['record_id']];sid=row['source_case_id'];require(source['source_case_id']==sid,'Canonical source join differs')
        origins[sid]={}
        for z in METALS:
            for q in ORDER3:
                v=source['matrix'][z][q];n=nc['tasks'][v['native_task_id']];d=descriptor(sid,z,q,v,n['result']);cached.append(d)
                if q=='origin':origins[sid][z]=d['original']
    alias=next(r for r in saved['rows'] if r['source_case_id']==REUSE)
    for z in METALS:
        for q in ORDER6[3:]:
            e=alias['evidence'][z][q];cached.append(descriptor(REUSE,z,q,{'xyz':e['xyz'],'MACE':e['MACE'],'MACE_eV':e['MACE_eV']},e['native_result']))
    require(len(cached)==156,'Exact150 canonical plus6 A0AC cached cells required')
    cache_keys={}
    for i,d in enumerate(cached):cache_keys.setdefault(identity(d),i)
    cells=[];targets=[];new_keys={};prepared={r['source_case_id']:r for r in prep['rows']}
    for row in sr['rows']:
        sid=row['source_case_id'];pr=prepared[sid];source=records[row['record_id']]
        require(row['prepared_mapping']==pr['mapping'] and read(pr['admission'])['pass'],'Final prepared mapping differs')
        for z in METALS:
            for q in ORDER6:
                v=row['matrix'][z][q];base=origins[sid][z]
                t={'task_id':sid+'__canonical__'+q+'__'+z,'source_case_id':sid,'case_id':row['case_id'],'candidate':q,'metal':z,
                   'xyz':v['xyz'],'MACE':v['receipt'],'MACE_eV':v['energy_eV'],**{k:base[k] for k in ('charge','multiplicity','atom_count','valence_electrons')},
                   'maps':source['maps'] if q in ORDER3 else {metal:pr['mapping'] for metal in METALS},'state_signature':source['state_signature'],'source_preparation':source['source_preparation'],
                   'prepared_admission':pr['admission'],'proposal_receipt':v.get('proposal_receipt'),'full_source_admission':v.get('full_source_admission')}
                for k in ('xyz','MACE','state_signature','source_preparation'):check(t[k])
                key=identity(t)
                if key in cache_keys:route={'kind':'cached','cache_index':cache_keys[key]}
                elif key in new_keys:route={'kind':'new','task_id':new_keys[key]}
                else:
                    require(q not in ORDER3 and sid!=REUSE,'Required original reuse missing')
                    cells.append(t);new_keys[key]=t['task_id'];route={'kind':'new','task_id':t['task_id']}
                targets.append({'cell':t,'route':route})
    require(len(targets)==300 and len(cells)<=144,'Finite endpoint scope differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('SOLVENT.py','RUN.sbatch','PLAN.md','RECOVERY_V2.md'):shutil.copy2(HERE/name,out/name)
    failed=read(FAILED/'preflight_5135/EXECUTION.json');require(failed['status']=='failed' and failed['reason']=="'reuse'" and all(failed[k]==0 for k in ('SCC_started','SCC_terminal','native_component_calls')),'Original zero-call schema failure differs')
    legacy=read(LEGACY/'manifest.json');legacy_pins={'manifest':pin(LEGACY/'manifest.json'),'driver':legacy['driver'],'collection':pin(LEGACY/'FINAL_COLLECTION.json')}
    for p in legacy_pins.values():check(p)
    c={k:prior[k] for k in KEYS};c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),wrapper=pin(out/'RUN.sbatch'),plan=pin(out/'PLAN.md'),
       template=pin(OLD/'CONFIG.json'),dispatcher=prior['prior_adapter'],completed_A0AC_collection=pin(OLD/'execution/FINAL_COLLECTION.json'),completed_A0AC_execution=pin(OLD/'execution/EXECUTION.json'),
       score_config=pin(score_config),score_collection=pin(score_collection),score_execution=pin(expath),score_adapter=sc['implementation'],mace_scorer=sc['old_scorer'],
       GPU_admission_adapter=sc.get('execution_adapter') if sc.get('GPU_transport')==GPU_KIND else None,
       preparation=sc['preparation'],preparation_manifest=sc['preparation_manifest'],canonical_collection=sc['canonical_collection'],calibration_inventory=cc['inventory'],canonical_native_collection=cc['canonical_native_collection'],
       calibration_plan=sc['calibration_plan'],search_plan=sc['plan'],generic_contact_helper=sc['generic_contact_helper'],
       cells=cells,targets=targets,cached=cached,rows=[{k:r[k] for k in ('source_case_id','case_id','record_id','protein','expected_class','biological_group','source_conditioning_metal')} for r in sr['rows']],
       protein='canonical25',CPUs=7,declared_sources=25,declared_added_matrix_cells=150,declared_new_unique_cells=len(cells),reused_original_cells=150,reused_A0AC_additional_cells=6,
       maximum_SCC_calls=4*len(cells),maximum_native_component_calls=2*len(cells),pool_orders=ORDERS,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0,automatic_submission=False,
       legacy_component=legacy_pins,schema_recovery_note=pin(out/'RECOVERY_V2.md'),original_failed_config=pin(FAILED/'CONFIG.json'),original_failed_preflight=pin(FAILED/'preflight_5135/EXECUTION.json'))
    old_config=read(c['original_failed_config'])
    require(failed['config']==c['original_failed_config'],'Original failure config differs')
    for k,value in old_config.items():
        if k not in ('implementation','wrapper','plan'):require(c[k]==value,'Schema repair altered scientific configuration: '+k)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,f,r,n,a,h,selector,search,mace):
    require(pin(sys.executable)==c['numerical_python'] and c['protocol_id']==PROTOCOL,'Interpreter/protocol differs')
    for k in ('plan','wrapper','search_plan','calibration_plan','generic_contact_helper','completed_A0AC_execution','schema_recovery_note','original_failed_config','original_failed_preflight'):check(c[k])
    prior=read(c['template'])
    for k in KEYS:require(c[k]==prior[k],'Qualified solvent dependency differs: '+k)
    require(c['dispatcher']==prior['prior_adapter'],'Executed solvent dispatcher differs')
    require(read(c['ORCA_correction_policy'])['candidate_sha256']==c['corrected_ORCA_parameter_provenance']['sha256'],'Full La parameter lineage differs');check(c['corrected_ORCA_parameter_provenance'])
    sc=read(c['score_config']);sr=read(c['score_collection']);gate(sc,sr,read(c['score_execution']),c['score_config'])
    require(c['score_adapter']==sc['implementation'] and c['mace_scorer']==sc['old_scorer'] and c['preparation']==sc['preparation'] and c['pool_orders']==ORDERS,'Canonical MACE protocol join differs')
    worker=read(sc['worker_manifest']);b=load('_canonical_solvent_worker_metadata',worker['implementation']);b.parent_env();b.modules(worker)
    import numpy as np
    from mace_site_kinematics import Kinematics
    mace.physical_result=search.physical_result
    if sc.get('GPU_transport'):
        require(c['GPU_admission_adapter']==sc['execution_adapter'],'Qualified GPU admission source differs')
        gpu=load('_canonical_solvent_GPU_admission',c['GPU_admission_adapter']);cpu_admission=mace.admitted_result
        require(gpu.KIND==GPU_KIND,'Unknown GPU admission contract')
        base=read(sc['base_CPU_config'])
        for key,value in base.items():
            if key not in ('output','wrapper'):require(sc[key]==value,'GPU bridge altered canonical scientific field: '+key)
        def admit(value,request,task,model,versions):
            if value.get('transport_kind')==GPU_KIND:
                require(value['adapter_config']==c['score_config'],'GPU receipt belongs to another execution')
                return gpu.gpu_admission(value,request,task,model,versions,np,search,b)
            return cpu_admission(value,request,task,model,versions)
        mace.admitted_result=admit
    else:require(c['GPU_admission_adapter'] is None,'Unexpected GPU admission adapter')
    pf,completion=search.preflight(sc,b,mace,np)
    replay=search.collect(sc,Path(c['score_collection']['path']).parent,mace,pf)
    if sc.get('GPU_transport'):replay['GPU_transport']=GPU_KIND
    require(replay==sr,'Actual full MACE common6 replay differs')
    require(worker['model']==c['MACE_model'],'MACE model lineage differs')
    runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for k in ('xtb','xtb_parameter_file','tblite_python'):check(runtime[k])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:check(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:check(p)
    for p in list(nm['implementation'].values())+nm['build']['pins']+build['linked_libraries']:check(p)
    check(build['binary']);require(build['manifest']==c['native_manifest'] and build['binary']==c['native_binary'] and nm['MACE_model']==c['MACE_model'],'Native build/model differs')
    prep=read(c['preparation']);prepared={p['source_case_id']:p for p in prep['rows']};helper=load('_canonical_solvent_contact',c['generic_contact_helper']);mr={};physical=[]
    for row in sr['rows']:
        sid=row['source_case_id'];pr=prepared[sid];contact=pr['contact'];kin=Kinematics(read(pr['mapping'])['context'])
        with np.load(check(contact['origin']),allow_pickle=False) as data:union=data['positions'].copy();symbols=data['symbols'].tolist()
        ids=read(contact['identities'])['centers']
        for q in ORDER6[3:]:
            coords=None
            if q!='amide_flip':
                owner=q.removeprefix('relaxed_flip_');proposal=read(row['proposals'][owner]);point=proposal['proposal'];pm=read(proposal['manifest']);task,=[t for t in pm['tasks'] if t['task_id']==proposal['task_id']]
                require(task['mapping']==pr['mapping'] and task['contact']==contact and pm['settings']==worker['settings'],'Final task/map/settings differ')
                original=next(t for t in pr['task_templates'] if t['metal']==owner)
                for k in ('active_indices','active_mode_ids','selector','mode_count'):require(task[k]==original[k],'Final original modes differ: '+k)
                full=completion.angular.full_q(task,np.asarray(point['active_q_radian']));require(np.array_equal(full,np.asarray(point['full_q'])),'Final q differs')
                coords=kin.evaluate(full)[1];fg=completion.angular.final_geometry(kin,task,np.asarray(point['active_q_radian']),[v[0] for v in b.xyz(check(point['coordinate']))]);require(fg['physical_feasible'],'Final physical replay failed')
                changed=union.copy();changed[:pr['atom_count']]=coords;cg=helper.contact_check(np,union,changed,symbols,contact['moving_indices'],ids,set(contact['caps']),pr['atom_count']);require(cg['pass'],'Final source contact replay failed')
            pair={}
            for z in METALS:
                target,=[t for t in c['targets'] if t['cell']['source_case_id']==sid and t['cell']['metal']==z and t['cell']['candidate']==q];cell=target['cell'];v=row['matrix'][z][q]
                require((cell['xyz'],cell['MACE'],cell['MACE_eV'])==(v['xyz'],v['receipt'],v['energy_eV']),'Current MACE cell differs')
                value=read(cell['MACE']);mace.admitted_result(value,read(value['request']),cell,c['MACE_model'],pf['versions'])
                pair[z]=b.xyz(check(cell['xyz']));require(len(pair[z])==pr['atom_count']==cell['atom_count'] and pair[z][0][0]==z,'Generic endpoint atom/metal differs')
                require(cell['charge']==pr['endpoint_charges'][z] and cell['multiplicity']==1 and sum(f.VALENCE[v[0]] for v in pair[z])-cell['charge']==cell['valence_electrons'],'Actual electron/state inventory differs')
                if q=='amide_flip':require(cell['xyz']['sha256']==next(t for t in pr['task_templates'] if t['metal']==z)['xyz']['sha256'],'Fixed flip bytes changed')
                else:
                    require(np.max(np.abs(np.asarray([v[1:] for v in pair[z]])-coords))<=1e-12,'Final frame/source mapping differs')
                    savedgate=read(v['full_source_admission']);require(savedgate['pass'] and savedgate['task_id']==task['task_id'] and savedgate['active_q_radian']==point['active_q_radian'] and savedgate['contact_inputs']==contact,'Final contact evidence differs')
                mr[cell['task_id']]={'status':'complete','receipt':cell['MACE'],'energy_eV':cell['MACE_eV']}
            require(pair['Ca'][1:]==pair['La'][1:] and pair['Ca'][0][1:]==pair['La'][0][1:],'Paired physical frame differs')
            physical.append({'source_case_id':sid,'candidate':q,'physical_admitted':True})
    require(len(c['targets'])==300 and len(c['cached'])==156 and len(c['cells'])==c['declared_new_unique_cells']<=144,'Cell denominator differs')
    legacy=LegacyComponentFacade(c,r);negative=legacy.negative_metadata_checks()
    cached=[a.reusable(c,d,n,legacy) for d in c['cached']]
    require(len(legacy.evidence)==12 and all(v['accepted'] for v in legacy.evidence),'Exact12 legacy components did not qualify')
    cc=read(c['canonical_collection']);inv=read(c['calibration_inventory']);residuals=[]
    require(cc['inventory']==c['calibration_inventory'] and c['canonical_native_collection']==cc['canonical_native_collection'] and selector.__file__==c['scorer']['path'],'Canonical baseline/scorer differs')
    require([row['record_id'] for row in c['rows']]==inv['calibration_rule']['designated_record_ids'],'Designated25 changed')
    for source,stored in zip(c['rows'],cc['rows']):
        require(source['record_id']==stored['record_id'],'Original canonical order differs');evidence={z:{} for z in METALS}
        for d,v in zip(c['cached'][:150],cached[:150]):
            if d['source_case_id']==source['source_case_id']:evidence[d['metal']][d['candidate']]=v
        matrix={z:{q:component(v) for q,v in evidence[z].items()} for z in METALS};h.reproduce(h.choose(selector,matrix,ORDER3),stored['pool'],residuals)
    declared={t['task_id']:t for t in c['cells']}
    for target in c['targets']:
        t,route=target['cell'],target['route']
        source=c['cached'][route['cache_index']] if route['kind']=='cached' else declared[route['task_id']]
        require(identity(t)==identity(source),'Exact same-source/state/XYZ cache routing differs')
    return {'status':'passed','physical_checks':physical,'reused_MACE':{t['task_id']:mr[t['task_id']] for t in c['cells']},'cached_evidence':cached,
       'baseline_replay':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},'new_molecular_calls':0,'reused_original_cells':150,'reused_A0AC_additional_cells':6,
       'legacy_schema_qualification':{'status':'passed','pins':c['legacy_component'],'real_component_replays':legacy.evidence,'negative_metadata_checks':negative,'scientific_formulas_and_tolerances_changed':False}}

def component(v):
    return {'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} if v['status']=='available' else None

def summarize(c,pf,sr,a,h,n,r,selector):
    native={v['task_id']:v for v in sr};new={};failed=[]
    for cell in c['cells']:
        v=native.get(cell['task_id'])
        try:
            require(v is not None and v['status']=='available','Required new native endpoint unavailable: '+str(v.get('reason') if v else 'not executed'))
            original={**cell,'native_frozen_ALPB_result':v['result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
            new[cell['task_id']]=a.reusable(c,{'original':original,**{k:cell[k] for k in ('maps','state_signature','source_preparation')}},n,r)
        except Exception as e:new[cell['task_id']]={'status':'unavailable','reason':str(e)};failed.append({'task_id':cell['task_id'],'reason':str(e),'native_record':v})
    rows=[];residuals=[]
    for source in c['rows']:
        sid=source['source_case_id'];evidence={z:{} for z in METALS}
        for target in (t for t in c['targets'] if t['cell']['source_case_id']==sid):
            cell,route=target['cell'],target['route'];z,q=cell['metal'],cell['candidate'];v=copy.deepcopy(pf['cached_evidence'][route['cache_index']] if route['kind']=='cached' else new[route['task_id']])
            if v['status']=='available':
                v['native_actual_MACE_use']={'MACE':v['MACE'],'MACE_eV':v['MACE_eV'],'xyz':v['xyz']};v.update(MACE=cell['MACE'],MACE_eV=cell['MACE_eV'],xyz=cell['xyz'],exact_geometry_reuse=route)
            evidence[z][q]=v
        matrix={z:{q:component(v) for q,v in evidence[z].items()} for z in METALS};pools={label:h.choose(selector,matrix,order) for label,order in ORDERS.items()}
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for before,after in [('common3','common4'),('common4','common6'),('common3','common6')]:
            changes[before+'_to_'+after]={}
            for policy in ('operational','mathematical','origin'):
                b,e=pools[before][policy],pools[after][policy];value=None
                if b is not None and e is not None:
                    dz={z:a.difference(works[z][e[z]['candidate']],works[z][b[z]['candidate']]) for z in METALS};delta=a.difference(dz['Ca'],dz['La']);raw=e['R_model_kcal_mol']-b['R_model_kcal_mol'];h.close(raw,delta['total'],residuals)
                    value={'per_metal_change_kcal_mol':dz,'differential_change_kcal_mol':delta,'R_model_change_kcal_mol':raw,'before':{z:b[z]['candidate'] for z in METALS},'after':{z:e[z]['candidate'] for z in METALS}}
                changes[before+'_to_'+after][policy]=value
        rows.append({**source,'pools':pools,'changes':changes,'all_origin_relative_works_kcal_mol':works,'evidence':evidence})
    old=read(c['canonical_collection']);calibrations={label:selector.calibrate([{**row,'pool':row['pools'][label]} for row in rows]) for label in ORDERS}
    for row in rows:
        row['decisions_by_pool']={label:selector.decisions({'pool':row['pools'][label],'expected_class':row['expected_class']},calibrations[label]) for label in ORDERS}
        row['old_band_transfer']={label:selector.decisions({'pool':row['pools'][label],'expected_class':row['expected_class']},old['calibrations']) for label in ORDERS}
    require(calibrations['common3']==old['calibrations'],'Original canonical calibration changed')
    return {'protocol_id':PROTOCOL,'declared_sources':25,'declared_new_unique_cells':len(c['cells']),'available_new_unique_cells':sum(v['status']=='available' for v in new.values()),
       'reused_original_cells':150,'reused_A0AC_additional_cells':6,'all25_common6_available':all(row['pools']['common6']['status']=='available' for row in rows),
       'rows':rows,'calibrations':calibrations,'old_calibrations':old['calibrations'],'failed_new_cells':failed,'pool_orders':ORDERS,
       'algebra_checks':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
       'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'production_changed':False,'automatic_follow_on':False,'development_references_repeatedly_inspected':True}

def execute(path,sha,action):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports');cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['implementation']==pin(__file__),'Wrong frozen adapter');matrix=load('_canonical_solvent_placement',c['matrix']);placement=matrix.placement()
    out=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID']);start=time.monotonic();pf=None;sr=[];receipt={'status':'failed','config':cp,'placement':placement,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'));dispatcher=load('_canonical_solvent_dispatch',c['dispatcher']);f=load('_canonical_factorial',c['factorial']);r=load('_canonical_component',c['component']);n=load('_canonical_native',c['native_adapter'])
            a=load('_canonical_reuse',c['solvent_adapter']);h=load('_canonical_pool',c['pool_adapter']);selector=load('_canonical_selector',c['scorer']);search=load('_canonical_search',c['score_adapter']);mace=load('_canonical_mace_predicate',c['mace_scorer'])
            pf=preflight(c,f,r,n,a,h,selector,search,mace);write(out/'PREFLIGHT.json',pf);write(out/'LEGACY_COMPONENT_SCHEMA_TEST.json',pf['legacy_schema_qualification'])
            if action=='execute':sr=dispatcher.solvent(c,out,f,r,n,pf['reused_MACE'])
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
        finally:
            try:
                if pf and action=='execute':
                    if not sr:
                        for cell in c['cells']:
                            p=out/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                            if p.exists():
                                v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    write(out/'PARTIAL_RESULTS.json',{'config':cp,'native_rows':sr,'declared_new_unique_cells':len(c['cells'])})
                    result={'config':cp,**summarize(c,pf,sr,a,h,n,r,selector)};write(out/'FINAL_COLLECTION.json',result)
                    if not result['all25_common6_available']:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
                try:receipt.update(SCC_started=sum(1 for _ in out.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in out.glob('solvent/points/*/*/EXECUTION.json')),native_component_calls=sum(1 for _ in out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
                    measured_resource_receipts={'SCC':[pin(p) for p in sorted(out.glob('solvent/points/*/*/EXECUTION.json'))],'native_ALPB':[pin(p) for p in sorted(out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]})
                except BaseException as e:receipt.update(status='failed',cost_reporting_error=str(e))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);q=s.add_parser('bind');q.add_argument('--score-config',required=True);q.add_argument('--score-collection',required=True);q.add_argument('--output',required=True)
    for action in ('preflight','execute'):
        q=s.add_parser(action);q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='bind':bind(a.score_config,a.score_collection,a.output)
    else:sys.exit(execute(a.config,a.config_sha,a.action))
