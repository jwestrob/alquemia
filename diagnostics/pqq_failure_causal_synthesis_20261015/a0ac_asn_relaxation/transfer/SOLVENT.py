#!/usr/bin/env python3
"""Transfer population join around the executed fixed-qvac solvent kernels."""
import argparse, collections, fcntl, os, shutil, time, traceback
from COMMON import *
PROTOCOL='uniform_terminal_amide_transfer_solvent_v1'

def score_gate(sc,sr,ex,cp):
    require(sr['config']==ex['config']==cp and ex['status'] in ('complete','partial') and sr['stage_complete'] and sr['all_contexts_terminal'],'Actual terminal full MACE stage required')
    require(sr['pool_orders']==ORDERS and sc['pool_orders']==ORDERS,'Common-pool definition differs')
    require(len(sr['rows'])==257 and {r['context_id'] for r in sr['rows']}==set(sc['all_context_ids']),'Full context denominator differs')
    require(sc['stage'] in ('all','remaining'),'Scout cannot release transfer solvent')

def bind(score_config,score_collection,output):
    cp=pin(score_config);sc=read(cp);rp=pin(score_collection);sr=read(rp);ep=pin(Path(score_collection).with_name('EXECUTION.json'));score_gate(sc,sr,read(ep),cp)
    prep=read(sc['preparation']);gm=read(sc['preparation_manifest']);gate=canonical_gate(gm['canonical_gate']['collection']);require(gate==gm['canonical_gate'],'Canonical gate changed')
    template=read(gate['config']);kernel=load('_transfer_solvent_metadata',template['implementation']);inv=read(gm['inventory']);ledger=read(gm['ledger']);a=ledger['artifacts']
    native={}
    for authority in ('canonical_native_collection','transfer_native_collection'):native.update(read(a[ledger['authority'][authority]])['tasks'])
    records={r['record_id']:r for r in inv['records']};sources={s['context_id']:s for s in gm['sources']};prepared={r['context_id']:r for r in prep['rows']}
    cached=[];origins={};targets=[];cells=[];seen={};new_keys={};rows=[]
    def cache(d):
        key=kernel.identity(d)
        if key not in seen:seen[key]=len(cached);cached.append(d)
        return seen[key]
    for row in sr['rows']:
        cid=row['context_id'];s=sources[cid];source=s['record'];sid=source['source_case_id'];pr=prepared[cid];origins[cid]={}
        require(row['source_case_id']==sid and row['case_id']==source['case_id'] and row['memberships']==s['memberships'],'MACE source/context membership changed')
        for z in METALS:
            for q in ORDER3:
                v=source['matrix'][z][q];nr=native[v['native_task_id']];d=kernel.descriptor(sid,z,q,v,nr['result']);cache(d)
                if q=='origin':origins[cid][z]=d['original']
        if 'reuse' in s:
            for z in METALS:
                for q in ORDER6[3:]:
                    e=s['reuse']['solvent_row']['evidence'][z][q];actual=e.get('native_actual_MACE_use',e)
                    cache(kernel.descriptor(sid,z,q,{k:actual[k] for k in ('xyz','MACE','MACE_eV')},e['native_result']))
        for z in METALS:
            for q in ORDER6:
                v=row['matrix'][z][q];paired=all(row['matrix'][zz][q]['status']=='complete' for zz in METALS)
                if not paired:
                    targets.append({'context_id':cid,'metal':z,'candidate':q,'route':{'kind':'unavailable','reason':'Required paired MACE geometry unavailable'},'MACE_cell':v});continue
                base=origins[cid][z]
                t={'task_id':cid+'__'+q+'__'+z,'context_id':cid,'source_case_id':sid,'case_id':source['case_id'],'candidate':q,'metal':z,'xyz':v['xyz'],'MACE':v['receipt'],'MACE_eV':v['energy_eV'],
                  **{k:base[k] for k in ('charge','multiplicity','atom_count','valence_electrons')},'maps':source['maps'] if q in ORDER3 else {zz:pr['mapping'] for zz in METALS},
                  'state_signature':source['state_signature'],'source_preparation':source['source_preparation'],'prepared_admission':pr.get('admission'),'proposal_receipt':v.get('proposal_receipt'),'full_source_admission':v.get('full_source_admission')}
                key=kernel.identity(t)
                if key in seen:route={'kind':'cached','cache_index':seen[key]}
                elif key in new_keys:route={'kind':'new','task_id':new_keys[key]}
                else:
                    require(q not in ORDER3 and 'reuse' not in s,'Exact mandatory cached endpoint missing')
                    cells.append(t);new_keys[key]=t['task_id'];route={'kind':'new','task_id':t['task_id']}
                targets.append({'context_id':cid,'metal':z,'candidate':q,'cell':t,'route':route})
        rows.append({**{k:source[k] for k in ('source_case_id','case_id','record_id','selection_id','protein','expected_class','biological_group','source_conditioning_metal','family')},'context_id':cid,'memberships':s['memberships']})
    require(len(rows)==257 and len(targets)==257*12 and len(cells)<=223*6,'Finite transfer denominator differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('COMMON.py','SOLVENT.py','SOLVENT.sbatch','PLAN.md'):shutil.copyfile(HERE/name,out/name)
    c={k:template[k] for k in kernel.KEYS};c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),common=pin(out/'COMMON.py'),wrapper=pin(out/'SOLVENT.sbatch'),plan=pin(out/'PLAN.md'),
      template=gate['config'],canonical_gate=gate,canonical_kernel=template['implementation'],dispatcher=template['dispatcher'],legacy_component=template['legacy_component'],
      score_config=cp,score_collection=rp,score_execution=ep,score_adapter=sc['implementation'],mace_scorer=sc['old_scorer'],preparation=sc['preparation'],preparation_manifest=sc['preparation_manifest'],
      inventory=gm['inventory'],ledger=gm['ledger'],historical_collection=a[ledger['authority']['frozen_qvac_final_collection']],generic_contact_helper=gm['generic_contact_helper'],
      rows=rows,cached=cached,targets=targets,cells=cells,CPUs=7,protein='uniform_terminal_amide_transfer',declared_contexts=257,declared_memberships=260,
      pool_orders=ORDERS,maximum_SCC_calls=4*len(cells),maximum_native_component_calls=2*len(cells),new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0,automatic_submission=False)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def physical(c,search,b,scorer,completion,np,pf,f):
    """Same final q, physical and full-source checks as canonical solvent preflight."""
    from mace_site_kinematics import Kinematics
    sr=read(c['score_collection']);prep={r['context_id']:r for r in read(c['preparation'])['rows']};worker=read(read(c['score_config'])['worker_manifest']);helper=load('_transfer_solvent_contacts',c['generic_contact_helper']);mr={};checks=[]
    for row in sr['rows']:
        cid=row['context_id'];pr=prep[cid]
        if pr['status']!='prepared':continue
        contact=pr['contact'];kin=Kinematics(read(pr['mapping'])['context'])
        with np.load(check(contact['origin']),allow_pickle=False) as data:union=data['positions'].copy();symbols=data['symbols'].tolist()
        ids=read(contact['identities'])['centers']
        for q in ORDER6[3:]:
            if not all(row['matrix'][z][q]['status']=='complete' for z in METALS):continue
            coords=None
            if q!='amide_flip':
                owner=q.removeprefix('relaxed_flip_');proposal=read(row['proposals'][owner]);point=proposal['proposal'];pm=read(proposal['manifest']);task,=[t for t in pm['tasks'] if t['task_id']==proposal['task_id']]
                require(task['mapping']==pr['mapping'] and task['contact']==contact and pm['settings']==worker['settings'],'Final task/map/settings differ')
                original=next(t for t in pr['task_templates'] if t['metal']==owner)
                for k in ('active_indices','active_mode_ids','selector','mode_count'):require(task[k]==original[k],'Original selected mode changed')
                full=completion.angular.full_q(task,np.asarray(point['active_q_radian']));require(np.array_equal(full,np.asarray(point['full_q'])),'Final full q differs')
                coords=kin.evaluate(full)[1];fg=completion.angular.final_geometry(kin,task,np.asarray(point['active_q_radian']),[v[0] for v in b.xyz(check(point['coordinate']))]);require(fg['physical_feasible'],'Final physical replay failed')
                changed=union.copy();changed[:pr['atom_count']]=coords;require(helper.contact_check(np,union,changed,symbols,contact['moving_indices'],ids,set(contact['caps']),pr['atom_count'])['pass'],'Final full-source contact failed')
            pair={}
            for z in METALS:
                target,=[v for v in c['targets'] if (v['context_id'],v['metal'],v['candidate'])==(cid,z,q)];cell=target['cell'];v=row['matrix'][z][q]
                require((cell['xyz'],cell['MACE'],cell['MACE_eV'])==(v['xyz'],v['receipt'],v['energy_eV']),'MACE endpoint changed')
                value=read(cell['MACE']);scorer.admitted_result(value,read(value['request']),cell,c['MACE_model'],pf['versions'])
                pair[z]=b.xyz(check(cell['xyz']));require(len(pair[z])==pr['atom_count']==cell['atom_count'] and pair[z][0][0]==z,'Endpoint count/metal differs')
                require(cell['charge']==pr['endpoint_charges'][z] and cell['multiplicity']==1 and sum(f.VALENCE[v[0]] for v in pair[z])-cell['charge']==cell['valence_electrons'],'Electron/state inventory differs')
                if q=='amide_flip':require(cell['xyz']['sha256']==next(t for t in pr['task_templates'] if t['metal']==z)['xyz']['sha256'],'Fixed flip differs')
                else:
                    require(np.max(np.abs(np.asarray([v[1:] for v in pair[z]])-coords))<=1e-12,'Final physical frame differs')
                    cg=read(v['full_source_admission']);require(cg['pass'] and cg['task_id']==task['task_id'] and cg['active_q_radian']==point['active_q_radian'] and cg['contact_inputs']==contact,'Saved final contact evidence differs')
                mr[cell['task_id']]={'status':'complete','receipt':cell['MACE'],'energy_eV':cell['MACE_eV']}
            require(pair['Ca'][1:]==pair['La'][1:] and pair['Ca'][0][1:]==pair['La'][0][1:],'Paired coordinates differ');checks.append({'context_id':cid,'candidate':q,'pass':True})
    return mr,checks

def preflight(c,kernel,f,r,n,a,h,selector):
    require(pin(sys.executable)==c['numerical_python'] and c['protocol_id']==PROTOCOL,'Numerical interpreter/protocol differs');template=read(c['template'])
    require(canonical_gate(c['canonical_gate']['collection'])==c['canonical_gate'] and c['canonical_gate']['config']==c['template'],'Canonical gate changed')
    for name in kernel.KEYS:require(c[name]==template[name],'Qualified solvent field differs: '+name)
    for name in ('common','plan','wrapper','score_adapter','generic_contact_helper','dispatcher'):check(c[name])
    require(c['dispatcher']==template['dispatcher'] and c['legacy_component']==template['legacy_component'],'Qualified legacy/dispatch route differs')
    require(read(c['ORCA_correction_policy'])['candidate_sha256']==c['corrected_ORCA_parameter_provenance']['sha256'],'Full La parameter lineage differs');check(c['corrected_ORCA_parameter_provenance'])
    sc=read(c['score_config']);sr=read(c['score_collection']);score_gate(sc,sr,read(c['score_execution']),c['score_config'])
    search=load('_transfer_MACE_replay',c['score_adapter']);worker=read(sc['worker_manifest']);b=load('_transfer_MACE_worker',worker['implementation']);scorer=load('_transfer_MACE_predicate',c['mace_scorer']);b.parent_env();b.modules(worker)
    import numpy as np
    pf,completion=search.preflight(sc,b,scorer,np)
    require(search.collect(sc,Path(c['score_collection']['path']).parent,scorer,pf)==sr,'Actual transfer MACE collection replay differs')
    search.install_admission(sc,b,scorer,np)
    mr,checks=physical(c,search,b,scorer,completion,np,pf,f)
    runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for name in ('xtb','xtb_parameter_file','tblite_python'):check(runtime[name])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:check(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:check(p)
    for p in list(nm['implementation'].values())+nm['build']['pins']+build['linked_libraries']:check(p)
    check(build['binary']);require(build['manifest']==c['native_manifest'] and build['binary']==c['native_binary'] and nm['MACE_model']==c['MACE_model']==worker['model'],'Native build/model differs')
    legacy=kernel.LegacyComponentFacade(c,r);negative=legacy.negative_metadata_checks();cached=[a.reusable(c,d,n,legacy) for d in c['cached']]
    require(all(v['accepted'] for v in legacy.evidence),'Legacy component replay failed')
    declared={v['task_id']:v for v in c['cells']};require(len(declared)==len(c['cells']),'Repeated molecular dispatch ID')
    for target in c['targets']:
        route=target['route']
        if route['kind']=='unavailable':continue
        actual=c['cached'][route['cache_index']] if route['kind']=='cached' else declared[route['task_id']]
        require(kernel.identity(target['cell'])==kernel.identity(actual),'Native source/state/coordinate route differs')
    historical={v['record_id']:v for v in read(c['historical_collection'])['rows']};residuals=[]
    for row in c['rows']:
        matrix={z:{} for z in METALS}
        for target in c['targets']:
            if target['context_id']!=row['context_id'] or target['candidate'] not in ORDER3:continue
            require(target['route']['kind']=='cached','Original common3 was dispatched');value=cached[target['route']['cache_index']];matrix[target['metal']][target['candidate']]=kernel.component(value)
        actual=h.choose(selector,matrix,ORDER3)
        for rid in row['memberships']:h.reproduce(actual,historical[rid]['pool'],residuals)
    return {'status':'passed','cached_evidence':cached,'physical_checks':checks,'reused_MACE':{cell['task_id']:mr[cell['task_id']] for cell in c['cells']},
      'baseline_replay':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},'legacy_schema_qualification':{'replays':legacy.evidence,'negative_metadata_checks':negative},'new_molecular_calls':0}

def summarize(c,pf,sr,kernel,a,h,n,r,selector):
    native={v['task_id']:v for v in sr};new={};failed=[];residuals=[]
    for cell in c['cells']:
        v=native.get(cell['task_id'])
        try:
            require(v is not None and v['status']=='available','Required new native endpoint unavailable')
            original={**cell,'native_frozen_ALPB_result':v['result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
            new[cell['task_id']]=a.reusable(c,{'original':original,**{k:cell[k] for k in ('maps','state_signature','source_preparation')}},n,r)
        except Exception as e:new[cell['task_id']]={'status':'unavailable','reason':str(e)};failed.append({'task_id':cell['task_id'],'reason':str(e),'native_record':v})
    unique=[];records={v['record_id']:v for v in read(c['inventory'])['records']};members=[];cal=c['canonical_gate']['calibrations'];old=c['canonical_gate']['old_calibrations']
    for source in c['rows']:
        evidence={z:{} for z in METALS}
        for target in (v for v in c['targets'] if v['context_id']==source['context_id']):
            route=target['route'];z,q=target['metal'],target['candidate']
            if route['kind']=='unavailable':value={'status':'unavailable','reason':route['reason'],'MACE_cell':target['MACE_cell']}
            else:
                value=copy.deepcopy(pf['cached_evidence'][route['cache_index']] if route['kind']=='cached' else new[route['task_id']]);cell=target['cell']
                if value['status']=='available':value['native_actual_MACE_use']={k:value[k] for k in ('MACE','MACE_eV','xyz')};value.update(MACE=cell['MACE'],MACE_eV=cell['MACE_eV'],xyz=cell['xyz'],exact_geometry_reuse=route)
            evidence[z][q]=value
        matrix={z:{q:kernel.component(v) for q,v in evidence[z].items()} for z in METALS};pools={label:h.choose(selector,matrix,order) for label,order in ORDERS.items()}
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for before,after in [('common3','common4'),('common4','common6'),('common3','common6')]:
            changes[before+'_to_'+after]={}
            for policy in ('operational','mathematical','origin'):
                b,e=pools[before][policy],pools[after][policy];delta=None
                if b is not None and e is not None:
                    dz={z:a.difference(works[z][e[z]['candidate']],works[z][b[z]['candidate']]) for z in METALS};dd=a.difference(dz['Ca'],dz['La']);raw=e['R_model_kcal_mol']-b['R_model_kcal_mol'];h.close(raw,dd['total'],residuals)
                    delta={'per_metal_change_kcal_mol':dz,'differential_change_kcal_mol':dd,'R_model_change_kcal_mol':raw,'before':{z:b[z]['candidate'] for z in METALS},'after':{z:e[z]['candidate'] for z in METALS}}
                changes[before+'_to_'+after][policy]=delta
        result={**source,'pools':pools,'changes':changes,'all_origin_relative_works_kcal_mol':works,'evidence':evidence};unique.append(result)
        for rid in source['memberships']:
            record=records[rid];row={**result,**{k:record[k] for k in ('record_id','source_case_id','case_id','selection_id','protein','biological_group','expected_class','source_conditioning_metal','family','role')}}
            row['decisions_by_pool']={label:selector.decisions({'pool':pools[label],'expected_class':row['expected_class']},cal[label]) for label in ORDERS}
            row['old_band_transfer']={label:selector.decisions({'pool':pools[label],'expected_class':row['expected_class']},old) for label in ORDERS};members.append(row)
    ledger=read(c['ledger']);ids_by_cohort=cohort_ids(ledger);expected=set(ids_by_cohort['expanded135'])|set(ids_by_cohort['current_ca125']);require(len(members)==len({v['record_id'] for v in members})==260 and {v['record_id'] for v in members}==expected,'Final cohort membership changed')
    cohorts={name:{label:{policy:dict(collections.Counter(v['decisions_by_pool'][label][policy]['outcome'] for v in members if v['record_id'] in ids)) for policy in selector.POLICIES} for label in ORDERS} for name,ids in ids_by_cohort.items()}
    triples={label:selector.triples(read(c['inventory']),[{**v,'pool':v['pools'][label]} for v in members],cal[label]) for label in ORDERS}
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_contexts':257,'declared_memberships':260,'unique_context_rows':unique,'rows':members,'cohorts':cohorts,'existing_triples':triples,
      'historical_ledger':c['ledger'],'superseded_primary_records':[v for v in ledger['records'] if v['record_id'] not in expected],
      'original_historical_collection':c['historical_collection'],'calibration_reference':c['canonical_gate'],'new_calibration':False,'failed_new_cells':failed,
      'declared_new_unique_cells':len(c['cells']),'available_new_unique_cells':sum(v['status']=='available' for v in new.values()),'pool_orders':ORDERS,
      'algebra_checks':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},'production_changed':False,'development_evidence':True}

def execute(path,sha,action):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports');cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['implementation']==pin(__file__) and c['common']==pin(Path(__file__).with_name('COMMON.py')),'Frozen transfer implementation differs')
    matrix=load('_transfer_solvent_placement',c['matrix']);placement=matrix.placement();out=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID']);start=time.monotonic();pf=None;sr=[]
    receipt={'status':'failed','config':cp,'placement':placement,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'));kernel=load('_transfer_canonical_solvent_kernel',c['canonical_kernel']);dispatcher=load('_transfer_solvent_dispatcher',c['dispatcher'])
            f=load('_transfer_factorial',c['factorial']);r=load('_transfer_component',c['component']);n=load('_transfer_native',c['native_adapter']);a=load('_transfer_reuse',c['solvent_adapter']);h=load('_transfer_pool',c['pool_adapter']);selector=load('_transfer_policy',c['scorer'])
            pf=preflight(c,kernel,f,r,n,a,h,selector);write(out/'PREFLIGHT.json',pf)
            if action=='execute':sr=dispatcher.solvent(c,out,f,r,n,pf['reused_MACE']) if c['cells'] else []
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
        finally:
            try:
                if pf and action=='execute':
                    if not sr:
                        for cell in c['cells']:
                            p=out/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                            if p.exists():v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    write(out/'PARTIAL_RESULTS.json',{'config':cp,'native_rows':sr});result=summarize(c,pf,sr,kernel,a,h,n,r,selector);write(out/'FINAL_COLLECTION.json',{'config':cp,**result})
                    if result['available_new_unique_cells']!=len(c['cells']):receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_error=str(e),collection_traceback=traceback.format_exc())
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
                try:receipt.update(SCC_started=sum(1 for _ in out.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in out.glob('solvent/points/*/*/EXECUTION.json')),native_component_calls=sum(1 for _ in out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),measured_resource_receipts={'SCC':[pin(p) for p in sorted(out.glob('solvent/points/*/*/EXECUTION.json'))],'native_ALPB':[pin(p) for p in sorted(out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]})
                except BaseException as e:receipt.update(status='failed',cost_reporting_error=str(e))
                write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);b=s.add_parser('bind');b.add_argument('--score-config',required=True);b.add_argument('--score-collection',required=True);b.add_argument('--output',required=True)
    for action in ('preflight','execute'):
        x=s.add_parser(action);x.add_argument('--config',required=True);x.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='bind':bind(a.score_config,a.score_collection,a.output)
    else:sys.exit(execute(a.config,a.config_sha,a.action))
