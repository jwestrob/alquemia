#!/usr/bin/env python3
"""Canonical25 geometry-only rebase. Metadata freeze; numerical work in Slurm."""
import argparse, copy, hashlib, importlib.util, json, os, shutil, socket, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_relaxation'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
GEOMETRY_SHA='5f69e8c7d57a71d6601d7ced65e8dcc419c8ad2ac154f755141b6faaf395c968'
REUSE_ID='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3'
METALS=('Ca','La'); MOVED=('OD1','ND2','HD21','HD22'); PROTOCOL='canonical25_fixed_flip_rebased_geometry_only_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    require(pin(p['path'])==p,'Changed input: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def require(v,why):
    if not v:raise ValueError(why)
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def key(s):return tuple(s.get(k,'') for k in ('chain','resnum','insertion_code','resname','atom','element'))

def freeze(output):
    gp=pin(WORK/'asn_canonical_geometry_v2/execution/PREPARATION.json');require(gp['sha256']==GEOMETRY_SHA,'Wrong actual canonical geometry')
    geometry=read(gp);gm=read(geometry['manifest']);scope=read(gm['scope_manifest'])
    require(geometry['status']=='complete' and geometry['prepared_cells']==50 and len(geometry['rows'])==25,'Actual complete canonical geometry required')
    require(read(WORK/'asn_canonical_geometry_v2/execution/A0AC_REGRESSION.json')['status']=='passed','Actual fixed-flip regression missing')
    priorp=pin(WORK/'asn_relaxation_v2/preparation/PREPARATION.json');prior=read(priorp);pc=read(prior['config'])
    oldscore=read(pc['old_config']);worker=read(oldscore['worker_manifest'])
    finalp=pin(WORK/'asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json');final=read(finalp)
    require(prior['status']=='complete' and final['all20_searches_and_common6_complete'],'Actual prior rebase/searches incomplete')
    sources=[]
    for r,sc in zip(geometry['rows'],scope['rows']):
        require(r['record_id']==sc['record_id'] and r['source_case_id']==sc['source_case_id'],'Canonical row/order differs')
        sm=read(sc['search_manifest']);require(sm['settings']==worker['settings'],'Source optimizer settings differ')
        require([t for t in sm['tasks'] if t['case_id']==sc['case_id']]==sc['original_tasks'],'Original paired task records differ')
        for p in (r['geometry_admission'],r['state_signature'],r['source_preparation'],r['protonated_source'],*r['maps'].values()):check(p)
        for z in METALS:check(r['matrix'][z]['amide_flip']['xyz']);check(r['matrix'][z]['amide_flip']['original_xyz'])
        sources.append({'row':r,'scope':sc,'contact_identities':pin(Path(r['geometry_admission']['path']).with_name('CONTACT_IDENTITIES.json'))})
    source,=[s for s in sources if s['row']['source_case_id']==REUSE_ID]
    oldrow,=[r for r in prior['rows'] if r['source_case_id']==REUSE_ID];executed,=[r for r in final['rows'] if r['source_case_id']==REUSE_ID]
    require(oldrow['status']=='prepared' and executed['status']=='complete','Real A0AC alias unavailable')
    require(all(source['row'][k]==oldrow['original_row'][k] for k in ('maps','state_signature','source_preparation','protonated_source')),'A0AC source/state aliases differ')
    for z in METALS:
        t,=[t for t in source['scope']['original_tasks'] if t['metal']==z]
        require(t==oldrow['old_searches'][z]['task'],'A0AC original task differs')
        require(source['row']['matrix'][z]['amide_flip']['xyz']['sha256']==oldrow['original_row']['matrix'][z]['amide_flip']['xyz']['sha256'],'A0AC flip bytes differ')
        for candidate in ('origin','adaptive_Ca','adaptive_La'):
            a=source['scope']['six_saved_cells'][z][candidate];b=executed['matrix'][z][candidate]
            require(a['MACE']==b['receipt'] and a['xyz']==b['xyz'],'A0AC original six-cell aliases differ');check(a['MACE']);check(a['xyz'])
        check(executed['proposals'][z])
        for cell in executed['matrix'][z].values():check(cell['receipt']);check(cell['xyz'])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir()
    implementation={}
    for name in ('PREPARE.py','RUN.sbatch','PLAN.md'):
        target=out/'implementation'/name;shutil.copyfile(HERE/name,target);implementation[name]=pin(target)
    c={'protocol':PROTOCOL,'implementation':implementation,'geometry':gp,'geometry_manifest':geometry['manifest'],
       'scope_manifest':gm['scope_manifest'],'canonical_collection':geometry['canonical_collection'],
       'prior_preparation':priorp,'prior_config':prior['config'],'prior_implementation':pc['implementation']['PREPARE.py'],
       'completed_search_collection':finalp,'completed_search_execution':pin(WORK/'asn_relaxation_remaining_v1/execution/EXECUTION.json'),
       'source_preparation_helper':pc['source_preparation_helper'],'source_parser':pc['source_parser'],
       'generic_contact_helper':gm['implementation']['PREPARE.py'],'python':pc['python'],
       'source_snapshot':worker['source_snapshot'],'source_modules':worker['source_modules'],'worker_manifest':oldscore['worker_manifest'],
       'settings':worker['settings'],'technical_recovery':pc['technical_recovery'],
       'sources':sources,'reuse':{'source_case_id':REUSE_ID,'prepared_row':oldrow,'executed_row':executed},
       'declared_sources':25,'new_energy_calls':0,'new_optimizations':0,'output':str(out/'execution')}
    write(out/'manifest.json',c);print(json.dumps(pin(out/'manifest.json')))

def base_row(s):
    r=s['row']
    return {**{k:r[k] for k in ('source_case_id','case_id','record_id','selection_id','protein','atom_count','role','endpoint_charges','water_inventory','terminal_mode_selected_by_metal','original_four_modes')},
       'status':'unavailable','reason':'not prepared','original_row':r,'geometry_only':True,'new_energy_calls':0,'new_optimizations':0}

def prepare_row(s,d,c,np,Kinematics,angular,source_helper,base,generic):
    r=s['row'];n=r['atom_count'];gate=read(r['geometry_admission']);cp=read(r['source_preparation']);state=read(r['state_signature'])
    require(r['status']=='prepared' and gate['pass'],'Fixed flip not admitted')
    require(cp['new_atom_count']==n and state['endpoint_charges']==r['endpoint_charges'] and state['water_inventory']==r['water_inventory'] and state['multiplicities']=={'Ca':1,'La':1},'Actual source state/count differs')
    maps={z:read(r['maps'][z]) for z in METALS};require(maps['Ca']==maps['La'],'Original paired maps differ');original=maps['Ca'];g=copy.deepcopy(original)
    before_atoms=base.xyz_atoms(check(r['matrix']['Ca']['amide_flip']['original_xyz']));flips={z:base.xyz_atoms(check(r['matrix'][z]['amide_flip']['xyz'])) for z in METALS}
    require(len(before_atoms)==n and all(len(flips[z])==n and flips[z][0]['element']==z for z in METALS),'Wrong paired atom counts/elements')
    require([a['xyz'] for a in flips['Ca']]==[a['xyz'] for a in flips['La']] and [a['element'] for a in flips['Ca'][1:]]==[a['element'] for a in flips['La'][1:]],'Paired fixed-flip frames differ')
    before=np.asarray([a['xyz'] for a in before_atoms]);after=np.asarray([a['xyz'] for a in flips['Ca']]);symbols=[a['element'] for a in flips['Ca']]
    entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']};require(set(entries)==set(range(1,n)),'Incomplete current context map')
    source_qm={key(x['source']):i for i,x in entries.items() if 'source' in x}
    named=r['named_Asn_indices_zero_based'];moved={named[a] for a in MOVED};require(moved==set(r['moved_indices_zero_based'])==set(gate['moved_indices_zero_based']),'Moved terminal identities differ')
    physical_qm={i:source_qm[key(x)] for i,x in enumerate(g['context']['source_atom_metadata']) if x and key(x) in source_qm}
    oldkin=Kinematics(original['context']);zero=np.zeros(len(oldkin.modes));old_core=oldkin.evaluate(zero)[1]
    require(np.max(np.abs(old_core-before))<=1e-12,'Original full q0 does not replay')
    changed=[]
    for i,qm in physical_qm.items():
        if qm in moved:g['context']['positions_A'][i]=after[qm].tolist();changed.append(i)
    require(len(changed)==4,'Four source sites not mapped');g['context']['core_positions_A']=after.tolist()
    kin=Kinematics(g['context']);p,core,_,_=kin.evaluate(zero)
    require(np.max(np.abs(core-after))<=1e-12,'Rebased full q0 does not replay')
    require(all(g['context'][k]==original['context'][k] for k in original['context'] if k not in ('positions_A','core_positions_A')),'Physical topology or cap Jacobian definition changed')
    fixed=[i for i in range(len(p)) if i not in changed];require(np.array_equal(kin.positions[fixed],oldkin.positions[fixed]),'Physical exterior changed')
    mapping_checks=kin.check(zero);require(mapping_checks['pass'],'Existing analytic Jacobian/covalent check failed');g['checks']=mapping_checks
    g['rebased_origin']={'original_mapping':r['maps'],'admitted_flip':r['matrix'],'geometry_admission':r['geometry_admission'],'changed_physical_indices':changed,'new_domain_origin':True}
    caps=set(gate['cap_indices']);require(caps=={i for i,x in entries.items() if x['kind']=='sigma_link_H'} and len(caps)==len(state['caps']),'Original cap inventory differs')
    cap_order=sorted(caps)
    require(np.array_equal(after[cap_order],before[cap_order]),'Raw before/after caps differ')
    require(np.array_equal(kin.core[cap_order],oldkin.core[cap_order]),'Raw map caps differ')
    require(np.array_equal(core[cap_order],old_core[cap_order]),'Computed before/after caps differ')
    cap_comparison={'raw_before_after_exact':True,'raw_map_caps_exact':True,'old_new_computed_caps_exact':True,
      'old_cap_q0_residual_A':float(np.max(np.abs(old_core[cap_order]-before[cap_order]))),
      'new_cap_q0_residual_A':float(np.max(np.abs(core[cap_order]-after[cap_order]))),
      'unchanged_full_q0_tolerance_A':1e-12,'diagnostic':c['technical_recovery']['diagnostic_result']}
    ids=read(s['contact_identities'])['centers'];require(len(ids)==gate['contact_union_atom_count'],'Fixed-flip source-union count differs')
    physical=source_helper.physical_source(check(r['protonated_source']),base);lookup={key(a['source']):a for a in physical}
    require(len(physical)==gate['full_source_atom_count'],'Actual full-source atom count differs')
    exterior=[lookup[key(a['source'])]['xyz'] for a in ids[n:]]
    union=np.vstack((after,np.asarray(exterior)));union_symbols=symbols+[lookup[key(a['source'])]['element'] for a in ids[n:]]
    require(len(union)==gate['contact_union_atom_count'] and len(exterior)==gate['appended_source_atoms'],'Actual source-union boundary differs')
    oldtasks={t['metal']:t for t in s['scope']['original_tasks']}
    require(set(oldtasks)==set(METALS) and oldtasks['Ca']['active_mode_ids']==oldtasks['La']['active_mode_ids'] and oldtasks['Ca']['active_indices']==oldtasks['La']['active_indices'],'Paired original modes differ')
    active=oldtasks['Ca']['active_indices'];require(len(active)==4 and [kin.modes[i]['id'] for i in active]==oldtasks['Ca']['active_mode_ids'],'Original four mode IDs/order differ')
    moving={i for k in active for i in kin.modes[k]['moving_indices']};require(moving<=set(physical_qm),'Selected moving source member missing from physical context')
    require(all(anchor not in kin.modes[k]['moving_indices'] or (kin.modes[k]['unit']=='radian' and anchor in kin.modes[k]['axis_indices'])
                for link in kin.links if link[1]=='cap' for anchor in link[2:4] for k in active),'Unsupported selected moving cap boundary')
    moved_qm=sorted({physical_qm[i] for i in moving});require(not caps.intersection(moved_qm),'Synthetic cap motion added')
    contacts=generic.contact_check(np,union,union,union_symbols,moved_qm,ids,caps,n);require(contacts['pass'],'Source-union q0 contact check failed')
    final_checks={}
    for z in METALS:
        t=oldtasks[z];require(t['mode_count']==len(kin.modes) and t['charge']==state['endpoint_charges'][z] and t['multiplicity']==state['multiplicities'][z] and t['mapping']==r['maps'][z] and t['xyz']==r['matrix'][z]['amide_flip']['original_xyz'],'Original task/state join differs')
        final_checks[z]=angular.final_geometry(kin,t,np.zeros(4),[a['element'] for a in flips[z]])
    write(d/'mapping_candidate.json',g)
    np.savez_compressed(d/'contact_origin_candidate.npz',positions=union,symbols=np.asarray(union_symbols));write(d/'contact_identities_candidate.json',{'centers':ids})
    mapping=pin(d/'mapping_candidate.json');contact={'origin':pin(d/'contact_origin_candidate.npz'),'identities':pin(d/'contact_identities_candidate.json'),'moving_indices':moved_qm,'caps':cap_order,'source':r['protonated_source']}
    admission={'status':'prepared','pass':True,'mapping_checks':mapping_checks,'q0_paired_replay':True,'caps_fixed':True,'cap_comparison':cap_comparison,
      'all_physical_H_retained':True,'source_contacts':contacts,'old_mode_ids':oldtasks['Ca']['active_mode_ids'],
      'source_case_id':r['source_case_id'],'original_geometry_admission':r['geometry_admission'],'local_domain_reference':'admitted_flip',
      'original_discrete_displacement_not_in_local_domain':True,'q0_geometry_checks':final_checks,'new_energy_calls':0}
    write(d/'ADMISSION.json',admission)
    templates=[]
    for z in METALS:
        t={k:oldtasks[z][k] for k in ('active_indices','active_mode_ids','active_roles','selection_id','selector','metal','charge','multiplicity','mode_count')}
        t.update(source_case_id=r['source_case_id'],xyz=r['matrix'][z]['amide_flip']['xyz'],mapping=mapping,source_preparation=r['source_preparation'],
          atom_count=n,contact=contact,q0_status='geometry_only_energy_and_force_gradient_unavailable',q0=None,origin_reuse=None,settings=c['settings'])
        templates.append(t)
    row={**base_row(s),'status':'prepared','reason':None,'mapping':mapping,'contact':contact,'admission':pin(d/'ADMISSION.json'),
      'task_templates':templates,'force_gradient_qualification':'UNRUN','optimization_executable':False}
    if r['source_case_id']==REUSE_ID:
        prior=c['reuse']['prepared_row'];pm=read(prior['mapping']);pa=read(prior['admission'])
        require({k:v for k,v in g.items() if k!='rebased_origin'}=={k:v for k,v in pm.items() if k!='rebased_origin'},'Actual A0AC rebased geometry/checks differ')
        for k in ('mapping_checks','q0_paired_replay','caps_fixed','cap_comparison','old_mode_ids'):require(admission[k]==pa[k],'Actual A0AC admission differs: '+k)
        previous=np.load(check(prior['contact']['origin']),allow_pickle=False)
        require(np.array_equal(previous['positions'],union) and np.array_equal(previous['symbols'],np.asarray(union_symbols)),'Actual A0AC source-union arrays differ')
        def labels(x):
            if isinstance(x,dict):return {k:labels(v) for k,v in x.items() if k!='omitted_membership'}
            if isinstance(x,list):return [labels(v) for v in x]
            return 'unchanged_carver_cofactor_H' if x=='unchanged_carver_PQQ_H' else x
        require(labels(read(prior['contact']['identities'])['centers'])==labels(ids) and labels(pa['source_contacts'])==labels(contacts),'Actual A0AC identities/contact checks differ')
        require(prior['contact']['moving_indices']==moved_qm and prior['contact']['caps']==cap_order,'Actual A0AC moving membership differs')
        for t in prior['tasks']:
            require(t['xyz']['sha256']==r['matrix'][t['metal']]['amide_flip']['xyz']['sha256'] and t['active_indices']==oldtasks[t['metal']]['active_indices'] and t['active_mode_ids']==oldtasks[t['metal']]['active_mode_ids'],'Actual A0AC origin alias differs')
        regression={'status':'passed','source_case_id':REUSE_ID,'exact_rebased_context_and_checks':True,'exact_source_union_arrays':True,
          'exact_cap_q0_and_contacts':True,'prepared_alias':prior['admission'],'executed_collection':c['completed_search_collection'],'new_energy_calls':0}
        write(d/'A0AC_REGRESSION.json',regression)
        row.update(mapping=prior['mapping'],contact=prior['contact'],prepared_state_alias=prior['admission'],regression=pin(d/'A0AC_REGRESSION.json'),
          reused_prepared_tasks=prior['tasks'],executed_result_alias=c['reuse']['executed_row'],repeat_search_permitted=False)
        for t in row['task_templates']:
            prior_task,=[p for p in prior['tasks'] if p['metal']==t['metal']]
            t.update(xyz=prior_task['xyz'],mapping=prior['mapping'],contact=prior['contact'],new_evaluation_permitted=False)
    return row

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before numerical imports')
    mp={'path':str(Path(path).resolve()),'sha256':sha};c=read(mp)
    require(c['protocol']==PROTOCOL and c['implementation']['PREPARE.py']==pin(__file__) and c['python']==pin(sys.executable),'Frozen implementation/interpreter differs')
    require(c['declared_sources']==25 and len(c['sources'])==25,'Canonical denominator differs')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    rows=[base_row(s) for s in c['sources']];receipt={'status':'failed','manifest':mp,'job_id':os.environ['SLURM_JOB_ID'],'new_energy_calls':0,'new_optimizations':0}
    try:
        for name in ('geometry','geometry_manifest','scope_manifest','canonical_collection','prior_preparation','prior_config','prior_implementation','completed_search_collection','completed_search_execution','source_preparation_helper','source_parser','generic_contact_helper','python','worker_manifest'):check(c[name])
        for p in c['implementation'].values():check(p)
        for k,p in c['technical_recovery'].items():
            if k!='kind':check(p)
        require(c['geometry']['sha256']==GEOMETRY_SHA and c['settings']==read(c['worker_manifest'])['settings'],'Actual geometry/unchanged settings differ')
        require([s['row'] for s in c['sources']]==read(c['geometry'])['rows'] and [s['scope'] for s in c['sources']]==read(c['scope_manifest'])['rows'],'Actual canonical inputs/order differ')
        require(read(c['completed_search_execution'])['status']=='complete','Actual reused search execution failed')
        base=load('_canonical_rebase_source_parser',c['source_parser']);receipt['placement']=base.allocation()
        source_helper=load('_canonical_rebase_original_flip',c['source_preparation_helper']);generic=load('_canonical_rebase_generic_contact',c['generic_contact_helper'])
        sys.path.insert(0,c['source_snapshot'])
        import numpy as np
        from mace_site_kinematics import Kinematics
        import adaptive_angular_proposals as angular
        for module in tuple(sys.modules.values()):
            path=getattr(module,'__file__',None)
            if path and Path(path).resolve().parent==Path(c['source_snapshot']).resolve():
                require(Path(path).name in c['source_modules'] and pin(path)==c['source_modules'][Path(path).name],'Actual imported source pin differs: '+str(path))
        for k in ('maximum_physical_heavy_displacement_A','angular_bound_radian','final_physical_tolerance_A','angular_tolerance_radian','physical_boundary_margin_A','angular_boundary_margin_radian'):
            require(angular.SETTINGS[k]==c['settings'][k],'Frozen geometry bound differs: '+k)
        first=next(i for i,s in enumerate(c['sources']) if s['row']['source_case_id']==REUSE_ID)
        for i in (first,*[j for j in range(25) if j!=first]):
            s=c['sources'][i];d=out/s['row']['source_case_id'];d.mkdir()
            try:rows[i]=prepare_row(s,d,c,np,Kinematics,angular,source_helper,base,generic)
            except Exception as e:
                rows[i].update(reason=f'{type(e).__name__}: {e}',traceback=traceback.format_exc())
                if (d/'ADMISSION.json').is_file():rows[i]['candidate_admission']=pin(d/'ADMISSION.json')
            write(d/'PREPARATION.json',rows[i])
            if i==first:
                require(rows[i]['status']=='prepared' and read(rows[i]['regression'])['status']=='passed','Exact A0AC rebased regression failed; remaining preparation held')
                receipt['A0AC_regression']=rows[i]['regression']
        receipt['status']='complete' if all(r['status']=='prepared' for r in rows) else 'partial'
    except BaseException as e:
        receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'status':receipt['status'],'declared_sources':25,
          'prepared_sources':sum(r['status']=='prepared' for r in rows),'rows':rows,'execution':pin(out/'EXECUTION.json'),
          'new_energy_calls':0,'new_optimizations':0,'force_gradient_qualification':'No new qualification; actual A0AC saved evidence aliased',
          'new_search_execution_authorized':False,'new_calibration_or_scores':False})
    return 0 if receipt['status']=='complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='op',required=True)
    a=sub.add_parser('freeze');a.add_argument('--output',required=True)
    a=sub.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--manifest-sha',required=True)
    args=p.parse_args()
    if args.op=='freeze':freeze(args.output)
    else:sys.exit(execute(args.manifest,args.manifest_sha))
