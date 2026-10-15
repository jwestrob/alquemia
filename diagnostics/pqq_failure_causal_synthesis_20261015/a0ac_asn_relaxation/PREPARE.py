#!/usr/bin/env python3
"""Metadata freeze, then Slurm-only rebased-map preparation. No model calls."""
import argparse, copy, hashlib, importlib.util, json, os, shutil, socket, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation'
OLD=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_scoring_v1'
METALS=('Ca','La'); MOVED=('OD1','ND2','HD21','HD22'); PROTOCOL='A0AC_admitted_amide_basin_frozen_four_mode_relaxation_v1'
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
    old=read(OLD/'CONFIG.json');prep=read(old['preparation']);pm=read(old['preparation_manifest']);inv=read(old['inventory'])
    collection=read(OLD/'execution/FINAL_COLLECTION.json');require(collection['complete_new_cells']==20,'Actual complete flip matrix required')
    sources=[]
    for row in prep['rows']:
        src,=[r for r in inv['rows'] if r['source_case_id']==row['source_case_id']]; searches={}
        for z in METALS:
            rp=src['optimizer_evidence']['adaptive_'+z]['receipt'];r=read(rp);m=read(r['manifest']);t,=[t for t in m['tasks'] if t['task_id']==r['task_id']]
            require(r['status']=='proposal_available','Original search unavailable');searches[z]={'result':rp,'manifest':r['manifest'],'task':t}
        ip=Path(row['geometry_admission']['path']).with_name('CONTACT_IDENTITIES.json')
        sources.append({'row':row,'source':src,'old_searches':searches,'contact_identities':pin(ip)})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir()
    for name in ('PREPARE.py','PREPARE.sbatch','PLAN.md','RECOVERY_V2.md'):shutil.copy2(HERE/name,out/'implementation'/name)
    c={'protocol':PROTOCOL,'implementation':{n:pin(out/'implementation'/n) for n in ('PREPARE.py','PREPARE.sbatch')},'plan':pin(out/'implementation/PLAN.md'),
       'old_config':pin(OLD/'CONFIG.json'),'old_collection':pin(OLD/'execution/FINAL_COLLECTION.json'),'old_scorer':old['implementation'],
       'old_solvent_collection':pin(ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1/execution/FINAL_COLLECTION.json'),
       'old_preparation':old['preparation'],'old_preparation_manifest':old['preparation_manifest'],'inventory':old['inventory'],'sources':sources,
       'source_preparation_helper':pm['implementation']['PREPARE.py'],'source_parser':pm['source_parser_reference'],
       'python':old['python'],'source_ids':[s['row']['source_case_id'] for s in sources],
       'scout_source_id':sources[0]['row']['source_case_id'],'output':str(out/'preparation'),'new_energy_calls':0}
    previous=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v1'
    diagnostic=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_caps_v1/execution'
    c['technical_recovery']={'kind':'cap_like_for_like_exact_comparison_v2','note':pin(out/'implementation/RECOVERY_V2.md'),
        'preserved_config':pin(previous/'CONFIG.json'),'preserved_preparation':pin(previous/'preparation/PREPARATION.json'),
        'diagnostic_result':pin(diagnostic/'RESULT.json'),'diagnostic_execution':pin(diagnostic/'EXECUTION.json')}
    require(c['scout_source_id']=='a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-0','Declared first source changed')
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def prepare_row(s,d,c,b,scorer,base,source_helper,m,np,Kinematics,angular,project,versions):
    r=s['row'];sid=r['source_case_id'];src=s['source'];gate=read(r['geometry_admission']);cp=read(r['source_preparation']);state=read(r['state_signature'])
    require(gate['pass'] and r['status']=='prepared','Flipped source unprepared');read(s['contact_identities']);check(r['protonated_source'])
    maps={z:read(r['maps'][z]) for z in METALS};require(maps['Ca']==maps['La'],'Original paired maps differ');original=maps['Ca'];g=copy.deepcopy(original)
    before=np.asarray([a[1:] for a in b.xyz(check(r['matrix']['Ca']['amide_flip']['original_xyz']))]);flips={z:b.xyz(check(r['matrix'][z]['amide_flip']['xyz'])) for z in METALS}
    require(flips['Ca'][1:]==flips['La'][1:] and flips['Ca'][0][1:]==flips['La'][0][1:],'Flip pair differs');after=np.asarray([a[1:] for a in flips['Ca']]);symbols=[a[0] for a in flips['Ca']]
    entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']};source_qm={key(x['source']):i for i,x in entries.items() if 'source' in x}
    named=r['named_Asn_indices_zero_based'];moved={named[n] for n in MOVED};require(moved==set(r['moved_indices_zero_based']),'Moved identity differs')
    physical_qm={i:source_qm[key(x)] for i,x in enumerate(g['context']['source_atom_metadata']) if x and key(x) in source_qm}
    oldkin=Kinematics(original['context']);zero=np.zeros(len(oldkin.modes));old_core=oldkin.evaluate(zero)[1];require(np.max(np.abs(old_core-before))<=1e-12,'Original map mismatch')
    changed=[]
    for i,qm in physical_qm.items():
        if qm in moved:g['context']['positions_A'][i]=after[qm].tolist();changed.append(i)
    require(len(changed)==4,'Four source sites not mapped');g['context']['core_positions_A']=after.tolist()
    kin=Kinematics(g['context']);p,core,j,cj=kin.evaluate(zero);require(np.max(np.abs(core-after))<=1e-12,'Rebased map does not reproduce actual flip')
    require(all(g['context'][k]==original['context'][k] for k in original['context'] if k not in ('positions_A','core_positions_A')),'Map topology or cap Jacobian changed')
    fixed=[i for i in range(len(p)) if i not in changed];require(np.array_equal(kin.positions[fixed],oldkin.positions[fixed]),'Physical exterior changed during rebase')
    mapping_checks=kin.check(zero);require(mapping_checks['pass'],'Rebased analytic Jacobian/covalent check failed');g['checks']=mapping_checks
    g['rebased_origin']={'original_mapping':r['maps'],'admitted_flip':r['matrix'],'geometry_admission':r['geometry_admission'],'changed_physical_indices':changed,'new_domain_origin':True}
    d.mkdir();write(d/'mapping.json',g);mapping=pin(d/'mapping.json');caps=set(gate['cap_indices'])
    cap_order=sorted(caps)
    require(np.array_equal(after[cap_order],before[cap_order]),'Raw admitted cap coordinates changed')
    require(np.array_equal(kin.core[cap_order],oldkin.core[cap_order]),'Raw map cap coordinates changed')
    require(np.array_equal(core[cap_order],old_core[cap_order]),'Rebased computed caps changed')
    cap_comparison={'raw_before_after_exact':True,'raw_map_caps_exact':True,'old_new_computed_caps_exact':True,
        'old_cap_q0_residual_A':float(np.max(np.abs(old_core[cap_order]-before[cap_order]))),
        'new_cap_q0_residual_A':float(np.max(np.abs(core[cap_order]-after[cap_order]))),
        'unchanged_full_q0_tolerance_A':1e-12,'diagnostic':c['technical_recovery']['diagnostic_result']}
    ids=read(s['contact_identities'])['centers'];require(len(ids)==gate['contact_union_atom_count'],'Contact identity count changed')
    physical=source_helper.physical_source(check(r['protonated_source']),base);lookup={key(a['source']):a for a in physical}
    exterior=[lookup[key(a['source'])]['xyz'] for a in ids[204:]];union=np.vstack((after,np.asarray(exterior)));union_symbols=symbols+[lookup[key(a['source'])]['element'] for a in ids[204:]]
    require(len(union)==9546 and len(physical)==9529,'Frozen full-source inventory differs')
    np.savez_compressed(d/'contact_origin.npz',positions=union,symbols=np.asarray(union_symbols));write(d/'contact_identities.json',{'centers':ids})
    oldtasks={z:s['old_searches'][z]['task'] for z in METALS};require(oldtasks['Ca']['active_mode_ids']==oldtasks['La']['active_mode_ids'],'Old paired modes differ')
    active=oldtasks['Ca']['active_indices'];require(len(active)==4 and [kin.modes[i]['id'] for i in active]==oldtasks['Ca']['active_mode_ids'],'Frozen selected modes changed')
    moving={i for k in active for i in kin.modes[k]['moving_indices']};require(moving<=set(physical_qm),'Moving source atom missing from physical context')
    require(all(not (link[2] in moving or link[3] in moving) for link in kin.links if link[1]=='cap'),'Unsupported moving cap boundary')
    moved_qm=sorted({physical_qm[i] for i in moving});require(not caps.intersection(moved_qm),'Synthetic motion added')
    contact={'origin':pin(d/'contact_origin.npz'),'identities':pin(d/'contact_identities.json'),'moving_indices':moved_qm,'caps':sorted(caps),'source':r['protonated_source']}
    checks=source_helper.contact_check(np,union,union,union_symbols,moved_qm,ids,caps);require(checks['pass'],'New source baseline contact check failed')
    matrix=next(x for x in read(c['old_collection'])['rows'] if x['source_case_id']==sid)['matrix'];tasks=[]
    for z in METALS:
        oldsearch=read(s['old_searches'][z]['result']);oldmanifest=read(s['old_searches'][z]['manifest']);require(oldmanifest['settings']==m['settings'],'Old optimizer settings differ')
        require(oldsearch['status']=='proposal_available' and oldsearch['proposal']['MACE']==src['matrix'][z]['adaptive_'+z]['MACE'],'Original proposal join changed')
        cell=matrix[z]['amide_flip'];value=read(cell['receipt']);req=read(value['request']);t=copy.deepcopy(oldtasks[z])
        t.update(task_id=sid+'__amide_flip__'+z,case_id=sid+'__amide_flip',source_case_id=sid,xyz=r['matrix'][z]['amide_flip']['xyz'],mapping=mapping,
                 source_preparation=r['source_preparation'],original_source_preparation=r['source_preparation'],rebased_mapping=mapping,contact=contact,atom_count=204)
        scorer.admitted_result(value,req,t,m['model'],versions);require(value['energy_eV']==cell['energy_eV'],'Cached flip energy changed')
        force=np.load(check(value['forces']),allow_pickle=False);require(force.shape==(204,3) and np.isfinite(force).all(),'Cached force invalid')
        coords,vectors,raw,normed,lengths=project(kin.data,zero,force);require(np.max(np.abs(coords-after))<=1e-12,'Rebased projected origin geometry differs')
        gradient=-np.einsum('mij,ij->m',cj[active],force)*scorer.EV;require(np.array_equal(gradient,raw[active]),'Force/Jacobian contraction differs')
        point={'status':'complete','coordinate':t['xyz'],'MACE':cell['receipt'],'MACE_eV':value['energy_eV'],'forces':value['forces'],
               'full_q':zero.tolist(),'active_q_radian':[0.]*4,'gradient_kcal_mol_rad':gradient.tolist(),'reused_scientific_origin':True}
        t.update(q0_status='native_only_composite_unavailable',q0={'components':{'MACE_eV':value['energy_eV']},'native_MACE_receipt':cell['receipt'],'low':None},
                 origin_reuse={'proposal_receipt':cell['receipt'],'point':point,'receipt_kind':'actual_completed_flip_with_rebased_projection'},
                 frozen_selector_origin=s['old_searches'][z]['manifest'],selector_policy='retain_exact_original_four_IDs_no_reselection')
        angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in flips[z]]);tasks.append(t)
    write(d/'ADMISSION.json',{'status':'prepared','pass':True,'source_case_id':sid,'mapping_checks':mapping_checks,'original_geometry_admission':r['geometry_admission'],
       'q0_paired_replay':True,'caps_fixed':True,'cap_comparison':cap_comparison,'all_physical_H_retained':True,'source_contacts':checks,'old_mode_ids':oldtasks['Ca']['active_mode_ids'],
       'original_discrete_displacement_not_in_local_domain':True,'local_domain_reference':'admitted_flip','old_searches':s['old_searches']})
    return {'status':'prepared','source_case_id':sid,'conditioning':r['conditioning'],'sample':r['sample'],'case_key':r['case_key'],'record_id':r['record_id'],
        'mapping':mapping,'contact':contact,'admission':pin(d/'ADMISSION.json'),'tasks':tasks,'baseline':matrix,'original_row':r,'old_searches':s['old_searches']}

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before numerical imports');require(int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])>=1,'No CPU allocation')
    mp={'path':str(Path(path).resolve()),'sha256':sha};c=read(mp);require(pin(__file__)==c['implementation']['PREPARE.py'] and pin(sys.executable)==c['python'],'Wrong frozen preparation/interpreter')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();rows=[];receipt={'status':'failed','config':mp,'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'new_energy_calls':0}
    try:
        check(c['plan']);check(c['implementation']['PREPARE.sbatch']);check(c['old_solvent_collection'])
        recovery=c['technical_recovery'];require(recovery['kind']=='cap_like_for_like_exact_comparison_v2','Unsupported recovery')
        for k,v in recovery.items():
            if k!='kind':check(v)
        old=read(c['old_config']);scorer=load('_asn_prior_scorer',c['old_scorer']);m=read(old['worker_manifest']);b=load('_asn_warm',m['implementation']);b.parent_env()
        pf=scorer.preflight(old,b);write(out/'OLD_PREFLIGHT.json',pf)
        import numpy as np
        from mace_site_kinematics import Kinematics
        import adaptive_angular_proposals as angular
        from adaptive_force_diagnostic import project
        import slsqp_precision
        require(slsqp_precision.engine()[0].completion.SETTINGS==m['settings'],'Frozen completion settings differ')
        source_helper=load('_asn_flip_preparer',c['source_preparation_helper']);base=load('_asn_source_parser',c['source_parser'])
        require([s['row']['source_case_id'] for s in c['sources']]==c['source_ids'] and len(c['sources'])==10,'Finite source inventory changed')
        for s in c['sources']:
            try:rows.append(prepare_row(s,out/s['row']['source_case_id'],c,b,scorer,base,source_helper,m,np,Kinematics,angular,project,pf['versions']))
            except Exception as exc:rows.append({'source_case_id':s['row']['source_case_id'],'status':'unavailable','reason':str(exc),'traceback':traceback.format_exc()})
        status='complete' if len(rows)==10 and all(r['status']=='prepared' for r in rows) else 'partial'
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'config':mp,'status':status,'declared_sources':10,'declared_new_searches':20,'reused_old_searches':20,'reused_endpoint_cells':80,'rows':rows,'versions':pf['versions'],'new_energy_calls':0})
        receipt.update(status=status,preparation=pin(out/'PREPARATION.json'))
    except BaseException:receipt['traceback']=traceback.format_exc();raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu);write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);f=s.add_parser('freeze');f.add_argument('--output',required=True)
    e=s.add_parser('execute');e.add_argument('--config',required=True);e.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='freeze':freeze(a.output)
    else:sys.exit(execute(a.config,a.config_sha))
