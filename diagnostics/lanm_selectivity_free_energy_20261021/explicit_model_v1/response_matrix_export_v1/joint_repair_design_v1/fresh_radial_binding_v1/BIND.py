#!/usr/bin/env python3
"""Four existing radial geometries; one future reviewed joint candidate."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
IDS=tuple('aquo9_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))
FUTURE=('fit_result','fit_review','regression_result','regression_report','root_acceptance')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,v):
    p=Path(p);t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');os.replace(t,p);return pin(p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sources():return {k:pin(HERE/(k+'.py')) for k in ('BIND','MODEL','REFERENCE','COLLECT')}

def prepare(out):
    tick,cpu=time.monotonic(),time.process_time();require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread')
    inp=read(pin(HERE/'SOURCE_INPUTS.json'));old=module(inp['old_reference'],'radial_old_reference')
    rc,ref,original,geometry,cells,centers=old.load(verify(inp['old_reference_config']))
    rm=read(inp['reference_manifest']);require(rm['configuration']==inp['old_reference_config'] and rm['native_panel']==inp['reference_panel'],'Frozen original reference preparation')
    ec=read(inp['evaluator_template']);ev=module(inp['evaluator'],'radial_evaluator');em=read(ec['manifest']);mp=read(ec['moment_preparation'])
    moment=module(ec['moment_adapter'],'radial_moment');field=module(inp['field_helper'],'radial_field');body=module(ec['body_helper'],'radial_body')
    require(ec['implementation']==inp['evaluator'] and mp['configuration']['field_helper']==inp['field_helper'],'Frozen helper lineage')
    rows={}
    for gid in IDS:
        g=geometry['rows'][gid];metal=g['geometry']['metal'];central=em['rows']['aquo9_central_'+metal]
        require(g['input']['types']==central['types'],'Original typed atom order')
        rr=[s.split() for s in verify(g['input']['xyz']).read_text().splitlines()[1:] if s.strip()]
        xyz=[[float(x) for x in a[2:5]] for a in rr]
        require(xyz==g['geometry']['positions_A'] and [int(a[5]) for a in rr]==central['types'],'Literal prepared geometry')
        row=copy.deepcopy(central)
        for k in ('baseline','baseline_authority','prepared_initial_input','initial_render','original_inventory_row'):row.pop(k,None)
        row.update(id=gid,positions_A=xyz,source_input=dict(central['source_input'],xyz=g['input']['xyz']),needs_export=False,needs_gradient=True,energy_replay_kcal_mol=None,baseline_u_available=False)
        masses=row['masses_amu'];row['COM_A']=[sum(m*p[j] for m,p in zip(masses,xyz))/sum(masses) for j in range(3)]
        multipoles=moment.multipoles(verify(row['source_input']['parameters']).read_text(),verify(row['source_input']['key']).read_text())
        frames=[moment.frame_signature(rr,i,multipoles) for i in range(28)]
        row['permanent_lab_Debye']=field.permanent_moment({'frame_signatures':frames})
        require(body.bodies(g)==g['bodies'],'Original body graph')
        rows[gid]=dict(geometry=g['geometry'],input=g['input'],bodies=g['bodies'],runtime=row,frames=frames,parent='aquo9_central_'+metal)
    m=dict(status='prepared_fresh_radial_binding_unbound',sources=sources(),inputs=pin(HERE/'SOURCE_INPUTS.json'),specification=pin(HERE/'SPEC.md'),
           rows=rows,geometry_manifest=rc['geometry_manifest'],reference_manifest=inp['reference_manifest'],reference_panel=inp['reference_panel'],
           anchor_reference=rm['anchor_reference'],evaluator_manifest=ec['manifest'],units=em['units'],native_programs=dict(reference=4,model=8),
           workers=4,reference_MPI_ranks=4,new_FD_energies=0,declared_checks=12,new_native_calls=0,fresh_outcomes_read=False)
    out=Path(out);out.mkdir(parents=True,exist_ok=False);manifest=save(out/'MANIFEST.json',m)
    cfg=dict(implementation=pin(__file__),manifest=manifest,**{k:None for k in FUTURE})
    save(out/'CONFIG_BIND_TEMPLATE.json',cfg)
    return save(out/'RESULT.json',dict(status=m['status'],manifest=manifest,new_native_calls=0,fresh_outcomes_read=False,local_threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu))

def load(cfg):
    require(cfg['implementation']==pin(__file__) and set(cfg)=={'implementation','manifest',*FUTURE},'Exact binder configuration')
    m=read(cfg['manifest']);require(m['sources']==sources() and m['inputs']==pin(HERE/'SOURCE_INPUTS.json') and m['specification']==pin(HERE/'SPEC.md'),'Frozen binding sources')
    require(set(m['rows'])==set(IDS) and m['declared_checks']==12,'Four/twelve denominator')
    inp=read(m['inputs']);ec0=read(inp['evaluator_template']);ev=module(inp['evaluator'],'radial_current_evaluator')
    return m,inp,ec0,ev

def candidate(cfg):
    m,inp,template,ev=load(cfg);require(all(isinstance(cfg[k],dict) for k in FUTURE),'Actual candidate/regression/root acceptance intentionally unbound')
    fr=read(cfg['fit_result']);review=read(cfg['fit_review']);reg=read(cfg['regression_result']);accept=read(cfg['root_acceptance']);verify(cfg['regression_report'])
    require(fr['status']=='bounded_joint_response_fit_complete' and fr['candidate_frozen'] is True and fr['error'] is None and fr['original_EF_checks_passed'] is True,'Terminal frozen all25-passing candidate required')
    fc=read(fr['configuration']);ec=read(fc['evaluator_config']);selected=read(fr['selected'])
    require(fc['evaluator']==inp['evaluator']==ec['implementation'] and ec['manifest']==m['evaluator_manifest'],'Selected evaluator lineage')
    for k in ev.SOURCES:require(ec[k]==template[k],'Changed evaluator source '+k);verify(ec[k])
    em=read(ec['manifest']);domain=read(ec['domain']);values=ev.checked(selected['parameters'],domain)
    require(fr['selected']==fr['vectors'][fr['selected_key']]['artifact'] and selected['status']=='complete' and len(selected['cells'])==46 and set(selected['cells'])==set(em['rows']),'Selected46 denominator')
    require(all(r['status']=='complete' for r in selected['cells'].values()) and len(selected['score']['weighted_residuals'])==531 and selected['score']['original_EF_checks_passed'] is True,'Selected531/all25 admission')
    require(selected['score']['objective']==fr['selected_objective'],'Selected objective join')
    require(review['status']=='passed_terminal_joint_fit_review_continuation' and review['aggregate_audit_status']=='passed_using_preserved_completed_phases_and_narrow_continuation' and review['fit_result']==cfg['fit_result'] and review['selected']==fr['selected'] and review['original_EF_checks_passed'] is True,'Actual aggregate selected raw review')
    require(review['original_audit_passed'] is False and review['original_completed_check_count']==327825 and review['selected_objective']==fr['selected_objective'],'Preserved failed audit/completed phases')
    original=read(review['original_configuration']);require(original['driver_config']==fr['configuration'] and original['fit_result']==cfg['fit_result'],'Original review/driver configuration')
    require(review['implementation']==inp['candidate_review_source'] and review['driver_implementation']==fc['implementation'] and review['evaluator_implementation']==inp['evaluator'],'Frozen aggregate-review sources')
    for k in ('implementation','failed_review','strict_continuation_failure','configuration'):verify(review[k])
    require(read(review['failed_review'])['status']=='unavailable_terminal_joint_fit_review' and read(review['strict_continuation_failure'])['status']=='unavailable_terminal_joint_fit_review_continuation','Both failed audit attempts preserved')
    raw=read(review['selected_raw_review']);require(raw['status']=='passed_selected_raw_replay' and raw['selected']==fr['selected'] and set(raw['rows'])==set(em['rows']) and len(raw['rows'])==46,'Full46 selected replay denominator')
    require(all(r['status'] in ('passed_selected_raw_replay','inherited_from_completed_initial_baseline_audit') for r in raw['rows'].values()),'Every selected raw row admitted')
    require(all(r.get('regression_energy_moment_subset_exact') is True for r in raw['rows'].values() if r['status']=='passed_selected_raw_replay') and review['regression_dependency']['all_selected_energy_moment_subsets_exact'] is True,'Exact selected energy/moment dependency replay')
    require(review['new_native_calls']==0 and review['new_parameter_evaluations']==0 and review['new_optimization_calls']==0 and review['reserved_outcomes_opened'] is False,'No new science or held outcomes in audit')
    require(reg['status']=='consumed_trimer_regression_complete' and reg['implementation']==inp['regression_source'] and reg['error'] is None and reg['fit_result']==cfg['fit_result'] and reg['selected']==fr['selected'] and reg['parameters']==values and len(reg['cycles'])==8 and all(r['status']=='complete' for r in reg['cycles'].values()),'Complete same-candidate consumed eight cycles');verify(reg['implementation'])
    require(accept['status']=='accepted_joint_candidate_for_fresh_W9_radial' and all(accept[k]==cfg[k] for k in ('fit_result','fit_review','regression_result','regression_report')) and accept['selected']==fr['selected'] and accept['all25_passed'] is True and accept['complete_consumed_cycles']==8 and accept['candidate_frozen_before_fresh_outcomes'] is True,'Root scientific progression acceptance')
    build=read(ec['build_receipt']);q=read(ec['qualification_review'])
    require(build['status']=='complete' and build['manifest']==ec['source_manifest'] and all(build['executables'][k]==v for k,v in ec['executables'].items()),'Completed exact native build')
    require(q['status']=='passed_local_response_native_qualification' and all(q[k]==ec[k] for k in ('source_manifest','build_receipt','executables')) and all(q[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Actual native qualification')
    for p in [ec['source_manifest'],q['implementation'],q['result'],*ec['executables'].values()]:verify(p)
    h={k:ev.module(ec[k],'radial_'+k) for k in ev.HELPERS};proxy=SimpleNamespace(h=h,m=em,config=ec)
    proxy.render=lambda vals,cid,exports=None:ev.render(vals,em['rows'][cid],domain,h['resolver'],exports)
    for metal in ('La','Dy'):
        cid='aquo9_central_'+metal;r=selected['cells'][cid];row=em['rows'][cid]
        if r['inherited_native']:
            require(ev.effective(values,metal,domain)==ev.effective(em['initial_vector'],metal,domain) and all(r[k]==v for k,v in em['baseline_cells'][cid].items()),'Inherited central effective-input authority')
        else:
            ev.Evaluator._saved_inputs(proxy,r,row,values,Path(r['receipts'][0]['cwd']).parent)
            parsed=ev.Evaluator._parse(proxy,r,row)
            require(parsed['model']==r['model'] and selected['model'][cid]==r['model'],'Selected central raw E/F replay')
    return m,inp,ev,ec,em,fr,selected,values,domain,proxy

def bind(config,out):
    cfg=json.loads(Path(config).read_text());m,inp,ev,ec,em,fr,selected,values,domain,proxy=candidate(cfg)
    out=Path(out);out.mkdir(parents=True,exist_ok=False);rows={}
    for gid,meta in m['rows'].items():
        row=meta['runtime'];key,render=ev.render(values,row,domain,proxy.h['resolver'],exports=False);central=selected['cells'][meta['parent']]
        expected=ev.render(values,em['rows'][meta['parent']],domain,proxy.h['resolver'],exports=False)[0]
        require(key==expected,'Displaced versus selected central rendered key')
        if not central['inherited_native']:
            require(all(verify(central['program_inputs'][p]['key']).read_text()==key for p in ('analyze','testgrad')),'Exact selected central key bytes')
            require(all(central['program_inputs'][p]['parameters']['sha256']==row['source_input']['parameters']['sha256'] for p in ('analyze','testgrad')),'Exact selected central parameter bytes')
        folder=out/'inputs'/gid;folder.mkdir(parents=True)
        for k,n in (('xyz','system.xyz'),('parameters','amoeba09.prm')):(folder/n).write_bytes(verify(row['source_input'][k]).read_bytes())
        (folder/'system.key').write_text(key);files={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        rows[gid]=dict(meta,files=files,render=render,central_inherited_zero=central['inherited_native'])
    b=dict(status='fresh_W9_radial_bound_no_calls',implementation=pin(__file__),sources=sources(),config=pin(config),manifest=cfg['manifest'],
           candidate_result=cfg['fit_result'],candidate_review=cfg['fit_review'],selected=fr['selected'],parameters=values,root_acceptance=cfg['root_acceptance'],
           evaluator_config=read(fr['configuration'])['evaluator_config'],geometry_manifest=m['geometry_manifest'],rows=rows,
           anchor_model={metal:selected['model']['aquo9_central_'+metal] for metal in ('La','Dy')},anchor_reference=m['anchor_reference'],
           new_native_calls=0,fresh_outcomes_read=False)
    save(out/'BOUND.json',b);return b

def bound_context(bp):
    b=read(bp);cfg=read(b['config']);ctx=candidate(cfg);m,inp,ev,ec,em,fr,s,v,domain,proxy=ctx
    require(b['status']=='fresh_W9_radial_bound_no_calls' and b['implementation']==pin(__file__) and b['sources']==sources() and b['manifest']==cfg['manifest'],'Bound source/manifest')
    require(b['candidate_result']==cfg['fit_result'] and b['candidate_review']==cfg['fit_review'] and b['root_acceptance']==cfg['root_acceptance'] and b['selected']==fr['selected'] and b['parameters']==v and b['evaluator_config']==read(fr['configuration'])['evaluator_config'],'Bound candidate joins')
    require(set(b['rows'])==set(IDS) and b['geometry_manifest']==m['geometry_manifest'] and b['anchor_reference']==m['anchor_reference'] and b['anchor_model']=={metal:s['model']['aquo9_central_'+metal] for metal in ('La','Dy')},'Bound denominator/anchors')
    for gid,row in b['rows'].items():
        require(all(row[k]==x for k,x in m['rows'][gid].items()),'Frozen radial metadata')
        key,render=ev.render(v,row['runtime'],domain,proxy.h['resolver'],exports=False)
        require(verify(row['files']['key']).read_text()==key and row['render']==render,'Exact candidate key')
        for k in ('xyz','parameters'):require(verify(row['files'][k]).read_bytes()==verify(row['runtime']['source_input'][k]).read_bytes(),'Literal radial input')
    return b,ctx

def release(bp,rp,kind):
    b,ctx=bound_context(bp);r=read(rp);m=ctx[0]
    require(r['status']=='approved_fresh_W9_radial_'+kind and r['bound']==bp and r['sources']==sources() and r['candidate_result']==b['candidate_result'] and r['selected']==b['selected'] and r['root_acceptance']==b['root_acceptance'],'Exact arm/candidate release')
    review=read(r['source_review']);require(review['status']=='passed_fresh_radial_binding_source_review' and review['sources']==sources() and review['manifest']==b['manifest'] and review['inputs']==m['inputs'],'Independent frozen binder review')
    require(r['cell_ids']==list(IDS) and r['automatic_retry'] is False,'Finite released denominator')
    if kind=='model':require(r['maximum_native_programs']==8 and r['workers']==4 and r['new_FD_energies']==0,'Eight analytic model programs')
    else:require(r['maximum_reference_calls']==4 and r['ranks']==4 and set(r['task_output_directories'])=={'ecp_pbe0__'+g for g in IDS} and len(set(r['task_output_directories'].values()))==4,'Four reference tasks')
    return b,ctx,r

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('prepare');a.add_argument('--out',required=True)
    a=s.add_parser('bind');a.add_argument('--config',required=True);a.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(prepare(a.out) if a.action=='prepare' else {'status':bind(a.config,a.out)['status']}))
