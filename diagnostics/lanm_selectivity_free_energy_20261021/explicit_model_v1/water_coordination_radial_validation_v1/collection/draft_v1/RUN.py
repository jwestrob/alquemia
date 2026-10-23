#!/usr/bin/env python3
"""Saved-only four-cell/twelve-check comparison; frozen helper arithmetic."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;BRANCH=HERE.parent;ROOT=HERE.parents[4]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/water_coordination_radial_validation_v1'
GIDS=tuple('aquo9_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))
FUTURE=('model_bound','model_result_path','reference_release','reference_result_inventory')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    return pin(path)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def default_config():
    r=json.loads((BRANCH/'reference/READY.json').read_text());m=json.loads((BRANCH/'model/CONFIG_SOURCE.json').read_text())
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(BRANCH/'ROOT_DECLARATION.md'),
        'reference_source':r['implementation'],'reference_configuration':r['configuration'],'reference_manifest':r['prepared_manifest'],
        'reference_source_review':pin(WORK/'reference/independent_review_v1/REVIEW.json'),
        'model_source':m['implementation'],'model_source_configuration':pin(BRANCH/'model/CONFIG_SOURCE.json'),
        'model_scorer_review':pin(WORK/'independent_model_scorer_review_v1/REVIEW.json'),
        'score_source':pin(BRANCH/'scorer/SCORE.py'),'score_specification':pin(BRANCH/'scorer/SPEC.md'),
        'geometry_manifest':m['geometry_manifest'],'body_helper':m['body_helper'],
        'cell_ids':list(GIDS),'declared_checks':12,'native_calls':0,**{k:None for k in FUTURE}}

def load(config):
    cfg=json.loads(Path(config).read_text());static=copy.deepcopy(cfg)
    for k in FUTURE:static[k]=None
    require(static==default_config(),'Exact frozen collector sources/configuration required')
    for value in static.values():
        if isinstance(value,dict) and set(value)=={'path','sha256'}:verify(value)
    rr=read(cfg['reference_source_review']);mr=read(cfg['model_scorer_review'])
    require(rr['status']=='source_preparation_review_passed' and mr['status']=='source_review_passed_unbound','Source reviews did not pass')
    require(mr['model_source']==cfg['model_source'] and mr['scorer_source']==cfg['score_source'],'Model/scorer review source differs')
    return cfg

def check_binding(cfg,bound,mcfg,release,geometry):
    require(set(bound['rows'])==set(geometry['rows'])==set(GIDS),'Exactly four bound geometry slots required')
    require(bound['candidate_result']==release['candidate_result']==mcfg['candidate_result'],'Model/reference candidate mismatch')
    require(bound['candidate_review']==release['selected_review']==mcfg['candidate_review'],'Model/reference independent selected review mismatch')
    require(bound['parameters']==release['parameters'],'Model/reference frozen parameters differ')
    require(bound['geometry_manifest']==cfg['geometry_manifest']==mcfg['geometry_manifest'],'Model/reference geometry manifest differs')
    require(release['prepared_manifest']==cfg['reference_manifest'] and release['configuration']==cfg['reference_configuration'] and release['implementation']==cfg['reference_source'],'Reference release source/preparation differs')
    require(bound['implementation']==cfg['model_source'] and mcfg['implementation']==cfg['model_source'] and mcfg['body_helper']==cfg['body_helper'],'Model/body source differs')
    for gid in GIDS:
        b=bound['rows'][gid];g=geometry['rows'][gid]
        require(b['geometry']==g['geometry'] and b['bodies']==g['bodies'],'Bound geometry/body differs: '+gid)
        require(b['files']['xyz']['sha256']==g['input']['xyz']['sha256'],'Bound coordinate bytes differ: '+gid)
        for key in ('types','bonds_zero_based'):require(b['input'][key]==g['input'][key],'Bound topology differs: '+gid)

def unavailable_model(reason,bound,authority=None):
    return {'status':'saved_model_admission_partial','rows':{g:{'status':'unavailable','model':None,'reason':reason} for g in GIDS},'complete':0,
        'anchor_model':bound['anchor_model'],'new_native_calls':0,'raw_result':authority,'reason':reason}

def admit_model(cfg,bound,mcfg,model):
    path=Path(cfg['model_result_path']).resolve();authority=None;raw=None
    if not path.exists():
        model.bound_inputs(mcfg,bound)
        return unavailable_model('missing_result',bound),{},None
    try:
        authority=pin(path);raw=read(authority)
        require(raw['bound']==cfg['model_bound'] and raw['candidate_result']==bound['candidate_result'] and raw['candidate_review']==bound['candidate_review'],'Actual model candidate/binding differs')
        admitted=model.admit_result(mcfg,bound,raw);admitted['raw_result']=authority
    except Exception as exc:
        # Invalid native/result JSON cannot erase the four-cell denominator.
        model.bound_inputs(mcfg,bound);admitted=unavailable_model(repr(exc),bound,authority)
    cost={k:raw.get(k) for k in ('job_id','native_program_invocations','new_FD_energies','elapsed_seconds','parent_CPU_seconds','aggregate_child_user_CPU_seconds','aggregate_child_system_CPU_seconds','per_program_child_CPU_may_overlap_not_summed')} if isinstance(raw,dict) else {}
    return admitted,cost,authority

def compare(bound,reference,admitted,score,body):
    require(set(reference['rows'])==set(admitted['rows'])==set(GIDS) and reference['declared']==4,'Admission denominator differs')
    data={'rows':{}};models={}
    for gid in GIDS:
        b=bound['rows'][gid];r=reference['rows'][gid]
        require(r['geometry']==b['geometry'] and r['bodies']==b['bodies'],'Reference/model geometry/body mismatch: '+gid)
        for key in ('types','bonds_zero_based'):require(r['input'][key]==b['input'][key],'Reference/model topology mismatch: '+gid)
        require(r['input']['xyz']['sha256']==b['files']['xyz']['sha256'],'Reference/model coordinate bytes differ: '+gid)
        data['rows'][gid]={'geometry':b['geometry'],'input':b['input'],'bodies':b['bodies'],'reference':r['reference'] if r['status']=='complete' else None}
        models[gid]=admitted['rows'][gid]['model'] if admitted['rows'][gid]['status']=='complete' else None
    for metal in ('La','Dy'):
        for key in ('energy','forces'):require(reference['anchor_reference'][metal][key]==bound['anchor_reference'][metal][key],'Central reference anchor differs')
        require(admitted['anchor_model'][metal]==bound['anchor_model'][metal],'Central model anchor differs')
    anchors={'model':bound['anchor_model'],'reference':reference['anchor_reference']}
    result=score.evaluate(data,models,body,anchors)
    require(result['declared_checks']==12 and (not result['all_validation_checks_passed'] or result['available_checks']==12),'Scorer denominator/pass contradiction')
    return result

def report(result):
    s=result['score'];lines=['# Fresh W9 radial comparison','',f"Status: **{result['status']}**. Available checks: {s['available_checks']}/12; passed: {s['passed_checks']}/12.",'',
        '| Comparison | Status | Error or RMS / maximum | Pass |','|---|---|---:|---|']
    groups=[('work '+k,v,'energy') for k,v in s['individual_radial_work'].items()]+[('paired work '+k,v,'energy') for k,v in s['paired_radial_work'].items()]+[(f'force {sign} {metal}',v,'force') for sign,d in s['rigid'].items() for metal,v in d.items()]
    for name,v,kind in groups:
        value='unavailable' if v['status']!='complete' else (f"{v['error_kcal_mol']:.6f} kcal/mol" if kind=='energy' else f"{v['RMS_kcal_mol_A']:.6f} / {v['maximum_component_kcal_mol_A']:.6f} kcal/mol/Å")
        lines.append(f"| {name} | {v['status']} | {value} | {v['passed']} |")
    stages=result['reference_admission']['rows'];changed=[g for g,r in stages.items() if r['reference'] is not None and r['reference']['changed_central_energy_stage_signature']]
    lines+=['',f"Changed central-to-displaced native stage signatures: {', '.join(changed) if changed else 'none among admitted references'}.",
        'Raw SCF, final-grid corrections, electronic-state fields and central stage metadata remain in RESULT.json and REFERENCE_ADMISSION.json.',
        '',f"New native calls in this collection: 0. Saved reference allocated CPU-seconds: {result['reference_admission']['allocated_native_CPU_seconds']}. Model measured costs are retained without summing overlapping per-program CPU measurements.",
        '', 'This compares fresh radial transfer against the unchanged provisional frozen-f ECP-PBE0-D4 reference. It does not qualify electronic states, analytic-force derivatives, equilibrium hydration, bulk/protein transfer or exchange free energies. Missing checks cannot pass. No fitting or production release follows.','']
    return '\n'.join(lines)

def collect(config,out):
    start,cpu=time.monotonic(),time.process_time();require(all(os.environ.get(k)=='1' for k in THREADS),'Saved collection requires one local CPU thread')
    cfg=load(config);require(all(cfg[k] is not None for k in FUTURE),'Actual root-bound inputs intentionally unavailable')
    bound=read(cfg['model_bound']);mcfg=read(bound['config']);release=read(cfg['reference_release']);geometry=read(cfg['geometry_manifest'])
    check_binding(cfg,bound,mcfg,release,geometry)
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    ref=module(cfg['reference_source'],'radial_collect_reference');model=module(cfg['model_source'],'radial_collect_model');score=module(cfg['score_source'],'radial_collect_score');body=module(cfg['body_helper'],'radial_collect_body')
    # The existing helper owns the native release gate and raw reference replay.
    rp=ref.collect(cfg['reference_configuration']['path'],cfg['reference_manifest']['path'],verify(cfg['reference_release']),verify(cfg['reference_result_inventory']),out/'REFERENCE_ADMISSION.json')
    refs=read(rp);admitted,cost,mp=admit_model(cfg,bound,mcfg,model);ap=save(out/'MODEL_ADMISSION.json',admitted)
    score_result=compare(bound,refs,admitted,score,body)
    status='validation_unavailable' if score_result['available_checks']<12 else ('validation_passed' if score_result['all_validation_checks_passed'] else 'validation_failed')
    result={'status':status,'implementation':pin(__file__),'configuration':pin(config),'candidate_result':bound['candidate_result'],'selected_review':bound['candidate_review'],'parameters':bound['parameters'],
        'declared_cells':4,'score':score_result,'reference_admission':refs,'reference_admission_file':rp,'model_admission':admitted,'model_admission_file':ap,'raw_model_result':mp,'model_actual_cost':cost,
        'new_native_calls':0,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'physical_accuracy_qualified':False,'automatic_fitting_release':False,'production_released':False}
    resultpin=save(out/'RESULT.json',result);(out/'REPORT.md').write_text(report(result))
    print(json.dumps({'result':resultpin,'report':pin(out/'REPORT.md'),'status':status}));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('config');a.add_argument('--out',required=True)
    a=s.add_parser('collect');a.add_argument('--config',required=True);a.add_argument('--out',required=True)
    args=p.parse_args()
    if args.action=='config':print(json.dumps(save(args.out,default_config())))
    else:collect(args.config,args.out)
