"""Frozen charge-feature ablation inventory and gates in the existing OMOL runner."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import accepted_attempt, check_atoms, EV_TO_KCAL
from mace_omol_ablation import ADAPTER, COMPONENT
from mace_omol_intact import EVALUATION, collect as collect_endpoints
from mace_omol_panel import ADAPTER_TASK

PROTOCOL='mace_omol_intact_charge_feature_ablation_descriptor_v1'
SEMANTICS='energy_like_descriptor_not_quantum_endpoint'
SOURCES={
    'core_collection':'020fb97916c75f9ea7f8f140f6ae4f9791fbcd7221e2436643111f2732a25dc7',
    'numerical_collection':'af17dcaead1584940bffbc666d5df811d909dd8f2565ebbfcc59f23abfb77961',
    'source_report':'7fdf787ccd81eeea798ad96f77d0febedb3d253b99b968dc8fa0c661b37ded82',
    'spectator_collection':'a5cb255530a270a3077adbb68627b3532f644984665dc8f5ddfb828645d58cbc'}
SETTINGS={'core_calls':8,'numerical_calls':14,'other_primary_calls':16,'spectator_calls':4,
          'numerical_tolerance_model_kcal':.01,'direction_margin_model_kcal':.02,
          'spectator_tolerance_model_kcal':.1,'output_semantics':SEMANTICS}


def source_task(result):
    mp=verify(result['manifest']);m=read_json(mp)
    for pin in m['implementation'].values():verify(pin)
    verify(m['software']);verify(m['model']['checkpoint'])
    found=[t for t in m['tasks'] if t['task_id']==result['task_id']]
    if len(found)!=1:raise InvalidArtifact('source task missing or duplicated')
    t=found[0]
    if not any(accepted_attempt(a,t,mp)==result for a in
               sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*'))):
        raise InvalidArtifact('source endpoint lacks its actual matching receipt')
    if not t.get('energy_only') or not t.get('capture_native_readout'):
        raise InvalidArtifact('pinned energy-only source required')
    return t,m


def expected(sources):
    data={}
    for name,sha in SOURCES.items():
        pin=sources[name]
        if pin['sha256']!=sha:raise InvalidArtifact('declared source selection changed: '+name)
        data[name]=read_json(verify(pin))
    tasks=[];parents=[]
    def add(result,group,mode='batched'):
        source,parent=source_task(result);parents.append(parent)
        t=copy.deepcopy(source)
        for k in ('cache_key','edge_adapter','execution_device'):t.pop(k,None)
        t.update(task_id=source['task_id']+'_mask_'+mode,ablation_group=group,
                 original_endpoint={'manifest':result['manifest'],'task_id':source['task_id']},
                 energy_component=COMPONENT,charge_feature_adapter=ADAPTER,
                 output_semantics=SEMANTICS)
        if mode=='batched':t['edge_adapter']=copy.deepcopy(ADAPTER_TASK)
        tasks.append(t)
    for result in data['core_collection']['rows'].values():
        for mode in ('native','batched'):add(result,'core',mode)
    for result in data['numerical_collection']['rows'].values():add(result,'numerical')
    source=data['source_report']
    for name,case in source['scores'].items():
        if name=='ALPHA_1F6S':continue
        for metal in ('La','Ca'):
            for position in ('bound','detached'):add(case['endpoints'][metal][position],'primary')
    for result in data['spectator_collection']['rows'].values():
        if result['task_id'].startswith('GGR_1GLG_'):add(result,'spectator')
    counts={group:sum(t['ablation_group']==group for t in tasks)
            for group in ('core','numerical','primary','spectator')}
    if counts!={'core':8,'numerical':14,'primary':16,'spectator':4}:
        raise InvalidArtifact('fixed 42-forward inventory changed')
    if len({t['task_id'] for t in tasks})!=42:raise InvalidArtifact('duplicate descriptor task')
    parent=parents[0]
    if any(any(p[k]!=parent[k] for k in ('software','inventory')) for p in parents):
        raise InvalidArtifact('source software/inventory differ')
    return tasks,parent,source


def descriptor_model(software):
    from mace_omol import model
    m=model(software)
    m.update(energy_component=COMPONENT,charge_feature_adapter=ADAPTER,output_semantics=SEMANTICS,
             execution_adapter={'selection':'native_core_and_exact_batched_tasks','batched':ADAPTER_TASK})
    return m


def component_check(path,implementation):
    c=read_json(path)
    for key in ('checkpoint','adapter','implementation','agreement','features'):verify(c[key])
    from mace_omol import CHECKPOINT_SHA
    if (c['status']!='pass' or c['checkpoint']['sha256']!=CHECKPOINT_SHA
            or c['adapter']['sha256']!=implementation['mace_omol_ablation.py']['sha256']
            or not all(c[k] is True for k in ('exact_reference_projection_match','all_201_charges_identical',
                                             'parameter_versions_unchanged'))
            or c['model_energy_forward_calls']!=0 or c['embedding_block_rows']!=201):
        raise InvalidArtifact('actual matching charge adapter component qualification required')
    return c


@cached_file_checks
def prepare(core_collection,numerical_collection,source_report,spectator_collection,component,agreement,output):
    from mace_omol import common,seal
    sources={name:record(path) for name,path in zip(SOURCES,(core_collection,numerical_collection,source_report,spectator_collection))}
    tasks,parent,_=expected(sources)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,'ablation_development')
    c=component_check(component,m['implementation'])
    if c['agreement']!=m['agreement']:raise InvalidArtifact('component check belongs to another plan')
    m.update(protocol_id=PROTOCOL,sources=sources,component_qualification=record(component),
             model=descriptor_model(verify(parent['software'])),settings=SETTINGS,
             energy_evaluation=EVALUATION,tasks=tasks,output_semantics=SEMANTICS)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL
    m=read_json(manifest);tasks,parent,_=expected(m['sources'])
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage']!='ablation_development'
            or m['settings']!=SETTINGS or m['energy_evaluation']!=EVALUATION or m['tolerances']!=TOL
            or m['output_semantics']!=SEMANTICS or m['reused'] or m['compute_budget'] is not None
            or m['model']!=descriptor_model(verify(m['software']))
            or m['software']!=parent['software'] or m['inventory']!=parent['inventory']):
        raise InvalidArtifact('charge ablation model/inventory policy changed')
    s=read_json(verify(m['software']))
    for pin in [m['agreement'],*m['implementation'].values(),s['python'],s['requirements'],s['checkpoint'],
                s['inspection'],s['download_receipt'],*read_json(verify(s['backend_source_inventory']))['files']]:verify(pin)
    c=component_check(verify(m['component_qualification']),m['implementation'])
    if c['agreement']!=m['agreement']:raise InvalidArtifact('ablation plan changed')
    if len(m['tasks'])!=42:raise InvalidArtifact('descriptor task count changed')
    for t,expected_task in zip(m['tasks'],tasks):
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if payload!=expected_task:raise InvalidArtifact('descriptor physical state or feature policy changed')
        if check_atoms(xyz(verify(t['xyz'])),t['charge'])!=t['state']:
            raise InvalidArtifact('descriptor electron state changed')
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('descriptor cache identity differs')
    return {'status':'pass','tasks':42,'manifest':record(manifest)}


def numerical_checks(m,rows):
    checks={'core':[],'numerical':[]}
    for group in checks:
        tasks=[t for t in m['tasks'] if t['ablation_group']==group]
        if not all(rows[t['task_id']]['status']=='computed' for t in tasks):continue
        for t in tasks:
            error=rows[t['task_id']]['native_readout']['component_sum_error_kcal_mol']
            checks[group].append({'name':t['task_id']+'_accounting','error_model_kcal':error,'pass':abs(error)<=.01})
        def value(case,metal,position,variant,mode='batched'):
            found=[t for t in tasks if t['case_id']==case and t['metal']==metal and t['position']==position
                   and t['variant']==variant and t['task_id'].endswith('_mask_'+mode)]
            if len(found)!=1:raise InvalidArtifact('ambiguous numerical endpoint')
            return rows[found[0]['task_id']]['energy_eV']
        def add(name,error):checks[group].append({'name':name,'error_model_kcal':error,'pass':abs(error)<=.01})
        if group=='core':
            for case in ('1H4I','4MAE'):
                errors={metal:(value(case,metal,'bound','primary')-value(case,metal,'bound','primary','native'))*EV_TO_KCAL for metal in ('La','Ca')}
                for metal,error in errors.items():add(case+'_'+metal+'_native_batched',error)
                add(case+'_paired_native_batched',errors['Ca']-errors['La'])
        else:
            for variant in ('repeat','rotate','farther'):
                errors={}
                for position in (('detached',) if variant=='farther' else ('bound','detached')):
                    for metal in ('La','Ca'):
                        error=(value('ALPHA_1F6S',metal,position,variant)-value('ALPHA_1F6S',metal,position,'primary'))*EV_TO_KCAL
                        errors[metal,position]=error;add(metal+'_'+position+'_'+variant,error)
                    add('paired_'+position+'_'+variant,errors['Ca',position]-errors['La',position])
                if variant!='farther':add('score_'+variant,errors['Ca','bound']-errors['La','bound']-errors['Ca','detached']+errors['La','detached'])
    return checks


@cached_file_checks
def execution_gate(manifest,gate):
    m=read_json(manifest);result=collect_endpoints(manifest)
    checks=numerical_checks(m,result['rows'])
    groups=('core',) if gate=='core' else ('core','numerical')
    return all(checks[g] and all(c['pass'] for c in checks[g]) for g in groups)


@cached_file_checks
def collect(manifest):
    m=read_json(manifest);result=collect_endpoints(manifest);rows=result['rows']
    checks=numerical_checks(m,rows);scores={}
    source=read_json(verify(m['sources']['source_report']))
    for case,original in source['scores'].items():
        terms={};endpoints={}
        for metal in ('La','Ca'):
            endpoints[metal]={}
            for position in ('bound','detached'):
                found=[t for t in m['tasks'] if t['case_id']==case and t['metal']==metal
                       and t['position']==position and t['variant']=='primary' and t['kind']=='full']
                if len(found)!=1:raise InvalidArtifact('missing/duplicate primary descriptor state')
                endpoints[metal][position]=rows[found[0]['task_id']]
            pair=endpoints[metal]
            terms[metal]=pair['bound']['energy_eV']-pair['detached']['energy_eV'] if all(r['status']=='computed' for r in pair.values()) else None
        value=(terms['Ca']-terms['La'])*EV_TO_KCAL if all(v is not None for v in terms.values()) else None
        scores[case]={'R_mask_model_kcal':value,'original_R_coord_kcal_mol':original['R_coord_kcal_mol'],
                      'bound_minus_detached_model_eV':terms,'endpoints':endpoints,'evidence':original['evidence'],
                      'preparation':original['preparation']}
    contrasts=[]
    for old in source['contrasts']:
        a,b=(scores[old[k]]['R_mask_model_kcal'] for k in ('positive_case','negative_case'))
        value=a-b if a is not None and b is not None else None
        contrasts.append({'positive_case':old['positive_case'],'negative_case':old['negative_case'],
                          'delta_model_kcal':value,'original_delta_R_kcal_mol':old['delta_R_kcal_mol'],
                          'pass':value is not None and value>SETTINGS['direction_margin_model_kcal']})
    spectator={}
    for metal in ('La','Ca'):
        pair={t['position']:rows[t['task_id']] for t in m['tasks'] if t['ablation_group']=='spectator' and t['metal']==metal}
        spectator[metal]=pair['bound']['energy_eV']-pair['detached']['energy_eV'] if all(r['status']=='computed' for r in pair.values()) else None
    svalue=(spectator['Ca']-spectator['La'])*EV_TO_KCAL if all(v is not None for v in spectator.values()) else None
    original=scores['GGR_1GLG']['R_mask_model_kcal'];shift=svalue-original if svalue is not None and original is not None else None
    numeric=result['numerical_gate_pass'] and all(checks[g] and all(c['pass'] for c in checks[g]) for g in checks)
    spectator_pass=shift is not None and abs(shift)<=SETTINGS['spectator_tolerance_model_kcal']
    result.update(protocol_id=PROTOCOL,output_semantics=SEMANTICS,score_unit='kcal_equivalent_model_units',
                  numerical_checks=checks,numerical_gate_pass=bool(numeric),scores=scores,contrasts=contrasts,
                  spectator={'R_mask_model_kcal':svalue,'change_model_kcal':shift,'pass':spectator_pass},
                  predictive_development_gate_pass=all(c['pass'] for c in contrasts),
                  canonical_extension_permitted=bool(numeric and spectator_pass and all(c['pass'] for c in contrasts)),
                  baseline_changed=False,production_promotion=False,broad_affinity_validated=False)
    return result


@cached_file_checks
def report(manifest,output):
    validate(manifest);result=collect(manifest)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    result['report_implementation']=record(__file__)
    write_new(out/'result.json',result)
    lines=['# Charge-feature ablation development result','',
           'Modified MACE descriptors in kcal-equivalent model units; not quantum energies or binding free energies.',
           'All cases were previously consumed. Baseline/default unchanged.','',
           '| Case | Masked score | Original native score (kcal/mol) |','|---|---:|---:|']
    for name,row in result['scores'].items():lines.append(f"|{name}|{row['R_mask_model_kcal']}|{row['original_R_coord_kcal_mol']}|")
    lines+=['',f"Numerical gate: {result['numerical_gate_pass']}. Three-direction gate: {result['predictive_development_gate_pass']}.",
            f"Distant-sodium change: {result['spectator']['change_model_kcal']}; consistency pass: {result['spectator']['pass']}.",
            f"Conditional canonical extension permitted: {result['canonical_extension_permitted']}.",'',
            'See result.json for unrounded components, actual receipts, unavailable endpoints and checks.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'status':result['status'],'result':record(out/'result.json')}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for key in (*SOURCES,'component','agreement','output'):p.add_argument('--'+key.replace('_','-'),required=True)
    p=sub.add_parser('report')
    for key in ('manifest','output'):p.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    print(json.dumps({'prepare':prepare,'report':report}[command](**args),indent=2))
