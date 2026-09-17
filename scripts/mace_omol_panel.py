"""Whole-chain canonical OMOL tasks using the existing immutable task executor."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, xyz
from mace_hybrid import accepted_attempt, check_atoms, write_xyz
from mace_omol_intact import EVALUATION, PROTOCOL as INTACT_PROTOCOL, energy_task, geometry, collect as collect_endpoints
from mace_omol_panel_prepare import validate as prepared
from mace_omol_edge_run import CONFIG as EDGE_CONFIG, PRODUCT_CONFIG
from mace_omol_products import ADAPTER
from mace_file_checks import cached_file_checks

PROTOCOL='mace_omol_intact_canonical_class_v1'
DECISION={'calibration_count':25,'transfer_count':3,'minimum_gap_kcal_mol':0.02,
          'Ca_band':'R_coord<=maximum_Ca_calibration', 'La_band':'R_coord>=minimum_La_calibration',
          'inside_gap':'inconclusive','charge_range_exclusions':False}
ADAPTER_TASK={'id':ADAPTER,'chunk_size':1024,'product_chunk_size':1024}


def expected_tasks(data):
    from mace_omol import COMPONENT
    tasks=[]
    for row in data['rows']:
        if row['status']!='prepared':continue
        p=read_json(verify(row['preparation']))
        for position in ('bound','detached'):
            for metal in ('La','Ca'):
                endpoint=p['endpoints'][metal]
                t={'task_id':f'{row["case_id"]}_{metal}_{position}_primary',
                   'case_id':row['case_id'],'metal':metal,'kind':'full','position':position,'variant':'primary',
                   'metal_index':len(p['physical_atoms'])-1,'source_xyz':endpoint['xyz'],
                   'preparation':row['preparation'],'assembly':p['assembly'],'microstate':p['microstate'],
                   'explicit_waters':p['explicit_waters'],'evidence':p['evidence'],
                   'evaluation_role':row['evaluation_role'],'sequence_group':row['sequence_accession_group'],
                   'expected_class':row['expected_class'],'require_isolated_metal':position=='detached',
                   'energy_component':COMPONENT,'edge_adapter':ADAPTER_TASK,
                   'outside_reported_training_charge_range':row['outside_reported_training_charge_range'],
                   'outside_reported_training_size_range':row['outside_reported_training_size_range'],
                   **{k:v for k,v in endpoint.items() if k!='xyz'}}
                tasks.append(energy_task(t))
    return tasks


def product_qualification(path):
    from mace_omol_edge_report import verified
    result=verified(path)
    manifest=read_json(verify(read_json(verify(result['edge_collection']))['manifest']))
    if manifest['stage']!='product_intact' or manifest['adapter']!=PRODUCT_CONFIG:
        raise InvalidArtifact('completed exact product qualification required')
    return result,manifest


def reuse(data,source_collection):
    """Verify eight actual earlier endpoints against exact desired physical states.

    The earlier edge-only execution is explicitly reused after both adapters'
    native equivalence checks. It never satisfies an ordinary product task cache.
    """
    from mace_omol import model
    saved=read_json(source_collection);mp=verify(saved['manifest']);m=read_json(mp)
    expected_model=model(verify(m['software']));expected_model['execution_adapter']=EDGE_CONFIG
    if (m['protocol_id']!=INTACT_PROTOCOL or m['stage']!='intact_benchmark'
            or m['energy_evaluation']!=EVALUATION or m['model']!=expected_model
            or saved['status']!='complete' or not saved['numerical_gate_pass']):
        raise InvalidArtifact('completed compatible intact source required')
    for pin in m['implementation'].values():verify(pin)
    original={t['task_id']:t for t in m['tasks']};refs={};values={}
    for desired in expected_tasks(data):
        if desired['case_id'] not in ('1H4I','4MAE'):continue
        key='PQQ_'+desired['task_id'];source=original[key]
        for field in ('source_xyz','preparation','charge','spin_multiplicity','state','metal_index',
                      'position','variant','assembly','microstate','explicit_waters','energy_component',
                      'require_isolated_metal','energy_only','capture_native_readout'):
            if source[field]!=desired[field]:raise InvalidArtifact('declared crystal reuse differs: '+field)
        payload={k:v for k,v in source.items() if k!='cache_key'}
        if source['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('source task cache identity changed')
        actual=xyz(verify(source['xyz']));wanted=geometry(desired)
        if ([r[0] for r in actual]!=[r[0] for r in wanted]
                or np.max(np.abs(np.asarray([r[1:] for r in actual])-np.asarray([r[1:] for r in wanted])))>1e-12):
            raise InvalidArtifact('source crystal geometry is not the desired bound/detached state')
        found=[]
        for attempt in saved['attempts']:
            if attempt['task_id']==key and attempt['accepted']:
                verify(attempt['receipt'])
                value=accepted_attempt(Path(attempt['path']),source,mp)
                if value==saved['rows'][key]:found.append(value)
        if not found:raise InvalidArtifact('saved crystal endpoint lacks matching actual execution receipt')
        refs[desired['task_id']]={'source_collection':record(source_collection),'source_task_id':key}
        values[desired['task_id']]=found[-1]
    if len(refs)!=8:raise InvalidArtifact('all eight exact crystal reuses required')
    return refs,values,m


@cached_file_checks
def prepare(preparation,product_equivalence,source_collection,agreement,output):
    from mace_omol import common,seal
    from mace_omol_backbone_audit import inspect
    data=prepared(preparation)
    integrity=inspect(data)
    if any(a['status']!='pass' and r['status']=='prepared' for a,r in zip(integrity['rows'],data['rows'])):
        raise InvalidArtifact('legacy preparation contains unsupported peptide geometry; prepare with the version2 connectivity guard')
    eq,parent=product_qualification(product_equivalence)
    refs,_,old=reuse(data,source_collection)
    if parent['software']!=old['software'] or parent['inventory']!=data['source_inventory']:
        raise InvalidArtifact('qualified model/software/source inventory differs')
    _,out,m=common(verify(data['source_inventory']),verify(parent['software']),agreement,output,'panel_canonical')
    for name in ('mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('execution adapter differs from full qualification')
    m.update(protocol_id=PROTOCOL,energy_evaluation=EVALUATION,decision_policy=DECISION,
             preparation=record(preparation),product_equivalence=record(product_equivalence),
             source_collection=record(source_collection),reused=refs,
             unprepared_cases=[r for r in data['rows'] if r['status']!='prepared'],
             reuse_policy='exact_physical_states_and_separately_native_qualified_execution_adapters')
    m['model']['execution_adapter']=PRODUCT_CONFIG
    m['tasks']=[t for t in expected_tasks(data) if t['task_id'] not in refs]
    for t in m['tasks']:
        if t['position']=='bound':t['xyz']=t['source_xyz']
        else:
            path=out/(t['task_id']+'.xyz');write_xyz(path,geometry(t));t['xyz']=record(path)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL,model
    m=read_json(manifest)
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage']!='panel_canonical'
            or m['energy_evaluation']!=EVALUATION or m['decision_policy']!=DECISION or m['tolerances']!=TOL
            or m['reuse_policy']!='exact_physical_states_and_separately_native_qualified_execution_adapters'):
        raise InvalidArtifact('intact canonical method/acceptance changed')
    data=prepared(verify(m['preparation']));_,parent=product_qualification(verify(m['product_equivalence']))
    refs,_,old=reuse(data,verify(m['source_collection']))
    wanted_model=model(verify(m['software']));wanted_model['execution_adapter']=PRODUCT_CONFIG
    if (m['inventory']!=data['source_inventory'] or m['model']!=wanted_model
            or m['software']!=parent['software'] or m['software']!=old['software']
            or m['reused']!=refs or m['unprepared_cases']!=[r for r in data['rows'] if r['status']!='prepared']):
        raise InvalidArtifact('panel source/model/reuse changed')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    for name in ('mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('panel execution differs from product qualification')
    wanted=[t for t in expected_tasks(data) if t['task_id'] not in refs]
    if len(wanted)!=len(m['tasks']) or len(wanted)>104:raise InvalidArtifact('finite panel task inventory changed')
    for t,e in zip(m['tasks'],wanted):
        if {k:v for k,v in t.items() if k not in ('xyz','cache_key')}!=e:
            raise InvalidArtifact('panel scientific task differs from preparation')
        actual=xyz(verify(t['xyz']));coords=geometry(t)
        if actual!=coords or check_atoms(actual,t['charge'])!=t['state']:
            raise InvalidArtifact('exact panel coordinates/state changed')
        if t['position']=='bound' and t['xyz']!=t['source_xyz']:
            raise InvalidArtifact('bound source bytes changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('panel cache identity changed')
    return {'status':'pass','tasks':len(wanted),'reused_endpoints':len(refs),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    m=read_json(manifest);data=prepared(verify(m['preparation']))
    result=collect_endpoints(manifest)
    refs,values,_=reuse(data,verify(m['source_collection']))
    if refs!=m['reused']:raise InvalidArtifact('panel reuse changed')
    result['reused_rows']={key:{**ref,'result':values[key]} for key,ref in refs.items()}
    result['unprepared_cases']=m['unprepared_cases']
    if m['unprepared_cases']:
        result['status']='incomplete';result['numerical_gate_pass']=False
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('preparation','product-equivalence','source-collection','agreement','output'):
        p.add_argument('--'+key,required=True)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
