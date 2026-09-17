"""Finite numerical qualification of the exact OMOL interaction batching adapter."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify
from mace_hybrid import EV_TO_KCAL
from mace_omol_edges import ADAPTER
from mace_omol_intact import qualified as native_qualified,validate as native_validate,collect as native_collect
from mace_file_checks import cached_file_checks

PROTOCOL='mace_omol_exact_edge_qualification_v1'
CONFIG={'id':ADAPTER,'edge_order':'original','allowed_chunk_sizes':[1024,2048],'gradients':False}
CPU_PROTOCOL='mace_omol_native_cpu_qualification_v1'
CPU_CONFIG={'id':'native_ScaleShiftMACE_CPU_float64','execution_adapter':None,'gradients':False}


def qualified(path):
    c=read_json(path);mp=verify(c['manifest']);m=read_json(mp)
    if m['stage'] not in ('edge_core','cpu_core'):raise InvalidArtifact('exact core qualification required')
    validate(mp);actual=collect(mp)
    if c!=actual or not c['numerical_gate_pass']:raise InvalidArtifact('edge core equivalence did not pass')
    return c,m


def tasks(native,stage):
    result=[]
    cpu=stage.startswith('cpu_')
    for chunk in ((None,) if cpu else (1024,2048) if stage=='edge_core' else (1024,)):
        for old in native['tasks']:
            t=copy.deepcopy(old);t.pop('cache_key',None)
            t['reference_task_id']=old['task_id']
            if cpu:t['execution_device']='cpu'
            else:t['edge_adapter']={'id':ADAPTER,'chunk_size':chunk}
            if stage.endswith('_core'):t['task_id']=old['task_id']+('_cpu' if cpu else '_chunk_'+str(chunk))
            result.append(t)
    return result


@cached_file_checks
def prepare(native_core_collection,intact_manifest,agreement,output,edge_qualification=None,backend='edge'):
    from mace_omol import common,seal
    _,core=native_qualified(native_core_collection,'intact_core')
    native_validate(intact_manifest);full=read_json(intact_manifest)
    if full['core_qualification']!=record(native_core_collection):raise InvalidArtifact('intact source used another native bridge')
    if backend not in ('edge','cpu'):raise InvalidArtifact('unsupported engineering backend')
    stage=backend+('_intact' if edge_qualification else '_core');cpu=backend=='cpu'
    _,out,m=common(verify(core['inventory']),verify(core['software']),agreement,output,stage)
    m.update(protocol_id=CPU_PROTOCOL if cpu else PROTOCOL,adapter=CPU_CONFIG if cpu else CONFIG,native_core_collection=record(native_core_collection),
             intact_manifest=record(intact_manifest),edge_qualification=record(edge_qualification) if edge_qualification else None)
    if cpu:m['model']['execution_device']='cpu'
    else:m['model']['execution_adapter']=CONFIG
    if edge_qualification:
        _,parent=qualified(edge_qualification)
        if (parent['native_core_collection']!=m['native_core_collection'] or parent['intact_manifest']!=m['intact_manifest']
                or parent['stage']!=backend+'_core'):
            raise InvalidArtifact('edge qualification sources differ')
    m['tasks']=tasks(full if edge_qualification else core,stage)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL,model
    m=read_json(manifest)
    cpu=m['stage'].startswith('cpu_');config=CPU_CONFIG if cpu else CONFIG
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=(CPU_PROTOCOL if cpu else PROTOCOL) or m['adapter']!=config
            or m['stage'] not in ('edge_core','edge_intact','cpu_core','cpu_intact') or m['tolerances']!=TOL or m['reused']):
        raise InvalidArtifact('edge engineering definition changed')
    _,core=native_qualified(verify(m['native_core_collection']),'intact_core')
    fullpath=verify(m['intact_manifest']);native_validate(fullpath);full=read_json(fullpath)
    if full['stage']!='intact_qualification' or full['core_qualification']!=m['native_core_collection']:
        raise InvalidArtifact('wrong native intact source')
    expected_model=model(verify(core['software']))
    if cpu:expected_model['execution_device']='cpu'
    else:expected_model['execution_adapter']=CONFIG
    if m['model']!=expected_model or m['software']!=core['software'] or m['inventory']!=core['inventory']:
        raise InvalidArtifact('edge adapter changed native scientific model/source')
    verify(m['agreement'])
    for ref in m['implementation'].values():verify(ref)
    if m['stage'].endswith('_intact'):
        _,parent=qualified(verify(m['edge_qualification']))
        if (parent['native_core_collection']!=m['native_core_collection'] or parent['intact_manifest']!=m['intact_manifest']
                or parent['stage']!=('cpu_core' if cpu else 'edge_core')):
            raise InvalidArtifact('edge parent belongs to other sources')
    elif m['edge_qualification'] is not None:raise InvalidArtifact('unexpected edge parent')
    is_core=m['stage'].endswith('_core');expected=tasks(core if is_core else full,m['stage'])
    if len(expected)!=len(m['tasks']) or len(expected)!=((4 if cpu else 8) if is_core else 14):
        raise InvalidArtifact('finite edge task count changed')
    for t,e in zip(m['tasks'],expected):
        if {k:v for k,v in t.items() if k!='cache_key'}!=e:raise InvalidArtifact('edge task differs from exact native input')
        verify(t['xyz'])
        if t['cache_key']!=cache_key({'task':e,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('edge adapter cache identity changed')
    return {'status':'pass','tasks':len(expected),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    from mace_omol import TOL
    m=read_json(manifest);c=native_collect(manifest)
    if c['status']=='complete' and m['stage'].endswith('_core'):
        original=read_json(verify(m['native_core_collection']));mp=verify(original['manifest'])
        if original!=native_collect(mp):raise InvalidArtifact('native bridge no longer matches actual receipts')
        errors={}
        for t in m['tasks']:
            error=(c['rows'][t['task_id']]['energy_eV']-original['rows'][t['reference_task_id']]['energy_eV'])*EV_TO_KCAL
            errors[(t['case_id'],t['metal'],t.get('edge_adapter',{}).get('chunk_size'))]=error
            c['checks'].append({'name':t['task_id']+'_native_equivalence','error_kcal_mol':error,
                                'pass':abs(error)<=TOL['energy_kcal_mol']})
        for case in ('1H4I','4MAE'):
            for chunk in ((None,) if m['stage']=='cpu_core' else (1024,2048)):
                error=errors[(case,'Ca',chunk)]-errors[(case,'La',chunk)]
                c['checks'].append({'name':case+'_paired_native_equivalence_'+str(chunk),
                                    'error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
    c['numerical_gate_pass']=c['status']=='complete' and bool(c['checks']) and all(r['pass'] for r in c['checks'])
    c['native_intact_equivalence_status']='separate_comparison_required' if m['stage'].endswith('_intact') else 'not_applicable'
    return c


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('native-core-collection','intact-manifest','agreement','output'):p.add_argument('--'+key,required=True)
    p.add_argument('--edge-qualification','--backend-qualification',dest='edge_qualification')
    p.add_argument('--backend',choices=('edge','cpu'),default='edge')
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
