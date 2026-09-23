"""Execute the declared triple-transfer manifest with existing scientific kernels."""
from __future__ import annotations
import argparse
import copy
import json
import math
import os
from pathlib import Path
import time
import numpy as np
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from affordable_workflow import execute as run_native
from compact_solvation_compare import native_endpoint
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
import adaptive_angular_proposals as angular
import consistent_context as context
import compact_solvation as solvent
import union_adaptive as original
import slsqp_precision as precision
from slsqp_precision_transfer import qualification,low_prepare
import union_triple_transfer as transfer
from union_triple_adaptive import PROTOCOL,POOL_PROTOCOL

SELECTOR_DIAGNOSTIC_COPY_TOLERANCE=1e-12
SELECTOR_DIAGNOSTIC_FIELDS=('Ca_kcal_mol_A','La_kcal_mol_A','differential_kcal_mol_A','independent_geometric_fraction')


def selector_equivalence(frozen,recomputed):
    """Only replayed diagnostic rounding may differ; retain the frozen selector."""
    a=copy.deepcopy(frozen);b=copy.deepcopy(recomputed);maximum=0.
    if len(a['selected'])!=len(b['selected']):raise InvalidArtifact('selected subspace differs')
    for x,y in zip(a['selected'],b['selected']):
        for key in SELECTOR_DIAGNOSTIC_FIELDS:
            u=x.pop(key);v=y.pop(key)
            if not isinstance(u,float) or not isinstance(v,float) or not math.isfinite(u) or not math.isfinite(v):
                raise InvalidArtifact('selector diagnostic type/nonfinite differs')
            delta=abs(u-v);maximum=max(maximum,delta)
            if delta>SELECTOR_DIAGNOSTIC_COPY_TOLERANCE:raise InvalidArtifact('selector diagnostic exceeds copy tolerance')
    if a!=b:raise InvalidArtifact('selected subspace or nondiagnostic selector data differs')
    return maximum


def execute_origins(manifest):
    transfer.validate_origins(manifest)
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=65536:
        raise InvalidArtifact('declared32CPU/64GiB allocation required')
    return run_native(manifest)


def collect_origins(stage,output):
    start=time.monotonic();stage=Path(stage);ready=read_json(stage/'READY.json');a=read_json(verify(ready['audit']))
    mp=verify(ready['origin_MACE_manifest']);mm=read_json(mp);lp=verify(ready['origin_GFN2_manifest'])
    transfer.validate_origins(lp);rows=[]
    for entry in a['rows']:
        if entry['pool_reuse']:continue
        cid=entry['pair_id'];c=entry['source'];prior=entry['origin_reuse']
        r={'case_id':cid,'source_case_id':entry['case_id'],'selection_id':entry['selection_id'],'source':c,
            'mapping':entry['mapping'],'native_endpoints':{},'solvent_endpoints':{}}
        for z in ('Ca','La'):
            ep=c['representations']['context']['endpoints'][z]
            try:
                reused=z in prior['native_endpoints']
                pin=prior['native_endpoints'][z]['receipt'] if reused else record(mp.parent/'results'/(cid+'__context__'+z)/'result.json')
                native=native_endpoint(read_json(verify(pin)),mm['model'])
                if not context.reusable_state(ep,native):raise InvalidArtifact('origin actual native state differs')
                original.origin(native,mm['model'])
                r['native_endpoints'][z]={**native,'status':'complete','reused':reused}
            except (OSError,InvalidArtifact,KeyError) as exc:r['native_endpoints'][z]={'status':'unavailable','reason':str(exc)}
            r['solvent_endpoints'][z]={}
            for medium in ('vacuum','alpb'):
                try:
                    pin=prior['solvent_endpoints'].get(z,{}).get(medium);reused=bool(pin)
                    if not pin:pin=solvent.completed(lp,cid+'__context__'+z+'__'+medium)
                    if pin is None:raise InvalidArtifact('origin solvent unavailable')
                    lm=read_json(verify(pin['manifest']));lt=next(t for t in lm.get('all_tasks',lm['tasks']) if t['task_id']==pin['task_id'])
                    if not context.reusable_state(ep,lt):raise InvalidArtifact('origin low-level state differs')
                    audit=solvent.diagnostics(pin,lt)
                    if audit['charge_sanity_status']!='pass':raise InvalidArtifact('actual native electronic state audit failed')
                    r['solvent_endpoints'][z][medium]={**pin,'status':'complete','reused':reused,'component_audit':audit}
                except (OSError,InvalidArtifact,KeyError) as exc:r['solvent_endpoints'][z][medium]={'status':'unavailable','reason':str(exc)}
        r['status']='complete' if all(r['native_endpoints'][z]['status']=='complete' and all(r['solvent_endpoints'][z][s]['status']=='complete' for s in ('vacuum','alpb')) for z in ('Ca','La')) else 'unavailable'
        r['source_preparation']=c['representations']['context']['preparation'];rows.append(r)
    result={'protocol_id':transfer.PROTOCOL,'ready':record(stage/'READY.json'),'audit':ready['audit'],'rows':rows,
        'denominator':55,'complete':sum(r['status']=='complete' for r in rows),'new_calls_in_collection':0,'wall_seconds':time.monotonic()-start}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def build_tasks(row,model):
    cid=row['case_id'];points={};projections={};new=[]
    for z in ('Ca','La'):
        ep=row['native_endpoints'][z];native,forces=original.origin(ep,model)
        kin=Kinematics(read_json(verify(row['mapping'][z]))['context']);zero=np.zeros(len(kin.modes))
        coords,v,raw,normed,_=project(kin.data,zero,forces);low=row['solvent_endpoints'][z]
        points[z]={'status':'complete','coordinate':ep['xyz'],'MACE':ep['native_MACE_receipt'],'MACE_eV':native['energy_eV'],
            'forces':native['forces'],'full_q':zero.tolist(),'active_q_radian':[0.]*4,
            'source_receipt_format':'actual_legacy_native_OMOL','reused_scientific_origin':True}
        projections[z]=(kin,v,raw,normed)
        new.append({'task_id':cid+'__'+z,'case_id':cid,'source_case_id':row['source_case_id'],'selection_id':row['selection_id'],
            'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],'mapping':row['mapping'][z],
            'source_preparation':row['source_preparation'],'mode_count':len(kin.modes),'q0_status':'available',
            'q0':{'components':{'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],
                'GFN2_ALPB_hartree':low['alpb']['energy_hartree']},'low':low,'native_MACE_receipt':ep['native_MACE_receipt']}})
    choice=preview(projections['Ca'][0].modes,projections['Ca'][1],projections['Ca'][3],projections['La'][3])
    ids=[x['id'] for x in choice['selected']]
    if len(ids)!=4:raise InvalidArtifact('four independent physical modes unavailable')
    indices=[[x['id'] for x in projections['Ca'][0].modes].index(mid) for mid in ids]
    for t in new:
        z=t['metal'];point={**points[z],'gradient_kcal_mol_rad':projections[z][2][indices].tolist()}
        t.update(active_indices=indices,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,selector=choice,
            origin_reuse={'proposal_receipt':points[z]['MACE'],'point':point,'receipt_kind':'actual_archived_native_origin_not_a_proposal'})
        angular.final_geometry(projections[z][0],t,np.zeros(4),[x[0] for x in xyz(verify(t['xyz']))])
    paired(verify(new[1]['xyz']),verify(new[0]['xyz']),new[1]['charge'],new[0]['charge'])
    if read_json(verify(new[0]['mapping']))!=read_json(verify(new[1]['mapping'])):raise InvalidArtifact('paired physical maps differ')
    return new,choice


def prepare(origins,template,agreement,output):
    data=read_json(origins);a=read_json(verify(data['audit']));base=read_json(template)
    if data['denominator']!=55 or base['model']!=a['config']['model'] or base['software']!=a['config']['software']:
        raise InvalidArtifact('frozen55 source/runtime scope differs')
    ready=read_json(verify(data['ready']));oi=read_json(verify(read_json(verify(ready['origin_GFN2_manifest']))['inputs']))
    qualification(verify(oi['rank_qualification']));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);pins=[]
    expected=[r['pair_id'] for r in a['rows'] if not r['pool_reuse']]
    if [r['case_id'] for r in data['rows']]!=expected:raise InvalidArtifact('source-pair order differs')
    for shard in range(2):
        chosen=data['rows'][shard::2];tasks=[];cases=[];unavailable=[]
        for r in chosen:
            if r['status']!='complete':unavailable.append({'case_id':r['case_id'],'status':'origin_unavailable','origin_row':r});continue
            try:ts,choice=build_tasks(r,base['model'])
            except InvalidArtifact as exc:unavailable.append({'case_id':r['case_id'],'status':'physical_mapping_unavailable','reason':str(exc),'origin_row':r});continue
            tasks.extend(ts);cases.append({'case_id':r['case_id'],'source_case_id':r['source_case_id'],'selection_id':r['selection_id'],
                'status':'prepared','origin_row':r,'union':r['source'],'selection':choice,'atom_count':r['source']['atom_count']})
        m={k:base[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
        dest=out/('shard_'+str(shard))/'proposals';dest.mkdir(parents=True)
        m.update(protocol_id=PROTOCOL,settings=precision.SETTINGS,agreement=record(agreement),origin_collection=record(origins),audit=data['audit'],
            source=record(template),reference=a['reference'],rank_qualification=oi['rank_qualification'],shard=shard,
            cases=cases,tasks=tasks,declared_case_ids=[c['case_id'] for c in cases],declared_source_pairs=[r['case_id'] for r in chosen],
            unavailable=unavailable,maximum_optimizer_starts=2*len(cases),declared_maximum_optimizer_starts=2*len(chosen),
            maximum_cross_MACE_calls=2*len(cases),maximum_GFN2_calls=8*len(cases),new_origin_calls=0,new_DFT_calls=0,production_changed=False,
            population='all100_triples_125_pairs_shard_'+str(shard),optimizer_software={'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
            execution_adapter=record(__file__),
            implementation=original.snapshot(dest/'implementation'))
        mp=dest/'manifest.json';write_new(mp,m);check=validate(mp);write_new(dest/'PREFLIGHT.json',check);pins.append(record(mp));print(json.dumps(check),flush=True)
    result={'audit':data['audit'],'reference':a['reference'],'origin_collection':record(origins),'shards':pins,'sources_per_shard':[28,27],
        'new_molecular_calls':0,'counts':a['counts']};write_new(out/'READY.json',result);return result


def validate(manifest):
    m=read_json(manifest);a=read_json(verify(m['audit']));d=read_json(verify(m['origin_collection']));base=read_json(verify(m['source']))
    qualification(verify(m['rank_qualification']));ref=read_json(verify(m['reference']));verify(m['agreement']);verify(m['execution_adapter'])
    expected=d['rows'][m['shard']::2];ids=[r['case_id'] for r in expected]
    if (m['shard'] not in (0,1) or len(expected)!=(28 if m['shard']==0 else 27) or m['declared_source_pairs']!=ids or
        m['protocol_id']!=PROTOCOL or m['settings']!=precision.SETTINGS or m['reference']!=a['reference'] or
        ref['reference_id']!=transfer.REFERENCE_ID or m['production_changed'] or m['maximum_optimizer_starts']!=len(m['tasks']) or
        len(m['tasks'])!=2*len(m['cases']) or m['declared_case_ids']!=[c['case_id'] for c in m['cases']] or
        set(m['declared_case_ids'])|{c['case_id'] for c in m['unavailable']}!=set(ids)):
        raise InvalidArtifact('finite fixed triple shard/settings differ')
    for key in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources'):
        if m[key]!=base[key]:raise InvalidArtifact('actual runtime/model differs')
    for p in m['implementation'].values():verify(p)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer implementation changed')
    rows={r['case_id']:r for r in expected};tasks=[]
    for c in m['cases']:
        r=rows[c['case_id']]
        if c['origin_row']!=r or c['union']!=r['source'] or r['status']!='complete':raise InvalidArtifact('actual source row changed')
        ts,choice=build_tasks(r,m['model']);tasks.extend(ts)
        selector_equivalence(c['selection'],choice)
    for expected_task,actual_task in zip(m['tasks'],tasks):
        selector_equivalence(expected_task['selector'],actual_task['selector'])
        actual_task['selector']=expected_task['selector']
    if tasks!=m['tasks']:raise InvalidArtifact('actual q0 force/state/projection task changed')
    return {'status':'validated','manifest':record(manifest),'declared_sources':len(expected),'prepared_sources':len(m['cases']),
        'unavailable':len(m['unavailable']),'searches':len(tasks),'maximum_cross_MACE':2*len(m['cases']),'maximum_GFN2':8*len(m['cases'])}


def engine():
    p=precision.private('slsqp_precision');p.validate=validate;u,shared=p.engine()
    u.PROTOCOL=PROTOCOL;u.POOL_PROTOCOL=POOL_PROTOCOL;shared.low_prepare=low_prepare
    return u,shared


def execute(manifest):return engine()[0].execute(manifest)
def collect(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    specs={'execute_origins':('manifest',),'collect_origins':('stage','output'),'prepare':('origins','template','agreement','output'),
        'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),
        'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output')}
    for op,fields in specs.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f.replace('_','-'),required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
