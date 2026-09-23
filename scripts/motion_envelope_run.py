"""Finite envelope native origins and bounded proposals; scalar components stay absent."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from compact_solvation_compare import native_endpoint
from adaptive_force_diagnostic import project,preview
from mace_site_kinematics import Kinematics
import adaptive_angular_proposals as angular
import consistent_context as context
import motion_envelope as envelope
import slsqp_precision as precision
import union_adaptive as original
from union_triple_transfer_run import selector_equivalence

PROTOCOL='Nikasha_motion_envelope4p3_native_only_four_mode_ftol1e8_proposals_v1'


def collect_origins(stage,output):
    stage=Path(stage);ready=read_json(stage/'READY.json');inputs=read_json(verify(ready['inputs']));mm=read_json(verify(ready['MACE_manifest']))
    rows=[]
    for source in inputs['cases']:
        tasks=[t for t in mm['tasks'] if (t['source_case_id'],t['selection_id'])==(source['case_id'],source['selection_id'])]
        row={'case_id':tasks[0]['case_id'],'source_case_id':source['case_id'],'selection_id':source['selection_id'],
            'source':source,'mapping':source['maps'],'source_preparation':source['representations']['context']['preparation'],
            'native_endpoints':{},'solvent_endpoints':None,'composite_status':'unavailable_not_evaluated'}
        for t in tasks:
            z=t['metal'];path=verify(ready['MACE_manifest']).parent/'results'/t['task_id']/'result.json'
            try:
                ep=native_endpoint(read_json(path),mm['model']);original.origin(ep,mm['model'])
                if not context.reusable_state(source['representations']['context']['endpoints'][z],ep):raise InvalidArtifact('actual origin physical state differs')
                row['native_endpoints'][z]={**ep,'status':'complete'}
            except (OSError,InvalidArtifact,KeyError) as exc:row['native_endpoints'][z]={'status':'unavailable','reason':str(exc)}
        row['status']='native_complete' if len(row['native_endpoints'])==2 and all(ep['status']=='complete'for ep in row['native_endpoints'].values()) else 'unavailable'
        rows.append(row)
    result={'protocol_id':envelope.PROTOCOL,'ready':record(stage/'READY.json'),'inputs':ready['inputs'],'inventory':inputs['inventory'],
        'model':mm['model'],'rows':rows,'denominator':34,'native_complete':sum(r['status']=='native_complete'for r in rows),
        'composite_complete':0,'new_calls_in_collection':0,'implementation':record(__file__)}
    write_new(output,result);return {'collection':record(output),'native_complete':result['native_complete'],'composite_complete':0}


def build_tasks(row,model):
    cid=row['case_id'];points={};projections={};tasks=[]
    for z in ('Ca','La'):
        ep=row['native_endpoints'][z];native,forces=original.origin(ep,model)
        kin=Kinematics(read_json(verify(row['mapping'][z]))['context']);zero=np.zeros(len(kin.modes))
        coords,v,raw,normed,_=project(kin.data,zero,forces)
        points[z]={'status':'complete','coordinate':ep['xyz'],'MACE':ep['native_MACE_receipt'],'MACE_eV':native['energy_eV'],
            'forces':native['forces'],'full_q':zero.tolist(),'active_q_radian':[0.]*4,
            'source_receipt_format':'actual_legacy_native_OMOL','reused_scientific_origin':True}
        projections[z]=(kin,v,raw,normed)
        tasks.append({'task_id':cid+'__'+z,'case_id':cid,'source_case_id':row['source_case_id'],'selection_id':row['selection_id'],
            'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],'mapping':row['mapping'][z],
            'source_preparation':row['source_preparation'],'mode_count':len(kin.modes),'q0_status':'native_only_composite_unavailable',
            'q0':{'components':{'MACE_eV':native['energy_eV']},'low':None,'native_MACE_receipt':ep['native_MACE_receipt']}})
    choice=preview(projections['Ca'][0].modes,projections['Ca'][1],projections['Ca'][3],projections['La'][3])
    ids=[x['id']for x in choice['selected']]
    if len(ids)!=4:raise InvalidArtifact('four independent physical modes unavailable')
    indices=[[m['id']for m in projections['Ca'][0].modes].index(mid)for mid in ids]
    for t in tasks:
        z=t['metal'];point={**points[z],'gradient_kcal_mol_rad':projections[z][2][indices].tolist()}
        t.update(active_indices=indices,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,selector=choice,
            origin_reuse={'proposal_receipt':points[z]['MACE'],'point':point,'receipt_kind':'actual_archived_native_origin_not_a_proposal'})
        angular.final_geometry(projections[z][0],t,np.zeros(4),[r[0]for r in xyz(verify(t['xyz']))])
    paired(verify(tasks[1]['xyz']),verify(tasks[0]['xyz']),tasks[1]['charge'],tasks[0]['charge'])
    if read_json(verify(tasks[0]['mapping']))!=read_json(verify(tasks[1]['mapping'])):raise InvalidArtifact('paired physical map changed')
    return tasks,choice


def prepare(origins,template,agreement,output):
    d=read_json(origins);base=read_json(template);audit=read_json(verify(d['inventory']))
    if d['denominator']!=34 or d['model']!=base['model'] or audit['counts']['pilot_denominator']!=34:raise InvalidArtifact('declared34/model differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[];missing=[]
    for r in d['rows']:
        if r['status']!='native_complete':missing.append({'case_id':r['case_id'],'reason':'actual native origins unavailable'});continue
        try:ts,choice=build_tasks(r,base['model'])
        except InvalidArtifact as exc:missing.append({'case_id':r['case_id'],'reason':str(exc)});continue
        tasks.extend(ts);cases.append({'case_id':r['case_id'],'status':'prepared','origin_row':r,'union':r['source'],'selection':choice})
    m={k:base[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    m.update(protocol_id=PROTOCOL,settings=precision.SETTINGS,agreement=record(agreement),origin_collection=record(origins),
        source=record(template),cases=cases,tasks=tasks,declared_case_ids=[c['case_id']for c in cases],
        declared_all_sources=[r['case_id']for r in d['rows']],unavailable=missing,maximum_optimizer_starts=len(tasks),
        declared_maximum_optimizer_starts=68,maximum_cross_MACE_calls=0,maximum_GFN2_calls=0,new_origin_calls=0,new_DFT_calls=0,
        production_changed=False,reference=None,composite_status='unavailable_not_evaluated',
        optimizer_software={'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
        execution_adapter=record(__file__),implementation=original.snapshot(out/'implementation'))
    write_new(out/'manifest.json',m);v=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest);d=read_json(verify(m['origin_collection']));base=read_json(verify(m['source']))
    if (m['protocol_id']!=PROTOCOL or m['settings']!=precision.SETTINGS or m['declared_all_sources']!=[r['case_id']for r in d['rows']]
        or len(m['declared_all_sources'])!=34 or m['declared_maximum_optimizer_starts']!=68 or m['maximum_GFN2_calls']!=0
        or m['maximum_cross_MACE_calls']!=0 or m['reference'] is not None or len(m['tasks'])!=2*len(m['cases'])
        or set(m['declared_case_ids'])|{r['case_id']for r in m['unavailable']}!=set(m['declared_all_sources'])):
        raise InvalidArtifact('finite native-only scope or policy differs')
    verify(m['agreement']);verify(m['execution_adapter'])
    for p in m['implementation'].values():verify(p)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software changed')
    for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources'):
        if m[k]!=base[k]:raise InvalidArtifact('runtime or checkpoint changed')
    rows={r['case_id']:r for r in d['rows']};actual=[]
    for c in m['cases']:
        r=rows[c['case_id']]
        if c['origin_row']!=r or c['union']!=r['source']:raise InvalidArtifact('actual source row differs')
        ts,selection=build_tasks(r,m['model']);selector_equivalence(c['selection'],selection);actual.extend(ts)
    for frozen,current in zip(m['tasks'],actual):
        selector_equivalence(frozen['selector'],current['selector']);current['selector']=frozen['selector']
        if set(frozen['q0']['components'])!={'MACE_eV'} or frozen['q0']['low'] is not None:raise InvalidArtifact('native proposal must not invent scalar components')
    if actual!=m['tasks']:raise InvalidArtifact('actual native force task differs')
    return {'status':'validated','manifest':record(manifest),'declared_sources':34,'prepared_sources':len(m['cases']),
        'bounded_searches':len(m['tasks']),'new_q0_calls':0,'GFN2_calls':0,'cross_MACE_calls':0,'composite_available':False}


def engine():
    u,_=precision.engine();u.PROTOCOL=PROTOCOL;u.validate=validate
    return u

def execute(manifest):return engine().execute(manifest)
def collect(manifest,output):return engine().collect(manifest,output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for op,keys in {'collect_origins':('stage','output'),'prepare':('origins','template','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps(globals()[cmd](**a),indent=2))
