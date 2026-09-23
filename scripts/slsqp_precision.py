"""One isolated SLSQP precision test on four pinned UNION contexts."""
from __future__ import annotations
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import scipy
import adaptive_completion as original_completion
import adaptive_angular_proposals as angular
import union_adaptive as original_union
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired,HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
from accommodation_folds_compare import decision

PROTOCOL='union_four_angular_native_OMOL_SLSQP_ftol1e8_v1'
POOL_PROTOCOL='union_four_angular_ftol1e8_common_native_OMOL_GFN2_v1'
CASES=('q89gy2-pqq-la_model','q9z4j7-pqq-la_model','1H4I',
       'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1')
SETTINGS={**original_completion.SETTINGS,'optimizer_ftol':1e-8}
GATES={'proposal_energy_kcal_mol':.001,'GFN2_cell_work_kcal_mol':.1,'pooled_R_kcal_mol':.2}


def private(name):
    path=Path(__file__).with_name(name+'.py');key='_precision_private_'+name
    spec=importlib.util.spec_from_file_location(key,path);module=importlib.util.module_from_spec(spec)
    sys.modules[key]=module;spec.loader.exec_module(module);return module


def engine():
    completion=private('adaptive_completion');completion.SETTINGS=copy.deepcopy(SETTINGS)
    union=private('union_adaptive');pool=private('nikasha_pool')
    union.completion=completion;union.PROTOCOL=PROTOCOL;union.POOL_PROTOCOL=POOL_PROTOCOL
    union.SETTINGS=copy.deepcopy(SETTINGS);union.validate=validate;union.pool=pool
    pool.validate=union.validate_pool
    return union,pool


def prepare(canonical,pilot,canonical_pool,pilot_pool,reference,agreement,output):
    parents=[read_json(p) for p in (canonical,pilot)];collections=[read_json(p) for p in (canonical_pool,pilot_pool)]
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    tasks=[];cases=[];sources=[]
    for cid in CASES:
        i=0 if cid==CASES[0] else 1;p=parents[i]
        if p['protocol_id']!=original_union.PROTOCOL or p['settings']!=original_completion.SETTINGS:
            raise InvalidArtifact('original numerical policy differs')
        c=next(c for c in p['cases'] if c['case_id']==cid);cases.append(c)
        prior=next(c for c in collections[i]['cases'] if c['case_id']==cid)
        if prior['pool']['status']!='available':raise InvalidArtifact('original complete pool missing')
        for z in ('Ca','La'):
            t=next(t for t in p['tasks'] if (t['case_id'],t['metal'])==(cid,z));tasks.append(t)
        sources.append({'case_id':cid,'parent':record((canonical,pilot)[i]),'old_pool':record((canonical_pool,pilot_pool)[i])})
    p=parents[0];m={k:p[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources','optimizer_software')}
    m.update(protocol_id=PROTOCOL,settings=SETTINGS,gates=GATES,agreement=record(agreement),reference=record(reference),sources=sources,
             cases=cases,tasks=tasks,declared_case_ids=list(CASES),population='numerical_four',maximum_optimizer_starts=8,
             new_origin_calls=0,maximum_cross_MACE_calls=8,maximum_GFN2_calls=32,production_changed=False,new_DFT_calls=0,
             implementation=original_union.snapshot(out/'implementation'))
    path=out/'manifest.json';write_new(path,m);v=validate(path);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['gates']!=GATES or m['declared_case_ids']!=list(CASES):raise InvalidArtifact('precision policy or scope changed')
    if len(m['tasks'])!=8 or m['maximum_optimizer_starts']!=8 or m['maximum_cross_MACE_calls']!=8 or m['maximum_GFN2_calls']!=32:raise InvalidArtifact('finite scope differs')
    for k in ('agreement','reference','software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software changed')
    if {k for k in SETTINGS if SETTINGS[k]!=original_completion.SETTINGS[k]}!={'optimizer_ftol'}:raise InvalidArtifact('more than tolerance changed')
    expected=[]
    for cid,source,c in zip(CASES,m['sources'],m['cases']):
        if source['case_id']!=cid or c['case_id']!=cid:raise InvalidArtifact('source order differs')
        parent=read_json(verify(source['parent']));old=read_json(verify(source['old_pool']));om=read_json(verify(old['manifest']))
        if om['source_manifest']!=source['parent'] or parent['settings']!=original_completion.SETTINGS or next(x for x in parent['cases'] if x['case_id']==cid)!=c:raise InvalidArtifact('original context changed')
        for k in ('model','software','orca','cpu_python','gpu_python','resources'):
            if m[k]!=parent[k]:raise InvalidArtifact('original electronic/runtime method differs')
        data=[];pair=[]
        for z in ('Ca','La'):
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z));expected.append(t['task_id']);pair.append(t)
            if t!=next(oldt for oldt in parent['tasks'] if oldt['task_id']==t['task_id']):raise InvalidArtifact('source task changed')
            ep=c['origin_row']['native_endpoints'][z];native,forces=original_union.origin(ep,m['model']);kin=Kinematics(read_json(verify(t['mapping']))['context'])
            coords,v,raw,normed,_=project(kin.data,np.zeros(len(kin.modes)),forces);data.append((kin,v,normed))
            point=t['origin_reuse']['point']
            if not np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad']) or native['energy_eV']!=t['q0']['components']['MACE_eV']:raise InvalidArtifact('origin energy/gradient changed')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('q0 mapping differs')
            angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
        if pair[0]['mapping']!=pair[1]['mapping'] and read_json(verify(pair[0]['mapping']))!=read_json(verify(pair[1]['mapping'])):raise InvalidArtifact('paired maps differ')
        choice=preview(data[0][0].modes,data[0][1],data[0][2],data[1][2])
        if choice!=c['selection'] or any(t['selector']!=choice for t in pair):raise InvalidArtifact('same paired selector required')
    if [t['task_id'] for t in m['tasks']]!=expected:raise InvalidArtifact('task order differs')
    return {'status':'validated','manifest':record(manifest),'cases':4,'optimizer_starts':8,'cross_MACE_maximum':8,'GFN2_maximum':32,'new_q0_calls':0,'model_calls_in_validation':0}


def execute(manifest):return engine()[0].execute(manifest)
def collect(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)


def compare(collection,output):
    fresh=read_json(collection);pm=read_json(verify(fresh['manifest']));m=read_json(verify(pm['source_manifest']));validate_pool(verify(fresh['manifest']))
    ref=read_json(verify(m['reference']));bands=ref['variants']['operational']['bands'];rows=[]
    for source in m['sources']:
        cid=source['case_id'];c=next(c for c in fresh['cases'] if c['case_id']==cid);old=next(c for c in read_json(verify(source['old_pool']))['cases'] if c['case_id']==cid)
        endpoints={};components=[]
        for z in ('Ca','La'):
            oldr=read_json(Path(verify(source['parent'])).parent/'proposals'/(cid+'__'+z)/'result.json');npth=Path(verify(pm['source_manifest'])).parent/'proposals'/(cid+'__'+z)/'result.json';nr=read_json(npth) if npth.exists() else None
            ok=bool(nr and nr['status']=='proposal_available' and oldr['status']=='proposal_available');row={'old_receipt':record(Path(verify(source['parent'])).parent/'proposals'/(cid+'__'+z)/'result.json'),'new_receipt':record(npth) if nr else None,'status':'available' if ok else 'unavailable'}
            if ok:
                a=np.array([x[1:] for x in xyz(verify(oldr['proposal']['coordinate']))]);b=np.array([x[1:] for x in xyz(verify(nr['proposal']['coordinate']))]);de=(nr['proposal']['MACE_eV']-oldr['proposal']['MACE_eV'])*EV_TO_KCAL
                row.update(native_proposal_delta_kcal_mol=de,native_energy_pass=abs(de)<=GATES['proposal_energy_kcal_mol'],maximum_coordinate_delta_A=float(np.linalg.norm(b-a,axis=1).max()),full_q_delta=(np.array(nr['proposal']['full_q'])-np.array(oldr['proposal']['full_q'])).tolist(),old_optimizer=oldr['optimizer'],new_optimizer=nr['optimizer'],old_wall_seconds=oldr['wall_seconds'],new_wall_seconds=nr['wall_seconds'],old_geometry=oldr['final_geometry'],new_geometry=nr['final_geometry'],old_boundary=oldr['boundary_flag'],new_boundary=nr['boundary_flag'])
            endpoints[z]=row
            if c['pool']['status']=='available':
                for name in ('adaptive_Ca','adaptive_La'):
                    oc=old['matrix'][z][old['aliases'][name]['representative']]['components'];nc=c['matrix'][z][c['aliases'][name]['representative']]['components']
                    dv=(nc['GFN2_vacuum_hartree']-oc['GFN2_vacuum_hartree'])*HA_TO_KCAL;da=(nc['GFN2_ALPB_hartree']-oc['GFN2_ALPB_hartree'])*HA_TO_KCAL
                    components.append({'metal':z,'candidate':name,'MACE_delta_kcal_mol':(nc['MACE_eV']-oc['MACE_eV'])*EV_TO_KCAL,'vacuum_delta_kcal_mol':dv,'ALPB_delta_kcal_mol':da,'transfer_delta_kcal_mol':da-dv,'vacuum_pass':abs(dv)<=.1,'ALPB_pass':abs(da)<=.1})
        available=c['pool']['status']=='available';newR=c['pool']['operational']['composite_R_model_kcal_mol'] if available else None;oldR=old['pool']['operational']['composite_R_model_kcal_mol'];dr=newR-oldR if available else None
        tail=endpoints['La'];tailpass=bool(tail['status']=='available' and tail['new_optimizer']['function_evaluations']<tail['old_optimizer']['function_evaluations'] and tail['new_wall_seconds']<tail['old_wall_seconds']) if cid==CASES[0] else None
        rows.append({'case_id':cid,'endpoints':endpoints,'components':components,'old_pool':old['pool'],'new_pool':c['pool'],'old_R':oldR,'new_R':newR,'delta_R':dr,'R_pass':abs(dr)<=.2 if dr is not None else False,'old_decision':decision(oldR,bands),'new_decision':decision(newR,bands),'Q89_tail_removed':tailpass,'numerical_pass':bool(available and all(v.get('native_energy_pass',False) for v in endpoints.values()) and len(components)==4 and all(v['vacuum_pass'] and v['ALPB_pass'] for v in components) and abs(dr)<=.2)})
    result={'protocol_id':PROTOCOL,'manifest':pm['source_manifest'],'collection':record(collection),'reference':m['reference'],'bands_used_only_for_transfer_check':bands,'gates':GATES,'rows':rows,'all_numerical_pass':all(x['numerical_pass'] for x in rows),'Q89_tail_removed':rows[0]['Q89_tail_removed'],'calibration_changed':False,'production_changed':False,'new_calls_in_comparison':0,'implementation':record(__file__)}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    specs={'prepare':('canonical','pilot','canonical_pool','pilot_pool','reference','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output'),'compare':('collection','output')}
    for op,args in specs.items():
        a=sub.add_parser(op)
        for k in args:a.add_argument('--'+k.replace('_','-'),required=True,type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
if __name__=='__main__':main()
