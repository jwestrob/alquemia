"""Supplemental scoring of one archived fixed-state hydrogen repair."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
import union_adaptive as original
import adaptive_completion
from preparation_coverage import CASE,PROTOCOL as PREPARATION_PROTOCOL

STATIC='archived_H_repair_local_union_native_OMOL_GFN2_static_v1'
PROTOCOL='archived_H_repair_union_four_angular_native_OMOL_SLSQP_v1'
POOL='archived_H_repair_union_adaptive_minimal_native_OMOL_GFN2_v1'
PREPARATION_SHA='5b3f68031e97c7d669b2389636dd5131b2172ee4c6dc13684359e08922cdcade'


def private(name):
    key='_coverage_private_'+name
    spec=importlib.util.spec_from_file_location(key,Path(__file__).with_name(name+'.py'))
    mod=importlib.util.module_from_spec(spec);sys.modules[key]=mod;spec.loader.exec_module(mod);return mod


def preparation(path):
    p=read_json(path)
    if record(path)['sha256']!=PREPARATION_SHA or p['protocol_id']!=PREPARATION_PROTOCOL or p['status']!='prepared':raise InvalidArtifact('exact repaired preparation required')
    for kind in ('local','union'):
        row=p[kind];rep=row['representations']['context'];ep=rep['endpoints']
        if row['case_id']!=CASE or row['status']!='prepared':raise InvalidArtifact('supplemental source differs')
        paired(verify(ep['La']['xyz']),verify(ep['Ca']['xyz']),ep['La']['charge'],ep['Ca']['charge'])
        if [(ep[z]['charge'],ep[z]['multiplicity']) for z in ('Ca','La')]!=[(-3,1),(-2,1)]:raise InvalidArtifact('chemical state differs')
        if len(xyz(verify(ep['Ca']['xyz'])))!={'local':160,'union':184}[kind]:raise InvalidArtifact('context membership differs')
        verify(rep['preparation'])
    return p


def common_checks(m):
    for name in ('agreement','source','software','orca'):verify(m[name])
    for pin in m['implementation'].values():verify(pin)
    source=read_json(verify(m['source']))
    for key in ('model','software','orca','cpu_python','gpu_python','resources'):
        if m[key]!=source[key]:raise InvalidArtifact('original electronic/runtime policy changed')


def prepare_origins(prepared,source,agreement,reference,output):
    p=preparation(prepared);old=read_json(source)
    if old['protocol_id']!=original.PROTOCOL or old['settings']!=adaptive_completion.SETTINGS:raise InvalidArtifact('original union optimizer required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);cases=[];tasks=[]
    for kind in ('local','union'):
        rep=p[kind]['representations']['context'];case={'case_id':kind,'status':'prepared','source':p[kind],'aliases':{},
          'candidates':[{'id':'origin','xyz':rep['endpoints']['Ca']['xyz']}],'matrix':{z:{} for z in ('Ca','La')}}
        for z in ('Ca','La'):
            ep=rep['endpoints'][z];tid=kind+'__'+z+'__origin'
            task={'task_id':tid,'case_id':kind,'metal':z,'candidate':'origin',**{k:ep[k] for k in ('xyz','charge','multiplicity')},'source_preparation':rep['preparation']}
            tasks.append(task);case['matrix'][z]['origin']={'status':'pending','task_id':tid,'xyz':ep['xyz'],'reused':False}
        cases.append(case)
    m={k:old[k] for k in ('model','software','orca','cpu_python','gpu_python','resources')}
    m.update(protocol_id=STATIC,agreement=record(agreement),source=record(source),preparation=record(prepared),reference=record(reference),
      cases=cases,tasks=tasks,declared_case_ids=['local','union'],settings={**original.pool.SETTINGS,'candidate_order':['origin']},
      GFN2_maxiter=500,numerical_policy_id='native_GFN2_MaxIter500_unchanged_convergence_v1',shard_count=1,new_MACE_cells=4,
      maximum_new_GFN2_calls=8,new_DFT_calls=0,new_optimizations=0,production_changed=False,implementation=original.snapshot(out/'implementation'))
    mp=out/'manifest.json';write_new(mp,m);check=validate_origins(mp);original.pool.low_prepare(mp);write_new(out/'PREFLIGHT.json',check);return check


def validate_origins(manifest):
    m=read_json(manifest);p=preparation(verify(m['preparation']));common_checks(m);verify(m['reference'])
    if (m['protocol_id']!=STATIC or m['GFN2_maxiter']!=500 or m['declared_case_ids']!=['local','union'] or
        len(m['tasks'])!=4 or m['new_MACE_cells']!=4 or m['maximum_new_GFN2_calls']!=8):raise InvalidArtifact('finite origin scope differs')
    expected=[]
    for kind,c in zip(('local','union'),m['cases']):
        if c['case_id']!=kind or c['source']!=p[kind]:raise InvalidArtifact('origin source differs')
        rep=p[kind]['representations']['context']
        for z in ('Ca','La'):
            tid=kind+'__'+z+'__origin';expected.append(tid);t=next(t for t in m['tasks'] if t['task_id']==tid);ep=rep['endpoints'][z]
            if any(t[k]!=ep[k] for k in ('xyz','charge','multiplicity')) or t['source_preparation']!=rep['preparation']:raise InvalidArtifact('origin coordinates/state changed')
            if c['matrix'][z]['origin']!={'status':'pending','task_id':tid,'xyz':ep['xyz'],'reused':False}:raise InvalidArtifact('origin result invented')
    if [t['task_id'] for t in m['tasks']]!=expected:raise InvalidArtifact('origin order differs')
    return {'status':'validated','manifest':record(manifest),'MACE_calls':4,'GFN2_calls':8,'new_calls_in_validation':0}


def execute_origins(manifest):
    pool=private('nikasha_pool');pool.validate=validate_origins;return pool.execute_mace(manifest)


def collect_origins(manifest,output):
    validate_origins(manifest);return original.pool.collect(manifest,output)


def warm_origin(ep,model):
    r=read_json(verify(ep['native_MACE_receipt']));req=read_json(verify(r['request']));m=read_json(verify(req['manifest']))
    if m['protocol_id']!=STATIC or r['status']!='complete' or r['model']!=model or r['energy_eV']!=ep['native_MACE_energy_eV']:raise InvalidArtifact('actual warm origin missing or incompatible')
    t=next(t for t in m['tasks'] if t['task_id']==req['task_id'])
    if any(req[k]!=ep[k] or req[k]!=t[k] for k in ('xyz','charge','multiplicity')):raise InvalidArtifact('origin request state differs')
    f=np.load(verify(r['forces']),allow_pickle=False)
    if f.shape!=(len(xyz(verify(req['xyz']))),3) or not np.isfinite(f).all():raise InvalidArtifact('actual origin forces missing')
    return r,f


def population(inputs,calibration):
    d=read_json(inputs)
    if (d['population']!='supplemental_archived_H_repair_one' or d['cases']!=[{'case_id':CASE}] or d['calibration']!=record(calibration)):raise InvalidArtifact('supplemental membership differs')
    preparation(verify(d['preparation']));return d['population'],1


def lookup(prepared,crystals,calibration,transfer,inputs):
    p=preparation(prepared);d=read_json(inputs);population(inputs,calibration)
    if record(prepared)!=d['preparation'] or record(crystals)!=record(prepared) or record(transfer)!=record(calibration):raise InvalidArtifact('source adapter pins differ')
    col=read_json(calibration);mp=verify(col['manifest']);validate_origins(mp)
    if read_json(mp)['preparation']!=record(prepared):raise InvalidArtifact('scored preparation differs')
    c=next(c for c in col['cases'] if c['case_id']=='union')
    if c['pool']['status']!='available':raise InvalidArtifact('actual union origin incomplete; search unavailable')
    rep=p['union']['representations']['context'];eps={};lows={}
    for z in ('Ca','La'):
        cell=c['matrix'][z]['origin'];eps[z]={**rep['endpoints'][z],'native_MACE_receipt':cell['MACE'],'native_MACE_energy_eV':cell['components']['MACE_eV']}
        lows[z]={s:{k:cell['low'][s][k] for k in ('energy_hartree','output','receipt','manifest','task_id')} for s in ('vacuum','alpb')}
    row={'case_id':CASE,'status':'complete','source_preparation':rep['preparation'],'native_endpoints':eps,'solvent_endpoints':lows}
    return read_json(verify(p['existing_union_preparation'])),[{'case_id':CASE,'actual_union_case_id':CASE,'base':{'case_id':CASE},
        'union':p['union'],'collection':record(calibration),'origin_row':row}]


def engine():
    union=private('union_adaptive');pool=private('nikasha_pool')
    union.PROTOCOL=PROTOCOL;union.POOL_PROTOCOL=POOL;union.lookup=lookup;union.declared_population=population;union.origin=warm_origin
    union.pool=pool;pool.validate=union.validate_pool
    return union,pool


def prepare_proposals(origins,output):
    from second_shell_context import parent_state
    from coordination_preparation_context import geometry
    from mace_site_kinematics import Kinematics
    from adaptive_force_diagnostic import project,preview
    col=read_json(origins);m=read_json(verify(col['manifest']));validate_origins(verify(col['manifest']));out=Path(output).resolve()
    ip=out.parent/'SUPPLEMENTAL_INPUTS.json';write_new(ip,{'population':'supplemental_archived_H_repair_one','cases':[{'case_id':CASE}],
       'preparation':m['preparation'],'calibration':record(origins),'reference':m['reference']})
    prep,rows=lookup(verify(m['preparation']),verify(m['preparation']),origins,origins,ip)
    sm=read_json(verify(m['source']));out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[]
    for row in rows:
        u=row['union'];old=row['origin_row'];context=read_json(verify(old['source_preparation']))
        state=parent_state(u['original_core'],prep['config']['topology'],require_endpoint_receipts=False);points={};projections={};new=[]
        for z in ('Ca','La'):
            ep=old['native_endpoints'][z];native,forces=warm_origin(ep,sm['model']);atoms=xyz(verify(ep['xyz']))
            geo=geometry(state,context,atoms,xyz(verify(u['original_core']['endpoints'][z]['xyz'])));gp=out/'maps'/(CASE+'__'+z+'.json');write_new(gp,geo)
            kin=Kinematics(geo['context']);zeros=np.zeros(len(kin.modes));_,vectors,raw,normalized,_=project(kin.data,zeros,forces);low=old['solvent_endpoints'][z]
            points[z]={'status':'complete','coordinate':ep['xyz'],'MACE':ep['native_MACE_receipt'],'MACE_eV':native['energy_eV'],'forces':native['forces'],
                'full_q':zeros.tolist(),'active_q_radian':[0.]*4,'source_receipt_format':'actual_warm_native_OMOL_request','reused_scientific_origin':True}
            projections[z]=(vectors,raw,normalized)
            new.append({'task_id':CASE+'__'+z,'case_id':CASE,'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],
                'mapping':record(gp),'source_preparation':old['source_preparation'],'mode_count':len(kin.modes),'q0_status':'available',
                'q0':{'components':{'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':low['alpb']['energy_hartree']},
                      'low':low,'native_MACE_receipt':ep['native_MACE_receipt']}})
        if read_json(verify(new[0]['mapping']))!=read_json(verify(new[1]['mapping'])):raise InvalidArtifact('paired physical maps differ')
        choice=preview(kin.modes,projections['Ca'][0],projections['Ca'][2],projections['La'][2]);ids=[v['id'] for v in choice['selected']]
        if len(ids)!=4:raise InvalidArtifact('four common modes unavailable')
        selected=[[v['id'] for v in kin.modes].index(i) for i in ids]
        for t in new:
            z=t['metal'];point={**points[z],'gradient_kcal_mol_rad':projections[z][1][selected].tolist()}
            t.update(active_indices=selected,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,selector=choice,
                origin_reuse={'proposal_receipt':point['MACE'],'point':point,'receipt_kind':'actual_archived_native_origin_not_a_proposal'})
            tasks.append(t)
        cases.append({**row,'status':'prepared','selection':choice,'atom_count':len(atoms)})
    manifest={k:sm[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    manifest.update(protocol_id=PROTOCOL,settings=adaptive_completion.SETTINGS,agreement=m['agreement'],inputs=record(ip),source=m['source'],
      preparation=m['preparation'],crystals=m['preparation'],calibration=record(origins),transfer=record(origins),cases=cases,tasks=tasks,
      declared_case_ids=[CASE],implementation=original.snapshot(out/'implementation'),population='supplemental_archived_H_repair_one',
      optimizer_software={'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),'kernel':record(scipy_kernel.__file__)},
      maximum_optimizer_starts=2,new_origin_calls=0,maximum_cross_MACE_calls=2,maximum_GFN2_calls=8,production_changed=False,new_DFT_calls=0,reference=m['reference'])
    mp=out/'manifest.json';write_new(mp,manifest);result=validate(mp);write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest):return engine()[0].validate(manifest)
def execute(manifest):return engine()[0].execute(manifest)
def collect(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)


def compare(origins,pool,local_reference,union_reference,output):
    from accommodation_folds_compare import decision
    from accommodation_fold_proposals import outcome
    from accommodation_nonlinear import relative_components
    o=read_json(origins);pm=verify(o['manifest']);validate_origins(pm);m=read_json(pm)
    c=read_json(pool);validate_pool(verify(c['manifest']));cm=read_json(verify(c['manifest']));am=read_json(verify(cm['source_manifest']))
    if am['calibration']!=record(origins) or am['preparation']!=m['preparation'] or c['denominator']!=1:raise InvalidArtifact('adaptive source differs')
    lr=read_json(local_reference);ur=read_json(union_reference);ar=read_json(verify(m['reference']))
    if ar['optimizer_settings']!=am['settings'] or ar['settings']!=cm['settings'] or ar['model']!=am['model']:raise InvalidArtifact('adaptive reference/method mismatch')
    if record(local_reference)['sha256']!='451f959d995d289cc246df5e937ef9d33be668ad227d23cdd5fd5bfa213b0d67':raise InvalidArtifact('exact released band authority required')
    lb=lr['frozen_bands']
    if lb['representation']!='context' or lb['solver']!='native':raise InvalidArtifact('released model representation differs')
    bands={'local_static':{'Ca_max':lb['Ca_supported_max_R_model_kcal_mol'],'La_min':lb['La_supported_min_R_model_kcal_mol']},'union_static':ur['bands'],
       'union_adaptive':ar['variants']['operational']['bands'],'union_adaptive_mathematical':ar['variants']['mathematical']['bands']}
    rows={}
    for kind in ('local','union'):
        row=next(x for x in o['cases'] if x['case_id']==kind);chosen=original.pool.choose_rows(row['matrix'],['origin'])
        if chosen!=row['pool']:raise InvalidArtifact('origin score algebra differs')
        name=kind+'_static';value=chosen['operational']['composite_R_model_kcal_mol'] if chosen['status']=='available' else None
        call=decision(value,bands[name]);rows[name]={'status':chosen['status'],'R_model_kcal_mol':value,'decision':call,'outcome':outcome(call,'La'),
          'components':chosen['operational'],'matrix':row['matrix'],'preparation_transfer_check':True}
    ac=c['cases'][0];chosen=original.pool.choose_rows(ac['matrix'],[x['id'] for x in ac['candidates']]) if ac['status']=='prepared' else ac['pool']
    if chosen!=ac['pool']:raise InvalidArtifact('adaptive pool algebra differs')
    for variant,name in (('operational','union_adaptive'),('mathematical','union_adaptive_mathematical')):
        value=chosen[variant]['composite_R_model_kcal_mol'] if chosen['status']=='available' else None;call=decision(value,bands[name])
        work={z:relative_components(ac['matrix'][z][chosen['rows'][z][variant+'_candidate']]['components'],ac['matrix'][z]['origin']['components']) for z in ('Ca','La')} if value is not None else None
        if work and abs((value-rows['union_static']['R_model_kcal_mol'])-(work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol']))>1e-6:raise InvalidArtifact('work/sign/unit closure failed')
        rows[name]={'status':chosen['status'],'R_model_kcal_mol':value,'decision':call,'outcome':outcome(call,'La'),'components':chosen[variant],
          'selected_work':work,'selected_candidates':{z:chosen['rows'][z][variant+'_candidate'] for z in ('Ca','La')} if work else None,'preparation_transfer_check':True}
    proposals=read_json(verify(cm['source']));optimizers={}
    for ep in proposals['endpoints']:
        r=read_json(verify(ep['proposal_receipt'])) if ep['proposal_receipt'] else None
        optimizers[ep['metal']]={'status':ep['status'],'reason':ep['reason'],'boundary_flag':ep['boundary_flag'],'receipt':ep['proposal_receipt'],
          'optimizer':r.get('optimizer') if r else None,'wall_seconds':r.get('wall_seconds') if r else None,'final_geometry':r.get('final_geometry') if r else None}
    result={'case_id':CASE,'expected_class':'La','label_scope':'consumed_PQQ_functional_class_structural_replicate',
       'origins':record(origins),'pool':record(pool),'preparation':m['preparation'],'references':{'local':record(local_reference),'union':record(union_reference),'adaptive':m['reference']},
       'bands':bands,'rows':rows,'optimizers':optimizers,'historical225_modified':False,'historical_preparation_failure_retained':True,
       'independent_biological_observations_added':0,'new_calls_in_comparison':0,'production_changed':False,'implementation':record(__file__)}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,args in {'prepare_origins':('prepared','source','agreement','reference','output'),'validate_origins':('manifest',),
       'execute_origins':('manifest',),'collect_origins':('manifest','output'),'prepare_proposals':('origins','output'),
       'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),
       'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output'),
       'compare':('origins','pool','local_reference','union_reference','output')}.items():
        q=sub.add_parser(op)
        for arg in args:q.add_argument('--'+arg.replace('_','-'),type=Path,required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
if __name__=='__main__':main()
