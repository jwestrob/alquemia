"""Thin approved triple-membership adapter for existing four-mode/native pool kernels."""
from __future__ import annotations
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from second_shell_context import parent_state
from coordination_preparation_context import geometry
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
import adaptive_angular_proposals as angular
import union_adaptive as original
import slsqp_precision as precision
from slsqp_precision_transfer import qualification,low_prepare
from compact_solvation import completed
from compact_solvation_compare import MINIMUM_CALIBRATION_GAP
from accommodation_folds_compare import decision
from union_triple_preparation import STRESS

PROTOCOL='Nikasha_three_La_union_four_mode_ftol1e8_proposals_v1'
POOL_PROTOCOL='Nikasha_three_La_union_minimal_native_OMOL_GFN2_rank1_v1'
CANDIDATES=['origin','adaptive_Ca','adaptive_La']


def prepare(origins,pool_reuse,template,old_reference,agreement,output):
    data=read_json(origins);inputs=read_json(verify(data['inputs']));reuse=read_json(pool_reuse);base=read_json(template)
    if inputs['pool_reuse']!=record(pool_reuse) or reuse['complete_pool_reuses']!=19 or data['denominator']!=29:
        raise InvalidArtifact('approved29 source/19 reused pools required')
    qualification(verify(inputs['rank_qualification']))
    needed=[r['case_id'] for r in reuse['rows'] if r['status']!='exact_new_ftol_pool_reuse']+[STRESS]
    if len(needed)!=10 or len(set(needed))!=10:raise InvalidArtifact('ten new source pools required')
    rows={r['case_id']:r for r in data['rows']};out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False);cases=[];tasks=[];unavailable=[]
    for cid in needed:
        row=rows[cid]
        if row['status']!='complete':unavailable.append({'case_id':cid,'status':'origin_unavailable','origin_row':row});continue
        prepared=row['source'];state=parent_state(prepared['original_core'],inputs['config']['topology'],require_endpoint_receipts=False)
        audit=read_json(verify(row['source_preparation']));points={};projections={};new=[]
        for z in ('Ca','La'):
            ep=row['native_endpoints'][z];native,forces=original.origin(ep,base['model'])
            atoms=xyz(verify(ep['xyz']));core=xyz(verify(prepared['original_core']['endpoints'][z]['xyz']))
            geo=geometry(state,audit,atoms,core);gp=out/'maps'/(cid+'__'+z+'.json');write_new(gp,geo)
            kin=Kinematics(geo['context']);zero=np.zeros(len(kin.modes));coords,v,raw,normed,lengths=project(kin.data,zero,forces)
            low=row['solvent_endpoints'][z]
            point={'status':'complete','coordinate':ep['xyz'],'MACE':ep['native_MACE_receipt'],'MACE_eV':native['energy_eV'],
                'forces':native['forces'],'full_q':zero.tolist(),'active_q_radian':[0.]*4,
                'source_receipt_format':'actual_legacy_native_OMOL','reused_scientific_origin':True}
            points[z]=point;projections[z]=(v,raw,normed)
            new.append({'task_id':cid+'__'+z,'case_id':cid,'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],
                'mapping':record(gp),'source_preparation':row['source_preparation'],'mode_count':len(kin.modes),'q0_status':'available',
                'q0':{'components':{'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],
                    'GFN2_ALPB_hartree':low['alpb']['energy_hartree']},'low':low,'native_MACE_receipt':ep['native_MACE_receipt']}})
        ca,la=new;paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if read_json(verify(ca['mapping']))!=read_json(verify(la['mapping'])):raise InvalidArtifact('paired physical maps differ')
        kin=Kinematics(read_json(verify(ca['mapping']))['context']);choice=preview(kin.modes,projections['Ca'][0],projections['Ca'][2],projections['La'][2])
        ids=[r['id'] for r in choice['selected']]
        if len(ids)!=4:raise InvalidArtifact('four independent physical modes unavailable')
        selected=[[m['id'] for m in kin.modes].index(x) for x in ids]
        for t in new:
            z=t['metal'];point={**points[z],'gradient_kcal_mol_rad':projections[z][1][selected].tolist()}
            t.update(active_indices=selected,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,selector=choice,
                origin_reuse={'proposal_receipt':points[z]['MACE'],'point':point,'receipt_kind':'actual_archived_native_origin_not_a_proposal'})
            angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))]);tasks.append(t)
        cases.append({'case_id':cid,'status':'prepared','origin_row':row,'union':prepared,'selection':choice,'atom_count':len(atoms)})
    m={k:base[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    m.update(protocol_id=PROTOCOL,settings=precision.SETTINGS,agreement=record(agreement),origin_collection=record(origins),
        source=record(template),pool_reuse=record(pool_reuse),old_reference=record(old_reference),reference=None,rank_qualification=inputs['rank_qualification'],
        cases=cases,tasks=tasks,declared_new_case_ids=needed,declared_case_ids=[c['case_id'] for c in cases],unavailable=unavailable,
        optimizer_software={'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
        population='canonical25_crystal3_separate_stress_three_La_membership_new10',maximum_optimizer_starts=2*len(cases),new_origin_calls=0,
        maximum_cross_MACE_calls=2*len(cases),maximum_GFN2_calls=8*len(cases),new_DFT_calls=0,production_changed=False,
        implementation=original.snapshot(out/'implementation'))
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest);data=read_json(verify(m['origin_collection']));reuse=read_json(verify(m['pool_reuse']));base=read_json(verify(m['source']))
    qualification(verify(m['rank_qualification']));verify(m['agreement']);verify(m['old_reference'])
    needed=[r['case_id'] for r in reuse['rows'] if r['status']!='exact_new_ftol_pool_reuse']+[STRESS]
    if (m['protocol_id']!=PROTOCOL or m['settings']!=precision.SETTINGS or m['declared_new_case_ids']!=needed or len(needed)!=10 or
        m['maximum_optimizer_starts']!=2*len(m['cases']) or len(m['tasks'])!=2*len(m['cases']) or m['reference'] is not None or m['production_changed']):
        raise InvalidArtifact('frozen triple-adaptive scope/settings differ')
    if set(m['declared_case_ids'])|{r['case_id'] for r in m['unavailable']}!=set(needed):raise InvalidArtifact('missing source denominator')
    if [c['case_id'] for c in m['cases']]!=m['declared_case_ids']:raise InvalidArtifact('source order differs')
    for pin in m['implementation'].values():verify(pin)
    for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources'):
        if m[k]!=base[k]:raise InvalidArtifact('runtime/model changed')
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer implementation differs')
    rows={r['case_id']:r for r in data['rows']}
    for c in m['cases']:
        if c['origin_row']!=rows[c['case_id']] or c['union']!=rows[c['case_id']]['source']:raise InvalidArtifact('prepared source changed')
        if c['origin_row']['status']!='complete':raise InvalidArtifact('incomplete original endpoint')
        projections=[];pair=[]
        for z in ('Ca','La'):
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],z));pair.append(t)
            ep=c['origin_row']['native_endpoints'][z];native,forces=original.origin(ep,m['model'])
            if (t['xyz'],t['charge'],t['multiplicity'],t['source_preparation'])!=(ep['xyz'],ep['charge'],ep['multiplicity'],c['origin_row']['source_preparation']):raise InvalidArtifact('source origin state changed')
            kin=Kinematics(read_json(verify(t['mapping']))['context']);coords,v,raw,normed,lengths=project(kin.data,np.zeros(len(kin.modes)),forces);projections.append((kin,v,normed))
            point=t['origin_reuse']['point']
            if (point['MACE']!=ep['native_MACE_receipt'] or point['forces']!=native['forces'] or point['MACE_eV']!=native['energy_eV'] or
                any(point['full_q']) or not np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad'])):raise InvalidArtifact('actual origin force projection differs')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],rtol=0,atol=1e-12):raise InvalidArtifact('map origin differs')
            for medium,key in [('vacuum','GFN2_vacuum_hartree'),('alpb','GFN2_ALPB_hartree')]:
                low=c['origin_row']['solvent_endpoints'][z][medium]
                if t['q0']['low'][medium]!=low or low['energy_hartree']!=t['q0']['components'][key]:raise InvalidArtifact('origin low components differ')
            if t['q0']['components']['MACE_eV']!=native['energy_eV']:raise InvalidArtifact('origin MACE component differs')
            angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
        if read_json(verify(pair[0]['mapping']))!=read_json(verify(pair[1]['mapping'])):raise InvalidArtifact('paired map differs')
        chosen=preview(projections[0][0].modes,projections[0][1],projections[0][2],projections[1][2])
        if chosen!=c['selection'] or any(t['selector']!=chosen for t in pair):raise InvalidArtifact('paired mode selector changed')
        ids=[x['id'] for x in chosen['selected']]
        indices=[[x['id'] for x in projections[0][0].modes].index(mid) for mid in ids]
        if len(ids)!=4 or any(t['active_indices']!=indices or t['active_mode_ids']!=ids for t in pair):
            raise InvalidArtifact('actual four selected coordinates differ')
    return {'status':'validated','manifest':record(manifest),'declared_new_sources':10,'prepared_sources':len(m['cases']),
        'new_searches':len(m['tasks']),'maximum_cross_MACE':2*len(m['cases']),'maximum_GFN2':8*len(m['cases']),'new_origin_calls':0,'new_DFT':0}


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


def reference(collection,output):
    d=read_json(collection);pm=read_json(verify(d['manifest']));m=read_json(verify(pm['source_manifest']));validate_pool(verify(d['manifest']))
    origins=read_json(verify(m['origin_collection']));inputs=read_json(verify(origins['inputs']));reuse=read_json(verify(m['pool_reuse']))
    available={c['case_id']:c for c in d['cases']};saved={r['case_id']:r for r in reuse['rows'] if r['status']=='exact_new_ftol_pool_reuse'}
    old=read_json(verify(m['old_reference']));rows=[]
    for source in inputs['cases']:
        cid=source['case_id'];reused=cid in saved;case=saved[cid]['pool'] if reused else available.get(cid)
        label=source['source']['expected_class'];canonical=bool(source['source'].get('canonical_coordinate_match',False))
        row={'case_id':cid,'root_case_id':source['source'].get('root_case_id',cid),'expected_class':label,
            'role':'canonical_calibration' if canonical else 'separate_stress' if cid==STRESS else 'consumed_crystal_transfer',
            'pool_reused':reused,'source_preparation':source,'status':case['pool']['status'] if case else 'unavailable','pool':case['pool'] if case else None,
            'matrix':case['matrix'] if case else None,'source_collection':saved[cid]['pool_collection'] if reused else record(collection),
            'old_tenfold_band_transfer':{}}
        for v in ('operational','mathematical'):
            value=case['pool'][v]['composite_R_model_kcal_mol'] if case and case['pool']['status']=='available' else None
            row['old_tenfold_band_transfer'][v]={'R':value,'decision':decision(value,old['variants'][v]['bands'])}
        rows.append(row)
    canonical=[r for r in rows if r['role']=='canonical_calibration']
    if len(rows)!=29 or len(canonical)!=25:raise InvalidArtifact('frozen calibration population differs')
    variants={}
    for v in ('operational','mathematical'):
        ready=[r for r in canonical if r['status']=='available'];r={'status':'incomplete','bands':None,'available':len(ready),'denominator':25,'gap_model_kcal_mol':None}
        if len(ready)==25:
            ca=[x['pool'][v]['composite_R_model_kcal_mol'] for x in ready if x['expected_class']=='Ca'];la=[x['pool'][v]['composite_R_model_kcal_mol'] for x in ready if x['expected_class']=='La']
            gap=min(la)-max(ca);r.update(status='available' if gap>MINIMUM_CALIBRATION_GAP else 'overlapping_or_unresolved',gap_model_kcal_mol=gap,
                class_extrema={'Ca_max':max(ca),'La_min':min(la)},class_spread={'Ca':max(ca)-min(ca),'La':max(la)-min(la)})
            if r['status']=='available':r['bands']=r['class_extrema']
        variants[v]=r
    for r in rows:
        r['own_reference_decisions']={v:decision(r['old_tenfold_band_transfer'][v]['R'],variants[v]['bands']) if variants[v]['bands'] else 'unavailable_reference' for v in variants}
    result={'protocol_id':POOL_PROTOCOL,'reference_id':'Nikasha_three_La_membership_ftol1e8_canonical25_v1','frozen_UTC':datetime.now(timezone.utc).isoformat(),
        'canonical_denominator':25,'crystal_denominator':3,'separate_stress_denominator':1,'rows':rows,'variants':variants,
        'collection':record(collection),'pool_reuse':m['pool_reuse'],'old_reference':m['old_reference'],'model':m['model'],'settings':pm['settings'],
        'optimizer_settings':m['settings'],'rank_qualification':m['rank_qualification'],'new_calls_in_reference':0,'production_changed':False,
        'crystals_or_stress_used_for_fit':False,'minimum_gap_model_kcal_mol':MINIMUM_CALIBRATION_GAP,'agreement':m['agreement']}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    commands={'prepare':('origins','pool_reuse','template','old_reference','agreement','output'),'validate':('manifest',),
        'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),'validate_pool':('manifest',),
        'execute_mace':('manifest',),'collect_pool':('manifest','output'),'reference':('collection','output')}
    for op,fields in commands.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f.replace('_','-'),required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))


if __name__=='__main__':main()
