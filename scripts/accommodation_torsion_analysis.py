"""Analyze finite physical work without assigning an affinity correction."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics
from compact_solvation import completed
from accommodation_torsion_profiles import CASES,ANGLES,PROTOCOL


def donor_distances(point):
    kin=Kinematics(read_json(verify(point['mapping']))['context']);positions=kin.evaluate(point['q'])[0]
    result={};roles=read_json(verify(point['mapping']))['context']['source_atom_metadata']
    for i,meta in enumerate(roles):
        if meta is None or meta['element'] not in ('O','N'):continue
        ident=f"{meta['chain']}/{meta['resnum']}/{meta['atom']}"
        distance=float(np.linalg.norm(positions[i]-positions[0]))
        if distance<4.5:result[ident]=distance
    return result


def analyze(manifest,mace_result,low_collection,output,dft_collection=None):
    m=read_json(manifest);d=read_json(verify(m['design']));mr=read_json(mace_result);low=read_json(low_collection)
    if mr['manifest']!=record(manifest) or low['manifest']!=m['low_manifest']:raise InvalidArtifact('collection source differs')
    mace={r['task_id']:r for r in mr['rows']};lr={r['task_id']:r for r in low['rows']}
    native={}
    if dft_collection:
        dc=read_json(dft_collection)
        if dc['manifest']!=m['dft_manifest']:raise InvalidArtifact('DFT source differs')
        native={r['point_id']:r for r in dc['rows']}
    native.update({k:{**v,'status':'complete','reused':True} for k,v in d['reused_DFT_centers'].items()})
    points=[]
    for p in d['points']:
        tid=p['task_id'];a=mace[tid];v,b=[lr[tid+'__'+s] for s in ('vacuum','alpb')]
        row={k:p[k] for k in ('task_id','case_id','metal','point','role','mode_id','angle_radian','geometry_checks')}
        row.update(status='unavailable',MACE_eV=a['energy_eV'],GFN2_vacuum_hartree=v['energy_hartree'],GFN2_ALPB_hartree=b['energy_hartree'],
            composite_energy_kcal_mol=None,solvent_transfer_kcal_mol=None,DFT_hartree=None,donor_distances_A=donor_distances(p))
        if all(x['status']=='complete' for x in (a,v,b)):
            solvent=(b['energy_hartree']-v['energy_hartree'])*HA_TO_KCAL
            row.update(status='complete',solvent_transfer_kcal_mol=solvent,composite_energy_kcal_mol=a['energy_eV']*EV_TO_KCAL+solvent)
        n=native.get(tid)
        if n and n['status']=='complete':row['DFT_hartree']=n['energy_hartree']
        points.append(row)
    pi={p['task_id']:p for p in points};profiles=[];validation=[];controls=[]
    for case in d['cases']:
        cid=case['case_id']
        origins={z:pi[cid+'__origin__'+z] for z in ('Ca','La')}
        for role,mode in case['modes'].items():
            rows=[]
            for angle in ANGLES:
                pair={}
                for metal in ('Ca','La'):
                    p=next(p for p in points if p['case_id']==cid and p['metal']==metal and
                        ((angle==0 and p['point']=='origin') or (p['role']==role and p['angle_radian']==angle)))
                    origin=origins[metal];work={}
                    if p['status']==origin['status']=='complete':
                        work={'MACE':(p['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL,
                            'solvent_transfer':p['solvent_transfer_kcal_mol']-origin['solvent_transfer_kcal_mol']}
                        work['composite']=work['MACE']+work['solvent_transfer']
                    native_work=(p['DFT_hartree']-origin['DFT_hartree'])*HA_TO_KCAL if p['DFT_hartree'] is not None and origin['DFT_hartree'] is not None else None
                    pair[metal]={'task_id':p['task_id'],'work_kcal_mol':work,'DFT_work_kcal_mol':native_work,
                        'donor_distances_A':p['donor_distances_A']}
                    if native_work is not None and angle:
                        validation.append({'case_id':cid,'role':role,'angle_radian':angle,'metal':metal,
                            'DFT_work_kcal_mol':native_work,'composite_work_kcal_mol':work.get('composite'),
                            'composite_minus_DFT_kcal_mol':work['composite']-native_work if work else None,
                            'same_work_sign':bool(work['composite']*native_work>0) if work else None})
                diff={term:pair['Ca']['work_kcal_mol'][term]-pair['La']['work_kcal_mol'][term] for term in ('MACE','solvent_transfer','composite')} if all(pair[z]['work_kcal_mol'] for z in pair) else None
                ndiff=pair['Ca']['DFT_work_kcal_mol']-pair['La']['DFT_work_kcal_mol'] if all(pair[z]['DFT_work_kcal_mol'] is not None for z in pair) else None
                rows.append({'angle_radian':angle,'endpoints':pair,'Ca_minus_La_work_kcal_mol':diff,'DFT_Ca_minus_La_work_kcal_mol':ndiff})
            derivatives={}
            for metal in ('Ca','La'):
                good=all(r['endpoints'][metal]['work_kcal_mol'] for r in rows)
                if good:
                    vals={r['angle_radian']:r['endpoints'][metal]['work_kcal_mol']['composite'] for r in rows}
                    derivatives[metal]={'central_slope_kcal_mol_rad':{str(h):(vals[h]-vals[-h])/(2*h) for h in (.2,.4)},
                        'central_curvature_kcal_mol_rad2':{str(h):(vals[h]+vals[-h])/(h*h) for h in (.2,.4)},
                        'lowest_sampled_angle_radian':min(vals,key=vals.get),'boundary_is_lowest_sample':min(vals,key=vals.get) in (-.4,.4)}
            profiles.append({'case_id':cid,'label_scope':case['label_scope'],'role':role,'mode_id':mode,'rows':rows,'finite_difference_diagnostics':derivatives})
        if cid in CASES[:2]:
            cm=read_json(verify(d['compact_source']));eps={}
            for metal in ('Ca','La'):
                a=case['archived_control_endpoints'][metal];de=(origins[metal]['MACE_eV']-a['native_MACE_energy_eV'])*EV_TO_KCAL
                old=[]
                for medium in ('vacuum','alpb'):
                    t=next(t for t in cm['all_tasks'] if (t['case_id'],t['representation'],t['metal'],t['medium'])==(cid,'context',metal,medium))
                    pin=cm['reused'].get(t['task_id']) or completed(verify(d['compact_source']),t['task_id'])
                    if pin is None:raise InvalidArtifact('archived control native source missing')
                    old.append(pin['energy_hartree'])
                shift=origins[metal]['solvent_transfer_kcal_mol']-(old[1]-old[0])*HA_TO_KCAL
                eps[metal]={'MACE_delta_kcal_mol':de,'transfer_delta_kcal_mol':shift,'composite_delta_kcal_mol':de+shift}
            pair=eps['Ca']['composite_delta_kcal_mol']-eps['La']['composite_delta_kcal_mol']
            controls.append({'case_id':cid,'endpoints':eps,'composite_contrast_delta_kcal_mol':pair,
                'pass':all(abs(v['MACE_delta_kcal_mol'])<=.01 and abs(v['transfer_delta_kcal_mol'])<=.1 for v in eps.values()) and abs(pair)<=.2})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'MACE_result':record(mace_result),'low_collection':record(low_collection),
        'DFT_collection':record(dft_collection) if dft_collection else None,'analyzer':record(__file__),'points':points,'profiles':profiles,
        'complete_points':sum(p['status']=='complete' for p in points),'point_denominator':64,'control_repeats':controls,'DFT_validation':validation,
        'DFT_validation_denominator':12,'affinity_correction':None,'entropy':None,'PLM_accuracy':None,'baseline_changed':False}
    write_new(output,result)
    return {'complete_points':result['complete_points'],'profiles':len(profiles),'DFT_validation':len(validation),'control_repeats':controls}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('manifest','mace_result','low_collection','output'):p.add_argument('--'+k.replace('_','-'),required=True)
    p.add_argument('--dft-collection');print(json.dumps(analyze(**vars(p.parse_args())),indent=2))
if __name__=='__main__':main()
