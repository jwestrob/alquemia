"""Validate only the fixed, eligible minima from a completed mechanics assessment."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,read_json,record,verify,write_new,xyz
from mace_curvature import actual_collection
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz
from mace_mechanics_run import preparation,TOL,UNAVAILABLE
from mace_short_engine import COMPONENT

SCHEMA='alquemia.mace_mechanics_minimum_short.v1'
ELIGIBLE='eligible_for_independent_DFT_minimum_validation'


def eligible(assessment):
    a=read_json(assessment)
    if a['status']!='complete' or not a['predicted_minima_are_not_validated_scores']:
        raise InvalidArtifact('complete pre-validation prediction record required')
    for ref in [a['implementation'],*a['sources'].values()]:verify(ref)
    selected={}
    for case,row in a['rows'].items():
        for metal,pred in row['stationary_points'].items():
            if pred['status']!=ELIGIBLE:continue
            z=np.array(pred['normalized_optimum'])
            if (not all(row['gates'].values()) or z.shape!=(2,) or not np.isfinite(z).all() or
                    np.max(np.abs(z))>1.+1e-12 or pred['predicted_relaxation_kcal_mol'] is None or
                    min(pred['fine_eigenvalues_kcal_mol'])<=0 or min(pred['coarse_eigenvalues_kcal_mol'])<=0):
                raise InvalidArtifact('minimum eligibility/state is inconsistent')
            selected[f'{case}_{metal}']={'case_id':case,'metal':metal,'prediction':pred}
    if not 0<len(selected)<=8:raise InvalidArtifact('no eligible minima, or conditional inventory exceeded')
    return a,selected


def prepare(assessment,output):
    # Preparation-only imports: GPU execution environment does not need gemmi.
    from mace_mechanics import physical_model,displaced,membership,H,POLICY
    from ggr_sensitivity import METHOD,ORCA
    a,selected=eligible(assessment);p=preparation(verify(a['sources']['prepared']))
    parent,pm=actual_collection(verify(a['sources']['short']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    qr=out/'quantum';sr=out/'short';qr.mkdir();sr.mkdir()
    pins=snapshot(sr,(*pm['implementation'],'mace_mechanics_minimum.py'))
    tasks=[];quantum=[];minimums={}
    for key,item in selected.items():
        name,metal=item['case_id'],item['metal'];c=read_json(verify(p['cases'][name]))
        repair=read_json(verify(c['source_preparation']))
        graph,coords,metalpos,mk,specs,selection=physical_model(repair)
        if specs!=c['specs'] or selection!=c['selected_coordinate']:raise InvalidArtifact('physical motion definition changed')
        physical=read_json(verify(c['physical_atoms']));z=item['prediction']['normalized_optimum']
        if not np.array_equal(H*np.array(z),item['prediction']['physical_optimum']):raise InvalidArtifact('coordinate units changed')
        core,full,state=displaced(repair,graph,coords,metalpos,mk,specs,z,physical,membership(graph,coords,metalpos,mk))
        directory=out/key;directory.mkdir()
        minimums[key]={**item,'physical_validation':state,'endpoints':{}}
        for kind,atoms in [('core',core),('full',full)]:
            atoms=[(metal if atom[0]=='La' else atom[0],*atom[1:]) for atom in atoms]
            xp=directory/f'{kind}.xyz';write_xyz(xp,atoms)
            center=(c['grids']['center']['endpoints'][metal] if kind=='core' else
                    p['physical_cases'][c['global_id']]['grids']['center']['endpoints'][metal])
            state0={'xyz':record(xp),'charge':center['charge'],'spin_multiplicity':center['spin_multiplicity'],
                    'state':check_atoms(atoms,center['charge'])}
            minimums[key]['endpoints'][kind]=state0
            tasks.append({'task_id':kind+'__'+key,'case_id':name,'metal':metal,'kind':kind,
                          'variant':'minimum_validation','energy_component':COMPONENT,
                          'prediction':record(assessment),'normalized_coordinates':z,**state0})
        td=qr/key;td.mkdir();xp=td/'core.xyz';xp.write_bytes(verify(minimums[key]['endpoints']['core']['xyz']).read_bytes())
        charge=minimums[key]['endpoints']['core']['charge'];ip=td/'endpoint.inp'
        ip.write_text(f'! {METHOD}\n* xyzfile {charge} 1 core.xyz\n')
        settings={'protocol_id':POLICY,'prediction':record(assessment),'source_preparation':c['source_preparation'],
                  'xyz':record(xp),'input':record(ip),'point':'minimum_validation',
                  'physical_coordinates':item['prediction']['physical_optimum'],
                  'charge':charge,'spin_multiplicity':1,'metal':metal,'case_id':name}
        quantum.append({'task_id':key,'case_id':name,'point':'minimum_validation','metal':metal,'charge':charge,'multiplicity':1,
                        'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),'engrad_path':None,
                        'task_type':'single_point','scientific_settings':settings,'cache_key':cache_key(settings)})
    top={'assessment':record(assessment),'original_preparation':a['sources']['prepared'],'minima':minimums,
         'implementation':pins,'new_DFT_calls':len(quantum),'new_short_calls':len(tasks),**UNAVAILABLE}
    write_new(out/'preparation.json',top)
    old_q=read_json(verify(p['DFT_tasks']));qm=copy.deepcopy(old_q)
    qm.update(schema_version='alquemia.mace_mechanics_minimum_DFT.v1',tasks=quantum,implementation=pins,
              minimum_preparation=record(out/'preparation.json'),orca=record(ORCA))
    write_new(qr/'manifest.json',qm)
    m={'schema_version':SCHEMA,'protocol_id':'mace_polar_1m_fixed_minimum_validation_v1',
       'agreement':p['agreement'],'minimum_preparation':record(out/'preparation.json'),
       'source_assessment':record(assessment),'parent_short_collection':a['sources']['short'],
       'numerical_reference':pm['numerical_reference'],'software':pm['software'],'model':copy.deepcopy(pm['model']),
       'implementation':pins,'tasks':tasks,'tolerances':TOL,**UNAVAILABLE}
    for task in tasks:task['cache_key']=cache_key({'task':task,'model':m['model'],'software':m['software'],'implementation':pins})
    write_new(sr/'manifest.json',m)
    return validate(sr/'manifest.json')


def validate(manifest):
    m=read_json(manifest);a,selected=eligible(verify(m['source_assessment']))
    top=read_json(verify(m['minimum_preparation']));parent,pm=actual_collection(verify(m['parent_short_collection']))
    if (m['schema_version']!=SCHEMA or m['model']!=pm['model'] or m['software']!=pm['software'] or
            m['tolerances']!=TOL or numerical_parent_gate(m)['status']!='pass' or
            top['assessment']!=m['source_assessment'] or top['original_preparation']!=a['sources']['prepared'] or
            set(top['minima'])!=set(selected) or len(m['tasks'])!=2*len(selected)):
        raise InvalidArtifact('minimum protocol/inventory differs from fixed prediction')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    expected={kind+'__'+key:(kind,key) for key in selected for kind in ('full','core')}
    if {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('minimum task identity changed')
    for task in m['tasks']:
        kind,key=expected[task['task_id']];item=selected[key];saved=top['minima'][key]
        if saved['prediction']!=item['prediction']:raise InvalidArtifact('frozen optimum was changed')
        for field,value in {'kind':kind,'case_id':item['case_id'],'metal':item['metal'],
                            'variant':'minimum_validation','energy_component':COMPONENT,
                            'prediction':m['source_assessment'],'normalized_coordinates':item['prediction']['normalized_optimum'],
                            **saved['endpoints'][kind]}.items():
            if task.get(field)!=value:raise InvalidArtifact('task differs from eligible fixed physical minimum')
        check_atoms(xyz(verify(task['xyz'])),task['charge'])
        payload={k:v for k,v in task.items() if k!='cache_key'}
        if task['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('minimum scientific cache key changed')
    return {'status':'pass','new_short_calls':len(m['tasks']),'new_DFT_calls':len(selected),'manifest':record(manifest)}


def report(prepared,dft,short,output):
    from mace_mechanics_assess import verified_new_dft,all_core_rows
    top=read_json(prepared);a,selected=eligible(verify(top['assessment']))
    dc,dm=verified_new_dft(dft);sc,sm=actual_collection(short);validate(verify(sc['manifest']))
    if dm['minimum_preparation']!=record(prepared) or sm['minimum_preparation']!=record(prepared):
        raise InvalidArtifact('validation calculations have different frozen predictions')
    mc,mm=actual_collection(verify(a['sources']['core']));p=preparation(verify(a['sources']['prepared']))
    core=all_core_rows(mc,p,verify(mm['reference_mace']));center_short,_=actual_collection(verify(a['sources']['short']))
    rows={};scores={}
    for key,item in selected.items():
        name,metal=item['case_id'],item['metal'];c=a['rows'][name];mapping=read_json(verify(p['cases'][name]))
        d0=c['grid_energies_relative_kcal_mol'][metal]['DFT_center_hartree']
        delta_dft=(dc['rows'][key]['energy_hartree']-d0)*HA_TO_KCAL
        full0=center_short['rows'][f'full__{mapping["global_id"]}_center_{metal}']['energy_eV']
        core0=core[f'{name}_center_{metal}']['energy_components_eV'][COMPONENT]
        delta_j=((sc['rows']['full__'+key]['energy_eV']-full0)-(sc['rows']['core__'+key]['energy_eV']-core0))*EV_TO_KCAL
        linear=float(np.array(c['gradient_components'][metal]['DFT_scaled_b'])@item['prediction']['normalized_optimum'])
        actual=delta_dft+delta_j;predicted=item['prediction']['predicted_relaxation_kcal_mol']
        tolerance=max(.02,.25*abs(delta_dft-linear));error=predicted-actual
        rows[key]={'DFT_change_kcal_mol':delta_dft,'J_change_kcal_mol':delta_j,
                   'predicted_relaxation_kcal_mol':predicted,'actual_energy_change_kcal_mol':actual,
                   'error_kcal_mol':error,'tolerance_kcal_mol':tolerance,'pass':abs(error)<=tolerance,
                   'actual_energy_lower_than_center':actual<0,
                   'actual_stationary_minimum_verified':False,
                   'validated_endpoint_response_kcal_mol':actual if abs(error)<=tolerance else None,
                   'prediction':item['prediction'],'DFT_receipt':dc['rows'][key]['receipt']}
    for name,c in a['rows'].items():
        states={metal:rows.get(f'{name}_{metal}',{}).get('validated_endpoint_response_kcal_mol') for metal in ('La','Ca')}
        delta=None if any(v is None for v in states.values()) else states['Ca']-states['La']
        scores[name]={'endpoint_response_kcal_mol':states,'response_correction_kcal_mol':delta,
                     'response_status':'both_endpoints_validated' if delta is not None else 'response_model_not_validated',
                     'DFT_R_kcal_mol':c['DFT_R_kcal_mol'],
                     'corrected_R_kcal_mol':None if delta is None else c['DFT_R_kcal_mol']+delta,
                     'calibrated_class':None,'entropy_correction_kcal_mol':None}
    # No paired claim can be made unless both GGR representations are supported.
    x,y=(scores[k]['response_correction_kcal_mol'] for k in ('GGR_extended','GGR_connected'))
    partition={'status':'unavailable','difference_kcal_mol':None,'pass':None}
    if x is not None and y is not None:
        partition={'status':'computed','difference_kcal_mol':y-x,'pass':abs(y-x)<=2.}
    result={'status':'complete','sources':{k:record(v) for k,v in [('prepared',prepared),('DFT',dft),('short',short)]},
            'implementation':record(__file__),'rows':rows,'scores':scores,'partition_check':partition,
            'validation_scope':'energy change at a fixed predicted minimum; no DFT gradient at that point',
            'predictive_improvement_claimed':False,'all_cases_consumed_development':True}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for key in ('assessment','output'):a.add_argument('--'+key,required=True)
    a=sub.add_parser('report')
    for key in ('prepared','dft','short','output'):a.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    result={'prepare':prepare,'report':report}[command](**args)
    print(json.dumps(result if command=='prepare' else {'status':result['status'],'endpoint_passes':{k:v['pass'] for k,v in result['rows'].items()},'scores':result['scores']},indent=2))
