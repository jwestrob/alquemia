"""Independent native DFT energy/gradient checks at fixed 3D-response predictions."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import (InvalidArtifact,HA_TO_KCAL,BOHR_TO_A,cache_key,read_json,record,verify,write_new,xyz)
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_metal_response import POLICY,TOL,LIMIT,UNAVAILABLE,shift,prediction,collect as collect_grid
SCHEMA='alquemia.mace_metal_minimum_short.v1'


def same_prediction(a,b):
    """Allow float64 solver roundoff across pinned NumPy/BLAS environments."""
    if isinstance(a,dict) and isinstance(b,dict):
        return set(a)==set(b) and all(same_prediction(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list):
        return len(a)==len(b) and all(same_prediction(x,y) for x,y in zip(a,b))
    if isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)) and not isinstance(b,bool):
        return bool(np.isfinite(a) and np.isfinite(b) and np.isclose(a,b,atol=1e-10,rtol=1e-10))
    return a==b


def eligible(assessment):
    a=read_json(assessment);p=read_json(verify(a['preparation']))
    if a['status']!='complete' or a['policy_id']!=POLICY or not a['predictions_are_not_validated_scores']:
        raise InvalidArtifact('complete fixed response predictions required')
    selected={}
    for name,row in a['rows'].items():
        for metal,state in row.items():
            pred=state['prediction']
            g=np.array(state['DFT_gradient_kcal_mol_A'])+state['J_gradient_kcal_mol_A']
            if not same_prediction(pred,prediction(g,state['relative_grid_energies_kcal_mol'])):
                raise InvalidArtifact('stored prediction differs from the fixed grid/gradient model')
            if pred['status']!='eligible_for_DFT_validation':continue
            u=np.array(pred['proposed_displacement_A'])
            if u.shape!=(3,) or not np.isfinite(u).all() or np.linalg.norm(u)>LIMIT or pred['predicted_response_kcal_mol'] is None:
                raise InvalidArtifact('invalid fixed prediction')
            selected[f'{name}_{metal}']={'case_id':name,'metal':metal,'prediction':pred}
    if len(selected)>8:raise InvalidArtifact('conditional inventory exceeded')
    return a,p,selected


def prepare(assessment,output):
    from mace_mechanics import physical_model
    from ggr_sensitivity import membership,METHOD,ORCA
    a,p,selected=eligible(assessment)
    old=read_json(verify(p['assessment']));old_p=read_json(verify(old['sources']['prepared']))
    parent=read_json(verify(old['sources']['short']));pm=read_json(verify(parent['manifest']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    qr=out/'quantum';sr=out/'short';qr.mkdir();sr.mkdir()
    pins=snapshot(sr,(*pm['implementation'],'mace_metal_response.py','mace_metal_minimum.py'))
    tasks=[];quantum=[];states={};excluded={}
    for key,item in selected.items():
        name,metal=item['case_id'],item['metal'];u=np.array(item['prediction']['proposed_displacement_A'])
        repair=read_json(verify(p['cases'][name]['source_preparation']))
        graph,coords,pos,mk,_,_=physical_model(repair)
        if membership(graph,coords,pos+u,mk)!=membership(graph,coords,pos,mk):
            excluded[key]={'status':'donor_membership_changed','prediction':item['prediction']};continue
        d=out/key;d.mkdir();states[key]={**item,'endpoints':{},'physical_validation':'only_selected_metal_moves; typed_donors_unchanged'}
        for kind in ('core','full'):
            src=p['centers'][name][metal][kind+'_state'];base=xyz(verify(src['xyz']));rows,i=shift(base,metal,u)
            xp=d/(kind+'.xyz');write_xyz(xp,rows)
            state={'xyz':record(xp),'charge':src['charge'],'spin_multiplicity':1,'state':check_atoms(rows,src['charge']),
                   'selected_metal_index':i}
            states[key]['endpoints'][kind]=state
            tasks.append({'task_id':kind+'__'+key,'case_id':name,'metal':metal,'kind':kind,'point':'fixed_minimum',
                          'variant':'minimum_validation','energy_component':'interaction_energy',
                          'prediction':record(assessment),'displacement_A':u.tolist(),**state})
        td=qr/key;td.mkdir();xp=td/'core.xyz';xp.write_bytes(verify(states[key]['endpoints']['core']['xyz']).read_bytes())
        charge=states[key]['endpoints']['core']['charge'];ip=td/'endpoint.inp'
        ip.write_text(f'! {METHOD} EnGrad\n* xyzfile {charge} 1 core.xyz\n')
        settings={'protocol_id':POLICY,'prediction':record(assessment),'source_preparation':p['cases'][name]['source_preparation'],
                  'xyz':record(xp),'input':record(ip),'point':'fixed_minimum','physical_coordinates_A':u.tolist(),
                  'charge':charge,'spin_multiplicity':1,'metal':metal,'case_id':name}
        quantum.append({'task_id':key,'case_id':name,'point':'fixed_minimum','metal':metal,'charge':charge,'multiplicity':1,
                        'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),
                        'engrad_path':str(td/'endpoint.engrad'),'task_type':'analytic_gradient',
                        'scientific_settings':settings,'cache_key':cache_key(settings)})
    top={'assessment':record(assessment),'original_preparation':a['preparation'],'states':states,'excluded':excluded,
         'implementation':pins,'new_DFT_calls':len(quantum),'new_short_calls':len(tasks),**UNAVAILABLE}
    write_new(out/'preparation.json',top)
    if not quantum:return {'status':'no_eligible_predictions','prepared':record(out/'preparation.json'),'excluded':excluded}
    qm=copy.deepcopy(read_json(verify(old_p['DFT_tasks'])))
    qm.update(schema_version='alquemia.mace_metal_minimum_DFT.v1',protocol_id=POLICY,tasks=quantum,implementation=pins,
              minimum_preparation=record(out/'preparation.json'),agreement=p['agreement'],orca=record(ORCA))
    write_new(qr/'manifest.json',qm)
    model=copy.deepcopy(pm['model']);model['preparation_policy']=POLICY
    m={'schema_version':SCHEMA,'protocol_id':POLICY+'_minimum_short','kind':'minimum_short',
       'agreement':p['agreement'],'minimum_preparation':record(out/'preparation.json'),'source_assessment':record(assessment),
       'parent_short_collection':old['sources']['short'],'numerical_reference':pm['numerical_reference'],
       'software':pm['software'],'model':model,'implementation':pins,'tasks':tasks,'tolerances':TOL,**UNAVAILABLE}
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':model,'software':m['software'],'implementation':pins})
    write_new(sr/'manifest.json',m);result=validate(sr/'manifest.json')
    from affordable_workflow import dry_run
    result.update(DFT=dry_run(qr/'manifest.json'),excluded=excluded)
    return result


def validate(manifest):
    m=read_json(manifest);a,p,selected=eligible(verify(m['source_assessment']))
    top=read_json(verify(m['minimum_preparation']));pc=read_json(verify(m['parent_short_collection']));pm=read_json(verify(pc['manifest']))
    model=copy.deepcopy(pm['model']);model['preparation_policy']=POLICY
    if (m['schema_version']!=SCHEMA or m['model']!=model or m['software']!=pm['software'] or m['tolerances']!=TOL or
            numerical_parent_gate(m)['status']!='pass' or top['assessment']!=m['source_assessment'] or
            set(top['states'])|set(top['excluded'])!=set(selected) or len(m['tasks'])!=2*len(top['states'])):
        raise InvalidArtifact('fixed-minimum preparation/method/inventory changed')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    expected={kind+'__'+key:(kind,key) for key in top['states'] for kind in ('core','full')}
    if {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('minimum tasks changed')
    for t in m['tasks']:
        kind,key=expected[t['task_id']];item=selected[key];name,metal=item['case_id'],item['metal']
        saved=top['states'][key];u=item['prediction']['proposed_displacement_A']
        if saved['prediction']!=item['prediction'] or t['displacement_A']!=u or any(t.get(k)!=v for k,v in saved['endpoints'][kind].items()):
            raise InvalidArtifact('frozen prediction or prepared state changed')
        expected_fields={'case_id':name,'metal':metal,'kind':kind,'point':'fixed_minimum',
                         'variant':'minimum_validation','energy_component':'interaction_energy',
                         'prediction':m['source_assessment']}
        if any(t.get(k)!=v for k,v in expected_fields.items()):
            raise InvalidArtifact('minimum task identity or energy scope changed')
        rows=xyz(verify(t['xyz']));base=xyz(verify(p['centers'][name][metal][kind+'_state']['xyz']));wanted,mi=shift(base,metal,u)
        if [r[0] for r in rows]!=[r[0] for r in wanted] or not np.allclose([r[1:] for r in rows],[r[1:] for r in wanted],atol=5e-10,rtol=0) or t['selected_metal_index']!=mi:
            raise InvalidArtifact('minimum physical displacement differs')
        check_atoms(rows,t['charge']);payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('minimum cache key changed')
    return {'status':'pass','new_short_calls':len(m['tasks']),'new_DFT_calls':len(top['states']),'manifest':record(manifest)}


def report(prepared,dft,short,output,allow_partial=False):
    from affordable_response import read_engrad
    from mace_mechanics_assess import verified_new_dft
    top=read_json(prepared);a,p,selected=eligible(verify(top['assessment']))
    if allow_partial:
        from ggr_sensitivity import executed
        dc=read_json(dft);dm,actual=executed(verify(dc['manifest']))
        for key,row in actual.items():
            if row['status']!='complete':continue
            if any(dc['rows'][key].get(k)!=v for k,v in row.items()):
                raise InvalidArtifact('partial DFT result differs from actual execution')
            for ref in dc['gradients'][key]['artifacts'].values():verify(ref)
    else:dc,dm=verified_new_dft(dft)
    sc=read_json(short);sm=read_json(verify(sc['manifest']));validate(verify(sc['manifest']))
    if ((not allow_partial and dc['status']!='complete') or sc['status']!='complete' or collect_grid(verify(sc['manifest']))['rows']!=sc['rows'] or
            dm['minimum_preparation']!=record(prepared) or sm['minimum_preparation']!=record(prepared)):
        raise InvalidArtifact('independent validation results incomplete/mismatched')
    rows={};scores={}
    for key,item in top['states'].items():
        if dc['rows'][key]['status']!='complete':continue
        name,metal=item['case_id'],item['metal'];center=p['centers'][name][metal]
        dft_change=(dc['rows'][key]['energy_hartree']-center['DFT_energy_Ha'])*HA_TO_KCAL
        core=sc['rows']['core__'+key];full=sc['rows']['full__'+key]
        j_change=((full['energy_eV']-center['full_short']['energy_eV'])-(core['energy_eV']-center['core_short']['energy_eV']))*EV_TO_KCAL
        actual=dft_change+j_change;pred=item['prediction']['predicted_response_kcal_mol'];u=np.array(item['prediction']['proposed_displacement_A'])
        nonlinear=dft_change-float(np.array(center['DFT_gradient_kcal_mol_A'])@u);tolerance=max(.05,.25*abs(nonlinear))
        # Analytic derivative of the actual native DFT+short energy, not the cheap model.
        g=read_engrad(verify(dc['gradients'][key]['artifacts']['engrad']))
        gd=np.array(g['gradient_Ha_per_bohr'][0])*HA_TO_KCAL/BOHR_TO_A
        gj=(-np.load(verify(full['forces']))[center['full_metal_index']]+np.load(verify(core['forces']))[0])*EV_TO_KCAL
        residual=float(np.linalg.norm(gd+gj));initial=float(np.linalg.norm(np.array(center['DFT_gradient_kcal_mol_A'])+center['J_gradient_kcal_mol_A']))
        force_tolerance=max(2.,.25*initial)
        checks={'energy_prediction':abs(pred-actual)<=tolerance,'lower_actual_energy':actual<0,'stationarity':residual<=force_tolerance}
        rows[key]={'prediction':item['prediction'],'DFT_change_kcal_mol':dft_change,'J_change_kcal_mol':j_change,
                   'actual_energy_change_kcal_mol':actual,'prediction_error_kcal_mol':pred-actual,
                   'energy_tolerance_kcal_mol':tolerance,'actual_gradient_kcal_mol_A':(gd+gj).tolist(),
                   'gradient_norm_kcal_mol_A':residual,'gradient_tolerance_kcal_mol_A':force_tolerance,
                   'checks':checks,'pass':all(checks.values()),'validated_response_kcal_mol':actual if all(checks.values()) else None,
                   'DFT_receipt':dc['rows'][key]['receipt']}
    for name,c in a['scores'].items():
        vals={m:rows.get(name+'_'+m,{}).get('validated_response_kcal_mol') for m in ('La','Ca')}
        raw={m:rows.get(name+'_'+m,{}).get('actual_energy_change_kcal_mol') for m in ('La','Ca')}
        delta=None if None in vals.values() else vals['Ca']-vals['La']
        exploratory=None if None in raw.values() else raw['Ca']-raw['La']
        scores[name]={'baseline_R_kcal_mol':c['baseline_R_kcal_mol'],'endpoint_validated_delta_R_kcal_mol':delta,
                     'endpoint_validated_R_kcal_mol':None if delta is None else c['baseline_R_kcal_mol']+delta,
                     'validated_delta_R_kcal_mol':None,'validated_R_kcal_mol':None,
                     'exploratory_actual_delta_R_kcal_mol':exploratory,
                     'exploratory_actual_R_kcal_mol':None if exploratory is None else c['baseline_R_kcal_mol']+exploratory,
                     'calibrated_class':None,'entropy_correction_kcal_mol':None}
    x,y=(scores[k]['endpoint_validated_delta_R_kcal_mol'] for k in ('GGR_extended','GGR_connected'))
    partition={'difference_kcal_mol':None if x is None or y is None else y-x,
               'pass':None if x is None or y is None else abs(y-x)<=2.}
    for row in scores.values():
        qualified=partition['pass'] is True and row['endpoint_validated_delta_R_kcal_mol'] is not None
        row['qualification_status']='qualified_development_response' if qualified else 'response_model_not_validated'
        if qualified:
            row['validated_delta_R_kcal_mol']=row['endpoint_validated_delta_R_kcal_mol']
            row['validated_R_kcal_mol']=row['endpoint_validated_R_kcal_mol']
    contrasts={}
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        contrasts[alpha]={field:None if any(scores[k][field] is None for k in (alpha,'GGR_extended')) else scores[alpha][field]-scores['GGR_extended'][field]
                          for field in ('baseline_R_kcal_mol','exploratory_actual_R_kcal_mol','validated_R_kcal_mol')}
    result={'status':'complete' if len(rows)==len(top['states']) else 'incomplete',
            'pending_DFT_endpoints':sorted(set(top['states'])-set(rows)),'sources':{k:record(v) for k,v in [('prepared',prepared),('DFT',dft),('short',short)]},
            'implementation':record(__file__),'rows':rows,'excluded':top['excluded'],'scores':scores,'contrasts':contrasts,
            'partition_check':partition,'production_score_available':False,'all_cases_consumed_development':True,
            'validation_scope':'energy and analytic metal-gradient at fixed predicted displacements; no native DFT Hessian or global basin proof'}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    reporter=out/'report_source.py';reporter.write_bytes(Path(__file__).read_bytes())
    result['implementation']=record(reporter);write_new(out/'result.json',result)
    return {'status':result['status'],'pending_DFT_endpoints':result['pending_DFT_endpoints'],'endpoint_checks':{k:v['checks'] for k,v in rows.items()},'scores':scores,'contrasts':contrasts,'partition':partition}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for command,keys in [('prepare',('assessment','output')),('validate',('manifest',)),('report',('prepared','dft','short','output'))]:
        q=sub.add_parser(command)
        for key in keys:q.add_argument('--'+key,required=True)
        if command=='report':q.add_argument('--allow-partial',action='store_true')
    args=vars(parser.parse_args());command=args.pop('command')
    print(json.dumps({'prepare':prepare,'validate':validate,'report':report}[command](**args),indent=2))
