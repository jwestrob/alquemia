"""Fixed-radius metal response, with actual GGR reuse and independent alpha DFT."""
from __future__ import annotations
import argparse
import copy
import json
import time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,BOHR_TO_A,read_json,record,verify,write_new,xyz,cache_key
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_metal_response import LIMIT,TOL,UNAVAILABLE,shift,prediction as free_prediction,collect as collect_short
from mace_metal_minimum import same_prediction

POLICY='DFT_anchored_MACE_GB_bounded_metal_response_v1'
SCHEMA='alquemia.mace_bounded_response_short.v1'


def ball(k,g,stored_interior=None):
    """Solve a positive quadratic on the declared sphere, in physical Angstroms."""
    k=np.asarray(k,dtype=float);g=np.asarray(g,dtype=float)
    if k.shape!=(3,3) or g.shape!=(3,) or not np.isfinite(k).all() or not np.isfinite(g).all() or not np.allclose(k,k.T,atol=1e-10,rtol=0):
        raise InvalidArtifact('finite symmetric physical3x3 model required')
    eig=np.linalg.eigvalsh(k)
    if eig.min()<=0:return {'status':'unstable_curvature','eigenvalues':eig.tolist()}
    u=-np.linalg.solve(k,g);lam=0.;calls=1;active=np.linalg.norm(u)>LIMIT
    if active:
        low,high=0.,1.
        while calls<200:
            trial=-np.linalg.solve(k+high*np.eye(3),g);calls+=1
            if np.linalg.norm(trial)<=LIMIT:break
            low,high=high,high*2
        else:raise InvalidArtifact('quadratic multiplier bracket failed')
        while calls<200:
            lam=(low+high)/2;u=-np.linalg.solve(k+lam*np.eye(3),g);calls+=1
            norm=float(np.linalg.norm(u))
            if abs(norm-LIMIT)<=1e-12:break
            if norm>LIMIT:low=lam
            else:high=lam
        if abs(np.linalg.norm(u)-LIMIT)>1e-10:raise InvalidArtifact('sphere equation did not converge')
    elif stored_interior is not None:
        old=np.asarray(stored_interior)
        if old.shape!=(3,) or not np.allclose(old,u,atol=1e-10,rtol=1e-10):
            raise InvalidArtifact('reused interior point differs from physical solution')
        u=old.copy()
    residual=float(np.linalg.norm((k+lam*np.eye(3))@u+g)/max(np.linalg.norm(g),np.finfo(float).tiny))
    if residual>1e-10:raise InvalidArtifact('quadratic equation residual failed')
    return {'status':'solved','displacement_A':u.tolist(),'displacement_norm_A':float(np.linalg.norm(u)),
            'boundary_active':bool(active),'lambda_kcal_mol_A2':lam,'linear_solves':calls,
            'equation_relative_residual':residual,'predicted_change_kcal_mol':float(g@u+.5*u@k@u),
            'eigenvalues_kcal_mol_A2':eig.tolist()}


def predict(state):
    g=np.array(state['DFT_gradient_kcal_mol_A'])+state['J_gradient_kcal_mol_A']
    old=state['prediction'];actual=free_prediction(g,state['relative_grid_energies_kcal_mol'])
    if not same_prediction(old,actual):raise InvalidArtifact('parent grid/gradient prediction changed')
    parts={key:{s:np.array(k) for s,k in value.items()} for key,value in old['matrices_kcal_mol_A2'].items()}
    fine=sum(v['fine'] for v in parts.values());coarse=sum(v['coarse'] for v in parts.values())
    result=ball(fine,g,old['proposed_displacement_A']);other=ball(coarse,g)
    if result['status']!='solved' or other['status']!='solved':
        return {'status':'unstable_curvature','fine':result,'coarse':other}
    u=np.array(result['displacement_A']);errors={key:float(.5*u@(v['fine']-v['coarse'])@u) for key,v in parts.items()}
    errors['combined']=float(.5*u@(fine-coarse)@u)
    delta=float(np.linalg.norm(u-other['displacement_A']))
    result.update(refinement_energy_changes_kcal_mol=errors,refinement_displacement_change_A=delta,
                  gradient_kcal_mol_A=g.tolist(),status='eligible_for_native_check' if max(map(abs,errors.values()))<=.02 and delta<=.01 else 'curvature_not_converged')
    return result


def sources(assessment,reuse_result,full_receipt_check=False):
    from mace_curvature import actual_collection
    a=read_json(assessment);p=read_json(verify(a['preparation']));reuse=read_json(reuse_result)
    verify(reuse['implementation']);rp=read_json(verify(reuse['sources']['prepared']))
    if (a['status']!='complete' or not a['predictions_are_not_validated_scores'] or
            rp['assessment']!=record(assessment) or reuse['status']!='complete' or
            set(reuse['rows'])!={'GGR_extended_Ca','GGR_extended_La','GGR_connected_Ca','GGR_connected_La'} or
            not reuse['partition_check']['pass']):raise InvalidArtifact('exact complete parent/GGR physical validation required')
    dc=read_json(verify(reuse['sources']['DFT']))
    verify(dc['manifest'])
    for row in reuse['rows'].values():
        receipt=read_json(verify(row['DFT_receipt']))
        for artifact in receipt['artifacts'].values():verify(artifact)
    if full_receipt_check:
        from mace_mechanics_assess import verified_new_dft
        verified_new_dft(verify(reuse['sources']['DFT']))
    actual_collection(verify(reuse['sources']['short']))
    for name,row in a['rows'].items():
        for metal,state in row.items():
            key=name+'_'+metal;pred=predict(state)
            if name.startswith('GGR'):
                original=rp['states'][key]
                if pred.get('status')!='eligible_for_native_check' or pred['boundary_active'] or pred['displacement_A']!=original['prediction']['proposed_displacement_A']:
                    raise InvalidArtifact('GGR interior solution changed; archive cannot satisfy this task')
                for kind in ('core','full'):
                    src=p['centers'][name][metal][kind+'_state'];wanted,_=shift(xyz(verify(src['xyz'])),metal,pred['displacement_A'])
                    archived=xyz(verify(original['endpoints'][kind]['xyz']))
                    if [x[0] for x in wanted]!=[x[0] for x in archived] or not np.allclose([x[1:] for x in wanted],[x[1:] for x in archived],atol=5e-10,rtol=0):
                        raise InvalidArtifact('reused GGR geometry differs')
                    if original['endpoints'][kind]['charge']!=src['charge'] or original['endpoints'][kind]['spin_multiplicity']!=src['spin_multiplicity']:
                        raise InvalidArtifact('reused GGR charge/spin differs')
    return a,p,reuse,rp


def prepare(assessment,reuse_result,agreement,output):
    from mace_mechanics import physical_model
    from ggr_sensitivity import membership,METHOD,ORCA
    start=time.monotonic();cpu=time.process_time();a,p,reuse,rp=sources(assessment,reuse_result,full_receipt_check=True)
    ancestor=read_json(verify(p['assessment']));oldp=read_json(verify(ancestor['sources']['prepared']))
    pc=read_json(verify(ancestor['sources']['short']));parent=read_json(verify(pc['manifest']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);qr=out/'quantum';sr=out/'short';qr.mkdir();sr.mkdir()
    pins=snapshot(sr,(*parent['implementation'],'mace_bounded_response.py','mace_metal_response.py','mace_metal_minimum.py'))
    predictions={};states={};excluded={};tasks=[];quantum=[];reused={}
    for name,row in a['rows'].items():
        for metal,state in row.items():
            key=name+'_'+metal;pred=predict(state);predictions[key]=pred
            if pred['status']!='eligible_for_native_check':excluded[key]={'reason':pred['status']};continue
            if name.startswith('GGR'):
                reused[key]={'validation':record(reuse_result),'state':rp['states'][key],'prediction':pred}
                continue
            repair=read_json(verify(p['cases'][name]['source_preparation']))
            graph,coords,pos,mk,_,_=physical_model(repair);u=np.array(pred['displacement_A'])
            before=membership(graph,coords,pos,mk);after=membership(graph,coords,pos+u,mk)
            if before!=after:
                excluded[key]={'reason':'typed_donor_membership_changed','before':before,'after':after};continue
            directory=out/key;directory.mkdir();states[key]={'case_id':name,'metal':metal,'prediction':pred,'endpoints':{}}
            for kind in ('core','full'):
                source=p['centers'][name][metal][kind+'_state'];rows,index=shift(xyz(verify(source['xyz'])),metal,u)
                if source['spin_multiplicity']!=1:raise InvalidArtifact('only unchanged singlet endpoints supported')
                xp=directory/(kind+'.xyz');write_xyz(xp,rows)
                endpoint={'xyz':record(xp),'charge':source['charge'],'spin_multiplicity':1,'selected_metal_index':index,
                          'state':check_atoms(rows,source['charge'])}
                states[key]['endpoints'][kind]=endpoint
                tasks.append({'task_id':kind+'__'+key,'case_id':name,'metal':metal,'kind':kind,'variant':'bounded_metal_response',
                              'energy_component':'interaction_energy','displacement_A':u.tolist(),**endpoint})
            td=qr/key;td.mkdir();xp=td/'core.xyz';xp.write_bytes(verify(states[key]['endpoints']['core']['xyz']).read_bytes())
            charge=states[key]['endpoints']['core']['charge'];ip=td/'endpoint.inp'
            ip.write_text(f'! {METHOD} EnGrad\n* xyzfile {charge} 1 core.xyz\n')
            settings={'protocol_id':POLICY,'parent_assessment':record(assessment),'agreement':record(agreement),
                      'source_preparation':p['cases'][name]['source_preparation'],'prediction':pred,
                      'xyz':record(xp),'input':record(ip),'charge':charge,'spin_multiplicity':1,'metal':metal,'case_id':name}
            quantum.append({'task_id':key,'case_id':name,'point':'bounded_metal_response','metal':metal,'charge':charge,'multiplicity':1,
                            'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),'engrad_path':str(td/'endpoint.engrad'),
                            'task_type':'analytic_gradient','scientific_settings':settings,'cache_key':cache_key(settings)})
    top={'policy_id':POLICY,'agreement':record(agreement),'parent_assessment':record(assessment),'reuse_result':record(reuse_result),
         'physical_preparation':a['preparation'],'predictions':predictions,'states':states,'excluded':excluded,'reused':reused,
         'implementation':pins,'new_DFT_calls':len(quantum),'new_short_calls':len(tasks),'radius_A':LIMIT,
         'preparation_wall_seconds':time.monotonic()-start,'preparation_CPU_seconds':time.process_time()-cpu,**UNAVAILABLE}
    write_new(out/'preparation.json',top)
    if not quantum:return {'status':'no_eligible_new_states','preparation':record(out/'preparation.json'),'excluded':excluded}
    qm=copy.deepcopy(read_json(verify(oldp['DFT_tasks'])))
    qm.update(schema_version='alquemia.mace_bounded_response_DFT.v1',protocol_id=POLICY,tasks=quantum,implementation=pins,
              bounded_preparation=record(out/'preparation.json'),agreement=record(agreement),orca=record(ORCA))
    write_new(qr/'manifest.json',qm)
    model=copy.deepcopy(parent['model']);model['preparation_policy']=POLICY
    m={'schema_version':SCHEMA,'protocol_id':POLICY+'_short','kind':'bounded_short','agreement':record(agreement),
       'bounded_preparation':record(out/'preparation.json'),'parent_short_collection':ancestor['sources']['short'],
       'numerical_reference':parent['numerical_reference'],'software':parent['software'],'model':model,
       'implementation':pins,'tasks':tasks,'tolerances':TOL,**UNAVAILABLE}
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':model,'software':m['software'],'implementation':pins})
    write_new(sr/'manifest.json',m)
    from affordable_workflow import dry_run
    return {'status':'prepared','short':validate(sr/'manifest.json'),'DFT':dry_run(qr/'manifest.json'),
            'excluded':excluded,'reused_GGR_endpoints':len(reused),'preparation':record(out/'preparation.json')}


def validate(manifest):
    m=read_json(manifest);top=read_json(verify(m['bounded_preparation']))
    a,p,reuse,rp=sources(verify(top['parent_assessment']),verify(top['reuse_result']))
    pc=read_json(verify(m['parent_short_collection']));parent=read_json(verify(pc['manifest']))
    model=copy.deepcopy(parent['model']);model['preparation_policy']=POLICY
    if m['schema_version']!=SCHEMA or top['policy_id']!=POLICY or top['radius_A']!=LIMIT or m['model']!=model or m['software']!=parent['software'] or m['tolerances']!=TOL or numerical_parent_gate(m)['status']!='pass':
        raise InvalidArtifact('bounded physical method/radius changed')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],*read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    wanted={name+'_'+metal:predict(state) for name,row in a['rows'].items() for metal,state in row.items()}
    if set(top['states'])|set(top['excluded'])|set(top['reused'])!=set(wanted) or not same_prediction(top['predictions'],wanted):
        raise InvalidArtifact('prediction population differs from actual parent model')
    expected={kind+'__'+key:(kind,key) for key in top['states'] for kind in ('core','full')}
    if {t['task_id'] for t in m['tasks']}!=set(expected) or len(m['tasks'])!=2*len(top['states']) or len(m['tasks'])>8:
        raise InvalidArtifact('finite alpha task inventory changed')
    for t in m['tasks']:
        kind,key=expected[t['task_id']];state=top['states'][key];name,metal=state['case_id'],state['metal'];pred=wanted[key]
        if name not in ('ALPHA_1F6S','ALPHA_6IP9') or pred['status']!='eligible_for_native_check' or not same_prediction(state['prediction'],pred):
            raise InvalidArtifact('unsupported alpha preparation/prediction')
        for k,v in {'case_id':name,'metal':metal,'kind':kind,'variant':'bounded_metal_response','energy_component':'interaction_energy',**state['endpoints'][kind]}.items():
            if t.get(k)!=v:raise InvalidArtifact('task identity/component/state differs')
        if not np.allclose(t['displacement_A'],pred['displacement_A'],atol=1e-10,rtol=1e-10):raise InvalidArtifact('frozen displacement changed')
        base=xyz(verify(p['centers'][name][metal][kind+'_state']['xyz']));rows=xyz(verify(t['xyz']));desired,index=shift(base,metal,t['displacement_A'])
        if index!=t['selected_metal_index'] or [r[0] for r in rows]!=[r[0] for r in desired] or not np.allclose([r[1:] for r in rows],[r[1:] for r in desired],atol=5e-10,rtol=0):
            raise InvalidArtifact('physical geometry differs from the fixed bounded point')
        check_atoms(rows,t['charge']);payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('bounded task cache key differs')
    return {'status':'pass','manifest':record(manifest),'new_short_calls':len(m['tasks']),'new_DFT_calls':len(top['states'])}


def endpoint_checks(pred,center,dft_change,j_change,actual_gradient):
    g=np.asarray(actual_gradient);u=np.array(pred['displacement_A']);actual=dft_change+j_change
    initial=np.linalg.norm(np.array(center['DFT_gradient_kcal_mol_A'])+center['J_gradient_kcal_mol_A'])
    nonlinear=dft_change-float(np.array(center['DFT_gradient_kcal_mol_A'])@u)
    error=pred['predicted_change_kcal_mol']-actual;etol=max(.05,.25*abs(nonlinear));gtol=max(2.,.25*initial)
    radial=float(u@g);norm=float(np.linalg.norm(g));tangent=g-u*radial/(u@u) if u@u>0 else g
    checks={'energy_prediction':abs(error)<=etol,'lower_actual_energy':actual<0}
    if pred['boundary_active']:
        checks.update(boundary_radial_sign=radial<=0,boundary_tangent_gradient=np.linalg.norm(tangent)<=gtol)
    else:checks['interior_gradient']=norm<=gtol
    return {'prediction':pred,'DFT_change_kcal_mol':dft_change,'J_change_kcal_mol':j_change,
            'actual_energy_change_kcal_mol':actual,'prediction_error_kcal_mol':error,'energy_tolerance_kcal_mol':etol,
            'actual_gradient_kcal_mol_A':g.tolist(),'gradient_norm_kcal_mol_A':norm,
            'tangent_gradient_norm_kcal_mol_A':float(np.linalg.norm(tangent)),'gradient_tolerance_kcal_mol_A':float(gtol),
            'radial_gradient_dot_displacement_kcal_mol':radial,'checks':{k:bool(v) for k,v in checks.items()},
            'pass':all(checks.values()),'endpoint_validated_response_kcal_mol':actual if all(checks.values()) else None}


def report(prepared,dft,short,output):
    from mace_mechanics_assess import verified_new_dft
    from affordable_response import read_engrad
    start=time.monotonic();cpu=time.process_time();top=read_json(prepared)
    a,p,reuse,rp=sources(verify(top['parent_assessment']),verify(top['reuse_result']),full_receipt_check=True)
    dc,dm=verified_new_dft(dft);sc=read_json(short);sm=read_json(verify(sc['manifest']));validate(verify(sc['manifest']))
    if dc['status']!='complete' or sc['status']!='complete' or collect_short(verify(sc['manifest']))['rows']!=sc['rows'] or dm['bounded_preparation']!=record(prepared) or sm['bounded_preparation']!=record(prepared):
        raise InvalidArtifact('actual bounded validation results incomplete/mismatched')
    rows={}
    for key,item in top['reused'].items():
        old=reuse['rows'][key];name=item['state']['case_id'];metal=item['state']['metal'];center=p['centers'][name][metal]
        row=endpoint_checks(item['prediction'],center,old['DFT_change_kcal_mol'],old['J_change_kcal_mol'],old['actual_gradient_kcal_mol_A'])
        row.update(execution_status='exact_archive_reuse',source_validation=top['reuse_result'],DFT_receipt=old['DFT_receipt']);rows[key]=row
    for key,item in top['states'].items():
        name,metal=item['case_id'],item['metal'];center=p['centers'][name][metal]
        dft_change=(dc['rows'][key]['energy_hartree']-center['DFT_energy_Ha'])*HA_TO_KCAL
        full=sc['rows']['full__'+key];core=sc['rows']['core__'+key]
        j_change=((full['energy_eV']-center['full_short']['energy_eV'])-(core['energy_eV']-center['core_short']['energy_eV']))*EV_TO_KCAL
        grad=read_engrad(verify(dc['gradients'][key]['artifacts']['engrad']))
        gd=np.array(grad['gradient_Ha_per_bohr'][0])*HA_TO_KCAL/BOHR_TO_A
        gj=(-np.load(verify(full['forces']))[center['full_metal_index']]+np.load(verify(core['forces']))[0])*EV_TO_KCAL
        row=endpoint_checks(item['prediction'],center,dft_change,j_change,gd+gj)
        row.update(execution_status='new_native_validation',DFT_receipt=dc['rows'][key]['receipt']);rows[key]=row
    scores={}
    for name,base in a['scores'].items():
        ends={m:rows.get(name+'_'+m,{}).get('endpoint_validated_response_kcal_mol') for m in ('Ca','La')}
        actual={m:rows.get(name+'_'+m,{}).get('actual_energy_change_kcal_mol') for m in ('Ca','La')}
        validated=None if None in ends.values() else ends['Ca']-ends['La']
        raw=None if None in actual.values() else actual['Ca']-actual['La']
        scores[name]={'baseline_R_kcal_mol':base['baseline_R_kcal_mol'],'endpoint_validated_delta_R_kcal_mol':validated,
                      'actual_exploratory_delta_R_kcal_mol':raw,'actual_exploratory_R_kcal_mol':None if raw is None else base['baseline_R_kcal_mol']+raw,
                      'qualified_R_kcal_mol':None,'qualified_delta_R_kcal_mol':None,'calibrated_class':None,'entropy_correction_kcal_mol':None}
    x,y=(scores[k]['endpoint_validated_delta_R_kcal_mol'] for k in ('GGR_extended','GGR_connected'))
    partition={'difference_kcal_mol':None if x is None or y is None else y-x,'pass':None if x is None or y is None else abs(y-x)<=2.}
    for row in scores.values():
        if partition['pass'] is True and row['endpoint_validated_delta_R_kcal_mol'] is not None:
            row['qualified_delta_R_kcal_mol']=row['endpoint_validated_delta_R_kcal_mol'];row['qualified_R_kcal_mol']=row['baseline_R_kcal_mol']+row['qualified_delta_R_kcal_mol']
    comparisons={}
    for name in ('ALPHA_1F6S','ALPHA_6IP9'):
        r={k:None if any(scores[n][k] is None for n in (name,'GGR_extended')) else scores[name][k]-scores['GGR_extended'][k]
           for k in ('baseline_R_kcal_mol','actual_exploratory_R_kcal_mol','qualified_R_kcal_mol')}
        r['qualified_direction_pass']=None if r['qualified_R_kcal_mol'] is None else r['qualified_R_kcal_mol']>0
        r['actual_margin_change_kcal_mol']=None if r['actual_exploratory_R_kcal_mol'] is None else r['actual_exploratory_R_kcal_mol']-r['baseline_R_kcal_mol']
        comparisons[name]=r
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);source=out/'report_source.py';source.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'sources':{k:record(v) for k,v in [('prepared',prepared),('DFT',dft),('short',short)]},
            'implementation':record(source),'rows':rows,'excluded':top['excluded'],'scores':scores,'partition_check':partition,
            'comparisons':comparisons,'validation_scope':'conditional bounded response with fixed exterior; boundary KKT is not unconstrained stationarity',
            'absolute_reference':None,'all_cases_consumed_development':True,'production_baseline_changed':False,
            'report_wall_seconds':time.monotonic()-start,'report_CPU_seconds':time.process_time()-cpu}
    write_new(out/'result.json',result)
    return {'status':'complete','endpoint_checks':{k:v['checks'] for k,v in rows.items()},'scores':scores,'comparisons':comparisons,'partition':partition}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for name,keys in [('prepare',('assessment','reuse-result','agreement','output')),('validate',('manifest',)),('report',('prepared','dft','short','output'))]:
        q=sub.add_parser(name)
        for key in keys:q.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    print(json.dumps({'prepare':prepare,'validate':validate,'report':report}[command](**args),indent=2))
