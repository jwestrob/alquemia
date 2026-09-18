"""Native validation at frozen matched-H vacuum hybrid response predictions."""
from __future__ import annotations
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,BOHR_TO_A,read_json,record,verify,write_new,xyz,cache_key
from affordable_response import read_engrad
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,check_atoms
from mace_omol_hybrid_response import POLICY,CASES,sources,task,collect as collect_model
from mace_omol_gradient_run import CONFIG,check_parent,model as gradient_model
from mace_omol_vacuum import METHOD,scientific_input,parse_endpoint
from mace_metal_minimum import same_prediction
from mace_bounded_response import ball,endpoint_checks
STAGE='matched_hybrid_response_validation'


@cached_file_checks
def inputs(assessment,full_check=False):
    a=read_json(assessment);verify(a['implementation']);prep=read_json(verify(a['prepared']))
    p,r,gr,gm,states,fulls,reused=sources(verify(prep['source_preparation']),verify(prep['static_report']),verify(prep['ggr_gradients']),full_check)
    if a['status']!='complete' or a['protocol_id']!=POLICY or not a['predictions_are_not_validated_scores'] or set(a['rows'])!=set(CASES):raise InvalidArtifact('complete declared predictions required')
    selected={}
    for name,row in a['rows'].items():
        for metal,s in row.items():
            pred=s['prediction'];k={key:np.array(v) for key,v in s['matrices_kcal_mol_A2'].items()};g=np.array(s['gradient_kcal_mol_A'])
            actual=ball(k['fine'],g);coarse=ball(k['coarse'],g)
            if actual['status']=='solved' and coarse['status']=='solved':
                u=np.array(actual['displacement_A']);error=float(.5*u@(k['fine']-k['coarse'])@u);delta=float(np.linalg.norm(u-coarse['displacement_A']))
                actual.update(refinement_error_kcal_mol=error,refinement_displacement_A=delta,status='eligible_for_native_check' if abs(error)<=.02 and delta<=.01 and all(v['pass'] for v in a['odd_checks'][states[name][metal]['global_id']][metal]) else 'numerical_response_not_qualified')
            if not same_prediction(actual,pred):raise InvalidArtifact('prediction differs from declared physical model')
            if pred['status']=='eligible_for_native_check':selected[name+'_'+metal]={'case_id':name,'metal':metal,'prediction':pred}
    return a,prep,p,r,gm,states,selected


@cached_file_checks
def prepare(assessment,output):
    from mace_omol import common,seal
    from mace_global_benchmark import snapshot
    from mace_mechanics import physical_model
    from ggr_sensitivity import membership,ORCA
    start=time.monotonic();cpu=time.process_time();a,initial,p,r,gm,centers,selected=inputs(assessment,True)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    _,md,m=common(verify(p['inventory']),verify(p['software']),verify(initial['agreement']),out/'mace',STAGE)
    states={};excluded={};tasks=[];quantum=[]
    unavailable={name+'_'+metal:{'reason':state['prediction']['status'],'prediction':state['prediction']}
                 for name,row in a['rows'].items() for metal,state in row.items()
                 if name+'_'+metal not in selected}
    qr=out/'quantum';qr.mkdir()
    parent_quantum=read_json(verify(r['quantum']['manifest']))
    qpins=snapshot(qr,(*m['implementation'],*parent_quantum['implementation']))
    for key,item in selected.items():
        name,metal=item['case_id'],item['metal'];c=read_json(verify(p['cases'][name]));original=read_json(verify(c['source_mapping']));repair=read_json(verify(original['source_preparation']))
        graph,coords,pos,mk,_,_=physical_model(repair);u=item['prediction']['displacement_A']
        if membership(graph,coords,pos,mk)!=membership(graph,coords,pos+u,mk):excluded[key]={'reason':'typed_donor_membership_changed','prediction':item['prediction']};continue
        s=centers[name][metal];states[key]={**item,'source_case':p['cases'][name],'endpoints':{}}
        for kind in ('core','full'):
            t=task(s[kind],name,metal,kind,'bounded',u,md);t['energy_only']=False;tasks.append(t)
            states[key]['endpoints'][kind]={'xyz':t['xyz'],'charge':t['charge'],'spin_multiplicity':1,'metal_index':t['metal_index'],'state':t['state']}
        d=qr/key;d.mkdir();xp=d/'core.xyz';xp.write_bytes(verify(states[key]['endpoints']['core']['xyz']).read_bytes());ip=d/'endpoint.inp';ip.write_text(scientific_input(s['core']['charge']))
        t={'task_id':key,'case_id':name,'metal':metal,'charge':s['core']['charge'],'multiplicity':1,
           'input':record(ip),'xyz':record(xp),'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),
           'task_type':'analytic_gradient','source_mapping':p['cases'][name],'displacement_A':u,'assessment':record(assessment)}
        t['cache_key']=cache_key({'task':t,'method':METHOD,'protocol':POLICY,'implementation':qpins});quantum.append(t)
    # Preserve the existing matched-H runner's largest-core-first task policy.
    quantum.sort(key=lambda t:(-len(xyz(verify(t['xyz']))),t['task_id']))
    top={'protocol_id':POLICY,'assessment':record(assessment),'states':states,'excluded':excluded,'agreement':initial['agreement'],
         'new_DFT_calls':len(quantum),'new_MACE_gradient_calls':len(tasks),'implementation':m['implementation'],
         'unavailable_predictions':unavailable,
         'preparation_wall_seconds':time.monotonic()-start,'preparation_CPU_seconds':time.process_time()-cpu}
    write_new(out/'preparation.json',top)
    if not states:return {'status':'no_supported_native_points','preparation':record(out/'preparation.json'),'excluded':excluded,'unavailable_predictions':unavailable}
    m.update(protocol_id=POLICY,model=gradient_model(verify(p['software'])),gradient_settings=CONFIG,core_report=gm['core_report'],
             minimum_preparation=record(out/'preparation.json'),tasks=tasks)
    check_parent(m);seal(md,m)
    q={'protocol_id':POLICY,'stage':STAGE,'minimum_preparation':record(out/'preparation.json'),'agreement':initial['agreement'],
       'method':METHOD,'energy_scope':'isolated_vacuum_endpoint','implementation':qpins,'tasks':quantum,'orca':record(ORCA),
       'execution_policy':{'task_runner':qpins['run_orca_task_manifest.py'],'runtime_renderer':qpins['render_orca_runtime_input.py']},
       'cost_tracking':{'compute_budget':None,'wall_time_limit':None}}
    write_new(qr/'manifest.json',q)
    return {'status':'prepared','DFT':validate_quantum(qr/'manifest.json'),'MACE':validate(md/'manifest.json'),'excluded':excluded}


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);top=read_json(verify(m['minimum_preparation']));a,initial,p,r,gm,centers,selected=inputs(verify(top['assessment']))
    if (m['protocol_id']!=POLICY or m['stage']!=STAGE or m['model']!=gradient_model(verify(p['software'])) or m['software']!=p['software'] or m['gradient_settings']!=CONFIG or
            set(top['states'])|set(top['excluded'])!=set(selected) or len(m['tasks'])!=2*len(top['states']) or len(m['tasks'])>16):raise InvalidArtifact('fixed validation method or inventory differs')
    check_parent(m)
    for ref in [m['agreement'],*m['implementation'].values()]:verify(ref)
    expected={kind+'__'+key+'__bounded':(kind,key) for key in top['states'] for kind in ('core','full')}
    if {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('validation task identity differs')
    for t in m['tasks']:
        kind,key=expected[t['task_id']];item=selected[key];s=centers[item['case_id']][item['metal']][kind];pred=top['states'][key]['prediction']
        if not same_prediction(pred,item['prediction']):raise InvalidArtifact('frozen displacement prediction differs')
        if t['energy_only'] or t['metal']!=item['metal'] or t['case_id']!=item['case_id'] or t['source_xyz']!=s['xyz'] or t['kind']!=kind or t['descriptor_gradient_experiment']!=POLICY or t['derivative_backend']!='checkpointed':raise InvalidArtifact('validation gradient/state identity differs')
        u=np.array(item['prediction']['displacement_A']);old=xyz(verify(s['xyz']));new=xyz(verify(t['xyz']));i=s['metal_index']
        want=np.array([row[1:] for row in old]);want[i]+=u
        if [x[0] for x in new]!=[x[0] for x in old] or not np.allclose([x[1:] for x in new],want,atol=5e-10,rtol=0) or t['charge']!=s['charge'] or t['spin_multiplicity']!=1 or t['metal_index']!=i:raise InvalidArtifact('physical geometry/charge changed')
        if check_atoms(new,t['charge'])!=t['state']:raise InvalidArtifact('parity or atom inventory differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('validation cache differs')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


@cached_file_checks
def validate_quantum(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);top=read_json(verify(m['minimum_preparation']));a,initial,p,r,gm,centers,selected=inputs(verify(top['assessment']))
    if m['protocol_id']!=POLICY or m['stage']!=STAGE or m['method']!=METHOD or len(m['tasks'])>8 or {t['task_id'] for t in m['tasks']}!=set(top['states']):raise InvalidArtifact('declared quantum inventory/method differs')
    for ref in [m['agreement'],m['orca'],*m['implementation'].values()]:verify(ref)
    for t in m['tasks']:
        item=top['states'][t['task_id']];state=item['endpoints']['core']
        if t['task_id'] not in selected or not same_prediction(item['prediction'],selected[t['task_id']]['prediction']):raise InvalidArtifact('quantum task is not the fixed eligible prediction')
        if t['charge']!=state['charge'] or t['xyz']['sha256']!=state['xyz']['sha256'] or verify(t['input']).read_text()!=scientific_input(t['charge']):raise InvalidArtifact('quantum input or physical state differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'method':METHOD,'protocol':POLICY,'implementation':m['implementation']}):raise InvalidArtifact('quantum cache differs')
    return dry_run(manifest)


@cached_file_checks
def collect_quantum(manifest):
    from ggr_sensitivity import executed
    validate_quantum(manifest);m,rows=executed(manifest)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:r.update(parse_endpoint(t,verify(r['output']),t['engrad_path']))
        except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    return {'status':'complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete','manifest':record(manifest),'rows':rows}


@cached_file_checks
def report(prepared,quantum,mace,output):
    start=time.monotonic();cpu=time.process_time();top=read_json(prepared);a,initial,p,r,gm,centers,selected=inputs(verify(top['assessment']),True)
    qc=collect_quantum(quantum);validate(mace);mc=collect_model(mace)
    if qc['status']!='complete' or mc['status']!='complete':raise InvalidArtifact('actual native validation incomplete')
    if read_json(quantum)['minimum_preparation']!=record(prepared) or read_json(mace)['minimum_preparation']!=record(prepared):raise InvalidArtifact('native validation belongs to another preparation')
    rows={};scores={}
    for key,item in top['states'].items():
        name,metal=item['case_id'],item['metal'];s=centers[name][metal];old=s['DFT_source'];new=qc['rows'][key]
        full=mc['rows']['full__'+key+'__bounded'];core=mc['rows']['core__'+key+'__bounded']
        dft=(new['energy_hartree']-old['energy_hartree'])*HA_TO_KCAL
        j=((full['energy_eV']-s['full']['archived_center_energy_eV'])-(core['energy_eV']-s['core']['archived_center_energy_eV']))*EV_TO_KCAL
        gd=np.array(new['gradient_kcal_mol_per_A'][0]);gj=(np.load(verify(full['gradient']))[s['full']['metal_index']]-np.load(verify(core['gradient']))[0])*EV_TO_KCAL
        center={'DFT_gradient_kcal_mol_A':a['rows'][name][metal]['DFT_gradient_kcal_mol_A'],'J_gradient_kcal_mol_A':a['rows'][name][metal]['context_gradient_kcal_mol_A']}
        row=endpoint_checks(item['prediction'],center,dft,j,gd+gj);row.update(DFT_receipt=new['receipt'],MACE_full=full['manifest'],MACE_core=core['manifest']);rows[key]=row
    for name in CASES:
        ends=[rows.get(name+'_'+metal,{}) for metal in ('Ca','La')];values=[v.get('actual_energy_change_kcal_mol') for v in ends]
        delta=None if None in values else values[0]-values[1];valid=delta if all(v.get('pass') for v in ends) else None;base=r['cases'][name]['hybrid_R_kcal_scale']
        scores[name]={'static_R_kcal_scale':base,'actual_delta_R_kcal_scale':delta,'actual_R_kcal_scale':None if delta is None else base+delta,
                      'endpoint_validated_delta_R_kcal_scale':valid,'qualified_R_kcal_scale':None,'calibrated_class':None}
    x,y=[scores[name] for name in ('GGR_extended','GGR_connected')]
    raw_response=None if x['actual_delta_R_kcal_scale'] is None or y['actual_delta_R_kcal_scale'] is None else y['actual_delta_R_kcal_scale']-x['actual_delta_R_kcal_scale']
    raw_final=None if x['actual_R_kcal_scale'] is None or y['actual_R_kcal_scale'] is None else y['actual_R_kcal_scale']-x['actual_R_kcal_scale']
    partition={'response_difference_kcal':None,'final_difference_kcal':None,'pass':False,
               'actual_response_difference_kcal':raw_response,'actual_final_difference_kcal':raw_final,
               'actual_response_threshold_pass':None if raw_response is None else abs(raw_response)<=2,
               'actual_final_threshold_pass':None if raw_final is None else abs(raw_final)<=2}
    if all(v['endpoint_validated_delta_R_kcal_scale'] is not None for v in (x,y)):
        d=y['endpoint_validated_delta_R_kcal_scale']-x['endpoint_validated_delta_R_kcal_scale'];final=y['actual_R_kcal_scale']-x['actual_R_kcal_scale']
        partition.update(response_difference_kcal=d,final_difference_kcal=final,pass_=abs(d)<=2 and abs(final)<=2)
        partition['pass']=partition.pop('pass_')
    for row in scores.values():
        if partition['pass'] and row['endpoint_validated_delta_R_kcal_scale'] is not None:row['qualified_R_kcal_scale']=row['actual_R_kcal_scale']
    contrasts=[]
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        for ggr in ('GGR_extended','GGR_connected'):
            x,y=scores[alpha],scores[ggr];original=x['static_R_kcal_scale']-y['static_R_kcal_scale'];raw=None if x['actual_R_kcal_scale'] is None or y['actual_R_kcal_scale'] is None else x['actual_R_kcal_scale']-y['actual_R_kcal_scale'];q=None if x['qualified_R_kcal_scale'] is None or y['qualified_R_kcal_scale'] is None else x['qualified_R_kcal_scale']-y['qualified_R_kcal_scale']
            contrasts.append({'alpha':alpha,'GGR':ggr,'static_margin_kcal':original,'actual_margin_kcal':raw,'actual_direction_pass':None if raw is None else raw>.02,'qualified_margin_kcal':q,'qualified_direction_pass':q is not None and q>.02})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);source=out/'report_source.py';source.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'prepared':record(prepared),'quantum':qc,'MACE':mc,'implementation':record(source),'rows':rows,'scores':scores,'partition':partition,'contrasts':contrasts,'all_four_actual_direction_gate':all(c['actual_direction_pass'] is True for c in contrasts),'all_four_direction_gate':all(c['qualified_direction_pass'] for c in contrasts),'excluded':top['excluded'],'unavailable_predictions':top['unavailable_predictions'],'baseline_changed':False,'all_cases_consumed_development':True,'aqueous_score':None,'entropy_correction':None,'report_wall_seconds':time.monotonic()-start,'report_CPU_seconds':time.process_time()-cpu}
    write_new(out/'result.json',result);return {k:result[k] for k in ('status','scores','partition','contrasts','all_four_actual_direction_gate','all_four_direction_gate')}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for name,keys in [('prepare',('assessment','output')),('validate',('manifest',)),('validate_quantum',('manifest',)),('collect_quantum',('manifest',)),('report',('prepared','quantum','mace','output'))]:
        q=sub.add_parser(name)
        for key in keys:q.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command');print(json.dumps(globals()[command](**args),indent=2))
