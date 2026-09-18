"""Fixed native checks and all-structure report for matched GGR transfer."""
from __future__ import annotations
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_bounded_response import ball,endpoint_checks
from mace_omol_hybrid_transfer import POLICY,STAGE,STRUCTURES,REPRESENTATIONS,prepared_states,quantum_manifest
from mace_omol_hybrid_response import task as response_task
from mace_omol_gradient_run import CONFIG,check_parent,model
from mace_omol_vacuum import METHOD,scientific_input


def inputs(assessment):
    a=read_json(assessment);p,cases=prepared_states(verify(a['preparation']))
    if a['status']!='complete' or a['protocol_id']!=POLICY or set(a['rows'])!=set(cases):raise InvalidArtifact('transfer assessment incomplete/different')
    for c in a['MACE_collections']:
        verify(c['manifest'])
        if c['status']!='complete' or not c['numerical_checks_pass']:raise InvalidArtifact('failed initial learned collection')
    verify(a['quantum']['manifest']);verify(a['implementation']);selected={}
    for name,row in a['rows'].items():
        for metal,s in row.items():
            pred=s['prediction'];k={key:np.array(value) for key,value in s['matrices_kcal_mol_A2'].items()};g=np.array(s['gradient_kcal_mol_A'])
            if not np.allclose(g,np.array(s['DFT_gradient_kcal_mol_A'])+s['J_gradient_kcal_mol_A'],atol=1e-10,rtol=0):raise InvalidArtifact('hybrid gradient accounting differs')
            fine=ball(k['fine'],g);coarse=ball(k['coarse'],g)
            if pred['status']!='eligible_for_native_check':continue
            if fine['status']!='solved' or coarse['status']!='solved':raise InvalidArtifact('unstable prediction marked eligible')
            u=np.array(fine['displacement_A']);err=float(.5*u@(k['fine']-k['coarse'])@u);delta=np.linalg.norm(u-coarse['displacement_A'])
            if abs(err)>.02 or delta>.01 or not all(c['pass'] for c in a['odd_checks'][cases[name]['global_id']+'_'+metal]):raise InvalidArtifact('numerically unsupported prediction marked eligible')
            if not np.allclose(u,pred['displacement_A'],atol=1e-10,rtol=0) or abs(fine['predicted_change_kcal_mol']-pred['predicted_change_kcal_mol'])>1e-10:raise InvalidArtifact('frozen predicted position/energy differs')
            selected[name+'_'+metal]={'case_id':name,'metal':metal,'prediction':pred}
    return a,p,cases,selected


@cached_file_checks
def prepare(assessment,output):
    from mace_omol import common,seal
    from mace_mechanics import physical_model
    from ggr_sensitivity import membership
    start=time.monotonic();cpu=time.process_time();a,p,cases,selected=inputs(assessment)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    _,md,m=common(verify(p['inventory']),verify(p['software']),verify(p['agreement']),out/'mace',STAGE)
    states={};excluded={};tasks=[]
    for key,item in selected.items():
        name,metal=item['case_id'],item['metal'];c=cases[name];repair=read_json(verify(c['source_preparation']))
        graph,coords,pos,mk,_,_=physical_model(repair);u=np.array(item['prediction']['displacement_A'])
        before=membership(graph,coords,pos,mk);after=membership(graph,coords,pos+u,mk)
        if before!=after:
            excluded[key]={'status':'donor_membership_changed','before':before,'after':after};continue
        ends={}
        for kind,e in [('core',c['endpoints'][metal]),('full',p['fulls'][c['global_id']][metal])]:
            t=response_task(e,name,metal,kind,'bounded',u.tolist(),md);t['energy_only']=False;tasks.append(t);ends[kind]={k:t[k] for k in ('xyz','charge','spin_multiplicity','state','metal_index')}
        states[key]={**item,'endpoints':ends,'donor_membership':before}
    top={'protocol_id':POLICY,'base_preparation':a['preparation'],'assessment':record(assessment),'states':states,'excluded':excluded,
         'unavailable_predictions':{n+'_'+metal:s['prediction'] for n,row in a['rows'].items() for metal,s in row.items() if n+'_'+metal not in selected},
         'new_DFT_calls':len(states),'new_MACE_calls':len(tasks),'preparation_wall_seconds':time.monotonic()-start,'preparation_CPU_seconds':time.process_time()-cpu}
    write_new(out/'preparation.json',top)
    if not states:return {'status':'no_supported_native_points','preparation':record(out/'preparation.json')}
    m.update(protocol_id=POLICY,phase='native',group='native',preparation=a['preparation'],native_preparation=record(out/'preparation.json'),core_report=p['core_report'],gradient_settings=CONFIG,model=model(verify(p['software'])),tasks=tasks)
    check_parent(m);seal(md,m)
    qs={key:{'case_id':s['case_id'],'metal':s['metal'],'core':s['endpoints']['core']} for key,s in states.items()}
    q=quantum_manifest(p,verify(a['preparation']),qs,out/'quantum','native',record(out/'preparation.json'))
    return {'status':'prepared','preparation':record(out/'preparation.json'),'quantum':q,'MACE':record(md/'manifest.json'),'excluded':excluded}


def expected_native(m):
    top=read_json(verify(m['native_preparation']));a,p,cases,selected=inputs(verify(top['assessment']))
    if m['preparation']!=top['base_preparation'] or top['base_preparation']!=a['preparation'] or set(top['states'])|set(top['excluded'])!=set(selected):raise InvalidArtifact('native preparation/eligible inventory differs')
    expected={}
    for key,s in top['states'].items():
        item=selected[key];name,metal=item['case_id'],item['metal'];u=item['prediction']['displacement_A'];c=cases[name]
        if s['case_id']!=name or s['metal']!=metal or s['prediction']!=item['prediction']:raise InvalidArtifact('native point differs from fixed prediction')
        for kind,e in [('core',c['endpoints'][metal]),('full',p['fulls'][c['global_id']][metal])]:
            old=xyz(verify(e['xyz']));new=xyz(verify(s['endpoints'][kind]['xyz']));want=np.array([x[1:] for x in old]);want[e['metal_index']]+=u
            if [x[0] for x in old]!=[x[0] for x in new] or not np.allclose([x[1:] for x in new],want,atol=5e-10,rtol=0) or s['endpoints'][kind]['charge']!=e['charge']:raise InvalidArtifact('fixed native coordinates/charge differ')
            expected[kind+'__'+key+'__bounded']=(e,name,metal,kind,'bounded',u)
    return p,expected


@cached_file_checks
def validate_quantum(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);p,expected=expected_native(m);top=read_json(verify(m['native_preparation']))
    if m['phase']!='native' or m['stage']!=STAGE or m['protocol_id']!=POLICY or m['method']!=METHOD or len(m['tasks'])>8 or len(m['tasks'])!=len(top['states']) or {t['task_id'] for t in m['tasks']}!=set(top['states']):raise InvalidArtifact('native quantum method/inventory differs')
    for pin in [m['agreement'],m['orca'],*m['implementation'].values()]:verify(pin)
    for t in m['tasks']:
        state=top['states'][t['task_id']];e=state['endpoints']['core'];payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['case_id']!=state['case_id'] or t['metal']!=state['metal'] or t['source_mapping']!=p['cases'][t['case_id']] or t['xyz']['sha256']!=e['xyz']['sha256'] or t['charge']!=e['charge'] or t['multiplicity']!=1 or verify(t['input']).read_text()!=scientific_input(e['charge']):raise InvalidArtifact('native quantum source/input differs')
        if t['cache_key']!=cache_key({'task':payload,'method':METHOD,'protocol':POLICY,'implementation':m['implementation']}):raise InvalidArtifact('native quantum cache differs')
    return dry_run(manifest)


@cached_file_checks
def report(prepared,quantum,mace,output):
    from mace_omol_hybrid_transfer import collect,collect_quantum
    start=time.monotonic();cpu=time.process_time();top=read_json(prepared);a,p,cases,selected=inputs(verify(top['assessment']))
    if read_json(quantum)['native_preparation']!=record(prepared) or read_json(mace)['native_preparation']!=record(prepared):raise InvalidArtifact('native manifests/preparation differ')
    q=collect_quantum(quantum);m=collect(mace)
    if q['status']!='complete' or m['status']!='complete' or not m['numerical_checks_pass']:raise InvalidArtifact('native actual results incomplete/invalid')
    rows={};scores={};partition={};reference=read_json(verify(p['reference_report']))
    for key,s in top['states'].items():
        name,metal=s['case_id'],s['metal'];initial=a['static_scores'][name]['endpoints'][metal];d=q['rows'][key];core=m['rows']['core__'+key+'__bounded'];full=m['rows']['full__'+key+'__bounded'];idx=s['endpoints']['full']['metal_index']
        dd=(d['energy_hartree']-initial['DFT_hartree'])*HA_TO_KCAL;dj=((full['energy_eV']-initial['full_eV'])-(core['energy_eV']-initial['core_eV']))*EV_TO_KCAL
        g=np.array(d['gradient_kcal_mol_per_A'][0])+(np.load(verify(full['gradient']))[idx]-np.load(verify(core['gradient']))[0])*EV_TO_KCAL
        rows[key]=endpoint_checks(s['prediction'],a['rows'][name][metal],dd,dj,g)
    for name,c in cases.items():
        ends=[rows.get(name+'_'+metal,{}) for metal in ('Ca','La')];values=[v.get('actual_energy_change_kcal_mol') for v in ends];delta=None if None in values else values[0]-values[1];base=a['static_scores'][name]['static_R_kcal_scale']
        scores[name]={'static_R_kcal_scale':base,'actual_delta_R_kcal_scale':delta,'actual_R_kcal_scale':None if delta is None else base+delta,'endpoint_validated_delta_R_kcal_scale':delta if all(v.get('pass') for v in ends) else None,'qualified_R_kcal_scale':None,'calibrated_class':None}
    for gid in STRUCTURES:
        x,y=[scores[gid+'_'+r] for r in ('extended','connected')];available=all(x.get(k) is not None and y.get(k) is not None for k in ('actual_R_kcal_scale','actual_delta_R_kcal_scale'))
        dr=y['actual_delta_R_kcal_scale']-x['actual_delta_R_kcal_scale'] if available else None;df=y['actual_R_kcal_scale']-x['actual_R_kcal_scale'] if available else None
        ok=available and abs(dr)<=2 and abs(df)<=2 and all(z['endpoint_validated_delta_R_kcal_scale'] is not None for z in (x,y))
        partition[gid]={'actual_response_difference_kcal':dr,'actual_final_difference_kcal':df,'actual_response_threshold_pass':None if dr is None else abs(dr)<=2,'actual_final_threshold_pass':None if df is None else abs(df)<=2,'pass':ok}
        if ok:
            for z in (x,y):z['qualified_R_kcal_scale']=z['actual_R_kcal_scale']
    combined=copy.deepcopy(scores)
    for rep in REPRESENTATIONS:combined['GGR_1GLG_'+rep]=reference['scores']['GGR_'+rep]
    contrasts=[]
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        x=reference['scores'][alpha]
        for gid in ('GGR_1GLG',*STRUCTURES):
            for rep in REPRESENTATIONS:
                name=gid+'_'+rep;y=combined[name];margin=None if y['actual_R_kcal_scale'] is None else x['actual_R_kcal_scale']-y['actual_R_kcal_scale'];qualified=None if x['qualified_R_kcal_scale'] is None or y['qualified_R_kcal_scale'] is None else x['qualified_R_kcal_scale']-y['qualified_R_kcal_scale']
                contrasts.append({'alpha':alpha,'GGR':name,'reused_comparison':gid=='GGR_1GLG','static_margin_kcal':x['static_R_kcal_scale']-y['static_R_kcal_scale'],'actual_margin_kcal':margin,'actual_direction_pass':None if margin is None else margin>.02,'qualified_margin_kcal':qualified,'qualified_direction_pass':qualified is not None and qualified>.02})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'report_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'prepared':record(prepared),'reference_report':p['reference_report'],'quantum':q,'MACE':m,'rows':rows,'scores':scores,'partition':partition,'contrasts':contrasts,'implementation':record(sp),'all_twelve_actual_directions_pass':all(x['actual_direction_pass'] is True for x in contrasts),'all_twelve_qualified_directions_pass':all(x['qualified_direction_pass'] for x in contrasts),'biological_groups':2,'all_cases_consumed_development':True,'excluded':top['excluded'],'unavailable_predictions':top['unavailable_predictions'],'baseline_changed':False,'aqueous_score':None,'entropy_correction':None,'report_wall_seconds':time.monotonic()-start,'report_CPU_seconds':time.process_time()-cpu}
    write_new(out/'result.json',result);return {k:result[k] for k in ('status','partition','contrasts','all_twelve_actual_directions_pass','all_twelve_qualified_directions_pass')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for op,keys in [('prepare',('assessment','output')),('report',('prepared','quantum','mace','output'))]:
        q=sub.add_parser(op)
        for k in keys:q.add_argument('--'+k,required=True)
    a=vars(p.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
