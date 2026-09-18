"""Native validation of the fixed coupled path, with all structural comparisons."""
from __future__ import annotations
import argparse,json,shutil,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_site_path import POLICY,STAGE,prepared,task,collect,select
from mace_site_kinematics import Kinematics
from mace_omol_gradient_run import CONFIG,check_parent,model
from mace_omol_vacuum import METHOD,scientific_input,parse_endpoint


def check_selection(a,p):
    """Recompute the declared choice from real grid results, without new inference."""
    actual={k:r for c in a['collections'] for k,r in c['rows'].items()}
    expected={k for k,s in p['states'].items() if s['status']=='prepared'}
    if set(a['selected'])!=expected or a['excluded']!=p['excluded']:
        raise InvalidArtifact('selected/unsupported state inventory changed')
    for key in expected:
        state=p['states'][key];row=read_json(verify(state['source_projection']));src=row['source']
        kin=Kinematics(read_json(verify(state['kinematics'])))
        dg=np.array(src['DFT']['gradient_kcal_mol_per_A'])-np.load(verify(src['MACE_core']['gradient']))*EV_TO_KCAL
        values={'center':0.}
        for point,v in state['points'].items():
            core=kin.evaluate(v['q'])[1];r=actual['full__'+key+'__'+point]
            values[point]=(r['energy_eV']-src['MACE_full']['energy_eV'])*EV_TO_KCAL+float(np.sum(dg*(core-kin.core)))
            if a['selected'][key]['all_points'][point]['MACE']!=r:
                raise InvalidArtifact('selected point does not retain its actual grid output')
        # Dicts loaded from JSON are key-sorted; explicitly retain smaller-lambda tie order.
        order=('center','p025','p050','p075','p100');chosen=min(order,key=lambda k:values[k])
        saved=a['selected'][key]
        if saved['point']!=chosen or abs(saved['prediction_kcal']-values[chosen])>1e-9:
            raise InvalidArtifact('selected point differs from fixed all-point energy rule')


@cached_file_checks
def prepare(selection,output):
    from mace_omol import common,seal
    from mace_global_benchmark import snapshot
    from ggr_sensitivity import ORCA
    a=read_json(selection);p=prepared(verify(a['preparation']))
    if a['status']!='complete_selection' or a['protocol_id']!=POLICY:raise InvalidArtifact('fixed completed selection required')
    # Selection must replay actual finite initial receipts before any DFT.
    for c in a['collections']:
        if collect(verify(c['manifest']))!=c:raise InvalidArtifact('selection collection differs from actual receipts')
    check_selection(a,p)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    _,md,m=common(verify(p['inventory']),verify(p['software']),verify(p['agreement']),out/'mace',STAGE)
    tasks=[task(k,p['states'][k],s['point'],kind,False) for k,s in a['selected'].items() if s['point']!='center' for kind in ('full','core')]
    m.update(protocol_id=POLICY,phase='native',group='native',preparation=a['preparation'],selection=record(selection),core_report=p['core_report'],gradient_settings=CONFIG,model=model(verify(p['software'])),tasks=tasks)
    check_parent(m);seal(md,m)
    qd=out/'quantum';qd.mkdir();parent=read_json(verify(p['quantum_parent']))
    names=(*parent['implementation'],'mace_site_path.py','mace_site_path_native.py','mace_site_kinematics.py')
    pins=snapshot(qd,names);qt=[]
    for key,s in a['selected'].items():
        if s['point']=='center':continue
        state=p['states'][key];e=state['points'][s['point']]['endpoints']['core'];d=qd/key;d.mkdir();xp=d/'core.xyz';shutil.copyfile(verify(e['xyz']),xp);ip=d/'endpoint.inp';ip.write_text(scientific_input(e['charge']))
        t={'task_id':key,'case_id':state['case_id'],'metal':state['metal'],'charge':e['charge'],'multiplicity':1,'input':record(ip),'xyz':record(xp),
           'source_mapping':state['source_projection'],'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),'task_type':'analytic_gradient'}
        t['cache_key']=cache_key({'task':t,'method':METHOD,'protocol':POLICY,'implementation':pins});qt.append(t)
    qt.sort(key=lambda t:(-len(xyz(verify(t['xyz']))),t['task_id']))
    qm={'protocol_id':POLICY,'stage':STAGE,'phase':'native','preparation':a['preparation'],'selection':record(selection),'agreement':p['agreement'],
        'implementation':pins,'tasks':qt,'orca':record(ORCA),'method':METHOD,'energy_scope':'isolated_vacuum_endpoint',
        'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},'cost_tracking':{'compute_budget':None,'wall_time_limit':None}}
    write_new(qd/'manifest.json',qm);validate_quantum(qd/'manifest.json')
    top={'selection':record(selection),'preparation':a['preparation'],'mace':record(md/'manifest.json'),'quantum':record(qd/'manifest.json'),'new_DFT_calls':len(qt),'new_MACE_calls':len(tasks)}
    write_new(out/'preparation.json',top);return top


@cached_file_checks
def validate_quantum(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);p=prepared(verify(m['preparation']));a=read_json(verify(m['selection']))
    check_selection(a,p)
    if m['protocol_id']!=POLICY or m['stage']!=STAGE or m['phase']!='native' or m['method']!=METHOD or a['preparation']!=m['preparation']:raise InvalidArtifact('coupled native method differs')
    expected={k:s for k,s in a['selected'].items() if s['point']!='center'}
    if len(m['tasks'])!=len(expected) or {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('native finite inventory differs')
    for ref in m['implementation'].values():verify(ref)
    for t in m['tasks']:
        state=p['states'][t['task_id']];point=expected[t['task_id']]['point'];e=state['points'][point]['endpoints']['core'];payload={k:v for k,v in t.items() if k!='cache_key'}
        if (t['xyz']['sha256']!=e['xyz']['sha256'] or t['charge']!=e['charge'] or t['multiplicity']!=1 or t['source_mapping']!=state['source_projection']
            or verify(t['input']).read_text()!=scientific_input(e['charge']) or t['cache_key']!=cache_key({'task':payload,'method':METHOD,'protocol':POLICY,'implementation':m['implementation']})):
            raise InvalidArtifact('native selected geometry/state/cache differs')
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
    return {'manifest':record(manifest),'rows':rows,'status':'complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete'}


@cached_file_checks
def report(preparation,output):
    start=time.monotonic();top=read_json(preparation);a=read_json(verify(top['selection']));p=prepared(verify(top['preparation']));q=collect_quantum(verify(top['quantum']));m=collect(verify(top['mace']))
    if q['status']!='complete' or m['status']!='complete':raise InvalidArtifact('native execution incomplete; collect missing/failed tasks before reporting')
    rows={};static={};scores={};parts={};contrasts=[]
    for key,state in p['states'].items():
        row=read_json(verify(state['source_projection']));s=row['source']
        static[key]=s['DFT']['energy_hartree']*HA_TO_KCAL+(s['MACE_full']['energy_eV']-s['MACE_core']['energy_eV'])*EV_TO_KCAL
        if key not in a['selected']:
            rows[key]={'status':'unsupported','reason':p['excluded'][key],'actual_change_kcal':None,'pass':False};continue
        selected=a['selected'][key];point=selected['point']
        if point=='center':
            rows[key]={'status':'reused_exact_center','actual_change_kcal':0.,'DFT_change_kcal':0.,'full_change_kcal':0.,'core_change_kcal':0.,'pass':True,'checks':{'no_proposed_displacement':True}};continue
        d=q['rows'][key];f,c=[m['rows'][kind+'__'+key+'__'+point] for kind in ('full','core')];kin=Kinematics(read_json(verify(state['kinematics'])));pos,core,j,cj=kin.evaluate(state['points'][point]['q'])
        dx=core-kin.core;g0=np.array(s['DFT']['gradient_kcal_mol_per_A']);dc=(d['energy_hartree']-s['DFT']['energy_hartree'])*HA_TO_KCAL
        fc=(f['energy_eV']-s['MACE_full']['energy_eV'])*EV_TO_KCAL;cc=(c['energy_eV']-s['MACE_core']['energy_eV'])*EV_TO_KCAL;actual=dc+fc-cc
        dg0=g0-np.load(verify(s['MACE_core']['gradient']))*EV_TO_KCAL
        actual_g=np.einsum('mij,ij->m',cj,np.array(d['gradient_kcal_mol_per_A'])-np.load(verify(c['gradient']))*EV_TO_KCAL)+np.einsum('mij,ij->m',j,np.load(verify(f['gradient']))*EV_TO_KCAL)
        pred_g=np.einsum('mij,ij->m',cj,dg0)+np.einsum('mij,ij->m',j,np.load(verify(f['gradient']))*EV_TO_KCAL)
        scales=np.array([mode['linear_probe_amplitude'] for mode in kin.modes]);gerr=(pred_g-actual_g)*scales;gtol=np.maximum(.04,.25*np.abs(actual_g*scales))
        nonlinear=dc-float(np.sum(g0*dx));etol=max(.05,.25*abs(nonlinear));error=selected['prediction_kcal']-actual
        replay=(f['energy_eV']-selected['all_points'][point]['MACE']['energy_eV'])*EV_TO_KCAL
        checks={'energy_prediction':abs(error)<=etol,'actual_descent':actual<0,'all_mode_gradient_prediction':bool(np.all(np.abs(gerr)<=gtol)),'grid_full_replay':abs(replay)<=.01}
        rows[key]={'status':'computed','point':point,'prediction_kcal':selected['prediction_kcal'],'actual_change_kcal':actual,'DFT_change_kcal':dc,'full_change_kcal':fc,'core_change_kcal':cc,
                   'DFT_nonlinear_change_kcal':nonlinear,'prediction_error_kcal':error,'energy_tolerance_kcal':etol,'projected_actual_gradient':actual_g.tolist(),'projected_predicted_gradient':pred_g.tolist(),
                   'probe_gradient_errors_kcal':gerr.tolist(),'probe_gradient_tolerances_kcal':gtol.tolist(),'full_replay_error_kcal':replay,'checks':{k:bool(v) for k,v in checks.items()},'pass':all(checks.values()),'DFT_receipt':d['receipt']}
    for name in sorted({s['case_id'] for s in p['states'].values()}):
        ca,la=[rows[name+'_'+met] for met in ('Ca','La')];center=static[name+'_Ca']-static[name+'_La'];delta=None if any(s['actual_change_kcal'] is None for s in (ca,la)) else ca['actual_change_kcal']-la['actual_change_kcal']
        scores[name]={'static_R_kcal':center,'response_delta_kcal':delta,'actual_R_kcal':None if delta is None else center+delta,'endpoints_pass':ca['pass'] and la['pass'],'qualified_R_kcal':None}
    for gid in ('GGR_1GLG','GGR_2FW0','GGR_2FVY'):
        x,y=[scores[gid+'_'+rep] for rep in ('connected','extended')];final=None if x['actual_R_kcal'] is None or y['actual_R_kcal'] is None else x['actual_R_kcal']-y['actual_R_kcal'];response=None if final is None else x['response_delta_kcal']-y['response_delta_kcal']
        parts[gid]={'final_difference_kcal':final,'response_difference_kcal':response,'pass':final is not None and abs(final)<=2 and abs(response)<=2}
    for name,s in scores.items():
        gate=s['endpoints_pass'] and (not name.startswith('GGR_') or parts[name.rsplit('_',1)[0]]['pass'])
        if gate:s['qualified_R_kcal']=s['actual_R_kcal']
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        for gid in ('GGR_1GLG','GGR_2FW0','GGR_2FVY'):
            for rep in ('extended','connected'):
                name=gid+'_'+rep;x,y=scores[alpha],scores[name];v=None if x['actual_R_kcal'] is None or y['actual_R_kcal'] is None else x['actual_R_kcal']-y['actual_R_kcal']
                contrasts.append({'alpha':alpha,'GGR':name,'static_margin_kcal':x['static_R_kcal']-y['static_R_kcal'],'actual_margin_kcal':v,'raw_pass':v is not None and v>.02,
                                  'qualified_pass':v is not None and v>.02 and x['qualified_R_kcal'] is not None and y['qualified_R_kcal'] is not None})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'report_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'native_preparation':record(preparation),'selection':top['selection'],'quantum':q,'MACE':m,'rows':rows,'scores':scores,'partition':parts,'contrasts':contrasts,
            'raw_pass_count':sum(v['raw_pass'] for v in contrasts),'qualified_pass_count':sum(v['qualified_pass'] for v in contrasts),'comparison_denominator':12,'biological_groups':2,'all_consumed_development':True,
            'baseline_changed':False,'calibrated_decision':None,'relaxation_free_energy':None,'entropy_correction':None,'is_stationary_minimum':False,'implementation':record(sp),'wall_seconds':time.monotonic()-start}
    write_new(out/'result.json',result);return {k:result[k] for k in ('status','raw_pass_count','qualified_pass_count','comparison_denominator','partition','contrasts')}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for op,keys in [('prepare',('selection','output')),('validate_quantum',('manifest',)),('collect_quantum',('manifest',)),('report',('preparation','output'))]:
        a=sub.add_parser(op)
        for key in keys:a.add_argument('--'+key,required=True)
    a=vars(parser.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
