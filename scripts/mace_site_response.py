"""Project real archived hybrid gradients onto source-defined donor coordinates."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,check_atoms
from mace_site_coordinates import SiteCoordinates,POLICY


@cached_file_checks
def sources(config):
    from mace_omol_hybrid_response import sources as old_sources,collect as old_collect
    from mace_omol_hybrid_transfer import prepared_states,collect,collect_quantum
    cfg=read_json(config);initial=read_json(verify(cfg['matched_initial']))
    p,r,gr,gm,states,fulls,reused=old_sources(verify(initial['source_preparation']),verify(initial['static_report']),verify(initial['ggr_gradients']),True)
    refs=[ref for ref in initial['jobs'] if read_json(verify(ref))['group']=='centers']
    if len(refs)!=1:raise InvalidArtifact('ambiguous original gradient centers')
    center=old_collect(verify(refs[0]))
    if center['status']!='complete' or not center['numerical_checks_pass']:raise InvalidArtifact('original gradient centers incomplete')
    result={};receipts=[center]
    for name,ends in states.items():
        c=read_json(verify(p['cases'][name]));orig=read_json(verify(c['source_mapping']));label='GGR_1GLG_'+name.split('_',1)[1] if name.startswith('GGR_') else name
        for metal,s in ends.items():
            full=(reused[s['global_id']+'_'+metal]['result'] if s['global_id']=='GGR_1GLG' else center['rows']['full__'+s['global_id']+'_'+metal+'__center'])
            core=center['rows']['core__'+name+'_'+metal+'__center']
            result[label+'_'+metal]={'case_id':label,'global_id':s['global_id'],'metal':metal,'repair':orig['source_preparation'],
                'whole_preparation':c['normalized_global_preparation'],'core_state':s['core'],'full_state':s['full'],
                'DFT':s['DFT_source'],'MACE_core':core,'MACE_full':full}
    a=read_json(verify(cfg['transfer_assessment']));tp,tc=prepared_states(verify(a['preparation']))
    actual_q=collect_quantum(verify(a['quantum']['manifest']))
    if actual_q['status']!='complete' or actual_q!=a['quantum']:raise InvalidArtifact('transfer native centers differ from receipts')
    refs=[c for c in a['MACE_collections'] if read_json(verify(c['manifest']))['group']=='centers']
    if len(refs)!=1:raise InvalidArtifact('ambiguous transfer gradient centers')
    cm=collect(verify(refs[0]['manifest']))
    if cm!=refs[0] or cm['status']!='complete' or not cm['numerical_checks_pass']:raise InvalidArtifact('transfer gradient centers differ from receipts')
    receipts.append(cm)
    for name,c in tc.items():
        for metal,e in c['endpoints'].items():
            key=name+'_'+metal
            result[key]={'case_id':name,'global_id':c['global_id'],'metal':metal,'repair':c['source_preparation'],'whole_preparation':c['whole_preparation'],
                'core_state':e,'full_state':tp['fulls'][c['global_id']][metal],'DFT':a['quantum']['rows'][key],
                'MACE_core':cm['rows']['core__'+key+'__center'],'MACE_full':cm['rows']['full__'+c['global_id']+'_'+metal+'__center']}
    if len(result)!=16:raise InvalidArtifact('fixed 16-state inventory differs')
    return result,receipts


@cached_file_checks
def prepare(config,agreement,output):
    start=time.monotonic();cpu=time.process_time();states,collections=sources(config)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);rows={};shared={};checks=[]
    for key,s in states.items():
        repair=read_json(verify(s['repair']));whole=read_json(verify(s['whole_preparation']));core=xyz(verify(s['core_state']['xyz']));full=xyz(verify(s['full_state']['xyz']))
        geom=SiteCoordinates(repair,whole,verify(s['core_state']['xyz']));check=geom.checks()
        if not check['pass']:raise InvalidArtifact('physical-coordinate check failed: '+key+' '+str(check))
        expected=16 if s['global_id'].startswith('GGR_') else 11
        if len(geom.modes)!=expected:raise InvalidArtifact('declared physical coordinate count differs')
        for kind,coords in [('core',core),('full',full)]:
            e=s[kind+'_state']
            if check_atoms(coords,e['charge'])!=e['state']:raise InvalidArtifact('actual paired electronic state differs')
        if not np.array_equal([a[1:] for a in full],geom.positions):raise InvalidArtifact('whole gradient physical coordinates differ')
        gd=np.array(s['DFT']['gradient_kcal_mol_per_A']);gc=np.load(verify(s['MACE_core']['gradient']))*EV_TO_KCAL;gf=np.load(verify(s['MACE_full']['gradient']))*EV_TO_KCAL
        if gd.shape!=gc.shape or gd.shape!=(len(core),3) or gf.shape!=(len(full),3) or not all(np.isfinite(g).all() for g in (gd,gc,gf)):raise InvalidArtifact('actual gradient shape/order differs')
        for r in (s['MACE_core'],s['MACE_full']):
            if r['gradient_unit']!='eV_equivalent_per_A' or r['gradient_definition']!='Cartesian_gradient_of_masked_descriptor':raise InvalidArtifact('MACE derivative sign/unit differs')
        projections={name:np.einsum('mij,ij->m',jac,g) for name,jac,g in [('DFT',geom.core_jacobian,gd),('MACE_core',geom.core_jacobian,gc),('MACE_full',geom.full_jacobian,gf)]}
        gh=projections['DFT']+projections['MACE_full']-projections['MACE_core']
        physical=geom.source_gradient(gd-gc)+gf
        replay=np.einsum('mij,ij->m',geom.full_jacobian,physical)
        dual_error=float(np.max(np.abs(gh-replay)))
        if dual_error>1e-9:raise InvalidArtifact('source/cap chain-rule projections disagree')
        group=s['global_id']
        if group in shared:
            old=shared[group]
            if old['modes']!=geom.modes or not np.array_equal(np.load(verify(old['full_jacobian'])),geom.full_jacobian):raise InvalidArtifact('paired/partition physical coordinate measure differs')
            fj=old['full_jacobian']
        else:
            p=out/(group+'_full_J.npy');np.save(p,geom.full_jacobian);fj=record(p);shared[group]={'modes':geom.modes,'full_jacobian':fj,'physical_ids':geom.ids,'whole_preparation':s['whole_preparation']}
        directory=out/key;directory.mkdir();p=directory/'core_J.npy';np.save(p,geom.core_jacobian);cp=record(p)
        p=directory/'hybrid_physical_gradient_kcal_A.npy';np.save(p,physical)
        row={'source':s,'modes':geom.modes,'full_jacobian':fj,'core_jacobian':cp,'hybrid_physical_gradient':record(p),
             'components':{k:v.tolist() for k,v in projections.items()},'hybrid_projected_gradient':gh.tolist(),
             'linearized_probe_change_kcal_equivalent':(gh*np.array([m['linear_probe_amplitude'] for m in geom.modes])).tolist(),
             'coordinate_checks':check,'source_projection_error':dual_error,'response_model_status':'response_model_not_validated','score_correction':None}
        write_new(directory/'projection.json',row);rows[key]=record(directory/'projection.json');checks.append({'state':key,'pass':True,'coordinate_check':check,'projection_error':dual_error})
    source=out/'coordinate_source.py';source.write_bytes((Path(__file__).parent/'mace_site_coordinates.py').read_bytes());driver=out/'analysis_source.py';driver.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete_gradient_projection','policy_id':POLICY,'config':record(config),'agreement':record(agreement),'rows':rows,'physical_systems':shared,'checks':checks,
            'implementation':{'coordinates':record(source),'driver':record(driver)},'new_DFT_calls':0,'new_MACE_calls':0,'actual_archived_electronic_states':16,
            'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'baseline_changed':False,'response_model_status':'response_model_not_validated','score_correction':None}
    write_new(out/'preparation.json',result);return {'status':result['status'],'preparation':record(out/'preparation.json'),'states':len(rows),'physical_systems':len(shared)}


@cached_file_checks
def report(prepared,output):
    start=time.monotonic();p=read_json(prepared)
    if p['status']!='complete_gradient_projection' or p['policy_id']!=POLICY or len(p['rows'])!=16 or not all(c['pass'] for c in p['checks']):raise InvalidArtifact('incomplete/failed gradient projections')
    pairs={};partitions={}
    for name in sorted({key.rsplit('_',1)[0] for key in p['rows']}):
        ca,la=[read_json(verify(p['rows'][name+'_'+metal])) for metal in ('Ca','La')]
        if ca['modes']!=la['modes'] or ca['full_jacobian']!=la['full_jacobian']:raise InvalidArtifact('paired physical coordinate measures differ')
        modes=ca['modes'];delta=np.array(ca['hybrid_projected_gradient'])-la['hybrid_projected_gradient'];probes=delta*np.array([m['linear_probe_amplitude'] for m in modes])
        pairs[name]={'modes':modes,'projected_gradient_R':delta.tolist(),'linearized_R_probe_kcal_equivalent':probes.tolist(),
                     'metal_max_abs_linearized_R_probe':float(np.max(np.abs(probes[:3]))),'donor_max_abs_linearized_R_probe':float(np.max(np.abs(probes[3:]))),
                     'Ca':{'gradient':ca['hybrid_projected_gradient'],'probe_change':ca['linearized_probe_change_kcal_equivalent']},'La':{'gradient':la['hybrid_projected_gradient'],'probe_change':la['linearized_probe_change_kcal_equivalent']},
                     'is_energy_evaluation':False,'score_correction':None}
    for gid in ('GGR_1GLG','GGR_2FW0','GGR_2FVY'):
        a,b=[pairs[gid+'_'+rep] for rep in ('extended','connected')]
        if a['modes']!=b['modes']:raise InvalidArtifact('GGR representations have different physical coordinates')
        d=np.array(b['linearized_R_probe_kcal_equivalent'])-a['linearized_R_probe_kcal_equivalent'];partitions[gid]={'linearized_probe_difference':d.tolist(),'maximum_abs':float(np.max(np.abs(d))),'not_an_energy_partition_test':True}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'report_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete_gradient_diagnostic','prepared':record(prepared),'pairs':pairs,'partition_derivatives':partitions,'implementation':record(sp),'wall_seconds':time.monotonic()-start,
            'new_DFT_calls':0,'new_MACE_calls':0,'biological_groups':2,'all_cases_consumed_development':True,'response_model_status':'response_model_not_validated','score_correction':None,'baseline_changed':False,'affinity_or_uncertainty_score':None}
    write_new(out/'result.json',result);return {'status':result['status'],'pairs':{n:{k:v[k] for k in ('metal_max_abs_linearized_R_probe','donor_max_abs_linearized_R_probe')} for n,v in pairs.items()},'partition_derivatives':partitions}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for op,keys in [('prepare',('config','agreement','output')),('report',('prepared','output'))]:
        q=sub.add_parser(op)
        for k in keys:q.add_argument('--'+k,required=True)
    a=vars(p.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
