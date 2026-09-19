"""Four-state masked-MACE accounting for context-prepared water hydrogens."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms,write_xyz
from mace_file_checks import cached_file_checks
from contextual_water_prepare import completed_result,isolated_dispatch

STAGE='contextual_water_scoring'
PROTOCOL='masked_omol_context_prepared_water_four_state_v1'


def reuse_score(pin,task,model,software):
    r,t,m=completed_result(pin)
    if m['model']!=model or m['software']!=software:raise InvalidArtifact('scorer reuse changes model/software')
    for key in ('charge','spin_multiplicity','state','metal_index','assembly','microstate','energy_component',
                'charge_feature_adapter','edge_adapter','energy_only','capture_native_readout','require_isolated_metal'):
        if t.get(key)!=task.get(key):raise InvalidArtifact('scorer reuse changes '+key)
    if xyz(verify(t['xyz']))!=xyz(verify(task['xyz'])):raise InvalidArtifact('scorer reuse changes physical coordinates')
    return r


@cached_file_checks
def prepare(collection,output):
    import mace_omol as omol
    from contextual_water_prepare import collect as collect_prepared
    from mace_omol_ablation_run import descriptor_model
    from mace_omol_prepared import tasks as original_tasks
    from mace_omol_intact import geometry
    c=read_json(collection)
    if c['status']!='complete' or c!=collect_prepared(verify(c['manifest'])):raise InvalidArtifact('complete actual contextual preparation required')
    source=read_json(verify(c['manifest']));req=read_json(verify(source['request']))
    parent=read_json(verify(req['scorer_source_manifest']))
    if parent['model']!=descriptor_model(verify(parent['software'])):raise InvalidArtifact('qualified masked scorer required')
    _,out,m=omol.common(verify(parent['inventory']),verify(parent['software']),verify(req['agreement']),output,STAGE)
    isolated_dispatch(out,m);m['model']=parent['model'];all_tasks=[];reused={}
    definitions={r['case_id']:r for r in req['cases']}
    for case in c['cases']:
        if case['operation']=='exact_zero_water_identity':continue
        for t in original_tasks(verify(case['source_preparation'])):
            metal=t['metal'];position=t['position'];tid=case['case_id']+'__'+metal+'_'+position
            t.update(task_id=tid,case_id=case['case_id'],source_xyz=case['endpoints'][metal]['xyz'],
                     contextual_water_collection=record(collection),contextual_case_id=case['case_id'])
            path=out/(tid+'.xyz');write_xyz(path,geometry(t));t['xyz']=record(path)
            t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
            pin=definitions[case['case_id']].get('reuse_scores',{}).get(metal+'_'+position)
            if pin:reuse_score(pin,t,m['model'],m['software']);reused[tid]=pin
            all_tasks.append(t)
    m.update(protocol_id=PROTOCOL,prepared_collection=record(collection),request=source['request'],
        scorer_tasks=all_tasks,reused=reused,tasks=[t for t in all_tasks if t['task_id'] not in reused],
        zero_water_cases=[r for r in c['cases'] if r['operation']=='exact_zero_water_identity'],
        baseline_changed=False,calibrated_decision=None,occupancy_probabilities=None)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


@cached_file_checks
def validate(manifest):
    from mace_omol_ablation_run import descriptor_model
    from mace_omol_intact import geometry
    m=read_json(manifest);c=read_json(verify(m['prepared_collection']));req=read_json(verify(m['request']))
    if m['stage']!=STAGE or m['protocol_id']!=PROTOCOL or m['model']!=descriptor_model(verify(m['software'])):
        raise InvalidArtifact('contextual scorer model changed')
    for pin in [m['agreement'],req['scorer_source_manifest'],*m['implementation'].values()]:verify(pin)
    cases={r['case_id']:r for r in c['cases']};expected={name+'__'+z+'_'+p for name,r in cases.items()
        if r['operation']!='exact_zero_water_identity' for z in ('Ca','La') for p in ('bound','detached')}
    if (len(m['scorer_tasks'])!=len(expected) or {t['task_id'] for t in m['scorer_tasks']}!=expected
            or len(m['tasks'])!=len(expected-set(m['reused'])) or {t['task_id'] for t in m['tasks']}!=expected-set(m['reused'])):
        raise InvalidArtifact('four-state coverage changed')
    for t in m['scorer_tasks']:
        case=cases[t['contextual_case_id']]
        if (t['source_xyz']!=case['endpoints'][t['metal']]['xyz'] or xyz(verify(t['xyz']))!=geometry(t)
                or t['charge']!=case['endpoints'][t['metal']]['charge'] or t['state']!=check_atoms(xyz(verify(t['source_xyz'])),t['charge'])):
            raise InvalidArtifact('prepared scoring coordinates/state changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scoring cache identity changed')
        if t['task_id'] in m['reused']:reuse_score(m['reused'][t['task_id']],t,m['model'],m['software'])
    return {'status':'pass','tasks':len(m['tasks']),'reused_endpoints':len(m['reused']),
        'zero_water_identities':len(m['zero_water_cases']),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest);m=read_json(manifest);mp=Path(manifest);rows={};cases=[]
    source=read_json(verify(m['prepared_collection']));req=read_json(verify(m['request']));definition={r['case_id']:r for r in req['cases']}
    for t in m['scorer_tasks']:
        pin=m['reused'].get(t['task_id']);r=None
        if pin:r=reuse_score(pin,t,m['model'],m['software'])
        else:
            for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
                candidate=accepted_attempt(a,t,mp)
                if candidate is not None:r=candidate;pin=record(a/'result.json')
        rows[t['task_id']]={'status':'computed','energy_eV':r['energy_eV'],'result':pin,'reused':t['task_id'] in m['reused']} if r else {'status':'unavailable'}
    for case in source['cases']:
        entry={'case_id':case['case_id'],'operation':case['operation'],'interaction_R_model_kcal':None,
            'calibrated_decision':None,'occupancy_probability':None,'source_preparation':case['source_preparation']}
        if case['operation']=='exact_zero_water_identity':
            entry.update(status='exact_identity_use_original_scorer_or_cache',endpoints=case['endpoints']);cases.append(entry);continue
        tasks=[t for t in m['scorer_tasks'] if t['case_id']==case['case_id']]
        if len(tasks)!=4 or any(rows[t['task_id']]['status']!='computed' for t in tasks):
            entry.update(status='unavailable',failure='four explicit bound/detached endpoints required');cases.append(entry);continue
        energies={(t['metal'],t['position']):rows[t['task_id']]['energy_eV'] for t in tasks}
        contrast=((energies['Ca','bound']-energies['Ca','detached'])-(energies['La','bound']-energies['La','detached']))*EV_TO_KCAL
        entry.update(status='computed',interaction_R_model_kcal=contrast,endpoint_energies_eV={z+'_'+p:e for (z,p),e in energies.items()},
            common_paired_geometry=case['common_paired_geometry'],bound_total_R_model_kcal=None,
            detached_environment_difference_model_kcal=None,component_closure_error_model_kcal=None)
        disconnected=definition[case['case_id']].get('disconnected_reference')
        if disconnected:
            from hydration_scanner import contrasts
            reference=read_json(verify(disconnected['inventory']));row=next(r for r in reference['rows'] if r['case_id']==disconnected['case_id'])
            entry.update(contrasts(energies,row['disconnected']));entry['disconnected_reference']=disconnected
        cases.append(entry)
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete','manifest':record(manifest),
        'rows':rows,'cases':cases,'baseline_changed':False,'occupancy_probabilities':None,'calibrated_decision':None}


def report(collection,output):
    c=read_json(collection)
    if c!=collect(verify(c['manifest'])):raise InvalidArtifact('scored collection changed')
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',c)
    lines=['# Context-prepared masked MACE components','',
        'Uncalibrated interaction descriptors; zero-water inputs retain the original scorer/cache.','',
        '| Case | Status | Ca−La interaction model-kcal | Detached environment term |','|---|---|---:|---:|']
    for r in c['cases']:lines.append(f"| {r['case_id']} | {r['status']} | {r['interaction_R_model_kcal']} | {r.get('detached_environment_difference_model_kcal')} |")
    (out/'TABLES.md').write_text('\n'.join(lines)+'\n');return {'status':c['status'],'cases':len(c['cases'])}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in [('prepare',('collection','output')),('dry-run',('manifest',)),('execute',('manifest',)),('collect',('manifest','output')),('report',('collection','output'))]:
        q=sub.add_parser(op)
        for name in fields:q.add_argument('--'+name,required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='dry-run':r=validate(**a)
    elif op=='execute':
        validate(a['manifest']);m=read_json(a['manifest'])
        if m['tasks']:
            run=subprocess.run([sys.executable,str(verify(m['implementation']['mace_hybrid.py'])),'execute','--manifest',str(Path(a['manifest']).resolve()),'--memory-mode','native'],check=True)
            r={'status':'executor_returned_collect_required','returncode':run.returncode}
        else:r={'status':'all_scoring_reused_or_zero_water_identity','new_model_calls':0}
    elif op=='collect':
        dest=a.pop('output');r=collect(**a);write_new(dest,r);r={'status':r['status'],'cases':len(r['cases'])}
    else:r=globals()[op](**a)
    print(json.dumps(r,indent=2))
