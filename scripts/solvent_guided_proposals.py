"""Archived ranking replay and contained solvent-guided proposal development."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import json
from pathlib import Path
import statistics
import numpy as np

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from accommodation_nonlinear import relative_components
from accommodation_proposals import score
from nikasha_pool import SCALED_PROTOCOL
from nikasha_pool_compare import load_results
import adaptive_angular_proposals as angular
from mace_site_kinematics import Kinematics
from mace_hybrid import write_xyz

PROTOCOL='solvent_guided_common_pool_ranking_replay_v1'
TOLERANCE=0.1
PROBE_PROTOCOL='solvent_selected_secant_angular_probes_v1'
PROBE_SETTINGS={'maximum_step_radian':.10,'geometry_only_step_sequence_radian':[.10,.05,.025,.0125],
    'direction_zero_tolerance_radian':1e-12,'coordinate_replay_tolerance_A':1e-12,
    'adapter_physical_tolerance_A':1e-10,
    'active_modes':'same_four_source_angular_modes','maximum_angles_radian':.8,
    'maximum_source_heavy_displacement_A':.8,'points_per_metal':2,'source_count':8,
    'maximum_new_geometries':32,'maximum_new_MACE_calls':64,'maximum_new_GFN2_calls':128,
    'direction_rule':'composite_minus_native_if_different_else_composite_minus_origin',
    'selection_tolerance_kcal_mol':.1,'geometry_checks':angular.SETTINGS}


def probe_direction(center,native,task,winners_differ):
    """Finite source-connected secant, never a claimed force or gradient."""
    center=np.asarray(center,dtype=float);native=np.asarray(native,dtype=float)
    if center.shape!=(task['mode_count'],) or native.shape!=center.shape:
        raise InvalidArtifact('full source-coordinate dimension differs')
    inactive=sorted(set(range(task['mode_count']))-set(task['active_indices']))
    if np.any(center[inactive]!=0) or np.any(native[inactive]!=0):
        raise InvalidArtifact('center or native winner uses a different physical subspace')
    d=center-native if winners_differ else center.copy()
    kind='native_to_composite_winner' if winners_differ else 'origin_to_agreed_winner'
    norm=float(np.max(np.abs(d)))
    if norm<=1e-12:return {'status':'unavailable','reason':'zero source-connected direction','direction_kind':kind,'direction':None}
    return {'status':'available','direction_kind':kind,'direction':(d/norm).tolist()}


def admissible_probe(kin,task,center,direction,sign,symbols):
    attempts=[]
    for step in PROBE_SETTINGS['geometry_only_step_sequence_radian']:
        q=np.asarray(center)+sign*step*np.asarray(direction)
        try:
            check=angular.final_geometry(kin,task,q[task['active_indices']],symbols)
            if check['maximum_heavy_displacement_A']>.8+PROBE_SETTINGS['adapter_physical_tolerance_A']:
                raise InvalidArtifact('finite adapter physical tolerance exceeded')
            attempts.append({'maximum_angle_step_radian':step,'status':'admissible'})
            return {'status':'available','full_q':q.tolist(),'geometry_checks':check,
                    'step_radian':step,'backoff_attempts':attempts}
        except InvalidArtifact as exc:attempts.append({'maximum_angle_step_radian':step,'status':'rejected_geometry','reason':str(exc)})
    return {'status':'unavailable','reason':'no admissible step in frozen geometry-only sequence','full_q':None,
            'step_radian':None,'backoff_attempts':attempts}


@lru_cache(maxsize=None)
def pinned(path,sha256):
    return read_json(verify({'path':path,'sha256':sha256}))


def pinread(pin):return pinned(pin['path'],pin['sha256'])


def actual_point(case,candidate_id,task,kin):
    cell=case['matrix']['Ca'][candidate_id]
    if candidate_id=='origin':
        q=np.zeros(task['mode_count']);point_pin=None
    else:
        # Existing evaluation record accompanies its actual generated XYZ.
        p=Path(verify(next(c for c in case['candidates'] if c['id']==candidate_id)['xyz'])).parent/'result.json'
        point_pin=record(p);point=read_json(p);q=np.asarray(point['full_q'],dtype=float)
        if point['status']!='complete':raise InvalidArtifact('archived source point incomplete')
        coords=xyz(verify(point['coordinate']))
        target=xyz(verify(cell['xyz']))
        if [a[0] for a in coords[1:]]!=[a[0] for a in target[1:]] or np.max(np.abs(np.asarray([a[1:] for a in coords])-np.asarray([a[1:] for a in target])))>1e-12:
            raise InvalidArtifact('point geometry differs from cross-scored geometry')
    coords=kin.evaluate(q)[1];target=np.asarray([a[1:] for a in xyz(verify(cell['xyz']))])
    difference=float(np.max(np.abs(coords-target)))
    if difference>1e-12:raise InvalidArtifact('actual source coordinates cannot replay in common physical map')
    return {'candidate_id':candidate_id,'full_q':q.tolist(),'point':point_pin,'scored_coordinate':cell['xyz'],
            'physical_coordinate_replay_difference_A':difference}


def prepare(inputs,rankings,agreement,probe_plan,output):
    given=read_json(inputs);rank=read_json(rankings)
    if rank['protocol_id']!=PROTOCOL or rank['settings']['origin_selection_tolerance_kcal_mol']!=.1:
        raise InvalidArtifact('archived ranking policy differs')
    if len(given['cases'])!=8 or len({c['case_id'] for c in given['cases']})!=8:
        raise InvalidArtifact('common8 source denominator differs')
    if given['agreement']!=record(agreement):raise InvalidArtifact('shared execution agreement differs')
    lookup={r['case_id']:r for r in rank['rows']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);rows=[]
    for s in given['cases']:
        col=pinread(s['pool_collection']);base=next(c for c in col['cases'] if c['case_id']==s['case_id'])
        prop=pinread(s['proposal_collection']);pm=pinread(prop['manifest'])
        tasks={z:next(t for t in pm['tasks'] if (t['case_id'],t['metal'])==(s['case_id'],z)) for z in ('Ca','La')}
        if col['manifest']!=s['pool_manifest'] or base['pool']['status']!='available':raise InvalidArtifact('complete actual source pool required')
        if s['pool_collection'] not in rank['source_collections']:raise InvalidArtifact('ranking and geometry pool differ')
        if (tasks['Ca']['active_mode_ids']!=tasks['La']['active_mode_ids'] or
                pinread(tasks['Ca']['mapping'])!=pinread(tasks['La']['mapping']) or
                tasks['Ca']['source_preparation']!=tasks['La']['source_preparation']):
            raise InvalidArtifact('paired physical representation differs')
        kin=Kinematics(pinread(tasks['Ca']['mapping'])['context']);atoms=xyz(verify(tasks['Ca']['xyz']));symbols=[a[0] for a in atoms]
        probes=[];candidates=[]
        for z in ('Ca','La'):
            ranks=row_rank(base,z)
            if ranks!=lookup[s['case_id']]['metals'][z]:raise InvalidArtifact('frozen ranking differs')
            t=tasks[z]
            center=actual_point(base,ranks['composite_operational'],t,kin)
            native=actual_point(base,ranks['native_operational'],t,kin)
            direction=probe_direction(center['full_q'],native['full_q'],t,ranks['operational_choices_differ'])
            for sign,name in [(-1,'minus'),(1,'plus')]:
                pid='solvent_'+z+'_'+name
                r={'id':pid,'generated_for_metal':z,'center':center,'native_winner':native,'direction':direction,'sign':sign}
                if direction['status']!='available':r.update(status='unavailable',reason=direction['reason'],full_q=None)
                else:r.update(admissible_probe(kin,t,center['full_q'],direction['direction'],sign,symbols))
                if r['status']=='available':
                    cd=out/'coordinates'/s['case_id'];cd.mkdir(parents=True,exist_ok=True)
                    xp=cd/(pid+'.xyz');coords=kin.evaluate(r['full_q'])[1]
                    write_xyz(xp,[(a,*map(float,p)) for a,p in zip(symbols,coords)])
                    r['coordinate']=record(xp)
                    candidates.append({'id':pid,'full_q':r['full_q'],'coordinate':r['coordinate']})
                probes.append(r)
        rows.append({'case_id':s['case_id'],'status':'prepared' if candidates else 'unavailable',
            'reason':None if candidates else 'no admissible geometry in frozen probe rule',
            'candidates':candidates,'declared_probe_count':4,'probe_coverage':{'available':len(candidates),'unsupported':4-len(candidates)},
            'probes':probes,'base_collection':s['pool_collection'],'source_proposal_manifest':prop['manifest'],
            'source_tasks':tasks,'known_class_report_only':s['known_class']})
    spec={'branch':'solvent_guided','protocol_id':PROBE_PROTOCOL,'settings':PROBE_SETTINGS,
        'inputs':record(inputs),'agreement':record(probe_plan),'execution_agreement':record(agreement),
        'rankings':record(rankings),'coordinate_limits':{'maximum_angle_radian':.8,'maximum_heavy_displacement_A':.8},
        'maximum_candidates_per_case':4,'cases':rows,'new_molecular_calls':0,'production_changed':False,
        'search_coverage_policy':'unsupported geometric probes retained; required molecular failures never removed',
        'implementation':record(__file__)}
    path=out/'specification.json';write_new(path,spec)
    summary={'specification':record(path),'sources':8,'declared_probes':32,
        'admissible_geometries':sum(len(r['candidates']) for r in rows),
        'unsupported_probes':sum(r['probe_coverage']['unsupported'] for r in rows),
        'new_molecular_calls':0,'maximum_MACE_before_dedup':2*sum(len(r['candidates']) for r in rows),
        'maximum_GFN2_before_dedup':4*sum(len(r['candidates']) for r in rows)}
    write_new(out/'PREPARATION.json',summary);return summary


def row_rank(case,metal):
    order=[q['id'] for q in case['candidates']]
    cells=case['matrix'].get(metal,{})
    absent=[q for q in order if cells.get(q,{}).get('status')!='complete']
    if not order or 'origin' not in order or absent:
        return {'status':'unavailable','reason':case.get('reason') or 'required candidate cell unavailable',
                'missing_candidates':absent,'candidate_order':order}
    works={q:relative_components(cells[q]['components'],cells['origin']['components']) for q in order}
    native=min(order,key=lambda q:works[q]['native_MACE_kcal_mol'])
    composite=min(order,key=lambda q:works[q]['composite_kcal_mol'])
    native_op=native if works[native]['native_MACE_kcal_mol'] < -TOLERANCE else 'origin'
    composite_op=composite if works[composite]['composite_kcal_mol'] < -TOLERANCE else 'origin'
    if case['pool']['status']=='available':
        old=case['pool']['rows'][metal]
        if (old['mathematical_candidate'],old['operational_candidate'],old['work_from_origin_kcal_mol'])!=(composite,composite_op,works):
            raise InvalidArtifact('saved composite choice/work replay differs')
    regret=works[native_op]['composite_kcal_mol']-works[composite]['composite_kcal_mol']
    if regret < -1e-8:raise InvalidArtifact('negative composite regret')
    return {'status':'available','candidate_order':order,'works_kcal_mol':works,
            'native_ranking':sorted(order,key=lambda q:works[q]['native_MACE_kcal_mol']),
            'composite_ranking':sorted(order,key=lambda q:works[q]['composite_kcal_mol']),
            'native_mathematical':native,'native_operational':native_op,
            'composite_mathematical':composite,'composite_operational':composite_op,
            'mathematical_choices_differ':native!=composite,
            'operational_choices_differ':native_op!=composite_op,
            'native_minimum_composite_regret_kcal_mol':works[native]['composite_kcal_mol']-works[composite]['composite_kcal_mol'],
            'native_selected_composite_regret_kcal_mol':regret,
            'selected_composite_energy_difference_kcal_mol':works[native_op]['composite_kcal_mol']-works[composite_op]['composite_kcal_mol'],
            'regret_exceeds_selection_tolerance':regret>TOLERANCE}


def stats(values):
    return {'count':len(values),'minimum':min(values) if values else None,
            'median':statistics.median(values) if values else None,'maximum':max(values) if values else None}


def summarize(rows):
    out={'source_denominator':len(rows),'available_pairs':sum(r['status']=='available' for r in rows),
         'unavailable_pairs':sum(r['status']!='available' for r in rows),'metal_rows':{}}
    for z in ('Ca','La'):
        rr=[r['metals'][z] for r in rows]; ok=[r for r in rr if r['status']=='available']
        out['metal_rows'][z]={'denominator':len(rr),'available':len(ok),'unavailable':len(rr)-len(ok),
            'mathematical_choices_differ':sum(r['mathematical_choices_differ'] for r in ok),
            'operational_choices_differ':sum(r['operational_choices_differ'] for r in ok),
            'regret_exceeds_0_1':sum(r['regret_exceeds_selection_tolerance'] for r in ok),
            'regret_kcal_mol':stats([r['native_selected_composite_regret_kcal_mol'] for r in ok]),
            'selected_native_candidates':dict(Counter(r['native_operational'] for r in ok)),
            'selected_composite_candidates':dict(Counter(r['composite_operational'] for r in ok))}
    valid=[r for r in rows if r['status']=='available']
    out['delta_R_composite_vs_native_selection_kcal_mol']=stats([r['delta_R_composite_vs_native_selection_kcal_mol'] for r in valid])
    out['abs_delta_R_exceeds_0_1']=sum(abs(r['delta_R_composite_vs_native_selection_kcal_mol'])>TOLERANCE for r in valid)
    out['unavailable_reasons']=dict(Counter(r['reason'] for r in rows if r['status']!='available'))
    return out


def replay(collections,agreement,output):
    bundle=load_results(collections)
    if bundle['protocol_id']!=SCALED_PROTOCOL:raise InvalidArtifact('not the declared completed scaled pool')
    if {m['population'] for m in bundle['manifests']}!={'scaled30','scaled225'} or len(bundle['rows'])!=255:
        raise InvalidArtifact('exact original30 and primary225 union required')
    membership={c['case_id']:m['population'] for m in bundle['manifests'] for c in m['cases']}
    # Canonical and folded records identify the same biological groups through
    # their actual recorded root_case_id, without treating them as new proteins.
    by_root={c['old_result']['root_case_id']:c['old_result']['biological_group']
             for c in bundle['rows'].values() if c['old_result'].get('root_case_id') and c['old_result'].get('biological_group')}
    rows=[]
    for cid,c in bundle['rows'].items():
        old=c['old_result']; metals={z:row_rank(c,z) for z in ('Ca','La')}
        good=all(r['status']=='available' for r in metals.values())
        r={'case_id':cid,'population':membership[cid],
           'biological_group':old.get('biological_group',by_root.get(cid,old.get('root_case_id',cid))),
           'source_conditioning_metal':old.get('source_conditioning_metal','original30'),
           'expected_class_report_only':old.get('expected_class'),'label_scope':old.get('label_scope'),
           'status':'available' if good else 'unavailable','reason':c['pool'].get('reason'),
           'metals':metals,'composite_R_at_native_selected':None,'composite_R_at_composite_selected':None,
           'delta_R_composite_vs_native_selection_kcal_mol':None}
        if good:
            chosen={mode:{z:c['matrix'][z][metals[z][mode+'_operational']]['components'] for z in ('Ca','La')} for mode in ('native','composite')}
            rn=score(chosen['native']['Ca'],chosen['native']['La'])['composite_R_model_kcal_mol']
            rc=score(chosen['composite']['Ca'],chosen['composite']['La'])['composite_R_model_kcal_mol']
            algebra=metals['La']['selected_composite_energy_difference_kcal_mol']-metals['Ca']['selected_composite_energy_difference_kcal_mol']
            if abs((rc-rn)-algebra)>1e-7:raise InvalidArtifact('selection score sign/algebra differs')
            r.update(composite_R_at_native_selected=rn,composite_R_at_composite_selected=rc,
                     delta_R_composite_vs_native_selection_kcal_mol=rc-rn)
        rows.append(r)
    groups={}
    for field in ('population','source_conditioning_metal','biological_group'):
        grouped=defaultdict(list)
        for r in rows:grouped[r[field]].append(r)
        groups[field]={k:summarize(v) for k,v in sorted(grouped.items())}
    result={'protocol_id':PROTOCOL,'agreement':record(agreement),'source_collections':bundle['collections'],
        'settings':{'origin_selection_tolerance_kcal_mol':TOLERANCE,'tie_policy':'archived_candidate_order'},
        'rows':rows,'summary':summarize(rows),'groups':groups,'new_molecular_calls':0,
        'no_label_fit_or_new_reference':True,'structural_replicates_are_not_independent_biology':True}
    write_new(output,result);return {'output':record(output),'summary':result['summary'],'populations':groups['population']}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='operation',required=True)
    a=sub.add_parser('replay');a.add_argument('--collections',nargs='+',required=True);a.add_argument('--agreement',required=True);a.add_argument('--output',required=True)
    a=sub.add_parser('prepare');a.add_argument('--inputs',required=True);a.add_argument('--rankings',required=True);a.add_argument('--agreement',required=True);a.add_argument('--probe-plan',required=True);a.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation');print(json.dumps(globals()[op](**args)))


if __name__=='__main__':main()
