"""Paired comparison of real frozen and responsive-trial density GK candidates."""
from __future__ import annotations
import argparse
import json
import math
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_density_gk_hybrid import PROTOCOL,TRIAL_PROTOCOL,TOL,validate


def compare(parent_path,trial_path,output):
    parent=read_json(parent_path);trial=read_json(trial_path)
    if parent['protocol']!=PROTOCOL or trial['protocol']!=TRIAL_PROTOCOL or not parent['complete'] or not trial['complete']:
        raise InvalidArtifact('complete real candidate pair required; no fallback comparison')
    pm=validate(verify(parent['manifest']));tm=validate(verify(trial['manifest']))
    audit=read_json(verify(tm['sources']['trial_source_audit']))
    if audit['parent']!=record(parent_path):raise InvalidArtifact('trial was not prepared from this frozen parent')
    parent_tasks={t['task_id']:t for t in pm['static_tasks']+pm['response_tasks']}
    for t in tm['static_tasks']+tm['response_tasks']:
        old=parent_tasks[t['task_id']]
        if any(t[k]['sha256']!=old[k]['sha256'] for k in ('xyz','key','mask')):
            raise InvalidArtifact('matched physical state/cavity/solver changed')
    pc=parent['variants']['primary']['cases'];tc=trial['variants']['primary']['cases'];cases={}
    if set(pc)!=set(tc) or len(pc)!=4:raise InvalidArtifact('candidate comparison inventory differs')
    ordinary=('direct_density_kcal','GK_permanent_transfer_kcal','environment_induction_transfer_kcal','short_context_kcal')
    for case in pc:
        a,b=pc[case],tc[case];pa=a['components_R_kcal'];pb=b['components_R_kcal']
        if a['evidence']!=b['evidence']:raise InvalidArtifact('evidence label/stratum changed')
        if set(pa)!={'DFT_vacuum_kcal',*ordinary} or set(pb)!={'DFT_intrinsic_trial_kcal',*ordinary}:
            raise InvalidArtifact('missing/incompatible candidate energy component')
        components={'intrinsic_core_response_kcal':pb['DFT_intrinsic_trial_kcal']-pa['DFT_vacuum_kcal']}
        components.update({k:pb[k]-pa[k] for k in ordinary})
        if components['short_context_kcal']!=0.:raise InvalidArtifact('learned context changed in density comparison')
        delta=b['R_kcal']-a['R_kcal'];error=delta-math.fsum(components.values())
        if abs(error)>TOL['algebra_kcal']:raise InvalidArtifact('paired score-change algebra does not close')
        cases[case]=dict(frozen_R_kcal=a['R_kcal'],trial_R_kcal=b['R_kcal'],change_R_kcal=delta,
            change_components_kcal=components,algebra_error_kcal=error,evidence=a['evidence'])
    contrasts=[]
    for a,b in zip(parent['variants']['primary']['contrasts'],trial['variants']['primary']['contrasts']):
        if (a['alpha'],a['GGR'])!=(b['alpha'],b['GGR']):raise InvalidArtifact('relative comparison pairing changed')
        contrasts.append(dict(alpha=a['alpha'],GGR=a['GGR'],frozen_difference_kcal=a['difference_kcal'],
            trial_difference_kcal=b['difference_kcal'],change_kcal=b['difference_kcal']-a['difference_kcal'],
            frozen_direction_pass=a['pass_'],trial_direction_pass=b['pass_']))
    result=dict(status='computed_development_comparison',parent=record(parent_path),trial=record(trial_path),
        implementation=record(__file__),physical_geometry_cavity_and_solver_matched=True,cases=cases,contrasts=contrasts,
        frozen_partition_kcal=parent['variants']['primary']['partition_kcal'],trial_partition_kcal=trial['variants']['primary']['partition_kcal'],
        numerical_pass=parent['numerical_pass'] and trial['numerical_pass'],trial_radius_sensitivity_pass=trial['radius_sensitivity_pass'],
        trial_partition_pass=trial['partition_pass'],trial_ordering_pass=trial['ordering_pass'],
        independent_biological_groups=sorted({c['evidence']['group'] for c in cases.values()}),
        evidence_use='Consumed method-development groups; no blind validation or refitted threshold',
        new_DFT_calls=0,new_MACE_calls=0,new_native_calls=0,reference=None,calibrated_class=None,baseline_changed=False)
    write_new(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('parent','trial','output'):p.add_argument('--'+name,required=True)
    a=p.parse_args();r=compare(a.parent,a.trial,a.output)
    print(json.dumps({k:r[k] for k in ('status','numerical_pass','trial_partition_pass','trial_ordering_pass')}))
