"""Explicit grouped comparisons against immutable same-model native collections."""
from __future__ import annotations
from affordable_common import InvalidArtifact,read_json,verify
from mace_density_gk_hybrid import TRIAL_PROTOCOL,TOL

PROTOCOL='declared_source_graph_responsive_density_GK_POLAR_panel_v2'
SPEC='frozen_density_hybrid_grouped_comparisons_v1'


def validate_spec(pin,case_ids,software):
    s=read_json(verify(pin));verify(s['plan'])
    if s['protocol']!=SPEC or s['ordering_tolerance_kcal']!=TOL['ordering_kcal']:
        raise InvalidArtifact('comparison model/tolerance differs')
    declared=[i for v in s['ordered_vectors'].values() for i in v]
    if sorted(declared)!=sorted(case_ids) or len(declared)!=len(set(declared)):
        raise InvalidArtifact('comparison vectors omit or duplicate physical sites')
    refs={}
    for name,p in s['references'].items():
        r=read_json(verify(p));m=read_json(verify(r['manifest']))
        if not r['complete'] or not r['numerical_pass'] or not r['radius_sensitivity_pass']:
            raise InvalidArtifact('reference calculation incomplete or numerically unsupported')
        if m.get('scientific_protocol',m['protocol'])!=TRIAL_PROTOCOL or m['software']!=software:
            raise InvalidArtifact('reference energy model/backend differs')
        refs[name]=r
    ids=set()
    for c in s['comparisons']:
        if c['id'] in ids or c['left_new_case'] not in case_ids or c['right_reference'] not in refs:
            raise InvalidArtifact('invalid or duplicate comparison')
        ids.add(c['id'])
        if c['right_case'] not in refs[c['right_reference']]['variants']['primary']['cases']:
            raise InvalidArtifact('reference case unavailable')
        if c['evidence_stratum']!='supporting_cross_study':raise InvalidArtifact('unsupported comparison evidence stratum')
    if not ids:raise InvalidArtifact('empty comparison inventory')
    return s,refs


def comparisons(m,cases,label):
    s,refs=validate_spec(m['sources']['comparison_spec'],m['case_ids'],m['software']);rows=[]
    for c in s['comparisons']:
        left=cases[c['left_new_case']]['R_kcal']
        right=refs[c['right_reference']]['variants'][label]['cases'][c['right_case']]['R_kcal']
        delta=left-right
        rows.append(dict(comparison_id=c['id'],left=c['left_new_case'],right=c['right_case'],
            evidence_stratum=c['evidence_stratum'],difference_kcal=delta,pass_=delta>TOL['ordering_kcal']))
    return dict(status='complete',cases=cases,contrasts=rows,ordering_pass=all(c['pass_'] for c in rows),
        partition_kcal=None,partition_pass=None,partition_status='not_tested_on_new_sites',
        ordered_vectors={group:[dict(case_id=i,R_kcal=cases[i]['R_kcal']) for i in ids] for group,ids in s['ordered_vectors'].items()},
        biological_classification=None)
