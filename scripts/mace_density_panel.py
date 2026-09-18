"""Declared source-graph panels using the frozen responsive-density GK hybrid."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new
from mace_density_gk_hybrid import TRIAL_PROTOCOL,TOL,VARIANTS,NULLS,key,derive,rotation
from mace_density_boundary_panel import validate as boundary_manifest
from mace_density_observations import validate as observations
from mace_density_short import collect as short_collect

PROTOCOL='declared_source_graph_responsive_density_GK_POLAR_panel_v1'


def inputs(boundary,environment,density,short_report,reference):
    br,er,dr,sr,rr=map(read_json,(boundary,environment,density,short_report,reference))
    bm=boundary_manifest(verify(br['manifest']));em=boundary_manifest(verify(er['manifest']))
    dm=observations(verify(dr['manifest']));sm=read_json(verify(sr['manifest']))
    if not all(r['complete'] and r['gates_pass'] for r in (br,er)) or bm['mode']!='source' or em['mode']!='zero_source':
        raise InvalidArtifact('qualified source and zero-source boundaries required')
    if bm['preparation']!=em['preparation'] or bm['frameworks']!=em['frameworks'] or bm['charges']!=record(density):
        raise InvalidArtifact('physical preparation/density mismatch')
    if dm['environment_boundary']!=record(environment) or sm['preparation']!=bm['preparation']:
        raise InvalidArtifact('density/short environment mismatch')
    qr=read_json(verify(dm['quantum']))
    if sm['quantum_manifest']!=qr['manifest']:raise InvalidArtifact('density/short quantum state mismatch')
    if dr['status']!='complete' or not all(dr[k] for k in ('fit_quality_gate_pass','projection_gate_pass','numerical_pass')):
        raise InvalidArtifact('density observation gates failed or unavailable')
    actual=short_collect(verify(sr['manifest']))
    if sr!=actual or sr['status']!='complete':raise InvalidArtifact('short report differs from actual execution')
    if rr['protocol']!=TRIAL_PROTOCOL or not rr['complete'] or not rr['numerical_pass'] or not rr['radius_sensitivity_pass']:
        raise InvalidArtifact('compatible frozen reference comparison required')
    # References are immutable raw contrasts, not an inherited classification band.
    rp=read_json(verify(rr['manifest']))
    for pin in rp['implementation'].values():verify(pin)
    if set(bm['cases'])!=set(em['cases']) or set(sr['cases'])!=set(bm['cases']):
        raise InvalidArtifact('case inventory differs across representations')
    terms={}
    for case in bm['cases']:
        meta=read_json(verify(bm['cases'][case]));other=read_json(verify(em['cases'][case]))
        if meta!=other:raise InvalidArtifact('source/reference physical boundary differs')
        for metal in ('Ca','La'):
            tid=case+'_'+metal;row=dr['rows'][tid];dt=next(t for t in dm['tasks'] if t['task_id']==tid)
            bt=next(t for t in bm['tasks'] if t['case_id']==case and t['metal']==metal)
            if verify(dt['xyz']).read_bytes()!=verify(bt['source_QM_xyz']).read_bytes() or dt['source_receipt']!=qr['rows'][tid]['receipt']:
                raise InvalidArtifact('density/boundary actual quantum source differs')
            value=dt['quantum_energy_hartree']*HA_TO_KCAL-row['old_generating_interaction_kcal']
            if value!=row['intrinsic_trial_core_kcal']:raise InvalidArtifact('old-field subtraction differs')
            short=sr['cases'][case]['endpoints'][metal]
            terms[tid]=dict(intrinsic_core_kcal=value,embedded_energy_hartree=dt['quantum_energy_hartree'],
                old_field_coupling_kcal=row['old_generating_interaction_kcal'],short_context_kcal=short['short_context_kcal'],
                short_core_eV=short['core_eV'],short_full_eV=short['full_eV'],source_output=dt['source_output'],
                source_receipt=dt['source_receipt'],evidence=meta['evidence'],evidence_use=meta['evidence_use'])
    return br,bm,er,em,dr,dm,rr,terms


def prepare(config,output):
    start=time.monotonic();cpu=time.process_time();c=read_json(config)
    if c['scientific_protocol']!=TRIAL_PROTOCOL or c['tolerances']!=TOL:raise InvalidArtifact('frozen model settings changed')
    args={k:verify(v) for k,v in c['inputs'].items()}
    br,bm,er,em,dr,dm,rr,terms=inputs(**args)
    cases=sorted(bm['cases'])
    if cases!=c['case_ids'] or c['identity_case'] not in cases:raise InvalidArtifact('declared case inventory differs')
    if 'comparison_spec' in c:
        from mace_density_comparisons import validate_spec
        validate_spec(c['comparison_spec'],cases,c['software'])
    elif set(c['reference_cases'])!={'ALPHA_1F6S','ALPHA_6IP9'}:raise InvalidArtifact('both frozen alpha references required')
    sw=read_json(verify(c['software']));parent=read_json(verify(rr['manifest']))
    if c['software']!=parent['software'] or sw['returncode']:raise InvalidArtifact('frozen native backend differs')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    paths={n:verify(p) for src in (bm,em,dm) for n,p in src['implementation'].items()}
    for n in ('mace_density_panel.py','mace_density_gk_hybrid.py','mace_density_short.py','mace_density_comparisons.py'):
        paths[n]=Path(__file__).with_name(n)
    for n,p in paths.items():shutil.copyfile(p,impl/n)
    write_new(root/'terms.json',terms)
    protocol=PROTOCOL
    if 'comparison_spec' in c:
        from mace_density_comparisons import PROTOCOL as comparison_protocol
        protocol=comparison_protocol
    m=dict(protocol=protocol,scientific_protocol=TRIAL_PROTOCOL,case_ids=cases,identity_case=c['identity_case'],
        reference_cases=c.get('reference_cases',[]),sources={**c['inputs'],'config':record(config)},software=c['software'],plan=c['plan'],
        implementation={p.name:record(p) for p in impl.glob('*.py')},terms=record(root/'terms.json'),threads=64,tolerances=TOL,
        static_tasks=[],response_tasks=[],new_energy_calls=12*len(cases)+3,new_field_queries=12*len(cases),
        new_response_solves=15*len(cases),new_DFT_calls=0,new_MACE_calls=0,**NULLS)
    if 'comparison_spec' in c:m['sources']['comparison_spec']=c['comparison_spec']
    def task(case,state,variant,mode,eps):
        metal='Ca' if state=='environment' else state
        src,report=(em,er) if state=='environment' else (bm,br)
        old=next(t for t in src['tasks'] if t['case_id']==case and t['metal']==metal)
        meta=read_json(verify(old['boundary']))
        rid=mode+'_'+variant+'_'+case+'_'+state+('_standard' if eps==1e-7 else '')
        d=root/'tasks'/rid;d.mkdir(parents=True)
        mat=rotation() if variant=='rigid' else np.eye(3);shift=np.array([17.,-23.,9.]) if variant=='rigid' else np.zeros(3)
        coords=np.array(meta['positions_A'])@mat.T+shift
        lines=verify(old['xyz']).read_text().splitlines();rows=[lines[0]]
        if len(lines)!=len(coords)+1:raise InvalidArtifact('native physical coordinate inventory differs')
        for line,pos in zip(lines[1:],coords):
            f=line.split();rows.append(' '.join(f[:2]+[format(v,'.17g') for v in pos]+f[5:]))
        (d/'framework.xyz').write_text('\n'.join(rows)+'\n')
        diameter=format(3.6497*VARIANTS.get(variant,1.),'.17g');opts=[]
        for line in verify(old['key']).read_text().splitlines():
            f=line.split()
            if f[0]=='polar-eps':line='polar-eps '+str(eps)
            if f[0]=='solute' and f[1] in ('4998','4999'):f[4]=diameter;line=' '.join(f)
            opts.append(line)
        (d/'framework.key').write_text('\n'.join(opts)+'\n')
        for k in ('mask','overrides'):shutil.copyfile(verify(old[k]),d/Path(old[k]['path']).name)
        dt=next(t for t in dm['tasks'] if t['task_id']==case+'_'+metal);drt=dr['rows'][dt['task_id']]
        t=dict(task_id=rid,case_id=case,state=state,metal=metal,variant=variant,mode=mode,poleps=eps,
            rotation=mat.tolist(),translation_A=shift.tolist(),metal_GK_diameter=diameter,
            xyz=record(d/'framework.xyz'),key=record(d/'framework.key'),mask=record(d/'framework.freeze'),overrides=record(d/'framework.charges'),
            prepared_result=report['tasks'][old['task_id']]['result'],boundary=old['boundary'],frozen_indices=old['frozen_indices'],
            density_task={**dt,'field_arrays':drt['arrays']},density_arrays=drt['arrays'])
        t['cache_key']=key(t,m);return t
    for var in VARIANTS:
        for case in cases:
            for state in ('environment','Ca','La'):
                m['static_tasks'].append(task(case,state,var,'static',1e-9))
                m['response_tasks'].append(task(case,state,var,'solve',1e-9))
                if var=='primary':m['response_tasks'].append(task(case,state,var,'solve',1e-7))
    for state in ('environment','Ca','La'):m['static_tasks'].append(task(c['identity_case'],state,'identity','identity',1e-9))
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol'] not in (PROTOCOL,'declared_source_graph_responsive_density_GK_POLAR_panel_v2') or m['scientific_protocol']!=TRIAL_PROTOCOL or m['tolerances']!=TOL or m['threads']!=64:
        raise InvalidArtifact('panel model/settings differ')
    for pin in [m['software'],m['plan'],m['terms'],*m['sources'].values(),*m['implementation'].values()]:verify(pin)
    c=read_json(verify(m['sources']['config']));n=len(m['case_ids'])
    if any(m[k]!=c[k] for k in ('case_ids','identity_case','scientific_protocol','software','plan','tolerances')) or m['reference_cases']!=c.get('reference_cases',[]):
        raise InvalidArtifact('declared panel scope differs')
    if 'comparison_spec' in c:
        from mace_density_comparisons import validate_spec,PROTOCOL as comparison_protocol
        if m['protocol']!=comparison_protocol or m['sources'].get('comparison_spec')!=c['comparison_spec']:raise InvalidArtifact('comparison scope differs')
        validate_spec(c['comparison_spec'],m['case_ids'],m['software'])
    elif m['protocol']!=PROTOCOL:raise InvalidArtifact('missing comparison specification')
    if {k:m['sources'][k] for k in c['inputs']}!=c['inputs']:raise InvalidArtifact('panel source selection differs')
    if [m['new_energy_calls'],m['new_field_queries'],m['new_response_solves']]!=[12*n+3,12*n,15*n] or m['new_DFT_calls'] or m['new_MACE_calls']:
        raise InvalidArtifact('panel calls differ')
    sw=read_json(verify(m['software']))
    for k in ('executable','library','original_source','original_routine','derived_routine','frontend','implementation','log'):verify(sw[k])
    a,b=derive(verify(sw['original_source']).read_text())
    if a!=verify(sw['original_routine']).read_text() or b!=verify(sw['derived_routine']).read_text():raise InvalidArtifact('native solver equations changed')
    expected={(case,state,var,'static',1e-9) for case in m['case_ids'] for state in ('environment','Ca','La') for var in VARIANTS}
    expected|={(case,state,var,'solve',1e-9) for case in m['case_ids'] for state in ('environment','Ca','La') for var in VARIANTS}
    expected|={(case,state,'primary','solve',1e-7) for case in m['case_ids'] for state in ('environment','Ca','La')}
    expected|={(m['identity_case'],state,'identity','identity',1e-9) for state in ('environment','Ca','La')}
    tasks=m['static_tasks']+m['response_tasks']
    if len(tasks)!=len(expected) or {(t['case_id'],t['state'],t['variant'],t['mode'],t['poleps']) for t in tasks}!=expected:
        raise InvalidArtifact('panel task inventory differs')
    if len(m['static_tasks'])!=12*n+3 or any(t['mode']=='solve' for t in m['static_tasks']) or any(t['mode']!='solve' for t in m['response_tasks']):
        raise InvalidArtifact('panel execution order differs')
    for t in tasks:
        if t['cache_key']!=key(t,m):raise InvalidArtifact('changed scientific task key')
        for k in ('xyz','key','mask','overrides','prepared_result','boundary','density_arrays'):verify(t[k])
        for k in ('probes','field_arrays','source_receipt','source_output','projection'):verify(t['density_task'][k])
        if 'reuse_component' in t:raise InvalidArtifact('undeclared panel reuse')
    return m


def comparisons(m,cases,label):
    if 'comparison_spec' in m['sources']:
        from mace_density_comparisons import comparisons as grouped_comparisons
        return grouped_comparisons(m,cases,label)
    reference=read_json(verify(m['sources']['reference']))['variants'][label]
    contrasts=[]
    for a in m['reference_cases']:
        for g,r in cases.items():
            value=reference['cases'][a]['R_kcal']-r['R_kcal']
            contrasts.append(dict(alpha=a,GGR=g,difference_kcal=value,pass_=value>TOL['ordering_kcal']))
    return dict(status='complete',cases=cases,contrasts=contrasts,ordering_pass=all(c['pass_'] for c in contrasts),
        partition_kcal=None,partition_pass=None,partition_status='not_tested_on_new_structures')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare');s.add_argument('--config',required=True);s.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect'):
        s=sub.add_parser(op);s.add_argument('--manifest',required=True)
        if op=='execute':s.add_argument('--retry-failed',action='store_true')
        if op=='collect':s.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output);print(json.dumps(dict(tasks=len(r['static_tasks'])+len(r['response_tasks']))))
    elif a.command=='dry-run':r=validate(a.manifest);print(json.dumps(dict(status='pass',static=len(r['static_tasks']),response=len(r['response_tasks']))))
    elif a.command=='execute':
        from mace_density_gk_hybrid import execute
        execute(a.manifest,a.retry_failed)
    else:
        from mace_density_gk_hybrid import collect
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass','radius_sensitivity_pass','ordering_pass')}))
