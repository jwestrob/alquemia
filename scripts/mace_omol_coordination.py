"""Finite-range matched coordination descriptor; not an ionic binding free energy."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, rotation, write_xyz, accepted_attempt, check_atoms
from mace_omol_readout import source, CASES

PROTOCOL='mace_omol_matched_coordination_descriptor_v1'
DEFINITION={'metal_x_beyond_ligand_max_A':30., 'farther_displacement_A':10.,
            'score':'(E_bound_Ca-E_detached_Ca)-(E_bound_La-E_detached_La)',
            'fragment_ionic_states_certified':False, 'solvent':None}


def bound_source(collection):
    from mace_omol import benchmark_tasks
    from mace_canonical_run import inventory
    c,m=source(collection); tasks=benchmark_tasks(inventory(verify(m['inventory'])),verify(m['mechanics']))
    rows=dict(c['rows']);rows.update({key:r['result'] for key,r in c['reused_rows'].items()})
    if len(tasks)!=64 or set(rows)!={t['task_id'] for t in tasks}:
        raise InvalidArtifact('all64 original bound endpoints required')
    return c,m,tasks,rows


def expected_tasks(parent_tasks, stage):
    tasks=[]
    for parent in parent_tasks:
        if stage=='coordination_qualification' and parent['case_id'] not in CASES:
            continue
        variants=('detached','repeat','farther','rotate') if stage=='coordination_qualification' and parent['case_id']==CASES[0] else ('detached',)
        for variant in variants:
            t=copy.deepcopy(parent);t.pop('cache_key',None)
            t.update(source_task_id=parent['task_id'],task_id=parent['task_id']+'_'+variant,
                     source_xyz=parent['xyz'],variant=variant,capture_native_readout=True,require_isolated_metal=True)
            t.pop('xyz');tasks.append(t)
    return tasks


def transformed(task):
    rows=xyz(verify(task['source_xyz']))
    if rows[0][0]!=task['metal'] or len(rows)<2:
        raise InvalidArtifact('selected metal must be atom zero in a real coordination core')
    pos=np.array([r[1:] for r in rows]);pos[0,0]=max(pos[1:,0])+DEFINITION['metal_x_beyond_ligand_max_A']
    if task['variant']=='farther':pos[0,0]+=DEFINITION['farther_displacement_A']
    if task['variant']=='rotate':pos=(pos-pos[0])@rotation().T+pos[0]
    if min(np.linalg.norm(pos[1:]-pos[0],axis=1))<=6.:
        raise InvalidArtifact('separated geometry has metal edges within native cutoff')
    return [(row[0],*p) for row,p in zip(rows,pos)]


def qualification(path):
    saved=read_json(path);mp=verify(saved['manifest']);m=read_json(mp)
    if m['stage']!='coordination_qualification':raise InvalidArtifact('wrong coordination numerical source')
    validate(mp);actual=collect(mp)
    if saved!=actual or not actual['numerical_gate_pass']:
        raise InvalidArtifact('passing actual detached-reference qualification required')
    return actual,m


def prepare(collection, agreement, output, qualification_collection=None):
    from mace_omol import common, seal
    _,parent,bt,_=bound_source(collection)
    stage='coordination_benchmark' if qualification_collection else 'coordination_qualification'
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,stage)
    m.update(protocol_id=PROTOCOL,source_collection=record(collection),definition=DEFINITION)
    wanted=expected_tasks(bt,stage)
    if qualification_collection:
        qc,qm=qualification(qualification_collection)
        if qm['source_collection']!=m['source_collection'] or qm['model']!=m['model'] or qm['software']!=m['software']:
            raise InvalidArtifact('qualification used a different bound benchmark or model')
        m['qualification']=record(qualification_collection)
        m['reused']={t['task_id']:{'source_collection':record(qualification_collection),'source_task_id':t['task_id']}
                     for t in qm['tasks'] if t['variant']=='detached'}
        wanted=[t for t in wanted if t['task_id'] not in m['reused']]
    for t in wanted:
        p=out/(t['task_id']+'.xyz');write_xyz(p,transformed(t));t['xyz']=record(p)
    m['tasks']=wanted
    return seal(out,m)


def validate(manifest):
    from mace_omol import SCHEMA,TOL,model
    m=read_json(manifest)
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['definition']!=DEFINITION
            or m['tolerances']!=TOL or m['stage'] not in ('coordination_qualification','coordination_benchmark')):
        raise InvalidArtifact('changed coordination descriptor definition/protocol')
    _,parent,bt,_=bound_source(verify(m['source_collection']))
    software=verify(m['software']);s=read_json(software)
    if m['software']!=parent['software'] or m['model']!=parent['model'] or m['model']!=model(software) or m['inventory']!=parent['inventory']:
        raise InvalidArtifact('coordination calculation changed the qualified model/source')
    for pin in [m['agreement'],*m['implementation'].values(),s['python'],s['requirements'],s['checkpoint'],
                s['inspection'],s['download_receipt'],*read_json(verify(s['backend_source_inventory']))['files']]:verify(pin)
    expected=expected_tasks(bt,m['stage'])
    if m['stage']=='coordination_benchmark':
        qc,qm=qualification(verify(m['qualification']))
        if qm['source_collection']!=m['source_collection'] or qm['model']!=m['model']:
            raise InvalidArtifact('qualification bound source differs')
        reuse={t['task_id']:{'source_collection':m['qualification'],'source_task_id':t['task_id']}
               for t in qm['tasks'] if t['variant']=='detached'}
        if m['reused']!=reuse or len(reuse)!=4:raise InvalidArtifact('qualified detached cache differs')
        expected=[t for t in expected if t['task_id'] not in reuse]
    elif m['reused']:raise InvalidArtifact('unexpected qualification cache reuse')
    if len(expected)!=(60 if m['stage']=='coordination_benchmark' else 10) or len(m['tasks'])!=len(expected):
        raise InvalidArtifact('finite detached-reference inventory changed')
    for t,e in zip(m['tasks'],expected):
        if any(t.get(k)!=v for k,v in e.items()):raise InvalidArtifact('detached task differs from frozen source')
        actual=xyz(verify(t['xyz']));wanted=transformed(t)
        if [r[0] for r in actual]!=[r[0] for r in wanted] or np.max(np.abs(np.array([r[1:] for r in actual])-np.array([r[1:] for r in wanted])))>1e-12:
            raise InvalidArtifact('unapproved detached geometry transformation')
        check_atoms(actual,t['charge'])
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('coordination scientific cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'reused':len(m['reused']),'manifest':record(manifest)}


def collect(manifest):
    from mace_omol import TOL,UNAVAILABLE
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[];reused={}
    for t in m['tasks']:
        found=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:found.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=found[-1] if found else {'status':'unavailable','energy_eV':None}
    complete=all(r['status']=='computed' for r in rows.values())
    if complete:
        for key,r in rows.items():
            force=float(np.max(np.abs(np.load(verify(r['forces']))[0])))
            checks.append({'name':key+'_isolated_metal','max_metal_force_eV_A':force,
                           'edge_count':r['input_state_check']['metal_neighbor_edge_count'],
                           'pass':force<=TOL['force_max_eV_A'] and r['input_state_check']['metal_neighbor_edge_count']==0})
    if complete and m['stage']=='coordination_qualification':
        for variant in ('repeat','farther','rotate'):
            de={}
            for metal in ('La','Ca'):
                a=rows[CASES[0]+'_'+metal+'_detached'];b=rows[CASES[0]+'_'+metal+'_'+variant]
                de[metal]=(b['energy_eV']-a['energy_eV'])*EV_TO_KCAL
                matrix=rotation() if variant=='rotate' else np.eye(3)
                df=float(np.max(np.abs(np.load(verify(b['forces']))@matrix-np.load(verify(a['forces'])))))
                checks.append({'name':metal+'_'+variant,'energy_error_kcal_mol':de[metal],
                               'force_error_eV_A':df,'pass':abs(de[metal])<=TOL['energy_kcal_mol'] and df<=TOL['force_max_eV_A']})
            checks.append({'name':'R_'+variant,'error_kcal_mol':de['Ca']-de['La'],
                           'pass':abs(de['Ca']-de['La'])<=TOL['energy_kcal_mol']})
    if m['stage']=='coordination_benchmark':
        qc,_=qualification(verify(m['qualification']))
        for key,ref in m['reused'].items():reused[key]={**ref,'result':qc['rows'][ref['source_task_id']]}
    return {'status':'complete' if complete else 'incomplete','manifest':record(mp),'rows':rows,'attempts':attempts,
            'checks':checks,'reused_rows':reused,'numerical_gate_pass':complete and bool(checks) and all(c['pass'] for c in checks),**UNAVAILABLE}


def report(collection, original_report, output):
    from mace_omol import TOL
    from mace_canonical_run import inventory
    from mace_canonical_report import fit_bands,decision
    saved=read_json(collection);mp=verify(saved['manifest']);m=read_json(mp);validate(mp)
    if saved!=collect(mp) or m['stage']!='coordination_benchmark':raise InvalidArtifact('actual coordination benchmark collection required')
    original=read_json(original_report)
    if original['sources']['OMOL']!=m['source_collection']:raise InvalidArtifact('original comparator uses another bound benchmark')
    for ref in original['sources'].values():verify(ref)
    _,_,bt,bound=bound_source(verify(m['source_collection']))
    detached=dict(saved['rows']);detached.update({key:r['result'] for key,r in saved['reused_rows'].items()})
    data=inventory(verify(m['inventory']));calibration_names={r['case_id']:r for r in data['rows']}
    baseline={r['case_id']:r for r in original['canonical_rows']};baseline.update(original['nonPQQ_rows'])
    if len(detached)!=64:raise InvalidArtifact('all64 separated endpoints must remain in report')
    scores={}
    for name in dict.fromkeys(t['case_id'] for t in bt):
        original_r=(bound[name+'_Ca']['energy_eV']-bound[name+'_La']['energy_eV'])*EV_TO_KCAL
        if abs(original_r-baseline[name]['R_kcal_mol'])>1e-8:
            raise InvalidArtifact('saved original OMOL contrast differs from executed bound endpoints')
        endpoints={};missing=[]
        for metal in ('La','Ca'):
            key=name+'_'+metal;b=bound[key];d=detached[key+'_detached']
            supported=(d['status']=='computed' and d['input_state_check']['metal_neighbor_edge_count']==0
                       and float(np.max(np.abs(np.load(verify(d['forces']))[0])))<=TOL['force_max_eV_A'])
            value=b['energy_eV']-d['energy_eV'] if b['status']=='computed' and supported else None
            if value is None:missing.append(metal)
            endpoints[metal]={'bound':b,'detached':d,'detached_numerical_support':supported,
                              'D_eV':value,'D_kcal_mol':None if value is None else value*EV_TO_KCAL}
        value=None if missing else (endpoints['Ca']['D_eV']-endpoints['La']['D_eV'])*EV_TO_KCAL
        source_row=calibration_names.get(name)
        scores[name]={'case_id':name,'status':'unavailable' if missing else 'computed','R_kcal_mol':value,
                      'missing_endpoints':missing,'endpoints':endpoints,'original_OMOL_R_kcal_mol':baseline[name]['R_kcal_mol'],
                      'baseline':baseline[name]['baseline'],'prospectively_blind':False,
                      'expected_class':source_row['expected_class'] if source_row else None,
                      'evaluation_role':source_row['evaluation_role'] if source_row else 'retrospective_nonPQQ_direction',
                      'sequence_accession_group':source_row['sequence_accession_group'] if source_row else name.split('_')[0],
                      'evidence_stratum':source_row['evidence_stratum'] if source_row else baseline[name]['evidence'],
                      'calibrated_class':None}
    canonical=[scores[r['case_id']] for r in data['rows']];bands=fit_bands(canonical)
    for r in canonical:r['calibrated_class']=decision(r['R_kcal_mol'],bands['decision_bands'])
    transfer=[r for r in canonical if r['evaluation_role']!='calibration']
    contrasts=[]
    for ggr in CASES:
        for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
            a,g=scores[alpha]['R_kcal_mol'],scores[ggr]['R_kcal_mol'];v=None if a is None or g is None else a-g
            contrasts.append({'positive_case':alpha,'negative_case':ggr,'delta_R_kcal_mol':v,
                              'role':'primary' if ggr==CASES[0] else 'representation_robustness','required_min_kcal_mol':.02,
                              'pass':v is not None and v>.02})
    g0,g1=(scores[name]['R_kcal_mol'] for name in CASES)
    result={'status':saved['status'],'protocol_id':PROTOCOL,'definition':DEFINITION,
            'sources':{'collection':record(collection),'original_report':record(original_report),'agreement':m['agreement']},
            'numerical_gate_pass':saved['numerical_gate_pass'],'calibration':bands,'canonical_rows':canonical,
            'canonical_gate_pass':saved['numerical_gate_pass'] and bands['status']=='calibratable' and all(r['calibrated_class']==r['expected_class']+'-supported' for r in transfer),
            'nonPQQ_rows':{name:r for name,r in scores.items() if name not in calibration_names},'nonPQQ_contrasts':contrasts,
            'nonPQQ_primary_gate_pass':all(r['pass'] for r in contrasts if r['role']=='primary'),
            'nonPQQ_robustness_gate_pass':all(r['pass'] for r in contrasts),
            'GGR_representation_shift_kcal_mol':None if g0 is None or g1 is None else g1-g0,
            'broad_affinity_validated':False,'baseline_changed':False,'binding_free_energy_kcal_mol':None}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/Path(__file__).name)
    result['implementation']=record(out/Path(__file__).name);write_new(out/'result.json',result)
    lines=['# Matched OMOL coordination descriptor','',f'Numerical gate: {result["numerical_gate_pass"]}. Canonical gate: {result["canonical_gate_pass"]}. Non-PQQ primary/robustness: {result["nonPQQ_primary_gate_pass"]}/{result["nonPQQ_robustness_gate_pass"]}.',
           '',f'Canonical calibration gap: {bands["gap_kcal_mol"]} kcal/mol.',
           '', '| Case | Coordination contrast, kcal/mol | Own calibration decision |','|---|---:|---|']
    for r in [*transfer,*result['nonPQQ_rows'].values()]:lines.append(f'| {r["case_id"]} | {r["R_kcal_mol"]} | {r["calibrated_class"]} |')
    lines+=['','| Non-PQQ comparison | Difference, kcal/mol | Pass |','|---|---:|---|']
    for r in contrasts:lines.append(f'| {r["positive_case"]} minus {r["negative_case"]} | {r["delta_R_kcal_mol"]} | {r["pass"]} |')
    lines+=['','These are consumed development comparisons. Both crystals overlap calibration accessions, two alpha structures are one biological group, and both GGR representations are another. PQQ composition confounding remains.',
            'The finite-range coordination descriptor cancels geometry-independent terms but does not certify fragment ionic states or physical dissociation energies. No solvent, aquo offset, entropy, training or fitted correction is added. Production remains unchanged.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for key in ('collection','agreement','output'):a.add_argument('--'+key,required=True)
    a.add_argument('--qualification-collection')
    a=sub.add_parser('report')
    for key in ('collection','original-report','output'):a.add_argument('--'+key,required=True)
    args=vars(p.parse_args());command=args.pop('command');r={'prepare':prepare,'report':report}[command](**args)
    print(json.dumps({k:r[k] for k in ('status','manifest','tasks','canonical_gate_pass','nonPQQ_primary_gate_pass','nonPQQ_robustness_gate_pass','GGR_representation_shift_kcal_mol') if k in r},indent=2))
