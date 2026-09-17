"""Declared charged-spectator consistency test using existing OMOL execution."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_global_prepare import CASES
from mace_hybrid import check_atoms, EV_TO_KCAL, write_xyz
from mace_omol_intact import PROTOCOL as SOURCE_PROTOCOL, EVALUATION, collect as collect_endpoints
from mace_omol_locality import endpoint
from mace_omol_panel import product_qualification, ADAPTER_TASK
from mace_omol_edge_run import PRODUCT_CONFIG

PROTOCOL = 'mace_omol_intact_disconnected_sodium_check_v1'
SETTINGS = {'spectator':'Na+', 'offset_from_maximum_x_A':10000.,
            'charge_increment':1, 'score_tolerance_kcal_mol':0.1,
            'maximum_direct_Coulomb_difference_kcal_mol':0.001,
            'Coulomb_kcal_A_per_mol_e2':332.063713299}


def expected(source_report):
    source = read_json(source_report)
    if (source['protocol_id']!=SOURCE_PROTOCOL or source['status']!='complete'
            or not source['numerical_gate_pass'] or set(source['scores'])!=set(CASES)):
        raise InvalidArtifact('completed original five-case report required')
    for key in ('collection','physical_preparation','agreement','report_implementation'):
        verify(source[key])
    tasks=[];coordinates={};reused={};models=[]
    for name in CASES:
        case=source['scores'][name];prep=read_json(verify(case['preparation']))
        original={};original_xyz={};results={}
        for metal in ('La','Ca'):
            for position in ('bound','detached'):
                r=case['endpoints'][metal][position];t,_=endpoint(r)
                key=(metal,position);original[key]=t;original_xyz[key]=xyz(verify(t['xyz']));results[key]=r
                models.append(read_json(verify(r['manifest'])))
                if t['charge']!=prep['endpoints'][metal]['charge']:
                    raise InvalidArtifact('original endpoint differs from protein state')
        x=max(a[1] for rows in original_xyz.values() for a in rows)+SETTINGS['offset_from_maximum_x_A']
        index=case['metal_index'];bound=np.array(original_xyz[('La','bound')][index][1:])
        spectator=np.array([x,bound[1],bound[2]])
        # Bound and detached states have identical spectator geometry for both metals.
        rb=float(np.linalg.norm(spectator-bound))
        rd=float(np.linalg.norm(spectator-np.array(original_xyz[('La','detached')][index][1:])))
        direct=SETTINGS['Coulomb_kcal_A_per_mol_e2']*abs(1/rb-1/rd)
        if direct>SETTINGS['maximum_direct_Coulomb_difference_kcal_mol']:
            raise InvalidArtifact('spectator construction has excessive direct Coulomb contrast')
        for (metal,position),old in original.items():
            rows=original_xyz[(metal,position)]+[('Na',*map(float,spectator))]
            q=old['charge']+1
            if not -10<=q<=10:raise InvalidArtifact('spectator state outside declared training charge range')
            task=copy.deepcopy(old)
            for key in ('xyz','cache_key'):task.pop(key)
            key=f'{name}_{metal}_{position}_sodium'
            task.update(task_id=key,variant='disconnected_sodium',source_xyz=old['xyz'],
                        charge=q,state=check_atoms(rows,q),edge_adapter=ADAPTER_TASK,
                        spectator_index=len(rows)-1,
                        spectator={'element':'Na','formal_charge':1,'xyz_A':spectator.tolist(),
                                   'origin':'declared_constructed_consistency_perturbation',
                                   'minimum_distance_A':float(min(np.linalg.norm(spectator-np.array(a[1:])) for a in rows[:-1])),
                                   'direct_Coulomb_score_bound_kcal_mol':direct},
                        original_endpoint={'manifest':results[(metal,position)]['manifest'],'task_id':old['task_id']},
                        microstate=old['microstate']+'; one disconnected Na+ spectator')
            tasks.append(task);coordinates[key]=rows;reused[key]=results[(metal,position)]
    if len(tasks)!=20:raise InvalidArtifact('declared spectator inventory changed')
    reference=models[0]
    if any(m['inventory']!=reference['inventory'] or m['software']!=reference['software']
           or m['model']['checkpoint']!=reference['model']['checkpoint'] for m in models):
        raise InvalidArtifact('source endpoint methods differ')
    return source,tasks,coordinates,reused,reference


@cached_file_checks
def prepare(source_report,product_equivalence,agreement,output):
    from mace_omol import common,seal
    _,tasks,coordinates,_,old=expected(source_report)
    _,qualified=product_qualification(product_equivalence)
    if qualified['software']!=old['software']:
        raise InvalidArtifact('source software differs from qualified adapter')
    _,out,m=common(verify(old['inventory']),verify(old['software']),agreement,output,'spectator_intact')
    m.update(protocol_id=PROTOCOL,settings=SETTINGS,energy_evaluation=EVALUATION,
             source_report=record(source_report),product_equivalence=record(product_equivalence),tasks=tasks)
    m['model']['execution_adapter']=PRODUCT_CONFIG
    for name in ('mace_omol_edges.py','mace_omol_products.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=qualified['implementation'][name]['sha256']:
            raise InvalidArtifact('qualified execution adapter changed')
    for task in tasks:
        path=out/(task['task_id']+'.xyz');write_xyz(path,coordinates[task['task_id']])
        # The standard writer's historical comment says real geometry; the
        # explicit extension must not be mistaken for a deposited sodium ion.
        lines=path.read_text().splitlines()
        lines[1]='Original protein geometry plus declared constructed distant Na+; see manifest'
        path.write_text('\n'.join(lines)+'\n')
        task['xyz']=record(path)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL,model
    m=read_json(manifest)
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage']!='spectator_intact'
            or m['settings']!=SETTINGS or m['energy_evaluation']!=EVALUATION or m['tolerances']!=TOL):
        raise InvalidArtifact('spectator method/settings changed')
    _,tasks,coordinates,_,old=expected(verify(m['source_report']))
    _,qualified=product_qualification(verify(m['product_equivalence']))
    wanted=model(verify(m['software']));wanted['execution_adapter']=PRODUCT_CONFIG
    if (m['model']!=wanted or m['software']!=old['software'] or m['software']!=qualified['software']
            or m['inventory']!=old['inventory'] or m['reused'] or len(m['tasks'])!=20):
        raise InvalidArtifact('source/model/inventory mismatch')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    for name in ('mace_omol_edges.py','mace_omol_products.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=qualified['implementation'][name]['sha256']:
            raise InvalidArtifact('spectator adapter differs from qualification')
    for actual,task in zip(m['tasks'],tasks):
        if {k:v for k,v in actual.items() if k not in ('xyz','cache_key')}!=task:
            raise InvalidArtifact('spectator task differs from fixed construction')
        rows=xyz(verify(actual['xyz']))
        if rows!=coordinates[task['task_id']] or check_atoms(rows,task['charge'])!=task['state']:
            raise InvalidArtifact('spectator coordinates/state differ')
        payload={k:v for k,v in actual.items() if k!='cache_key'}
        if actual['cache_key']!=cache_key({'task':payload,'model':m['model'],
                                         'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('spectator cache identity mismatch')
    return {'status':'pass','tasks':20,'reused_reference_endpoints':20,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    m=read_json(manifest)
    source,_,_,reused,_=expected(verify(m['source_report']))
    result=collect_endpoints(manifest);scores={};checks=[]
    for name in CASES:
        terms={};components={}
        for metal in ('La','Ca'):
            pair={p:result['rows'][f'{name}_{metal}_{p}_sodium'] for p in ('bound','detached')}
            terms[metal]=(pair['bound']['energy_eV']-pair['detached']['energy_eV']) if all(
                r['status']=='computed' for r in pair.values()) else None
            if terms[metal] is not None:
                from mace_omol_readout import checked_arrays
                tasks={t['task_id']:t for t in m['tasks']}
                arrays={p:checked_arrays(r,tasks[r['task_id']]) for p,r in pair.items()}
                components[metal]={k:float((arrays['bound'][k]-arrays['detached'][k]).sum()) for k in arrays['bound']}
        value=(terms['Ca']-terms['La'])*EV_TO_KCAL if all(v is not None for v in terms.values()) else None
        baseline=source['scores'][name]['R_coord_kcal_mol']
        change=None if value is None else value-baseline
        scores[name]={'R_coord_kcal_mol':value,'original_R_coord_kcal_mol':baseline,
                      'change_kcal_mol':change,'bound_minus_detached_eV':terms,'components_eV':components,
                      'pass':change is not None and abs(change)<=SETTINGS['score_tolerance_kcal_mol']}
        checks.append({'name':name,'change_kcal_mol':change,'pass':scores[name]['pass']})
    contrasts=[]
    for contrast in source['contrasts']:
        a,b=(scores[contrast[k]]['R_coord_kcal_mol'] for k in ('positive_case','negative_case'))
        value=None if a is None or b is None else a-b
        change=None if value is None else value-contrast['delta_R_kcal_mol']
        contrasts.append({**contrast,'modified_delta_R_kcal_mol':value,'change_kcal_mol':change,
                          'consistency_pass':change is not None and abs(change)<=SETTINGS['score_tolerance_kcal_mol']})
        checks.append({'name':contrast['positive_case']+'_minus_'+contrast['negative_case'],
                       'change_kcal_mol':change,'pass':contrasts[-1]['consistency_pass']})
    result.update(protocol_id=PROTOCOL,source_report=m['source_report'],scores=scores,contrasts=contrasts,
                  reused_reference_endpoints=reused,consistency_checks=checks,
                  spectator_consistency_pass=result['numerical_gate_pass'] and all(c['pass'] for c in checks),
                  baseline_changed=False,broad_affinity_validated=False,production_promotion=False)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('source-report','product-equivalence','agreement','output'):
        parser.add_argument('--'+key,required=True)
    print(json.dumps(prepare(**vars(parser.parse_args())),indent=2))
