"""Pinned native OMOL sources and paired compact-solvation comparisons."""
from __future__ import annotations
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms
from mace_omol import COMPONENT, accepted_state

SCHEMA = 'compact_native_OMOL_solvation_source_inventory_v1'
COMPOSITE_PROTOCOL = 'native_OMOL_plus_native_GFN2_ALPB_transfer_v1'
MINIMUM_CALIBRATION_GAP = 0.02


@lru_cache(maxsize=None)
def pinned_manifest(path, sha256):
    return read_json(verify({'path':path,'sha256':sha256}))


def native_endpoint(result, expected_model):
    """Read the actual endpoint receipt and its actual prepared task."""
    pin=record(verify(result['forces']).parent/'result.json')
    actual=read_json(verify(pin))
    if actual!=result:raise InvalidArtifact('embedded native endpoint differs from original receipt')
    source=result['manifest'];m=pinned_manifest(source['path'],source['sha256'])
    task=next(t for t in m['tasks'] if t['task_id']==result['task_id'])
    if (m['model']!=expected_model or result['status']!='computed' or
            result['memory_mode']!='native' or result['energy_component']!=COMPONENT or
            not accepted_state(result,task)):
        raise InvalidArtifact('not a compatible native vacuum OMOL endpoint')
    rows=xyz(verify(task['xyz']));check_atoms(rows,task['charge'])
    return {'xyz':task['xyz'],'charge':task['charge'],'multiplicity':task['spin_multiplicity'],
            'native_MACE_energy_eV':result['energy_eV'],'native_MACE_receipt':pin,
            'source_manifest':source,'source_task_id':task['task_id'],
            'source_protocol_id':m['protocol_id'],'energy_component':COMPONENT,
            'atoms':len(rows),'metal_index':task.get('metal_index',0),'status':'available'}


def representation(endpoints, preparation, protocol_id):
    ca,la=(xyz(verify(endpoints[m]['xyz'])) for m in ('Ca','La'))
    if len(ca)!=len(la) or ca[0][0]!='Ca' or la[0][0]!='La' or any(a[0]!=b[0] for a,b in zip(ca[1:],la[1:])):
        raise InvalidArtifact('paired endpoint composition differs')
    if endpoints['La']['charge']-endpoints['Ca']['charge']!=1:
        raise InvalidArtifact('paired formal charge difference differs')
    return {'protocol_id':protocol_id,'preparation':preparation,'endpoints':endpoints,
            'native_R_model_kcal_mol':(endpoints['Ca']['native_MACE_energy_eV']-endpoints['La']['native_MACE_energy_eV'])*EV_TO_KCAL,
            'paired_nonmetal_coordinate_difference_indices':[i for i,(a,b) in enumerate(zip(ca,la)) if i and a!=b],
            'new_geometry_generated':False}


def inventory(root, output):
    root=Path(root).resolve()
    core_path=root/'diagnostics/mace_omol_20260917/result.json';core=read_json(core_path)
    context_dir=root/'workspaces/environment_pqq_20260919/prepared_v1'
    context_manifest_path=context_dir/'manifest.json';context=read_json(context_manifest_path)
    context_result=read_json(context_dir/'result_1202474.json')
    collection=read_json(context_dir/'collection_1202474.json')
    model=core['model']
    if context['model']!=model:raise InvalidArtifact('core/context checkpoint differs')
    expanded={r['task_id']:r for r in [read_json(verify(p)) for p in collection['results']]}
    expanded.update({r['task']['task_id']:r['accepted'] for r in context['reused']})
    cases=[]
    for original in core['canonical_rows']:
        cid=original['case_id'];prepared=next(c for c in context['cases'] if c['source']['case']==cid)
        reps={}
        for rep in ('core','context'):
            eps={metal:native_endpoint(original['endpoints'][metal] if rep=='core' else expanded[cid+'__expanded__'+metal],model)
                 for metal in ('Ca','La')}
            prep=prepared['source']['parent'] if rep=='core' else prepared['preparation']
            protocol=core['protocol_id'] if rep=='core' else context['protocol_id']
            reps[rep]=representation(eps,prep,protocol)
        expected=next(r for r in context_result['rows'] if r['case']==cid)
        if reps['core']['native_R_model_kcal_mol']!=original['R_kcal_mol'] or reps['context']['native_R_model_kcal_mol']!=expected['context_R_model_kcal_mol']:
            raise InvalidArtifact('native PQQ score replay differs')
        cases.append({'case_id':cid,'biological_group':original['sequence_accession_group'],
                      'sequence_accession_group':original['sequence_accession_group'],
                      'role':original['evaluation_role'],'evidence_stratum':original['evidence_stratum'],
                      'expected_class':original['expected_class'],'label_scope':'canonical_PQQ_functional_class',
                      'all_evidence_consumed':True,'representations':reps})
    sources=[record(core_path),record(context_manifest_path),record(context_dir/'result_1202474.json'),record(context_dir/'collection_1202474.json')]
    for directory,ids,collection_name in [
        (root/'workspaces/second_shell_20260919/prepared_v2',('1F6S','6IP9','1GLG'),'mace_collection_1202083.json'),
        (root/'workspaces/environment_replicas_20260919/prepared_v1',('2FW0','2FVY'),'mace_collection_1202454.json')]:
        cm=read_json(directory/collection_name);dm=read_json(directory/'manifest.json')
        results=[r['accepted'] for r in cm['rows']] if 'rows' in cm else [read_json(verify(p)) for p in cm['results']]
        index={r['task_id']:r for r in results};sources.extend([record(directory/collection_name),record(directory/'manifest.json')])
        for cid in ids:
            prepared=next(s for s in dm['states'] if s['case']==cid);reps={}
            for rep,variant in [('core','core'),('context','expanded')]:
                eps={metal:native_endpoint(index[cid+'__'+variant+'__'+metal],model) for metal in ('Ca','La')}
                reps[rep]=representation(eps,prepared['source']['parent'] if rep=='core' else prepared['preparation'],
                                         eps['Ca']['source_protocol_id'])
            group=prepared['source']['group'];label='La' if group=='alpha_lactalbumin' else 'Ca'
            if group not in ('alpha_lactalbumin','GGR'):raise InvalidArtifact('unexpected direct-direction group')
            cases.append({'case_id':cid,'biological_group':group,'sequence_accession_group':None,
                          'role':'direct_direction_development','evidence_stratum':prepared['source']['evidence_stratum'],
                          'expected_class':label,'label_scope':'relative_alpha_lactalbumin_minus_GGR_direction_only',
                          'all_evidence_consumed':True,'representations':reps})
    if len(cases)!=33 or len({c['case_id'] for c in cases})!=33:raise InvalidArtifact('33-case scope differs')
    result={'schema_version':SCHEMA,'source_records':sources,'model':model,'software':core['software'],'cases':cases,
            'native_PQQ_bands':{'core':core['calibration'],
                                'context':{'Ca_supported_max_R_kcal_mol':context_result['Ca_band_max'],
                                           'La_supported_min_R_kcal_mol':context_result['La_band_min'],
                                           'gap_kcal_mol':context_result['gap_model_kcal_mol'],'calibration_count':25,
                                           'transfer_rows_used':False}},
            'comparison_definitions':{'calibration_cases':25,'consumed_crystal_transfers':3,
                                      'alpha_GGR_pairs':[{'La_like':a,'Ca_like':b,'biological_comparison_group':'alpha_lactalbumin_vs_GGR'}
                                                         for a in ('1F6S','6IP9') for b in ('1GLG','2FW0','2FVY')],
                                      'alpha_GGR_structural_directions':6,'alpha_GGR_biological_comparisons':1,
                                      'PQQ_group_note':'Accession groups retained; crystal replicas overlap calibration. Homologs are not automatically independent biology.'},
            'units':{'native_MACE_endpoint':'eV','low_level_endpoint':'hartree','comparison':'kcal_per_mol_or_model_kcal',
                     'eV_to_kcal_mol':EV_TO_KCAL,'hartree_to_kcal_mol':HA_TO_KCAL},
            'total_case_representations':66,'available_native_MACE_endpoints':132,
            'all_evidence_consumed':True,'absolute_aquo_reference':None,'composite_bands':None,
            'new_scientific_calls':0,'new_geometry_generated':False,'baseline_changed':False}
    write_new(output,result);return result


def mix_pair(mace_eV, vacuum_hartree, alpb_hartree):
    """Pure unit/sign algebra; no unavailable contribution is replaced by zero."""
    for values in (mace_eV,vacuum_hartree,alpb_hartree):
        if set(values)!= {'Ca','La'} or any(not isinstance(v,(float,int)) or not math.isfinite(v) for v in values.values()):
            raise InvalidArtifact('both finite endpoint energies required for every method')
    native=(mace_eV['Ca']-mace_eV['La'])*EV_TO_KCAL
    transfers={m:(alpb_hartree[m]-vacuum_hartree[m])*HA_TO_KCAL for m in ('Ca','La')}
    correction=((alpb_hartree['Ca']-vacuum_hartree['Ca'])-(alpb_hartree['La']-vacuum_hartree['La']))*HA_TO_KCAL
    return {'native_R_model_kcal_mol':native,'endpoint_solvation_kcal_mol':transfers,
            'solvation_delta_R_kcal_mol':correction,'composite_R_model_kcal_mol':native+correction,
            'GFN2_vacuum_R_kcal_mol':(vacuum_hartree['Ca']-vacuum_hartree['La'])*HA_TO_KCAL,
            'GFN2_ALPB_R_kcal_mol':(alpb_hartree['Ca']-alpb_hartree['La'])*HA_TO_KCAL}


def checked_transfer(row, source, inventory_pin):
    from compact_solvation import completed, METHOD
    if row['status']!='complete':return None
    energies={}
    for medium in ('vacuum','alpb'):
        e=row['endpoints'][medium]
        if e['status']!='complete':raise InvalidArtifact('complete transfer contains unavailable endpoint')
        mp=verify(e['manifest']);m=pinned_manifest(str(mp),e['manifest']['sha256'])
        if m['method_id']!=METHOD or m['protocol_id']!=COMPOSITE_PROTOCOL or m['inventory']!=inventory_pin:
            raise InvalidArtifact('low-level endpoint method/inventory differs')
        task=next(t for t in m['all_tasks'] if t['task_id']==e['task_id'])
        if (task['source_endpoint']!=source or task['case_id']!=row['case_id'] or
                task['representation']!=row['representation'] or task['metal']!=row['metal'] or
                task['medium']!=medium or task['solver']!=row['solver']):
            raise InvalidArtifact('low-level endpoint state/source differs')
        actual=completed(mp,e['task_id'])
        if actual is None or any(actual[k]!=e[k] for k in ('energy_hartree','output','receipt','manifest','task_id')):
            raise InvalidArtifact('low-level endpoint receipt/energy does not replay')
        if actual['energy_hartree']!=row[medium+'_hartree']:
            raise InvalidArtifact('row energy differs from original endpoint')
        energies[medium]=actual['energy_hartree']
    if row['delta_solv_hartree']!=energies['alpb']-energies['vacuum']:
        raise InvalidArtifact('reported transfer algebra differs')
    return energies


def numerical_check(paths, primary_lookup, sources, inventory_pin):
    result={'status':'separate_parent_audit_required','collections':[],'rows':[],
            'endpoint_transfer_tolerance_kcal_mol':0.10,'score_correction_tolerance_kcal_mol':0.20,
            'endpoint_transfer_differences_kcal_mol':None,'score_correction_differences_kcal_mol':None}
    if not paths:return result
    expected={(case,'context',metal) for case in ('1H4I','q9z4j7-pqq-la_model') for metal in ('Ca','La')}
    numerical={}
    for path in paths:
        col=read_json(path);result['collections'].append(record(path))
        if col['inventory']!=inventory_pin or col['protocol_id']!=COMPOSITE_PROTOCOL:
            raise InvalidArtifact('numerical collection source differs')
        for row in col['rows']:
            key=(row['case_id'],row['representation'],row['metal'])
            if row['solver']!='ordinary_tight' or key not in expected or key in numerical:
                raise InvalidArtifact('numerical check differs from frozen pilot scope')
            for endpoint in row['endpoints'].values():
                for pin in endpoint.get('available_artifacts',[]):verify(pin)
            if row['status']=='complete':checked_transfer(row,sources[key],inventory_pin)
            numerical[key]=row;result['rows'].append(row)
    if set(numerical)!=expected or any(r['status']!='complete' for r in numerical.values()):
        result['status']='numerical_crosscheck_unavailable';return result
    if any(k not in primary_lookup or primary_lookup[k]['status']!='complete' for k in expected):
        result['status']='numerical_crosscheck_unavailable';return result
    differences={k:(numerical[k]['delta_solv_hartree']-primary_lookup[k]['delta_solv_hartree'])*HA_TO_KCAL for k in expected}
    score={case:differences[(case,'context','Ca')]-differences[(case,'context','La')] for case in ('1H4I','q9z4j7-pqq-la_model')}
    result.update(status='pass' if max(abs(x) for x in differences.values())<=0.10 and max(abs(x) for x in score.values())<=0.20 else 'fail',
                  endpoint_transfer_differences_kcal_mol=[{'case_id':k[0],'representation':k[1],'metal':k[2],'difference':v} for k,v in sorted(differences.items())],
                  score_correction_differences_kcal_mol=score)
    return result


def comparison(inventory_path, collections, solver, numerical_collections=()):
    inv=read_json(inventory_path);inv_pin=record(inventory_path)
    if inv['schema_version']!=SCHEMA or solver not in ('native','ordinary_tight'):
        raise InvalidArtifact('unsupported source inventory or numerical solver')
    if numerical_collections and solver!='native':
        raise InvalidArtifact('the numerical crosscheck must compare against the primary native solver')
    lookup={};col_pins=[];failure_history=[];component_audits=[]
    sources={(c['case_id'],rep,metal):e for c in inv['cases'] for rep,r in c['representations'].items()
             for metal,e in r['endpoints'].items()}
    for path in collections:
        col=read_json(path);col_pins.append(record(path))
        if col['inventory']!=inv_pin or col['protocol_id']!=COMPOSITE_PROTOCOL:
            raise InvalidArtifact('collection source or protocol differs')
        component_audits.append(col.get('component_audit','unavailable'))
        for row in col['rows']:
            if row['solver']!=solver:continue
            key=(row['case_id'],row['representation'],row['metal'])
            if key not in sources:raise InvalidArtifact('collection contains undeclared case/representation')
            if row['status']!='complete':
                failure_history.append({'collection':record(path),'row':row})
                lookup.setdefault(key,row);continue
            checked_transfer(row,sources[key],inv_pin)
            if key in lookup and lookup[key]['status']=='complete':
                old=lookup[key]
                if any(old[k]!=row[k] for k in ('vacuum_hartree','alpb_hartree','delta_solv_hartree')) or any(
                        old['endpoints'][medium]['receipt']!=row['endpoints'][medium]['receipt'] for medium in ('vacuum','alpb')):
                    raise InvalidArtifact('conflicting duplicate actual transfer; select explicit collection')
            lookup[key]=row
    rows=[]
    for c in inv['cases']:
        for rep,r in c['representations'].items():
            eps={m:lookup.get((c['case_id'],rep,m)) for m in ('Ca','La')}
            available=all(e is not None and e['status']=='complete' for e in eps.values())
            row={k:c[k] for k in ('case_id','biological_group','sequence_accession_group','role','evidence_stratum','expected_class','label_scope')}
            row.update(representation=rep,native_status='available',composite_status='available' if available else 'unavailable',
                       native_R_model_kcal_mol=r['native_R_model_kcal_mol'],composite_R_model_kcal_mol=None,
                       solvation_delta_R_kcal_mol=None,endpoint_solvation_kcal_mol=None,
                       GFN2_vacuum_R_kcal_mol=None,GFN2_ALPB_R_kcal_mol=None,
                       low_level_endpoints=eps,source_preparation=r['preparation'],composite_decision=None,
                       unavailable_reason=None if available else {m:'not_run_or_not_collected' if e is None else e['status'] for m,e in eps.items() if e is None or e['status']!='complete'})
            if available:
                energies={m:r['endpoints'][m]['native_MACE_energy_eV'] for m in ('Ca','La')}
                row.update(mix_pair(energies,{m:eps[m]['vacuum_hartree'] for m in ('Ca','La')},
                                    {m:eps[m]['alpb_hartree'] for m in ('Ca','La')}))
            rows.append(row)
    calibration={};pair_matrix=[];context_effects=[]
    for rep in ('core','context'):
        members=[r for r in rows if r['representation']==rep and r['role']=='calibration']
        valid=[r for r in members if r['composite_status']=='available']
        cal={'denominator':25,'available':len(valid),'class_extrema':None,'gap_model_kcal_mol':None,'bands':None,
             'all_cross_class_pairs_ordered':None,'ordered_cross_class_pairs':None,'total_cross_class_pairs':None,
             'class_spread':None,'minimum_gap_model_kcal_mol':MINIMUM_CALIBRATION_GAP,
             'transfer_rows_used':False,'threshold_coefficients_fitted':False}
        if len(valid)==25:
            ca=[r['composite_R_model_kcal_mol'] for r in valid if r['expected_class']=='Ca']
            la=[r['composite_R_model_kcal_mol'] for r in valid if r['expected_class']=='La']
            gap=min(la)-max(ca);ordered=sum(a>b for a in la for b in ca);total=len(ca)*len(la)
            cal.update(class_extrema={'Ca_max':max(ca),'La_min':min(la)},gap_model_kcal_mol=gap,
                       class_spread={'Ca':max(ca)-min(ca),'La':max(la)-min(la)},
                       all_cross_class_pairs_ordered=ordered==total,ordered_cross_class_pairs=ordered,total_cross_class_pairs=total)
            if gap>MINIMUM_CALIBRATION_GAP:
                cal['bands']={'Ca_supported_max_R_model_kcal_mol':max(ca),'La_supported_min_R_model_kcal_mol':min(la),
                              'protocol_id':COMPOSITE_PROTOCOL,'representation':rep,'solver':solver}
                for row in rows:
                    if row['representation']!=rep or row['role']=='direct_direction_development' or row['composite_status']!='available':continue
                    value=row['composite_R_model_kcal_mol']
                    row['composite_decision']='Ca-supported' if value<=max(ca) else 'La-supported' if value>=min(la) else 'inconclusive'
                    row['composite_expected_region_pass']=row['composite_decision']==row['expected_class']+'-supported'
        calibration[rep]=cal
        by_case={r['case_id']:r for r in rows if r['representation']==rep}
        for pair in inv['comparison_definitions']['alpha_GGR_pairs']:
            a,b=by_case[pair['La_like']],by_case[pair['Ca_like']]
            ready=a['composite_status']==b['composite_status']=='available'
            pair_matrix.append(dict(pair,representation=rep,native_margin_model_kcal_mol=a['native_R_model_kcal_mol']-b['native_R_model_kcal_mol'],
                                    composite_margin_model_kcal_mol=a['composite_R_model_kcal_mol']-b['composite_R_model_kcal_mol'] if ready else None,
                                    composite_status='available' if ready else 'unavailable'))
    for c in inv['cases']:
        both={rep:next(r for r in rows if r['case_id']==c['case_id'] and r['representation']==rep) for rep in ('core','context')}
        ready=all(r['composite_status']=='available' for r in both.values())
        context_effects.append({'case_id':c['case_id'],'native_context_minus_core_R_model_kcal_mol':both['context']['native_R_model_kcal_mol']-both['core']['native_R_model_kcal_mol'],
                               'composite_context_minus_core_R_model_kcal_mol':both['context']['composite_R_model_kcal_mol']-both['core']['composite_R_model_kcal_mol'] if ready else None,
                               'status':'available' if ready else 'unavailable'})
    qualification=numerical_check(numerical_collections,lookup,sources,inv_pin)
    return {'protocol_id':COMPOSITE_PROTOCOL,'inventory':inv_pin,'collections':col_pins,'solver':solver,'implementation':record(__file__),
            'rows':rows,'calibration':calibration,'alpha_GGR_matrix':pair_matrix,'context_effects':context_effects,
            'native_PQQ_bands':inv['native_PQQ_bands'],'failed_or_incomplete_row_history':failure_history,
            'available_composite_case_representations':sum(r['composite_status']=='available' for r in rows),
            'case_representation_denominator':66,'component_audits':component_audits,
            'numerical_qualification':qualification['status'],'numerical_check':qualification,
            'interpretation':'preliminary utility separate from numerical qualification',
            'absolute_aquo_reference':None,'S_kcal_mol':None,'baseline_changed':False,'all_evidence_consumed':True,
            'alpha_GGR_independent_biological_comparisons':1,'new_scientific_calls':0}


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('inventory');q.add_argument('--root',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('compare');q.add_argument('--inventory',required=True);q.add_argument('--collection',action='append',required=True)
    q.add_argument('--solver',choices=('native','ordinary_tight'),required=True);q.add_argument('--output',required=True)
    q.add_argument('--numerical-collection',action='append',default=[])
    a=vars(p.parse_args());op=a.pop('command')
    if op=='inventory':
        result=inventory(**a);summary={k:result[k] for k in ('total_case_representations','available_native_MACE_endpoints','new_scientific_calls')}
    else:
        result=comparison(a['inventory'],a['collection'],a['solver'],a['numerical_collection']);write_new(a['output'],result)
        summary={k:result[k] for k in ('available_composite_case_representations','case_representation_denominator','new_scientific_calls')}
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
