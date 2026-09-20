#!/usr/bin/env python3
"""Clone existing fast-preparer roles for every archived reference sample.

Read-only structural inspection plus manifest creation. No normalization,
protonation, folding or energy execution. Ca-conditioned inputs explicitly
require the named residue-only adapter before they can pass the old preparer.
"""
import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import gemmi

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from detect_pqq import inspect_pqq_candidates


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    diag=ROOT/'diagnostics/accommodation_controls_20260920'
    ws=ROOT/'workspaces/accommodation_controls_20260920'
    ip=ws/'pqq_fold_samples_inventory.json'
    original=ROOT/'diagnostics/pqq_fast_release_20260920/SOURCES.json'
    inv=json.loads(ip.read_text()); template=json.loads(original.read_text())
    originals={c['case_id']:c for c in template['cases']}
    cases=[]; meta=[]; errors=[]
    for row in inv['rows']:
        root_case=originals[row['panel_id']]
        assert root_case['normalization']=='protenix_generic_PQQ'
        for conditioning,mode in row['modes'].items():
            for sample in mode['samples']:
                case=copy.deepcopy(root_case)
                case['case_id']=row['panel_id']+'__conditioned_'+conditioning+'__'+sample['sample_id']
                assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',case['case_id'])
                case['source_structure']={k:sample['cif'][k] for k in ('path','sha256')}
                case['role']='structural_robustness_development'
                case['root_case_id']=row['panel_id']
                case['source_conditioning_metal']=conditioning
                case['canonical_coordinate_match']=sample['matches_canonical_coordinate_fingerprint']
                case['primary_evaluation_pool']=not case['canonical_coordinate_match']
                if not case['canonical_coordinate_match']:
                    for key in list(case):
                        if key.startswith('archived_'):
                            del case[key]
                structure=gemmi.read_structure(case['source_structure']['path'])
                model=structure[0]
                metals=[(c,r,a) for c in model for r in c for a in r if a.element.name in ('Ca','La')]
                assert len(metals)==1
                chain,residue,atom=metals[0]
                assert atom.element.name==conditioning
                raw={'chain':chain.name,'resnum':residue.seqid.num,'icode':residue.seqid.icode.strip(),
                     'resname':residue.name,'atom':atom.name,'element':atom.element.name,
                     'xyz_A':list(atom.pos)}
                normalized={k:raw[k] for k in ('chain','resnum','icode','atom')}
                normalized['resname']='LA'
                # Preserve actual native atom name and element. The standard
                # protonator's LA residue alias does not make Ca into La.
                case['metal']=normalized
                case['raw_source_metal']=raw
                case['required_metal_normalization']={
                    'operation':'residue_name_only_to_LA','source':raw,'target_selector':normalized,
                    'preserve_native_element':True,'preserve_atom_name':True,'preserve_coordinates':True,
                    'existing_crystal_precedent':'write_selected_normalized_structure retains native element and atom name while aliasing residue to LA'}
                inspection_error=None
                try:
                    pqq=inspect_pqq_candidates(model)
                    assert len(pqq)==1
                    locator=pqq[0].locator
                    assert locator.chain==case['pqq']['chain'] and locator.resnum==case['pqq']['resnum']
                    pqq_status='one_complete_source_PQQ_matches_canonical_selector'
                except Exception as exc:
                    inspection_error=f'{type(exc).__name__}: {exc}'
                    pqq_status='unsupported_source_PQQ'
                    errors.append({'case_id':case['case_id'],'error':inspection_error})
                readiness=('source_PQQ_failed' if inspection_error else
                           'ready_for_existing_preparer_not_yet_protonated' if conditioning=='La' else
                           'requires_Ca_source_residue_only_adapter_before_existing_protonator')
                case['source_readiness']=readiness
                cases.append(case)
                meta.append({'case_id':case['case_id'],'root_case_id':row['panel_id'],
                             'biological_group':case['biological_group'],'expected_class':case['expected_class'],
                             'conditioning':conditioning,'seed':sample['seed'],'sample':sample['sample'],
                             'canonical_coordinate_match':case['canonical_coordinate_match'],
                             'canonical_file_hash_match':sample['matches_canonical_file_sha256'],
                             'primary_evaluation_pool':case['primary_evaluation_pool'],
                             'source':case['source_structure'],'summary_confidences':sample['summary_confidences'],
                             'confidences':sample['confidences'],'input_data':mode['input_data'],
                             'raw_source_metal':raw,'normalization_target':normalized,
                             'PQQ_readiness':pqq_status,'readiness':readiness,
                             'error':inspection_error})
    assert len(cases)==250 and len({c['case_id'] for c in cases})==250
    assert sum(c['canonical_coordinate_match'] for c in cases)==25
    assert sum(c['primary_evaluation_pool'] for c in cases)==225
    result={'schema_version':template['schema_version'],'config':template['config'],
            'cases':cases,'source_template':pin(original),'full_inventory':pin(ip),
            'scope':'All 250 archived samples. Primary pools are100 unselected La-conditioned and125 Ca-conditioned;25 canonical geometries are replay only.',
            'baseline_changed':False,'new_calls_executed':0,
            'source_readiness':{'La_raw_inputs':125,'Ca_raw_inputs_requiring_residue_alias_adapter':125,
                                'PQQ_inspection_failures':len(errors)},
            'normalization_warning':'Do not pass125 Ca-conditioned cases through the unmodified La-only AF3 normalizer. Implement the declared residue-only alias while retaining nativeCa element/CA1 atom/coordinates; never cloneLA1 by accident.',
            'errors':errors}
    target=diag/'PQQ_ALL250_SOURCES.json'
    target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    meta_path=ws/'pqq_all250_metadata.json'
    meta_path.write_text(json.dumps({'source_manifest':pin(target),'cases':meta,'counts':{'all':250,'canonical_replay':25,'primary_La':100,'primary_Ca':125},'errors':errors},indent=2,sort_keys=True)+'\n')
    print(json.dumps({'manifest':pin(target),'metadata':pin(meta_path),'cases':250,'canonical_replay':25,'primary_La':100,'primary_Ca':125,'PQQ_failures':len(errors)},indent=2))


if __name__=='__main__':
    main()
