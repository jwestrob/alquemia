"""Exact canonical PQQ source inventory for a separately calibrated MACE scorer."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from mace_hybrid import check_atoms

SOURCE_PROTOCOL='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'


def source_scorer(calibration):
    path=Path(calibration).resolve().parent/'score.py'
    spec=importlib.util.spec_from_file_location('canonical_archive_scorer',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module,path


def endpoint_states(manifest,charges,directory):
    states={}
    for metal in ('La','Ca'):
        def resolved(ref):
            return {**ref,'path':str((Path(directory)/ref['path']).resolve())}
        xp=verify(resolved(manifest['outputs'][metal+'_xyz']));ip=verify(resolved(manifest['outputs'][metal+'_input']))
        lines=[line.strip() for line in ip.read_text().splitlines() if line.strip()]
        declarations=[line.split() for line in lines if line.lower().startswith('* xyzfile')]
        if len(declarations)!=1 or declarations[0][2:4]!=[str(charges[metal]),'1']:
            raise InvalidArtifact('frozen endpoint charge or singlet declaration differs')
        atoms=xyz(xp)
        if atoms[0][0]!=metal:raise InvalidArtifact('canonical metal index/element differs')
        states[metal]={'xyz':record(xp),'charge':charges[metal],'spin_multiplicity':1,
                       'state':check_atoms(atoms,charges[metal]),'baseline_input':record(ip)}
    paired(verify(states['La']['xyz']),verify(states['Ca']['xyz']),charges['La'],charges['Ca'])
    return states


def archived_pair(target,row,score,pins):
    mp=verify(target['manifest']);m=read_json(mp)
    if m['protocol_id']!=SOURCE_PROTOCOL:raise InvalidArtifact('source protocol changed')
    artifacts={}
    for metal in ('La','Ca'):
        task=next(t for t in m['tasks'] if t['task_id']==metal)
        value,art=score.score_task(task_id=metal,task=task,manifest_path=mp,
                                   manifest_hash=target['manifest']['sha256'],pins=pins)
        if art!=row['artifacts'][metal] or value!=row['energies_hartree'][metal]:
            raise InvalidArtifact('released DFT result does not match verified archived execution')
        execution=read_json(verify(art['execution']))
        for ref in execution['artifacts'].values():verify(ref)
        artifacts[metal]=art
    contrast=(artifacts['Ca']['energy_hartree']-artifacts['La']['energy_hartree'])*score.HA2KCAL
    if abs(contrast-row['R_kcal_mol'])>1e-8:raise InvalidArtifact('released baseline contrast changed')
    return artifacts


def target_record(case_id,label,role,group,target,baseline):
    mp=verify(target['manifest']);m=read_json(mp)
    if m['protocol_id']!=SOURCE_PROTOCOL:raise InvalidArtifact('wrong canonical source protocol')
    if (m['fixed_core']['water_policy']!='dry_exclude_all_source_and_synthetic_waters' or
            m['pqq']['microstate_id']!='pqq_ox_3minus_v1' or m['pqq']['formal_charge']!=-3 or
            not m['pqq']['present'] or m['pqq']['multiplicity']!=1 or
            any(f['kind'] not in ('fixed_core_pqq','fixed_core_protein_sidechain','fixed_core_cationic_sidechain') for f in m['qm_fragments'])):
        raise InvalidArtifact('canonical cofactor, water inventory or fragment chemistry differs')
    for key in ('heavy_coordinate_check','normalization_manifest','normalized','protonated','protonation_manifest'):
        verify(target[key])
    source=target.get('source_cif',target.get('raw_source'));verify(source)
    if not target['nonmetal_coordinates_byte_identical']:raise InvalidArtifact('source paired coordinates not certified')
    endpoints=endpoint_states(m,target['charges'],mp.parent)
    for metal,q in [('La',3),('Ca',2)]:
        if sum(f['formal_charge'] for f in m['qm_fragments'])+q!=endpoints[metal]['charge']:
            raise InvalidArtifact('fragment formal charge ledger does not close')
    return {'case_id':case_id,'expected_class':label,'evaluation_role':role,
            'evidence_stratum':'canonical_PQQ_functional_class','prospectively_blind':False,
            'sequence_accession_group':group,'homology_independence_claimed':False,
            'source_geometry':source,'source_manifest':target['manifest'],
            'preparation_provenance':{k:target[k] for k in ('heavy_coordinate_check','normalization_manifest','protonation_manifest')},
            'assembly':{'scope':'unchanged_frozen_core','selected_site':m['selected_site']},
            'fixed_core':m['fixed_core'],'cofactor':m['pqq'],'explicit_water_inventory':[],
            'source_charge_ledger':m['charge_ledger'],'source_fragments':m['qm_fragments'],
            'descriptors':{'coordination':m['coordination'],'acidic_Dplus2':target.get('acidic_Dplus2'),
                           'atom_composition':dict(Counter(a[0] for a in xyz(verify(endpoints['La']['xyz']))[1:])),
                           'charges':target['charges']},
            'endpoints':endpoints,'baseline':baseline,'preparation_status':'exact_archived_inputs_verified'}


def audit(calibration,holdouts,external_preparation,external_collection,panel,agreement,output):
    # This function reads existing scientific artifacts; it runs no model/ORCA.
    from ggr_sensitivity import executed
    score,scorer_path=source_scorer(calibration)
    cal=read_json(calibration);hold=read_json(holdouts);external=read_json(external_preparation)
    if cal['protocol_id']!=SOURCE_PROTOCOL or hold['protocol_id']!=SOURCE_PROTOCOL:
        raise InvalidArtifact('canonical released protocols required')
    prep=read_json(verify(cal['preparation']));hp=read_json(verify(hold['preparation']))
    pins=read_json(verify(cal['implementation_pins']))
    frozen=list(csv.DictReader(Path(panel).open(),delimiter='\t'))
    if any(r['metal_label'] not in ('Ln','Ca') for r in frozen):raise InvalidArtifact('unsupported frozen label')
    labels={r['candidate_id'].casefold():('La' if r['metal_label']=='Ln' else 'Ca') for r in frozen}
    if len(labels)!=25 or len(cal['scores'])!=25 or len(prep['targets'])!=25:
        raise InvalidArtifact('frozen25 calibration inventory required')
    sources={k:record(v) for k,v in [('calibration',calibration),('holdouts',holdouts),('external_preparation',external_preparation),
                                      ('external_collection',external_collection),('panel',panel),('agreement',agreement)]}
    rows=[];targets={t['panel_id']:t for t in prep['targets']}
    for row in cal['scores']:
        name=row['panel_id'];group=name.split('-pqq-')[0].casefold()
        if group not in labels or row['class']!=labels[group]:raise InvalidArtifact('frozen calibration ID/label differs')
        target=targets[name];art=archived_pair(target,row,score,pins)
        baseline={'artifacts':art,'published_R_kcal_mol':row['R_kcal_mol'],
                  'published_S_kcal_mol':row['S_aquo_gauge_kcal_mol'],'class':row['class'],'release':record(calibration)}
        rows.append(target_record(name,row['class'],'calibration',group,target,baseline))
    if len({r['sequence_accession_group'] for r in rows})!=25 or Counter(r['expected_class'] for r in rows)!=Counter({'La':11,'Ca':14}):
        raise InvalidArtifact('calibration membership/counts changed')
    ht={t['pdb_id']:t for t in hp['targets']}
    if {r['pdb_id'] for r in hold['scores']}!={'1H4I','4MAE'}:raise InvalidArtifact('exact consumed crystal pair required')
    group_evidence={}
    for row in hold['scores']:
        name=row['pdb_id'];target=ht[name];art=archived_pair(target,row,score,pins)
        source=verify(target['raw_source']);accessions=[]
        for line in source.read_text().splitlines():
            fields=line.split()
            if fields[:3]==['DBREF',name,'A'] and len(fields)>7 and fields[5]=='UNP':accessions.append(fields[6])
        if len(set(accessions))!=1:raise InvalidArtifact('unambiguous chainA UniProt group required')
        group=accessions[0].casefold();group_evidence[name]={'source':target['raw_source'],'chain':'A','accession':accessions[0],
                                                        'also_in_calibration':group in labels}
        label={'1H4I':'Ca','4MAE':'La'}[name]
        baseline={'artifacts':art,'published_R_kcal_mol':row['R_kcal_mol'],'published_S_kcal_mol':row['S_aquo_gauge_kcal_mol'],
                  'class':row['frozen_band_call'],'release':record(holdouts)}
        rows.append(target_record(name,label,'retrospective_structural_transfer',group,target,baseline))
    ex=read_json(external_collection);em,actual=executed(verify(ex['manifest']))
    er=next(r for r in ex['rows'] if r['case']=='pqq_1kb0' and r['lane']=='fixed_core')
    target=external['target'];mp=verify(target['manifest']);manifest=read_json(mp)
    if target['pdb_id']!='1KB0' or external['expected_band']!='Ca-supported' or er['status']!='complete':
        raise InvalidArtifact('completed frozen1KB0 control required')
    for metal in ('La','Ca'):
        task=next(t for t in em['tasks'] if t['task_id']==f'pqq_1kb0_fixed_core_{metal}')
        if task['source_manifest']!=target['manifest'] or task['source_xyz']!=manifest['outputs'][metal+'_xyz']:
            raise InvalidArtifact('external execution does not use exact frozen core')
        ar=actual[task['task_id']]
        if any(er['endpoints'][metal].get(k)!=v for k,v in ar.items()):
            raise InvalidArtifact('external result differs from successful execution receipt')
    rows.append(target_record('1KB0','Ca','retrospective_external_class_transfer','q46444',target,
                              {'artifacts':er['endpoints'],'published_R_kcal_mol':er['score']['R_kcal_mol'],
                               'published_S_kcal_mol':er['score']['S_kcal_mol'],'class':er['decision'],'release':record(external_collection)}))
    result={'status':'complete','schema_version':'alquemia.mace_canonical_inventory.v1','sources':sources,
            'source_protocol_id':SOURCE_PROTOCOL,'source_verifier':record(scorer_path),
            'implementation':record(__file__),'rows':rows,'sequence_group_evidence':group_evidence,
            'calibration_count':25,'transfer_count':3,'logical_endpoints_per_checkpoint':56,
            'interpretation':'functional class and retrospective structural transfer; no independent-affinity or incremental-information claim',
            'baseline_bands':cal['calibration']['released_supported_bands']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    implementation=out/'implementation';implementation.mkdir()
    preserved=implementation/'mace_canonical.py';shutil.copyfile(__file__,preserved)
    result['implementation']=record(preserved)
    write_new(out/'inventory.json',result)
    return {'status':'complete','calibration_count':25,'transfer_count':3,'inventory':record(out/'inventory.json'),
            'sequence_group_evidence':group_evidence}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('audit')
    for key in ('calibration','holdouts','external-preparation','external-collection','panel','agreement','output'):
        a.add_argument('--'+key,required=True)
    args=vars(p.parse_args());command=args.pop('command');print(json.dumps(audit(**args),indent=2))
