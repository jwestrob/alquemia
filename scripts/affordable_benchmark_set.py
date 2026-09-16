"""Inventory real benchmark structures and assemble evidence without scoring them."""
from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import json
import re
import shutil
from affordable_common import read_json,record,verify,write_new,InvalidArtifact


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def structure_inventory(manifest,output):
    import gemmi
    if output.exists():raise InvalidArtifact('inventory exists; choose a new version')
    source=read_json(manifest);records=[];fasta=[]
    for item in source['records']:
        path=verify(item);sid=item['pdb_id'];st=gemmi.read_structure(str(path))
        if not st:raise InvalidArtifact(f'no model: {sid}')
        block=gemmi.cif.read_file(str(path)).sole_block()
        polymers=[]
        for row in block.find(['_entity_poly.entity_id','_entity_poly.pdbx_seq_one_letter_code_can','_entity_poly.pdbx_strand_id']):
            seq=re.sub(r'\s+','',gemmi.cif.as_string(row[1]))
            chains=gemmi.cif.as_string(row[2]).split(',')
            polymers.append({'entity_id':row[0],'author_chains':chains,'sequence':seq,'sequence_sha256':sequence_hash(seq),'sequence_length':len(seq),'nonstandard_sequence_tokens':bool(re.search('[^A-Z]',seq))})
            fasta.append(f'>{sid}|entity={row[0]}|chains={",".join(chains)}|sequence_sha256={sequence_hash(seq)}\n{seq}\n')
        atoms=[(c,r,a) for c in st[0] for r in c for a in r]
        waters=sum(r.name in ('HOH','WAT','DOD') for c in st[0] for r in c)
        hetero=dict(Counter(r.name for c in st[0] for r in c if r.het_flag=='H' and r.name not in ('HOH','WAT','DOD')))
        def selector(c,r,a):
            return {'chain':c.name,'resname':r.name,'resnum':r.seqid.num,'icode':r.seqid.icode.strip(),'atom':a.name,'altloc':a.altloc.strip('\x00 ').strip(),'element':a.element.name,'occupancy':a.occ,'xyz_A':[a.pos.x,a.pos.y,a.pos.z]}
        sites=[]
        for c,r,a in atoms:
            if not a.element.is_metal:continue
            neighbors=[]
            for cc,rr,aa in atoms:
                if aa.element.name not in ('O','N','S'):continue
                dist=a.pos.dist(aa.pos)
                if dist<=3.3:
                    neighbors.append(dict(selector(cc,rr,aa),distance_A=dist,is_water=rr.name in ('HOH','WAT','DOD'),qualifies_at_3p1_A=dist<=3.1))
            sites.append({'metal':selector(c,r,a),'raw_neighbors_3p3_A':sorted(neighbors,key=lambda x:x['distance_A']),
                          'interpretation':'geometric inventory only; alternate conformers unselected; no affinity-site assignment inferred'})
        records.append({'pdb_id':sid,'source':record(path),'model_count':len(st),'inventoried_model':st[0].num,'assembly_policy':'deposited_coordinates_inventory_not_biological_assembly_selection','resolution_A':st.resolution or None,'experimental_methods':list(block.find_values('_exptl.method')),'polymers':polymers,'hetero_residues_first_model':hetero,'water_count_first_model':waters,'metal_sites_first_model':sites})
    groups={}
    for row in records:
        for poly in row['polymers']:groups.setdefault(poly['sequence_sha256'],[]).append({'pdb_id':row['pdb_id'],'entity_id':poly['entity_id'],'chains':poly['author_chains']})
    payload={'schema_version':'alquemia.benchmark_structure_inventory.v1','source_manifest':record(manifest),'implementation':record(__file__),'gemmi_version':gemmi.__version__,'records':records,'exact_sequence_replicate_groups':[v for v in groups.values() if len(v)>1]}
    write_new(output,payload)
    with output.with_suffix('.fasta').open('x') as f:f.write(''.join(fasta))
    return payload


def evidence_rows(directory):
    """Read only explicit evidence_rows.json packages, never arbitrary source JSON."""
    rows=[];packages=[]
    for path in sorted(directory.glob('*/evidence_rows.json')):
        value=read_json(path);batch=value if isinstance(value,list) else value['rows']
        packages.append(record(path))
        for row in batch:rows.append(dict(row,evidence_package=record(path)))
    return rows,packages


def historical_rows(root,output):
    """Carry already inspected controls as historical records, never new tests."""
    base=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
    release=read_json(base/'result.json');prep=read_json(verify(release['preparation']))
    targets={t['panel_id']:t for t in prep['targets']};rows=[]
    for score in release['scores']:
        tid=score['panel_id'];target=targets[tid]
        for endpoint in score['artifacts'].values():
            for key in ('input','xyz','output','execution'):verify(endpoint[key])
        rows.append({'target_id':'CALIBRATION_'+tid,'protein_name':tid,'direction':score['class'],
            'label_type':'historical_PQQ_functional_calibration','biological_group':tid.split('-pqq-')[0],
            'fold_group':'PQQ_eight_blade_ADH','primary_sources':[{'kind':'frozen_historical_label_record','record':release['preregistration']}],
            'measurements':[],'construct':'see pinned original panel and preparation','conditions':'see pinned original panel',
            'structure_ids':[],'structure_paths':[target['source_cif']],
            'usable_for':['historical_calibration_regression'],'mapping_status':'existing_frozen_preparation',
            'limitations':['Already used for calibration; not independent new evaluation.','Experimental labels inherited from frozen source record, not reaudited in this construction phase.','Motif, charge and composition already separate this panel.'],
            'score_exposure':'calibration_scores_previously_inspected','existing_scores':[dict(score,protocol_id=release['protocol_id'],release=record(base/'result.json'))],
            'prepared_manifest':target['manifest']})
    hp=base/'reserved_crystal_holdout/result/holdout_result.json';holdout=read_json(hp)
    for score in holdout['scores']:
        sid=score['pdb_id']
        rows.append({'target_id':'PQQ_CRYSTAL_'+sid,'protein_name':sid,'direction':'Ca' if score['expected_band']=='Ca-supported' else 'La',
            'label_type':'historical_PQQ_class_transfer','biological_group':sid,'fold_group':'PQQ_eight_blade_ADH',
            'primary_sources':[{'kind':'frozen_transfer_record','record':record(hp)}],'measurements':[],
            'structure_ids':[sid],'usable_for':['historical_transfer_regression'],'mapping_status':'existing_frozen_preparation',
            'limitations':['Class-transfer result, not a direct affinity label.','Sequence overlap with calibration must stay grouped; not a new biological-family test.'],
            'score_exposure':'transfer_scores_previously_inspected','existing_scores':[dict(score,protocol_id=holdout['protocol_id'],release=record(hp))]})
    current=root/'diagnostics/baseline_benchmark_20260915/RESULT.json';data=read_json(current)
    score=next(s for s in data['scores'] if s['case']=='pqq_1kb0')
    rows.append({'target_id':'PQQ_CRYSTAL_1KB0','protein_name':'Q46444 quinoprotein alcohol dehydrogenase','direction':'Ca',
        'label_type':'historical_PQQ_class_transfer','biological_group':'Q46444','fold_group':'PQQ_eight_blade_ADH',
        'primary_sources':[{'kind':'frozen_transfer_record','record':record(root/'diagnostics/pqq_q46444_1kb0_external_validation_20260915/EXPERIMENT.md')}],
        'measurements':[],'structure_ids':['1KB0'],'usable_for':['historical_transfer_regression'],'mapping_status':'existing_frozen_preparation',
        'limitations':['Functional/structural Ca class, not direct La/Ca affinity.','Score already inspected.'],
        'score_exposure':'transfer_scores_previously_inspected','existing_scores':[score]})
    rows.append({'target_id':'PQQ_CRYSTAL_6OC6','protein_name':'6OC6 secondary crystal','direction':'unresolved',
        'label_type':'historical_secondary_not_run','biological_group':'C5B120','fold_group':'PQQ_eight_blade_ADH',
        'primary_sources':[{'kind':'frozen_transfer_record','record':record(hp)}],'measurements':[],
        'structure_ids':['6OC6'],'usable_for':['structural_replicate_only'],'mapping_status':'secondary_not_run',
        'limitations':['Sequence already represented by C5B120 in calibration.','Not part of primary transfer verdict; no success or failure inferred.'],
        'score_exposure':'not_run','existing_scores':[]})
    write_new(output,rows);return {'rows':rows}


def preserve_implementation(output):
    from affordable_common import snapshot_implementation
    inventory=snapshot_implementation(output)
    records=read_json(verify(inventory))
    source=Path(__file__).resolve().parent/'protonate_cif.py'
    destination=output/source.name;shutil.copyfile(source,destination)
    # Preserved bytes remain authoritative if the live checkout changes later.
    copies={name:item['preserved_copy'] for name,item in records.items()}
    copies[source.name]=record(destination)
    write_new(output/'preserved.json',copies)
    return record(output/'preserved.json')


def prepare_generic(config,root,output,protonation_report=None):
    from protonate_cif import protonate
    from carve_generic import carve
    from affordable_peptide import repair
    cfg=read_json(config);source=verify(cfg['source']);topology=verify(cfg['topology'])
    if cfg['protocol_id']!='generic_peptide_amide_vertical_native_r2scan3c_v3':raise InvalidArtifact('unsupported preparation protocol')
    if output.exists():raise InvalidArtifact('refusing existing preparation')
    if not output.resolve().is_relative_to((root/'workspaces').resolve()):raise InvalidArtifact('workspace required')
    output.mkdir(parents=True);protonated=output/'source_protonated.pdb'
    snapshot=preserve_implementation(output/'implementation')
    result={'config':record(config),'source':record(source),'implementation_snapshot':snapshot,'sites':[], 'new_high_level_endpoint_evaluations':0}
    try:
        if protonation_report:
            prior=read_json(protonation_report)
            verify(prior['config']);verify(prior['source'])
            if prior['config']!=record(config) or prior['source']!=record(source):raise InvalidArtifact('protonation replay configuration mismatch')
            p=prior['protonation']
            if p['ph']!=cfg['protonation_pH'] or p['add_missing_residues']:raise InvalidArtifact('protonation replay policy mismatch')
            shutil.copyfile(verify(p['output']),protonated)
            result['protonation']=dict(p,output=record(protonated))
            result['protonation_reuse']={'report':record(protonation_report),'scope':'exact archived protonated coordinates; no new hydrogen placement; original staging builder bytes were not retained'}
        else:
            result['protonation']=protonate(source,protonated,ph=cfg['protonation_pH'],add_missing_residues=False)
        for site in cfg['sites']:
            row={'site':site,'status':'attempted'};stem=cfg['target_id']+'_'+site['site_id'];d=output/site['site_id']
            try:
                carve(protonated,d/'historical_v2_intermediate',stem,site_chain=site['chain'],site_resnum=site['resnum'],site_icode=site['icode'],site_atom=site['atom'],qm_inclusion_cut=3.3)
                mp=d/'historical_v2_intermediate'/f'{stem}_carve_manifest.json';row['carve_manifest']=record(mp)
                repaired=repair(mp,d/'amide_v3',topology)
                row.update(status='prepared',repair_manifest=record(d/'amide_v3/repair_manifest.json'),atom_count=repaired['paired_invariants']['atom_count'])
            except Exception as exc:row.update(status='unsupported',reason=str(exc))
            result['sites'].append(row)
    except Exception as exc:result.update(status='preparation_failed',reason=str(exc))
    write_new(output/'preparation_report.json',result);return result


def task_manifest(reports,root,output,agreement,evidence):
    import gemmi
    from affordable_common import paired,xyz,cache_key
    from affordable_benchmark import audit_reference
    root=root.resolve();output=output.resolve()
    if output.exists():raise InvalidArtifact('refusing existing task bundle')
    if not output.resolve().is_relative_to((root/'workspaces').resolve()):raise InvalidArtifact('workspace required')
    snapshot=preserve_implementation(output/'implementation')
    evidence_records,packages=evidence_rows(evidence)
    evidence_lookup={row['target_id']:row for row in evidence_records}
    tasks=[];cases=[];audits=[]
    for report_path in reports:
        report=read_json(report_path);source=verify(report['source']);st=gemmi.read_structure(str(source))
        raw={}
        for chain in st[0]:
            for residue in chain:
                for atom in residue:
                    if atom.element.is_hydrogen:continue
                    key=(chain.name,residue.seqid.num,residue.seqid.icode.strip(),atom.name.strip(),atom.altloc.strip('\x00 '))
                    if key in raw:raise InvalidArtifact('ambiguous original source atom')
                    raw[key]=atom.pos
        cfg=read_json(verify(report['config']));target=cfg.get('evidence_target_id',cfg['target_id'])
        if target=='HANS_LANM_WT':target='hans_lanm_wt'
        if target not in evidence_lookup:raise InvalidArtifact('prepared target missing from evidence: '+target)
        for site in report['sites']:
            if site['status']!='prepared':
                audits.append({'target_id':target,'site':site['site'],'status':site['status'],'reason':site.get('reason')});continue
            mp=verify(site['repair_manifest']);m=read_json(mp);positions=xyz(verify(m['outputs']['La']['xyz']));displacements=[]
            for atom in m['atom_graph']['source_to_qm']:
                if atom['kind']!='source' or atom['source']['element'] in ('H','D'):continue
                a=atom['source'];key=(a['chain'],a['resnum'],a['insertion_code'],a['atom'],a['altloc'])
                if key not in raw:raise InvalidArtifact(f'unmapped source heavy atom {key}')
                p=positions[atom['qm_index']];q=raw[key]
                displacements.append(((p[1]-q.x)**2+(p[2]-q.y)**2+(p[3]-q.z)**2)**0.5)
            maximum=max(displacements,default=0.)
            if maximum>0.002:raise InvalidArtifact('source coordinates changed during protonation/preparation')
            name=report_path.parent.name+'_'+site['site']['site_id'];ids=[];copies={}
            for metal in ('La','Ca'):
                endpoint=m['outputs'][metal];ip=verify(endpoint['input']);xp=verify(endpoint['xyz']);d=output/name/metal;d.mkdir(parents=True)
                i=d/ip.name;x=d/xp.name;shutil.copyfile(ip,i);shutil.copyfile(xp,x);copies[metal]=x;tid=name+'_'+metal;ids.append(tid)
                tasks.append({'task_id':tid,'case':name,'lane':'amide_v3','metal':metal,'charge':endpoint['charge'],'input':record(i),'xyz':record(x),'output_path':str(d/(ip.stem+'.out')),'source_input':record(ip),'source_xyz':record(xp),'source_manifest':record(mp),'protocol_id':m['protocol_id']})
            row={'case':name,'target_id':target,'lane':'amide_v3','protocol_id':m['protocol_id'],'source_manifest':record(mp),'preparation_config':report['config'],'tasks':ids,'biological_group':target,'evidence_stratum':'protein_level_or_condition_limited_support','evaluation_role':'unscored_frozen_supporting_case','paired_invariants':paired(copies['La'],copies['Ca'],m['outputs']['La']['charge'],m['outputs']['Ca']['charge']),'source_heavy_max_displacement_A':maximum,'explicit_water_inventory':m['explicit_water_inventory'],'descriptors':{'coordination':m['coordination'],'charge_La':m['outputs']['La']['charge'],'charge_Ca':m['outputs']['Ca']['charge']}}
            row.update(biological_group=evidence_lookup[target]['biological_group'],observation_group=target,evidence_stratum=evidence_lookup[target]['label_type'])
            row['cache_key']=cache_key(row);cases.append(row);audits.append({'case':name,'status':'pass','source_heavy_max_displacement_A':maximum,'source_heavy_atom_count':len(displacements),'paired':row['paired_invariants']})
    if not tasks:raise InvalidArtifact('no prepared sites')
    payload={'schema_version':'alquemia.benchmark_ready_tasks.v1','status':'prepared_not_executed','cases':cases,'tasks':tasks,'preparation_reports':[record(p) for p in reports],'preparation_audits':audits,'agreement':record(agreement),'reference':audit_reference(root),'release':record(root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'),'implementation_snapshot':snapshot,'execution_policy':{n:record(root/'scripts'/f) for n,f in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},'orca':record('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca'),'new_high_level_endpoint_evaluations':0,'decision_policy':'raw contrasts/common gauge only; no PQQ bands, universal zero, individual-site labels or threshold fitting','group_policy':'Hans EF1-EF2-EF3 ordered vector; alpha-lactalbumin both source-conditioned geometries retained with different water inventories'}
    payload['evidence_packages']=packages
    write_new(output/'manifest.json',payload);return payload


def assemble(evidence,inventory,output,completed_ledger=None,prepared_tasks=None):
    if output.exists():raise InvalidArtifact('refusing existing benchmark release')
    rows,packages=evidence_rows(evidence)
    if not rows:raise InvalidArtifact('no evidence packages')
    structures=read_json(inventory);lookup={r['pdb_id']:r for r in structures['records']}
    for structure in structures['records']:verify(structure['source'])
    def check_pins(value):
        if isinstance(value,dict):
            expected=value.get('sha256',value.get('file_sha256'))
            if isinstance(value.get('path'),str) and isinstance(expected,str):verify({'path':value['path'],'sha256':expected})
            for child in value.values():check_pins(child)
        elif isinstance(value,list):
            for child in value:check_pins(child)
    completed=read_json(completed_ledger)['scores'] if completed_ledger else []
    prepared=read_json(prepared_tasks)['cases'] if prepared_tasks else []
    score_case_map={'GGR_1GLG':['ggr_1glg_GGR'],'AEQUORIN_1SL8_VECTOR':['aequorin_1sl8_EF1','aequorin_1sl8_EF3','aequorin_1sl8_EF4'],
                    'CARP_PV_CD':['carp_parvalbumin_4cpv_CD'],'CARP_PV_EF':['carp_parvalbumin_4cpv_EF']}
    seen=set()
    for row in rows:
        if row['target_id'] in seen:raise InvalidArtifact('duplicate target_id '+row['target_id'])
        seen.add(row['target_id'])
        check_pins(row)
        local_pins=[]
        for source in row['primary_sources']:
            if isinstance(source.get('local_capture'),str):local_pins.append(record(source['local_capture']))
        row['primary_local_capture_pins']=local_pins
        if row['target_id'] in score_case_map and completed_ledger:
            row['existing_scores']=[r for r in completed if r['case'] in score_case_map[row['target_id']]]
            if not row['existing_scores']:raise InvalidArtifact('missing completed score join')
            row['completed_release']=record(completed_ledger)
            row['score_exposure']='original_v2_and_repaired_v3_scores_previously_inspected'
        ready=[c for c in prepared if c['target_id']==row['target_id']]
        if ready:row['prepared_site_cases']=ready;row['prepared_tasks_manifest']=record(prepared_tasks)
        for key in ('label_type','direction','biological_group','primary_sources','mapping_status','limitations'):
            if key not in row:raise InvalidArtifact(f"{row['target_id']} missing {key}")
        ids=row.get('structure_ids',[])
        row['available_structure_records']=[{'pdb_id':sid,'source':lookup[sid]['source'],'metal_selectors':[s['metal'] for s in lookup[sid]['metal_sites_first_model']]} for sid in ids if sid in lookup]
        row['missing_inventory_structure_ids']=[sid for sid in ids if sid not in lookup]
        row.setdefault('score_exposure',row.get('score_visibility','prior_exposure_not_established'))
        row['new_calculation_performed']=False
        if row.get('existing_scores'):disposition='reuse_existing_scores'
        elif ready:disposition='prepared_unscored'
        elif row['direction']=='unresolved':disposition='unlabeled_or_unresolved_challenge'
        elif row['available_structure_records']:disposition='structure_available_preparation_or_interpretation_gate'
        elif row.get('structure_paths'):disposition='model_available_holo_preparation_gate'
        elif row.get('sequence'):disposition='sequence_available_structure_or_identity_gate'
        else:disposition='evidence_available_input_gate'
        row['benchmark_disposition']=disposition
        row['family_holdout_group']=row.get('fold_group',row['biological_group'])
        row['family_grouping_basis']='conservative existing fold/family annotation; no sequence-identity clustering or independence claim'
        check_pins(row)
    missing={c['target_id'] for c in prepared}-seen
    if missing:raise InvalidArtifact('prepared targets absent from evidence: '+str(sorted(missing)))
    payload={'schema_version':'alquemia.existing_evidence_benchmark.v1','status':'curated_not_newly_scored','implementation':record(__file__),'evidence_packages':packages,'structure_inventory':record(inventory),'record_count':len(rows),'biological_groups':sorted({r['biological_group'] for r in rows}),'label_type_counts':dict(Counter(r['label_type'] for r in rows)),'rows':rows,'interpretation':'Counts are evidence records, not independent biological observations. Unsupported mappings remain in the ledger. No pooled accuracy or new energy is inferred.'}
    output.mkdir(parents=True)
    write_new(output/'benchmark_manifest.json',payload)
    with (output/'benchmark_manifest.tsv').open('x') as f:
        keys=['target_id','protein_name','label_type','direction','biological_group','family_holdout_group','benchmark_disposition','usable_for','mapping_status','structure_ids','score_exposure','score_visibility']
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(row.get(k),ensure_ascii=False) if isinstance(row.get(k),(dict,list)) else row.get(k,'') for k in keys})
    return payload


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('inventory');q.add_argument('--manifest',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=s.add_parser('assemble');q.add_argument('--evidence',type=Path,required=True);q.add_argument('--inventory',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--completed-ledger',type=Path);q.add_argument('--prepared-tasks',type=Path)
    q=s.add_parser('history');q.add_argument('--root',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=s.add_parser('prepare-generic');q.add_argument('--root',type=Path,required=True);q.add_argument('--config',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--protonation-report',type=Path)
    q=s.add_parser('tasks');q.add_argument('--root',type=Path,required=True);q.add_argument('--reports',type=Path,nargs='+',required=True);q.add_argument('--agreement',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--evidence',type=Path,required=True)
    a=p.parse_args()
    if a.op=='inventory':r=structure_inventory(a.manifest,a.output)
    elif a.op=='history':r=historical_rows(a.root,a.output)
    elif a.op=='prepare-generic':r=prepare_generic(a.config,a.root,a.output,a.protonation_report)
    elif a.op=='tasks':r=task_manifest(a.reports,a.root,a.output,a.agreement,a.evidence)
    else:r=assemble(a.evidence,a.inventory,a.output,a.completed_ledger,a.prepared_tasks)
    print(json.dumps({'records':len(r.get('records',r.get('rows',[]))),'tasks':len(r.get('tasks',[])),'prepared_sites':len(r.get('sites',[]))}))
if __name__=='__main__':main()
