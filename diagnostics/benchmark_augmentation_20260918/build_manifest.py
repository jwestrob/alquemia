#!/usr/bin/env python3
"""Inventory acquired real structures and release the bounded evidence extension."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import gemmi

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'workspaces/benchmark_augmentation_20260918'
AA=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def build(output):
    output.mkdir(parents=True,exist_ok=False)
    baseline_path=ROOT/'workspaces/site_classifier_20260918/features_v2/features.json'
    baseline=[r for r in json.loads(baseline_path.read_text())['rows'] if r['target']=='PQQ_functional_class']
    scoring=gemmi.AlignmentScoring();scoring.match=4;scoring.mismatch=-2;scoring.gapo=-20;scoring.gape=-1
    # Dispositions are evidence/preparation decisions made without new score inspection.
    policy={
      '8GY2':('Ca','Ca_associated_structural_functional_control','prepared_DFT_unscored','Matched whole-protein MACE unavailable: heme-bearing assembly is outside frozen preparation.'),
      '4CVB':('Ca','Ca_associated_structural_functional_control','source_ready_protocol_gate','Different donor topology: GLN220/ASP333/GLU335/TYR476; no canonical Glu/Asn anchors. Requires agreed noncanonical preparation; no inherited bands.'),
      '7WMK':(None,'Ca_PQQ_structural_candidate_functional_label_not_reaudited','source_ready_evidence_and_protocol_gate','DepA enzyme study and Ca-PQQ structure identified; no direct La/Ca comparison verified. Noncanonical ligand mapping needs declared preparation.'),
      '4MH1':(None,'Ca_PQQ_structural_candidate_functional_label_not_reaudited','source_ready_evidence_and_protocol_gate','No direct La/Ca comparison verified; source Ca coordination differs from canonical core. Exact assay/metal-state crosswalk remains required.'),
      '3DAS':(None,'Ca_PQQ_structural_candidate_functional_label_not_reaudited','source_ready_evidence_and_protocol_gate','Different sugar-dehydrogenase core includes backbone carbonyls/waters; canonical protocol cannot be inherited.'),
      '3A9H':(None,'PQQ_structure_without_local_Ca_site','source_ready_site_gate','No deposited metal within 4 A of the first PQQ; do not add or move a metal.'),
      '2D0V':(None,'homologous_Ca_PQQ_structure','same_existing_homology_group','90.2% nearest current sequence; not a new independent group. Functional assay crosswalk not reaudited.'),
      '5XM3':(None,'PQQ_Mg_source_structure','source_metal_and_evidence_gate','Deposited metal is Mg, not Ca; do not relabel crystal occupancy as a Ca measurement.'),
      '1KV9':('Ca','existing_calibration_accession_structure','duplicate_biological_group','Q8GR64 already in calibration; structural replicate only.'),
      '1YIQ':('Ca','existing_calibration_accession_structure','duplicate_biological_group','Q4W6G0 already in calibration; structural replicate only.'),
      '6DAM':('La','existing_calibration_accession_structure','duplicate_biological_group','A0A3F2YLY8 already in calibration; structural replicate only.'),
      '6ZCV':('La','existing_calibration_accession_Pr_structure','duplicate_biological_group','Q88JH0 already in calibration. Pr geometry is not a new La/Ca affinity observation.'),
      '6ZCW':('La','existing_calibration_accession_Pr_structure','duplicate_biological_group','Q88JH0 already in calibration; same protein as 6ZCV.'),
      '7O6Z':(None,'existing_calibration_accession_Nd_structure','duplicate_biological_group','94.8% to canonical A0ACD6B9F2 over full sequence; Nd activity alone does not supply a new strict La label.'),
      '9M2J':(None,'unpublished_PQQ_structure','source_ready_evidence_gate','Deposited citation says to be published; no independently verified functional La/Ca label.'),
      '9M2K':(None,'unpublished_PQQ_structure','source_ready_evidence_gate','Same protein/variant family as 9M2J; no independently verified functional La/Ca label.'),
      '9OOZ':(None,'PqqT_Y161W_PQQ_structure','source_ready_metal_site_gate','No bound La/Ca/Gd in deposition. PqqT variant does not resolve the existing metal-site mapping gate.'),
      '9X0Q':(None,'8GY2_membrane_region_variant','duplicate_biological_group','Same ADH family as 8GY2; modified construct, not independent biological evidence.'),
      '9X0R':(None,'8GY2_membrane_region_variant','duplicate_biological_group','Second form of the same variant as 9X0Q; not independent biological evidence.'),
    }
    rows=[]
    for pdb in sorted(policy):
        f=WORK/'sources'/f'{pdb}.cif';b=gemmi.cif.read(str(f)).sole_block();st=gemmi.make_structure_from_block(b)
        polymer=[]
        for c in st[0]:
            seq=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in c if r.name in AA)
            if seq:polymer.append((c.name,seq))
        chain,seq=polymer[0]
        comparisons=[]
        for old in baseline:
            a=gemmi.align_string_sequences(list(seq),list(old['sequence']),[],scoring)
            comparisons.append({'case_id':old['case_id'],'identity_relative_longer':a.match_count/max(len(seq),len(old['sequence'])),'cigar':a.cigar_str()})
        comparisons.sort(key=lambda x:x['identity_relative_longer'],reverse=True)
        sites=[]
        for c in st[0]:
            for r in c:
                if r.name!='PQQ':continue
                nearby=[]
                for cc in st[0]:
                    for rr in cc:
                        for a in rr:
                            if not a.element.is_metal:continue
                            d=min(a.pos.dist(v.pos) for v in r)
                            if d<4:nearby.append({'selector':f'{cc.name}:{rr.name}{rr.seqid}/{a.name}','element':a.element.name,'nearest_PQQ_distance_A':d})
                sites.append({'PQQ':f'{c.name}:PQQ{r.seqid}','nearby_metals':nearby})
        label,kind,status,gate=policy[pdb]
        citation=[]
        for r in b.find('_citation.',['id','title','pdbx_database_id_DOI']):
            citation.append({'id':r[0],'title':gemmi.cif.as_string(r[1]),'doi':gemmi.cif.as_string(r[2])})
        row={'case_id':'EXT_'+pdb,'pdb_id':pdb,'label':label,'evidence_stratum':kind,'status':status,'remaining_gate':gate,
             'source':pin(f),'source_url':f'https://files.rcsb.org/download/{pdb}.cif','primary_citations':citation,
             'accessions':list(b.find_values('_struct_ref.pdbx_db_accession')),'deposited_assembly_ids':list(b.find_values('_pdbx_struct_assembly.id')),
             'coordinate_model':1,'sequence_author_chain':chain,'sequence':seq,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
             'sequence_length':len(seq),'nearest_existing_sequence':comparisons[0],
             'matches_existing_group_at_50pct':comparisons[0]['identity_relative_longer']>=.5,
             'sequence_comparisons':comparisons,'PQQ_sites':sites,
             'eligible_direct_affinity':False,'eligible_PQQ_association_transfer':pdb=='8GY2',
             'broad_fold':'PQQ_eight_blade_related_or_noncanonical_not_independent_fold_claim',
             'new_energy_evaluations':0,'new_score_inspected':False,
             'score_exposure':'No scores computed or inspected here. Existing diagnostic JSON/TSV/Markdown identifier scan found no 8GY2/O05542, 4CVB/Q93RE9, 7WMK or 4MH1/E3F069 records; no exhaustive historical blindness claim.'}
        if pdb=='8GY2':
            prep=WORK/'prepared_8gy2_v1/result.json';row['preparation']=pin(prep)
            row['DFT_protocol']='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
            row['MACE_status']='unsupported_heme_assembly'
            row['label_limitations']=['Ca-associated enzyme/structure; no matched La/Ca affinity experiment or demonstrated La exclusion.','Primary linked biochemical article is identified; publisher full text returned HTTP 403 in this curation.','New group under the frozen whole-sequence rule, not proof of independent fold or a balanced additional La/Ca family panel.']
        rows.append(row)
    # Candidate-to-candidate grouping uses the same frozen criterion, no new thresholds.
    links=[]
    for i,a in enumerate(rows):
        for b in rows[:i]:
            alignment=gemmi.align_string_sequences(list(a['sequence']),list(b['sequence']),[],scoring)
            ident=alignment.match_count/max(len(a['sequence']),len(b['sequence']))
            if ident>=.5:links.append({'a':a['case_id'],'b':b['case_id'],'identity_relative_longer':ident})
    old_ledger=ROOT/'diagnostics/benchmark_set_20260915/EVIDENCE.json'
    old=json.loads(old_ledger.read_text())['rows']
    inherited=[]
    names={'H19_FAM1_ExaF','H19_Rkho_XoxF5','H19_Tcon_XoxF5','H19_Gmar_XoxF5','5GB1C_MxaF','LW13_XoxF','LW13_MxaF','PQQT_WT','PQQT_K142A','PQQT_K142D','AQUALYSIN_LOW_CA','SIXB_B2','SIXB_C5'}
    for r in old:
        if r['target_id'] in names:
            inherited.append({'case_id':r['target_id'],'direction':r['direction'],'evidence_stratum':r['label_type'],
                              'status':r['benchmark_disposition'],'mapping_status':r['mapping_status'],
                              'primary_sources':r['primary_sources'],'source_ledger':pin(old_ledger),
                              'prepared_here':False,'gate_unchanged':True})
    result={'schema':'alquemia.existing_evidence_benchmark_augmentation.v1','agreement':pin(Path(__file__).with_name('AGREEMENT.md')),
            'implementation':pin(Path(__file__)),'reference_features':pin(baseline_path),'parent_ledger':pin(old_ledger),
            'new_structure_inventory':rows,'candidate_sequence_links_at_50pct':links,'inherited_gated_candidates':inherited,
            'ready_DFT_site_count':1,'ready_DFT_endpoint_count':2,'ready_matched_DFT_MACE_count':0,'new_direct_affinity_groups':0,
            'new_energy_evaluations':0,'new_classifier_fits':0,'new_cluster_jobs':0,
            'interpretation':'Extension of source/preparation inventory; no added accuracy estimate. Ca association is distinct from La/Ca affinity. Only 8GY2 is prepared. Candidates lacking evidence or compatible chemistry remain gated.'}
    (output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    fields=['case_id','label','evidence_stratum','status','sequence_length','matches_existing_group_at_50pct','remaining_gate']
    with (output/'READINESS.tsv').open('x') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps({'manifest':pin(output/'manifest.json'),'structures':len(rows),'ready_DFT_sites':1,'ready_matched_DFT_MACE':0,'new_direct_affinity_groups':0},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();build(args.output.resolve())
