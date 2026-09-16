#!/usr/bin/env python3
"""Select the agreed existing AF3 XoxF models and freeze explicit core mappings."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import sys
import gemmi

A=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
EQ=Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue')
AF=EQ/'coordination_review/local_geometry/af3_comparison'
OUT=A/'workspaces/plm_xoxf_fixed_core_20260916'
HERE=Path(__file__).resolve().parent
APPROVAL=EQ/'xoxf_fixed_core/authorization.json'
EXPECTED={'PQQSEQ_242fa05e3ffc20087d42':('PLM0_60_coex_jun17_scaffold_38_54',[190,274,316,318]),
          'PQQSEQ_faa97386262eec4316fc':('PLM2_30_coex_sep16_scaffold_2928_2',[186,252,302,304])}
read=lambda p:json.loads(Path(p).read_text())

def record(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def write(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def wrapper():
    p=A/'diagnostics/plm_adh9_fixed_core_20260916/prepare_candidates.py'
    if record(p)['sha256']!='8d988aa7bd5e10924f879b464b58f4ad3403829e0b21b7db6e55bfefb304ec5a':raise ValueError('Frozen wrapper changed')
    return load('xoxf_frozen_preparer',p)

def select(rows):
    return sorted(rows,key=lambda r:(-r['CN'],-r['protein_La_iptm'],-r['protein_CN'],r['sample']))[0]

def main():
    if not APPROVAL.is_file():raise ValueError('Current approval missing')
    w=wrapper();fixed=w.fixed;cp=A/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/implementation_pins.json'
    original=fixed.verify_pins(cp)
    sources=read(AF/'manifest.json');target_sources={t['target_id']:t for t in sources['targets']}
    saved=EQ/'coordination_review/saved_residue_observations.tsv'
    saved_rows=list(csv.DictReader(saved.open(),delimiter='\t'))
    OUT.mkdir(exist_ok=False,parents=True);(OUT/'adapter_inputs').mkdir()
    targets=[];all_models=[];reviews=[]
    for uid,(gene,positions) in EXPECTED.items():
        source=target_sources[uid];inp=read(source['input_path']);seq=inp['sequences'][0]['protein']['sequence']
        if gene!=source['gene_id'] or hashlib.sha256(seq.encode()).hexdigest()!=source['protein_sha256']:raise ValueError('Source sequence/gene provenance mismatch')
        role_names=['anchor_glutamate','anchor_asparagine','catalytic_aspartate','extra_acidic_ligand_homolog']
        resnames=['GLU','ASN','ASP','ASP'];roles={role:{'chain':'A','resname':name,'resnum':pos,'icode':''} for role,name,pos in zip(role_names,resnames,positions)}
        mapped=[r for r in saved_rows if uid in r.values()]
        target_base={'protein_chain':'A','metal':{'chain':'B','resnum':1,'resname':'LA','atom':'LA','icode':''},'pqq':{'chain':'C','resnum':1,'resname':'PQQ','icode':''}}
        models=[]
        for sample in range(3):
            folder=AF/f'output/{uid}_AF3/seed-101_sample-{sample}'
            cif=next(folder.glob('*_model.cif'));summary=next(folder.glob('*_summary_confidences.json'))
            structure=gemmi.read_structure(str(cif));w.heavy_map(structure)
            chain=structure[0]['A'];res=list(chain);observed=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in res)
            if observed!=seq or [r.seqid.num for r in res]!=list(range(1,len(seq)+1)):raise ValueError('Raw fold sequence/numbering mismatch')
            if any(a.element.name in {'H','D'} for c in structure[0] for r in c for a in r):raise ValueError('Raw fold contains H')
            if set(c.name for c in structure[0])!={'A','B','C'}:raise ValueError('Not one protein/La/PQQ')
            s=read(summary);confidence=float(s['chain_pair_iptm'][0][1])
            site,pkey,pqq,prepared,contacts,sulfur,untyped=w.site_state(structure[0],target_base)
            directN=sum(c.element.upper()=='N' for c in contacts)
            models.append({'target_id':uid,'gene_id':gene,'sample':sample,'CN':len(contacts),'protein_CN':sum(c.residue.chain=='A' for c in contacts),
              'protein_La_iptm':confidence,'protein_PQQ_iptm':s['chain_pair_iptm'][0][2],'iptm':s['iptm'],'direct_N':directN,
              'admission_pass':len(contacts)>=7 and directN<=2 and math.isfinite(confidence) and confidence>=.9,
              'source_cif':record(cif),'summary_confidence':record(summary),'protein_sequence_sha256':source['protein_sha256'],
              'typed_contacts':[fixed.contact_metadata(c,frozenset()) for c in contacts]})
        chosen=select(models);all_models.extend(models)
        if not chosen['admission_pass']:raise ValueError('Selected model fails frozen admission gate; no substitution')
        structure=gemmi.read_structure(chosen['source_cif']['path']);chain=structure[0]['A']
        bynum={r.seqid.num:r for r in chain};cat=bynum[positions[2]]
        if seq[positions[2]-2:positions[2]]!='WD' or seq[positions[3]-1]!='D':raise ValueError('Catalytic WD/D+2 mapping mismatch')
        for role,key in roles.items():
            if bynum[key['resnum']].name!=key['resname']:raise ValueError('Role residue identity mismatch')
        candidates=[]
        for r in chain:
            if r.name not in {'ARG','LYS'}:continue
            ds=sorted((a.pos.dist(b.pos),a.name,b.name) for a in cat if a.name in {'OD1','OD2'} for b in r if b.name in fixed.CATION_HBOND_ATOMS[r.name])
            candidates.append({'resname':r.name,'resnum':r.seqid.num,'distance_A':ds[0][0],'asp_atom':ds[0][1],'partner_atom':ds[0][2]})
        candidates.sort(key=lambda x:x['distance_A']);near=[x for x in candidates if x['distance_A']<=3.5]
        if len(near)!=1:raise ValueError('Selected model lacks unique mapped catalytic cationic partner')
        partner=near[0];roles['catalytic_asp_cationic_partner']={'chain':'A','resname':partner['resname'],'resnum':partner['resnum'],'icode':''}
        case=f'{uid}_AF3_sample{chosen["sample"]}';dst=OUT/'adapter_inputs'/f'{case}.cif';score=OUT/'adapter_inputs'/f'{case}_summary.json'
        shutil.copyfile(chosen['source_cif']['path'],dst);shutil.copyfile(chosen['summary_confidence']['path'],score)
        target={'case_id':case,'target_id':uid,'gene_id':gene,'rank':chosen['sample'],'source_cif':record(dst),'summary_confidence':record(score),
          'source_provenance':{'raw_AF3_model':chosen['source_cif'],'raw_AF3_summary':chosen['summary_confidence'],'fold_input':record(source['input_path']),
              'fold_manifest':record(AF/'manifest.json'),'saved_homology_residue_observations':record(saved)},'roles':roles,**target_base}
        w.role_state(structure[0],target)
        targets.append(target);reviews.append({'target_id':uid,'gene_id':gene,'selected_sample':chosen['sample'],'selection_metrics':chosen,
          'roles':roles,'cationic_partners_by_distance':candidates,'saved_homology_rows':mapped,
          'sequence_context':{k:seq[max(0,v['resnum']-9):v['resnum']+8] for k,v in roles.items()},
          'partner_mapping':'Unique side-chain Arg/Lys contact within3.5A to mapped catalytic Asp; D+27 Arg in both homologs; no model substitution'})
    candidate=OUT/'candidate_manifest.json';write(candidate,{'schema_version':'plm.adh9.candidate_manifest.v1','protocol_id':w.PROTOCOL,'approval':record(APPROVAL),'targets':targets})
    pins={'schema_version':'plm.adh9.fixed_core_pins.v1','protocol_id':w.PROTOCOL,'wrapper':record(w.__file__),
      'candidate_manifest':record(candidate),'calibration_implementation_pins':record(cp),
      'calibration_result':record(A/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'),
      'candidate_gates':{'minimum_confidence':.9,'minimum_CN':7,'cutoff_A':3.1,'maximum_direct_N':2,'maximum_Asp_partner_A':3.5},
      'authorities':{'approval':record(APPROVAL),'selection_script':record(__file__)}}
    write(OUT/'candidate_implementation_pins.json',pins)
    write(OUT/'selection_review.json',{'status':'PASS','target_count':2,'models_checked':6,'models':all_models,'selected':reviews,
      'selection_rule':'CN descending; protein-La iPTM descending; proteinCN descending; sampleindex ascending; then gates and selected-model partner check',
      'approval':record(APPROVAL),'script':record(__file__)})
    for x in reviews:print(x['gene_id'],'sample',x['selected_sample'],'CN',x['selection_metrics']['CN'],'partner',x['cationic_partners_by_distance'][0],flush=True)

if __name__=='__main__':main()
