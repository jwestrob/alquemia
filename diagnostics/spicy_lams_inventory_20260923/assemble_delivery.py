"""Join the completed metadata inventory; inspect file availability only."""
import argparse,collections,csv,json
from pathlib import Path
from build_inventory import ROOT,PROJECT,BINDING,SENTINEL,pin,write,table,filelist

def main(inventory,output):
 src=Path(inventory);out=Path(output);out.mkdir(parents=True,exist_ok=False)
 summary=json.loads((src/'SUMMARY.json').read_text());rows=list(csv.DictReader((src/'ORTHOLOGS.tsv').open(),delimiter='\t'));structures=json.loads((src/'STRUCTURES.json').read_text())
 coverage=[]
 for r in rows:
  hits=[m for m in structures if any(c['sequence_sha256']==r['sequence_sha256']for c in m['chains'])]
  coverage.append({'ortholog_id':r['ortholog_id'],'accession':r['accession'],'sequence_sha256':r['sequence_sha256'],'reserved_panel':r['reserved_panel'],
    'exact_mature_chain_structure_files':len(hits),'exact_structure_paths':[m['file']['path']for m in hits],
    'related_construct':'Mex: workbook113aa;8FNS105aa omitsMAPTTTTK;6MI5=workbook minusMA plusHis6'if r['id6']=='621'else None,
    'full_chain_coverage':'exact_mature_chain'if hits else'related_construct_only'if r['id6']=='621'else'not_located_in_declared_inventory_roots'})
 cov=table(out,'STRUCTURE_COVERAGE.tsv',coverage)
 links=[]
 for category in ('genomes','proteomes','annotations'):
  d=PROJECT/'gtdb_core_families'/category;files=list(d.iterdir())
  links.append({'category':category,'entries':len(files),'physical_files':sum(p.is_file()and not p.is_symlink()for p in files),
   'symlinks':sum(p.is_symlink()for p in files),'broken_symlinks':sum(p.is_symlink()and not p.exists()for p in files)})
 linkpin=write(out,'GTDB_FILESYSTEM_COUNTS.json',links)
 resources=[
 ('sealed_AF2',SENTINEL/'folding_runs/af2_template_v1','full mature chains','32models/16proteins;protein-only;2ranks;EF1–3 primary protocol','unscored sentinel;outcomes sealed',SENTINEL/'STRUCTURE_PROTOCOL.md'),
 ('Mex_FEP',BINDING/'fep','full117aa protein plus one selected metal, explicit solvent','three separate EF1/EF2/EF3 arms, not a jointly3/4occupied state','historical wrong series direction/magnitude;not reusable predictive validation',BINDING/'lanm_benchmark/SESSION_NOTES.md'),
 ('Hans_Amber',ROOT/'benchmarks/hans_lanm_amber_1264_v1','archived parameterized atomistic models/trajectories','20trajectories;La/Dy EF2 and EF4 challenge','completed negative coordination-discrimination result',ROOT/'HANS_LANM_AMBER_1264_CAPABILITY_RESULT_2026-08-04.md'),
 ('Hans_outer_water',ROOT/'benchmarks/hans_lanm_dy_outer_water_v1','archived atomistic Dy model/metadynamics','four15ns walkers;fixed local coordination states','NO_CALL;not wholeprotein folding free energy',ROOT/'benchmarks/hans_lanm_dy_outer_water_v1/STATUS.md'),
 ('Hans_QMMM',ROOT/'benchmarks/hans_lanm_dy_qmmm_correction_v1','local quantum-region/protein model machinery','historical correction chain','no qualified final production transfer established',ROOT/'diagnostics/lanm_series_followup_20260923/READINESS.md'),
 ('Hans_Protenix',ROOT/'hans_lanm_protenix_benchmark','full monomer/dimer predictions','16predictions;four paired seeds perassembly','negative metal-placement/coordination result',ROOT/'hans_lanm_protenix_benchmark/scoring/HANS_LANM_CAPABILITY_VERDICT.md'),
 ('FSRD_4La',BINDING/'folding_runs/lanm_4La_multimer_test/predict_out','110aa mature chain with4La','five same-sequence predictions;not fouroccupied-site evidence','assay sequence match only;no occupancy validation',BINDING/'folding_runs/lanm_4La_multimer_test/manifest.json'),
 ('Hans_Mex_MACE_local',ROOT/'workspaces/lanm_series_followup_20260923','44–50atom local cores','EF1/2/3 separately;La/Dy physical/effective spins explicit','initial signal fails Dy-source structural transfer',ROOT/'diagnostics/lanm_series_followup_20260923/DY_TRANSFER_REPORT.md')]
 rr=[];files=[]
 for name,p,representation,scope,status,evidence in resources:
  rr.append({'resource_id':name,'path':str(p),'exists':p.exists(),'representation':representation,'scope':scope,'historical_status':status,'evidence':pin(evidence)})
  if name in ('Hans_Amber','Hans_outer_water','Hans_QMMM','Hans_Mex_MACE_local'):
   for f in filelist(p):
    if f.suffix.lower()in('.dcd','.xtc','.prmtop','.parm7','.xml','.pdb','.cif','.xyz'):
     files.append({'resource_id':name,'path':str(f),'bytes':f.stat().st_size,'type':f.suffix,'trajectory_or_energy_contents_read':False})
 resourcepin=table(out,'GLOBAL_AND_LOCAL_RESOURCES.tsv',rr);filepin=table(out,'ADDITIONAL_REPRESENTATION_FILES.tsv',files)
 h=collections.defaultdict(list)
 for m in structures:h[m['file']['sha256']].append(m['file']['path'])
 counts={'source_inventory':pin(src/'SUMMARY.json'),'source_artifacts':summary['artifacts'],
  'structure_coverage':cov,'resource_summary':resourcepin,'additional_resource_files':filepin,'GTDB_counts':linkpin,
  'exact_mature_structure_accessions':sum(int(r['exact_mature_chain_structure_files'])>0 for r in coverage),
  'related_construct_only_accessions':sum(r['full_chain_coverage']=='related_construct_only'for r in coverage),
  'no_structure_located_accessions':sum(r['full_chain_coverage']=='not_located_in_declared_inventory_roots'for r in coverage),
  'identical_structure_file_groups':[v for v in h.values()if len(v)>1],
  'additional_representation_counts':{k:dict(collections.Counter(r['type']for r in files if r['resource_id']==k))for k in sorted({r['resource_id']for r in files})},
  'workbook_ID6_min_max':[min(int(r['id6'])for r in rows),max(int(r['id6'])for r in rows)],
  'absent_ID6_values_within_0_621':sorted(set(range(622))-{int(r['id6'])for r in rows}),
  'outcomes_read':False,'new_molecular_calls':0,'new_jobs':0,'originals_modified':False}
 write(out,'DELIVERY.json',counts);print(json.dumps({k:v for k,v in counts.items()if k not in ['source_artifacts','identical_structure_file_groups']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inventory',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(a.inventory,a.output)
