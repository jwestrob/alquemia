"""Read-only known-reference inventory. No coordinate preparation or model calls."""
from pathlib import Path
import collections,hashlib,json
ROOT=Path(__file__).resolve().parents[3]
def pin(path):
 p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(path):return json.loads(Path(path).read_text())
def checked(p):
 assert pin(p['path'])==p
 return read(p['path'])
def fragments(items):return sorted((r['kind'],tuple(sorted(r['source'].items()))) for r in items)
def main():
 selection_path=ROOT/'workspaces/motion_envelope_20260923/run_v1/SELECTION.json'
 selection=read(selection_path);inputs=checked(selection['source_inputs']);source=checked(inputs['source_preparation'])
 prepared_path=ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json';expanded=read(prepared_path)
 own_path=ROOT/'workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/manifest.json';own=read(own_path)
 ten_path=ROOT/'workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json';ten=read(ten_path)
 canon={r['case_id'].split('__conditioned_')[0]:r for r in selection['pilot'] if r['role']=='canonical_calibration'}
 groups={g['selection_id']:g for g in selection['groups']}
 expanded_rows={(r['case_id'],r['selection_id']):r for r in expanded['cases']}
 pools={(r['source_case_id'],r['selection_id']):r for r in own['cases']}
 ten_rows={r['case_id']:r for r in ten['cases']}
 rows=[]
 for r in source['cases']:
  if r['source']['source_conditioning_metal']!='Ca':continue
  c=r['source'];cid=r['case_id'];root=c['root_case_id'];sid=canon[root]['selection_id'];pair=(cid,sid)
  assert c['expected_class'] in ('Ca','La') and c['all_evidence_consumed']
  for p in [c['source_structure'],c.get('original_source_structure',c['source_structure'])]:assert pin(p['path'])==p
  old=ten_rows[cid];same=False
  if old['status']=='prepared':same=fragments(checked(old['union'])['fragments'])==fragments(groups[sid]['fragments'])
  e=expanded_rows.get(pair);pool=pools.get(pair)
  if pool:
   assert e and e['status']=='prepared' and e['mapping_eligible']
   for p in list(e['maps'].values())+[e['state_signature']]:checked(p)
   for z in ('Ca','La'):
    for cell in pool['matrix'][z].values():
     assert pin(cell['xyz']['path'])==cell['xyz'];checked(cell['MACE'])
  rows.append({'case_id':cid,'protein':root,'biological_group':c['biological_group'],'expected_class':c['expected_class'],
   'conditioning':'Ca','source_structure':c.get('original_source_structure',c['source_structure']),
   'source_preparation_status':r['status'],'source_preparation_reason':r.get('reason'),
   'existing_target_selection_id':sid,'target_union':groups[sid]['union'],
   'exact_expanded_preparation':e is not None,'exact_saved_ownthree_pool':pool['case_id'] if pool else None,
   'existing_tenfold_status':old['status'],'tenfold_fragment_list_equal_target':same,
   'required_work':'reuse_saved_ownthree' if pool else 'source_preparation_excluded' if r['status']!='prepared' else 'expanded_preparation_then_proposals_and_components'})
 assert len(rows)==125 and len(canon)==25
 counts=collections.Counter(r['required_work'] for r in rows)
 assert dict(counts)=={'expanded_preparation_then_proposals_and_components':107,'reuse_saved_ownthree':3,'source_preparation_excluded':15}
 assert not any(r['tenfold_fragment_list_equal_target'] for r in rows)
 n=counts['expanded_preparation_then_proposals_and_components']
 result={'protocol_id':'Ca_conditioned_expanded_existing_policy_inventory_v1','status':'inventory_only_no_scientific_calls',
  'sources':{'source_preparation':inputs['source_preparation'],'expanded_selection':pin(selection_path),'expanded_preparation':pin(prepared_path),'ownthree_revaluation':pin(own_path),'tenfold_preparation':pin(ten_path)},
  'denominator':125,'biological_groups':25,'available_source_groups':len({r['protein'] for r in rows if r['source_preparation_status']=='prepared'}),
  'source_classes':dict(collections.Counter(r['expected_class'] for r in rows)),
  'supported_source_classes':dict(collections.Counter(r['expected_class'] for r in rows if r['source_preparation_status']=='prepared')),
  'counts':dict(counts),'additional_exact_saved_expanded_pools_beyond_probes':0,
  'exact_tenfold_fragment_matches':0,'existing_probes_MACE_cells':18,'existing_probes_GFN_inputs_per_backend':36,
  'conditional_upper_work':{'expanded_preparation_requests':n,'bounded_MACE_searches':2*n,'MACE_origin_evaluations':2*n,'MACE_cross_proposal_evaluations':2*n,'MACE_final_matrix_cells_including_search_endpoints':6*n,'GFN_inputs_per_backend':12*n,'primary_plus_independent_GFN_calls':24*n,'optimizer_trial_calls':'not determined by inventory; existing bounded optimizer policy, measure actual'},
  'scope_distinction':'Apply each protein existing canonical La-three-source4.3A fragment union to each Ca source, exactly as declared probes. New Ca-only or mixed union discovery changes preparation rule and is not this extension.',
  'reuse_qualification':'Existing three ownthree geometries/MACE are verified. Corrected GFN results are reusable only after root current transfer finishes and each actual raw component passes; no missing energy presumed.',
  'gates':['root interprets current corrected transfer and measured throughput','preserve15 original source PQQ exclusions','source graph/cap/state and paired endpoint checks for each new prepared source','same bounded optimizer settings and common ownthree pool, no old-tenfold proposal substitution','same actual corrected canonical bands, no Ca-fold refit'],
  'new_folds':0,'new_prepared_coordinates':0,'new_molecular_calls':0,'production_changed':False,'rows':rows}
 out=Path(__file__).with_name('INVENTORY.json');out.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('rows','sources')},indent=2))
if __name__=='__main__':main()
