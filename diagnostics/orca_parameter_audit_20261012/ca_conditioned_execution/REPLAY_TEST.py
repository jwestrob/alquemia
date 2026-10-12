"""Real archived parser/algebra checks. Never launches a scientific executable."""
from pathlib import Path
import argparse,copy,importlib.util,json

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ca_solvent_adapter',HERE/'GFN2.py')
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
ROOT=g.ROOT

def main(output):
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
 scope_path=ROOT/'workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json'
 scope=g.f.read(scope_path)
 archive_path=ROOT/'workspaces/motion_envelope_20260923/pilot34_pool_v1/COLLECTION_FINAL.json'
 archive=g.f.read(archive_path)
 selected={r['case_id']:r for r in scope['rows'] if r['execution_status']=='reuse_actual_corrected_pool'}
 cases=[]
 for old in archive['cases']:
  if old['case_id'] not in selected:continue
  assert old['pool']['status']=='available'
  source=old['source']['origin_row']
  case={'case_id':old['case_id'],'source_case_id':source['source_case_id'],'selection_id':source['selection_id'],
        'source':old['source'],'candidates':{q['id']:{**q,'status':'available'} for q in old['candidates']},
        'native_pool_status':'available','matrix':{z:{} for z in ('Ca','La')}}
  for z in ('Ca','La'):
   for q in g.ORDER:
    cell=old['matrix'][z][q];receipt=g.f.read(cell['MACE'])
    state=receipt.get('state',receipt.get('input_state_check'))
    case['matrix'][z][q]={k:cell[k] for k in ('xyz','MACE','status')}
    case['matrix'][z][q].update(energy_eV=receipt['energy_eV'],charge=int(state['charge']),multiplicity=int(state['spin_multiplicity']))
  cases.append(case)
 assert len(cases)==3
 normalized={'evaluation_kind':'contract_normalization_of_three_real_archived_pools_not_new_science',
             'archive':g.f.pin(archive_path),'manifest':archive['manifest'],'terminal':True,'cases':cases}
 fixture=out/'ARCHIVE_CONTRACT.json';g.f.write(fixture,normalized)
 full=out/'full';summary=g.prepare(scope_path,fixture,full,fixture=True)
 assert summary['unique_inputs_per_backend']==36 and summary['actual_reused_inputs_per_backend']==36 and summary['fresh_inputs_per_backend']==0
 manifest=g.f.read(full/'manifest.json')
 for source in selected.values():assert g.verify_prior(source,manifest)==source['prior_row']
 result=g.collect(full,'COLLECTION.json');residuals=[]
 for row in result['rows']:
  expected=selected[row['case_id']]['prior_row']
  assert row['decisions']==expected['decisions']
  for policy in ('operational','mathematical','origin'):
   residual=row['pool'][policy]['R_model_kcal_mol']-expected['pool'][policy]['R_model_kcal_mol']
   assert residual==0.;residuals.append(residual)
 assert result['cost']['new_primary_calls']==result['cost']['new_independent_calls']==0
 broken=copy.deepcopy(normalized)
 broken['evaluation_kind']='deliberately_missing_cell_copy_of_real_fixture_parser_test_only'
 broken['cases'][0]['matrix']['Ca']['adaptive_La']={'status':'unavailable','reason':'deliberately_missing_parser_fixture'}
 broken['cases'][0]['native_pool_status']='unavailable'
 bp=out/'MISSING_CELL_CONTRACT.json';g.f.write(bp,broken)
 bd=out/'missing_cell';g.prepare(scope_path,bp,bd,fixture=True);missing=g.collect(bd,'COLLECTION.json')
 tested=next(row for row in missing['rows'] if row['case_id']==broken['cases'][0]['case_id'])
 assert tested['pool']['operational'] is None and tested['pool']['mathematical'] is None
 assert tested['pool']['origin']==selected[tested['case_id']]['prior_row']['pool']['origin']
 invalid=copy.deepcopy(cases[0]['matrix']['Ca']['origin']);invalid['charge']+=1
 try:g.verified_cell(invalid,'Ca',scope)
 except AssertionError:pass
 else:raise AssertionError('mismatched real-fixture state was accepted')
 def forbidden(*args,**kwargs):raise RuntimeError('scientific execution was invoked by parser test')
 original=g.f.execute_point;g.f.execute_point=forbidden
 try:
  try:g.execute(full,14)
  except AssertionError as exc:assert 'parser replay cannot launch science' in str(exc)
  else:raise AssertionError('fixture execution was not disabled')
 finally:g.f.execute_point=original
 report={'status':'passed','driver':g.f.pin(HERE/'GFN2.py'),'test':g.f.pin(__file__),
         'fixture':g.f.pin(fixture),'replayed_actual_pools':3,'reparsed_actual_inputs_per_backend':36,
         'score_residuals_kcal_mol':residuals,'checks':['three_actual_pools_exact_scores_and_decisions',
         'three_reused_pool_rows_and_raw_outputs_verified','missing_required_cell_preserves_origin_only',
         'wrong_state_rejected','parser_fixture_execution_disabled'],
         'scientific_calls':0,'new_biological_evidence':False,
         'complete_replay':g.f.pin(full/'COLLECTION.json'),'missing_cell_replay':g.f.pin(bd/'COLLECTION.json')}
 g.f.write(out/'RESULT.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True)
 main(p.parse_args().output)
