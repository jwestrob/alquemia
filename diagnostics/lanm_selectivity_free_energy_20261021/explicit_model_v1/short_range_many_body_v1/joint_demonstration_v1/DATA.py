"""Saved-data assembly for one joint demonstration. No dispatch, fit, or tensor solve.

load(inputs, donor_endpoint_exports=None) returns exact native fixtures and raw
method/channel-labelled targets. Call from the declared Slurm demonstration.
Optional endpoint_exports maps the eight La endpoint IDs to complete fixture pins.
No loss weights, accuracy tolerances or fitted coefficients are selected here.
"""
import hashlib,json,math,os
from pathlib import Path

DONOR='matched_dry_acetate_mono_central_La'
ENDPOINTS=tuple(f'pbe0__La__{direction}__{shift}__complex' for direction in ('radial','rotation') for shift in ('minus_010','minus_005','plus_005','plus_010'))

def pin(path):
 p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
 assert pin(p['path'])==p,('changed saved input',p);return json.loads(Path(p['path']).read_text())
def target(terms,refs,units):
 """Electronic common-parent-basis target; never insert D4 into water cycles."""
 return {'energy_kcal_mol':math.fsum(c*refs[k]['energy_electronic_hartree'] for k,c in terms.items())*units['hartree_kcal_mol'],
         'mu_e_bohr':[math.fsum(c*refs[k]['dipole_total_origin0_au'][j] for k,c in terms.items()) for j in range(3)],
         'method':'common-parent-basis electronic ECP PBE0; D4 excluded','origin_A':[0.,0.,0.],
         'moment_channel':'signed energy-conjugate total moment; native physical-u moment remains separate',
         'reference_terms':terms,'consumed_development':True,'physical_accuracy_qualified':False}

def load(inputs,donor_endpoint_exports=None):
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('saved target assembly is part of the declared Slurm campaign')
 if isinstance(inputs,(str,Path)):inputs=json.loads(Path(inputs).read_text())
 dataset=read(inputs['water_dataset']);index=read(dataset['fixtures']);water={id:read(p) for id,p in index['fixtures'].items()};assert len(water)==54
 static=read(dataset['static']);refs=read(dataset['references']);refjoins=read(dataset['reference_model_joins']);cycles=read(dataset['cycles']);fields=read(dataset['fields']);forces=read(dataset['forces']);base=read(dataset['baseline']);units=base['units']
 shared=read(inputs['shared_baseline']);union=read(inputs['water_union']);donor_union=read(inputs['donor_union']);donor=read(inputs['donor_center']);displaced=read(inputs['donor_displacements']);dy=read(inputs['donor_Dy_report_only'])
 assert shared['old_h_H_zero'] and union['shared_baseline']==inputs['shared_baseline'] and donor_union['shared_baseline']==inputs['shared_baseline'] and displaced['baseline']==inputs['shared_baseline']
 for id,f in water.items():assert f['compatible'] and f['native_operator'] and union['checks'][id]['same_applicable_pair_records']
 assert donor['compatible'] and donor_union['checks']['donor/'+DONOR]['compatible'];fixtures={**water,DONOR:donor};fixture_pins={**index['fixtures'],DONOR:inputs['donor_center']}
 dyid='matched_dry_acetate_mono_central_Dy';assert dy['compatible'] and donor_union['checks']['donor/'+dyid]['compatible'];fixtures[dyid]=dy;fixture_pins[dyid]=inputs['donor_Dy_report_only']
 # Full signed moments use individual raw references, never magnitudes or paired-only residuals.
 moments=[]
 for id,r in refs.items():
  join=refjoins[id];mid=join['fixture_id'];assert mid==r['model_fixture_id'] and join['real_geometry_exact'] and mid in fixtures|static
  f=fixtures.get(mid,static.get(mid));sensitive=f['geometry']['metal_atom'] is not None and bool(f['geometry']['ligands'])
  moments.append({'id':id,'model_id':mid,'target_mu_e_bohr':r['dipole_total_origin0_au'],'origin_A':r['dipole_origin_A'],'method':r['energy_recipe'],'parameter_sensitive':sensitive,'source':dataset['references'],'source_record_pointer':'/'+id,'raw_record':r['raw_moment_record'],'moment_channel':'conjugate','consumed_development':True})
 pair=[]
 for c in cycles['pair']:
  assert all(k in fixtures|static for k in c['model_terms']);row={'id':c['id'],'model_terms':c['model_terms'],**target(c['reference_terms'],refs,units),'stored_signed_mu_Debye':c['reference_pair_mu_Debye'],'print_bound_Debye_per_component':c['retained_reference_print_bound_Debye_per_component'],'source':dataset['cycles']};pair.append(row)
 trimer=[]
 for id,c in cycles['trimer'].items():
  assert all(k in fixtures|static for k in c['model_terms']);t=target(c['reference_terms'],refs,units)
  # Keep the existing admitted energy scalar as the primary target; raw terms remain available.
  trimer.append({'id':id,'model_terms':c['model_terms'],**t,'energy_from_raw_terms_kcal_mol':t['energy_kcal_mol'],'energy_kcal_mol':c['reference_electronic_E3_kcal_mol'],'stored_signed_mu_Debye':c['reference_mu3_Debye'],'print_bound_Debye_per_component':c['reference_print_bound_Debye_per_component'],'selection':c['selection'],'source':dataset['cycles']})
 assert len(pair)==34 and len(trimer)==8
 chi=[];cs=cycles['conditional_susceptibility']
 for name,rid in cs['state_ids'].items():
  estimates=[x for x in cs['reported_rows'] if x['quantity']==name and x['kind']=='state'];chosen=[x for x in estimates if x['step']=='half' and x['estimator']=='energy'];assert len(chosen)==1
  chi.append({'id':'water/'+name,'model_id':refjoins[rid]['fixture_id'],'direction_lab':[0.,0.,1.],'target_source_au':chosen[0]['reference_a0_cubed'],'reference_print_bound':chosen[0]['reference_print_bound_a0_cubed'],'selected_estimator':'energy curvature at saved half field','all_saved_estimators':estimates,'source':dataset['cycles'],'independent_scalar_modes':1,'not_seven_term_chi3':True,'physical_accuracy_tolerance':None})
 dt=donor['reference_targets'];ds=dt['susceptibility'];assert ds['qualified']
 chi.append({'id':'donor/La_acetate','model_id':DONOR,'direction_lab':ds['direction_lab'],'target_source_au':ds['value']['estimates']['half']['chi_energy'],'reference_print_bound':ds['value']['estimates']['half']['roundoff_bounds']['chi_energy'],'selected_estimator':'energy curvature at saved half field','all_saved_estimators':ds['value']['estimates'],'source':ds['source'],'independent_scalar_modes':1,'physical_accuracy_tolerance':None,'limitation':ds['limitation']})
 moments.append({'id':'donor/La_acetate','model_id':DONOR,'target_mu_e_bohr':dt['moment']['total_dipole_lab_au'],'origin_A':dt['moment']['common_origin_A'],'method':'Original ECP PBE0-D4; D4 field-independent in qualified response','parameter_sensitive':True,'source':inputs['donor_center'],'source_record_pointer':'/reference_targets/moment','raw_record':dt['moment'],'moment_channel':'conjugate; unequal native Ed/Ep retained','consumed_development':True})
 moments.append({'id':'donor/Dy_acetate','model_id':dyid,'target_mu_e_bohr':dy['reference_targets']['moment']['total_dipole_lab_au'],'origin_A':dy['reference_targets']['moment']['common_origin_A'],'method':'Original ECP PBE0-D4; report-only Dy transfer consistency','parameter_sensitive':True,'source':inputs['donor_Dy_report_only'],'source_record_pointer':'/reference_targets/moment','raw_record':dy['reference_targets']['moment'],'moment_channel':'conjugate; unequal native Ed/Ep retained','consumed_development':True,'fit_primary':False})
 endpoint_exports={} if donor_endpoint_exports is None else donor_endpoint_exports
 assert set(endpoint_exports)<=set(ENDPOINTS)
 endpoints={};missing=[]
 for id in ENDPOINTS:
  r=displaced['rows'][id];assert r['baseline_compatibility']['compatible'];reference=r['reference'];row={'id':id,'geometry':r['geometry'],'coordinate':r['coordinate'],'original_inputs':r['inputs'],'baseline_compatibility':r['baseline_compatibility'],'saved_baseline_energy_kcal_mol':r['model']['energy_kcal_mol'],'saved_components':r['model']['components'],'reference':reference,'target_relative_energy_kcal_mol':(reference['energy_hartree']-dt['energy_hartree'])*units['hartree_kcal_mol'],'reference_center_energy_hartree':dt['energy_hartree'],'relative_energy_method':'same original ECP PBE0-D4 complex; includes D4 at both geometries','source':inputs['donor_displacements'],'native_operator_fixture':None,'consumed_development':True}
  if id in endpoint_exports:
   e=read(endpoint_exports[id]);assert e['geometry']==r['geometry'] and e['baseline']==inputs['shared_baseline'] and e['compatible'] and e['native_operator'] is not None
   fixtures[id]=e;fixture_pins[id]=endpoint_exports[id];row['native_operator_fixture']=endpoint_exports[id]
  else:missing.append(id)
  endpoints[id]=row
 donor_targets={'center_model_id':DONOR,'reference':dt,'center_reference_method':'Original ECP PBE0-D4','mechanical_modes':dt['mechanical_modes'],'canonical_rotation_corroboration_not_mixed_loss':dt['canonical_rotation'],'endpoints':endpoints,'absolute_complex_energy_is_not_a_pair_interaction_target':True,'native_physical_mu_channel':'permanent+B^T u','energy_conjugate_mu_channel':'permanent+B^T(u+v)/2','saved_conjugate_mu_lab_au':donor['baseline_observables']['conjugate_mu_lab_au']}
 # Stable public interface for the solver; all raw source records remain attached.
 for row in chi:
  row['fixture_id']=row['model_id'];row['direction']=row['direction_lab'];row['target_a0_cubed']=row['target_source_au']
 donor_targets.update(center_id=DONOR,moment_target_au=dt['moment']['total_dipole_lab_au'],field_direction=ds['direction_lab'],endpoint_rows=endpoints)
 trimer_dict={r['id']:dict(r,primary_fit='aquo9_' in r['id']) for r in trimer}
 return {'schema':'joint_demonstration_saved_data_v1','provenance':inputs,'shared_baseline':inputs['shared_baseline'],'units':units,'operators':fixtures,'fixture_pins':fixture_pins,'static':static,'references':refs,'reference_model_joins':refjoins,'pair_cycles':pair,'trimer_cycles':trimer_dict,'individual_moment_targets':moments,'chi_targets':chi,'donor':donor_targets,'donor_Dy_report_only':{'center_id':dyid,'reference':dy['reference_targets'],'baseline_observables':dy['baseline_observables'],'source':inputs['donor_Dy_report_only'],'fit_primary':False,'failed_Dy_susceptibility_not_a_target':True,'canonical_Dy_torque_not_available':True},'water_fields':fields,'water_force_records':forces,'water_parent_components_and_cycles':cycles['parent'],'water_common_basis_components':cycles['common_basis_decomposition'],'missing_native_endpoint_operators':missing,'complete_for_joint_demonstration':not missing,'target_weights':None,'coefficient_bounds':None,'coefficient_values':None,'consumption':'All records consumed development; no holdout claim.','method_channel_rules':['Water energy targets use electronic common-parent-basis PBE0 only; no D4 insertion.','Donor relative energy/slopes use original ECP PBE0-D4; canonicalCC rotation remains separate corroboration.','Use kernel.response_to_source_au for literal native unit conversion; report geometric chi separately.','Signed moments are origin0 and conjugate; report physical native-u moment separately when Ed!=Ep.','Pair/trimer coefficients act before the coupled solve; evaluate full term dictionaries, not isolated corrections.','One scalar per five susceptibility modes; alternative steps/estimators are consistency evidence.','Tighter W9far references remain separate sensitivity; no replacement of original targets.','Static metal-free terms have identically zero metal-centered correction; retain their saved observables.','Primary fit:34pair cycles,4W9 trimer cycles,5chi modes,La donor moment/two slopes. All8trimers andDy drycenter reported.']}
