"""One saved-artifact inventory; no fitting, raw re-admission or native calls."""
import hashlib, json, math, time
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def read(path):return json.loads(Path(path).read_text())

def rank(rows):
    a=[[Fraction(x) for x in row] for row in rows];r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];v=a[r][col];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][col];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def main():
    wall,cpu=time.monotonic(),time.process_time()
    out=HERE/'TARGET_INVENTORY.json';assert not out.exists()
    paths={
        'development_data':W/'explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1/preparation_v2/DATA.json',
        'development_objective':D/'explicit_model_v1/paired_water_coordination_joint_v1/OBJECTIVE.py',
        'body_helper':D/'explicit_model_v1/matched_dry_v1/OBJECTIVE.py',
        'development_declaration':D/'explicit_model_v1/paired_water_coordination_joint_v1/ROOT_DECLARATION.md',
        'dipole_comparison':W/'redistribution_comparison_v1/actual_v1/RESULT.json',
        'dipole_arithmetic_review':W/'redistribution_comparison_v1/independent_arithmetic_review_v1/REVIEW.json',
        'quantum_full_precision_review':W/'reference_chemistry_v1/saved_electronic_redistribution_v1/independent_actual_review_v1/actual_v3/REVIEW.json',
        'model_moments':W/'explicit_model_v1/saved_redistribution_counterpart_v1/run_v2/execution/RESULT.json',
        'trimer_comparison':W/'explicit_model_v1/collective_mechanics_options_v1/collection_v1/actual_v1/RESULT.json',
        'trimer_native_admission':W/'explicit_model_v1/collective_mechanics_options_v1/actual_review_v1/actual_v1/RESULT.json',
        'field_comparison':W/'explicit_model_v1/conditional_susceptibility_v1/comparison/run_v1/execution/RESULT.json',
        'field_specification':D/'reference_chemistry_v1/conditional_susceptibility_v1/SPEC.md',
        'field_zero_inventory':D/'reference_chemistry_v1/conditional_susceptibility_v1/SOURCE_INVENTORY.json',
        'reserved_metadata':W/'explicit_model_v1/matched_dry_v1/reserved_validation_v1/RESERVED_PANEL.json',
        'fresh_radial_declaration':D/'explicit_model_v1/water_coordination_radial_validation_v1/ROOT_DECLARATION.md',
        'fresh_radial_preparation_review':W/'explicit_model_v1/water_coordination_radial_validation_v1/ROOT_PREPARATION_REVIEW.json',
        'donor_validation_declaration':D/'explicit_model_v1/paired_donor_validation_stageA_v1/ROOT_DECLARATION.md',
        'donor_terminal_report':D/'explicit_model_v1/paired_donor_fit_v1/FINAL_REPORT.md',
    }
    artifacts={k:pin(v) for k,v in paths.items()}
    data=read(paths['development_data']);dip=read(paths['dipole_comparison']);qr=read(paths['quantum_full_precision_review'])
    model=read(paths['model_moments']);trimer=read(paths['trimer_comparison']);field=read(paths['field_comparison'])
    assert dip['inputs']['quantum_review']==artifacts['quantum_full_precision_review'] and dip['inputs']['model_result']==artifacts['model_moments']
    assert len(trimer['cycles'])==8 and len(field['table'])==28 and field['numerical_checks_passed']
    assert model['candidate_result']==data['provenance']['candidate_result']
    artifacts['development_selected_review']=pin(data['provenance']['selected_review']['path'])
    assert artifacts['development_selected_review']==data['provenance']['selected_review']
    artifacts['trimer_geometry']=pin(trimer['inputs']['manifest']['path']);assert artifacts['trimer_geometry']==trimer['inputs']['manifest']
    artifacts['field_reference_collection']=pin(field['configuration']['reference_collection']['path'])
    assert artifacts['field_reference_collection']==field['configuration']['reference_collection']
    geoms=['water_minus','water_central','water_plus','aquo8_central','aquo9_central']
    ids=[q+'_'+m for q in geoms for m in ('La','Dy')];assert set(ids)==set(data['rows'])
    energies=[]
    for q in geoms[:3]+[geoms[-1]]:
        energies.append(dict(id='exchange__'+q,coefficients={q+'_Dy':1,q+'_La':-1,'aquo8_central_Dy':-1,'aquo8_central_La':1}))
    for m in ('La','Dy'):
        for sign in ('minus','plus'):
            energies.append(dict(id='work__'+m+'_'+sign,coefficients={'water_'+sign+'_'+m:1,'water_central_'+m:-1}))
    for sign in ('minus','plus'):
        energies.append(dict(id='paired_work__'+sign,coefficients={'water_'+sign+'_Dy':1,'water_central_Dy':-1,'water_'+sign+'_La':-1,'water_central_La':1}))
    for e in energies:
        e.update(reference_kcal_mol=math.fsum(v*data['rows'][k]['reference']['energy'] for k,v in e['coefficients'].items()),
            selected17676_model_kcal_mol=math.fsum(v*data['baseline_model'][k]['energy'] for k,v in e['coefficients'].items()),absolute_error_limit_kcal_mol=1)
    erank=rank([[e['coefficients'].get(k,0) for k in ids] for e in energies]);assert erank==6
    force=[]
    for kind in ('La','Dy','paired'):
        for geom in geoms:
            q=data['baseline_score']['per_fixture'][kind][geom]
            force.append(dict(id=geom+'__'+kind,reference_pointer='/baseline_score/per_fixture/'+kind+'/'+geom+'/reference_body_response',
                endpoint_coefficients={geom+'_Dy':1,geom+'_La':-1} if kind=='paired' else {geom+'_'+kind:1},
                RMS_limit_kcal_mol_A=2,maximum_component_limit_kcal_mol_A=5,selected17676_passed=q['passed']))
    dev=dict(source='development_data',consumption='consumed_development_fit',cells={k:dict(reference_pointer='/rows/'+k+'/reference',
        reference_native_result=data['rows'][k]['reference_result'],geometry_pointer='/rows/'+k+'/geometry',topology_bodies_pointer='/rows/'+k,
        model_input_pointer='/baseline_model/'+k,atoms=data['rows'][k]['geometry']['atoms']) for k in ids},
        energy_checks=energies,force_checks=force,checks=25,energy_checks_count=10,force_checks_count=15,energy_linear_rank=erank,
        units=dict(energy='kcal/mol',positions='angstrom',forces='kcal/mol/angstrom',torque_response='centroid torque divided by1angstrom'),
        conventions=dict(reference_force='negative energy gradient',errors='model minus reference',paired='Dy minus La',
            reference_recipe='Uncorrected original ECP-PBE0-D4 energy/gradient; do not substitute CP or CC fragment energies',
            body_weight='Per category equal; within category equal bodies; within body equal response types and xyz. Ion net force only; water net force and centroid torque/1A.'),
        objective=dict(source='development_objective',residual_count=413,energy_residuals=8,force_residuals=405,
            block_weights=dict(energy=.25,La=.25,Dy=.25,paired=.25),energy_subblock_weights=dict(exchange=.5,individual_radial_work=.5),
            paired_radial_works_in_loss=False,force_components_per_metal=135),
        redundancy=['Paired force responses are exactly Dy minus La; five paired fixtures add no new independent vector targets.',
            'Paired radial work = Dy work minus La work = exchange(displaced)-exchange(central).',
            'The ten scalar energy gates span six endpoint-energy combinations. All ten original acceptance gates remain required.',
            'Atomic force residuals and body projections derive from the same ten atomic force sets; no independent-likelihood interpretation.'])
    conversion=dip['Debye_per_atomic_unit'];pair=[]
    for n in (8,9):
        for metal in ('La','Dy'):
            gid=f'aquo{n}_central_{metal}';ion=f'pbe0__aquo{n}__{metal}__ion_ghost_shell'
            for i in range(1,n+1):
                stem=f'pbe0__aquo{n}__{metal}__water_{i:02d}';terms={stem+'__metal_water':1,ion:-1,stem+'__water_only':-1}
                modelterms={gid+f'__dimer_{i:02d}':1,'ion_'+metal:-1,gid+f'__water_{i:02d}':-1}
                mu=[math.fsum(a*qr['rows'][k]['total_dipole_au'][j] for k,a in terms.items())*conversion for j in range(3)]
                mm=[math.fsum(a*model['rows'][k]['moments']['induced_sum_Debye'][j] for k,a in modelterms.items()) for j in range(3)]
                bound=math.fsum(abs(a)*qr['rows'][k]['dipole_roundoff_bounds_au']['total_dipole_au'] for k,a in terms.items())*conversion
                pair.append(dict(id=gid+f'__pair_{i:02d}',parent=gid,water_index=i,reference_terms=terms,model_terms=modelterms,
                    reference_pair_mu_Debye=mu,selected17676_pair_mu_Debye=mm,retained_reference_print_bound_Debye_per_component=bound))
            summed=[math.fsum(x['reference_pair_mu_Debye'][j] for x in pair if x['parent']==gid) for j in range(3)]
            assert max(abs(summed[j]-dip['parents'][gid]['P']['quantum_full_precision_check_Debye'][j]) for j in range(3))<1e-12
    assert len(pair)==34
    dipoles=dict(consumption='consumed_mechanistic_development',reference_source='quantum_full_precision_review',model_source='model_moments',
        comparison_source='dipole_comparison',units='Debye',Debye_per_atomic_unit=conversion,frame='original parent laboratory frame; origin0',
        reference_convention='Total nuclear+electronic dipoles after common-origin conversion; ECP effective charge11, ghost charge0. Full parent orbital+aux basis.',
        model_convention='Charge-balanced fixed-coordinate cycles cancel permanent moments; induced moment sum equals total response for this unchanged model. Do not compare absolute quantum total dipole directly to model induced-only dipole.',
        pair_cycles=pair,parents={g:{k:dict(reference_mu_Debye=v['quantum_full_precision_check_Debye'],selected17676_model_mu_Debye=v['model_Debye'],
            pointer='/parents/'+g+'/'+k,reference_print_bound_Debye_per_component=v['quantum_print_bound_Debye_per_component'],
            model_print_bound_Debye_per_component=v['model_print_bound_Debye_per_component']) for k,v in pp.items()} for g,pp in dip['parents'].items()},
        formulas=dict(pair='mu(MWi)-mu(M)-mu(Wi)',C='mu(MWn)-mu(M)-mu(Wn)',P='sum_i pair_i; ion subtracted n times',B='C-P'),
        redundancy=['P is the sum of the34 literal parent-specific pair cycles.','Each parent C=P+B; do not treat all three as independent observations.',
            'W8 net vectors nearly cancel by symmetry; individual pair vectors remain nonzero. Do not fit near-zero net-vector relative errors.',
            'Norms, xyz components, paired-metal differences and cross-coordination differences reuse the same vectors.'],
        physical_acceptance_threshold=None,prospective_loss=None,print_bounds_are_accuracy_bounds=False)
    tri=dict(source='trimer_comparison',geometry_source='trimer_geometry',consumption='consumed_mechanistic_development',count=8,
        units=dict(E3='kcal/mol',mu3='Debye',positions='angstrom'),reference_recipe='Electronic PBE0 in full parent basis, no D4',
        formula='M12-M1-M2-12+M+1+2 for both energy and total common-origin dipole',
        cycles={k:{a:v for a,v in row.items() if a in ('selection','reference_terms','model_terms','reference_electronic_E3_kcal_mol',
            'model_E3_kcal_mol','reference_mu3_Debye','model_mu3_Debye','reference_print_bound_Debye_per_component','model_print_bound_Debye_per_component')} for k,row in trimer['cycles'].items()},
        redundancy=['All cycles share admitted ion and individual-water/dimer terms with pair/C/P/B observables.',
            'Do not sum the eight geometry-selected E3 values to the full-shell B; unmeasured triplets and higher orders remain.',
            'mu3=-dE3/dF at zero field; it is a new derivative observable relative to E3, but an odd finite-field slope on the same cycle remeasures it.'],
        physical_acceptance_threshold=None,prospective_loss=None)
    responses=dict(source='field_comparison',specification='field_specification',consumption='consumed_mechanistic_development',
        state_ids=field['matched_state_ids'],state_count=4,derived_increment_count=2,derived_paired_count=1,reported_rows=field['table'],
        units=dict(field='atomic units',susceptibility='a0^3',dipole='e*a0',energy='hartree'),field_axis='laboratory z',field_steps_au=[.001,.0005],
        shared_geometry='W9 far waters1,3 trimer and matching M+water1 dimer, both La/Dy; full same-parent basis in reference',
        definitions=dict(chi_energy='-(E(+s)+E(-s)-2E0)/s^2',chi_moment='(mu_z(+s)-mu_z(-s))/(2s)',
            odd_energy_moment='-(E(+s)-E(-s))/(2s)',increment='chi(trimer)-chi(dimer)',paired='increment(Dy)-increment(La)'),
        numerical_consistency_gates=dict(state_a0_cubed=.5,increment_and_paired_a0_cubed=1,odd_slope_e_a0=1e-4,passed=field['numerical_checks_passed']),
        independent_physical_quantity_count=4,reference_response_pointer='/reference_response',model_response_pointer='/model_response',
        redundancy=['Two steps and two estimators cross-check each of four state susceptibilities; they are not16 independent state targets.',
            'The two increments and their paired contrast derive from those four states; all28 report rows are correlated.',
            'The four zero-field dipoles reuse pair/trimer component moments; odd slopes are consistency checks, not independent new targets.',
            'This is response to adding a second water, not alpha3: added-water intrinsic, pair and metal-response contributions remain.',
            'Uniform-z response samples one global mode; it does not determine the response to an inhomogeneous metal field.'],
        physical_acceptance_threshold=None,prospective_loss=None)
    reserved=read(paths['reserved_metadata']);g=reserved['geometries'];gids=list(g) if isinstance(g,dict) else [x['geometry_id'] for x in g]
    water=[k for k in gids if 'aquo8' in k];donor=[k for k in gids if 'aquo8' not in k and any(x in k for x in ('acetate','acetamide'))];protein=[k for k in gids if k not in water+donor]
    assert (len(water),len(donor),len(protein))==(4,12,8)
    exclusion=dict(metadata_source='reserved_metadata',reserved_numerical_outcomes_read_for_this_inventory=False,
        historical_water_orientation=dict(ids=water,status='previously_inspected_validation_regression_not_fresh',included_in_repair_targets=False,
            authority='fresh_radial_declaration explicitly classifies earlier four orientations as historical regression evidence'),
        donor_orientation=dict(ids=donor,status='reserved_unopened_held_after_failed_donor_development',included_in_repair_targets=False,
            authority=['donor_terminal_report','donor_validation_declaration']),
        protein_stage_B=dict(ids=protein,status='reserved_unexecuted_separate_stage_B',included_in_repair_targets=False,authority='donor_validation_declaration'),
        fresh_W9_radial=dict(ids=[f'aquo9_{s}_{m}' for s in ('minus','plus') for m in ('La','Dy')],status='reserved_held_unbound_after_failed_joint_development',
            included_in_repair_targets=False,original_controls=12,authority=['fresh_radial_declaration','fresh_radial_preparation_review']))
    inventory=dict(status='saved_consumed_target_inventory',implementation=pin(__file__),artifacts=artifacts,
        fixed_model_candidate=data['provenance']['candidate_result'],candidate_pin_authority='development_data provenance plus selected raw review; no fit trajectory reread',
        fixed_six_parameters=data['baseline_parameters'],groups=dict(original_energy_force=dev,pair_and_collective_dipoles=dipoles,trimer_energy_and_dipoles=tri,conditional_susceptibility=responses),
        exclusions=exclusion,policy=dict(fitting_performed=False,new_loss_defined=False,new_weights_defined=False,existing_thresholds_changed=False,
            new_observables_descriptive_until_prospective_loss=True,native_calls=0,physical_state_or_reference_accuracy_qualified=False,
            quantum_limitations='Provisional same frozen-f ECP/PBE0 representation, finite basis/SCF limitations; CP and correlated sensitivities do not silently replace the original energy/force target.'),
        cost=dict(wall_seconds=time.monotonic()-wall,CPU_seconds=time.process_time()-cpu,local_threads=1))
    out.write_text(json.dumps(inventory,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(inventory=pin(out),bytes=out.stat().st_size,counts=dict(energy=10,force=15,pair_dipoles=34,parent_C_P_B=12,trimer_cycles=8,field_states=4,field_rows=28),cost=inventory['cost'])))

if __name__=='__main__':main()
