# Two fixed native witnesses: scientific interpretation

The nonlinear model reproduces the local tradeoff: **B-half improves component errors while preserving all25 physical gates and all three current electrical norms; D-half gains much more component accuracy but fails13 physical gates.** This is evidence about these two consumed directions, not a capacity theorem, model selection or production admission. Root and the independent actual auditor own admission separately.

Pinned [producer result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/two_direction_native_validation_v1/run_v1/execution/RESULT.json) (`1604c8ff`), [root completion](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/two_direction_native_validation_v1/run_v1/ROOT_COMPLETION.json) and [fixed plan](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/two_direction_native_validation_v1/PLAN.json); complete direct pins are in [INTERPRETATION_INPUTS.json](INTERPRETATION_INPUTS.json).

## Observed versus predicted

P and B are paired Dy−La pair and collective energy errors across W8/W9. q is the fixed-origin minimax squared-error ratio.

| Point | q | P RMS, kcal/mol | B RMS, kcal/mol | Original gates |
|---|---:|---:|---:|---:|
| Actual final center | .227271822 |12.720316|10.441828|25/25|
| B-half predicted | .208871747 |12.040271|10.193096|25/25|
| B-half native | .208830475 |12.040218|10.192088|25/25|
| D-half predicted | .056818176 |6.360170|5.221169|12/25|
| D-half native | .056837931 |6.361276|5.193072|12/25|

Maximum absolute differences in the four component-error predictions are.001054 kcal/mol for B-half and.029267 for D-half. Both reproduce every predicted gate pass/fail flag. B-half still lies close to boundaries: W9 exchange error−.996841 kcal/mol and paired water-plus force RMS1.996697 kcal/mol/Å against limits1 and2. D-half fails four energy gates and all nine La/Dy/paired single-water force groups; paired water-minus force RMS is7.066414 kcal/mol/Å. Its improved shell energies are not a usable-model result.

B-half's old E/F objective also improves from.386590 to.365160. Its saved `fixed_ceiling_pass` flags refer to historical origin24 ceilings, not the current-center electrical ceilings used in case B; those historical flags remain descriptive.

## Electrical response

These source-defined residual RMS values are descriptive, with no new physical acceptance tolerance:

| Point | Pair dipole, D | Collective dipole, D | Susceptibility, bohr³ |
|---|---:|---:|---:|
| Actual center |.639432|.085846|.680444|
| B-half predicted |.639430|.085812|.679596|
| B-half native |.639428|.085802|.679601|
| D-half predicted |.640616|.089201|.652126|
| D-half native |.640619|.089200|.652029|

The W9 collective dipole Dy−La z-contrast remains positive: B-half native+.128832819 D (predicted+.128840145), D-half native+.148448368 D (predicted+.148449000), versus reference−.104617953 D. B-half susceptibility errors remain La−.883783/−.902995 and Dy+.367083/+.340890 bohr³ for dimer/trimer. Conditional increment errors are smaller (La−.019213,Dy−.026193; difference-of-increments−.006981 bohr³), so a close increment does not cure the absolute or signed-response disagreement.

## Mechanics decision

The large gate-relaxed gain and corresponding dimer-force failures survive nonlinear evaluation; solver termination and linearization error do not explain this tested tradeoff. The earlier four-case comparison found little gain from removing box faces and only a modest penalty from preserving electrical norms. The next design decision should therefore examine a specific missing or untested freedom in the existing radial/coordination response, while retaining physical gates. The current functional already has distinct pair and collective drive/stiffness coefficients and vdW R/epsilon; this evidence does not justify merely adding a generic many-body term or choosing a unique microscopic cause. A separate bounded source-based design recommendation follows; no new fit or term is executed here.

B-half demonstrates remaining feasible improvement in the present model, but its12.04/10.19 kcal/mol component errors, sign disagreement and narrow gate margins remain material. Both points are inside the original parameter box. The targets retain their consumed common-basis/PBE0 reference limitations; no new reference, held-out outcome or free-energy evidence was obtained.

## Accounting

Exactly2 vectors,92 native cell records and112 program invocations completed; no solver or reference call. Producer cost:23.841330 wall seconds and74.781135 CPU-seconds (25.027738 controller+22.087874 child user+27.665523 child system), truthful login execution. No inherited or concurrent per-program child counters were added. This report only reads saved metrics.

## All25 frozen physical gates

The native pass/fail pattern matches the frozen prediction for both points. Energy values are signed kcal/mol; force values are RMS / maximum body component in kcal/mol/Å. Limits are |energy|≤1, force RMS≤2 and maximum≤5.

| Gate | B-half native | D-half native | B / D pass |
|---|---:|---:|---|
|energy/balanced_exchange/water_minus|-0.240471|-3.125077|Pass / Fail|
|energy/balanced_exchange/water_central|-0.396076|-1.318842|Pass / Fail|
|energy/balanced_exchange/water_plus|+0.062697|+0.173275|Pass / Pass|
|energy/balanced_exchange/aquo9_central|-0.996841|-0.321024|Pass / Pass|
|energy/individual_radial_work/La_minus|+0.029490|+0.807465|Pass / Pass|
|energy/individual_radial_work/La_plus|-0.259691|-0.715602|Pass / Pass|
|energy/individual_radial_work/Dy_minus|+0.185095|-0.998769|Pass / Pass|
|energy/individual_radial_work/Dy_plus|+0.199083|+0.776514|Pass / Pass|
|energy/paired_radial_work/minus|+0.155605|-1.806235|Pass / Fail|
|energy/paired_radial_work/plus|+0.458774|+1.492117|Pass / Fail|
|force/La/water_minus|0.327709 / 0.632855|3.597637 / 7.141810|Pass / Fail|
|force/La/water_central|0.516571 / 1.025902|2.497206 / 4.954231|Pass / Fail|
|force/La/water_plus|1.101744 / 2.184442|2.222824 / 4.405564|Pass / Fail|
|force/La/aquo8_central|1.606120 / 4.951625|1.387532 / 4.270861|Pass / Pass|
|force/La/aquo9_central|1.630916 / 3.632798|1.694157 / 3.632401|Pass / Pass|
|force/Dy/water_minus|1.973169 / 3.917508|3.470890 / 6.863942|Pass / Fail|
|force/Dy/water_central|0.283110 / 0.554857|3.006399 / 5.951974|Pass / Fail|
|force/Dy/water_plus|0.894979 / 1.773582|2.178974 / 4.317311|Pass / Fail|
|force/Dy/aquo8_central|1.124831 / 3.371436|1.205927 / 3.631194|Pass / Pass|
|force/Dy/aquo9_central|1.995798 / 3.884181|1.944333 / 3.884719|Pass / Pass|
|force/paired/water_minus|1.661798 / 3.284653|7.066414 / 14.005752|Pass / Fail|
|force/paired/water_central|0.796204 / 1.580759|5.502977 / 10.906205|Pass / Fail|
|force/paired/water_plus|1.996697 / 3.958024|4.401781 / 8.722874|Pass / Fail|
|force/paired/aquo8_central|0.500349 / 1.580225|0.194008 / 0.639702|Pass / Pass|
|force/paired/aquo9_central|1.158987 / 2.915792|0.792550 / 2.403742|Pass / Pass|
