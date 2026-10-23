# Directional quantum susceptibility of the consumed dry donor centers

## Question and exact ten-call scope

The original selected donor model underestimates the central dipole projected on the fixed rotation tangent: La0.114375 versus PBE0 0.580438 ea0 and Dy0.100684 versus0.635663 ea0. Its directional model susceptibilities are about54.12575/49.10297 a0³. The completed model40 result71d83185 reports all31 numerical checks passed; its independent actual review is separate. Quantum field curvature will test whether the moment mismatch is accompanied by a susceptibility mismatch in this same projection.

Use exactly the **two existing central dry mono-acetate geometries**. Each receives a compatible SP zero and four fields **−.001, −.0005, +.0005, +.001 au**: **ten new PBE0 single points**. No rotations, new geometry, fit, parameter point, new model call, central gradient, donor holdout or additional field size is included. Reuse the already computed model center fields and model zero energies.

Eight nonzero SPs would suffice for the previously defined dipole-slope estimator χ_mu at two steps. They would not supply the compatible E(0) required for both established energy-curvature checks. A four-energy replacement would change that estimator. This ten-cell design retains both existing estimators and odd-energy/zero-dipole controls; it does not assume the old EnGrad zeros interchangeable with SP zeros. Old center energies and common-origin dipoles remain contextual comparisons only, without an identity/accuracy gate against the new SPs.

## Frozen coordinates, field and native recipe

Use the exact central XYZ files, atom order, source input and basis assets bound by SOURCE_INVENTORY.json. Both metals retain the same coordinates and the metal at laboratory zero. The fixed unit vector is

```
t = (0.7562812564509392, -0.15594410634600103, -0.6353897204369424)
```

It is the original geometry-defined rotation tangent, shared with the completed model40 design. The field is F=s*t, with h=.001 and half=.0005 au. It is not chosen from response values. Preserve the full rendered vector precision and use the nominal s consistently in all estimators. Direction/rounding readback must agree with the frozen vector to floating-point/input-print precision.

Start with each actual central PBE0-D4 input. Remove only `EnGrad`; retain PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF, PModel, MaxIter300, the exact metal GTO/ECP/Aux-J blocks and all remaining scientific settings. Add `%scf EField Fx,Fy,Fz` and `EFieldOrigin 0,0,0`, and `%elprop Dipole true` / `Origin 0,0,0`. Use these same directives with an explicit zero vector in the zero SP. No automatic restart from an earlier field, new grid, solver or convergence threshold is permitted. Resource fields alone follow the actual grant.

All eight centers are real. Charge+2; effective restricted singlet; 40 explicit electrons;189 orbital functions, with counts[47,31,31,31,31,6,6,6]. Actual native Aux-J counts are La580 and Dy598. Real LaECP46/DyECP55 leave metal effective nuclear charge11. Physical parent La-singlet/Dy-sextet labels stay explicit and distinct from this frozen-f effective singlet. No ghost, local/canonical correlation or changed D4 convention is introduced.

ORCA field units are atomic units, with energy convention E(F)=E(0)−μ·F−½FᵀχF+…. Explicit property origins follow geometry units; zero is unambiguous. The native printed electric origin is in Bohr. Require actual field and both origin echoes, native property origin zero, full property coordinates, charge/electron/AO/Aux-J counts and finite dipoles. Use the qualified ORCA6.1.1 executable and inherited native parser/admission, not input keywords alone as evidence. Source/manual/native units evidence and the independently passed prior16-field quantum experiment are pinned in SOURCE_INVENTORY.json.

## Charge, energy stages and D4

Use full-precision property-file total dipoles at the laboratory origin. Nuclear and electronic vectors must close to the total and to the printed raw output within their token precision. The effective nuclear charges are [11,6,6,8,8,1,1,1], sum42;40 electrons give net+2. At origin0 require μ_nuclear=sum Z_eff*r. Do not add ECP core electrons back. For contextual old COM-origin values use μ_lab=μ_COM+2O, with O in Bohr, as already independently reconstructed in the dry-donor inventory. A common fixed origin shift changes all total dipoles by the same charge term, while susceptibility is origin-independent.

Retain native SCF electronic energy, final integration/COSX stages, D4 and final total separately. The energy estimator uses the full-precision **SCF_Energy.totalEnergy electronic energy** from the property file, joined to the native SCF total-energy token. Require property FinalEnergy to match native final total and SCF+D4 to close to final within the actual printed precision and the inherited native closure gate. The old D4-off field helper's equality FinalEnergy=SCF must be adapted explicitly, not copied into this D4-on panel.

At fixed geometry and charge the retained D4 term is field-independent. Require the raw D4 values within each metal's five calculations to match at their printed precision; report total-energy curvatures and their closure as a check. Do not double-subtract D4 or add a field-dependent dispersion correction. Preserve every energy stage and a failed stage join; no subtraction of observed COSX changes to make a gate pass.

## Estimators and sixteen numerical checks

Let E be the matched electronic SP energy in Hartree and μ=t·μ_total(0), in ea0. For each metal and s∈{h,half}:

```
chi_E(s)  = -[E(+s)+E(-s)-2E(0)]/s²
chi_mu(s) = [mu(+s)-mu(-s)]/(2s)
mu_E(s)   = -[E(+s)-E(-s)]/(2s)
```

Use the prior declared numerical resolutions unchanged:

- Per metal, |chi_E(s)−chi_mu(s)|≤**0.5 a0³** at both steps; |chi_E(h)−chi_E(half)|≤0.5 and |chi_mu(h)−chi_mu(half)|≤0.5 a0³; |mu_E(s)−mu(0)|≤**1e−4 ea0** at both steps. Six checks per metal.
- Apply the four curvature comparisons directly to Dy−La with **1 a0³**, retaining both signed metal values and all coefficients. Four paired checks. Report paired odd slopes and bounds descriptively.

All **ten native slots and sixteen checks** remain explicit. Missing/nonfinite/unadmitted data never pass. Retain both raw steps and both estimators; no averaging, extrapolation, new field step, retuned threshold or omitted failed check. Report full dipole vectors, projected estimates, even-moment response and old-EnGrad/new-SP differences.

Propagate actual decimal half-last-unit and binary conversion bounds. Projection uses sum |t_i|*b_mu_i; linear combinations use absolute coefficients. For example b_chiE=(b_plus+b_minus+2b_zero)/s² and b_chimu=(b_muplus+b_muminus)/(2|s|). These are printing/arithmetic bounds only. SCF tolerance is not an energy-error bound. The .5/1 a0³ and1e−4 ea0 resolutions are consistency checks, not physical accuracy margins or confidence intervals. Failure can reflect field-step nonlinearity, numerical convergence or energy/property-stage inconsistency; its cause remains unresolved and continuation stops.

## Interpretation and staging

Compare qualified quantum susceptibilities with the already saved original-model center chi_E at both steps, plus their Dy−La differences. No new cross-method accuracy cutoff is added. A susceptibility discrepancy establishes an additional electrical-response defect beyond the known moment offset in this direction. Similar susceptibility with a moment offset is consistent with a driving/permanent/equilibrium-response mismatch, but does not uniquely establish any one cause. This test cannot uniquely identify permanent versus induced contributions, charge transfer, a local parameter defect, or the source of the zero-field rotational energy error. A uniform field probes a global mode, unlike the spatially varying metal field.

PBE0 method, basis/diffuse, frozen-f/physical-state and single-reference limitations remain. Native Energy Check with other failed printed SCF criteria, and field-dependent final-integration stages, must remain visible. Numerical consistency does not certify physical reference accuracy. The separate canonical rotation test addresses a different method-sensitivity question and is not replaced by this panel.

First run the complete **La five-cell block**, including its compatible zero; require actual native/property admission, all six numerical checks and cost review before the exact Dy five-cell block. Root owns execution/releases. Four MPI ranks, one node, default scheduler memory/time, MaxCore derived from75% of actual grant per rank. No automatic retry or duplicate zero. Existing eight nearby rotation SPs took141.97–166.24s each at four ranks; the mean estimates about**1.70 allocated CPU-hours** for ten calls. This is an estimate, not a guarantee; actual field cost and behavior govern continuation. Design/source preparation makes zero molecular calls.
