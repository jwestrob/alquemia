# Do fixed water well depths limit the radial force response?

## Current decision and boundary

**Deferred option: continue the same four parameters first.** Root selected a separately declared, finite 32-additional-step continuation after the actual radial result showed continuing improvement at its iteration cap. The result does not establish that the current parameter family failed or stalled. No epsilon calculation is selected or released.

If epsilon investigation is subsequently warranted, use **one matched six-column derivative panel at one frozen point**, reusing exact existing calculations. This tests whether the two water well depths add useful response directions beyond the four radius/damping parameters. Maximum cost is **12 parameter-probe vectors, 48 new physical cells and 96 native programs**. If a finite six-parameter refinement is separately authorized, use this as its initial derivative batch and reuse it, rather than repeat the same derivatives after a separate diagnostic.

The proposed probe design and interpretation below were written before inspecting the terminal radial outcome; the terminal addendum records the subsequent exposure and current decision. This file selects no starting vector and releases no calculation or optimizer. Root must name and pin any future diagnostic point, its eight native records, exact keys and unchanged objective. Preserve the completed four-parameter protocol and all results. No threshold, weight, reference, geometry, solver, donor parameter or reserved-validation change is proposed here.

## What epsilon changes in this Hamiltonian

The pinned native buffered14-7 interaction can be written as

`U(r) = epsilon * ((1+delta)/(rho+delta))^7 * ((1+gamma)/(rho^7+gamma) - 2)`,

where `rho=r/Rmin`, `delta=0.07`, and `gamma=0.12`. Native `ehal1.f` multiplies both the energy and its radial derivative by the pair epsilon. At fixed coordinates and Rmin, epsilon therefore scales this pair contribution linearly. It changes the balance of repulsion/dispersion against the unchanged electrostatic/polarization forces; it does **not** independently reshape the buffered14-7 function or change its exponents.

The current effective metal–water-O epsilons are La **0.24429636467571708** and Dy **0.37332330078584813 kcal/mol**. Perturb these explicit pair values only: `VDWPR 34 500 Rmin epsilon` for La and `VDWPR 34 501 Rmin epsilon` for Dy. Native `kvdw.f` installs special pair Rmin/epsilon directly after ordinary combining rules. Do not change homoatomic ion `VDW` values: that would also change other combined pairs. Keep full pair Rmin as printed; do not halve it.

`POLPAIR` is keyed by atom types, so retain `POLPAIR 36 500/501 damping` unchanged during epsilon probes. Water-O type36 is class34. Native `kpolar.f` and `damping.f` use the Thole value, polarizabilities and distances; they do not take this VDWPR epsilon as input. At fixed coordinates, changing only epsilon should leave induced-dipole/polarization behavior unchanged up to solver/output precision. Water-H, donor pairs, permanent multipoles, intrinsic polarizabilities and all covalent terms remain untouched.

Use the existing native engine, parser and saved authorities. These equations explain the diagnostic; they are **not** a proposal to implement a replacement analytic energy engine.

## Matched probes

Retain exactly `water_{minus,central,plus}_{La,Dy}` and `aquo8_central_{La,Dy}`. Every vector has a complete eight-cell result, including explicit reuse of the opposite metal. At the frozen point:

1. Reuse the eight complete central E/force records; no central calls.
2. For the original four variables, retain the existing parameter finite-difference recipe: plus/minus0.005 of each original allowed interval, clipped to its original bounds, divided by the actual clipped span. These produce eight vectors.
3. For each metal's epsilon, use symmetric **±5% of its fixed baseline epsilon**, with all other values held fixed. These produce four vectors. Define dimensionless `eta_M=(epsilon_M−epsilon_M,0)/epsilon_M,0`; the new columns are divided by the actual eta span0.10. Store derivatives in both physical epsilon units and this dimensionless convention so later coordinate rescaling needs no repeated native calls.

Each vector changes one metal only: four affected physical cells, each needing one analyze plus one analytic `testgrad Y N`. Thus12×4×2=96 programs before exact reuse. These are parameter finite differences; they do not run Cartesian numerical differentiation and do not introduce its internal displaced-energy count. Reuse requires exact effective input, engine, state, raw-output and native-admission identity, not approximate parameter agreement or a relabeled cache key.

The last radial iteration's Jacobian belongs to its **iteration origin**. It is reusable here only when that exact parameter vector equals the frozen diagnostic point and all probe inputs/step widths match. Do not present it as a final-point Jacobian merely because it came from the last iteration.

A four-vector epsilon-only panel costs at most32 native programs. It is sufficient for the same six-column comparison **only if** an exact matching four-column panel already exists. Otherwise it is a weaker directional test at the chosen point; comparison to a nearby historical Jacobian must be labeled approximate and cannot establish a matched six-parameter local span. The recommended96-program ceiling avoids that ambiguity.

At most eight CPU workers, one native process/thread each, root-owned Slurm submission. Stop new starts on the first native/nonfinite/polarization failure, retain running siblings and all eight-cell denominators, and do not retry or expand the probe set automatically.

## Analysis without a fit or a selected new candidate

Use the exact frozen radial weighted-residual vector `r`: the same balanced-energy residuals, individual La/Dy and paired rigid responses, four0.25 blocks and per-fixture reporting. No raw unlike-element offset, unbalanced free fragment, new weight or statistical-independence assumption is added.

Build `J4` from the four original variables and `Je` from the two epsilon variables, at the **same point**. Preserve signed derivatives and actual spans. Report:

- Each epsilon column's effect on balanced exchanges, individual and paired radial works, every per-fixture rigid response, and outward force slopes. For original water side fixtures the moving body is the metal and `dE/dq=−F_metal·u`; do not replace the actual path with rigid-water translation.
- The size of each epsilon response outside the old span: `W=(I−P4)Je`, plus `norm(W_j)/norm(Je_j)` when the denominator is nonzero. Report raw norms and singular values as well; a normalized angle alone can conceal a negligible response.
- Residual alignment with these additional directions. Construct the old retained basis using the existing relative SVD cutoff1e−6. Use the corresponding **same absolute singular threshold** for the orthogonal increment W, not a new threshold relative to a tiny W. Report the nested old/extended projections and residual components by all four blocks and individual fixtures. This is a local linear projection, not a fitted parameter vector or achieved improvement. Retain discarded singular values and the coordinate normalization explicitly.
- Actual scores/control flags of all symmetric probes for transparency, without selecting a winner or changing parameters. Keep per-fixture energy≤1 kcal/mol, rigid RMS≤2 kcal/mol/Å and maximum body-component≤5 kcal/mol/Å as the existing development resolutions.

For epsilon-only probes the midpoint energy and force should reproduce the frozen center because only a linear VDW pair coefficient changed. Compare `(plus+minus)/2` to center using the existing replay scales:0.0001 kcal/mol for energy and0.001 kcal/mol/Å per force component. Inspect native energy-component and polarization receipts as available. A failure indicates an input/admission/numerical inconsistency to resolve; it does not establish that epsilon physics is nonlinear. This parameter-linearity check does not replace the separate Cartesian force qualification.

Do not infer a final-point rank from an earlier iteration. Do not solve for or execute a newly proposed parameter vector in this diagnostic. A later six-parameter optimizer would need its own frozen positive-epsilon bounds, trust/scaling conventions, finite budget and root release; it can reuse these exact native derivatives with explicit coordinate conversion.

## Interpretation

| Observation | Supported conclusion |
|---|---|
| Epsilon adds resolved directions aligned with currently unexplained residuals | Fixed well depths are a plausible local restriction; a bounded joint test may be useful |
| Epsilon directions largely lie in the old span | Little additional local freedom at this point; not proof of global redundancy |
| Additional directions exist but barely align with remaining errors | Epsilon freedom is unlikely to address these residuals locally |
| Aggregate linear projection improves but individual force errors remain | The stated mechanics target remains unmet; do not hide failed fixtures in the loss |
| Actual probes improve some errors or pass development checks | Finite observed response only; no automatic new-candidate selection, transfer or physical qualification |

Even a useful epsilon direction does not identify the well depth as the unique physical defect. Missing short-range polarization behavior, the pair functional form, other parameters or reference uncertainty remain alternatives. Conversely, a budget-limited four-parameter fit is not proof that its family cannot work.

## Primary local evidence and pins

All paths below are relative to the project root. No new external literature was needed because the executed source and its primary manual are available locally.

Radial SPEC: `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_repair_v1/SPEC.md`, SHA `0e937e4ef9ddf3d768ff93509c5e8e705a5384956895fa3e46efc89e8bda7f41`; objective `5b83cf5f3e77aeb96ec0f4d0a0c3f765883f406e8f4110648972e5bac0996f6b`; evaluator `1cf6025e27291523ac7c7242a07245da22c8927d67f4323c87392974133eda56`.

Native paths share prefix `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source/`:

| File | Relevant definition | SHA256 |
|---|---|---|
| `source/ehal1.f` | Pair energy/derivative, especially lines237–260 | `b98da2d7ad01d7888966ddb2b4281ce0b4266306ff71147b7a3286af6f9eb197` |
| `source/kvdw.f` | Ordinary diameter conversion and later direct pair override | `c9646b0e1bdf79b29164a152a304f16ecbcfce4670eb4d1a6fed4e3af59024a1` |
| `source/kpolar.f` | Type-indexed Thole pair overrides | `01c5cba4c5f7b9a14544b8084f57227aa5d8317f283ac550da653b619f024d8b` |
| `source/damping.f` | Thole damping depends on distance/polarizability/pair damping | `8e7b99d28a30548b4c7ffaff78794848c6c0bdbc6066d09857379445c232e4ee` |
| `source/initprm.f` | Gamma0.12,delta0.07 defaults | `4c48bfb322169acad353bc1b390db16a3923eef684c652ce4a547f94f254496a` |
| `source/prmkey.f` | Explicit overrides of those defaults | `56207706e8294e4b592a29c47b7b0794ca7655d59ac05778aba4f05b144afef4` |
| `doc/sphinx/tinker/text/keywords.rst` | VDWPAIR/POLPAIR and buffered14-7 documentation | `a9eed423d4f4edf14dda4b888066e471b75227f652840da254933000cb121b65` |

## Terminal grounding and deferred execution

The saved result from job17589 was inspected only after root reported it terminal. Authority: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_repair_v1/run_v1/execution/RESULT.json`, SHA `31cb5caa7f96a3f312cc3367964149923665faa16b41d8993e664ad22eeb259c`.

It completed12 outer steps,132 new parameter vectors and1344 native programs without an error, stopping at `outer_step_limit`. The selected objective is0.4970127555469721, improved from the last iteration origin's0.5295442466796308. Energy checks, every paired-force check and both central-geometry force checks pass. Four individual side-fixture RMS checks remain above the unchanged2 kcal/mol/Å threshold: La minus2.458380, La plus2.140032, Dy minus2.454490 and Dy plus2.320198. Every maximum body-component error is below5. These are remaining development discrepancies, not evidence that extra epsilon variables are already needed.

The selected four-parameter vector is `[3.800640492160888, 0.27995568058029474, 3.490554961303394, 0.2649340807140388]`. The last Jacobian was evaluated at `[3.7919180208770316, 0.27208105028179513, 3.481024470557103, 0.2536840807140388]`. These differ, and no recorded iteration origin equals the final selected point. Thus the saved last Jacobian cannot supply a matched final-point four-column span. A future endpoint may differ again after the separately declared continuation; derivative reuse must be assessed against that exact future authority.

Root's present choice is the unchanged four-parameter continuation. This document remains a future option only; the final17589 vector is recorded as evidence and is **not** assigned as an epsilon-test starting point. Probe inputs, runnable config and execution release remain **unassigned**. This file performs no molecular work and changes no completed protocol.
