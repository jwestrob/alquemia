# Coupled scaffold ledger implemented; finite-step force discrepancy diagnosed

**The selected classical model and physical maps work, and the apparent derivative discrepancy is dominated by finite-step LJ curvature. Both predeclared finite-difference gates nevertheless remain failed.** No electronic component or full hybrid response was evaluated; this is reusable preparation and classical derivative evidence, not new La/Dy discrimination.

## Actual model and physical motion

The repaired Hans source is represented by1,891 physical atoms:1,887 native protein/water atoms plus actual targetEF3 and La/La/Na spectators. There are191 realQM and1,700MM atoms. Four fixed-length caps exist only in the electronic mapping.

The classical energy is retained native bonded + MM–MM LJ + shifted-charge MM–MM Coulomb + realQM–MM LJ. It keeps all MM bonded terms, four native cut-bond stretches,12crossangles,21crosstorsion entries and4crossCMAP. It excludes whollyQM bonded terms and the8angles/11torsion entries represented by the actual caps. The latter include the doubly cappedCA65–C65–N66–CA66 torsion. Native MM exceptions use charge products rebuilt from the fixed redistributed charges and their verified original1–4 scale. Native LJ exceptions and Lorentz–Berthelot mixing are retained.

Actual installed pure12–6 La/Dy IOD/TIP3P and Na HFE/TIP3P values are pinned in `LEDGER.json`; noC4, capFFparticles, classicalQM–MM Coulomb or invented ion parameters. These parameter choices are available candidates, not a validated embedded-QM cross interaction fit. The omitted quantum component is **unavailable**, not zero. Full electronic+classical energy/force qualification remains false.

The real N83–CA83 graph bridge generates the prescribed backbone rotation:814protein atoms move, including176realQM/638MM. Of these,637 actual field rows move, maximum.02238Å at+.05°. Actual waters, all metals and upstream atoms remain fixed. Bond-length/angle changes are≤2.58e−14Å/4.31e−14rad. Core/field files and complete physical cap mappings exist at all three original configurations.

## Fixed-link checks

All four cap Jacobians agree with coordinate finite differences within6.60e−10. Their radial-MM action is≤3.74e−17, confirming that fixed caps do not supply real cut-bond stretch stiffness. Native retained radial stiffnesses are634 and674kcal/mol/Å²; direct harmonic curvature checks agree within1.87e−10kcal/mol/Å². The phi cap-tangent residual is≤7.97e−7Å/rad. These pass their frozen thresholds.

## Numerical outcome: retain both failures

Per-component La derivatives (Dy gives the same MM terms and nearly identical cross-LJ results), units kcal/mol/radian:

| Component | Analytic | Error at±.05° | Error at±.025° | Coarse/half error ratio | Richardson error |
|---|---:|---:|---:|---:|---:|
| Retained bonded |~0|1.04e−9|3.39e−9|roundoff|4.17e−9|
| MM LJ |−350.583155|+.207168|+.051768|4.001844|−3.18e−5|
| MM Coulomb |−94.986220|+.000144|+.0000364|3.95208|+5.81e−7|
| RealQM–MM LJ |+25.939235|+.012656|+.003164|4.00021|−2.24e−7|

MM-LJ exceeds the unchanged.027012 tolerance at both steps. Every other component passes its original tolerance. The separately named Richardson derivative `(4Dhalf−Dfull)/3` was requested for diagnosis; **it is not a retrospectively invented acceptance gate**. No smaller step or new angle was evaluated.

Independent pair algebra on the saved coarse configurations reproduces the MM-LJ analytic sum and finite difference. The largest error contribution is Met132HB2–water334O: origin2.46037Å/.88179kcal repulsion, swept2.47639→2.44435Å; derivative error+.15246kcal/rad. H40HB3–H113HB2 at2.03811Å adds+.13174, partly offset by other pairs. These are steep but finite contacts, not overlapping nuclei. The~fourfold error reduction and small Richardson discrepancy support O(h²) truncation, not an observed force-sign or unit error. They do not establish the missing complete electronic derivative.

## Execution, failure and writer repair

| Job | Work | Actual wall/allocatedCPU seconds |
|---|---|---:|
|1220320|Three configurations×two targetLJ assignments×four component queries|6/6|
|1220321|Direct pair diagnosis on same saved configurations|6/6|
|1220322|Prescribed halfstep, failed result serialization after16queries|3/3|
|1220325|Identical halfstep receipt recovery, no physics change|2/2|

Total17allocatedCPU-seconds,1sharedCPU per job, mem0, noGPU.56actual component energy/force queries including the repeated halfstep; pair diagnosis adds algebraic pair evaluations but no OpenMM Context. Original run process peakRSS174,728KiB. Slurm's zero MaxRSS samples are not interpreted as zero usage. Serial map generation/reporting adds unmeasured small CPU work. NoQM calls or optimization occurred.

The halfstep's NumPy boolean serialization failure is preserved with its partial output. The identical replay used explicit native-JSON conversion. A new **scoped research-only** `result_io.py` now serializes fully before exclusive file creation, safely converts NumPy scalars/arrays, rejects nonfinite data and preserves no-overwrite behavior. `write_configuration()` accepts actual energy/force arrays for immediate per-configuration receipts before aggregate checks. Four tests on actual saved molecular artifacts pass in.64s. Shared `affordable_common.py` and production writers remain unchanged.

The executed scripts remain immutable historical evidence. Do not rerun their old serializer. Future research execution must use the scoped writer and persist each full energy/force result before calculating aggregate gates. The completed diagnostic stores component energies and projected derivatives; full Cartesian component force arrays were not retained and are not reconstructed or zero-filled here.

## Artifacts and next decision

Workspace `workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1/` contains the pinned `LEDGER.json`, particles/native exceptions/selected bonded terms, exact cross-term decisions, original three physical/core/field configurations, complete cap mappings, `RESULT.json` and `PAIR_DIAGNOSIS.json`. The corrected halfstep result is `halfstep_receipt_recovery_v2/RESULT.json`; the failed first halfstep remains in `halfstep_v1/`. Reports retain all gates and unrounded numbers; `ACCOUNTING.txt` retains actual jobs.

This branch has reached its bounded stopping point. The next substantive test, if root chooses it after current electronic results, is the complete finite energy and its physical derivative on this same motion. Bulk solvent, responsive spectators and QM/MM short-range transfer remain limitations. No optimization, stiffness/entropy correction, production promotion or biological conclusion is justified by this classical test alone.
