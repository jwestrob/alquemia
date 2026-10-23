# Hydrated acetamide: finite-step curvature or derivative inconsistency

## Frozen question and inventory

The consumed hydrated-acetamide La center has analytic ligand-translation slope14.852538 versus the ±0.15Å secant13.139391kcal/mol/Å. Its1.713148 difference exceeds the existing coarse screen1.485254; that flag remains preserved. This finite test asks whether much smaller displacements support the saved analytic gradient, distinguish finite-step curvature, or leave a numerical energy/gradient inconsistency unresolved. Dy is the matched control. No reference accuracy, force-field fit or physical state is selected here.

Reuse exactly the two31-atom `hydrated_acetamide_central_{La,Dy}` geometries, central energies and full analytic gradients from consumed training result SHA256 `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`. The intact9-atom ligand is indices1–9 (zero based); metal0 and seven waters at10–30 remain fixed. Translate every ligand atom by q times the saved outward unit vector `[0.8660254037844289,0,0.5000000000000167]`, for q=−0.010,−0.005,+0.005,+0.010Å. Both metals use the same coordinates and direction. Preserve all covalent distances, atom order and water/metal coordinates. Save12-decimal XYZ values and check actual displacements within1e−12Å.

Exactly **eight new uncorrected PBE0-D4 single-point energies**:2metals×2step sizes×2signs. No center rerun, gradient recomputation, CC/CP, new geometry search, optimization, parameter fitting or reserved24 outcomes. Parent native results are La`77b64d307a4f6609d43bbd3efc72cfc227668a32965fbc16d9fb518833dee81d` and Dy`459e63b573a78fc68b7fe812ddf834a580057002be78c9fedfd126342e083553`.

## Exact native recipe and resources

Use the hydrated centers' **actual executed** deck, not a substituted Water48 recipe: PBE0 D4 def2-TZVP def2/J RIJCOSX; NoAutostart, DefGrid3, VeryTightSCF; PModel guess and SCF MaxIter300; the same pinned external La/Dy orbital/ECP/auxiliary assets and ORCA6.1.1. Charge+3, effective singlet, frozen-f coresLa46/Dy55 remain unchanged; physical La1/Dy6 labels are separate. Every center has110 explicit electrons and31 real atoms, with no ghost centers. Read and bind the full basis dimension from the original assets/strict native replay rather than treating a label as sufficient evidence.

Preparation must prove the original frozen renderer reproduces each actual executed center deck exactly at its recorded resource values. New single-point text differs only by removing `EnGrad`; runtime PAL/MaxCore are allocation controls. Any scientific setting mismatch blocks launch rather than creating an implicit method comparison. No seed/restart orbitals are reused. Energy-stage differences that appear despite identical scientific settings remain reportable outcomes.

Root submits eight singleton jobs, each4MPI×1thread, standard/standard_expand, scheduler-default memory and time. MaxCore is floor(0.75×actual granted MiB/4). Preserve actual allocation, input, geometry, output, stderr, exit code, elapsed and allocation cost. No job retry, extra step or changed solver is included. Source/preparation does not submit or make molecular calls.

## Admission, observables and unchanged criteria

Reuse the proven Water48/donor strict SP parser for native completion, SCF stop evidence, charge/multiplicity/electron/ECP inventory, all atom coordinates, per-center and total basis dimensions, restricted PBE0/25% exchange, RIJCOSX, real-molecule D4 and electronic+dispersion final-energy closure. No EnGrad file is required for these energy-only tasks. Retain every printed SCF residual separately from the native SOSCF stopping mode; do not retroactively require stale DIIS/density tolerances that were not original admission criteria.

At each h=0.010 and0.005Å report `[E(+h)−E(−h)]/(2h)` in kcal/mol/Å, separately for electronic, D4 and final PBE0-D4 energy. Compare final-energy finite differences to the saved central derivative **−sum(F_ligand)·u**, never a metal-only projection. Report La, Dy and Dy−La results, both step estimates, central analytic values, all differences and the preserved coarse result.

Use Water48's unchanged **0.1kcal/mol/Å** resolution independently for (a) the two-step finite-difference difference and (b) analytic-versus-FD at each step. The smaller-step result is the prespecified summary only when its two-step check passes. No extrapolation, third step, tolerance change or deletion of unfavorable rows is allowed.

Reuse Water48's full energy-stage metadata: COSX grid signatures, final-exchange recomputation, final-integration/count signatures and native excerpts. The center has EnGrad and displaced points have SP by design; do not assume these imply the same effective energy surface. Flag changed signatures; output-stage equality alone does not prove surface equivalence. Missing/admission-failed cells leave a complete eight-slot unavailable denominator and no scientific derivative conclusion.

If both steps agree with the analytic value and each other, the old coarse anomaly is consistent with finite-step curvature on this sampled coordinate; the original flag remains historical evidence. A resolved FD slope disagreeing with the analytic value indicates sampled derivative inconsistency requiring stage/method diagnosis. Unresolved step dependence permits neither conclusion. A successful numerical check does not establish the frozen-f reference's physical accuracy, all-coordinate derivative accuracy or eligibility of unrelated donor states.
