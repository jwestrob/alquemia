# Explicit electrostatics with the previously qualified MACE short readout

Declared after FULL_BOUNDARY_GB_REPORT.md and the saved-state coupling audit,
before new native-potential or short-readout outputs. Active-goal authorization
applies. This is a distinct candidate, not a rescue/relabeling of the failed
OMOL-plus-reaction-field model.

## Why this next step

The classical direct field and solvent cross term nearly cancel for alpha,
but the entangled OMOL context has no verified matching electrostatic part.
A beyond36A tail cannot repair that roughly120kcal failure: alpha has no such
atoms and GGR's direct tail is only about-13kcal. Do not add arbitrary Coulomb
energy to OMOL or fit an electrostatic weight.

MACE-POLAR already exposes a trained interaction_energy readout before global
charge restoration/field updates. This is NOT new: the old fixed-field screen
failed a PQQ boundary gate (medium6.60/large5.56kcal), and the short component
alone failed alpha/GGR directions. Preserve those failures. Here, use the
same previously qualified medium readout on matched normalized cores, actual
QM-derived charges, a full fixed solvent boundary and explicit direct coupling.
The short/field split was learned jointly; it is not a unique physical
non-electrostatic decomposition. This candidate must earn its own evidence.

## Fixed energy expression

A_M = E_DFT,vacuum(core,M)
      + C_exact[core_density_M, q_env]
      + G_full[Pq_QM,M + q_env]
      + T_short(full,M) - T_short(core,M).

C_exact=sum_j q_env,j * phi_DFT,vacuum,M(r_j), from the native saved-wavefunction
potential utility. Choose exact density coupling uniformly BEFORE its new
outputs, not whichever charge approximation favors a protein. Keep the fitted,
projected point-charge C from the prior audit as an approximation diagnostic;
it is not a second predictive model or a tunable contribution.

G_full is the unchanged completed OBC2 reaction energy. No coreCPCM, bare
QM-QM Coulomb, MACE total electrostatic/electron scalar, induced-dipole energy,
relaxation or entropy is added. DFT already contains intrinsic capped-core
energy; environment-only vacuum/solvent terms cancel within each Ca/La contrast.
Retain all noncancelling self and cross reaction components. There is no
self-consistent field/density response. Exact direct density plus projected
monopole reaction energy is an explicit approximation; quantify their direct
coupling mismatch rather than claiming a common exact density functional.

R_A=A_Ca-A_La. Convert Hartree/eV/atomic-unit potentials exactly once using
recorded project constants. Larger relative R is more La-like only as a
predeclared direction test. No compatible aquo reference, absolute zero,
baseline band or calibrated affinity class exists for this candidate.

## Data and execution inventory

Use all four existing normalized representations: GGR extended/connected,
alpha1F6S/6IP9. Same exact DFT endpoints, mapped heavy/H coordinates, synthetic
core caps, physical full atoms, peptide connectivity, formal charges, assembly,
protonation,2/3alpha waters and local ff19SB boundary ledger as the completed
full_boundary_GB_v1 manifest. No re-preparation by score or additional cofactor.

Primary MACE checkpoint is medium MACE-POLAR-1-M.model SHA
fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a,
float64, existing exact early-exit adapter polar_scale_shift_exact_readout_v1.
Reuse its qualified native scalar/gradient machinery and numerical reference.
No large-checkpoint comparison or charge-category sweep in this trial.

The read-only field_short_reuse_audit_v1.json verifies six exact whole-protein
short-component receipts from global benchmark job1200701. Their physical
coordinates, charges, spin and checkpoint match.536 named archived core-task
geometry comparisons found zero exact normalized-core matches. Therefore:
-8 new short-component core energy/analytic-force calls (fourCa/La pairs).
-8 native orca_vpot calls at every nonzero q_env site of each actual frozen
 full-boundary state, with paired-identical probe coordinates/weights.
-0new DFT/charge-fit/whole-MACE/GB/training calls. Reuse all saved endpoints,
 eight GB full terms and six whole short terms exactly.

Potential evaluation uses isolated copies of the saved GBW/densities/index;
no new SCF or field-polarized density. Check point/atom ordering and finite
outputs. No probe may coincide with a QM nucleus. Record nearest separation,
including near-boundary probes. This directly tests the region that matters
for coupling, beyond the earlier exterior-only potential test.

## Frozen checks and comparison

Require inherited actual-source/charge/solver gates and exact mapped states.
Compare C_projected versus C_exact for each endpoint and Ca-minus-La contrast.
The paired coupling-error gate is <=1kcal for EACH representation and <=1kcal
for its GGR partition change: half the existing2kcal boundary scale, independent
of labels. Passing the old exterior RMS gate alone cannot satisfy this check.
Do not silently replace the exact term by the point term or ignore a failing
representation gate. Report the raw diagnostic and a failed validity status.

Scalar component/direct score algebra closes within1e-7kcal-scale; all sources
must have accepted receipts. The unchanged short adapter already passed native
component and derivative tests; do not rerun a large qualification panel.
New forces are derivatives of the short component only, not combined forces.

Retain abs(GGR connected-minus-extended R_A)<=2 and all four alpha-minus-GGR
contrasts>0.02. Numerical/charge-representation, partition and biological
ordering judgments remain separate. Report all denominators and failures.
No reweighting, dielectric/radius/cap/water change, favorable checkpoint/core
selection or new threshold after results. These are two consumed biological
groups with qualified cross-study directions, not prospective validation.

Existing allocation/runner policies: oneA5000/16CPU/64474MiB for eight small
short calls; eight one-thread CPU utility workers/16GB for eight potentials.
Expected minutes; no application time/CPU budget. Pin implementations/software,
use locks, preserve partial attempts and record actual costs. Prepare/dry-run
and real-fixture tests first. Baseline/default and concurrent jobs unchanged.
If unsuccessful, retain the result and reconsider the model; do not expand this
candidate's benchmark or automatically switch its parameters to make it win.
