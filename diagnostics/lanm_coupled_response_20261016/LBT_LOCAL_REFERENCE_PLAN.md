# Local source-response reference: distinguish electronic models

Declared17October2026UTC after saved comparison12823, before new reference calls.
Protocol `lbt3_local_source_response_reference_v1`. Root owns execution/integration.
This is a finite electronic-response test, not a binding-affinity campaign.

## Question and decision

At identical source geometries and mechanical/solvent terms, replacing MACE with
GFN2 reduces the old full-pool source gap from28.262 to1.173kcal/mol (0.609 at
unchanged selected geometries). This is model sensitivity, not proof GFN2 is right.
The earlier native MACE whole/core control already had similar source-contrast
gaps (72.049/73.690kcal/mol); adding the scaffold alone did not remove the effect.
Before extending structural search, test whether an independent electronic model
supports the large local geometric penalty, the flatter GFN2 response, or neither.

Primary quantities for each method and metal are
`S_M = E_M(7CCN) - E_M(7CCO)` and `S_Dy - S_La` in kcal/mol.
Use same-metal native differences before conversion. No unlike-element absolute
totals, reference offsets or preference thresholds. Report continuous errors from
PBE0 and the already-used1kcal/mol interpretive resolution; it is not a statistical
error bar or biological acceptance threshold. Keep local and whole results separate.
A PBE0 disagreement diagnoses this conditional local surface only. Neither positive
exchange signs nor a flatter source response by themselves qualify a classifier.

## Frozen inputs and reuse

Use all four exact48-atom matched-core XYZs from preparation12469:
`workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/`
`{7CCO,7CCN}/matched_core_{La,Dy}.xyz`.
These contain the target and complete mapped donor-fragment chemistry with the
existing caps, charge-1, targetindex0. Preserve every coordinate, proton inventory,
cap, donor membership and source identity. No new fold, waters or optimization.
Physical states are La(III) singlet and Dy(III)4f9/sextet; record backend effective
states separately. Pin complete XYZ/state/source/implementation/runtime before calls.

Reuse all four matched-core MACE12471 energies and both accepted standalone7CCO
GFN2 vacuum energies from qualification12488. Add only the missing7CCN La/Dy
vacuum GFN2 singlepoints: same published xTB6.7.1 executable/author parameters,
UHF0, charge-1,300K electronic temperature, accuracy0.0001,maxiter500,no restart,
JSON/full populations and exact native convergence/charge gates. No nativeORCA
GFN2 defaults, solvent or regenerated original cells. This two-cell component can
run independently on two shared Slurm CPUs, one thread per actual endpoint.

The historical14-ion LBT calculations cannot supply these DFT references: they
used r2SCAN-3c/CPCM and a different DyECP28 explicit-f treatment; oneDy failed.
Retain all results/failures and do not rerun that campaign. No compatible completed
vacuum PBE0-D4/ECP55 result was found in the bounded inventory.

## Four independent reference endpoints, staged by actual execution

Use the already executable compact frozen-f reference recipe, not a new functional
or electronic-state search: nativeORCA6.1.1 PBE0-D4, CHNOdef2-TZVP/def2-J, metal
lcecp1-TZVP with its pinned ECP and AutoAuxJ, RIJCOSX, DefGrid3, VeryTightSCF,
PModel, NoAutostart, EnGrad. LaECP46, DyECP55; effective restricted valence singlet.
No solvent, external point charges, gCP, geometry optimization or wavefunction seed.
No numerical gradients. Analytic gradients are returned by the established
reference/parser path but primary comparison is energy; no force qualification
is implied by producing a gradient.

Reuse verified primitives `metal_environment_frozen_f.state/input_text/parse`
and the pinned compact execution assets from
`workspaces/metal_environment_response_20260926/lady_compact_exchange_v1/manifest.json`
(SHA25618ed7c0f4141b95cba20191bf36867bc839ae8722945b81dfb74fe3eeeda7f25).
Its old hardcoded Hans/Mex matrix must not be submitted. A new four-row LBT manifest
binds this exact recipe and inputs. State/electron counts derive from each actual
XYZ and declared ECP before launch, and must match observed native output.

First execute La/7CCO on12 actual MPI slots under current shared Slurm policy.
After real convergence/state/gradient admission, execute the remaining La/7CCN
and Dy source pair using the same settings. Numerical admission, not favorable
energy, releases those independent cells. Two endpoints per metal, four total.
Use the existing native convergence declaration and component/parser admission,
finite energies, ECP/restricted-state evidence and matching analytic gradient
energy/coordinates. Record all printed SCF residuals, their criteria, the actual
ConvCheckMode/ConvForced settings and a separate all-residuals-pass flag. The
unchanged compact input uses native total/one-electron-energy convergence; do
not quietly replace its gate with a forced-all-density/orbital protocol. A false
all-residuals flag remains a numerical limitation, never a hidden pass or a claim
of gradient/stability qualification. Preserve warnings and unavailable cells.
Energy convergence is not electronic stability or spin-orbit qualification. No automatic
iteration increase, changed guess, tighter/looser gate or failed-SCF resubmission.

The frozen4f Dy model cannot validate physical sextet/specific spin-orbit effects;
MACE's physical multiplicity input and GFN2's f-in-core treatment are distinct.
PBE0 is an independent approximate reference, not an oracle. A decisive unresolved
method dependence may justify a separately declared small second treatment later;
this plan does not launch a functional sweep or explicit-f recovery campaign.

## Resources and completion

All numerical work goes through Slurm. Normal shared allocations, at most12CPUs/job
and current60/user limit; no memory/time flags, no priority changes. Derive ORCA
MaxCore from the actual grant and match MPI ranks to slots. Record wall/allocated
CPU time, observed memory, every attempted call and all failures. Arm completion
and alert-only SCF-health wakes. No GPU needed. Existing50-atom references took
549s on40CPUs (oneDy) and1809s on40CPUs (five independent compact endpoints);
these are historical measurements, not runtime predictions for the new jobs.

Deliver allfour MACE/GFN/PBE0 rows or explicit missing cells, source strains,
metal-dependent strain differences and numerical/state limitations. Stop this
finite reference question after interpretation; choose the next modeling step
from the evidence rather than automatically enlarging the same search. No new
calibration, production promotion, whole-protein quantum run or library rescore.
