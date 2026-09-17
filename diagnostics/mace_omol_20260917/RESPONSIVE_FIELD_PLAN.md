# Responsive quantum core, fixed classical environment, local MACE context

Declared after EXPLICIT_FIELD_SHORT_REPORT.md and before any new embedded
endpoints. Active-goal authorization applies. Distinct version; frozen-density
candidate remains failed. No new threshold, favorable structure selection,
force-field/radius/water modification or checkpoint comparison.

## Question and precedent

Does allowing the core density to respond to the actual permanent protein field
supply useful missing information after the direct-charge representation has
passed its relevant coupling test? Recompute the reaction term from the new
density instead of attaching an old frozen-density term to a new quantum state.

This is not a newly discovered method. The older1H4I qm33/qm36 native field
experiment is in diagnostics/density_embedding_20260916/REPORT.md. Its endpoint
responses were about-31kcal, but differential partition effect only-0.974kcal;
it did not solve the remaining boundary error. That different PQQ preparation
cannot answer the present GGR/alpha case. Reuse its native point-charge input,
energy-accounting and receipt infrastructure, not its scientific outputs.

## Expression and fixed inputs

For each M:

`B_M = E_DFT,embedded(core,M; q_env) + G_full(Pq_embedded,M + q_env) + T_short(full,M) - T_short(core,M)`.

The native embedded energy must include core internal and direct permanent
core/environment interactions exactly once. Verify that expression against
installed ORCA documentation and actual native input/output before launch.
Do not add C_exact again. Fixed classical environment-only vacuum terms cancel
between Ca and La; retain the full reaction term and its self/cross components.
The quantum density responds to permanent protein charges; it does not respond
self-consistently to the added GB reaction field. This is an explicitly
approximate one-way coupling, not self-consistent QM/GB or a free energy.

Reuse all four normalized representations and the exact nonzero q_env from
explicit_field_short_v1. Same native ORCA6.1.1 r2SCAN-3c, DefGrid3, TightSCF,
NoAutostart, vacuum; only add the declared point-charge file. Paired identical
coordinates, charges, assembly, caps, water inventory and multiplicity. Preserve
native basis/ECP and composite corrections. No CPCM, geometry optimization,
trajectory, numerical gradient, response clamp or added empirical energy.
Use native single-point analytic-gradient endpoints as supported by the existing
runner; resulting derivatives belong to the embedded quantum energy only.
No combined-model gradient or mechanical correction is supported.

CHELPG stays the previously frozen native scheme and defaults, with the same
source/cap projection. Full ff19SB environment and charge redistribution stay
fixed; do not reassign them to suit the polarized density. No charge
renormalization. Retain endpoint charge precision and closed-shell/ECP checks.
Full physical GB boundary, radii, OBC2, dielectrics and salt unchanged. Reuse
both full and core MEDIUM short components from the completed explicit-field
trial: geometry/model has not changed.

## Finite execution inventory and resources

-8 new native embedded high-level endpoints, oneCa/La pair perrepresentation.
 Reuse eight vacuum energies and exact vacuum density couplings. Existing
 quantum runner/allocation policy, up to64CPUs; no applicationCPU/time budget.
 Recent matched eight-vacuum-endpoint job took751s and48064allocatedcore-s;
 older field job took1259s for4 endpoints on different hardware. These are
 cost evidence, not a guaranteed estimate. Record every attempt.
-8 native CHELPG utilities from the new saved densities.
-8 native potential utilities, each at the union of the existing exterior
 quality probes and all nonzero environmental charge sites. Keep both probe
 sets identified; no probe/nucleus overlap. One8CPU/16GB utility allocation.
-44 new cheap full-boundary GB calls:8full,8QM,8identity,8rotation,8translation,
 2GGRconnected repeats,2GGRconnected Reference. Reuse4environment-only terms
 with identical coordinates/charges/cavity. Current full48call job took57GPU-s;
 preserve existing GPU16CPU/64474MiB runner policy.
-0newMACE/whole inference/training/geometry-search calls.

Prepare/test/dry-run before submission; preserve wavefunction/density/index
and copied attempts. Failed convergence must remain explicit. Do not expand
this inventory automatically or secretly attach vacuum-fit GB when embedded
charge extraction fails. Prepare any later broader benchmark separately.

## Checks frozen before outputs

All inherited real-source, coordinate, charge, state and execution checks apply.
Keep CHELPG exterior RMS<=.005au OR relativeRMS<=10%, for endpoints and paired
potentials. Projection conserves charge<=1e-9e and dipole<=1e-8eA. Native printed
charge closure<=5e-5e; no renormalization.
Projected-minus-exact direct Ca/La coupling error<=1kcal EACH representation,
and absolute GGR partition error<=1kcal, as in the prior near-boundary test.
These assess the updated distribution; don't reuse a passed vacuum fit gate.
GB identity<=1e-6kcal, rotation/translation/repeat/Reference<=.01kcal;
component algebra<=1e-7kcal-scale. No score-dependent tolerance changes.

Record the electronic-response diagnostic
`P_M=E_embedded-M - E_vacuum-M - C_exact[vacuum_density_M,q_env]`.
Each P_M should be<=+.05kcal within the same stable variational state. A larger
positive value requires Hamiltonian/SCF accounting review; do not clamp it.
Also retain intrinsic and direct changes using the embedded native density
potential, and the change in each GB component. No unique causal attribution
from total energies alone; explicit environment polarization remains absent.

Same predictive tests: abs(GGRconnected-minus-extended contrast)<=2kcal;
allfour alpha-minus-GGR contrasts>.02kcal. Report physical/representation,
partition and ordering outcomes separately. Two already consumed biological
groups; alpha cross-study conditions unresolved and twoalpha water counts differ.
No fresh-blind claim, calibrated reference, old band, universal zero, fitted
weight or selected successful representation. Absolute class, combined gradient
and relaxation remain null even if these development gates pass.

If this model fails, retain its outcome. Do not widen this candidate or tune
its boundary/solvent on these outputs. Any different physical approach needs a
new version and predeclared question; the larger authorized goal continues.
