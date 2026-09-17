# Matched vacuum diagnostic: separate solvent response from learned error

Declared before preparing or executing new quantum outputs. Jacob's active-goal
authorization applies; no additional permission or project compute limit.
Production baseline, references and failed candidate records remain unchanged.

## Question and scope

The failed subtractive expressions combine native r2SCAN-3c/CPCM core energies
with a vacuum-trained learned term. Even an exact vacuum low-level model would
leave a core-dependent solvent contribution in that expression. Quantify this
confound before inventing another learned charge feature or expanding inference.
This is a mechanism diagnostic, not a new aqueous discriminator or rescue of
the rejected static candidates. Removing solvent does not make r2SCAN-3c and
OMOL's omegaB97M-V/def2-TZVPD training Hamiltonian identical.

Reuse the exact eight centers in `shared_neutral_core_v2/manifest.json`
(SHA256 `4b3420ef72fc742b25cc407b1992073e3a55da2d28f25f15c9055781bb7db814`):
GGR_extended, GGR_connected, ALPHA_1F6S, ALPHA_6IP9, each Ca/La. These are
consumed development cases. Preserve every XYZ byte, cap/source mapping,
physical charge (Ca -1, La 0), singlet multiplicity and explicit water.

Run eight analytic-gradient endpoints with ORCA6.1.1 native
`r2SCAN-3c NoAutostart DefGrid3 TightSCF EnGrad`: remove only `CPCM(Water)`
from the actual reference scientific inputs. No numerical DFT gradients,
optimization, Hessian, fitted charge extraction, new functional/basis/ECP,
additional dispersion, solvent solver, model inference or biological scoring.
Protocol: `native_r2scan3c_matched_vacuum_response_diagnostic_v1`.

The existing fixed-field vacuum archive is 1H4I/Asp303 with different coordinates
and partition, and cannot replace these endpoints. Inspect source manifests
for compatible exact states; retain a bounded archive audit rather than claim
an exhaustive filesystem search. Existing CPCM, masked and shared-neutral
energies/gradients will all be reused with their original pins and receipts.

## Output and interpretation frozen before execution

For each endpoint report E_vac and grad(E_vac), E_CPCM-E_vac and its gradient,
and differences versus BOTH already declared learned representations. Retain
physical projections onto the archived metal translation and peptide rotation;
do not interpret independent cap motions as physical coordinates. Export
grad(E_Ca-E_La) without inventing a covariance or response energy.

For GGR compute connected-minus-extended R=E_Ca-E_La in vacuum and CPCM:

    J_CPCM - J_learned = (J_CPCM - J_vacuum) + (J_vacuum - J_learned).

Require algebraic closure within1e-7kcal-scale. Compare the residual vacuum
partition disagreement to the existing2kcal-scale diagnostic tolerance for
both learned descriptors; this does not reopen either failed candidate's
conditional whole calculations. Report every component, including when solvent
worsens agreement. No newly fit threshold or aqueous reference. Eight center
gradients alone do not validate curvature, relaxation, forces against a changed
Hamiltonian or biological accuracy. Missing/nonconverged results stay unavailable.

## Execution and expected scale

Reuse the existing manifest ORCA runner, locks, MPI/thread fixes and allocation
policy. One standard/memory allocation of64CPUs, four concurrent16-rank endpoints;
use all allocated CPUs subject to the finite remaining inventory. Actual matched
CPCM center receipts provide per-task costs; record them in preparation. Expect
minutes per endpoint and order tens of thousands of allocated core-seconds for
this one-time diagnostic, not a promised production cost. No GPU is required.
There is no application time/CPU budget; scheduler limits still apply.

Every output gets an execution receipt and explicit vacuum Hamiltonian check,
including native D4/gCP and La46-electron ECP accounting. Do not pass vacuum
gradients through a parser that labels them `isolated_CPCM_endpoint`. Use the
existing scalar-gradient reader with an explicit vacuum scope. Retain failed
attempts and scratch needed to establish the method; retries require fresh
workspaces with identical scientific input. Do not overwrite completed work.

Tests use pinned real inputs, archived CPCM/vacuum outputs for rejection checks,
and actual new outputs once available. Actual integration remains explicitly
unrun until the scientific executable finishes. No fabricated outputs.
