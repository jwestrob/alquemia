# PedE solvent response — declared 13 October 2026

## Authorization and question

Jacob: “Continue your investigation however you see fit. Let's fix our solvent
correction.” This authorizes contained calculations and supersedes the preceding
saved-artifact-only pause. Root owns submission and integration. Production,
thresholds, earlier failures and unrelated work remain intact.

Does the large structure-dependent GFN2 solvent response reflect a credible
electronic response, a solvent approximation, or a preparation/state mismatch?
The completed diagnosis found a 21.1062 kcal/mol origin-score difference between
PedE sample3 and sample1: −0.1504 MACE plus +21.2566 solvent. This motivates a
matched reference, not subtraction of an inconvenient term or label fitting.

## Fixed inputs and parallel work

1. Decompose the saved solvent transfers for all five two-Ca PedE origins and
   all four compatible La-conditioned controls. Retain electrostatic, surface,
   hydrogen-bond and reference contributions plus the intrinsic remainder.
   Missing printed entropy stays missing. These are saved outputs, not new
   molecular evaluations.
2. Audit PQQ connectivity/proton inventory/source mapping and the actual fold
   template. Compare archived crystal geometries with their state limitations.
   A neutral fold template versus an oxidized PQQ3− scoring state is a question,
   not proof of an incorrect structure. Do not change protonation or redox state.
3. Prepare eight singlepoints: saved PedE two-Ca sample1 and sample3 origins ×
   Ca/La × vacuum/CPCM water. Preserve all 192 atoms, 14 caps, coordinates,
   zero explicit waters, Ca−3/La−2 physical singlets. Exact inputs and source
   maps are pinned in the reference manifest before execution. No optimization,
   new folds, numerical gradients or new candidate geometry.

## Independent reference and execution gate

Use installed native ORCA 6.1.1 r2SCAN-3c / DefGrid3 / TightSCF, built-in composite
basis/ECP/corrections, with CPCM water only in the water partner. No manual extra
dispersion/gCP or xTB parameter import applies to this DFT Hamiltonian. Verify
actual basis, La ECP46, electron inventory, successful SCF and normal termination
from output. This is a new diagnostic protocol; historical NormalSCF results are
not overwritten. Check final TRAH residuals where used, not stale DIIS rows.

Start sample1 La water alone, measure convergence/cost, then its vacuum partner.
Expand to the other six declared cells only after the scout is executable and
chemically interpretable. Preserve GBWs and populations for bounded recovery;
failed energies cannot substitute for accepted endpoints. Additional numerical
repairs need a diagnosed cause and separate record, not repeated blind retries.

All scientific code/tests use Slurm. Current shared QoS permits 12 CPUs/job and
60 CPUs/user. Use 12 MPI for the scout, native MPI fixes and allocation-derived
MaxCore with 25% operational headroom. No memory/time overrides, exclusivity
bypass, priority changes or production dependencies. Cheap audits use one CPU.
Do not cancel existing jobs2669/2832/2833 on the inconsistent node. Arm explicit
completion wakes on the current root thread and record jobs/commands/receipts.

## Comparison and decisions frozen before new energies

For each method, U_M(q) = E_water,M(q) − E_vacuum,M(q).
Primary diagnostic D = [U_Ca(1) − U_La(1)] − [U_Ca(3) − U_La(3)].
The current GFN2 value is +21.2566283 kcal/mol. Also report each metal's source
response, vacuum response and within-method group-charge changes. Do not compare
absolute charges from different population schemes as identical observables.

DFT/CPCM and GFN2/ALPB differ in both electronic and solvent models; neither is
an oracle. Compare saved ALPB non-electrostatic terms explicitly. A discrepancy
above 2 kcal/mol is consequential relative to the existing 3.4722 kcal/mol
calibration gap; this is a mechanistic triage scale, not a fitted acceptance rule
or a claim of DFT accuracy. Below that scale, preserve method/numerical uncertainty
and avoid diagnosing a unique cause. If numerical uncertainty could change this
judgment, qualify a small decisive subset rather than asserting certainty.

No affinity or classifier decision follows from eight reference endpoints. No
old aquo offset/bands transfer to a changed Hamiltonian automatically. Choose an
affordable corrective experiment from the physical findings, then use the
designated canonical calibration and all declared structural controls. Do not
remove terms, select structures, clip charges or alter thresholds to rescue PedE.
An unresolved or negative finding stays visible; no automatic promotion.

## Ownership

- Root: plan, scheduler, comparison, shared status and any promotion proposal.
- expanded_launch_review: saved solvent accounting/decomposition only.
- claude_task_review: source/cofactor geometry and state review only.
- la_parameter_scope: pinned eight-cell reference preparation/runner only.

Each agent owns its diagnostic subdirectory. Existing runners are reused through
immutable snapshots where needed, not edited in another agent's working files.
