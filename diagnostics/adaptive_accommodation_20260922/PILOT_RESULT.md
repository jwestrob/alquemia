# Adaptive angular accommodation: physical movement, modest score changes

**Four force-selected donor motions produce valid new candidates, but this pilot
does not demonstrate better classification.** Both known crystal controls retain
their classes under explicit old-band transfer. The two unlabeled PLM examples
remain Ca-supported and inconclusive. Their additional contrast changes are much
smaller than the relief already found by the earlier terminal-donor search.

Both metals competed over the same five geometries per source: original, earlier
Ca/La proposals, and new Ca/La proposals. Chemistry and the solvent-containing
energy expression stayed fixed. No quantum validation, entropy, occupancy change,
new reference calibration or production promotion belongs to this result.

| Source | Additional R versus prior shared pool (model kcal/mol) | Selected Ca / La geometry | Old-band transfer |
|---|---:|---|---|
| 1H4I | +1.339976614 | adaptive La / adaptive La | Ca-supported |
| 4MAE | −1.393052760 | adaptive Ca / adaptive La | La-supported |
| PLM83440678cbbd658047c9 | +2.516050966 | adaptive Ca / adaptive La | Ca-supported; label unknown |
| PLM07ab500e3df76b30d71c | +0.305200530 | adaptive Ca / adaptive La | Inconclusive; label unknown |

Positive changes are more La-like. These are changes of the calibrated
electronic descriptor, not affinity or binding-free-energy measurements. The
two-control separation narrows by2.733029374 model kcal/mol; this is not evidence
of failure across an untested broader population, or an accuracy gain.

## What physically changed

All eight optimizations returned valid candidates. Seven were interior to the
declared final displacement region. PLM8344 La reached0.8Angstrom at an extra-Asp
oxygen and retained a load in the selected space. Omitted motions also retain
forces. This identifies a limitation of the selected local freedom; it does not
establish that scaffold compliance is the unique missing effect.

The broader rotations lower native energies for both metals; most of that extra
relief cancels in their difference. Solvent selection matters: 1H4I's own adaptive
Ca proposal raises its composite energy, while Ca at the La-proposed geometry
lowers it. This actual cross-selection supports retaining the common-pool
primitive even though the pilot does not add a correct biological call.

Intermediate SLSQP trials can violate the final displacement constraint, as
declared before execution. All31 completed infeasible requests are retained;
the largest displacement was3.804Angstrom. Final candidates satisfy the bound,
angle limits, physical atom/cap mappings and overlap checks. No unconstrained
composite minimum is claimed. See [independent physical review](EXPANSION_REVIEW.md).

## Solver failure and explicit numerical completion

The original expanded collection remains **3/4 available**: La at the 4MAE
adaptive-Ca geometry failed the primary vacuum SCF after125iterations. It is not
reported as a successful pool or substituted with an old score.

A separately declared three-call numerical diagnostic changed only the iteration
ceiling to500. The failed point converged at129iterations. Two exact archived
controls reproduced within1.9e−9kcal/mol, with identical parameters/state and
unchanged convergence tolerances. [Diagnostic report](../adaptive_maxiter_20260922/REPORT.md).

`common_pool_v1/recovered_collection_v1.json` is a separately named4/4 result with
numerical policy `primary_native_GFN2_with_explicit_MaxIter500_completion_v1`.
It pins the original failed matrix, successful diagnostic and previous failed
attempt within the replaced component. The three other cases are unchanged.
No source geometry or Hamiltonian was adjusted to obtain a desired class.

## Actual cost and artifacts

- Eight optimizations:114 fresh native MACE energy/force calls,16 cache reuses.
- Common-pool additions:8 fresh MACE calls,8 own-candidate MACE reuses,32GFN2
  attempts (31 initially complete). Numerical diagnostic:3 furtherGFN2 calls.
- Jobs1209857/1209859/1209860/1209876: **9,976 allocated core-seconds and71
  GPU-seconds**, including the failure and both numerical controls. This is
  incremental to existing origins/proposals, not full source-to-score timing.
- Primary collection: `workspaces/adaptive_accommodation_20260922/common_pool_v1/collection_with_reasons.json`.
- Explicit completion: same directory, `recovered_collection_v1.json`.
- Actual proposals: `workspaces/adaptive_accommodation_20260922/proposals_v1/after_proposals_1209857.json`.
- Seventeen shared-pool/recovery tests pass with zero skips; they include the real
  failed result and real numerical completion. [Receipt](RECOVERY_TESTS_v1.txt).

**Next decision:** apply the same fixed selection rule to the designated original
reference population before making a fidelity or promotion claim. The remaining26
sources are being prepared. Retain the released scorer meanwhile; no PLM cohort
rescore is part of this pilot.
