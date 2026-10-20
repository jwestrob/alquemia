# Two native directional checks submitted

The completed whole-peptide reference supports the recent La energy lowering,
but its remaining forces disagree with the model-converged endpoint. This
finite test checks the PBE0 reference's own energy/force consistency before
using those forces to diagnose the mechanics model.

## Current execution

| Stage | Job | State at submission | Resources |
| --- | --- | --- | --- |
| Saved-input preparation | 15845 | Completed and reviewed | 1 CPU,6 seconds |
| Negative displacement | 15847 | Running | 12 MPI CPUs |
| Positive displacement | 15848 | Running | 12 MPI CPUs |

Both native jobs use `standard_expand`, scheduler default memory/time and the
original allocation-derived MaxCore policy. Each has its own completion wake
and alert-only health watcher. Both startup checks pass with48000MiB and MaxCore3100MB/rank. Completion and
health watcher processes are live. Final native admission remains pending;
startup is not a scientific result.

The two positions are exactly `x±0.01*u` Å around step725, where `u` is the
normalized full256-atom MACE-minus-PBE0 vacuum force discrepancy. The PBE0
method, basis/ECP, charge−2, singlet, SCF/grid and unseeded guess remain fixed.
Only EnGrad is removed from the native deck because the new calls need energies.
The existing central analytic gradient is reused.

Preparation15845 passed the source, coordinate, state and zero-call checks.
It consumed6 allocated CPU-seconds and made no molecular calls. Both independent
reviews pass. Stored PBE0/MACE projections along `u` are−66.4531999104 and
−9.02704252585kcal/mol/Å. These are saved analytic projections, not a new
finite-difference result.

## After both endpoints finish

Inspect each actual native result/SCF diagnostics, preserve failures and costs,
and bind their exact result hashes for oneCPU saved collection. The collector
compares the centered energy difference against the central PBE0 projection,
using Decimal subtraction of printed energies and the frozen numerical rule.
Do not run a waiting collector allocation or duplicate either endpoint.

A pass supports this one directional derivative at this step size. A failure
cannot uniquely identify finite-step, grid, SCF or electronic-state effects.
Neither outcome establishes the physically correct Hamiltonian, a complete
force qualification, solvent accuracy or metal preference. MACE's training
Hamiltonian differs from this comparator. No automatic retry, grid sweep,
central rerun, additional optimizer or other-metal calculation follows.

## Evidence

- [Declared protocol](SPEC.md), SHA `c79d6737`.
- [Method comparison](HAMILTONIAN_CONTEXT.md).
- [Independent source review](STATIC_REVIEW.md), SHA `64951613`.
- [Independent actual preparation review](ACTUAL_PREPARATION_REVIEW.md), SHA `ad55da10`.
- [Preparation archive](preparation_result_v1/ARCHIVE.json), SHA `f8d31a51`.
- Preparation RESULT `779412e9`, MANIFEST `585fef98`, DIRECTION `f20be914`.
- Source `cd4396c0`, prepare config `cec2f6e3`, run config `dd4875fa`.
- Root execution decision `62cd2ead`; only the two named singleton calls released.
- [Completed central reference report](../lbt_7ccn_converged_reference_v1/FINAL_REPORT.md).
