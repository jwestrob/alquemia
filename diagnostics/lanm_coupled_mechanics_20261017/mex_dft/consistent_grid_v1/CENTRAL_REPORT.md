# Mex consistent-grid central result

17 October 2026 UTC. Central job **13281 is complete and admitted under the
original native rules**. Root reviewed its measured cost and approved the two
previously declared displaced singlepoints. Directional force consistency is
still untested until both finish and the three-cell comparison passes.

Latest18October: **minus13597 is complete/admitted; plus13596 is still running**.
Minus used32,661s on12CPUs. Its energy-check stopping route, five failing printed
residuals and COSX scalar-closure flagfalse remain recorded in
`displacement_submission_v1/minus/result_v1/`. No paired derivative result yet.

Both singlepoints launched with12MPI each and completion/health notifications.
Config84a47050 and the actual
launches are archived in `displacement_submission_v1/`; no waiting collector
allocation was submitted.

## Actual result

- Normal ORCA termination, 46 SCF cycles, stopping through the gradient check.
- Complete 209-atom QM gradient, 1,840-row external-charge gradient and mapped
  2,049-atom physical gradient. Exact source/state/field joins pass.
- Grid2 and grid3 match in every parsed fingerprint, including 1,026,059 points;
  the original XC and grid1 settings remain unchanged.
- La ECP46, charge0, effective singlet and832 explicit electrons remain unchanged.
  Execution used12 MPI ranks,48,000MiB allocation and MaxCore3100MB per rank.
- The final printed COSX correction is −0.000000027Eh. Both recorded scalar
  closures pass their printed-rounding checks. This supplies no force-error bound.
- Saved O97 electronic gradient projection is −3.6736455776498964kcal/mol/Å.
  This is dE/dR, not force; it has not yet been tested against the new energy slope.

Only the orbital-gradient residual meets its printed threshold. Energy change,
both density changes, orbital rotation and DIIS error do not. These observations
remain explicit under the original native admission; they are not new criteria
introduced after inspecting the result. No electronic-state stability, global
force accuracy or finite-step truncation qualification is established.

Result SHA256:
`dbf29cb3d3b8948630f71b77d924f39630e6a399ba8a21f59bfdd5f118acb27f`.
Full result, assembly, root/independent reviews and actual terminal-health delivery
are preserved in [central_result_v1](central_result_v1/ROOT_REVIEW.json).

## Cost decision and bounded continuation

Actual scheduler elapsed was41,475seconds on12CPUs: **497,700allocatedCPU-seconds**.
This is11hours31minutes15seconds of wall time. Allocation accounting is distinct
from measured utilization. The native runner's parent-only CPU measurement does
not include its MPI children.

Root accepts two further12MPI singlepoints to complete this finite test. Both use
the exact originally prepared ±0.001Å external-O97 displacements, unchanged core,
state, method, SCF settings and consistent-grid recipe. The original acceptance
rule remains `0.005 + 0.0001*max(abs(FD),abs(analytic))`kcal/mol/Å. Existing MM and
ALPB terms are reused; no new such calls or optimizer steps are allowed.

At most three fresh molecular calls in total: the completed central and these
two singlepoints, with no retries. Scheduler memory/time defaults apply. Current
standard policy permits12CPUs/job,120/user,10jobs/user and a3-day maximum; no
limit, priority or other-user allocation is changed. Both runs must reproduce
the central allocation-derived numerical settings. A source, state, runtime or
native admission failure prevents a qualified comparison.

The exact decision and planned commands are in
[displacement_submission_v1](displacement_submission_v1/ROOT_CONTINUATION_DECISION.json).
Read CURRENT.md and actual submissions for launch state; preparation is not
evidence of execution. Original failed qualification and its remaining six/eight
dependent cells remain stopped on either outcome. No production or preference
claim follows from this central result.
