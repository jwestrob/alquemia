# Field-aware Nikasha: response measurable; candidate not yet qualified

27 September 2026. **The native reference detects a stable metal-dependent
response to a real protein-environment perturbation. The intended MACEPOL-EF
comparison is blocked by unavailable exact weights/interface, and native rigid
qualification remains failed.** Do not adopt this path as a preference scorer
or launch coupled LanM relaxation. Preserve the released PQQ workflow.

## Expanded-region result

The completed Thr159 partition test retains the response direction: +0.3528 to
+0.5208kcal/mol. This shows finite representation dependence without a large
jump on this case; it is not partition convergence or improved classification.
[Actual comparison](PARTITION_REPORT.md). Cumulative molecular development work
now42outputs/864128allocatedCPU-s,zeroGPU. All owned molecular jobs are terminal;
current independent transfer-preparation and legacy hybrid checks are in CURRENT.md.

## 27 September continuation

Strict-convergence worker1219450 added4outputs,746s/256624allocatedCPU-s;
all4remain invalid under its frozen printed-residual admission rules. Actual
TRAH switching limits interpretation of those summary density fields. Raw
rotation errors are essentially unchanged: the remedy is rejected, not repeated.
See [SCF diagnosis](SCF_CHECK_REPORT.md). Completed molecular-worker cumulative
cost is now708640allocatedCPU-s across38outputs, zeroGPU; this is development
cost, not routine inference cost.

An independent [partition diagnostic](PARTITION_DIAGNOSTIC_PLAN.md) is being
prepared to compare the same physical hydroxyl perturbation with Thr159 treated
electronically. Known numerical failures remain visible, and no passed
qualification, optimization or classifier follows from this test. Exact ML
weights remain blocked after the public-branch recheck. Check CURRENT.md for
current submission ownership; earlier terminal-state statements below are dated.

## What actually ran

On consumed1H4I, six scout endpoints, twenty force-check endpoints and eight
numerical-diagnostic endpoints completed:34 new native r2SCAN-3c analytic-gradient
calculations. Same54-atom electronic core, Ca-3/La-2 singlets, native La46-core
ECP,9087 fixed protein charges. EnvironmentB rotates actual Thr159 HG1 by10degrees;
this is a modeled perturbation with inherited source O-H geometry, not sampling.

The double difference (LaB-LaA)-(CaB-CaA) is+0.352796kcal/mol originally and
+0.352640 after quadrature refinement. Its stability supports interpreting this
specific electronic response. It does not establish affinity or discrimination.
Selected external-charge and mapped-boundary force finite differences and repeats
pass. Joint rotation fails the declared energy tolerance on both numerical grids;
translation passes. Denser quadrature did not solve that issue. Exact residuals
and all unavailable statuses remain in the linked result records.

## What is and is not supported

1. Public model source supports atomwise potential and field; exact paper weights
and matching engine implementation are unavailable. No alternative model was
silently substituted. Metal/spin field-response quality is untested.
2. Native selected analytic derivatives agree with energy differences. Rigid
qualification fails; constant-potential gauge has algebra-only coverage, and
refined-grid finite differences remain unrun. This is not complete coupling validation.
3. Only one consumed Ca/La core was evaluated. No ML/reference accuracy comparison,
expanded-region qualification,4MAE/non-PQQ transfer, or new La/Dy result exists.
4. Small rotational sensitivity is unresolved. Full hybrid cross repulsion,
dispersion and boundary mechanics remain missing; they were never filled with zero.
5. Next justified step is obtaining the exact model/interface and resolving the
native numerical failure through a specific diagnosis, not further generic grid
increases or whole-protein optimization. Author question is drafted, unsent.

## Costs and usable artifacts

Completed molecular workers:276+455+583=1314 allocation wall-seconds,
452016 allocated CPU-seconds, zero GPU. These include the original scout's
underused344CPU allocation; they are not measured CPU utilization or production
throughput. Collectors and notification tests are separately recorded, outside
this molecular-worker sum. MaxRSS is batch accounting, not verified aggregate MPI
peak memory. Latest worker used8x43MPI ranks with verified exclusive-node RAM.

- [Scout result](RESULT_1219207.json)
- [Directional forces](FORCE_CHECK_REPORT.md)
- [Grid diagnosis](GRID_CHECK_REPORT.md)
- [Upstream capabilities and unsent author question](upstream/AUDIT.md)
- [Checkpoint and recovery](CURRENT.md)

All owned science and collectors are terminal. Actual wake events worked for both
new jobs. No production promotion, remote push, email, library run or old LanM
vacuum retry occurred. The26September queued report is archived separately.
