# Direct donor-basin breadth: complete, no scanner gain demonstrated

**Result:** two of four one-coordinate integrals pass the frozen numerical and
domain checks. Their breadth contributions are only −0.107 and −0.084 kcal/mol.
The larger structural effect is in the minimum energy. Adding all grid geometries
to the prior adaptive pools leaves the three available pool minima exactly
unchanged; the fourth pool is explicitly unavailable. No new classification repair,
production change, calibrated threshold or whole-pocket entropy is delivered.

## Actual scope and execution

[PLAN.md](PLAN.md) was fixed before energies: four consumed PQQ sources; one real
source-supported donor torsion selected from the existing common force selector;
65 points,33 nested coarse points, a fixed inner75%domain,300K and identical
Ca/La physical dq/(2π) measure. Protein spectators, metal, PQQ, waters, protonation
and charge remained fixed. All260 geometries passed mapping, physical-bond,
fixed-atom, cap, overlap and0.8Å heavy-displacement checks. A physical atomwise
Jacobian check verified the constant rotor metric; caps are not thermal modes.

Executed512/512 successful native float64 MACE evaluations and1024 native GFN2
ALPB/vacuum attempts, of which1022 converged. Eight exact origin metal-cells were
reused. NoDFT,MD,optimization,Hessian,extra start,training or scientific retry.
The MaxIter500 policy was the already-qualified source-pool numerical protocol;
no tolerance or orbital-start policy changed to repair failures.

The two failures are1H4I La vacuum at basin_56 and basin_59. Actual ORCA outputs
report failure after500SCFcycles. Both remain in the denominator; no interpolation,
artificial barrier, favorable-point deletion or baseline fallback.

## Minimum and width are distinct

All entries below are Ca-minus-La changes relative to the source q0, in model
kcal/mol. Width is `−RT ln(W_Ca/W_La)` and belongs with the minimum of the *same*
65-node curve. It must not be added to an unrelated adaptive-pool minimum.

| Source | Sampled-minimum contribution | Breadth contribution | Total conditional change | Frozen checks |
|---|---:|---:|---:|---|
| 1H4I | unavailable | unavailable | unavailable | Two missing La solvent cells |
| 4MAE | −0.389541 | −0.106738 | −0.496279 | Pass |
| Q9Z4J7 | +4.786034 raw | −1.050726 raw | +3.735308 raw | Fail; correction unavailable |
| A0ACD6B9F2 Ca sample4 | +16.161591 | −0.083615 | +16.077975 | Pass |

4MAE effective widths are0.298247/0.249355rad (Ca/La). A0ACD6B9F2 widths are
0.120466/0.104701rad. Their paired fine/coarse effects are−0.000801 and−0.001387
kcal/mol, and paired outer/inner effects−0.000030 and+0.000695. Maximum outer
weight is0.00512% for4MAE and0.11662% forA0ACD6B9F2, below the frozen1%limit.

Q9Z4J7 fails with0.644739kcal/mol paired fine/coarse change,0.690298 extent change,
and68.6552%La weight outside the inner domain. Its La minimum lies outside that
inner domain. The raw integral remains visible, but the numerical correction is
null. Passing rows mean only **quadrature/domain checks pass**, not an independently
qualified electronic surface, global basin convergence or equilibrium populations.

## Numerical electronic-surface warning

The curves are unsmoothed. Q9Z4J7 La has a5.600864kcal/mol adjacent jump in the
ALPB-minus-vacuum contribution, whereas the native MACE surface varies smoothly.
A0ACD6B9F2 La has a24.713458kcal/mol solvent jump around−0.11rad, outside its
weighted well. Its integral can pass numerical quadrature despite that distant
spike. The attribution is to the observed solvent component; this pilot alone
does not establish the unique electronic state or diagnose its cause.

This is consequential: converged endpoint receipts and a converged integral do
not by themselves validate the physical electronic surface. No jump was fitted,
smoothed or removed. Native, ALPB, vacuum and composite values remain in the
actual collection/result; successful MACE values are retained at failed solvent
points. The figure also makes the incomplete1H4I grid explicit.

## Utility and cost

The two numerically resolved breadth terms are smaller than the existing0.2kcal/mol
composite numerical scale. No demonstrated discriminatory benefit comes from
adding this integral. The largeA0ACD6B9F2 source response is already explained by
lower-energy arrangement, and the existing adaptive pool has an even lower
admitted arrangement: all three available expanded pools retain exactly their
old adaptive minima. These are consumed development examples, not blind affinity
validation.

| Job | Purpose | Allocated wall s | Allocated CPU-s |
|---|---|---:|---:|
| 1210107 | 512MACE cells | 221 | 7,072 |
| 1210122 | GFN shard0 | 918 | 58,752 |
| 1210123 | GFN shard1 | 1,104 | 70,656 |
| 1210124 | GFN shard2 | 611 | 39,104 |
| 1210125 | GFN shard3 | 937 | 59,968 |
| **Total** | | | **235,552** |

GPU allocation:221requestedGPU-seconds. The MACE worker phase, including its
startup, took89.202s;131.798allocated seconds were outside that phase. For the
four solver jobs,99.499/86.239/129.946/137.738s were outside the span between first
native-task start and last native-task completion. These intervals include
preparation, startup and the repeated shared collector's geometry validation;
they are not all uniquely attributable to validation. Allocations were32CPU/
1GPU/200000MiB and64CPU/128GiB per solvent shard, each running eight8-rank tasks.
Reported Slurm per-step MaxRSS is retained; a combined simultaneous-process peak
is unavailable. Initial local preparation/validation, tests and final analysis
are additional unmetered local work; preparation itself recorded264.040s.

Original allocations1210108–1210111 never started and cost zero allocatedCPU-s.
They inherited standard-partition exclusive-node behavior that persisted after
partition/node updates. Following parent coordination, only those confirmed
pending zero-runtime jobs were canceled; the identical manifests were submitted
directly to gpu with explicit nodes. Fresh requests used OverSubscribe=OK and
started on shared CPU capacity. No other jobs, priority settings or running
chemistry were changed. Guarded pre/post/submission receipts are preserved.

## Deliverables and recommendation

- Actual final collection: `workspaces/local_basin_breadth_20260922/scoring_v1/after_solvent_3_1210125.json` (reused directly, no recollection).
- Integral and unrounded components: `workspaces/local_basin_breadth_20260922/RESULT_v2.json`.
- Costs and full accounting: `workspaces/local_basin_breadth_20260922/COST_v1.json`.
- Editable curves: `workspaces/local_basin_breadth_20260922/figures_v2/basin_profiles.svg` and `.pdf`; caption/source pins alongside.
- Nine real-artifact tests pass; see TESTS_v7.txt. These include archived-energy algebra and actual final-grid integration replay, separate from the executed molecular pilot.
- Runnable operations: [COMMANDS.md](COMMANDS.md). Protocol: `physical_single_rotor_conditional_breadth_v1`; shared scoring protocol and production remain distinct.

A reporting-only NumPy boolean serialization bug interrupted RESULT_v1.json;
that partial file is preserved. Analysis implementation_v4 repaired conversion
without changing any value, gate or scientific call. figures_v2 fixes only
presentation; figures_v1 remains preserved.

**Recommendation: stop expanding this breadth pilot.** Keep the reusable physical
mapping/integration diagnostic, prioritize useful geometry accommodation and the
observed solvent numerical problem, and retain the established scorer. The data
do not justify a routine entropy correction or a full-cohort integral campaign.
