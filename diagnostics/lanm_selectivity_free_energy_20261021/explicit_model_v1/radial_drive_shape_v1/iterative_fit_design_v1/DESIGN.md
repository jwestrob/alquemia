# One bounded radial14 calibration — decision draft

**Source-only advice.** No executable campaign, optimizer call, native call, new reference or heldout access is created or released. Root owns the declaration and resource envelope. Exact authorities are in `INPUTS.json`.

## Evidence and scientific purpose

The admitted local comparison (`6163eb69`, affine `81a04e9d`) and two measured half-points (`aada89fd`, admission `594e721d`) support testing a bounded joint calibration:

| Measured state | q | P RMS | B RMS | pair dipole RMS (D) | C dipole RMS (D) | susceptibility RMS (a0³) |
|---|---:|---:|---:|---:|---:|---:|
| feca, e=0 | .22727182 | 12.72032 | 10.44183 | .63943155 | .08584593 | .68044364 |
| original12 half | .20920128 | 12.04232 | 10.20113 | .63942965 | .08562692 | .67259677 |
| radial14 half | .18343700 | 11.06138 | 9.55234 | .63942730 | .08276691 | .65820174 |

P/B RMS units are kcal/mol. Both measured half-points pass all25 physical gates, three norm safeguards and signed-W9 safeguard. The radial half's q differs from its affine prediction by +.00060310; the original12 difference is −.00002950. Both local solver returns remain **AlmostSolved**, not established optima. This supports a useful sampled direction in the augmented family; it does not establish global superiority, convergence or physical adequacy.

The two added component-Jacobian columns are exactly zero in this saved central panel. At r0 the new radial drive vanishes but changes radial forces; its joint movement with the original parameters can relieve the force/component conflict. H is unchanged, so e cannot change susceptibility at fixed original coefficients. The W9 z contrast remains wrong-sign: radial half model +.12786475 D versus reference −.10461795 D, residual +.23248270 D. Nonincrease is a safeguard, not electrical acceptance.

## Fixed policy

Start at exact original `feca1810` artifact `d7ce78d0` plus `[0,0]` (14-view key `3c7547e9`). Do not substitute the measured half as the initial incumbent. Keep the original12 domain, append La/Dy e in [−1,1], unchanged radial-qualified vacuum Hamiltonian, targets, 46 dynamic cells and48 static terms.

Keep q = max(J_P / 711.9511599971815, J_B / 497.4305961564669). These original normalization authorities are distinct from the new start and safeguards. Do not reset them at feca or subsequent incumbents. Decreasing max-q does not require both component scores individually to decrease; report both.

Every accepted incumbent must satisfy the unchanged25 actual physical gates **and** all three actual electrical scores at or below their immutable feca values: pair_mu .40887270287953875, C_mu .00736952423619518, chi .46300354382546827. Also require absolute signed W9 Dy−La lab-z model-minus-reference error ≤ .23251750527198523 D. Report old_EF and all electrical losses, signed values and units; no old_EF aggregate ceiling and no new physical electrical tolerance.

For exact saved first-conic reuse, bind the literal admitted LINEARIZATION denominators: pair_mu .10221817571988469, C_mu .001842381059048795, chi .11575088595636707. They derive from half-weighted531 residual squares (score/4), with possible last-bit arithmetic differences from stored scores. Keep these numerical denominators fixed at every subsequent linearization; keep the actual score tests above separately exact. Keep the signed-W9 reference −.10461795254288386 D and ceiling fixed. Neither `normalized_problem`'s current-norm reset nor `MATH.problem`'s current-W9 reset may survive in the iterative variant.

At each of at most8 steps, use minus/plus probes for all14 normalized coordinates at ±.005, clipped to the original box, divided by their **actual** spans. Radial physical probes are e±.01 away from boundaries. Preserve the initial Dy_v clipped span .005170593745372654. Trust radius is .025 in each normalized coordinate. Fractions remain [1,.5,.25,.125].

## Minimal implementation changes

1. **Numerical wrapper:** reuse radial MATH dimensional views and radial SOLVE case B (three electrical cones plus signed-W9 sides), with an explicit general14 probe helper: the present helper only accepts e=0. Construct the immutable context from original normalization authority plus feca safety authority; existing NEW.make_fixed cannot simply receive feca because it validates origin24. Reuse current linearize/scalar certificate machinery. Radial14 retains20 SOCs,890 affine sides,1437 conic rows/15 unknowns; its full scalar check has910 margins. Do not reuse NEW's four-cone removal formulation.
2. **Proposal import:** bind the saved radial14 x15, ordered28 stencil, exact matrix/layout, coefficients, denominators, target and settings pins. The saved full/half witnesses remain immutable; construct all four fractions using the same scalar logic. Import the first solution with0 new solves only if these identities match exactly; fail preparation on a mismatch rather than silently rerun or approximately reuse it. A root-declared fresh-first-solve variant would require its own explicit policy. Subsequent steps allow at most1 Clarabel solve each: unchanged settings `7794882e`,200 iterations,60s,1 thread, no retry or alternate solver.
3. **Witness and selection:** a finite AlmostSolved return can supply a candidate direction only when independently reconstructed retained scalar margins, literal box/trust, all25 predicted gates, fixed norms and signed-W9 sides are nonnegative without the descriptive1e−8 feasibility relaxation, and predicted raw component q improves. Preserve status, residuals, duals and solver_success=false. Interpolate q from t², not t; preserve the existing auxiliary q=predicted max-component ratio+1e−8 construction, which changes no physical parameter or gate. Select only after every eligible declared fraction is completely native-admitted. The actual all25+three scores+W9 tests and raw q reduction decide acceptance; retain current-based epsilon1e−10*max(1,abs(q_current)), current-wins ties, then earlier fraction. Probes never become candidates except through a separately certified exact trial role.
4. **Execution/persistence:** adapt NEW's complete-batch driver and journal, radial render/parser/effective-per-metal identity, and v2 terminal persistence. Replace fixed four/two-point allowlists and112 limits with declared probe/trial roles and limits below. Preserve truthful resources, per-attempt fsynced intent before launch, return/cell evidence, lock ordering, whole-vector budget precheck and vector/failure/final ledger materialization. Carry v2 primary-error preservation and secondary snapshot/stop error reporting; do not restore the stale-terminal-result failure. No native kernel, method, parser tolerance or induction policy changes.

## Exact reuse and maximum new work

Import one center plus28 derivative vectors: the24 exact old zero views under bridge `d965f122`, and4 already admitted radial probes. Preserve old producers, raw records, costs and actual spans. The bridge does not authorize importing the entire historical old-binary registry. Same-metal cache identities must include all six original physical values, e, input bytes, native source/executables and parser authority. Other-metal cells are reused only under exact scientific identity.

| New work within8 steps | Vectors | Maximum native programs |
|---|---:|---:|
| Steps2–8 derivative stencils | 7×28 =196 |196×28 =5488 |
| Four trial slots at each step |8×4 =32 |32×56 =1792 |
| **Hard maximum** | **228** | **7280** |

There are at most8 logical proposal records,1 imported solve and7 new solver calls under exact first-proposal reuse. No new first stencil or center is needed. The measured radial half `869cc32e` can serve its explicit first half-trial slot after exact identity/admission checks, saving1 new vector/56 programs (all-other-work upper bound227/7224); savings never expand either hard cap. The measured original12 half is a comparator unless an exact separately declared radial14 trial coincides. Record inherited full vectors, inherited per-metal cells, requested roles and actual native attempts separately. No broad historical re-admission or repeated completed programs.

Stop on no eligible witness, no actual accepted improvement, a hard cap, native/admission failure, or8 steps. Preserve the last fully accepted incumbent and all failed/incomplete fractions; no partial-batch selection, automatic retry, extension, or second campaign. An incomplete derivative panel cannot support a proposal. Existing development milestones remain descriptive unless the root declaration explicitly adopts a stopping criterion.

## Decision and necessary preparation checks

Recommend this one finite development calibration: measured radial flexibility improved the same consumed objective while retaining the proposed safeguards. Its purpose is to test whether that improvement persists after relinearization. No additional reference job is needed. The consumed reference definitions and uncertainty limits remain: fixed common-basis electronic paired targets, electrical state/axis conventions and conditional susceptibilities; cluster CC and D4 are separate evidence and must not replace target rows. No result qualifies binding free energies, transfer, periodic dynamics or a physical electrical model.

Before any release, test only the policy changes: probes at nonzero e and boundaries/actual spans; frozen versus current ceilings; exact first-proposal import; strict AlmostSolved witness versus infeasible return; q-cushion isolation; actual safeguard rejection; incomplete-batch retention/ties;14/per-metal cache identities and zero-bridge scope;228/7280 intent accounting; v2 failure persistence. Reuse existing journal controls and native qualification rather than repeat them. Root must still declare resource/release bindings and terminal review; a consumed-trimer followup remains separately bound to the final admitted candidate.
