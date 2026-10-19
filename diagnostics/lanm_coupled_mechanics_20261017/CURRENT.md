# Coupled La/Dy mechanics: current checkpoint

17 October2026UTC. The PLAN remains authorized and unmet. Root owns all launches and integration. PQQ production stays unchanged.

## Exact checkpoint continuation submitted after actual preflight

19 October23:39UTC: preflight15192 passed exact original checkpoint/cache,
physical-source guards and version replay; zero native calls/model loads.
Result2dd9db66 is archived in lbt_7ccn_continuation_scout_v1/preflight_result_v1;
actual cost10CPU-s. Independent actual review16b664a9 confirms all bindings.

Root decisiona21d2c53/config7127a7fc release only the declared7CCN-La scout:
12 additional accepted steps/40 new distinct evaluations, cumulative92/120,
unchanged original80/80 history/trust/source/energy/force criteria. Job15254
is submitted with explicit standard_expand,12CPUs/oneH200, scheduler defaults.
Completion watcher3497887 and exact command are in the run workspace and
[archived submission](lbt_7ccn_continuation_scout_v1/run_submission_v1/run/SUBMISSION.json).
Read actual run results before subsequent work. No converged minimum, corrected
force accuracy or metal preference is established. Original artifacts and all
closed branches remain preserved. Goal remains unmet.

## Alternative QoS accepted for the pending preflight

19 October23:22UTC: Jacob relayed Rohan's instruction to try standard_expand.
Live policy permits24CPUs/job,500CPUs/user,20jobs/user,5days; user association
allows it. Root successfully updated owned pending preflight15192 in place:
QOS=standard_expand, unchanged1CPU/qualified host/job ID/watcher. Slurm
changed the default time limit3→5days. Both QoS priorities are0; scheduling
availability is not guaranteed. [Exact record](lbt_7ccn_continuation_scout_v1/qos_expand_v1/SCHEDULING.md).
Use explicit `sbatch --qos=standard_expand` for the already-declared12CPU/oneH200
scout after actual preflight admission. Pinned science, old launchers and original
failures remain unchanged. No new molecular run has been submitted yet.

## Checkpoint continuation preflight submitted

19 October22:57UTC: independent static review1813ff0a passes the exact
7CCN-La restart adapter0c14a03e/configad0d0996. Root submitted saved-only
preflight15192 on1CPU/noGPU to the qualified H200 host. Actual submission,
review and completion watcher3437814 are in
workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_continuation_scout_v1/preflight.
Read [the specification](lbt_7ccn_continuation_scout_v1/SPEC.md).
No new molecular execution has been submitted at this checkpoint. The actual
passed preflight and pinned root decision must precede the single12-step/40-call
scout. Existing inputs, original80-step result and production remain preserved.

## Saved direction supports a bounded continuation scout

19 October22:44UTC: saved audit15176 completes on1CPU/2s, zero molecular
calls. [Report](lbt_7ccn_direction_audit_v1/REPORT.md), resultcc102d11.
The direction the unchanged optimizer would form from its80-step7CCN-La
checkpoint, if continuation were permitted, has nominal dE/dalpha
MACE−0.647635/PBE0−0.744240kcal/mol. No fallback or trust scaling occurs.
This differs from the older deformation chord and supports local descent
under both saved gradients; finite-step lowering/force accuracy remain unproven.
The original optimizer stopped before forming that direction.

Root is preparing one separate7CCN-La scout:12 additional accepted steps or
40 new distinct evaluations, preserving exact history/trust/initial point,
all physical guards, energy terms and original force criteria. Cumulative
limits92/120 bind actual80/80, without changing the original80/200 artifacts.
Versioned exact-checkpoint budget admission and actual zero-call Slurm preflight
must pass before that single run. No other source/metal campaign, new DFT,
production change or natural-transfer conclusion is released. Distributed force
errors and full reference caveats remain. The broader goal remains unmet.

## Four 7CCN references and comparison complete

19 October22:37UTC: after-La14880 completed and collector15173 re-admits all
four actual results, zero new molecular calls. Read
[the final report](lbt_whole_7ccn_reference_v1/FINAL_REPORT.md).
WholePBE0+ALPB deformation work is La−131.838380/Dy−102.738676kcal/mol,
versus MACE−137.647592/−103.506621. Prior broad accommodation is supported,
but La error−5.809212 and differential error+5.041267 exceed the original
1kcal reporting resolution. No preference or converged minimum is established.

Terminal La chord force is PBE0−1.716526 versus MACE+1.580140; full force-error
RMS3.574975 remains. The chord is not necessarily the next optimizer direction.
Root is preparing a separate saved-only exact checkpoint-direction/force-error
analysis before any continuation. No new molecular calculation or optimization
follows automatically. Native residual/state/force-qualification limits and both
missing Dy model-force comparisons remain. Full resultfff95c2f and evidence are
in lbt_whole_7ccn_reference_v1/full_result_v1; total1,179,035allocatedCPU-s
including preflight/collection. All four native allocations are terminal.
The broader mechanics goal remains unmet. Closed branches and production stay
preserved. This supersedes the pending endpoint statuses below.

### Earlier three-cell checkpoint

19 October21:54UTC: after-Dy14882 is COMPLETE/admitted under original native
gates, result1aca406f/native receipt061ffabf. Full256 gradient retained;
native17-cycle convergence still has nonuniform printed residuals. Physical-state
and force qualification remain false; Dy MACE force comparisons remain unavailable.
Actual24,232s/12CPUs =290,784allocatedCPU-s. Exact evidence/cost/review:
[after-Dy record](lbt_whole_7ccn_reference_v1/endpoint_reviews_v1/after_Dy/ROOT_REVIEW.json).

Three of four cells are admitted. After-La14880 remains RUNNING in analytic
gradient calculation. No recovery is needed. Preserve its completion/health
wakes; collect the four actual results only after it finishes. No partial
paired-work result, optimization continuation or preference conclusion follows.

### Earlier before-Dy completion

19 October21:48UTC: before-Dy14881 completed and passes original native admission,
result391e1d99/native receiptd1ca4271. Full256 gradient retained. Native convergence
after17 cycles has nonuniform printed residual flags; physical-state and force
qualification remain false, and archived Dy MACE forces remain unavailable.
Actual23,858s/12CPUs =286,296allocatedCPU-s. Exact raw output, gradient, result,
cost and review are preserved in
[endpoint_reviews_v1/before_Dy](lbt_whole_7ccn_reference_v1/endpoint_reviews_v1/before_Dy/ROOT_REVIEW.json).

Together with14795, two of four reference cells are admitted. After-La14880 and
after-Dy14882 remain RUNNING in analytic-gradient calculation at this checkpoint;
no recovery or duplicate calculation is needed. The first immediate result-file
check preceded filesystem visibility; reread found the complete successful
result. Preserve that observation without inventing a wrapper failure. Wait for
both actual after-state results before the saved-only four-cell collection.
No deformation-work/preference conclusion or goal completion yet.

## Earlier continuation launch

19 October: [scout report](lbt_whole_7ccn_reference_v1/SCOUT_REPORT.md).
Job14795 completed one native before-La cell, result5fe74c9b. Original native
admission passes; non-energy printed residual flags and small-gap warning remain.
Physical-state and independent force qualification remain false. Actual23,440s
on12CPUs (281,280CPU-s) is within the prior reference range. Root decision74c11c68
binds the actual scout, native receipt, original preflight and measured cost.

Unchanged after-La14880, before-Dy14881 and after-Dy14882 are RUNNING, each12MPI
and48,000MiB with allocation-derived MaxCore3100. Startup replay passed; native
launch reached and completion/alert-only health wakes armed. Exact submissions
are in lbt_whole_7ccn_reference_v1/continuation_v1/. Do not repeat14795 or submit
a waiting collector. Collect all four admitted results before deformation-work
comparison; no current La/Dy preference conclusion. Candidate14806 stays closed.
Production and failed branches remain preserved; the broader goal remains unmet.

## Whole-Mex directional check complete

19 October: [actual report](mex_whole_force_v1/REPORT.md). Preflight14810 and
H200 run14811 complete; all32 state/direction/step comparisons pass for MACE,
fixed-charge ALPB and total under unchanged gates. Both aggregate inputs remain.
Actual64MACE/16DyALPB calls; four MACE centers, two ALPB centers and16La
displacements reused, no central reruns. Resultc8c61664/manifest27c91c6c,
sourcebf32bc23/config94718e39. Cost980CPU-s/81GPU-s including preflight.

This qualifies only the sampled BEFORE derivatives. No optimizer, physical-spin,
accuracy or transfer qualification follows. Preserve all native outputs and exact
requests in mex_whole_force_run_v1; delayed callbacks require acknowledgement,
not repeats. Independent7CCN DFT14795 remains RUNNING at10:39UTC,12CPUs;
its actual-result/cost gate remains. Candidate14806 stays closed. Goal unmet.

## Natural-source state contract prepared

19 October: [source-only preparation](natural_transfer_contract_v1/REPORT.md)
binds Mex8FNS and both Hans sources, their physical/topology/charge/hydration
artifacts and deposited assembly records. Inventory `e7d7ce48`; root metadata
check `5946410f` passes 298 assertions and 107 unique pins. Zero molecular calls
or scientific allocations. Configuration A/B stays distinct from biological
chains; local electronic states stay distinct from unassigned whole-system
states. Both Hans extracted-subunit comparisons remain separate.

Metadata fidelity passes; **transfer scoring readiness remains false**. No
Hamiltonian, aggregate-spin or partner-state default is supplied by this contract.
The mechanics/response gates remain unmet. Failed electronic-base candidate
stays closed; independent7CCN reference14795 retains its actual scout-admission
and measured-cost gate before the remaining three cells. Read current receipts.

## Electronic-base candidate closed; independent reference continues

19 October: GPU14804, saved assembly14805 and comparison14806 complete.
Read [the decision report](whole_base_delta_comparison_v1/FINAL_REPORT.md).
Candidate deformation costs La9.870817/Dy9.200125 versus wholePBE0
12.279593/11.722656 kcal/mol. Both frozen1kcal gates fail. WholeMACE remains
closer on this consumed fixture. Result `ae68645d`, independent review `35cce123`.

All six declared development calls complete; four conditional7CCN candidate
calls are **closed and unrun**. No optimizer, retuned gate or production change.
Raw pending-slot labels remain as produced; ROOT_DECISION explicitly closes them.
Preserved14803 failed before any molecular call; reviewed saved-authority recovery
14804 passed. Total260CPU-s/14GPU-s including failure. Exact raw results, forces,
parameters, source and costs are archived. Independent7CCN wholePBE0 scout14795
continues unchanged; review actual admission/cost before its remaining3 cells.
The overall mechanics goal remains active and unmet.

### Historical launch and recovery record

### Passed scout

Actual preflight **14799** passes all eight states, seven regional reuse pairs
and all eight ALPB receipts, zero molecular calls, 6 CPU-seconds. Native scout
**14800** also passes unchanged admission: one whole256 GFN call, 65 SCC
iterations, 14.72 s native wall time versus the declared 600 s limit; 20 allocated
CPU-seconds. Exact raw result and root review are preserved under
`whole_base_delta_v1/actual_scout_v1/`. No candidate work/accuracy result yet.

The separate continuation adapter is being reviewed before the remaining three
7CCO whole-GFN calls, one missing regional-GFN and one regional-MACE call.
Reuse the scout exactly; first failure stops new scheduling. Both-metal work
triage still controls the four prebound 7CCN calls. Independent DFT14795 continues.

### Prior submission record

One-CPU saved-data preflight **14799** is submitted, completion watcher armed.
Source `f760cd79` / config `8dc6cdf0` passed independent static review. This stage
replays all eight declarations and saved admissions with zero molecular calls.
Read actual workspace `whole_base_delta_preflight_v1` before the one-call scout.
No candidate molecular calculation has been launched at this checkpoint.

19 October 2026 UTC: live scheduler/native output confirms reference scout14795
running with 12 MPI ranks and SCF initialization reached. Its completion gate
remains unchanged. The preceding goal turn made scientific progress; this turn
also verifies the live wait and addresses a separate mechanics limitation.

[Exact native-source review](whole_model_physics_review_v1/SOURCE_PHYSICS_REVIEW.md)
shows finite-range MACE atom-energy interactions and broadcast total charge/spin
labels, without an explicit long-range vacuum electrostatic or charge-equilibration
solve. Fixed-charge ALPB does not feed back into MACE. These facts identify a
model limitation, not the cause or magnitude of the observed error.

Root declared a [bounded subtractive test](whole_base_delta_design_v1/SPEC.md):
whole GFN2 + regional MACE − the same regional GFN2 + existing ALPB. Include
whole-GFN2 + ALPB alone as a no-extra-call comparator. Preserve both fixed LBT
reference pairs, the baseline172/two-cap/charge−2 region, full cap derivatives
and exact qualified La/Dy parameters. GFN-level self-consistency does not make
the regional correction mutually polarizable.

Ceiling: eight whole-GFN2 calls and one missing regional MACE/GFN2 pair,
before exact reuse. Seven regional pairs and all eight ALPB cells already exist;
reuse the certified 7CCN after-La SCC component, never repeat its failed default.
Zero-call Slurm preflight and one CPU scout come first. Complete the four consumed
7CCO cells, then allow the four prebound 7CCN challenge cells only if both metals'
candidate work errors are within the frozen 1 kcal/mol development triage rule.
No new DFT, optimization or transfer scoring follows from preparation.

The exact specification and reuse inventory passed independent static review;
the declaration is preserved in pushed commit `48bb2f6`. Actual saved-data
preflight and scout results are still required. Linnaeus owns the native adapter,
Harvey the separate saved comparator, Arendt independent review, root submissions.
Read actual newer receipts before launching. The
[original-plan audit](plan_evidence_20261019_v1/REQUIREMENT_MATRIX.md) retains
0/6 converged broad searches, failed Mex derivative qualification, the challenged
hybrid scaffold and unfulfilled natural transfer. No smaller goal substitutes.

## Whole-system Mex panel complete; independent LBT reference running

19 October 2026 UTC: Mex collector 14794 completes all eight fixed whole-system
MACE + unchanged ALPB cells. One H200 scout and seven remaining forwards, no
optimization or repeated completed chemistry. Read the
[result](mex_whole_system_v2/REPORT.md) and independent collection review.
Result `4c0c8b73…`, source `dd611ca0…`, manifest `acdaa462…`.

For the same job12970 La-proposed motion, total work is La −1937.990241 /
Dy −1942.176256 kcal/mol under aggregate multiplicities 10/15, and
La −1944.455151 / Dy −1948.078345 under 4/9. Both favor Dy more along this
motion, reversing the old hybrid ordering. This is model sensitivity, not a
validated preference. Spin changes individual work by −6.464910 / −5.902089,
above the frozen 1 kcal/mol resolution; differential change +0.562821 is below
it. Do not pick one aggregate input as a unique physical prediction. No Mex
optimizer, new spin sweep or reopening of the failed derivative branch follows.

Independent 7CCN reference preflight 14793 passes all four fixed states, saved
ALPB and available force replay with zero molecular calls. Before-La scout
**14795 is RUNNING**, 12 MPI ranks / 48,000 MiB allocation / MaxCore 3,100 MB
per rank, existing native PBE0 method. Completion and alert-only health wakes
are armed. Read `lbt_whole_7ccn_reference_v1/` and workspace
`lbt_whole_7ccn_before_La_v1/`. Its START.json uses an initial failure sentinel;
native output and scheduler confirm active execution, with no terminal RESULT
at this checkpoint. Do not mistake that initial field for a completed failure.

Root reviews the actual native scout result and allocation cost before the
remaining three declared DFT cells. No duplicates or waiting collector allocation.
Whole-peptide accuracy and natural mixed-ion Mex state remain separate questions.
The goal remains active and unmet; production and historical failures are intact.

## Earlier declaration: coherent whole-system intervention

19 October 2026 UTC: root declared [the finite replacement test](whole_system_intervention_v1/SPEC.md),
independently reviewed. Candidate: whole-system native MACE plus the same fixed-FF
ALPB, with no added MM or GFN embedding term. LBT already uses this model; reuse
its completed qualification/searches. The inspected 7CCO DFT fixture is consumed
development evidence, not prospective validation.

Mex will use the exact job12970 La start/terminal, all2,049 atoms and four reused
ALPB receipts. Eight fixed gradient cells cover both metals and two declared
aggregate-spin hypotheses (La10/Dy15 and La4/Dy9). First before-La10 H200 scout;
remaining seven depend on finite/state and measured600s/90%memory gates. Neither
hypothesis identifies a localized spin state. No optimizer or preference claim.

A separate whole256 PBE0 challenge freezes job12968's 7CCN-La start/terminal,
both metals, reusing four existing MACE energies/ALPB receipts. First before-La
scout; remaining three depend on actual admission/cost review. Old native method,
effective/physical states and reference limitations remain. Missing Dy MACE
forces stay missing. Source preparation only at this checkpoint; no new jobs yet.
Root owns execution; Harvey owns Mex adapter, Linnaeus reference adapter, Arendt
independent review. Closed Mex derivative and challenged hybrid remain closed.

## Matched LBT scaffold transition complete

19 October 2026 UTC: audit14787 completed all four matched states and passed
all original closure gates:36contexts/816queries,20seconds on1CPU. Including
preserved zero-query failures14785/14786, cost25allocatedCPU-seconds. No new
electronic evaluations, geometries or optimization. Source d736a589/configbf88807d,
result1b320766 and independent review5d0e62aa are preserved under
[lbt_scaffold_transition_v3](lbt_scaffold_transition_v3/REPORT.md).

Classical work is−23.281560La/−23.280382Dy at the same shared coordinates.
Retained bonded work is−10.682960; its terminal-touching subset is−9.532816.
All three nonbonded groups also lower energy. The terminal force discrepancy
has a persistent offset and changes along the motion: La hybrid error-change
RMS4.463637 versus whole-MACE0.530091kcal/mol/Å. No single-term repair follows.
Full Dy hybrid work/transition and before-Dy whole-MACE forces remain missing.

Root reviewed and archived the complete result, raw query provenance and actual
costs; Arendt's independent result review passes. All jobs are terminal. The
next mechanics design should test consistent scaffold/coupling treatment on
this independent saved deformation, preserving the denominator and reference
limits. No fitted one-term correction or restart of the challenged optimizer.
The broader goal remains active and unmet; Mex remains closed.

## Mex consistent-grid test complete: derivative gate still fails

18 October 2026 UTC: all three native endpoints are admitted. Saved collector
14107 completes the comparison and returns `failed_derivative_gate`; its exit 1
records the scientific gate failure. **Close this finite candidate.**

Embedded PBE0 gives analytic slope −3.673646 versus energy-difference slope
−3.797017 kcal/mol/Å: residual −0.123372 exceeds the original 0.005380 tolerance.
MM and fixed-charge ALPB pass their individual checks; the total fails.
The new matched-grid recipe shows a smaller discrepancy than the historical
0.414201 result but still fails consistency. No unique cause or physical accuracy follows. Read the
[final report](mex_dft/consistent_grid_v1/final_collection_v1/REPORT.md).

The exact result `4aa50af0…`, fresh central assembly, all three native pins,
unchanged source `1600a0a5…`, config `ecee028c…`, costs and independent review are
preserved under `mex_dft/consistent_grid_v1/final_collection_v1/`.
Collector 14107 made zero molecular calls. All candidate jobs are terminal.
Molecular allocations total 1,363,404 CPU-seconds; listed support adds 33.

Retain every SCF residual and stopping-route caveat, including the minus COSX
scalar-closure flag. These are not error bounds or post hoc rejection criteria.
No retries, changed thresholds/steps, grid/solver sweep or production change.
The original remaining six qualification and eight response cells stay stopped.
The completed LBT reference still challenges the hybrid scaffold force balance;
that is the next mechanics target, not an excuse to change the input structures.
The broader goal remains unmet. Earlier execution checkpoints are historical.

## Whole LBT four-cell reference complete

All four native endpoints and the full saved collector13424 are complete/admitted.
The same fixed deformation costs PBE0+ALPB La+12.279593/Dy+11.722656 kcal/mol,
compared with whole-MACE+ALPB +11.290076/+11.082401 and hybrid La−28.244859.
This challenges the hybrid benefit on the tested fixture. After-state whole-DFT
and whole-MACE terminal forces strongly oppose the motion; initial projections
are positive, so the result does not rule out every small displacement.
Read [the full report](lbt_whole_reference_continuation/FINAL_REPORT.md).

All four retain electronic-state and force-qualification limitations. The
Dy−La work contrast −0.556937 is below the declared reporting resolution;
no preference, universal MACE validation or unique MM parameter defect follows.
The complete result9229f7be, full256 force reports, native pins, reviews and costs
are archived in `lbt_whole_reference_continuation/full_result_v1/`.
Four molecular allocations total1154004CPU-s; listed support47CPU-s.
After-La13290 actual completion delivery is separately acknowledged.

Do not continue or promote the challenged LBT hybrid response.

## Saved partition comparison complete

Job13463 completed the four admitted archived-gradient comparisons on1CPU/1s,
with zero molecular calls. Expanding the electronic region leaves the Phe1CA
force projection excessive: La51.063→47.742 versus whole-PBE0+ALPB3.341;
Dy51.065→47.765 versus3.346 kcal/mol/Å. The MM contribution remains dominant.
Phe1CA stays outside the electronic-real region in both partitions. Global RMS
difference falls, but maximum atom error increases slightly at the Ala17 terminal
carboxylate carbon. No expanded after-state or finite-path validation exists.
Read [the partition report](lbt_partition_force_reference_v1/REPORT.md).

Result e25076f5, four full256 tables and independent actual review are preserved
under `lbt_partition_force_reference_v1/archived_comparison_result_v1/`.
The separate raw reconstruction failures13446/13455/13460 remain failed:
first an incorrect arithmetic association in the analysis code, then cap-scatter
last-digit differences up to3.55e-15 kcal/mol/Å. No tolerance was changed.
The successful comparison explicitly uses the previously admitted POINT forces,
requiring exact internal component closure and original baseline replay.
All four support jobs total3allocatedCPU-s; no molecular reruns.

Moving the tested boundary alone does not repair the scaffold force balance.
No unique parameter fault or La/Dy preference follows. The separate Mex three-cell consistency test now fails as described above. The mechanics goal remains unmet. Earlier checkpoints below are historical.

## Before-Dy reference complete; two LBT continuations remain

Job13291 terminated normally at18:43UTC and passes unchanged singleton admission:
256-atom gradient,17SCF cycles, charge−2/NEL1030, DyECP55, effective singlet and
physical multiplicity6,12MPI/MaxCore3100. Actual cost23058seconds ×12CPUs =
276696allocatedCPU-seconds. The result is pinned as `fc6f6ad2…` in
`lbt_whole_reference_continuation/activation_v1/before_Dy_result_v1/`.
Root and independent review agree. Only the orbital-gradient residual passes
its printed threshold; the other five flags remain false. Native admission
does not establish physical-state stability or numerical force accuracy.

After-La13290 and after-Dy13292 are still calculating gradients; Mex13281 remains
in SCF. Preserve the completed before-La scout and this before-Dy result, then
collect all four cells only after both remaining jobs terminate. No partial
deformation-work comparison, new molecular submission or retry follows.
The completion notification has now arrived. Its separate `ROOT_ACK.json`
records actual delivery and links the earlier proactive result review; the watcher
receipt and completed scientific result remain unchanged.

## Execution progress at 17:07 UTC

All three whole-LBT continuations print SCF convergence: after-La13290 after18
cycles, before-Dy13291 after17, and after-Dy13292 after19. They are recomputing
final exchange energies; final gradients, normal termination and native admission
remain pending. Mex13281 remains in SOSCF, last printed iteration20. All four
jobs are RUNNING. Continue unchanged, retaining the full four-cell LBT collection
and Mex central-admission/cost gates. No new scientific calls or submissions.
Exact native tails and scheduler observations are preserved in
`lbt_whole_reference_continuation/activation_v1/scf_milestone_v1/`.
This supersedes the earlier startup-only progress below, without changing any
scientific result, input, failure or execution setting.

## Latest: whole-peptide DFT challenges the extreme hybrid terminal force

Whole256 before-La13051 and original saved collector13283 are complete/admitted.
Along the frozen initial hybrid Phe1CA force direction, PBE0+sameALPB gives
3.341027 versus hybrid51.063006 and whole-MACE0.798547kcal/mol/A. Global force
difference RMS is10.141916 for hybrid versus3.511987 for whole-MACE. The large
hybrid terminal response is poorly supported by this conditional reference;
MACE remains quantitatively different and is not validated by this comparison.
Read [the actual scout report](lbt_whole_reference/scout_result_v1/REPORT.md).

Native18-cycle termination/full256 gradient and exact collector replay pass.
Four printed residual flags and all electronic-state/force-accuracy limitations
remain explicit. Cost284892CPU-s plus3CPU-s for collection. Root approved the
original remaining three fixed cells after this measured-cost/force review,
with no repeated scout, retries or optimization. Exact decision and config are
in workspaces/lanm_coupled_mechanics_20261017/lbt_whole_reference_continuation_v1/.
Post-scout source/config review, actual Slurm preflight13288 and all seven
negative controls13289 pass (5CPU-s total). The original remaining cells now
run as13290afterLa/13291beforeDy/13292afterDy,12MPI each, with completion and
actual singleton-manifest health wakes. Native output creation is verified;
all three printed molecular state/ECP/MaxIter settings and actual12MPI/MaxCore3100
match the frozen inputs. Independent startup review passes; terminal admission
remains pending. No waiting collector is allocated.
Read [the execution checkpoint](lbt_whole_reference_continuation/activation_v1/REPORT.md).
Before-Dy13291's first-build silence alert was reviewed against advancing CPU
accounting and the completed scout's1949.7-second first iteration. Continue
unchanged; convergence remains unestablished. The checkpoint preserves the
actual alert, scheduler observations and independent review.
After-Dy13292's matching first-build alert was separately reviewed with its own
output and increasing CPU accounting; continue unchanged, with no recovery.
After-La13290 has the same reviewed first-build status; all three continuation
jobs remain active without a justified cancellation or retry.
All three initial silence alerts have now resolved: first SCF iterations are
printed (after-La1967.8s, before-Dy1948.6s, after-Dy1969.1s). Exact rows and live
accounting are in the execution checkpoint's `first_iteration_review_v2/`.
Mex central13281 has reached iteration6. All four remain active without terminal
admission; no new mechanics or preference result is established by this progress.

## Latest: Mex mismatch localized to final integration; original gate still fails

Saved diagnosis13277 completed with exact original replay and zero molecular
calls. Printed final COSX contributes -0.401606 of the -0.414201 kcal/mol/A
derivative mismatch. Native point-charge gradient mapping reproduces exactly.
Subtracting the printed correction leaves -0.012595; subtracting the full
post-SCF shift leaves -0.011124. Both exceed the old0.005409 tolerance.
The post-SCF shift includes an explicit remainder beyond printed rounding;
neither diagnostic energy definition is promoted. Read
[the saved diagnosis](mex_dft/gridx3_diagnosis_v1/REPORT.md).

Cost2CPU-s, no new molecular calls. The original remaining6 qualification and8
response cells remain stopped. A separate consistent-grid qualification is
specified under mex_dft/consistent_grid_v1/EXECUTION_SPEC_V1.md: fresh central
gradient first, then at most two displaced energies after actual admission/cost
review. Only SCF grid2 is brought to the original final grid3 construction;
original derivative tolerance and state/SCF settings remain. Root accepted this bounded design. Actual saved-only preflight13280 passes all
source, shared-component and three legacy-parser controls with zero molecular
calls (10CPU-s). Failed13279 remains preserved (10CPU-s); versioned adapter1600a0a5
corrects only a native text delimiter. Root and independent review admit the
prepared inputs. Fresh central13281 is submitted on12MPI, config28abed13,
completion watcher2759506 and alert-only health armer2759507. Actual native startup
passes full grid/state/runtime checks: grid2 equals original/final grid3; XC/grid1/
grid3 stay unchanged. Charge0/multiplicity1/NEL832, MaxIter125 and12MPI/MaxCore3100
match. Terminal gradient admission and measured cost still gate both displacements.
No plus/minus or original-six/response-eight cells are submitted. Exact receipts:
mex_consistent_grid_central_v1/SUBMISSION.json; see the consistent-grid report.
Whole13051 and its original collector13283 have completed as recorded above.

## Mex derivative scout failure

Both Mex endpoints12977/12978 terminate normally and pass native admission.
Original saved collector13271 fails PBE0/total directional consistency: embedded
PBE0 analytic -3.673754 versus finite difference -4.087955 kcal/mol/A;
residual -0.414201 exceeds 0.005409. MM and fixed-charge ALPB pass.
Read [the failed scout report](mex_dft/scout_pair_failure_v1/REPORT.md).
The original stop rule closes the remaining six qualification points and
blocks the dependent eight-cell response panel. No solver/step retry is launched.

Independent source review finds no mapping/sign/unit defect. The saved
COSX/gridx3 decomposition is complete under mex_dft/gridx3_diagnosis_v1/;
it makes no molecular calls and cannot override the original failed gate.
Whole-peptide reference13051 subsequently completed; its result and continuation
are recorded above. The goal remains unmet.

## What we learned

1. **Accessible motion matters.** Six broad-coordinate searches pass physical admission but none converges. The conditional Mex/LBT source discrepancy falls24.120991→6.209455kcal/mol. This compares unchanged models with more motion; scaffold Hamiltonians still differ. See broad_response/REPORT.md and BROAD_RESPONSE_CHECKPOINT.json.
2. **Matching scaffold treatment alone makes a modest difference.** The complete104-cell hybrid matrix retains the same four selected geometries. Source discrepancies are5.034374kcal/mol for the predeclared baseline region and5.125352 for the expanded region. Larger component changes compensate; small net region dependence does not validate each term. See lbt_hybrid_pool_recovery/result/REPORT.md and POOL_RECOVERY_CHECKPOINT.json.
3. **The new model changes terminal bonded mechanics.** Saved force audit13032 passes all four states. The largest force change is Phe1 CA, outside the electronic region, with hybrid maximum forces47–53 kcal/mol/Å at the old selected geometries. Saved MM groups make this a bonded-led change; no obvious duplicate interaction or charge-shift error was found. This does not establish which model is accurate. See hybrid_response_force_audit/result_v2/.
4. **Original response jobs are terminal.** Three finish80 accepted steps and remain nonconverged;7CCN/La fails embedded GFN after2 steps. The three admitted works are−28.245/−28.397/−31.301 kcal/mol. There is no complete paired response pool. See lbt_hybrid_response/original_four_state_result/.
5. **The finite numerical recovery failed at a later geometry.** Four-call13036 passes its three frozen agreement gates, and preflight13047 exactly reconstructs the original failed trial. Continuation13048 resumes correctly but fails embedded GFN2 again after10 accepted steps/11 evaluations. The new failed trial reaches500 SCC iterations with large residuals; isolated GFN2 converges. The frozen second-failure rule closes this branch. No four-state pool or partial-survivor substitute is available. See lbt_hybrid_response_recovery/failure_result/REPORT.md.

6. **The terminal force is specifically bond-and-angle dominated.** MM-only audit13122 passes all four original starts and all208 local terms, with exact native parameter/term checks and fullforce closure. The largest term in every case is N–CA stretching, followed by N–CA–C bending. Along the total MM force, bonds contribute24.96–30.62 and angles16.70–22.52kcal/mol/A; torsions/nonbonded are smaller. This identifies the modeled cost, not the accurate model. No copying/counting defect found. Audit+saved report13130 cost10allocatedCPU-s, zero other molecular calls. Read lbt_mm_term_audit/result/REPORT.md and ALL_TERMS.tsv.

7. **The two models disagree on the cost of exactly the same motion.** Whole-model completion 13140–13142 passes all native and saved replay gates. The 13017 La deformation lowers hybrid energy by 28.244859 kcal/mol but raises whole-MACE plus fixed-FF ALPB by 11.290076. Shared solvent cannot explain this disagreement; differing electronic regions and MM mechanics prevent unique attribution to stiffness. Whole-model Dy-minus-La deformation work is -0.207676 kcal/mol, below the existing 1 kcal reporting resolution. Read [the report](lbt_whole_model_completion/REPORT.md). Four energies and three force sets are available; hybrid after-Dy remains missing. No accuracy or preference claim follows.

LBT is a peptide control, not a second natural LanM.7CCN is Lu-conditioned. Mex has an EF3 target, three Nd spectators and167waters. No binding free energy, robust preference classifier or experimental transfer is established.

## Owned work

| Jobs / task | Work | Resources / gate |
| --- | --- | --- |
|12977/12978 complete;13271 failed derivative gate|Embedded PBE0 external-charge La± derivative scout|Both endpoints available; PBE0/total derivatives fail, MM/ALPB pass. Remaining six stopped.|
|13051/13283 complete|Whole256/q−2 PBE0 before-La scout and saved collector|Native admission and exact saved replay pass; one new reference cell.284895allocatedCPU-s total.|
|13281|Mex central La gradient with consistent SCF/final COSX grids|Actual zero-call preflight13280 and native startup checks pass.12MPI; terminal gradient admission and measured cost gate the two displaced calculations. Completion and alert-only health wakes armed.|
|13290/13291/13292 running|Three remaining whole-LBT reference cells|Original source/settings unchanged. Actual13051/13283 admitted; preflight13288 and seven controls13289 pass.12MPI per singleton, completion/health wakes armed.|
|Stopped at prerequisite|Mex eight-cell DFT response adapter|Prepared source/tests remain preserved. Original scout13271 fails derivative qualification; no response cells submitted.|
|13122/13130 complete|Four-start MM term attribution and saved report|4/4 states,208 terms, all closures pass. N–CA stretching and N–CA–C bending dominate.10allocatedCPU-s total. See lbt_mm_term_audit/result/REPORT.md.|
|13140/13141/13142 complete|Whole-LBT fixed-deformation model completion|Native admission and saved replay pass. Two MACE plus one ALPB calls; four energies/three force sets. 140 CPU-s/11 GPU-s. DFT gates unchanged.|

Read ACTIVE_CHECKPOINT.json for exact active submissions and watcher records. All four original response jobs13017–13020 and recovery13048 are terminal. Recovery costs151s on12CPUs/oneH200; partial values remain diagnostic only. Solver qualification13033/13036 and saved force audit13032 are terminal. No paired pool or figure is submitted. Root owns submissions; agents own adapters/reviews. All numerical science uses Slurm.

Historical single-endpoint checkpoint, before plus completion: Mex minus 12978 completed on 17 October at 09:38 UTC: 32 SCF cycles, normal native termination, adapter status complete/DFT available. Its final gridx3 exchange step took 1,920.673 s; the long silence resolved without intervention. Scheduler cost: 22,132 s × 12 CPUs = 265,584 allocated CPU-s. Native energy-based convergence does not satisfy every printed density/gradient criterion; all flags and physical-state limitations remain. Read [the single-endpoint checkpoint](mex_dft/scout_minus_result_v1/REPORT.md). No paired derivative or mechanics accuracy result is established. Plus 12977 and whole-peptide 13051 remain active; no dependent reference is launched.

**Next:** collect the whole-LBT four-cell comparison after actual13290–13292 terminal results. Preserve all missing/failed slots; no repeated scout or retries. For Mex13281, review the actual terminal central gradient and measured cost before either displaced calculation. Preflight13280 and central submission are already complete; do not repeat them. No displaced Mex or original response cells are submitted.

The saved [startup timing review](mex_dft/consistent_grid_v1/startup_review/PHASE_TIMING_REVIEW.json)
records that the prior Mex first SCF iterations took about20minutes after startup
and guess. Output silence during the new first build alone does not establish a
stall. Both jobs were RUNNING with advancing CPU accounting at12:00UTC; whole13051
was calculating the COSX gradient. These are activity observations, not convergence
or force-accuracy evidence. Existing health/completion wakes remain active.

The saved audit's first CPU-host attempt13031 failed an upstream exact derivative replay before producing rows. Unchanged source/config/gates pass on the original qualified host/thread setup as13032. Cross-host exact-arithmetic portability remains unresolved; use the qualified replay environment, and preserve the first failure. See hybrid_response_force_audit/failure_v1/ and result_v2/.

## Numerical qualification and failures

Hybrid total/component derivatives pass the declared boundary and field directions at both steps and both regions (12986/12989). Numerical force consistency is distinct from electronic accuracy. The original12975 exact-map failure is preserved;12984 demonstrated only derived-Jacobian roundoff, with all coordinates/state exact. Versioned686eba84 permits bounded derived-coefficient roundoff only.

Original pool12992 preserves66 successful cells, one genuine isolated-GFN SCC failure and37 untouched cells. Four-call13009 qualifies one isolated-component numerical recovery without changing geometry/state/temperature/parameters/iterations/gates. Saved13012 validates its explicit certificate;13013 completes only the37 untouched cells with the original solver. All104 cells now pass. The original failure is immutable; this is not a global solver fallback or a unique-ground-state claim.

The primary baseline electronic region (complete residues3–14) was frozen before full-pool outcomes. Expanded2–15 remains sensitivity only. Existing actual raw data, states, caps, MM ledger, native parameter lineage and source-local geometry pools are retained.

## Limits and restart

Experimental LBT terminal chemistry has already had a complete bounded primary-source review: [paper, full supplement and deposition audit](../lanm_coupled_response_20261016/LBT_TERMINAL_CHEMISTRY.md). Those records do not establish the synthesis terminal groups and retain a charge-declaration discrepancy. The free NH3+/COO− model remains conditional. Repeating those same sources cannot resolve the missing experimental specification; this uncertainty does not change the ongoing fixed-state force comparison or authorize a new terminal-state campaign.

The prior broad searches and current scout are not minima. No automatic optimizer-budget extension, new source selection, clamped-core restart, parameter fitting or production promotion follows from success. A new failed required cell leaves its dependent result unavailable.

Mex outside-QM motion is distributed: the largest3-residue window contributes only7.42% of peptide squared displacement. One local expanded-DFT region cannot validate global MM lowering; that branch is not launched. See mex_motion/FEASIBILITY.md.

Standard scheduling currently permits12CPU/job and120CPU/user; GPU partitions have separate limits. Use actual allocations and current policy, no memory/time/priority overrides. Existing unrelated jobs remain untouched. Read actual receipts before acting; delayed events never request repeated chemistry. Root owns status/submission. Agents own only disjoint source/review tasks.

Saved preparation13059 failed an upstream exact terminal physical/force replay before molecular calls. Its Python binary matched, but NumPy dependency metadata differed (1.26.4 versus original2.4.6). Recovery13060 exposed a missing Gemmi dependency before replay.13063 keeps the original environment and supplies only the exact existing Gemmi package/native libraries in an owned overlay. Neither original environment is modified. The specific terminal numerical field remains unresolved until actual evidence. See mex_dft_response/preparation_failure_v1/.

Mex response adapter delivery: mex_dft_response/READINESS_REPORT.md. All six preparation/diagnosis/preflight allocations are terminal and preserved; no molecular calls were made. The original full physical/source/force checks and native acceptance rules remain unchanged. Mex12977/12978 are complete and scout collection13271 failed its derivative gate. The separate consistent-grid central scout13281 and whole-LBT remaining13290–13292 are active; the original response panel remains stopped.

## Fixed-deformation model completion — 17 October2026UTC

A bounded cache audit found missing whole-MACE energies at the exact 13017 La terminal for both metals, plus its Dy fixed-FF ALPB. Root froze lbt_whole_model_completion/SPEC.md before execution; independent source review passed. Preflight 13140 passed with zero molecular calls. Run 13141 completed exactly those three evaluations on 12 CPUs/one H200. Saved collector 13142 replayed actual native outputs, input states, full force arrays and the comparison. Both before points and after-La solvent were reused. Energy coverage is 4/4; whole-MACE force coverage is 3/4 because archived before-Dy was energy-only. The missing hybrid after-Dy point stays unavailable.

The same La displacement has hybrid work -28.244859 versus whole-model work +11.290076 kcal/mol. This establishes Hamiltonian dependence, not which model is accurate or a unique MM defect. Total cost is 140 allocated CPU-seconds and 11 GPU-seconds. Full results, source and actual commands are preserved under lbt_whole_model_completion/; exact native pins are in ACTIVE_CHECKPOINT.json. Independent integration review agrees.

The ongoing before-La DFT scout tests the initial forces and cost. It alone cannot validate the finite deformation work; the existing admission/cost gates still control further reference cells. No closed-response recovery, new geometry or DFT continuation was launched by this adjunct. All three completion watchers are recorded; delayed callbacks require acknowledgement, never repeated chemistry.

## Historical source-preparation checkpoint — superseded by activation above

The original runner permits only before-La. A separate [continuation draft](lbt_whole_reference_continuation/CONTINUATION_DRAFT.py) now implements the three already-declared remaining cells as singleton runs and a saved four-cell collector. Read [SPEC_DRAFT.md](lbt_whole_reference_continuation/SPEC_DRAFT.md) and STATIC_DRAFT_REVIEW.json. Author AST inspection and independent static source/provenance review pass; no scientific imports, numerical tests, runnable configs, continuation decision or submissions were made.

The draft replays original13049 unchanged, then separately admits the completed13142 after-Dy solvent. It forbids repeating before-La, retains unavailable comparator slots and uses the original method/state/admission rules. Actual13051 completion and its original saved collection, root force/cost review, a bound continuation decision and real-artifact Slurm checks are still required before execution. No physics result or execution readiness is implied by source review. Root owns those gates and all future launches; source preparation is complete.

The separate saved-artifact validation harness is also prepared and independently
reviewed: `VALIDATE_DRAFT.py`, interface note and `STATIC_VALIDATION_REVIEW.json`
in the same continuation directory. It requires a genuine positive replay before
seven rejection checks; supplemental geometry/metal checks explicitly exercise
the unchanged native admission helper. AST/source review only: no scientific
imports, executed tests, bound configuration, decision or submission. Run the
actual checks through Slurm only after the original scout/collector and cost
review. Native continuation execution and the later collector need their own
actual validation; this source review does not establish either.
