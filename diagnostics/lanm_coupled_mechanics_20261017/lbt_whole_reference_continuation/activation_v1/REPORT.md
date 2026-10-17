# Whole-LBT remaining references: execution checkpoint

17 October 2026 UTC. The original before-La scout13051 and saved collector13283
are complete and admitted. Their conditional force comparison does not support
the extreme hybrid Phe1CA force. Root approved completion of the original four
fixed cells after reviewing the measured23741-second/12CPU scout cost.

The previously reviewed continuation source06e40d37 and specificatione23187e0
are bound unchanged. The historical draft names/protocol are retained for exact
reproducibility; actual authorization and evidence are in
`ROOT_CONTINUATION_DECISION.json` and `EXECUTION_REVIEW.json`.

## Actual checks and launches

| Job | Work | Result/resources |
| --- | --- | --- |
|13288|Zero-call continuation preflight|Passed,3seconds/1CPU|
|13289|Real-input replay and seven rejection controls|Passed,2seconds/1CPU|
|13290|after-La whole256 PBE0 EnGrad|Started,12MPI|
|13291|before-Dy whole256 PBE0 EnGrad|Started,12MPI|
|13292|after-Dy whole256 PBE0 EnGrad|Started,12MPI|

The positive replay precedes controls for missing/changed scout evidence,
changed retry permission, wrong cost receipt, repeating before-La, and incorrect
supplemental-solvent geometry/metal. All seven reject at their intended gates;
the forbidden scout creates no native directory. No external process or
molecular call occurred in validation. Original13049 replay precedes separate
admission of the completed after-Dy ALPB. Independent review verified465
referenced pins and exact source/config associations without a mismatch.

Each new native job contains exactly one remaining task. The runtime config
adds only the actual passed preflight pin. Charge/state, physical coordinates,
PBE0 recipe,12MPI and allocation-derived memory remain unchanged. No time,
memory or priority override is added. Actual source/config/command and
completion/health-watch records are in the three `*_SUBMISSION.json` files.
`STARTUP_RECORDS.json` records native output creation and armed health wakes.
The START file's initial failure-status placeholder is not a terminal result.
Subsequent native startup review in `native_startup/` confirms all three printed
charge−2/NEL1030/effective-singlet/ECP46-or55/MaxIter125 settings,12MPI and
MaxCore3100. Independent review agrees. Dy physical multiplicity6 and effective
molecular multiplicity1 remain distinct; atomic-guess multiplicities do not
override molecular state. No terminal native admission is implied.

## Before-Dy terminal result

Job13291 completed at18:43UTC with native normal termination and a full256-atom
gradient. Unchanged singleton admission passes, with17SCF cycles and the original
state/runtime settings. Actual scheduler cost is23058seconds on12CPUs,
276696allocatedCPU-seconds. Result `fc6f6ad2…`, native receipt, force report,
accounting and independent review are archived in `before_Dy_result_v1/`.
Only the orbital-gradient residual meets its printed threshold; energy, density,
rotation and DIIS flags remain false. Physical-state and force-consistency
qualification remain false. No full deformation-work comparison is available.

After-La13290 and after-Dy13292 still calculate gradients. Retain all four planned
cells and run the original saved collector after both remaining terminal results.
No molecular retry or new call. Earlier running/startup observations below are
historical; `ROOT_REVIEW` records proactive review, not delivered-event ACK.

## Interpretation and next action

Before-Dy13291 emitted a900-second output-silence alert during its first Fock
build. Actual scheduler observations show RUNNING with advancing aggregate CPU;
the completed before-La scout needed1949.7seconds for its first SCF iteration.
Root and independent review therefore retain the unchanged calculation and
watchers. This does not prove convergence, useful numerical progress or per-rank
utilization. Evidence and delivered-event acknowledgement are preserved in
`before_Dy_silence_review_v1/`. No cancellation, retry or new calculation follows.
After-Dy13292 subsequently emitted the same first-build alert. Its own raw
output, empty scheduler stderr, RUNNING state and advancing aggregate CPU were
reviewed separately; continue unchanged. The actual delivered-event snapshot
and decision are in `after_Dy_silence_review_v1/`. Neither Dy job has a terminal
result or a demonstrated SCF convergence trend at this review.
After-La13290's first-build alert was also reviewed against its own output,
empty scheduler stderr and increasing aggregate CPU. Continue unchanged; its
separate snapshot and acknowledgement are in `after_La_silence_review_v1/`.
All three alerts precede a printed first SCF iteration and establish no stall.

Subsequent native output resolves all three initial silence alerts: after-La,
before-Dy and after-Dy print their first iterations taking1967.8,1948.6 and1969.1
seconds respectively. Raw rows and live scheduler evidence are retained in
`first_iteration_review_v2/`. The v1 metadata field-count assertion is preserved;
no native calculation changed. All remain active; no convergence trend, terminal
admission or finite-work result follows from this startup progress.

At17:07UTC all three continuations print SCF convergence: after-La18 cycles,
before-Dy17 and after-Dy19. Each has reached final exchange-energy recomputation.
The `scf_milestone_v1/` snapshot preserves raw output tails, source manifests,
the unchanged runtime configuration and live scheduler accounting. Final energies,
gradients and normal termination remain pending; SCF convergence alone does not
admit a completed reference or establish electronic-state stability. Mex13281
remains in SOSCF at iteration20. Continue all four existing jobs unchanged.

These jobs test the energy and force changes along the same fixed La-proposed
deformation, for both metals. No optimization, geometry repair, original-scout
repeat or molecular retry is allowed. The complete four-cell collector retains
failed/unavailable cells and cannot report a survivor-only work contrast.
Missing before-Dy whole-MACE forces and after-Dy hybrid evidence remain explicit.

Wait for actual terminal native admission, then run the saved collector with
the original scout plus all three actual singleton results. No waiting collector
is allocated. PBE0 stability, numerical force accuracy, free termini, vacuum
dianion, frozen-f/SOC and solvent limitations remain; there is no new La/Dy
preference result or production change.
