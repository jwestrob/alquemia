# v5: scope DIIS parsing to the SCF iteration table

Verified bug in saved1220300health output: the unscoped eight-number regex treated paired orbital rows as DIIS iterations. Native energies2.0/0.0 and orbital indices90–112 entered residual fields. Successful convergence already suppressed pathology alerts for these four completed cases, but the recorded evidence was wrong and could mislead review.

The parser now opens only on the actual DIIS header with Energy(Eh),Delta-E,RMSDP,MaxDP,DIISErr,Damp,Time(sec). It closes on solver transitions, alternative iteration headers, convergence/failure, energy/orbital/gradient sections and termination. It retains the latest real DIIS table while marking it inactive afterSOSCF/TRAH. Subsequent orbital numbers never reopen it.

All four real1220300outputs now correctly endDIIS atiteration7, beforeSOSCF; their full-output parsedDIISrows exactly match the actual prefix beforeSOSCF. The realpost-SCFsection alone yields noDIISrows. Fourteenwatcher tests pass, including priorearlyDyalerts, near-rootstagnation andsuccessdebounce. See HEALTH_POLICY_V5_AND_TESTS.json for exact hashes/before-after evidence.

Future-only `native_scf_health_v5`: all scientific alert thresholds andv4debounce behavior remain unchanged. No existing watchers, frozen receipts, queues or jobs were modified. Use a freshv5receipt for future submissions. No molecular evaluation or scheduler submission.
