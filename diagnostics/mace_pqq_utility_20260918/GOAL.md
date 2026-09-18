# Agreed goal: useful PQQ discrimination with the working masked MACE model

Jacob: “we can't test on a set that has unproven labels. we need to test on
the set i have already prepared that allows us to discriminate. it needs to
be able to discriminate with the same fidelity, more or less, and faster”.
Then: “well. big dawg. set a new goal and proceed.” Reconfirmed: “continue”.

Determine whether the existing frozen, two-call masked-MACE PQQ discriminator
preserves the DFT baseline's classification fidelity on that prepared labelled
set and delivers a meaningful measured scoring speed advantage. Deliver the
paired results, runnable workflow, actual execution receipts and use/no-use
recommendation. Keep the production baseline and each model's own calibration
unchanged. New model experiments remain paused.

This supersedes the proposed PLM/extension analyses in DISCUSSION.md. Unlabelled
PLM predictions are not accuracy labels. No new external controls, scientific
models, fitted thresholds, geometry changes or non-PQQ experiments belong here.

Jacob cleared the earlier platform goal and requested this replacement. The
new goal is now active in the platform tracker, with no token budget. The older
goal was not falsely marked complete. This file records the execution scope.

## Execution plan, declared before fresh benchmark outputs

- Accuracy inventory: all 25 frozen canonical PQQ references, plus the three
  consumed crystal transfers, using each method's original published bands.
  Preserve calibration versus transfer roles, sequence groups and denominators.
  Reuse real scores; 1KB0 is unsupported by MACE and remains in the report.
- Timing: all 25 canonical cases, in the original fixed order, each once with
  fresh outputs for both methods: 50 masked-MACE bound endpoint evaluations
  and 50 native r2SCAN-3c/CPCM endpoint evaluations. No detached MACE evaluations
  are needed with the already qualified exact factorization reference.
- Use existing immutable prepared inputs and frozen MACE implementation. Time
  prepared-input-to-score operations, including validation, model loading,
  orchestration and reporting. Report historical preparation costs separately;
  this is not a measured raw-structure-to-score benchmark.
- Match CPU host hardware (node-128-512g-8gpu-1). MACE uses one A5000 and 16 CPUs,
  64474 MiB; DFT uses 32 CPUs, two concurrent endpoints with 16 MPI ranks each,
  no GPU. This follows existing endpoint policies instead of handicapping DFT
  to one core. Record actual allocation and overhead, including failures.
- Finite task inventory, no project CPU/time cutoff, no new scientific retry
  variants. Existing scheduler limits apply. Technical retries preserve inputs
  and remain visible. No historical cache hits count as fresh timing.
- Fidelity: require the existing 25/25 classifications and both supported
  transfers. Repeated scores must reproduce the originals within 0.01 kcal
  (model kcal for MACE). Report literal frozen-band classifications and any
  numerical boundary crossings separately; do not move thresholds or disguise
  an inconclusive result. This is reference fidelity, not independent accuracy.
- Practical speed target: at least 1.5x faster median prepared-to-score time
  and lower total time across all 25. Also report the full per-case distribution
  and resource differences. This engineering target is not a biological rule.
- Scope ends in an evidence-based recommendation. A negative result does not
  authorize returning to open-ended model development.
