# One finite 7CCN-La completion attempt

## Question and selection

Can the unchanged whole-MACE plus fixed-FF ALPB model reach its original full-force criteria with one substantial finite continuation? Root SCIENCE_DECISION.json declares this attempt. Actual15456 supports the prior step80-to92 lowering within the existing1kcal/mol reporting resolution (reference−4.726291 versus model−4.226978). Full force disagreement persists. This supports testing further model relaxation, not force accuracy, reference convergence or metal preference.

Exactly one state and restart: actual15254 RESULTc445cf4b, final checkpoint317c2b24, terminal POINTeb96df29. Preserve original initial POINT50a597c8 and original prepared7CCN source; no terminal rebasing. The checkpoint contains92accepted/92distinct evaluations, eight L-BFGS history pairs and trust0.1Å. Its canonical JSON SHA is981ec084b57fa86599125619b9ecce7613c883597b23e38d3521be1eaa5c0da2.

## Unchanged mechanics and finite allowance

Additional limits:80 accepted steps and200 distinct new composite evaluations. Cumulative limits:172accepted/292evaluations. Original optimizer maximum atom-vector force<=0.10 and atom-vector RMS<=0.03kcal/mol/Å must both pass for model convergence. All256 atoms remain movable; original connectivity, chirality, peptide, contact and bound-state rules remain unchanged. Preserve wholecharge−2/La multiplicity1, MACE weights/float64/runtime, fixed-FF charges/ALPB, native gradients and exact cache identities.

CARTESIAN_OPTIMIZER_COMPLETION.py is copied from the reviewed scout optimizer3badac29. OPTIMIZER_DIFF.txt changes only entry messages, the exact old92/120→new172/292 limit admission, and exact permitted checkpoint digest. Every numerical kernel, setting, history/trust update, rigid-mode handling, Armijo rule, trial-rejection behavior and force criterion is unchanged. It accepts no other checkpoint, later retry, fresh search or generic budget extension. The original force/trust/search code remains immutable.

CONTINUE.py calls pinned original CONTINUE.prepare0c14a03e with exact actual previous CONFIG7127a7fc. This performs original saved source/gradient/version/derivative qualification and seeds the original cache without model initialization or molecular evaluation. Reuse that same evaluator instance and its empty output directory. Verify the actual15254 execution/preflight chain, then admit the union of original and previous saved records:113 records,109 unique coordinate/cache identities, and92 historical requests. Every admitted old cache key must remain represented; replace generated saved aliases with authoritative archived pins. New evaluator requests/costs start at zero. Historical92 requests and history remain in the exact resume checkpoint and are preserved by the unchanged optimizer checkpoint writer.

Re-admit original initial, step80 and step92 full gradients/components and coordinates using original helpers; exact step80-to92 work and original-source terminal guard must replay. Actual15456 comparison is pinned rationale; its agreement sign/resolution is not an execution accuracy gate. No force correction, fitted interaction or changed convergence tolerance.

## Stages and resources

`--stage preflight --config CONFIG_PREFLIGHT.json --output FRESH` performs zero-call saved admission under Slurm,1CPU/noGPU on qualified H200 host node-224-2t-8gpu-1. Original MACE venv and native asset checks remain unchanged. No trial geometry, evaluator.evaluate or model initialization. Future preflight_result/run_decision must be null.

`--stage run --config ROOT_BOUND_CONFIG --output FRESH` requires the exact actual passed preflight and separately pinned root run decision, same scientific binding and exact preflight evidence. Use12CPUs/oneH200 on the same qualified host, gpu_h200 partition, standard_expand QoS, scheduler default memory/time, no requeue. Original model loads lazily at the first new request and stays warm. At most200 new MACE and200 ALPB calls, one of each per attempted new composite; preserve failures/costs. No newDFT/GFN/MM, other metal/source, cross-score, intermediate reference or extra search.

Root run decision schema: protocol `lbt_7ccn_completion_attempt_v1`, status `approved_single_completion_attempt`, exact `scientific_binding`, `preflight_result`, `checkpoint`, state `["7CCN","La"]`, and cumulative_limits `{"accepted_steps":172,"evaluations":292}`. Root owns launchers/submission/watchers/runtime configs. No stage launches another.

## Stopping and output

Use original search terminations: both full-force criteria, either finite budget, native execution failure, accepted state-changing departure, zero rigid-free direction, or no admissible decreasing trial under the original finite backtracking rule. Physical trial rejections may backtrack exactly as before. Preserve every trial/rejection and admitted checkpoint. A completed execution with state change remains a recorded outcome excluded from bound-state conclusions. Budget exhaustion is incomplete relaxation; no automatic restart, expansion or further reference.

Record additional/cumulative counts, original-start and step92-relative native-component work, all raw point/gradient/request pins, full terminal force metrics, physical/state status, model load/native costs and final checkpoint. All unavailable outcomes remain explicit. No minimum beyond declared model convergence, force-accuracy, stable-electronic-state, affinity or preference claim. Preserve vacuum-dianion, conditional terminal chemistry, frozen-f and reference SCF/gradient limitations from the actual comparison. All scientific imports/replay/tests run only on Slurm; preparation is source/metadata only.
