# Single7CCN-La checkpoint continuation scout

## Question and frozen state

Does a short continuation of unfinished whole-MACE plus fixed-FF ALPB relaxation reduce the remaining forces while preserving physical admission? Reference15173 supports the large prior deformation but reports La/differential work errors above1kcal/mol. Saved audit15176 finds nominal local descent under both admitted gradients along the conditional next optimizer direction. These are rationale, not new accuracy gates or proof of finite-step reference descent.

Exactly one state: original7CCN-La run12968, terminal POINTcf59ef7a and final checkpoint4e3a6123. Preserve its original initial point50a597c8, original prepared source and source-relative PhysicalGuard, eight L-BFGS pairs, trust0.1Å, 80 accepted steps/80 distinct evaluations, and all99 archived point records (97 unique cache keys). No source rebasing, reset, repeated central calculation, other search or opposite-metal score.

All256 atoms remain movable under original topology/chirality/peptide/singular-contact and bound-state rules. Whole q−2/La multiplicity1, native float64 MACE weights/versions, fixed-FF charges/ALPB executable/build/libraries and exact source identities remain unchanged. No MM/GFN/DFT, new direction selection, fitting or solvent change.

## Narrow budget revision

Additional allowance: at most12 accepted steps and40 distinct new composite energy/gradient evaluations; cumulative caps92accepted/120evaluations from actual80/80. Each new composite attempts at most one MACE and one ALPB call. Exact cache hits do not spend distinct-evaluation budget. Retain failed attempts and all accepted checkpoints; native failure stops without retry.

CARTESIAN_OPTIMIZER_SCOUT.py is an isolated copy of originald93821da. OPTIMIZER_DIFF.txt contains only hashlib import and budget-admission changes: require a non-null exact original checkpoint canonical digest2be9f7cd, permit only original80/200limits to new92/120, and compare consumed counters against both budgets rather than requiring unused old200limit<=new120. All numerical direction, trust, line search, force convergence and physical-guard bodies remain unchanged. No later checkpoint or other state is admitted by this copy.

Convergence stays max atom-vector force<=0.10 and atom-vector RMS<=0.03kcal/mol/Å over the full system. Budget exhaustion is incomplete relaxation. A state-changing accepted endpoint remains an explicit outcome excluded from bound-state conclusions.

## Stages and resources

`preflight --config CONFIG_PREFLIGHT.json --output FRESH`: saved-only Slurm, one CPU/no GPU on qualified node-224-2t-8gpu-1, original MACE venv. Reuse unchanged original loader/admission/evaluator initialization; re-admit initial and terminal full gradients/components, source guards, original derivative qualification, every prior cache record and native package versions. No evaluator.evaluate, model initialization or molecular call occurs. Historical build/library checks also run, requiring the qualified host. Complete original4-state saved input loading is inherited; only one state may execute.

`run --config ROOT_BOUND_CONFIG --output FRESH`: same host,12CPUs/oneH200, standard scheduler defaults, no memory/time override or requeue. Original NativeMACE initializes lazily on the first genuinely new point and remains warm. The config must pin actual passed preflight and root's review before run. Its scientific binding must exactly match preflight; no automatic run follows preflight.

Root review JSON: protocol `lbt_7ccn_La_checkpoint_continuation_scout_v1`, status `approved_single_continuation_scout`, exact `scientific_binding` from preflight, `preflight_result`, `checkpoint`, state `["7CCN","La"]`, cumulative_limits `{"accepted_steps":92,"evaluations":120}`. Future pins remain null until actual evidence exists. No scheduler submissions belong to this adapter.

## Required outputs and limits

Record new requests/points/costs, additional and cumulative counters, original-start and checkpoint-relative component work, full terminal gradient pins/force metrics, state labels and immutable final checkpoint. Failures preserve partial outputs and one-state denominator. Root reviews actual outcome before any other continuation or reference. Missing Dy reference-comparator model forces remain missing; no pool/preference/source-robustness or improved-accuracy claim follows. Preserve native reference SCF, frozen-f, vacuum-dianion and terminal-chemistry caveats. All numerical imports/tests and execution occur on Slurm.
