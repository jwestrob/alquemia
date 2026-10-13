# Two-Ca recovery: inputs and ExaA reference audit complete

CPU preparation job 2879 completed on node-224-2t-8gpu-1. All three exact cached inputs passed the actual AF3 upstream schema: seed 1, chains A/B/C/D, every original field unchanged, with only a separate `[Ca+2]` ligand D appended from B. Root has submitted Q9Z4J7's five-model scout as **job 2881**, twelve CPUs/one H200 GPU. Inference has now completed and the separate [generated-source audit](audit/REPORT.md) passed all five unique-PQQ-Ca site checks. The remaining ten models also completed and passed the same source-geometry checks; see [C5/Q88 audits](audit/CONTINUATION_REPORT.md). Full protonation/context/scoring admission is separate. The preparation-state details below retain their original scope.

## Why two Ca is physically motivated for ExaA

The actual [1FLG](https://www.rcsb.org/structure/1FLG) coordinate audit maps Q9Z4J7 Asp45/Thr48/Asp51 to crystal Asp11/Thr14/Asp17 in both chains, with 100% reported identity over the aligned residues. Each crystal N-terminal Ca has exactly the six donors observed around the original predicted Ca: Asp O/OD1, Thr O/OG1, Asp O/OD1. This is the same residue/donor pattern, not a coordinate-superposition claim.

| Crystal Ca | Existing typed CN (≤3.1 Å) | Nearest PQQ atom (Å) | Site |
|---|---:|---:|---|
| A:751 | 7 | 2.731787 | PQQ-associated |
| A:761 | 6 | 25.282441 | N-terminal structural |
| B:752 | 6 | 2.794068 | PQQ-associated |
| B:762 | 6 | 25.431111 | N-terminal structural |

The [primary paper](https://doi.org/10.1006/jmbi.2000.3603) describes the additional N-terminal Ca as contributing to stability. The earlier single-Ca Q9 predictions therefore occupy an experimentally observed Ca site while leaving the PQQ site empty. The finite two-Ca ExaA input addresses that occupancy ambiguity. Two-Ca occupancy in PedE Q88JH5 and ExaF C5AXV8 remains a hypothesis; this audit does not establish it experimentally for those proteins.

## Frozen execution and scope

Active manifest: `workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/manifest.json` (SHA256 1d34a717b2d97fe4c59852d2f894eddab24c4e69f6bdb0d86839ada694f88fb7). Immutable implementation: `diagnostics/ca_coverage_recovery_20261012/refold/v3/`.

Prepared inputs and complete maps: `run_v3_h200/prepared/PREPARED.json`, `1FLG_AUDIT.json`, downloaded/hash-pinned `1FLG.cif`, and three actual schema logs. The manifest declares three proteins/fifteen samples; Q9's first five require an actual complete executable scout before either other protein can be separately dispatched. All generated samples and failures are retained; no ranking selection, extra seeds or automatic retries.

New AF3 execution explicitly uses 3.0.4, seed 1, five diffusion samples and ten recycles, with the original cached MSAs/templates and no data-pipeline search. The historical source-fold version/recycle count remains unverified, so old/new predictions do not isolate Ca count as the only possible execution difference. The SIF and weights remain pinned to the working ADH9 artifacts. GPU thread limits and affinity use the actual twelve-CPU allocation; no memory/time requests were introduced.

Prepared inputs are not scoring endpoints. Both ions remain in full sources. Unique PQQ-associated active Ca, spectator location and unchanged chemical/context gates must pass the separate `spectator_preparation/PLAN.md` before scoring. A spectator within 4.3 Å of any retained physical heavy atom in the full core/union is explicitly unsupported by that plan. Original 125 outcomes and corrected canonical bands remain frozen; there is no new fit or production promotion.

## Actual runtime evidence and preserved failures

- `run_v1`: metadata-only unexecuted draft, preserved and marked superseded.
- Job 2877 / `run_v2`: failed before molecular imports because the node-local Apptainer executable differs from login. No prepared input/model call.
- Job 2878 / `run_v3` on node 5: pinned shared runtime passed; exact Q9 input was written. Apptainer 1.1.3 failed namespace startup with `Could not write info to setgroups: Permission denied`. **The upstream schema never ran** despite the generic exception saying schema rejected. Actual log/receipt preserved; zero folds/energies or mapping.
- Job 2879 / `run_v3_h200`: unchanged v3 implementation completed using node-local Apptainer 1.5.4. All three schema validations and the 1FLG audit passed. This establishes a working approved H200 transport; it does not diagnose the precise node 5 kernel/runtime cause or bypass any configuration.

Job 2879 scheduler elapsed 15 s; receipt reports 11.876681 CPU-s, 7.303611 s inside the prepared-output phase, 199924 KiB parent / 185528 KiB child peak RSS. Scheduler MaxRSS separately reports 5106820K; these measurements have different accounting scopes and are retained as observed. No energy or fold calls occurred in preparation. Job 2881's completed inference cost/result and audit2884 are now recorded separately in audit/REPORT.md. Root owns its watcher 3530162, event `ca-coverage-ExaA-two-Ca-AF3-scout-v1`.

This agent performed code/metadata work and read actual receipts on login; all input/schema/structure operations executed inside the real Slurm job. Python AST checks passed; no independent scientific tests or submissions were run by this agent.
