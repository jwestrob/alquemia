# Whole-protein MACE fits on RTX A5000

Completed 2026-09-16, job **1200381**. All 12 approved calls completed:
four original cores and eight full-protein primary/repeat/translation/rotation
calls. The full 9,141-atom 1H4I protein runs with MACE-POLAR-1 medium, unchanged
float64 weights and vacuum realspace physics. No new DFT was run.

**Memory feasibility passes. Predictive improvement is unestablished.**
Repeat, translation and charge checks pass; rotation and the previously observed
hybrid partition check fail their frozen tolerances. The production baseline
and its reference/decision bands are unchanged.

## What changed

The original backend retained huge all-pair and local neural-network arrays.
The opt-in implementation sums electrostatic pairs in 512-site blocks, with an
analytic backward that recomputes blocks. It also checkpoints 4,096-edge local
interactions and 256-atom symmetric/sparse tensor products. An indices-only
accumulation derivative avoids retaining edge messages. All neighbors still
contribute before node normalization; global fields and charge constraints remain.

Four complete GPU core comparisons match the original implementation: maximum
absolute energy difference 7.276e-12 eV, force difference 5.941e-12 eV/Angstrom,
and density-coefficient difference 7.317e-14. All 17 real-artifact/parser/kernel
tests pass. Original weights and buffers remain bitwise identical. Original full
execution could not fit, so full-system equivalence is supported by these kernel
and core comparisons, not a successful unblocked full-protein reference.

Implementation commit: `801111a`. Kernel IDs:
`graph044_blocked_analytic_backward_v1` and
`mace0316_edge_product_checkpoint_v5`. The scientific protocol remains
`mace_polar_1m_vacuum_r2scan3c_subtractive_pilot_v1`.
See [implementation details](BLOCKED_KERNEL.md) and [approved scope](MEMORY_AGREEMENT.md).

## Measured runtime and memory

One RTX A5000, 16 allocated CPUs, 64,474 MiB requested host RAM; native GPU
execution, no host offload. The primary pair includes energies and analytic forces.

| Measurement | Actual value |
|---|---:|
| La / Ca evaluation | 58.137043 / 58.324431 seconds |
| Primary pair evaluation | 116.461473 seconds |
| Primary pair, including subprocess startup | 130.417239 seconds |
| Primary pair allocated CPU time | 2,086.675825 core-seconds |
| Peak GPU tensor allocation, all full calls | 10,250,326,016 bytes = 9.546360 GiB |
| Peak GPU allocator reservation | 11,379,146,752 bytes = 10.597656 GiB |
| Maximum worker host RSS, all full calls | 1,725,688 KiB = 1.645744 GiB |
| Successful 12-call allocation | 564 GPU-allocation seconds; 9,024 core-seconds |
| All memory-development jobs, including failures | 898 GPU-allocation seconds; 14,368 core-seconds |

The memory-development jobs reported 1,000.692 actual CPU seconds. These are
allocated GPU times, not GPU utilization integrals. The preceding interface/core
phase separately used 69 GPU-allocation seconds and 1,104 allocated core-seconds.
CPU test/preparation work was not fully profiled. These measurements establish
single-case feasibility, not whole-pipeline speed or a baseline-relative cost ratio.

All failed attempts remain under `blocked_v1` through `blocked_v5`.
Jobs 1200368/69, 1200370, 1200371/72, 1200379 and 1200380 exposed edge weights,
TorchScript checkpoint handling, retained messages and per-atom tensor products.
The offload attempt 1200372 reached 81.017 GiB host RSS and failed on GPU memory;
its requested RAM was not an observed hard process cap. It had already failed
when cancellation was requested. The final implementation needs no offload.

## Frozen physical checks

| Check | Observation | Gate |
|---|---|---|
| Four-core original implementation agreement | All energy/force/density checks pass | PASS |
| Full-system charge closure | Maximum absolute error 1.187e-12 e | PASS, 1e-5 e |
| Repeat | Zero reported energy/contrast change; force maximum 4.49e-11 eV/Angstrom | PASS |
| Translation | Zero reported energy/contrast change; force maximum 3.69e-11 eV/Angstrom | PASS |
| Rotation: endpoint energies | La -1.867694; Ca -1.811157 kcal/mol | FAIL, 0.01 kcal/mol |
| Rotation: La/Ca contrast | +0.056536496 kcal/mol | FAIL, 0.01 kcal/mol |
| Rotation: forces after mapping back | Maximum 0.059996344 eV/Angstrom | FAIL, 0.001 eV/Angstrom |
| Hybrid qm36-minus-qm33 partition shift | -4.998674456 kcal/mol | FAIL, 2 kcal/mol |

The unchanged backend uses Cartesian displaced charges to represent dipoles,
with finite offsets (0.1 Angstrom features; 0.02 Angstrom energy). This is a
candidate explanation for rotation sensitivity, not an established attribution.
No offset, precision, tolerance or model was adjusted after inspecting results.
Original-vs-blocked rotated-core and representation checks would be a separate,
proposed investigation; they have not run.

The hybrid expression is `E_full,MACE + E_core,DFT - E_core,MACE`, in vacuum.
The partition difference cancels the shared full-system MACE term. Its roughly
5 kcal/mol residual is smaller than the archived 61.739 kcal/mol DFT-core shift,
but this is not evidence of improved biological classification. Raw elemental
energy offsets are large; no compatible aquo reference, calibrated S or class
is available. No baseline threshold is inherited.

## Reuse and remaining work

The memory machinery is reusable for compatible larger checkpoints; **MACE-POLAR
large has not been loaded or evaluated here**. Its architecture and original-core
agreement must be checked before asserting compatibility or predicting its memory.
Training, Hessians, periodic systems and unsupported architectures remain excluded.
Quadratic pair work remains even though its saved memory is bounded.

Recommendation: retain the baseline for predictions and keep this memory backend
for the next MACE experiment. Numerical credibility is partial, predictive value
is untested, and single-protein medium-model inference is affordable on A5000.
No memory-work jobs or H200 continuation remain active.

## Authoritative artifacts and next command

Workspace prefix: `workspaces/mace_hybrid_20260916/`.

- `blocked_v6/manifest.json`, SHA256 `8d5fc630c11b50457cd30927f9bcb827f99bd1c67191b307ccb01bc0712c7189`.
- `blocked_v6/collection_job_1200381.json`, SHA256 `025b5ae120368354f0368e656646ccbc2386e33f398b724b6ddfd3dc3f5276bb`.
- `blocked_v6/REPORT.md`: all unrounded endpoint energies, components and checks.
- `blocked_v6/core_comparison.json`: original-core equivalence; per-task execution directories contain forces, densities and receipts.
- `memory_sacct_final.tsv` and `memory_cost_final.json`: terminal accounting including every failed memory-development job.
- `blocked_regression_tests_v6.log`: 17 passing tests; `blocked_v6/job_1200381_accounting.json`: accounting watcher receipt.

Read-only verification, from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  scripts/mace_hybrid.py dry-run \
  --manifest workspaces/mace_hybrid_20260916/blocked_v6/manifest.json
```
