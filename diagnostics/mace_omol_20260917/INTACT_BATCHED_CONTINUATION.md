# Complete the original intact-chain comparison with a verified memory adapter

Declared2026-09-17 before evaluating the remaining four intact proteins. This
retains the scientific states, five-case panel, score expression and three
relative-order criteria in INTACT_CHAIN_PLAN.md. The physical protocol remains
`mace_omol_intact_chain_matched_coordination_v1`; the execution adapter has its
own identifier, source hashes and cache identity. No threshold is added.

Exact interaction batching has passed20 actual core checks at two batch sizes
and31 intact ALPHA_1F6S checks at1024edges per batch. The intact memory peak is
4,957,639,680bytes and inference approximately5.3seconds. This is feasibility
and numerical self-consistency, not yet full native equivalence or prediction.
H200 native reference1200809 remains queued. NATIVE_CPU_PLAN.md supplies an
independent unmodified float64 native reference with the same physical inputs.

Only after the native reference's full qualification and the24 endpoint/paired
comparisons in mace_omol_edge_report.py pass0.01kcal/mol, execute the remaining
16 endpoints on oneA5000/16CPUs/64474MiB host with1024edges per batch. Reuse the
four actual batched ALPHA_1F6S primary endpoints, with explicit source receipts.
Do not substitute native results for a failed batched endpoint. Verify actual
native/adapter equality again from saved receipts when preparing/collecting.

The16 calls are bound/detached x La/Ca for GGR_1GLG, ALPHA_6IP9, PQQ_1H4I and
PQQ_4MAE. All five full-chain scores are then assessed together. Required
differences remain:4MAE-1H4I, alpha1F6S-GGR and alpha6IP9-GGR, each>0.02kcal/mol.
No favorable site, core, charge, water inventory or decision band is selected.
All cases are consumed development evidence; alpha structures share one group.

Use the existing intact prepare command with explicit --edge-equivalence,
--core-qualification and --full-qualification collections. The source reference
may be verified native CPU or native H200; hardware stays explicit in receipts.
The energy_evaluation.forward field identifies the native ScaleShiftMACE energy
expression; model.execution_adapter records the batched interaction execution.
The original unbatched implementation and all failed attempts remain intact.
No newDFT, solver, training, gradients or relaxation; report total real costs,
unavailable cases and baseline protection. No production promotion is implied.
