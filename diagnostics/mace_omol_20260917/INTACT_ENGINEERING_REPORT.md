# Intact-chain OMOL: memory implementation checkpoint

The intact 1,932-atom ALPHA_1F6S model now fits on an A5000 using exact interaction
batching: peak allocated GPU memory4.958GB, approximately5.3seconds inference
per endpoint. The first unbatched call exhausted the same GPU. This is a
working memory improvement; no predictive improvement is claimed yet.

| Actual experiment | New successful calls | Numerical checks | Result |
|---|---:|---:|---|
| Native energy-only core bridge1200807 |4|10|All pass; exactly matches archived force-producing energies|
| Native intact A5000 attempt1200808 |0|0|One CUDA OOM; no energy; failed attempt retained|
| Batched core bridge1200810 |8|20|Both chunk sizes match native energies exactly|
| Batched intact qualification1200811 |14|31|All pass; max error2.12620406899e-6kcal/mol|
| Native CPU core bridge1200812 |4|10|All pass; exactly matches native GPU energies|
| Native CPU intact reference1200814 |14|31|All pass; all24 native/batched comparisons also agree exactly|

The adapter sums every original neighbor message and density before applying
the original nonlinear layer. Each full call processes134,760edges per layer
in the bound geometry. Original model parameters and buffers remain unchanged.
All physical atoms, charge/spin, geometry, cutoff and energy terms are retained.
Synthetic carve caps are absent. Energy-only results explicitly contain no
forces; this work does not validate gradients or relaxation.

Full-size native equivalence now passes: native CPU job1200814 completed14
calls and31 numerical checks; all24 direct native/batched endpoint and paired
comparisons agree exactly. H200 job1200809 remains queued. The remaining
16-call comparison is running as1200815 under INTACT_BATCHED_CONTINUATION.md.
No threshold, label or water state changes; predictive results remain pending.

## Actual cost at this checkpoint

44 successful forward evaluations and one failed attempt; zero newDFT, solver
or training calls. Completed allocations total497GPU-seconds,71,184allocated
core-seconds and18,440.933actualCPU-seconds. Native CPU qualification is one-time
validation, not a claim of cheap routine CPU scoring. Local preparation timing
receipts and exact job records are linked in INTACT_ENGINEERING_STATUS.json;
other local work is not fully profiled. Pending jobs have no completed cost
assigned. Model inference timing excludes import, preparation and collection.

Baseline and old experiments remain unchanged. All cases are consumed research
evidence. Intact inputs are prepared chainA; they do not imply complete assembly,
long-range electrostatics, certified separated ionic states or binding free
energies. Broad La/Ca affinity and production promotion remain unestablished.

The live report-source snapshot timing issue is documented in INTACT_REPORT_SOURCE_NOTE.md. Independent recomputation verified every saved comparison; inference always used frozen source snapshots. Further reporting uses an immutable implementation directory.

Nested checksum verification now has an opt-in operation-local cache. Four corruption/restoration tests pass, including same-size rewrites with restored mtime. Fresh files bypass reuse to handle coarse timestamp resolution. The complete comparison still agrees: replay took258.88wall/156.06CPU seconds versus370.45wall/311.72CPU seconds before caching. The scientific computation is unchanged; other validation overhead remains.

## Partial benchmark checkpoint, 2026-09-17

Job1200815 produced eight new accepted endpoints, then exhausted GPU memory at
the first PQQ endpoint. With the four qualified alpha1F6S reuses, intact-chain
R_coord values are GGR62.7148448127, alpha1F6S99.8684074321 and
alpha6IP990.8031326925kcal/mol. Both frozen alpha-minus-GGR criteria pass:
+37.1535626194 and+28.0882878798. This is one consumed, qualified biological
comparison on two alpha structural forms; no absolute calibrated class or broad
affinity validation follows. Unlike the primary core model, this representation
gets both relative directions right. Representation and total-charge conditioning
changed together, so this result does not isolate the physical cause.

The PQQ comparison remains unavailable. Unchanged-manifest allocator recovery
1200816 is running on the same A5000, reusing all eight valid tasks and retaining
the OOM attempt. Only the documented PyTorch expandable-segments allocator
setting changed; see ALLOCATOR_RECOVERY.md. No scientific protocol/criteria edits.
Final reporting will use the completed recovery collection and frozen reporter.
