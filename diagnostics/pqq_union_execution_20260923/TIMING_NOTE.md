# Timing interpretation correction before results

The frozen plan describes the released path as a warm two-endpoint batch per
source. The actual pinned implementation is more specific: one process per
source, but `mace_omol.worker` constructs a calculator independently for each
endpoint. This gives up to20 model initializations across ten processes. Its
`model_load_seconds` includes initialization and input setup, not a separately
isolated checkpoint-read timer. The candidate uses three persistent stage workers
(origins, search, cross-scoring) and records model initialization in their actual
worker summaries.

No executor, physics, batching or numerical setting has changed for this test.
Report the actual initialization/evaluation/stage receipts and complete sequential
same-host allocations. Any observed end-to-end difference includes preparation,
batching, initialization, rank/concurrency and scoring/collection overhead. It is
not a controlled speedup attributed solely to accommodation or MPI count, and is
not a comparison with a hypothetically optimized static implementation.
