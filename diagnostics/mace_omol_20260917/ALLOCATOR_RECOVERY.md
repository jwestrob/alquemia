# Unchanged-model allocator recovery — 2026-09-17

Job1200815 completed eight endpoints, then failed on the first 9088-atom PQQ
endpoint. Its outputs and exit collection are retained. Both alpha structures
order above GGR under the frozen relative criteria; PQQ remains unavailable.

Retry the same immutable16-task manifest and1024-edge adapter on the same A5000
allocation. The existing executor reuses eight accepted tasks and preserves the
failed attempt. Only PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True changes.
No geometry, method, precision, task inventory, score or criterion changes.
Record the exact environment and job in submission_allocator_recovery.json.

The failed allocation requested5.55GiB with14.34GiB allocated and4.32GiB
reserved but unused. This suggests allocator fragmentation may contribute;
it does not establish that the full forward will fit. PyTorch2.8 documents
expandable_segments as an experimental allocator option to reduce unusable
memory slices. If it still fails, retain that result and investigate exact
atom-local batching or the already authorized larger GPU.

Source: https://docs.pytorch.org/docs/2.8/notes/cuda.html#optimizing-memory-usage-with-pytorch-cuda-alloc-conf
Native/batched energy equivalence remains the required criterion for any new
arithmetic implementation. An allocator setting does not modify that arithmetic.
