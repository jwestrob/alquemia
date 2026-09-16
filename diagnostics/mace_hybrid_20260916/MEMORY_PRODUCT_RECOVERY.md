# Remaining neural-memory recovery

Same approved memory work, scientific model and numerical gates as
MEMORY_AGREEMENT.md and MEMORY_LOCAL_RECOVERY.md. No analysis expansion.

Full job1200372 revealed a further native operation with three large live
tensors: an8.93GiB allocation in the per-atom symmetric product contraction.
These operations act independently on each atom after neighbor aggregation.
Wrap the unchanged product module in checkpointed256-atom batches; core
validation uses17-atom batches to exercise multiple boundaries. The full
neighbor sums remain intact before these per-atom operations.

The edge accumulation adjoint now saves only receiver indices. Original
torch.index_add also retains the source edge messages, defeating activation
checkpointing; the derivative requires only gathering the outgoing gradient.
Tests on the real core neighbor graph compare both values and derivatives.

Nonreentrant checkpoints recompute the entire block: TorchScript converts
PyTorch's private early-stop exception into an opaque RuntimeError. This
technical compatibility fix retains all arithmetic and analytic derivatives.

Host offload in1200372 reached81.017GiB RSS despite the64474MiB request and
then failed on GPU memory. A cancellation was requested after the excessive
live RSS observation, but the job and watcher had already terminated. Do not
claim Slurm enforces this requested RAM as a measured process-memory cap.
No automatic repetition of that known offload failure is scheduled. Retry the
same real cores and full-system tasks natively with the complete blocking.
Only our own failed/pending jobs and watchers were managed; baseline preserved.
