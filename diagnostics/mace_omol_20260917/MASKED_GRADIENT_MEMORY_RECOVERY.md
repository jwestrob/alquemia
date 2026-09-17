# Whole-protein derivative memory recovery

Declared after actual full-stage job 1200884 failed on its first GGR Ca
endpoint, before a replacement adapter is evaluated. The original gradient
plan's structures, coordinates, weights, mask, spin, precision, signed
0.01-A displacements and all numerical acceptance criteria remain fixed.
Production energy-only adapters and scores are unchanged.

## Observed failure and specific technical change

Adapter v2 passed all 23 real-core derivative checks. Whole GGR (4698 atoms)
then exhausted its A5000 allocation during an edge tensor-product forward:
24,969,669,632 peak CUDA allocated bytes. This was one failed model call;
no whole-protein derivative or displacement result was produced.

Installed PyTorch 2.8 inspection using the real 73-atom core coordinate tensor
shows IndexAddBackward saves the source tensor; ScatterAddBackward saves only
the index. The current per-edge checkpoint recomputes each message, but the
subsequent index_add retains those messages outside the checkpoint. This is a
specific memory-retention mechanism; actual recovery must establish whether
removing it suffices.

Adapter v3 replaces the two index_add accumulations with scatter_add using a
broadcast view of exactly the same receiver indices. Every message, edge order,
receiver and density normalization is retained. No tensor needed for a
coordinate derivative is detached, no precision is reduced, and no source is
removed. No custom derivative is introduced.

## Qualification and execution

Requalify the new adapter with the same ten real-core tasks (native versus
adapted gradients/energies, rotation, signed displacements). Only after those
pass, rerun the same six-task whole-GGR manifest. New versioned directories,
hashes, receipts and costs; previous failures remain. These sixteen forwards
are additional recovery work, not hidden in the original sixteen-call count.
One A5000, 16 CPUs, 64474 MiB host memory, existing runner/environment. Actual
memory and time determine feasibility; no physical-force or relaxation claim.
