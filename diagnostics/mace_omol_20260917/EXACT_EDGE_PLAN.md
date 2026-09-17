# Exact interaction batching for OMOL energy-only inference

Declared2026-09-17 before adapted inference. The native intact ALPHA_1F6S
call1200808 failed CUDA OOM onA5000; unchanged H200 recovery1200809 is queued.
Jacob's blanket pilot/resource authorization applies. This is engineering
equivalence work, with no new physical model or altered benchmark decision.

The actual checkpoint has three RealAgnosticResidualNonLinearInteractionBlock
layers. Installed MACE0.3.16 blocks.py computes learned radial weights and
edge density, tensor products for every neighbor edge, then sums into receiving
atoms. The failed native tensor product tried an additional7.20GiB with18.41GiB
already allocated. No allocator setting could remove that live-memory demand.

Process edges in original order,1024 or2048 edges at a time. Use the original
embedding, radial MLP, tensor-product weights and cutoff on every edge. Sum
both messages and density into the full original atom inventory before applying
the original normalization, residual, nonlinear activation and linear output.
No edge/atom is discarded; no new carve, cap, graph partition, charge choice,
precision, learned parameter, energy correction or long-range term is added.
All global charge/spin conditioning and graph/node readouts remain untouched.
Only floating-point summation/batch order can differ. Unsupported interaction,
fused/LAMMPS execution or enabled gradients must fail explicitly.

Use typed execution-adapter metadata, immutable implementation snapshots and
cache keys. Never let an adapted result satisfy a native task silently. Model
parameters and buffers must retain their versions and identity; actual layer
and edge counts/chunk size must be recorded. This adapter supports energy-only
inference; no force or response claim.

Finite qualification inventory:

1. Four exact existing1H4I/4MAEcore endpoints at each of two chunk sizes:8calls.
   Compare every energy and paired Ca-minus-La contrast with actual native
   energy-only bridge1200807. Also retain native component sum checks.
2. Conditional on step1, replay all14 existing intact ALPHA_1F6S qualification
   geometries at1024edges/chunk. Same repeat/rotation/farther-distance and
   zero-detached-edge checks as INTACT_CHAIN_PLAN. Compare each to matching
   native H200 output when available, retaining unavailable status meanwhile.

Require all endpoint and paired-contrast differences<=0.01kcal/mol, frozen
before results. No comparison-driven tuning of chunk size or tolerance.22new
forwards planned;failed attempts/retries remain visible;zeroDFT/training/solver.
First use oneA5000/16CPU/64474MiBhost with measured actual costs. Keep native
H200 job and all native results. Report equivalence separately from predictive
usefulness; passing this work does not complete the discriminator goal.

No altered scientific score will be accepted as a memory fix. Larger intact
benchmark execution with this adapter requires successful core qualification,
intact numerical checks and native intact equivalence on this declared case.
