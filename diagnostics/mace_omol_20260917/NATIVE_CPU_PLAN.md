# Native CPU equivalence reference

Declared2026-09-17 before native CPU inference. The H200 native intact reference
job1200809 remains queued. Exact edge batching passed20 native core checks;
intact GPU qualification is being prepared. Native CPU evaluation is an
independent way to obtain the same unmodified full-system reference without
waiting for the large GPU. Keep all existing jobs and receipts.

Use the exact pinned checkpoint, native float64 ScaleShiftMACE forward under
no_grad, original full neighbor graph, charge/spin and energy terms. No edge
adapter, solvent, model fitting, force computation or physical input change.
Record CPU device explicitly; native CUDA and CPU execution cannot silently
satisfy each other's tasks. Retain the same checks on state and readout sums.

Finite inventory:4 native core calls (1H4I/4MAE x La/Ca), compared with actual
energy-only GPU bridge1200807. Require all four endpoint and two paired errors
<=0.01kcal/mol. Conditional on that gate,14 ALPHA_1F6S calls on the exact
existing bound/detached/repeated/rotated/farther geometries. Apply the original
intact numerical checks, then compare all14 energies and their paired contrasts
against the exact-edge A5000 outputs, using the same0.01kcal/mol tolerance.
18 new forwards, zeroDFT/training/solver; failures/retries remain recorded.

Allocate one standard node with64 CPUs and128GiB RAM. Native GPU failure needed
approximately26GiB of live tensors for this case; host peak and CPU scaling
remain measurements to make, not performance guarantees. Use the existing
task executor with a CPU allocation wrapper. All actual allocation and failed
attempt costs count. This is one-time engineering validation, not a claim that
routine scoring is affordable by moving GPU work onto CPU nodes. Production
remains unchanged, and no predictive claim follows from numerical equivalence.
