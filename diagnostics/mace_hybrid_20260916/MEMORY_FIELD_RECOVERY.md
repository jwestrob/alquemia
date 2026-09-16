# Field-response tensor-product workspace

Continuation of the approved memory work and preceding technical recovery
notes. Full job1200380 passed all four core comparisons and reached the learned
field update, where `SparseUvuTensorProduct` attempted a17.85GiB per-atom outer
product while5.97GiB was already allocated. Earlier neighbor/product bottlenecks
were traversed successfully. No scientific output was available for this full
endpoint; the failure is retained in blocked_v5.

Apply the same256-atom checkpointed batching to the original sparse tensor
product module, including its uses in field updates and the energy readout.
These products act independently on the leading atom dimension. Their global
input fields are still evaluated over all atoms; global charge equilibration
still takes place after collecting all node outputs. No interactions, charges,
weights, precision, self terms or physical settings change. Core verification
uses17-atom blocks and the original fixed tolerances before full execution.
New version blocked_v6 retains all prior outputs and implementation snapshots.
