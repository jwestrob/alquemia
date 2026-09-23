# CPU-host metadata comparison recovery

All110 proposals and110 cross-MACE calls completed under frozen inputs. Native
scalar jobs1211130/1211133 then failed in preflight on node-128-512g-8gpu-1,
each after3seconds, before any of the440 candidate GFN calls started. Their
outputs/allocations remain preserved. The failure was exact equality of
recomputed selector diagnostic metadata, not a changed proposal or failed SCF.

Read-only job1211145 replayed all55 source pairs on that host: zero molecular
energy/force calls,16.045326109975576s executor wall. Actual report:
`run_v2/SELECTOR_HOST_DIAGNOSTIC.json`. All selected mode IDs, source maps,
task forces/gradients, electronic states and physical settings match. All474
reported differences are repeated copies of four selected-mode diagnostic
floats, with maximum absolute difference1.0658141036401503e−14. No other task
field differs. Those fields are Ca/La/differential normalized projections and
independent geometric fraction.

Parent authorized the narrow technical comparison repair before scalar restart.
Execution adapterv2 permits absolute1e−12 rounding only in those four diagnostic
fields. Mode IDs, order, criteria, rank threshold, redundant modes, other metadata
and all physical/electronic task fields still compare exactly. The frozen selector
is retained. No optimization, MACE, coordinate, reference or physics is changed.
This is not an energy-convergence tolerance, altered force or selected-subspace
policy. Actual host replay tests distinguish this repair from a hypothesized
hardware explanation; the specific BLAS/hardware source of rounding is not claimed.

The recovery executes the original two scalar manifests only, with the same
rank1×32/32CPU/64GiB resources and existing runner. It does not repeat completed
chemistry. Both original executorv1 and new validatorv2 remain pinned. Final costs
include the two failed allocations and read-only replay allocation.
