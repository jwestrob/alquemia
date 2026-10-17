# Electronic substitution: preserved first attempt

Job **12770** left the full saved comparison unavailable. It admitted 20 of
36 LBT cells and all 10 Mex cells, and reproduced the previous comparison.
The 16 unavailable LBT cells all failed at `KeyError('charge_closure_e')`:
their archived energy-only receipts store signed charge-closure evidence in
`solver_observation` instead of the gradient receipt's top-level field.
This failure occurred while reading saved evidence; no molecular call failed
or was rerun.

All 15 declared menus remain in the original output. Only three origin menus
are available, so the partial output does not establish the full comparison
of electronic models on the same geometry pools. No favorable subset is used.

The separately versioned recovery verifies the actual archived population
vectors and applies the unchanged absolute charge-closure limit of 1e-5 e.
It preserves the Hamiltonians, states, geometry pools, selectors and numerical
thresholds. The original source, config and results remain unchanged.

## Actual cost

Slurm reports **FAILED, exit 1:0**, 10 seconds on 1 CPU and no GPU, with
6.676 seconds of measured CPU time. Application time was 9.809 seconds;
reported process CPU time was 6.576 seconds. There were **zero new molecular
calls, geometries or optimizer steps**. The raw scheduler timestamps are
retained without assuming a timezone conversion.

Exact input/output pins and scheduler records are in
[ATTEMPT_V1_ACCOUNTING.json](ATTEMPT_V1_ACCOUNTING.json) and
[ATTEMPT_V1_SACCT.txt](ATTEMPT_V1_SACCT.txt). The immutable scientific output
is [the original RESULT](../../../workspaces/lanm_coupled_response_20261016/electronic_substitution_v1/execution/RESULT.json).
