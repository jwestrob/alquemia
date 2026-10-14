# Joint-search reuse: technically feasible, not yet executed

The existing `adaptive_metal_proposals.MixedProposal` consumes the same
`transport.evaluate(request_pin) -> (native_result, result_pin)` contract as the
qualified seven-thread `CPU_BATCH.WarmCPU`. Its `optimize(manifest, task, transport)`
uses the saved physical mapping, three metal translations, four specified donor
angles, analytic cap-aware gradients and existing source/overlap/displacement
checks. No model retraining or second mechanics implementation is needed for
that interface.

This is a static compatibility finding. A new research manifest must bind the
same physical source, state, model and worker, and explicitly declare the mixed
coordinates and numerical settings. The worker checks those manifest settings
against its own; simply replacing a function in an old batch is not sufficient.
The historical joint `execute()` entrypoint hardcodes a32CPU/GPU allocation and
must not be used for the proposed seven-CPU transport.

There is a numerical distinction to retain: the archived joint kernel minimizes
in eV with `ftol_eV`, whereas the current donor-only `adaptive_completion` kernel
uses a consistently scaled Hartree-equivalent objective, gradient and tolerance.
Consequently an unchanged call to the old joint kernel is its own method
version. It cannot silently be described as the current optimizer with three
coordinates appended. Any future integration should specify the scaling and
test it explicitly on the actual source before launching a finite search.

Relevant source paths:

- `scripts/adaptive_metal_proposals.py`: `MixedProposal`, `optimize`, `constraints`,
  `physical_status` and `final_geometry`.
- `workspaces/motion_envelope_20260923/pilot34_searches_v1/implementation/`:
  exact archived joint/donor kernels and `mace_site_kinematics.py` used by the
  current mappings and physical preparation.
- `diagnostics/orca_parameter_audit_20261012/cpu_proposal_batch/CPU_BATCH.py`:
  `authorize_request`, `worker`, `WarmCPU`.
- `diagnostics/exaa_metal_position_20261014/scoring/SCORE.py`: actual reuse of
  that qualified seven-thread worker for the eight completed transplant calls.

This branch adds no optimizer, worker, molecular calls or future-run manifest.
The saved-force result remains the decision input.
