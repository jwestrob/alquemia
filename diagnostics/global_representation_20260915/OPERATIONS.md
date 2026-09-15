# Runnable operations

From `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs`. The following are preparation/integrity operations; use a fresh output path for regeneration.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_global.py audit --solver-result workspaces/affordable_challenger_20260915/solver_completion/solver_result.json --output workspaces/global_representation_20260915/audit_rerun.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_global.py prepare --solver-result workspaces/affordable_challenger_20260915/solver_completion/solver_result.json --carve-manifest diagnostics/pqq_boundary_pair_20260914/mxaf_1H4I_qm33/mxaf_1H4I_qm33_carve_manifest.json --orca /groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca --output workspaces/global_representation_20260915/preparation_rerun
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p 'test_affordable_global.py' -v
```

Prepared scientific inputs: `workspaces/global_representation_20260915/prepared_v1/global_manifest.json`. **New Hamiltonian proposed, not executed.** Once the scientific-method decision is agreed, this uses the existing runner inside a suitable allocation:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/run_orca_task_manifest.py workspaces/global_representation_20260915/prepared_v1/global_manifest.json --orca /groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca --workers 2 --nprocs 16
```

The runner refuses login-node execution and writes per-attempt receipts without overwriting failed attempts. Sixteen MPI ranks per endpoint follows the existing runner policy; actual native-xTB parallel behavior must be checked at launch. Choose an appropriate node for this two-endpoint workload, not a 344-CPU allocation with one serial worker. No new job was submitted, no walltime requested, and no watcher was started during this investigation. Execution must retain parameter files, outputs and charge/wavefunction artifacts required to assess the method.

## Approved execution

Jacob agreed to continue after the explicit method proposal. Job **1198999** runs the pair through `run_global.sbatch`, which divides the allocated CPUs equally between the two endpoints (172 MPI ranks each on node-344-8t-1). This overrides the small-carve 16-rank default for this whole-protein feasibility test; actual MPI activity was verified. Both exported native parameter sets contain the expected La/Ca entries and give 25,764 active electrons and 23,259 orbitals. No converged energies yet.

Collect completed/failed endpoint receipts without hiding partial failures:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_global_collect.py --manifest workspaces/global_representation_20260915/prepared_v1/global_manifest.json --output workspaces/global_representation_20260915/collection_1198999.json
```

Owned terminal watcher PID 450588 writes `workspaces/global_representation_20260915/terminal_1198999.json`. It monitors only this job and launches no new analyses.
