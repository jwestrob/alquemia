# Prepared; scientific execution awaits root review

2026-09-19 20:54 UTC. No new scientific energy calls or Slurm submissions.

`PLAN.md` records the exact five-case/ten-endpoint proposal sent to the parent
agent. The parent requested review before execution; no review response has yet
arrived. This is a delegated coordination checkpoint, not a reinstated standing
user approval gate. The baseline/default and existing outputs remain unchanged.

Prepared immutable manifest:
`workspaces/coordination_preparation_20260919/prepared_v2/manifest.json`.
The ten original-core and compact-context mappings pass covalent-bond and
analytic-Jacobian checks. Physical coordinate dimensions are 10/12/16/11/11 for
1H4I/4MAE/1GLG/1F6S/6IP9, respectively, in both endpoints.

Six tests passed, zero skips, in 6.015 seconds:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
-m unittest discover -s tests -p 'test_coordination_preparation_context.py' -v
```

These are actual pinned geometry/provenance/parser tests, not scientific energy
integration results. They cover ten paired maps, real archived starting energies,
constraint Jacobians, fixed waters/source atoms, and explicitly corrupted cap and
out-of-domain copies. No manufactured energies or labels were used.

Unexecuted `prepared_v1` is preserved: preparation initially required bitwise
identity through floating-point zero-angle rotation, which differs at ~machine
precision. The corrected check uses absolute tolerance 1e-10 Å, while source inputs
remain exact archived bytes. No scientific job was attempted for either version.

Next command, after root reviews the already submitted exact plan (from repo root):

```
sbatch --parsable \
  --output=diagnostics/coordination_preparation_20260919/mace_%j.out \
  --error=diagnostics/coordination_preparation_20260919/mace_%j.err \
  diagnostics/coordination_preparation_20260919/run_mace.sbatch \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/coordination_preparation_20260919/prepared_v2/manifest.json
```

The current live script and frozen snapshot match. After actual MACE proposals,
`prepare-dft` creates the finite original-core SP manifest, the existing
`affordable_workflow` executes it, and `collect` reports before/after raw contrasts
and group differences. Do not inherit baseline bands or call lower energies an
accuracy improvement. All files remain uncommitted for parent coordination.
