# Frozen-density/GK hybrid operations

The three recorded stages preserve the default scorer. Read the corresponding
plans/reports before interpreting any numerical output. All paths below are
explicit; run from any directory. Existing completed calculations are reused.
New collection paths are write-once: use another explicit output path if the
shown replay filename already exists. No source edits are required.

## Audit and collect the actual density/boundary inputs

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_multipole_coupling_v1/implementation/mace_density_multipoles.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_multipole_coupling_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_boundary_v1/implementation/mace_density_gk_boundary.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_boundary_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/implementation/mace_density_gk_hybrid.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/manifest.json
```

Density job1201063 already completed all8saved-density utilities. Boundary
initialization completed all16states locally. Neither needs a scientific rerun.
These validators inspect pinned software, method, geometry, prepared charges,
source observations and exact transformation of the native solver routine.

## Execute/recover and collect the declared hybrid inventory

```bash
sbatch --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/slurm_%j.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/slurm_%j.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/mace_omol_20260917/run_density_gk_hybrid.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/implementation/mace_density_gk_hybrid.py collect --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/manifest.json --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_v1/collection_replay_01.json
```

Submission uses64CPUs/64GB and the existing standard/memory policy; it acquires
an execution lock. Do not submit while another executor holds that lock.
Successful tasks are reused after hash checks. Failed attempts remain visible;
`execute --retry-failed` is an explicit recovery operation inside the same
allocation and recorded manifest, not permission to alter scientific settings.
Collection includes separate primary, rigid, radius and convergence results,
all component energies, charge/state provenance, missing statuses and checks.
No reference or calibrated classification is supplied by this development pilot.

## Rebuild and prepare without launching science

The following fresh output names deliberately preserve the completed build and
manifest. Preparation itself runs no energies, response solves, DFT or MACE.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/scripts/mace_density_gk_hybrid.py build --parent-software /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_boundary_software_v1/receipt.json --native-source /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/tinker_sources_v1/tinker_git/source/induce.f --frontend /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/scripts/mace_tinker_density_gk.f90 --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_rebuild_01
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/scripts/mace_density_gk_hybrid.py prepare --boundary /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_boundary_v1/result.json --density /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_multipole_coupling_v1/report_job_1201063/result.json --short-report /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/explicit_field_short_result_v1.json --software /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_rebuild_01/receipt.json --plan /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/mace_omol_20260917/DENSITY_GK_HYBRID_PLAN.md --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/density_gk_hybrid_reprepare_01
```

The isolated frontend preserves the installed native library. Original and
derived native source, compiler command/version, executable and build receipt
are retained. The derivation only supplies four direct-field arrays and exports
the final unrounded native residual; solver equations remain unchanged.

## Tests on real artifacts

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_density_multipoles.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_density_gk_boundary.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_density_gk_hybrid.py -v
```

Actual native integration tests skip explicitly when their real outputs are
unavailable. Algebra/source/cache/error tests do not manufacture successful
scientific records. Costs for new solver work must be combined with the
separately recorded upstream quantum, charge, field and MACE costs before any
production feasibility claim.
