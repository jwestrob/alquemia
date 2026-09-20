# Fixed torsion intervention operations

From the repository root. All outputs are write-once. Source/chemistry/mapping
preparation, 128 native GFN2 tasks and 16 native DFT validation tasks are frozen
under `workspaces/accommodation_torsion_20260920/prepared_v3/`.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_accommodation_torsion_profiles.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/accommodation_torsion_profiles.py validate --manifest workspaces/accommodation_torsion_20260920/prepared_v3/manifest.json
```

`prepare --help` lists explicit source-root, parent-plan and agreement paths for
a new copy; the approved input cases and physical angles are fixed in this
experiment version. Use a new version for a different scientific scope.
The CPU and MACE sbatch scripts require an absolute manifest path; actual
submitted commands/IDs are recorded in `SUBMISSIONS_PRIMARY.json`. Reuse valid
completed outputs rather than submitting them again.

The independent analyzer accepts `--manifest`, `--mace-result`,
`--low-collection`, optional `--dft-collection` and a fresh `--output`. It exports
all finite work, differential work, donor distances and DFT comparison errors.
Missing work stays null. No fitted affinity correction is produced.
