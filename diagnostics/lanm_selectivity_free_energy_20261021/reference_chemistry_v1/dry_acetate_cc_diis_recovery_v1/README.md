# One failed Dy CC cell: MaxDIIS12 recovery handoff

Source-only preparation is complete. The sole eligible cell is `cc__Dy__rotation__plus_005__complex`
(original failed job18866). Native inputs differ only by ` MaxDIIS 12`
inside `%mdci`; see [DECK.diff](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/preparation_v1/DECK.diff). Actual original readbacks
are MaxDIIS7, MaxIter100, residual tolerance1e−8 and amplitude level shift0.2.
The original CC trace ends at1.38e−6 after100iterations, with no final CCSD(T).
No original file, gate, method, geometry or failure record was replaced.

The exact current v2 canonical parser and deck are reused; the recovery parser
adds printed solver readback and retains failed native SCF/CC evidence. Only
one immutable fresh output is allowed. A successful scalar energy cannot
by itself complete or qualify the original eight-cell rotation comparison.
Original printed SCF criteria, frozen-f/single-reference and basis limits remain.
The source makes no assumption about DIIS vector storage or a qualified restart.

## Existing zero-call preparation

11 focused checks pass, including raw failure/source joins, one-line input
delta, unauthorized cell/path rejection, and changed/missing solver readback.
Synthetic solver text is labeled as synthetic and does not establish native
keyword support or convergence. Actual cost: 2.6139s wall,
2.4599CPU-s, one thread, zero native calls.

[Manifest](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/preparation_v1/MANIFEST.json) and [preflight](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/preparation_v1/PREFLIGHT.json) are frozen.
The old job's native wall×ranks proxy is 19735.592CPU-s;
this is historical failure cost, not a prediction for DIIS12 or scheduler usage.
Root must inspect actual native/runtime/memory evidence before any continuation.

## Practical commands — unexecuted

Optional new zero-call preflight (the saved preparation_v1 already passed):

```bash
env OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/REFERENCE.py prepare --out /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/preparation_v2
```

Root binds the [release template](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/ROOT_RELEASE_TEMPLATE.json) to
`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/run_v1/ROOT_RELEASE.json` and creates that parent directory. The
template status is null; no active release is prepared. Then the exact scout
command, if root elects execution, is:

```bash
sbatch --job-name=dy_cc_diis12 --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/run_v1/scheduler.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/run_v1/scheduler.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/donor_cp_v1/RUN.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/REFERENCE.py run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/preparation_v1/MANIFEST.json --cell cc__Dy__rotation__plus_005__complex --release /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/run_v1/ROOT_RELEASE.json --out /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/run_v1/execution
```

The reused wrapper supplies4 single-thread MPI ranks/one node, standard_expand,
and no memory/time override. Runtime MaxCore is floor(0.75×actual memory_MiB/4).
The adapter has no loop/retry/remaining-cell path and launches no watcher.
[COMMANDS.json](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_cc_diis_recovery_v1/COMMANDS.json) retains the structured argv and future
saved admission command. No user approval step is introduced.
