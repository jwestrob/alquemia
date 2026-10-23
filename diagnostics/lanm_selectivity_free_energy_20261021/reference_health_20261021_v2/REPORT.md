# Bounded reference-health review —21October2026,13:28UTC

**No new reference result is admitted.** One scheduler/native/resource observation per job; the donor fit was untouched. Running jobs remain unchanged.

| Job | Authoritative state / elapsed | Native evidence | Assessment |
|---|---|---|---|
|17260 Dy–water, explicit electrons|FAILED1:0;2:47:21,4CPUs|SCF abort after5291cycles; last macro299 error norm9.875006e−5; no engrad|Electronic convergence failure. The outer missing-file collection error is secondary.|
|17258 Dy–acetate, explicit electrons|RUNNING;4:59:14,4CPUs|Macro94→199; error norm8.143458e−5→1.957806e−4; repeated rejected steps|Computation advances, but convergence has not improved over the prior checkpoint.|
|16680 whole-peptide La scout|RUNNING;9:25:10,24allocatedCPUs/12MPI|SCF3→5; latest RMSDP231,MaxDP50400,DIISErr0.0752; iteration5575.2s|Slow concrete progress, still far from an admitted SCF/gradient result.|

## Resource evidence

17258 retains the same four native rank PIDs and substantial CPU accumulation since the corrected prior review. All four were in uninterruptible-sleep (D) state at the single new process snapshot; this is consistent with I/O waiting at that instant and does not establish a persistent stall. Slurm reports MaxRSS586500K, about234.65GB cumulative reads and19.04GB writes. The source output advanced by105macro iterations.

16680 retains all12 native ranks, each running (R), with increased accumulated CPU and zero observed swap. Slurm MaxRSS76668300K remains below its declared96000MiB allocation. The output added two complete SCF rows. CPU activity and output advancement do not establish convergence or a physically useful energy. No explicit OOM marker appears in the captured outputs.

17260 used10:41:52 actual scheduler CPU and4×10041=40164allocatedCPU-seconds. The native output explicitly terminates in LEANSCF; the adapter reports failed_collection/FileNotFoundError withreturncode0. Neither the return code nor intermediate energy overrides the native failure.

## Scientific decision value

The most relevant missing evidence is the matched explicit-electron Dy water/donor comparison. An admitted17258 result could test the acetate endpoint, but17260’s failed water anchor prevents the intended balanced water-to-acetate cycle. A separately declared root-owned electronic-state/solver recovery would be needed to resolve that missing endpoint; this review performs none. Even a completed AE/ECP comparison changes basis and relativistic representation as well as core treatment, and does not isolate f-electron or spin–orbit uncertainty. Existing La pilot/native-state caveats remain.

16680 addresses the older whole-peptide La method/model discrepancy under a version-limited omegaB97M-V/SCNL Hamiltonian. It cannot answer the Dy frozen-f or paired exchange question and is not a prerequisite for the active free-energy program. The present evidence supports neither a new accuracy claim nor treating all live work as idle.

## Exact evidence

- [Health receipt](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/RESULT.json); [scheduler snapshot](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/SCHEDULER.json); [resource snapshot](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/SSTAT.json).
- [17260 native abort](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/17260/endpoint.snapshot.out:6197); [actual failed result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/pilot_v1/tasks/ae_pbe0__water_central_Dy/execution/RESULT.json).
- [17258 latest macro/error rows](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/17258/endpoint.snapshot.out:4275).
- [16680 latest SCF rows](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_health_20261021_v2/16680/endpoint.snapshot.out:2214).
- [Corrected earlier review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/active_reference_health_v2/RESULT.json); its earlier reversed-job-label measurements remain invalid.
- [Active reference requirements](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/PLAN.md:80).
