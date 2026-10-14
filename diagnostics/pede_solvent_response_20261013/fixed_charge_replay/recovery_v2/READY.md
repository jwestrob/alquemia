# CMake integration recovery v2 — prepared, unexecuted

Slurm4000 failed during CMake generation: exporting xTB's libraries through an
`ALIAS` exposed the local `mctc-lib` target outside its export set. Actual
`configure.err` and terminal collection are pinned in the successor manifest.
There were zero native component evaluations and zero SCC/DFT calls. The
libcurl version warning was nonfatal; the export error stopped generation.

The sole build fix replaces the mctc-lib alias with an imported interface linked
to the same target, matching xTB's own `cmake/modules/xtb-utils.cmake` subproject
pattern. Native source archives, compiler/library pins, Fortran wrapper, driver,
eight endpoint inputs, parameters, charge precision, tolerances and accounting
are unchanged. Failed v1 source and execution artifacts remain intact.

Manifest:
`workspaces/pede_solvent_response_20261013/fixed_charge_replay_v2/manifest.json`

SHA256 `d1392c78cf3c0446a73ad4337db105f9ab371d10d2932b0586023ba223820689`.

```bash
REPLAY_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/fixed_charge_replay_v2
REPLAY_SHA=d1392c78cf3c0446a73ad4337db105f9ab371d10d2932b0586023ba223820689
sbatch --nodelist=node-64-768g-5 "$REPLAY_WORK/implementation/RUN.sbatch" \
  build-test-execute "$REPLAY_WORK/manifest.json" "$REPLAY_SHA" \
  "$REPLAY_WORK/execution_v1"
```

Same12-task shared request, no memory/time override; at most12 compiler workers
and one numerical thread per component evaluation. Root submits and arms the
completion wake on `execution_v1/FINAL_COLLECTION.json`. No build/test/replay
has run for v2 yet. The four native water-charge component gates remain mandatory
before the four vacuum-charge evaluations. No successful calculation is repeated.
