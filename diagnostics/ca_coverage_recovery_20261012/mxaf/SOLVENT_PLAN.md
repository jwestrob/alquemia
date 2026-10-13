# Corrected solvent completion for the two MxaF restarts

**Implemented, not executed.** This dependent stage waits for the terminal native
recovery collection from job 2869. It does not rerun the original 106-source batch
or recalibrate the existing decision bands.

## Finite work and reuse

`SOLVENT.py` verifies the unchanged 125-source scope, primary result, frozen
reference, both source states/maps, and the eight unchanged origin/Ca-candidate
native cells. The primary GFN2 collection already contains eight accepted raw
inputs per source (two metals × origin/Ca proposal × vacuum/water), sixteen total.
All exact old input keys enter the reuse cache; their raw output, runtime,
parameters, charges and solver evidence are checked using the frozen helper.

Only Ca/La at a newly admitted La proposal can require new components: at most
eight new inputs across both sources, each evaluated once by primary xTB and once
by independent tblite. Coordinate aliases can reduce that count. A missing new
proposal or failed component remains unavailable; valid origin scores remain
separate. The same complete origin/Ca/La pool, selection tolerance, raw-score
algebra, canonical bands and evidence labels are retained.

Before any new calls, the stage collects reused components and verifies exact
agreement of eight prior solvent components and both origin scores/decisions
against the actual primary collection. This is a real-artifact parser/algebra
check, performed under Slurm. Qualification keeps the existing per-metal
ALPB-minus-vacuum agreement ≤1e-6 Eh, atomic charge agreement ≤1e-5 e and existing
SCF/state/charge-sum checks. Primary xTB energies remain authoritative; a known
absolute-energy offset in tblite is not subtracted from them.

## Existing implementations retained

The adapter reuses `GFN2.py` cell validation and collection, `REEVALUATE.py`
component and pool algebra, `context_factorial_gfn2.py` process execution, and
`GFN2_STAGE.py` allocation-aware execution. Their pinned files remain unchanged.
Twelve allocated CPUs support up to twelve independent one-thread input workers;
each input's two backends run sequentially. No memory request, time limit,
account switch, new retry, model inference or optimization is introduced here.
All scientific imports, preparation, replay checks and execution require an
actual twelve-CPU Slurm job before running.

`FINAL_COLLECTION.json` contains only the two recovery sources. The separate
`PRIMARY_RECOVERY_OVERLAY.json` retains each of the original 125 rows under
`primary` and appends `technical_recovery` for these two sources. Original
counts—including eighteen unavailable—remain explicit alongside the recovery
counts. The two records are never silently merged into the historical result.
Backend process costs and full stage allocation time are reported separately.

## Root launch command

After job 2869 (root must attach the dependency and completion watcher):

```bash
sbatch --dependency=afterok:2869 \
  diagnostics/ca_coverage_recovery_20261012/mxaf/solvent.sbatch \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/mxaf_v2/FINAL_COLLECTION.json \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_slurm_v1/FINAL_COLLECTION.json \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/mxaf_solvent_v1
```

Only AST/Bash syntax and independent static review occur before submission.
Scientific fixture replay and actual new component calls are unrun at this
checkpoint. No speedup or recovered classification is claimed in advance.
