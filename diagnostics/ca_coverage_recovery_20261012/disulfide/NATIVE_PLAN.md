# One repaired MMOL_2048 source: unchanged native accommodation

Slurm2870 repaired only the malformed existing Cys128 HB2 coordinate and passed
the original envelope/state/map gates. Both real controls replayed. This stage
asks whether the recovered input can now produce the same finite native geometry
pool used by the corrected Ca-conditioned benchmark. It does not score a new
biological label, recalibrate, or change production.

One seven-CPU warm model worker uses the already qualified CPU MACE-OMOL weights,
head, charge/spin convention, force and precision settings. The adapter invokes
the pinned CPU_BATCH.run_source unchanged: two origin energies/forces, the
original force-selected four-mode selector, two original SLSQP searches,
and at most two opposite-metal candidate evaluations. Existing original
ftol/maxiter/bounds/admission/tie settings stay unchanged. Both metals share the
complete candidate pool. Failed cells/searches stay unavailable. A single
exclusive-created attempt prevents accidental repeated partial searches.

Before new source evaluations, replay the same two actual archived CPU/GPU
qualification fixtures under the existing0.01 kcal/mol energy and0.001 eV/Å
maximum-force tolerances. These are the existing runtime fixtures, not the two
same-protein preparation controls; those are not rescored. Recheck actual source
pair, q0 physical map, charges, original optimizer binary and code pins on Slurm.

Task manifests contain only the repaired target and two runtime fixtures.
The existing run_source slices its first212 task slots but filters by the
target case ID; the fixture IDs differ and are not admitted to the source pool.
The worker authorizes requests against the new exact two-origin manifest.

Finite work: two qualification calls, two new origins, two unchanged bounded
optimizer starts and at most two new cross calls. No DFT, xTB, new waters,
protons, source folds, training or calibration. GFN2 composite selection remains
unavailable until a separately collected corrected solvent stage. Preserve
original125 counts and the original unavailable row in every later overlay.

Record all actual model forwards, terminal and interrupted calls, wall/CPU
costs and placement. SevenCPU Slurm shared allocation, no memory/time flags.
Root submits and arms the existing completion watcher.

```bash
sbatch --parsable diagnostics/ca_coverage_recovery_20261012/disulfide/run_native.sbatch "$PWD/workspaces/ca_coverage_recovery_20261012/disulfide/native_v2/manifest.json"
```
