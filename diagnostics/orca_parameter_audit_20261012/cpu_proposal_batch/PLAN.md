# Finite corrected-context Ca-conditioned CPU proposal batch

Preparation only until root reviews the actual manifest and dry-run receipt.
Root authorized this engineering stage on12October2026 after the completed
CPU force/proposal scout. No batch model calls have run from this plan.

## Scientific scope

Use exactly Curie's106 newly usable Ca-conditioned preparations under the frozen
existing4.3Å union policy. Preserve the one new preparation failure and all
historical exclusions in the source inventory. No reserved outcomes, new
protonation/geometry policy, label-selected subset or calibration tuning.

For every prescribed usable source, evaluate actual Ca/La native-vacuum OMOL
origin energy/forces, select the same four common physical angular modes with
the historical paired-force rule, and run exactly one original SLSQP start for
each metal. Preserve original checkpoint, state, atom order, caps, kinematics,
selector, optimizer binary, ftol1e-8, maxiter200, angular±0.8rad and heavy0.8Å
bounds/admission. No new call/time limit, fallback, extra start or altered
failure handling. Failed/oversized trials and unavailable proposals remain
explicit. Boundary proposals keep their actual flags, without a minimum claim.

Build each source's common origin/adaptive_Ca/adaptive_La matrix from its actual
admitted coordinates. Reuse exact own-metal results; compute only missing
opposite-metal cells, with source-specific charge/multiplicity. Use existing
1e-12Å geometry aliasing. Missing required proposals leave that source's full
pool unavailable; retain any valid origin/proposal/cross cells independently.
All MACE receipts remain native-only: solvent/DFT/composite values are absent.
Root performs the separate corrected-GFN2 stage after the MACE batch completes.

## Resource and qualification policy

Two original optimizer-parent processes, each with one numerical thread, and
two persistent model children, each with seven Torch threads. For each lane
the parent blocks during child inference and the model waits during parent
numerics: at most14 active numerical threads across both lanes, leaving two
of Jacob's16 for other light work. Single graph, float64, unchanged native
checkpoint/head/charge/spin embeddings, no GPU and no fake Slurm environment.

Every newly loaded worker first evaluates the archived Ca and La origin
fixtures from the completed202-atom CPU scout at seven threads: four total
qualification calls per execution attempt. Energy0.01kcal/mol and maximum
Cartesian force0.001eV/Å gates stay frozen. A barrier requires both workers'
two checks to pass before either evaluates a new source. A failed gate stops
the batch before new-source points. No automatic rerun of a failed attempt.

106 sources imply212 origin calls,212 declared optimizer starts and at most212
new opposite-metal MACE cells, plus four qualification calls on initial launch.
Optimizer trial counts are unknown and measured. Geometry deduplication and
actual unavailable sources can reduce later cells, never conceal exclusions.
Record every attempted/finished forward, failed call, worker/model startup,
elapsed and processCPU time, peakRSS and actual terminal source statuses.
There is no assumed thread-speedup or guaranteed campaign completion time.

## Restart and evidence policy

One foreground command owns both lanes and all descendants; a lock prevents
concurrent execution of the same immutable batch. Root arms its completion wake
when it later launches the command. Snapshot inputs and implementation before
execution. Persist owned PID/process-group identities and attempt records.

On explicit restart, verify manifests/receipts and skip all terminal sources.
Reuse completed origin/cross receipts; do not recompute them. Recover an actual
complete recorded output when postprocessing was interrupted. A started model
cell without a complete valid result remains unavailable. An interrupted
optimizer without a terminal proposal result remains unavailable, with retained
partial requests, rather than silently adding a fresh q=0 start. Continue only
unstarted stages/sources under the same finite manifest. Never overwrite
historical/partial scientific files or manufacture a successful receipt.

Root reviews the concrete manifest, schema and dry-run before any molecular
launch. Production and shared scientific kernels remain unchanged.

Each execution attempt retains its immutable final collection under attempts/.
The workspace-level FINAL_COLLECTION.json is a replaceable latest-summary index
for the existing completion watcher. It records terminal=false when interrupted
with unstarted sources. Scientific receipts and earlier attempt summaries remain
unchanged. No model is loaded when an explicit restart finds all106 terminal.
