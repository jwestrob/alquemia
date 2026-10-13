# MxaF Ca-conditioned technical recovery

## Question and status

Can two P12293 La searches interrupted by Slurm `NODE_FAIL` complete under the
unchanged proposal policy? Jacob authorized these separately named restarts as
part of recovering the 18 unavailable Ca-conditioned inputs. They are technical
recovery, not new structures, labels, calibration, or a production promotion.

**Root submitted job 2869 on node-64-768g-5; actual Slurm scientific preflight passed.**
Exact task/selector reconstruction and all four q0 maps passed (maximum coordinate
residual 8.9e-16 Å). Runtime qualification and recovery outcomes remain pending
at this checkpoint. The ready-to-launch preparation record remains immutable.
Root owns submission, the completion wake, corrected GFN2, frozen-band scoring,
and integration with the unchanged 125-source primary result. This runner owns
only samples 1 and 4 and returns a separate two-source native collection.

## Actual reusable evidence

Both targets have complete Ca/La origins, an admitted Ca proposal, Ca at that
proposal, and La at that proposal: eight native cells and two Ca proposals total.
Their La `INTERRUPTED.json` records explicitly preserve incomplete SLSQP searches
without claiming optimizer continuation. All old source, request, result, force,
interruption, and partial-trial artifacts remain read-only.

The recovery reuses both exact La q0 energies/forces. It restarts SLSQP from the
original q0 in fresh paths; interrupted trial evaluations are not treated as a
restartable optimizer state. Any repeated trial is a newly measured call. The
primary collection and its prior costs remain pinned alongside the new results.
Incomplete old calls retain unknown cost where no terminal receipt exists.

## Frozen method and finite work

Protocol: `MxaF_Ca_conditioned_same_policy_La_q0_restart_v1`.

- Same source/context graph, charges, singlets, atom order, mapping, selected four
  angular modes, and origin coordinates as the original prepared sources.
- Same native vacuum OMOL checkpoint `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`,
  CPU float64 calculator, SciPy 1.17.1 wrapper/kernel, and archived implementation.
- Same SLSQP: 200 iterations, ±0.8 rad, heavy-atom displacement ≤0.8 Å,
  Hartree-equivalent objective scaling and `optimizer_ftol=1e-8`. No added starts.
- Two original real CPU qualification controls must pass the existing
  0.01 kcal/mol energy and 0.001 eV/Å force tolerances before either search starts.
- Two new La optimizations; zero Ca optimizations and zero origin evaluations.
  At most two new Ca-at-La candidate evaluations complete the common pool.
  No DFT, GFN2, folding, state changes, or new geometry policy in this runner.
- Origin/Ca/La candidates retain the original 1e-12 Å coordinate alias policy.
  A missing cell keeps the common pool unavailable while valid origins remain.

One warm seven-thread model worker serves both sources sequentially. Its
one-thread optimizer parent and seven-thread model execute numerical work in
turn. Actual Slurm allocation and affinity are checked before scientific imports;
no memory request or wall-time limit is introduced. Worker timings, actual model
forwards, process CPU time, peak RSS, and total allocation time are retained.
There is no predicted successful runtime or hidden automatic retry.

## Inputs, commands, outputs

Frozen runnable plan:
`workspaces/ca_coverage_recovery_20261012/mxaf_v2/manifest.json`
(SHA-256 `14a14cbd8acb561fd255e582043928bed378a2676d64f831e7b7b8deedd26037`).

Root submission command, from repository root:

```bash
sbatch diagnostics/ca_coverage_recovery_20261012/mxaf/run.sbatch \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/mxaf_v2/manifest.json
```

The final collection will be
`workspaces/ca_coverage_recovery_20261012/mxaf_v2/FINAL_COLLECTION.json`.
It contains only these two recovery sources, the pinned primary collection,
per-cell raw receipts, pending/failure states, and incremental cost. It is an
input to root's recovery overlay; it does not replace the primary 106-source
native collection or the 125-source result.

`mxaf_v1` is preserved metadata-only preparation before the final receipt
hardening. It was not executed and must not be submitted. The v2 driver and
SBATCH pins are in `PREPARATION_RECORD.json`.

## Checks at preparation and submission

Stdlib-only preparation verified the original source results against the terminal
collection, exact task/settings/map pins, eight raw native cells, two Ca proposal
receipts, and both interrupted La searches. Python AST and Bash syntax checks
passed. An independent static review found no launch blocker. No model, NumPy,
SciPy, force projection, geometry admission test, or scientific parser test ran
on login. Exact scientific reconstruction and q0 geometry checks subsequently passed
inside job 2869; calculator qualification also runs inside that job. Preparation is not scientific
acceptance; failed qualification cannot silently substitute an origin score.

Submission/wake receipts: `mxaf_v2/SUBMISSION.json` and `WATCH.json` beneath the
workspace above. Root armed event `ca-coverage-mxaf-restarts-v2` (watcher 3490926).
Do not resubmit the ready command while that job or its results exist.
