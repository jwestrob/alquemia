# 106-source CPU proposal batch is prepared

The finite Ca-conditioned MACE batch passed its file-only preflight and bounded
independent code review. No batch model call was made during preparation.
Root launched the frozen snapshot as worker1620438 with completion watcher1620439.
Its startup preflight and seven-thread controls precede new-source scoring.
The immutable workspace LAUNCH.json records the actual launch command/ownership.

## Method and scope

The manifest pins Curie's106 new usable sources,212 paired origin inputs,
frozen4.3 Å context policy, exact model and scientific scope. The three completed
Ca probes are reused separately; the full125 denominator keeps15 prior exclusions
and the one unsupported disulfide preparation. No new calibration is introduced.

The runner imports the348-file historical snapshot and calls its unchanged
`motion_envelope_run.build_tasks` and private precision-configured
`CompletionProposal`/`optimize`/`Kinematics` implementation. The original SciPy
optimizer runs in its original environment. A separate persistent model process
provides actual CPU native OMOL energy and Cartesian forces with the unchanged
checkpoint, head, dtype and electronic-state inputs. Two sequential parent/model
lanes use at most14 active numerical threads.

Before new points, both seven-thread workers must independently pass archived
Ca and La origin energy/force controls: four finite calls per attempt. Gates
remain0.01 kcal/mol energy and0.001 eV/Å maximum force error. The earlier completed
eight-thread scout passed five replays and both exact archived optimizer starts;
that evidence does not replace the actual seven-thread startup checks.

Each source uses two real origins, the existing four-mode selector and one
original SLSQP start per metal. Own-metal energies are reused and at most two
opposite-metal calls complete the common origin/adaptive_Ca/adaptive_La matrix.
Only genuinely admitted proposals enter it; failed optimizer trials remain
preserved. Native-only output explicitly lacks solvent/composite scores.

## Verification

- All106 source checks passed; all212 prepared state/map/coordinate inputs were
  verified. Contexts contain169–215 atoms. Maximum reconstructed q0 mismatch:
  1.7763568394002505e-15 Å, below the unchanged1e-12 Å copy tolerance.
- All348 historical module pins, model/optimizer/runtime pins and exact scientific
  settings verified. Selection remains deferred until real paired origin forces.
- Five file-only checks passed: wrong state/geometry rejection, archived receipt
  reuse, begun incomplete model-cell retention, interrupted optimizer no-restart,
  and rebuilding an absent/stale final index from106 terminal sources without
  launching any process or model.
- Independent review identified a missing-index recovery branch; the frozen
  implementation fixes it. No remaining code-review blocker was reported.

The first preflight wrapper stopped because PREPARATION and ORIGIN_INPUTS order
differed, although their106-member sets were identical. The corrected check
requires the exact prepared membership and authoritative ORIGIN_INPUTS order.
The failed preflight_v0 and successful preflight_v1 remain separate and intact.
Root selected the latter as the execution directory, avoiding duplicate checks.

## Restart and limits

Terminal sources and actual completed origin/cross cells are reused. A previously
started optimizer lacking a terminal proposal is unavailable on explicit restart;
there is no hidden extra q0 start. Partial scientific receipts stay intact.
The latest collection index can change between attempts; every attempt's final
collection remains immutable. All-terminal restart rebuilds that index with zero
model loading/calls.

Trial counts and seven-thread throughput remain unmeasured at preparation.
There is no new optimizer call/time ceiling or assumed CPU speedup. Interrupted
calls lacking a terminal receipt have unknown completed-forward cost, recorded
explicitly alongside all observed receipts. Root's separate corrected-GFN2
adapter handles later scalar cells and fixed-reference classification.

See COMMANDS.md and READY.json for the executable command and exact artifact pins.
