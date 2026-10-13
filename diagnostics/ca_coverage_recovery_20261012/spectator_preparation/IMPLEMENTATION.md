# Explicit spectator preparation: implementation ready

The narrow adapter is implemented and frozen for Q9's five real two-Ca sources.
Root submitted job **2969**, six shared CPUs on node-64-768g-5; its completion
wake is armed. Preparation results are pending. No energy/model evaluation is
part of this job, and no production or original source result is changed.

The source geometry audits already show a distinct PQQ Ca and remote structural
Ca in all 15 new Q9/C5/Q88 folds. This adapter tests whether those sources remain
chemically compatible with the existing local scorer. It preserves both ions
and all source heavy coordinates, then admits only a fixed local context whose
nearest retained physical heavy atom is farther than 4.3 Å from the spectator.
The spectator is explicitly absent from the local energy model, not silently
assigned a charge or parameter.

The implementation reuses the exact core, context, canonical union, physical
mapping and state checks. Only the standard-only protonation wrapper's allowed
source inventory is extended. Its old RNG, pH, OpenMM CPU-one-thread setting,
forcefield=None and nonstandard-H policy remain. Actual per-sample selectors
handle Q9 sample3's reversed active/spectator chains correctly.

First, a fresh P12293 sample0 one-Ca replay must recover the archived core and
context coordinates and expanded maps/state within the existing 1e−12 Å rule.
Then five independent spawned source workers prepare Q9. Failed regression
blocks new preparation; individual source failures remain in the denominator.
All scientific imports, preparation and tests are guarded by the actual Slurm
allocation. Only metadata freeze, Python AST and shell syntax checks ran on login.
Root and Gauss reviewed the code statically; this is not a passed molecular test.

## Exact state

- Frozen manifest: `workspaces/ca_coverage_recovery_20261012/spectator_preparation_q9_v2/manifest.json`.
- SHA-256: `b79a1b21f9e7d348154e34ea8a63966429f9511fef36537bedcb2cf55aa177dc`.
- Expected result: its `execution/FINAL_PREPARATION.json`.
- `READY.json` contains all driver/script/plan pins and the runnable submission.
- The unexecuted v1 metadata is preserved; v2 incorporates the strict `<4.0 Å`
  association convention, explicit spawned workers, and actual standard partition.
- Original `PLAN.md` and `AUDIT.json` remain the pre-implementation record.

The same code can prepare the ten actual C5/Q88 sources after root interprets
this first result, reusing the successful regression receipt rather than repeating
it. Subsequent scoring belongs to the separate two_ca_scoring adapter. A new
occupancy source never silently replaces its original one-Ca benchmark row.
