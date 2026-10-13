# Actual3354 collection failure repaired without molecular repetition

The native La/water scout finished normally with a converged101-cycle SCF.
Job3354 was marked failed because its collection stage could not locate an
optional summary label. **Job3370 completed the separate collection-only
recovery with15/15 regression checks, zero new molecular calls.** The original
failed execution/collection and all molecular inputs/outputs remain unchanged.

The recovered result is
`workspaces/pede_solvent_response_20261013/reference_collection_recovery_v2/execution/COLLECTION.json`,
SHA `b9c597207d4689d02476de51e0ded5e7c49b3635dee915547a6fc3b95527a782`.
It contains one completed diagnostic and seven unexecuted cells. It is not a
complete solvent comparison or evidence of electronic-state stability.

## Exact defect and unchanged criterion

The original collector obtained the TRAH threshold from a `Last Orbital Gradient`
summary row. That row appeared in its previous SOSCF→TRAH example, but is absent
in this DIIS→TRAH output. The actual scout prints both:

```text
Converg. threshold (grad. norm) .... 1.000e-05
Orbital Gradient TolG          .... 1.000e-05
```

Its last NR macro residual is5.663000e−7. The new external collector reads the
two actual active thresholds and requires them to agree. If the optional summary
tolerance exists, that must agree too. Missing, conflicting, nonfinite or
nonpositive thresholds fail closed. The scientific criterion remains the same:
nonnegative finite final TRAH gradient norm at or below its printed threshold.
No density/DIIS summary is relabeled as final TRAH density convergence.

All ten other checks passed originally. Regression confirms every other
collection field/check is unchanged, including the energy,828-electron state,
192-atom populations,46-electron La ECP, retained wavefunction and diagnostics.
The normal output total is−6072.963646749778Eh; the macro SCF energy is a distinct
quantity before the full composite corrections. The original large density/DIIS
summary values and small-gap warnings remain available and do not establish
state stability. No new physical inference follows from correcting the parser.

## Actual regression evidence

Slurm3370 ran for approximately2s on one CPU. `REGRESSION.json`, SHA
`8216ba0cf3f7b20d7072909123360ff379f91573d06a34d2d406baf60e098553`, records:

- Current missing-summary output: exact printed threshold1e−5 and final
  residual5.663e−7; original failure reproduced and new collection admitted.
- Previous summary-present real output: matching threshold1e−5 and residual
 5.528534e−7.
- Missing/conflicting threshold rejection and above-threshold/negative-norm
  rejection. Modified negative fixtures existed only in parser memory.
- All other original checks and cell fields preserved; seven untouched cells
  remain explicitly not completed.

The raw scout output SHA is
`41fdcc7e535abdb8f6d07a599cb8884c9ce00a8970ca8e0bc16a3c7a34551834`;
its native execution receipt SHA is
`3ac872db5f874713838fb711dcecd469a7231b5b781c5b8e6bf9a9cff3ea4ac7`.
The measured native runtime was2h39m32.091s at12MPI. Reported maximum module
memory was732.6MB; this is not aggregateRSS. Root retains the actual allocation
and accounting receipts. This measured scout replaces any optimistic runtime
extrapolation from smaller historical contexts.

## Corrected external execution interface

`EXECUTE_V2.py` wraps the original frozen `preflight` and
`run_orca_task_manifest.run_manifest`. It requires the successful recovery
receipt, its tested corrected-collector hash, original protocol and all eight
scientific template/XYZ/state identities. It retains the per-workspace lock,
12MPI, allocation-derived memory, immutable outputs and explicit selected IDs.
The only changed collection behavior is the threshold-source repair above.
Each invocation writes fresh external preflight/collection/execution receipts;
the original v1 source and snapshots are not edited.

Invocation through the new12-CPU wrapper is:

```text
RUN_NATIVE.sbatch ABS_DRIVER ABS_MANIFEST MANIFEST_SHA
  ABS_RECOVERY_EXECUTION RECOVERY_EXECUTION_SHA DRIVER_SHA
  ABS_FRESH_OUTPUT TASK_ID [TASK_ID ...]
```

`ABS_RECOVERY_EXECUTION` is the recovery workspace's `execution/EXECUTION.json`,
SHA `04ef1aeb2bd09f552eda704416b0de52841322744038f1ee1bc451b851dc2577`.
Driver SHA is
`abaf1466089714dfa3fa6177155585116896bc221bd30a8134b6227a82ea2a36`;
tested `COLLECT_V2.py` SHA is
`24182c847f9fb76359a93b784eb1731a236b7306922fdfbcba54a6700509f412`.
Root owns the actual frozen dispatch, resource aggregate, launches and wakes.
The existing scout is reused, and seven remaining tasks are split across the
four declared workspaces. No molecule has been submitted by this agent.
