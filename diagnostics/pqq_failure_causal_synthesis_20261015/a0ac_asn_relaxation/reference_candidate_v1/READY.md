# Conditional adapter readiness — unbound and unrun

Source preparation only. No workspace/manifest, root decision, numerical preflight or native endpoint has been created by this task. The [plan](PLAN.md) fixes four selected old/new Ca4 energies; the ongoing uniform transfer remains the primary experiment. Nothing here changes production or makes DFT a new transfer qualification requirement.

## Future metadata binding

Only after root chooses to proceed, from the repository root:

```bash
a0ac_reference_source="$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/reference_candidate_v1"
a0ac_reference_work="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_selected_reference_v1"
/usr/bin/python3 "$a0ac_reference_source/PREPARE.py" --output "$a0ac_reference_work"
```

This optional command hashes/copies the pinned saved sources and unchanged native helpers; it does not import scientific modules or compute geometry/energies. It refuses an existing output directory and prints the new manifest path/SHA. Set `a0ac_reference_sha` to that actual emitted SHA, not an anticipated value. Four isolated `tasks/sample4_METAL_CANDIDATE_vacuum/` directories receive exact byte-copied `core.xyz`, fixed native `endpoint.inp` and future `endpoint.out` paths. Scientific state: 204 atoms, 14 caps, Ca−4/La−3 singlets, original H and no solvent.

## Future Slurm checks and four endpoint commands

Root owns every submission and completion wake. All three actions use one standard node, 12 MPI slots, no memory/time/priority overrides. No action below has been run:

```bash
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" preflight
```

Actual preflight replays source/paired/cap/selection joins, electronic counts, runtime qualification, inherited real parser fixtures, eight ExaA collections and six MMOL collections. No ORCA endpoint runs in `preflight`.

For any `execute`, supply the actual root-authored decision path/SHA with:

- `decision: "execute_four_selected_Ca4_native_endpoints"`;
- `manifest: {path, sha256}` for this actual new manifest;
- `transfer_collection` and `transfer_execution`: exact terminal full uniform-transfer pins, retaining all 257 contexts/260 memberships and failed/unavailable cells;
- `interpretation`: root's scientific reason to run this fixed reference after reading that result.

The adapter does not create that decision or infer it from a favorable label/score. Call each old-selected scout explicitly:

```bash
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" execute \
  sample4_Ca_adaptive_Ca_vacuum \
  --decision "$a0ac_decision_path" --decision-sha "$a0ac_decision_sha"
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" execute \
  sample4_La_adaptive_La_vacuum \
  --decision "$a0ac_decision_path" --decision-sha "$a0ac_decision_sha"
```

New-selected endpoints additionally require a root scout review containing `decision: "continue_two_new_selected"`, the same manifest pin, both actual admitted one-cell `scout_collections` pins, and `interpretation` of their state/convergence/runtime. Then the two explicit commands are:

```bash
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" execute \
  sample4_Ca_relaxed_flip_Ca_vacuum \
  --decision "$a0ac_decision_path" --decision-sha "$a0ac_decision_sha" \
  --scout-review "$a0ac_scout_review_path" --scout-review-sha "$a0ac_scout_review_sha"
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" execute \
  sample4_La_relaxed_flip_La_vacuum \
  --decision "$a0ac_decision_path" --decision-sha "$a0ac_decision_sha" \
  --scout-review "$a0ac_scout_review_path" --scout-review-sha "$a0ac_scout_review_sha"
```

No automatic retries or iteration changes occur. A negative/missing decision or failed scout blocks the dependent molecular action. These are internal scientific gates, not fresh user-permission requests.

## Future collection and limits

After terminal outcomes, collect all four declared cells, including missing/failed cells:

```bash
sbatch "$a0ac_reference_work/implementation/RUN.sbatch" \
  "$a0ac_reference_work/manifest.json" "$a0ac_reference_sha" collect
```

Attempts are immutable `attempts/JOBID_ACTION_TASK/` directories (all-four collection ends in `all_four`). Each has `PREFLIGHT.json`, applicable root/scout review receipts, `COLLECTION.json` and `EXECUTION.json`. Native task directories retain raw ORCA outputs, execution/resource receipts and GBWs. Collection reuses the count-aware/native corrected collector, preserving actual solver criteria, density residuals, warnings and failures. Per-metal work is new-selected minus old-selected; paired change is Ca work minus La work. Missing cells remain null. The 1 kcal/mol interpretation scale is not an acceptance test/error bar, and a positive paired change alone is not a validation claim. No DFT selection, solvent or class decision is calculated.

The unchanged parent is `mmol_reference_v1/manifest.json`, SHA `10c23bbfcdfacb6a5ff3ea5dbaa6b8661f6ef8cc79681289b4d83b1252b58dc8`. Binding copies its exact executor/render/parser/collector snapshots and pins the actual source/method/runtime fixtures; no mutable repository defaults supply native execution. The actual 12-MPI cost precedents and current input/cache distinctions are in [PLAN.md](PLAN.md). Readiness is not evidence of runtime or DFT success.
