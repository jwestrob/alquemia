# Conditional readiness — relaxed Asn solvent

Code is prepared; no solvent configuration is bound and no job is submitted here. Actual complete remaining-stage results and the frozen sign-blind gate are required. Root reviews them and the contact results before permitting this stage. Dependency pins are recorded in `DEPENDENCIES.json`; bind pins the final implementation and actual input collections without copying existing source trees.

## 1. Metadata binding, only after the actual gate passes

From the repository root:

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/solvent/SOLVENT.py bind \
  --score-config workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_remaining_v1/CONFIG.json \
  --score-collection workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_v1
```

Binding rejects a missing/partial collection, noncredible proposal, changed population, failed full-source admission or false recorded gate. It does not evaluate coordinates or energies. The immutable config pins all40 actual final cells, the80 prior cells through their completed collection/config, both actual MACE stage receipts, source maps/state and unchanged solvent dependencies.

## 2. Root submits a preflight-only job first

```bash
asn_relaxed_solvent_dir="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_v1"
asn_relaxed_solvent_sha="$(sha256sum "$asn_relaxed_solvent_dir/CONFIG.json" | cut -d ' ' -f 1)"
sbatch --parsable \
  --output="$asn_relaxed_solvent_dir/preflight_%j.out" \
  --error="$asn_relaxed_solvent_dir/preflight_%j.err" \
  "$asn_relaxed_solvent_dir/RUN.sbatch" "$asn_relaxed_solvent_dir/CONFIG.json" \
  "$asn_relaxed_solvent_sha" preflight
```

This performs saved-data validation only: replays80 accepted common4 cells, both stage/gate joins, all20 final physical/source/cap mappings and paired40 frames. It writes `preflight_JOBID/PREFLIGHT.json` and `EXECUTION.json`; no SCC, native replay, MACE or optimization starts. Root records the actual command/job and completion wake, then inspects these receipts. Existing outputs are exclusive-created.

## 3. Separate root-reviewed molecular execution

Only after the actual preflight passes and root accepts the physical/contact evidence, use the same immutable config/wrapper with final action `execute`:

```bash
sbatch --parsable \
  --output="$asn_relaxed_solvent_dir/slurm_%j.out" \
  --error="$asn_relaxed_solvent_dir/slurm_%j.err" \
  "$asn_relaxed_solvent_dir/RUN.sbatch" "$asn_relaxed_solvent_dir/CONFIG.json" \
  "$asn_relaxed_solvent_sha" execute
```

Execution repeats input validation, then calls the existing seven-worker component-v2 dispatcher/native replay. No automatic submission, retries or resumed partially executed points. Final collection preserves all ten rows/common3/common4/common6, including failed required cells; the actual execution receipt records backend attempts and costs. There are no class decisions, calibration changes or production edits.

Validation so far: static Python AST and Bash syntax only. Gauss independently cleared the actual source joins, gate, complete-pool handling and unchanged native recipe at `SOLVENT.py` SHA256 `9c121df6020041e278e191c88effea036ffaa861b705572c1026d987690273e2`. Wrapper SHA256: `1837ec9802529e95ed92d28e16ef6ac37de03bcecee83741640e5b717dc30e6e`. The complete-data preflight has not run.
