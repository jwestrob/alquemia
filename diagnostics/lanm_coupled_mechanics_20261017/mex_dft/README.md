# Execution handoff

All commands below run inside root-owned Slurm jobs. No molecular or numerical stage is permitted on the login node. `QUALIFY.py` does not submit jobs.

The frozen `CONFIG.json` is sufficient for the zero-energy preparation:

```bash
sbatch RUN_PREPARE.sbatch ABS_QUALIFY_PY ABS_CONFIG FRESH_OUTPUT prepare
```

Root has staged that preparation under `workspaces/lanm_coupled_mechanics_20261017/mex_dft_preparation_v1` as job12976. Inspect its actual result; do not duplicate it. The runtime copies in `native/` are byte-identical to the admitted native runner and allocation-aware renderer. Scientific imports use the pinned repository helpers, avoiding the historical runtime directory's unrelated scientific modules.

After successful preparation and input review, derive a run config by adding only `preparation_result: {path, sha256}` to the same frozen configuration. Run each scout endpoint in a fresh allocation/output:

```bash
sbatch RUN_SINGLEPOINT.sbatch ABS_QUALIFY_PY ABS_RUN_CONFIG FRESH_PLUS_OUTPUT MM_O97_La_plus
sbatch RUN_SINGLEPOINT.sbatch ABS_QUALIFY_PY ABS_RUN_CONFIG FRESH_MINUS_OUTPUT MM_O97_La_minus
```

Each execution `MANIFEST.json` names exactly one actual endpoint. Root attaches the existing pinned `ARM_HEALTH.py` and completion watcher. Each endpoint calls the existing MM and fixed-charge ALPB primitives once before its native DFT singlepoint. No MACE, GFN2 or new central gradient is called.

A saved scout collector config adds `collection_scope: "scout_pair"` and `endpoint_results`, mapping the two exact task IDs above to their actual result pins. Use `RUN_PREPARE.sbatch ... collect` on one CPU. Actual pair admission, its frozen derivative gate, and root's measured cost review precede the remaining six. Their run config additionally pins `scout_collection`.

The six remaining task IDs are `MM_O97_Dy_plus`, `MM_O97_Dy_minus`, `boundary_CA83_La_plus`, `boundary_CA83_La_minus`, `boundary_CA83_Dy_plus`, and `boundary_CA83_Dy_minus`. The final saved collector uses `collection_scope: "all_eight"` and all eight endpoint pins, retaining the original pair. Missing/failed states stay unavailable; no retries are encoded.

The source passed AST parsing and both wrappers passed shell syntax checks only before submission. The code and config are frozen for actual preparation; any subsequent correction needs a versioned recovery. The derivative result is not an electronic-state or physical-accuracy qualification, and the single fixed step does not establish truncation convergence.
