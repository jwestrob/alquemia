# Shutdown recovery — 7 October 2026

The remaining Hans8FNR pair did not finish. Saved Slurm accounting records job 1220342 cancelled by uid 0 at September 28 10:59:32, after 15,943 seconds on 224 allocated CPUs (3,571,232 allocated CPU-seconds). Both output files stop during TRAH, without SCF convergence or normal termination. Neither energy is usable. This is an interrupted numerical recovery, not a completed negative biological result.

Re-running the existing comparison on saved artifacts recovers four of six accepted endpoints. Hans8DQ2-minus-Mex remains +11.689857 kcal/mol with the matched finite classical terms; Hans8FNR remains unavailable. The separate Glu91 experiment supports a physical differential response but shifts the common-pool contrast by only −0.173940 kcal/mol. A source-robust LanM discriminator is still unproven.

The old completion receipt still says RUNNING because its watcher stopped before cancellation. Saved terminal accounting supersedes that stale receipt. The expected final collector and aggregate files were not produced. No chemistry was rerun, and no old output or receipt was overwritten. Compute remains unavailable according to Jacob; no submissions or watcher restarts were attempted.

Next: inspect saved orbital files for a defensible, separately recorded same-state continuation when compute returns. Unconverged orbitals can only be a numerical starting hypothesis, never an accepted endpoint or automatically admitted seed. Do not increase iteration limits or repeat PModel initialization without a specific numerical diagnosis. Preserve all four accepted endpoints. If the missing pair eventually converges, compare both Hans sources before deciding on matched structural accommodation; do not select the favorable source.

Reproduce the file-only comparison from the repository root (choose a fresh output filename):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARE.py --config diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARISON_RECOVERY_v1.json --output workspaces/metal_environment_response_20260926/recovery_20261007/comparison_repeat.json
```

RESULT.json pins the terminal evidence and regenerated comparison. PQQ production is unchanged.

## Saved orbital assessment

Both approximately116MiB GBW files survived and are pinned in ORBITAL_INVENTORY.json. The native ORCA6.1.1 exporter reports successful GBW reads for both and reads the matching2008 point charges. Its default export then reconstructs integrals; both owned utility processes were deliberately terminated before full export to avoid continuing computation on the login host during the outage. One earlier La read failed because environment.pc was absent from the working directory; that log is preserved. No SCF was run. Binary completeness, occupation/state integrity and successful MORead continuation remain unqualified. Existing accepted-source seed validation remains unchanged and correctly rejects these unconverged sources. Raw utility logs and scratch copies are under workspaces/metal_environment_response_20260926/recovery_20261007/orbital_read/.
