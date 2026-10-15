# Factorial execution complete

Root submitted the frozen one-CPU stage as **4813**, with completion watcher
PID3512684. It completed all8 new native off-diagonal evaluations and reused
all8 actual diagonal results; all four metal/source factorials are available.
See [REPORT.md](REPORT.md) and [RESULT.json](RESULT.json) for actual results.
Do not rerun this completed stage.

Frozen manifest:
`workspaces/pqq_failure_causal_synthesis_20261015/asn_solvent_factorial_v1/manifest.json`,
SHA256 `e3bd63c88a56fadb81687ae637cfc06518087c4c67bdc48880b4290a77c9bed4`.

Reviewed source pins:

- FACTORIAL.py: `9d670298861ae355721d851f974224962c7c6d8a411494c9334cc05a25af6c24`.
- RUN.sbatch: `6bc0ca633c5c91fe249623e596bb7ea5df7607c5d689f71fc19d55301aac050c`.
- PLAN.md: `b8056a623aef4b60af814bc59e50a2c705dabd1300a696f4a6492c3a8f397853`.

Independent final static review is in [REVIEW.md](REVIEW.md). Only AST/Bash
syntax and metadata binding ran locally. All scientific preflight, native
evaluation and factorial algebra ran in Slurm4813.

## Commands recorded for reproducibility

Metadata binding, already completed:

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_solvent_factorial/FACTORIAL.py prepare \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_solvent_factorial_v1
```

Root's submitted command, expressed from the repository root, is recorded below.
This is provenance, **not a retry instruction**:

```bash
asn_factorial_dir="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_solvent_factorial_v1"
sbatch --parsable \
  --output="$asn_factorial_dir/slurm_%j.out" \
  --error="$asn_factorial_dir/slurm_%j.err" \
  "$asn_factorial_dir/RUN.sbatch" "$asn_factorial_dir/manifest.json" \
  e3bd63c88a56fadb81687ae637cfc06518087c4c67bdc48880b4290a77c9bed4
```

The workspace `SUBMISSION.json` preserves the exact absolute argv and watcher
command. `execution/PREFLIGHT.json`, `FINAL_COLLECTION.json`, `PARTIAL_RESULTS.json`
and `EXECUTION.json` retain actual scientific admission and costs. No new diagonal,
SCC, MACE, DFT or optimizer calls occurred. All earlier branch records remain intact.
