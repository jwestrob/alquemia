# Saved report readiness

Single stdlib CLI; no molecular calls or scoring imports. Unrun until the actual transfer solvent collection is available. Root submits from the repository root on one Slurm CPU:

```bash
sbatch diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/report/RUN.sbatch \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/report/REPORT.py" \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_solvent_v1/execution/FINAL_COLLECTION.json" \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_solvent_v1/execution/EXECUTION.json" \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_report_v1/execution"
```

The supplied `RUN.sbatch` accepts the explicit reporter path followed by collection, execution and output paths; it never resolves code relative to Slurm's spooled script. Its single Python CLI is `REPORT.py --collection PATH --execution PATH --output DIR`. Use the actual final solver workspace if its version changes. Output must be new. Outputs: REPORT.md, RESULT.json, INPUTS.json, EXECUTION.json and COHORTS, MEMBERS, TRANSITIONS, REGRESSIONS, CANDIDATES, SOURCE_SPREADS, ORDERING_DIAGNOSTIC, TRIPLES, COSTS and SUPERSEDED_PRIMARY TSVs.

This is an exporter over recorded fields. Scientific criteria remain owned by SOLVENT.py. Partial results retain all rows; malformed/missing input fails with FAILURE.md/EXECUTION.json. No transfer fixture is synthesized and no work has been submitted.
