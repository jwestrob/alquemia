# Saved-summary readiness

Prepared only. Root binds after job5158 has written terminal EXECUTION.json.
No molecular calls or numerical tests were run during preparation. Python AST
and shell syntax checks passed; actual archived schema checks run in Slurm.

Metadata-only binding:

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_solvent/report/SUMMARY.py freeze \
  --source-config workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_v2/CONFIG.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_report_v1
```

The command prints the new CONFIG.json pin. Submit that exact bound file:

```bash
CANONICAL_REPORT_CONFIG="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_report_v1/CONFIG.json"
CANONICAL_REPORT_SHA=$(sha256sum "$CANONICAL_REPORT_CONFIG" | cut -d ' ' -f 1)
sbatch workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_report_v1/implementation/RUN.sbatch \
  "$CANONICAL_REPORT_CONFIG" "$CANONICAL_REPORT_SHA"
```

One shared CPU, default scheduler memory/time, no molecular evaluations. Root
owns submission and completion wake. Outputs are RESULT.json, CANONICAL.tsv,
A0AC.tsv, A0AC_RANGES.tsv and EXECUTION.json. A partial/failed source is reported
without fabricating its missing rows or calibration. The real5135 missing-result
fixture and actual5030 score/range regression execute in the same saved-only job.
