# Frozen expanded-candidate commands

Run from the repository root. Root already launched canonical v3; inspect its
launch and execution receipts before invoking any execution command again.

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py dry-run --directory workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/TEST_ARTIFACTS.py
```

Canonical execution (already launched by root, 14 total numerical threads):

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py execute --directory workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3 --stage canonical --workers 14
```

Only after root reads the actual `CANONICAL_COLLECTION.json` and confirms valid
operational calibration, run the separately gated transfer stage:

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py execute --directory workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3 --stage transfer --workers 14
```

Read-only partial collection, using a new unused output name:

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py collect --directory workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3 --output-name MANUAL_PARTIAL_001.json
```

Do not edit the frozen driver/plan during execution. A changed scientific rule
requires a new protocol/preparation. Completion notification is root-owned via
`WATCH_STAGE.py`; exact PID/starttick and queue command are in `CANONICAL_LAUNCH.json`.
