# Fixed-grid SCF diagnostic implementation

`scripts/metal_environment_scf_check.py` prepares exactly four real-source
Ca/La A/rigid endpoints. It preserves exact refined geometry and field bytes,
native electronic recipe and fixed quadrature. Root's SCF_CHECK_PLAN.md is the
agreement; root owns submission. No molecular work ran during this subtask.

Optional named parser profile `tight_all_active_criteria_forced_125_v1` adds
ConvCheckMode 0, ConvForced 1 and explicit unchanged MaxIter 125/TightSCF
thresholds. Old parser/default and grid-only inputs remain unchanged. Scientific
templates are exact and versioned; no arbitrary tolerance overrides are accepted.

Collection requires actual all-criteria/forced header evidence and independently
checks final energy, RMS/maximum density, orbital gradient and rotation residuals
against the declared thresholds and printed tolerances. DIIS error is retained;
when the actual solver switches to SOSCF, its historical DIIS value is not
misreported as the terminal active error. Otherwise its DIIS residual is checked.
Output residual observations survive convergence/receipt/parser failures. Rigid
energy and core/MM gradient tolerances are unchanged. Passing these four cells
would not qualify new finite differences or an A/B response: B is absent.

CLI:

```
python scripts/metal_environment_scf_check.py prepare \
  --source-collection workspaces/metal_environment_response_20260926/grid_check_v1/FINAL_COLLECTION.json \
  --agreement diagnostics/metal_environment_response_20260926/SCF_CHECK_PLAN.md \
  --output workspaces/metal_environment_response_20260926/scf_check_v1 \
  --workers 4 --mpi-ranks 86
```

`dry-run`, `execute`, `collect` use `--manifest`; collection optionally writes
`--output`. The source manifest is resolved and verified from the actual source
collection. Prepared snapshots retain existing execution receipts and memory
policy; no scientific result is inherited into a new endpoint.

Three real-artifact tests pass, zero skips, 3.18 seconds:

```
python -m unittest discover -s tests -p test_metal_environment_scf_check.py -v
```

Tests inspect actual residuals, retain old refined-output parsing, reject those
outputs as strict-SCF evidence, and prepare/validate a temporary four-cell
manifest with exact source bytes. Removed-cell corruption fails. There is no
fabricated successful strict-SCF output; actual new parsing remains unrun until
calculations complete.
