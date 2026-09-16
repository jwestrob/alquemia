# Assessor implementation and historical preliminary checks

**Final update:** all 25 outputs collected; all 40 tests passed with no skips.
The physical gate failed. See [REPORT.md](REPORT.md) and [RESULT.json](RESULT.json).
The preliminary records below remain historical.

2026-09-16. `scripts/global_electrostatic_assess.py` collects the frozen
25-task campaign, assesses its physical gates, and writes a short report.
It launches no calculations or conditional accuracy trial.

Use the full `surfaces_v1/campaign_manifest.json` for assessment even when
execution is divided into disjoint scheduling slices. Individual solver
receipts remain tied to their own immutable TABI manifests.

## Technical validation records

Under `workspaces/global_electrostatic_20260916/`:

- `assessment_preexecution_collection_v1.json`: rejected collector attempt.
  Exact charge-vector/provenance comparison passed, but whole-dictionary
  equality falsely rejected last-bit Python summation differences. The final
  check keeps exact atomic charges/provenance and permits at most 1e−12 e in
  the redundant `sum_e` scalar; it does not normalize any charge.
- `assessment_preexecution_collection_v2.json`: rejected collector attempt.
  Two rotated inputs had one 9.9995e−11 Å reconstruction difference from
  platform rounding. Final reconstruction requires exact atom order,
  charges, and radii, with at most one unit in the tenth coordinate decimal
  (1.01e−10 Å allowing binary arithmetic). Actual input hashes and executed
  mesh equality remain strict. No prepared coordinate was modified.
- `assessment_preexecution_collection_v3.json`,
  `assessment_preexecution_v3.json`, and `.md`: source/input validation passes;
  0 computed, 0 failed, 25 unavailable solver tasks at collection time.
  The physical gate is **unavailable**, not a scientific failure.

The failures in the first two files are collector-development diagnostics,
not failed quantum/solver executions. They must not enter the execution-cost
or failed-calculation denominator.

## Tests

```bash
python -m unittest discover -s tests -p 'test_global_electrostatic_assess.py' -v
```

Eight tests ran in 59.951 s: seven passed; the completed 25-task solver-component
integration test was explicitly skipped because its actual outputs were not
yet available. Tests cover real campaign provenance, paired/gate bookkeeping,
explicitly corrupted copies of actual inputs/collections, and archived direct
Coulomb arithmetic. Archived CPCM charges are used for arithmetic only.

Final real-output validation remains pending. No synthetic scientific result
stands in for the unrun integration test.

## Collection commands when outputs are ready

The commands below read the same frozen campaign and create new files without
overwriting preliminary records:

```bash
python scripts/global_electrostatic_assess.py collect --campaign workspaces/global_electrostatic_20260916/surfaces_v1/campaign_manifest.json --output workspaces/global_electrostatic_20260916/assessment_complete_collection_v1.json
python scripts/global_electrostatic_assess.py assess --collection workspaces/global_electrostatic_20260916/assessment_complete_collection_v1.json --output workspaces/global_electrostatic_20260916/assessment_complete_v1.json
python scripts/global_electrostatic_assess.py report --assessment workspaces/global_electrostatic_20260916/assessment_complete_v1.json --output workspaces/global_electrostatic_20260916/assessment_complete_v1.md
```

File naming does not confer completion: the assessor retains failed, partial,
or missing tasks and cannot pass the gate without every required check.
