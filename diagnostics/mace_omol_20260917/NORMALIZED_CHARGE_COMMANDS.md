# Normalized endpoint charge operations

Run from the repository root, with fresh output paths. Job1200970 was submitted
once; inspect its state before attempting execution. This diagnostic returns
no affinity score or solvent correction. See [plan](NORMALIZED_CHARGE_PLAN.md).

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/normalized_charge_v1/implementation/mace_omol_charges.py dry-run --manifest workspaces/mace_omol_20260917/normalized_charge_v1/manifest.json

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/normalized_charge_v1/implementation/mace_omol_charges.py report --manifest workspaces/mace_omol_20260917/normalized_charge_v1/manifest.json --output workspaces/mace_omol_20260917/normalized_charge_report_replay_v1

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_omol_charges.py -v
```

Preparation replay without cluster submission:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_charges.py prepare --quantum workspaces/mace_omol_20260917/matched_H_quantum_v1/manifest.json --plan diagnostics/mace_omol_20260917/NORMALIZED_CHARGE_PLAN.md --output workspaces/mace_omol_20260917/normalized_charge_preparation_replay_v1
```

The exact submission argv and manifest hash are in
`normalized_charge_v1/submission.json`. Wrapper `run_normalized_charges.sbatch`
uses eight one-thread utility workers on standard-shared with16GB host memory.
No application time/CPU budget. It executes the frozen module alongside its
manifest, with a nonblocking exclusive lock and distinct attempt directories.
Successful task receipts are reused on recovery; failed attempts require the
explicit `execute --retry-failed` flag inside an allocation. Preserve prior
attempts and record the retry launch. Never rerun a complete job to produce
a missing report: collection is a separate, read-only operation.

Each task copies the native GBW/densities/index, then runs CHELPG and exact
potential utilities sequentially. Utility resource receipts include process
startup. All source copies must retain their hashes after execution. Printed
charge precision is retained; there is no charge renormalization. Unavailable
outputs remain explicit, and actual integration tests skip until real results
exist. Software test success does not substitute for actual charge-quality gates.
