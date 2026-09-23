# Explicit operations

Run from the repository root with the existing driver:

```bash
CPCM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$CPCM_PY" scripts/solvent_cpcm_pilot.py validate --manifest workspaces/solvent_cpcm_20260922/pilot_v1/manifest.json
"$CPCM_PY" scripts/solvent_cpcm_matched.py validate --manifest workspaces/solvent_cpcm_20260922/matched_vacuum_v1/manifest.json
"$CPCM_PY" -m unittest discover -s tests -p test_solvent_cpcm_pilot.py -v
```

Jobs1209970/1209980 are terminal and the final collection exists. Do not resubmit
or rerun the completed observer. Submission receipts, allocation,
preflight and immutable implementation copies are alongside each manifest. Do
not resubmit completed/partial tasks in place. `execute` uses the existing locked
ORCA executor and explicitly rejects partial-attempt overwrites.

The observer wrote the final matched result. The following is a recovery command
only if that output were absent; it must not overwrite the existing result:

```bash
"$CPCM_PY" workspaces/solvent_cpcm_20260922/matched_vacuum_v1/implementation/solvent_cpcm_matched.py collect \
  --manifest workspaces/solvent_cpcm_20260922/matched_vacuum_v1/manifest.json \
  --output workspaces/solvent_cpcm_20260922/matched_vacuum_v1/final_collection.json
```

The paired result exposes MACE vacuum, native and ordinary GFN2 vacuum, ordinary
CPCM, matched transfer and separate model contrasts. Missing endpoints and absent
calibration remain null. All eight original CPCM attempts are reused from actual
receipts only; neither class labels nor convenient scores select a solver.
