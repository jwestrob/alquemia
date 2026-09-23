# Explicit operations

Run from the repository root with the existing driver:

```bash
CPCM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$CPCM_PY" scripts/solvent_cpcm_pilot.py validate --manifest workspaces/solvent_cpcm_20260922/pilot_v1/manifest.json
"$CPCM_PY" scripts/solvent_cpcm_matched.py validate --manifest workspaces/solvent_cpcm_20260922/matched_vacuum_v1/manifest.json
"$CPCM_PY" -m unittest discover -s tests -p test_solvent_cpcm_pilot.py -v
```

Jobs1209970/1209980 are already submitted. Submission receipts, allocation,
preflight and immutable implementation copies are alongside each manifest. Do
not resubmit completed/partial tasks in place. `execute` uses the existing locked
ORCA executor and explicitly rejects partial-attempt overwrites.

Once both jobs terminate, the installed observer writes the final matched result.
Only if that observer has failed and the final result is absent:

```bash
"$CPCM_PY" workspaces/solvent_cpcm_20260922/matched_vacuum_v1/implementation/solvent_cpcm_matched.py collect \
  --manifest workspaces/solvent_cpcm_20260922/matched_vacuum_v1/manifest.json \
  --output workspaces/solvent_cpcm_20260922/matched_vacuum_v1/final_collection.json
```

The paired result exposes MACE vacuum, native and ordinary GFN2 vacuum, ordinary
CPCM, matched transfer and separate model contrasts. Missing endpoints and absent
calibration remain null. All eight original CPCM attempts are reused from actual
receipts only; neither class labels nor convenient scores select a solver.
