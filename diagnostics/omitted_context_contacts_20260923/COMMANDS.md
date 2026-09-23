# Reproduce the fixed geometry-only audit

All paths are explicit; this command calls no molecular executable, reconstructs
no hydrogen and changes no source/candidate geometry. Use a new result filename.
The frozen input contains all55 newly scored threefold source/context pairs and
the same sources' actual tenfold-precision counterparts.

```bash
CONTACT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$CONTACT_PY" scripts/omitted_context_contacts.py \
  --inputs workspaces/omitted_context_contacts_20260923/INPUTS.json \
  --graph-policy diagnostics/omitted_context_contacts_20260923/GRAPH_POLICY.md \
  --output workspaces/omitted_context_contacts_20260923/RESULT_replay.json
```

The saved final audit is RESULT_v1.json. Complete per-source/context/geometry
rows are `workspaces/omitted_context_contacts_20260923/export_v1/CONTACTS.tsv`;
full atom-pair lists and source/mapping pins are in RESULT_v1.json. To regenerate
only the compact table/statistics from those saved results:

```bash
"$CONTACT_PY" diagnostics/omitted_context_contacts_20260923/export_contacts.py \
  --result workspaces/omitted_context_contacts_20260923/RESULT_v1.json \
  --output workspaces/omitted_context_contacts_20260923/export_replay
"$CONTACT_PY" -m unittest discover -s tests -p test_omitted_context_contacts.py -v
```

The2.0Å flags and3.5Å summaries are fixed diagnostics. They do not reject candidates,
modify energies, supply an omitted-environment correction or alter classifications.
