# Source-status clarification before model execution

The compact handoff and PLAN.md misidentified the existing 1KB0 whole-protein
block as TRO parameterization. The actual V2 preparation manifest rejects a
peptide connection across missing structure: A/511/C to A/513/N is 4.750089 A.
The earlier audit also found the 18.969263 A connection from C573 to N579.
This is the authoritative reason for keeping 1KB0 unavailable here. See
`diagnostics/mace_omol_20260917/INTACT_CANONICAL_INTEGRITY_ADDENDUM.md` and
`workspaces/mace_omol_20260917/intact_panel_prepared_v2/preparation_manifest.json`.

This corrects explanatory provenance only. The same 25 calibration structures,
two cached crystal transfers, unavailable third transfer, model, 50-call scope
and decision rule remain fixed. No new molecular evaluation had started when
this clarification was recorded. Preserve the original plan and its hash.
