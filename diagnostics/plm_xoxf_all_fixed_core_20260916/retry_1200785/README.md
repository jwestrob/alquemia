# Workflow-status repair after jobs1200785–1200788

The initial container launcher validation failed before any new MSA or AF3
calculation. All174 new proteins therefore have three missing model files.
The first preparation adapter mistakenly labeled these as `unsupported`, and
the publisher consequently reported zero workflow failures. It generated no new
chemical preparations or quantum tasks. The two prior completed results were
correctly reused.

This retry adapter preserves the first run's frozen code and output files.
It changes outcome classification only:

- Missing fold records, missing/invalid output files and source integrity
  failures are `preparation_failed`, with the stage and upstream reason retained.
- A valid selected model failing the existing CN/confidence/direct-N gates, or
  lacking a supported homologous core mapping, remains `unsupported`.
- `unsupported_targets`/`unsupported_count` contain only scientific exclusions.
  Separate `failed_targets`/`failed_count` retain workflow/preparation failures.

All176 proteins, two reused results, model selection, chemical helpers,
hydrogen policy, endpoint method and decision bands are unchanged. The CLI and
ready-only `prepared_pairs.json` interface are unchanged. Use this directory's
`prepare_batch.py` for the explicit retry; its three sibling Python files must
be pinned together. No automatic retry, model substitution, new protonation or
ORCA execution was added or performed here.

Sixteen lightweight tests passed in3.314s, including the actual retained missing
fold fixture and the distinction from a genuine selected-model gate failure.
A read-only replay of all174 missing target records yields **174 workflow
failures, zero chemistry exclusions and two reused results**. Replay products:
`workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/missing_output_status_replay/report.json`.
Detailed evidence and final hashes: `validation.json`.

The separate folding agent owns container launcher portability; the parent owns
the retry job chain and watcher. This repair does not reinterpret the missing
models as evidence about coordination or metal preference. The original
residual-H-force and computational-calibration limitations remain in force.
