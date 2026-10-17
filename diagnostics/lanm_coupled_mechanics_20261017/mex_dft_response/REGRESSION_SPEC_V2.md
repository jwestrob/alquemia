# Saved-fixture failure-gate regression V2

This tests the reviewed `RESPONSE_REFERENCE_V2.py` source `3acf0887…` and template
`a35bde86…`. It does not prepare molecules, admit a reference result, or execute
an engine. The failed preparation 13059 remains unchanged and cannot be used as
a passing fixture. Root runs this only after real preparation **and** preflight
have succeeded. A versioned adapter repair requires a separately reviewed test
revision; the runner does not accept another source silently.

V1 regression source/spec remain unchanged, with byte-identical copies in
`archive_v1/`. This version changes only the pinned adapter/template, test
protocol ID, and saved-replay hook signature for the added diagnosis argument.
The five controls and their acceptance rules are unchanged. The independently
reviewed adapter repair binds actual zero-call diagnosis 13065 and only admits
the six demonstrated origin/pilot wrapper aliases; the eight DFT task identities
are unchanged. Actual successful V2 preparation/preflight are still required.

## Five controls

| Control | What it proves |
| --- | --- |
| Saved preflight, no supplied qualification | The isolated harness reaches the intended successful preflight path using actual passing receipts. Preparation alone still does not authorize a native endpoint. |
| Endpoint, missing qualification | The actual qualification gate rejects `null` before a native request. |
| Endpoint, real preparation supplied as qualification | An existing, correctly hashed but scientifically incompatible receipt cannot stand in for the full eight-fixture derivative qualification. |
| Preflight, the same incompatible qualification | A late gate exception replaces the provisional `complete` status with `execution_failure`, retains `status_before_failure: complete`, and exits unsuccessfully. |
| Nonfinite approximate gradient | The actual force-report finite-array check rejects one deliberately corrupted in-memory element of a real approximate origin gradient; the outer main handler records failure. The matching reference origin and shared MM/ALPB arrays are real saved evidence. The NaN is a labeled test fault, never a molecular observation. |

Each negative case must emit `execution_failure`, a nonzero return code, the
specific expected failure, and zero requested/attempted component calls and
optimizer steps. A tripwire rejects any attempted engine or full qualification
collector entry. The original sources and all molecular artifacts stay intact.

## Scope of isolation

The runner copies exact adapter function code objects into a private namespace.
Four test-only hooks supply already admitted input context, task inventory,
common-pool replay and shared components from the successful saved fixtures.
This avoids repeating the complete preparation. The actual `main`, `endpoint`,
`preflight_gate`, `qualification_gate` and `force_report` code remains unchanged.
The nonfinite case additionally replaces only the in-memory read of its named
approximate-gradient artifact. No original module globals are modified.

These are focused negative-path regressions, not another native integration,
scientific tolerance, eight-cell qualification or measured-cost approval.
They do not establish that a future positive qualification or cost decision
passes; those remain real adapter execution gates.

## Execution

Use the interpreter pinned by the original recipe config, in a root-owned
one-CPU Slurm allocation. No memory/time override or GPU is needed. All paths
below must name actual artifacts; output must be fresh.

```text
PYTHON REGRESSION_V2.py --template CONFIG_V2.template.json \
  --preparation-result ACTUAL_PREPARATION/RESULT.json \
  --preflight-result ACTUAL_PREFLIGHT/RESULT.json \
  --output FRESH_REGRESSION
```

`RESULT.json` pins the runner, frozen adapter/template, actual fixture receipts
and every case result. It records test status and elapsed application time.
Scheduler allocation/cost evidence is attached by root separately. Source-only
AST and metadata review precede execution; no test result is claimed here.
