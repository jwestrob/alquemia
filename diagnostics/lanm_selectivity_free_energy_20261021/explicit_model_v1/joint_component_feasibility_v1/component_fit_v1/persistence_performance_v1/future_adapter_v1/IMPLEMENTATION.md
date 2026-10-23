# Future persistence prototype — unexecuted

This is an isolated source prototype. The live component campaign, original ADAPTER/DRIVER, science, native kernels and current selection remain unchanged. No runnable campaign configuration or release is created here. It is not a continuation mechanism.

## Contained delta

The copied DRIVER.py is byte-identical to the frozen driver. The future ADAPTER.py changes only persistence, incremental accounting of already-observed returncodes, and the new journal-source/policy binding. `_cell`, `_gate`, `authorize`, `required_programs` and the original native slot/global-lock structure are AST-identical. Rendering, parsing, assembly and native invocation still use the original parent evaluator. Scientific selection, budgets, source/build admission, invocation order and stop rules remain intact.

`JOURNAL.py` writes a compact event while holding the same controller RLock. Attempt intent is fsynced before the unchanged native helper call. Terminal events preserve the exact new native receipt and input-pin evidence. Cell-request events reference the durable per-cell RESULT or admitted immutable-import/baseline authority. There is no growing history copy or scan at attempt entry/return, and no full-ledger write at ordinary cell completion. Event cost is constant in historical attempt count; per-cell scientific artifact sizes remain those of the fixed panel.

The in-memory ledger remains compatible with the old schema. Native-attempt history now references the same mutable attempt objects under the existing lock rather than copying all prior objects twice per program. Returncode count updates once from the changed attempt's flag. Journal state is separately copied/reduced, so later mutation cannot rewrite old events.

Initialization, vector completion/failure, stop handling and snapshots materialize the complete legacy ledger atomically. Snapshot/finalization compares that ledger with incremental journal replay. A sidecar records its committed journal sequence/hash, making an older compatibility ledger distinguishable from the active journal. DRIVER calls snapshots at vector/step boundaries rather than per native/cell event. This prototype does not remove the lower-frequency growing-result/snapshot work.

Atomic JSON writes use file fsync, replace and parent-directory fsync. Journal writes handle short writes, fsync each event and chain sequence/checksum records. This requests filesystem durability through POSIX interfaces; it does not prove power-loss behavior of the shared filesystem or supply distributed consensus. No durability/performance claim is based on the denied py-spy dump.

## Crash and failure semantics

A durable intent consumes the original attempt cap before launch. If a process crashes before a return receipt, the event remains unresolved and charged; a launch cannot be inferred either way solely from that intent. Such cells cannot be reused or automatically rerun. Replay accepts a valid committed prefix plus an explicitly reported incomplete final suffix without editing either; committed-prefix corruption, duplicate IDs/returns, cap violations and events after stop are rejected as appropriate. Replay is read-only and does not resume execution.

A returned in-flight attempt may still publish its receipt after global stop. Queued slots remain untouched. A write/fsync failure prevents the helper from launching when intent durability has not completed; the journal refuses subsequent events after an I/O failure. Existing files and the committed prefix remain available for diagnosis. The prototype cannot guarantee a complete terminal JSON if the underlying filesystem itself refuses writes; it never conceals that as a successful run.

The original no-failed-cell-retry and fresh-cache guards remain. Imported scientific records still enter only through the unchanged immutable lineage joins. Completed reuse retains exact source/input/receipt identities. An exact duplicate vector still needs the original explicit probe/trial authorization; the persistence layer creates no selection rights.

## Configuration and review before any future use

A future configuration must bind the copied adapter/driver paths, plus:

```
persistence_policy = append_only_component_ledger_v1
journal_implementation = exact JOURNAL.py pin
```

These fields enter the configuration fingerprint, so an old release/source-review/preparation cannot authorize this version. The entire future configuration, immutable scientific lineage, source review and root release must be rebound prospectively. No live configuration is changed and no automatic migration/import of its mutable ledger is provided. Terminal reviewers should check the journal sidecar/replay in addition to the compatible legacy ledger before admitting this future implementation.

## Artificial controls and limitations

Thirty-two initial controls verify exact replay/materialization, unresolved intents, preserved partial suffixes, corruption/duplicate/cap rejection, stopped/untouched cells, native intent before a pure mock helper, immutable source joins and unchanged driver/science/lock structure. Eight focused controls exercise the actual future `_invoke` failure path, terminal receipt after stop, queued-slot rejection and a mocked journal-write failure that prevents any helper call. No subprocess/native call, fit proposal or scientific evaluation occurs.

All controls use one OS/numerical thread. Three exhaustive valid event interleavings with sequentially assigned attempt IDs cover two logical writers, including out-of-order returns. This verifies reducer order semantics under the existing lock; it is explicitly not a real multi-thread scheduler stress test. Filesystem crash/power-loss testing and independent source review remain future prerequisites, not performed qualifications. The first generator failed its exact-line replacement guard before creating any adapter; that source/failure is preserved, and the version2 generator only fixes the line anchor.

No speedup has been measured. The current observations justify testing a smaller persistence path, not a claimed percentage of runtime saved or a scientific change.
