# Operational persistence diagnostic

The single attempted nonblocking py-spy dump was denied permission (exit1). No stack was captured. Sampling stopped immediately; there was no retry, elevated access, pausing attach or workaround. The empty stdout, exact stderr and both own-UID process snapshots are preserved under `samples_v1`. No scientific result, candidate, Python local or mutable fit artifact was read.

Before/after observations at2026-10-22 22:10:07 UTC show the same controllerPID3944417/start tick61001214,13threads, one `analyze` descendant then none. Controller CPU increased0.29seconds across approximately0.285seconds between the snapshot timestamps. Cumulative process write_bytes were25,192,235,008 and wchar25,073,717,012; neither advanced during this short sample. These are process-wide counters, not measured ledger-only bytes or proof of physical disk throughput. RSS was198,749pages. The snapshots are non-atomic and instantaneous; they cannot establish sustained native utilization.

Sampling cost was0.549seconds wall,0.519seconds sampler CPU and0.020seconds child CPU. Only one py-spy invocation was made. The short observation supports controller activity and substantial cumulative writes, but **does not locate the executing frame or quantify persistence's runtime share**.

## Source-backed bottleneck hypothesis

Frozen ADAPTER.py `save` indents/serializes the whole object to a temporary file and replaces the target. Native attempt entry and return each deep-copy the complete growing attempt history into the ledger and rewrite the whole ledger while holding the global lock. Each completed cell also rewrites that ledger under the same lock. Return handling additionally scans all historical attempts to count observed returncodes. With P attempts, these history operations have quadratic aggregate growth; the lock can serialize workers that are ready to start or report a native program.

The parent evaluator's snapshot deep-copies requests and hashes the full ledger; DRIVER.persist includes that snapshot in a growing result at vector/step boundaries. Those boundary operations are separate from the more frequent per-attempt/per-cell hot paths. Parsing, assembly and solver work are other possible controller costs; no sampled stack distinguishes them. No speedup is asserted from the denied profile.

The current campaign and every pinned source remain unchanged. Performance evidence is not scientific evaluation or a reason to alter its candidate selection.

## Minimal future persistence design

Use one append-only journal owned by the controller, with a small synchronous record under the existing global lock at each event. No additional background worker or change to native scheduling/locking is needed. Events bind configuration/source/release, monotonically increasing sequence number, previous-record checksum and current checksum. Encode one bounded event, write all its bytes, and fsync before acknowledging durability. An I/O error sets the original stop condition; it must not permit the next native launch.

Record types: initialization, request creation, native-attempt intent, native return/failure, cell-request status, request completion/failure and global stop. Keep immutable per-cell artifacts and exact native receipts/input pins. Publish those with temporary-file fsync, atomic rename and directory fsync before journaling a reference to them. Do not copy all past attempts into every event. Maintain counters incrementally, with attempt IDs and cell/program identities unique and checked.

The exact attempt intent must be durable **before** calling the original native helper. A crash after intent and before launch is intrinsically ambiguous without a durable launcher handshake; charge that unresolved intent against the cap and prohibit automatic retry. Do not label it a proven completed program. Native return evidence must be written durably before its completion event. Completed cache reuse requires the same full scientific identity, admitted immutable cell artifact and receipt pins; an intent or partial output never authorizes reuse.

Replay checks a contiguous sequence and checksum chain, rejects duplicate attempt IDs or terminal transitions, enforces request/cell identities and original caps, and preserves stop/failure history. A final incomplete journal record is a crash suffix, retained and explicitly reported; corruption in the committed prefix is an error. Do not silently trim or append across an unresolved suffix. There is no automatic campaign resumption or failed-call retry.

Materialize the original complete ledger at vector boundaries, failure and finalization. Compare its semantic contents with journal replay, then atomically publish it. Preserve the legacy terminal schema and all accounting fields so existing saved reviewers can inspect the same record. During hot paths, the append-only journal is the authoritative durable state; a full ledger is explicitly a checkpoint at its recorded journal sequence. A future run needs a separately reviewed persistence version/configuration; this does not substitute into the live campaign.

Prospective controls are purely artificial: exact final ledger reconstruction, crash suffix and unresolved-intent handling, duplicate rejection, stopped/capped attempts, partial-cell failure and interleaved writer ordering. No optimizer, scientific evaluation or native program is needed. The next source prototype will retain the frozen scientific render/parse/assembly, invocation order, locks, budgets and selection code.
