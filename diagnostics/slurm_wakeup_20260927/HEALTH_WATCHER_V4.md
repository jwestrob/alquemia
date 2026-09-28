# v4: coalesce successful partial completion notifications

Future submissions only. Current processes and frozenv1/v2/v3 receipts are untouched. `native_scf_health_v4` changes only success-progress event selection. Existing pathology thresholds, failures and final scheduler-terminal events remain immediate. No automatic cancellation or runtime cap is introduced.

Successful endpoints are tracked in every saved snapshot but do not individually enqueue messages. If the job remainsRUNNING with both successful and unfinished endpoints for120seconds, emit one aggregate `successful_partial_persistent` alert perjob. Counts may increase during that interval without resetting it; completion or leavingRUNNING clears pending partial state. This is allocation-review notification, not inference of CPU utilization. A final terminalevent bypasses the timer and subsumes any unqueued successful progress.

A successful endpoint requires normaltermination, SCFconvergence, no explicitfailure and no nonzero returncode. Any terminal endpoint not meeting those conditions still receives an immediate reviewevent. Numerical pathology alerts always bypass the debounce, including while successful peers finish.

## Actual replay

Saved job1220300 receipts report first endpoint completion10:39:07UTC and last10:40:13UTC (about67seconds). The old watcher enqueued four messages from10:39:08through10:40:41, three carrying partial-success progress and one final. Replaying exactly those saved observation timestamps throughv4 emits only the finalterminalevent. No messages were actually queued by this test. Administrative completion states were reconstructed from genuine execution receipts; no molecular output or numerical energy was fabricated.

Additional tests use the retained earlier two-La-complete/two-Dy-unfinished scout to verify the120-second timer and one-event deduplication, and actual failedPModel output to verify failure/pathology bypass. All11Slurmwatcher tests pass. These replay tests prove eventselection, not a new end-to-end service test; the underlying supportedCodexqueue wake path is unchanged.

Use the same command flags with a freshv4receipt on the next submission. Do not restart activev3watchers, edit their receipt identities, or remove alreadyqueued messages to apply this change. The manifest/policy identity check intentionally rejects attempts to reuse an oldreceipt withv4.

Run: `python -m unittest discover -s tests -p 'test_slurm*' -v`. Actualreceipt pins andpolicy: HEALTH_POLICY_V4_AND_TESTS.json.
