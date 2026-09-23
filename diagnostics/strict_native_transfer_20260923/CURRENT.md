# Strict225 checkpoint

2026-09-23: four finite CPU-only jobs1211337–1211340 are running on
node-224-2t-8gpu-1,32CPU/64GiB each. Each runs600 one-rank native scalar calls.
No MACE, DFT, geometry or preparation calculations are added.

Exact prepared scope:225 original sources,208 supported;200 fresh pools and
eight exact shared reuses.17 old preparation exclusions remain. All2496 source
cells passed actual coordinate/state/receipt/MACE audit. Four prelaunch tests
pass in3.911s, zero skips. Preparation150.387s local wall; no molecular output
from the new jobs is claimed here.

Workspace: `workspaces/strict_native_transfer_20260923/run_v1/`.
Inventory SHA6e4608a5db68b5b07b34ea7680b101a73a49d76299565005cf2b8f7f83652fcd.
Actual `SUBMISSION.json` pins the four manifests, runner and commands.
Each wrapper always collects its own terminal/failed attempt into
`shard_N/COLLECTION.json`. Do not duplicate running collection or resubmit jobs.
The exact final compare command is in COMMANDS.md. The frozen fresh strict32
reference is shared and unchanged; no calibration runs in this transfer.

Root reviewed scope/code and explicitly authorized these four jobs. Historical
old precision, released static, DFT and all other methods remain in the joined
ledger alongside the strict candidate. All75 strict group summaries and100
correlated La triples are preserved. Production stays unchanged.

Durable report-only observer PID1618192 is active, with exact scope/pins in
`run_v1/OBSERVER.json`. It waits for these four jobs, recovers only missing
terminal collections, then writes COMPARISON.json, SUMMARY.json, COSTS.json,
REPORT.md, FINAL_TESTS.log and FINISH.json. No molecular retries or new jobs.
