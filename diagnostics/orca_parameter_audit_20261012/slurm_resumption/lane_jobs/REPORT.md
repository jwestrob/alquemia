# Frozen CPU proposal lanes on Slurm

The unchanged Ca-conditioned MACE batch now has a QOS-compliant Slurm transport.
Root submitted lane jobs **1028/1029**, collector **1030**, and downstream scalar
stage **1031** on 12 October 2026. This report records implementation and launch;
it does not claim new molecular results. Inspect the actual receipts for progress.

## Preserved scientific work

The clean drain preserved 51 of 106 terminal sources, leaving 55 sources for the
same two disjoint lanes. The original pinned `CPU_BATCH.py` supplies the model,
four-mode selector, optimizer, native geometry admission, six-cell common pool,
source reuse and collection. No scientific settings or decision bands changed.
Each new warm worker runs the original Ca/La qualification fixtures; all four
checks must pass before either lane starts new source optimization.

The older 14-CPU submission was rejected by the actual 12-CPU/job QOS and created
no job. The original source partitions now run as two independent seven-CPU jobs,
14 CPUs in aggregate within the 60-CPU/user limit. This changes placement only.
The account/QOS and scientific kernels remain unchanged.

## Exact execution record

Current plan: `workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/lane_jobs_v2/PLAN.json`.
It pins wrapper, original implementation, original manifest and actual drain,
and reserves original batch `attempts/003`. The earlier unused metadata plan v1
and attempt002 remain preserved. **Do not launch either plan again.**

`CHAIN.json`, `SUBMISSION.json`, per-job submission files and watcher records in
the current plan directory contain the actual commands, job IDs and wake setup.
Both source lanes precede the collector through `afterany`; scalar GFN2 follows
the collector through `afterok`.

`lane0.sbatch` and `lane1.sbatch` request seven tasks on distinct healthy nodes
`node-64-768g-1` and `node-64-768g-2`. The collector requests one task on the first
node. No script supplies memory, time, account or QOS overrides. The partition's
actual default is unlimited memory per node, so this is **not a memory-sharing
claim**. Distinct nodes prevent the two waiting qualification workers from
blocking one another through whole-node default memory reservations. Runtime
receipts retain actual scheduler RAM, task counts, CPU affinity and host.

## Ownership and failure behavior

The metadata-only login coordinator holds the existing `RUN.lock`, reads actual
Slurm states, and releases the common qualification barrier. It cancels only its
two explicitly registered startup jobs if qualification fails before release.
After release, a failed lane does not cancel its independent sibling. Root owns
submission, completion watches and any subsequent recovery decisions.

The collector runs scientific collection under Slurm after both jobs terminate.
It preserves partial sources and explicit unavailable scientific pools. Incomplete
source coverage or a failed lane returns failure and blocks dependent scoring;
an explicitly unavailable pool remains a valid terminal record. Original complete
source results and endpoint calls are reused, and interrupted searches are not
silently restarted. The frozen worker's historical placement field is supplemented
by the actual new scheduler receipts.

## Verification and limitations

Python compilation, shell syntax checks, exact pin checks, actual drain checks,
and read-only healthy-node/partition queries passed. Metadata preparation made
zero molecular calls. Root reviewed the complete adapter; a further read-only
review is recorded separately if received. These checks establish orchestration
readiness, not successful molecular execution or classifier improvement.

The coordinator and completion watchers must survive for automatic progression.
After interruption, inspect actual jobs, plan and existing receipts before any
recovery; do not duplicate jobs or replay completed science.
