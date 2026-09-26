# Field-aware metal response — restart checkpoint

2026-09-26. Root owns integration/execution. New user handoff authorizes the
staged Ca/La environmental-response scout and supersedes the previous discussion
pause. Historical whole-protein LanM vacuum subtraction remains closed.

## Current state (updated 2026-09-26 23:33 UTC)

- **Reference job1218751 queued afterany final discriminator collector1217591.**
  Six native analytic-gradient endpoints,64CPU task slots/256GiB, four16-rank
  workers, no GPU. All six cells currently unavailable/pending execution.
- **Collector1218752** runs afterany1218751, independently of this chat/login.
  It writes FINAL_COLLECTION.json, final_accounting.txt and AUTO_REPORT.md in
  the reference workspace. It launches no dependent scientific experiment.
- Upstream audit complete: code supports atomwise potential and field, but
  exact trained checkpoint and paper-compatible engine interface unavailable.
  Author question drafted in upstream/AUDIT.md and **not sent**.
- Real 1H4I preparation complete:54atoms,9087permanent charges, source-matched
  Ca−3/La−2 singlets; only Thr159HG1 rotates+10degrees in environmentB.
- Twelve unit/parser/allocation tests and five real preparation checks pass.
  No new molecular evaluations, native force qualification or predictive result.
- All three agents finished and released ownership; no duplicate job submission.
- Root owns interaction ledger, reference harness, numerical gates and submission.
- Production PQQ/PLM belongs to another session; preserve its jobs and dirty files.
- Reading anchor: git `87f1134`; existing unrelated tracked edits were present.

## Scheduling and interruption

Jacob says biotite shuts down Monday 28 September at 10 AM PST. Plan for the
earlier interpretation, 10 AM local Pacific, and do not rely on a surviving
login process. Research must queue **behind discriminator jobs**. Current
discriminator terminal collector is 1217591; other active chain IDs are
1217426,1217440,1217453,1217454,1217588,1217589,1217590. Requery queue and
dependency graph before submission; use explicit dependencies, not priority
changes to another session's jobs. Pending1218751 holds no running allocation.

Every submitted manifest, job ID, command, receipt and collection status will
be recorded here. No automatic rerun after restart: inspect scheduler and
receipts first. Keep partial and failed attempts. No email or remote push is
requested by this phase.

## Exact artifacts and next work

Workspace: `workspaces/metal_environment_response_20260926/reference_scout_v1/`.
Manifest SHA256: `c822b82da3857e45ead36a345d6130d8138833c56e065113dc83ed4ad3526696`.
Read SUBMISSION.json, COLLECTOR_SUBMISSION.json and DRY_RUN.json. Actual Python
dependencies and allocation-aware renderer are frozen under implementation/;
production scripts/environments were not altered. Native parser was exercised
on an actual archived embedded GGR output, not fabricated output.

Full additive QM/MM remains blocked by cross metal/PQQ repulsion/dispersion and
boundary-reference coverage. Scout measures only the declared **electronic
embedding component**, not total hybrid energy or affinity. PLAN.md freezes
its Hamiltonian, A/B rule and numerical tolerances. The complete model gate is
not passed by running a component test.

After1218751/1218752 finish, read actual collection/accounting. If complete,
declare finite native directional checks (MM plus mapped boundary), repeat and
rigid-transform cases, and expanded-region qualification before interpretation
or extension. No ML call without an exact accessible checkpoint.
No automatic4MAE/nonPQQ/LanM expansion, whole-protein optimization or promotion.
New molecular cost is zero while dependency-pending; preparation/code audit
time is separate and not a molecular timing benchmark.

Recovery commands (from repository root):

```bash
cat diagnostics/metal_environment_response_20260926/CURRENT.md
squeue -u jwestrob -o '%.18i %.50j %.10T %.30P %.30R'
git status --short --untracked-files=no
sacct -j 1218751,1218752 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,MaxRSS,ExitCode
cat workspaces/metal_environment_response_20260926/reference_scout_v1/AUTO_REPORT.md
```

If the collector is interrupted, first inspect jobs/receipts, then collect to
a new filename without rerunning chemistry:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/metal_environment_response_20260926/reference_scout_v1/implementation/metal_environment_reference.py collect \
  --manifest workspaces/metal_environment_response_20260926/reference_scout_v1/manifest.json \
  --output workspaces/metal_environment_response_20260926/reference_scout_v1/RECOVERED_COLLECTION_1.json
```

If that filename exists, inspect it rather than overwrite it. Partial/failed
endpoints require an explicit fresh attempt, never an automatic rerun.
Startup exit75 means new discriminator work or Monday's shutdown cutoff
prevented molecular execution; read the job log. The cutoff is a shutdown
safeguard, not a project compute/time budget.
