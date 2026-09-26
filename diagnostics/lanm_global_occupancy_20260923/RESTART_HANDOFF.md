# Nikasha restart handoff — 2026-09-25

> September26 update:1217219 is terminal; both seeded vacuum endpoints failed
> after500cycles. No new calculation is queued by this session. Read
> [the current Pro review](../../docs/PRO_REVIEW_LANM_20260926.md) and CURRENT.md.
> The queued-job/resubmission instructions below are historical and must not
> trigger an automatic relaunch of this failed branch.

**Read this first after a restart.** This is the compact state of Jacob's long
session, including the last approved action. Recover actual scheduler/artifact
state before doing anything: the snapshot below can become stale.

## Mission and permissions

Nikasha (formerly Alquemia) provides affordable structure-sensitive La/Ca
discrimination for the integrated PLM hillslope manuscript. The PQQ round is
packaged. Current research asks whether full-chain LanM modeling can distinguish
lanthanides, beginning with a few Hans/Mex sources and two-/four-ion occupancies.
Jacob prioritizes accuracy for LanM and authorized contained pilots, whole-node
CPU/RAM use, completion watchers and email to `jacobwestroberts@gmail.com`.
He removed project-wide compute budgets/per-analysis approval loops, but preserve
scientific scope and current local rules; do not bypass partition limits or
touch others' jobs. Keep responses short. Latest instruction: preserve all context
and vault documentation before biotite restarts.

## Where everything lives

Repository/root:
`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs`

Relative paths below are under that root. CPU Python:
`/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python`.
Do not install over environments or reset/stash the working tree: many unrelated
dirty files belong to other work. Recent scoped commits through this checkpoint:
`9f26861`, `dc24ed9`, `04ad3f8`, `68532c2`; use live HEAD, never reset to an anchor.
The final handoff itself is committed afterward. No remote push in this phase.

## Scientific state, without the conversational history

- PQQ: preserve the working production/static and explicit DFT paths. The strongest
  consumed225-structure accommodation comparison was207correct/0wrong/1inconclusive/
  17unavailable; same-solver static199/1/8/17. These are development results over
  correlated structures, not225 independent proteins. Practical three-source
  validation is separately91/100correct, or94/100 with named recovery; do not
  conflate its protocol with tenfold context. Read `docs/PLM_PQQ_SOP.md` and
  `diagnostics/plm_pqq_delivery_20260923/{REPORT,MANUSCRIPT}.md`.
- Figure/paper text: `diagnostics/discrimination_single_panel_20260923/`.
  Single-panel artifact is under
  `workspaces/discrimination_transfer_figures_20260923/accommodation_only_v1/`.
  It was emailed. User accepted Ca/La-conditioned source-fold distinctions.
- Spicy-Lams: `diagnostics/spicy_lams_inventory_20260923/REPORT.md` locates616assay
  sequences and16reserved-panel complete-chain models. Reserved numerical
  outcomes remain unopened. Throughput is secondary to accuracy for this branch.
- Prepared LanM scope: Hans8DQ2, Hans8FNR, Mex8FNS complete chain-A monomers, each
  EF12,EF23,EF1234. Fixed source-specific waters/protons. No cross-source raw
  energy minimization or equilibrium occupancy inference. The model does not
  represent unfolding/dimerization. Physical Dy spin and native effective
  singlet are explicitly distinct. See `PLAN.md`, pinned `prepared_v1` and
  `scoring_v2` manifests under `workspaces/lanm_global_occupancy_20260923/`.
- MACE first source Hans8DQ2 EF12:129 real evaluations in451.42worker seconds;
  full/native force comparison passed. **Both relaxed candidates broke covalent
  geometry gates** (carboxylate C–C stretches to~2.6–2.9Å) and were rejected.
  No valid accommodation or within-series preference has been demonstrated.
  Do not loosen gates or restart the old eight-system continuation1213040.
- Native origin job1216564:5h52m21s on112CPUs. Both ALPB cells passed native/state/
  strict-SCF checks; La−3040.435409827431Eh, Dy−3039.193803581625Eh. Both vacuum
  cells failed after500iterations. Those final unconverged energies are unusable.
  Earlier1213018/1216461/1216547 were retained MPI/memory/environment startup
  failures. All technical startup fixes are now implemented.

## Current exact work

**1217219**, `nikasha-lanm-vacuum-seeded`, was **PENDING,
QOSMaxJobsPerUserLimit**, at handoff. No runtime manifest existed yet.
It performs exactly two same-metal ALPB-seeded vacuum restarts at the unchanged
Hans8DQ2 EF12 origin coordinates. No successful solvent call repeats.
`SEEDED_VACUUM_PLAN.md` is the scope. Hamiltonian stays vacuum, nativeGFN2,
TolE1e−10,MaxIter500,300K. Matching native GBW+xtbw initial guesses are retained
and require actual `INITIAL GUESS: XTBRESTART` output. Single-seed convergence
does not qualify a unique ground state or a biological preference.

Executable: `run_seed_vacuum.sbatch`; prepares/executes/collects with
`seed_vacuum_recovery.py` and the pinned existing native executor.
Requests one exclusive memory node, full RAM, at least112task slots; actual
allocated CPUs are divided between two workers (56ranks each on112CPUs,
112each on224CPUs). MaxCore uses75% of RAM across ranks. `native_allocation.py`
reads scheduler `MinMemoryNode=0`/node `RealMemory`, **not** the optional missing
`SLURM_MEM_PER_NODE` variable. No oversubscription or altered scheduler settings.

Runtime root: `workspaces/lanm_global_occupancy_20260923/seeded_vacuum_v1/`.
Look for `manifest.json`, `SEEDS_BEFORE.json`, execution receipts and
`COLLECTION.json`. `seeded_vacuum_preflight_v1` is dry-run only, not computation.
Seed sources and valid solvent results remain under `native_feasibility_retry_v2`.
Never erase or overwrite them, the immutable seed copies or failed attempts.

## Monitoring survives chat, but not necessarily a host reboot

- Detached reporter for1217219 hadPID4160763. Its lock prevents duplicates.
  `DELIVERY_WATCHER_1217219.json` records the exact command. PIDs may be reused
  after reboot: check command identity, do not kill by an old PID.
- Reporter writes `NATIVE_RESULT_1217219.md`, `DELIVERY_1217219.json`, a vault
  capture and a summary email after **any** terminal outcome. Slurm END/FAIL
  mail is also configured. Previous1216564 delivery succeeded at local relay.
- `/root/lanm_completion_watch` was the agent observer. Agent messages do not
  reliably awaken idle root and cannot overcome exhausted usage. Do not claim
  an email wakes the assistant. Detached reporting needs no active assistant;
  restart it after a host reboot if absent. Agents may also need recreation.

### First read-only commands

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
git status --short
squeue -h -j 1217219 -o '%i %T %M %R'
sacct -j 1217219 -X --format=JobID,State,ExitCode,Elapsed,AllocCPUS,NodeList -P
ls -l diagnostics/lanm_global_occupancy_20260923/*1217219*
pgrep -af 'watch_native_delivery.py.*--job 1217219'
```

If delivery is absent and the reporter is not running, restart it (safe also for
a completed job; it collects/reports the terminal result, never submits work):

```bash
nohup /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/lanm_global_occupancy_20260923/watch_native_delivery.py --job 1217219 --manifest workspaces/lanm_global_occupancy_20260923/seeded_vacuum_v1/manifest.json --seeded-vacuum >> diagnostics/lanm_global_occupancy_20260923/delivery_1217219.log 2>&1 < /dev/null &
```

If1217219 is pending/running/requeued, **do not submit another job**. Slurm jobs
may survive controller maintenance. If it completed, inspect existing results.
If interrupted, preserve partial files. The executor intentionally refuses to
overwrite incomplete attempts; blindly requeueing a prepared directory will fail.
After confirming the old job is terminal and recovery is still needed, prepare
a fresh attempt via the batch's explicit output argument. This command repeats
both vacuum tasks, so use it only if neither has a validated completion receipt;
if one succeeded, prepare a one-task manifest using that actual success as reuse
instead of repeating it:

```bash
sbatch --parsable --output=diagnostics/lanm_global_occupancy_20260923/seed_vacuum_%j.out --error=diagnostics/lanm_global_occupancy_20260923/seed_vacuum_%j.err diagnostics/lanm_global_occupancy_20260923/run_seed_vacuum.sbatch workspaces/lanm_global_occupancy_20260923/seeded_vacuum_resume_20260925
```

Use a new fresh path if that named resume directory already exists. Attach the
reporter to the **returned new job ID and corresponding new manifest**, record
both in CURRENT.md. Do not reuse the1217219 watcher for another job.
No automatic new guess, increased cycle limit, relaxed convergence threshold,
geometry/protonation change or larger campaign is authorized by this recovery.
If seeded SCF fails, diagnose/propose the next scientific choice. If it succeeds,
report native numerical feasibility separately from the still-invalid relaxed
MACE candidates; revisit covalent-preserving global motions before expansion.

## Concurrent PLM work: recover its own current state

Another session owns the large PLM folding/scoring campaign. Do not relaunch it
or change its jobs from this handoff. Authoritative live state:
`/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-23_Nikasha_PQQ/CURRENT_STATE.json`.
Existing cohort workspace: `workspaces/plm_pqq_cohort_20260923/`.
Older counts in this chat are historical; use per-target receipts/latest state.
PQQ SOP, references, historical DFT references for Colin's paper, and original
classifications must remain accessible and unchanged.
