# Corrected Ca-source solvent stage on Slurm

Prepared only; no calculation or submission by this adapter's author. This adds
Slurm placement to the frozen `ca_conditioned_execution/GFN2.py` driver (SHA256
`c8c92032ca0e08541f8e4501da0d92a724b492f297d773d2c6dfbc3083d26a20`).
It does not change that driver, its manifest, parser, Hamiltonian, component gates,
125-source denominator, three corrected reuses, exclusions or frozen calibration.
The existing `context_factorial_gfn2.py` and `REEVALUATE.py` helper hashes are also
checked before import, so delayed execution cannot silently use edited helpers.

The stage runs prepare → dry-run → execute → collect inside the allocation.
It requires the actual terminal MACE collection. Submit it with an `afterok`
dependency on the integrator's Slurm MACE continuation job; a running/partial
MACE collection fails preparation. Watch both predecessor and solvent jobs so an
upstream failure cannot leave an unnoticed dependency waiting forever.

## Execution and recovery

The wrapper requires real `SLURM_JOB_ID` and `SLURM_CPUS_ON_NODE`, supports one
node, and derives concurrency from allocated CPUs and the number of pending
inputs. Root explicitly requested 64 CPU slots on a healthy 64-CPU node; the
wrapper rejects a different allocation or a CPU affinity mask smaller than that
allocation. The script requests `--ntasks=64 --nodes=1 -w node-64-768g-2` on
`standard`, because the current scheduler does not automatically allocate every
CPU with partition defaults. Read-only `scontrol` verified that node IDLE with
64 effective CPUs and no existing allocation. It is outside the known bad-node
list. This permits 64 independent one-thread inputs when enough remain.
Each input still runs the pinned xTB then independent tblite
backend sequentially through unchanged `context_factorial_gfn2.execute_point`.
Each executor thread is bound to one distinct inherited CPU; child processes
inherit this affinity and all numerical thread environment settings equal one.
No molecular timeout, new iteration limit, RAM request, CPU-per-task flag or
walltime is added. Partition defaults apply.

Both preparation and execution have nonblocking ownership locks. Existing raw
RESULT files are verified/reused; an existing incomplete point directory is
left unavailable and never automatically repeated. The inherited pending-task
check refuses live unfinished process groups. A partial preparation without a
manifest also fails visibly rather than overwriting its directory. Subsequent
attempts get new receipts/collections; previous results remain immutable.

Before calls, a full preexecution collection records every unavailable result.
After calls, the same frozen collector computes components, paired decisions and
measured backend wall/CPU time and peak process RSS. Wrapper receipts separately
record actual allocation, CPU affinity, worker count, elapsed time before and
after final collection and corresponding allocated CPU-seconds; the latter are
not utilization. The shell saves a Slurm
accounting snapshot on exit; root should obtain terminal accounting as usual.
A killed job may leave no final collection; preserved START/raw records plus the
preexecution collection make that interruption visible. Resumption computes only
never-started inputs. Failed native pools cannot silently become origin scores.

## Integrator command

The integrator supplies its actual predecessor job ID through `GFN_PREDECESSOR`.
Run from the canonical repository root; do not run the Python stage on login.

```bash
test -n "$GFN_PREDECESSOR"
sbatch --dependency="afterok:$GFN_PREDECESSOR" \
  diagnostics/orca_parameter_audit_20261012/slurm_resumption/GFN2_STAGE.sbatch \
  "$PWD/workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json" \
  "$PWD/workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/FINAL_COLLECTION.json" \
  "$PWD/workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_slurm_v1"
```

The input collection path must match the actual MACE continuation output; root
owns that choice and the completion watcher. No new optimization or production
mutation belongs to this stage. All scientific fixture/integration tests remain
for execution through Slurm; local verification is syntax/read-only review only.
