# Deviations from PLAN.md — second-shell test, stage A

Recorded by the lead agent at submission, 9 November 2026 (local) / 10 November UTC. The scientific
design in PLAN.md is unchanged.

1. **MPI layout of the job wrappers.** PLAN.md's Files section says the wrappers were changed to
   `-c 24` with ranks from `$SLURM_CPUS_ON_NODE`. The preparing subagent found that this gives OpenMPI
   one slot and added an untested oversubscribe setting. The lead agent restored the layout that
   completed the La partner jobs: `-N 1 --ntasks=24 -c 1`, ranks from `$SLURM_NTASKS`. The allocation is
   still 24 CPUs under `standard_expand`, with 24 ORCA ranks using all of them. The oversubscribe
   setting was removed. `restart_job.sh` keeps the STOP-file check.
2. **Watcher.** `watch_stage_a.sh` was added after the plan was frozen. It runs on the login node,
   enforces the plan's cost rule (1,000 allocated CPU-h summed over all four job IDs, or 30 h wall for
   any job), and runs the frozen `second_shell_analyze.py` when every job is terminal. The cap can be
   raised only by writing a number into `workspaces/lady_second_shell_20261109/CAP_CPUH`, which the
   lead agent does on Jacob's word alone.
3. **Report authorship.** Reports for this diagnostic are written by the lead agent, because the
   harness returns subagent findings as text.
