# Second-shell test — CURRENT (2026-11-10 01:50 UTC)

**Stage A is running.** Submitted 2026-11-10 01:17 UTC by the lead agent after Jacob authorized the test.

| Seed | Main job | Restart job (afterany) | Node |
|---|---|---|---|
| Dy_ssaq_w21_s1 | 69800 | 69801 | node-64-768g-6 |
| Dy_ssaq_w21_s2 | 69802 | 69803 | node-64-768g-7 |

- Both main jobs started at once on 24 MPI ranks; ORCA reports 24 parallel processes, `job.err` is
  empty, and the first SCF is converging. Expected 15–25 h each.
- **Watcher.** `watch_stage_a.sh` runs on the login node (PID in
  `workspaces/lady_second_shell_20261109/WATCH_PID`, log `watch.log`). It reports job state, runs the
  frozen `second_shell_analyze.py` when all four jobs are terminal, and then writes `WATCH_DONE`. The
  gate verdict lands in `GATE.json`.
- **If the watcher is gone** (login-node reboot): restart it with
  `nohup bash diagnostics/lady_second_shell_20261109/watch_stage_a.sh >> workspaces/lady_second_shell_20261109/watch.log 2>&1 < /dev/null &`.
  It reads `JOBS.txt`, so a restart is safe.
- **To finish by hand** once all four jobs are terminal:
  `/groups/banfield/users/jwestrob/envs/lnsel/bin/python diagnostics/lady_second_shell_20261109/second_shell_analyze.py`.
- **No cost cap.** Jacob removed the CPU-hour cap and the 30 h stop on 2026-11-10 UTC (DEVIATIONS.md
  item 4). The jobs run until ORCA finishes.
- **Stage B** (acetate, no cost cap) is released only by a PASS in `GATE.json`. Jacob pre-authorized
  it on that condition on 2026-11-09. FAIL or INCONCLUSIVE ends the test.

See `PLAN.md` (frozen, sha256 `31a7381d…640a`) and `DEVIATIONS.md`.
