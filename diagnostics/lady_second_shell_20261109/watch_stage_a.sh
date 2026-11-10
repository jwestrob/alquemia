#!/bin/bash
# Stage A completion watcher. Run on the login node:
#   nohup bash watch_stage_a.sh >> <workspace>/watch.log 2>&1 < /dev/null &
# Cost limits were removed on Jacob's instruction (DEVIATIONS.md item 4): the jobs run until ORCA
# finishes. This watcher only reports state and, once every job is terminal, runs the frozen gate
# analysis. Single-threaded; it polls sacct and does no compute.
R=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
D=$R/diagnostics/lady_second_shell_20261109
W=$R/workspaces/lady_second_shell_20261109
PY=/groups/banfield/users/jwestrob/envs/lnsel/bin/python
POLL=300
field () { sacct -j "$1" -X --format="$2" -n 2>/dev/null | head -1 | tr -d ' '; }
terminal () { case "$1" in COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE|PREEMPTED) return 0;; *) return 1;; esac; }
echo "$(date -u +%FT%TZ) watcher started (no cost limits) pid $$"
while true; do
  all_done=1; total=0; line=""
  while read -r n j1 j2; do
    for j in $j1 $j2; do
      s=$(field "$j" State); raw=$(field "$j" CPUTimeRAW)
      [ -z "$s" ] && s=UNKNOWN
      terminal "$s" || all_done=0
      total=$(( total + ${raw:-0} ))
      line="$line $j:$s"
    done
  done < "$W/JOBS.txt"
  cpuh=$(( total / 3600 ))
  echo "$(date -u +%FT%TZ) cpuh=$cpuh all_done=$all_done$line"
  [ "$all_done" -eq 1 ] && break
  sleep "$POLL"
done
cd "$R" && "$PY" "$D/second_shell_analyze.py" > "$W/analyze_final.log" 2>&1
echo "analysis rc=$? (see analyze_final.log, GATE.json)"
echo "ALL_TERMINAL cpuh=$cpuh" > "$W/WATCH_STATUS.txt"
date -u +%FT%TZ > "$W/WATCH_DONE"
