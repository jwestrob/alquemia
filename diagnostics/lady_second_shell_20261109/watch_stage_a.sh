#!/bin/bash
# Stage A completion watcher and cost-cap enforcer (PLAN.md, "Cost cap and runtime action").
# Run on the login node:  nohup bash watch_stage_a.sh > <workspace>/watch.log 2>&1 &
# It polls sacct, stops the jobs if the summed allocation reaches the cap or a job passes 30 h wall,
# and runs the frozen gate analysis once every job is terminal. Single-threaded; no compute.
R=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
D=$R/diagnostics/lady_second_shell_20261109
W=$R/workspaces/lady_second_shell_20261109
PY=/groups/banfield/users/jwestrob/envs/lnsel/bin/python
MAX_WALL_S=$((30*3600)); POLL=300
field () { sacct -j "$1" -X --format="$2" -n 2>/dev/null | head -1 | tr -d ' '; }
terminal () { case "$1" in COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE|PREEMPTED) return 0;; *) return 1;; esac; }
while true; do
  CAP=1000; [ -s "$W/CAP_CPUH" ] && CAP=$(head -1 "$W/CAP_CPUH" | tr -cd '0-9')   # Jacob may raise the cap by writing a number here
  all_done=1; total=0; over_wall=0; line=""
  while read -r n j1 j2; do
    for j in $j1 $j2; do
      s=$(field "$j" State); raw=$(field "$j" CPUTimeRAW); el=$(field "$j" ElapsedRaw)
      [ -z "$s" ] && s=UNKNOWN
      terminal "$s" || all_done=0
      total=$(( total + ${raw:-0} ))
      if ! terminal "$s" && [ "${el:-0}" -ge "$MAX_WALL_S" ]; then over_wall=1; fi
      line="$line $j:$s"
    done
  done < "$W/JOBS.txt"
  cpuh=$(( total / 3600 ))
  echo "$(date -u +%FT%TZ) cpuh=$cpuh cap=$CAP all_done=$all_done over_wall=$over_wall$line"
  [ "$all_done" -eq 1 ] && break
  if [ "$cpuh" -ge "$CAP" ] || [ "$over_wall" -eq 1 ]; then
    echo "LIMIT REACHED (cpuh=$cpuh cap=$CAP over_wall=$over_wall): stopping per PLAN.md"
    while read -r n j1 j2; do touch "$W/runs/$n/STOP"; scancel "$j2"; scancel "$j1"; done < "$W/JOBS.txt"
    echo "LIMIT_REACHED cpuh=$cpuh" > "$W/WATCH_STATUS.txt"
    sleep 90; break
  fi
  sleep "$POLL"
done
cd "$R" && "$PY" "$D/second_shell_analyze.py" > "$W/analyze_final.log" 2>&1
echo "analysis rc=$? (see analyze_final.log, GATE.json)"
[ -f "$W/WATCH_STATUS.txt" ] || echo "ALL_TERMINAL cpuh=$cpuh" > "$W/WATCH_STATUS.txt"
date -u +%FT%TZ > "$W/WATCH_DONE"
