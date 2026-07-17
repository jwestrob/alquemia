#!/bin/bash
# Long-running recarve watcher: every 60 min, rebuild the recarve queue and
# (if non-empty) invoke process_recarve_queue.py to re-submit any newly-arrived
# empty-carve workspaces.
#
# Lifecycle: write PID to /tmp/recarve_watcher.lock. Stop by deleting the lockfile.
# Concurrent runs are no-ops (PID check).
#
# Cadence rationale: 60 min matches typical SLURM job wallclock for a Ca/La
# pair on the memory partition. The upstream fold_daemon → discriminator
# pipeline drops empties at irregular intervals, so a hourly sweep keeps the
# recarve queue from accumulating without hammering the cluster.
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
LOCK=/tmp/recarve_watcher.lock
LOG=$ALCH/results/recarve_watcher.log
SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SLEEP_SECONDS=3600  # 60 min

echo $$ > "$LOCK"
echo "$(date -Iseconds) | recarve_watcher started (pid=$$)" | tee -a "$LOG"

while [ -f "$LOCK" ] && [ "$(cat "$LOCK" 2>/dev/null)" = "$$" ]; do
  # 1) Rebuild the queue
  n_queue=$("$SCAN_PYTHON" "$ALCH/scripts/rebuild_recarve_queue.py" 2>&1 \
            | tee -a "$LOG" \
            | grep -oE "Rebuilt .*: [0-9]+ new" \
            | awk '{print $(NF-2)}')
  n_queue=${n_queue:-0}

  # 2) If queue has entries, run the processor
  if [ "$n_queue" -gt 0 ]; then
    echo "$(date -Iseconds) | $n_queue entries in queue → running process_recarve_queue.py" | tee -a "$LOG"
    "$SCAN_PYTHON" "$ALCH/scripts/process_recarve_queue.py" 2>&1 \
        | tee -a "$LOG" \
        | grep -E "RECARVED|CARVE_AMBIGUOUS|SOLVENT_EXCLUDE|sbatched|outcome" \
        | tail -20
  else
    echo "$(date -Iseconds) | queue empty — sleeping ${SLEEP_SECONDS}s" | tee -a "$LOG"
  fi

  sleep "$SLEEP_SECONDS"
done

echo "$(date -Iseconds) | recarve_watcher stopped (lockfile gone or PID mismatch)" | tee -a "$LOG"
rm -f "$LOCK"
