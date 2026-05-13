#!/bin/bash
# Long-running inbox watcher: every 5 min, check inbox; if any CIFs present
# AND queue has slots free, run process_inbox.sh to consume.
# Stop on SIGTERM or by deleting the lockfile.
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
LOCK=/tmp/inbox_watcher.lock
echo $$ > $LOCK

while [ -f "$LOCK" ] && [ "$(cat $LOCK 2>/dev/null)" = "$$" ]; do
  cifs=$(ls $ALCH/inbox/*.cif 2>/dev/null | wc -l)
  total=$(squeue -u jwestrob -h -o "%i" 2>/dev/null | grep -vE "alphap|taxon_scan" | wc -l)
  free_slots=$((200 - total))

  if [ "$cifs" -gt 0 ] && [ "$free_slots" -gt 5 ]; then
    echo "$(date -Iseconds) | inbox=$cifs queue=$total free=$free_slots → running process_inbox.sh"
    bash $ALCH/scripts/process_inbox.sh 2>&1 | grep -E "OK\]|FAIL" | head -50 || true
  else
    echo "$(date -Iseconds) | inbox=$cifs queue=$total free=$free_slots — waiting"
  fi

  sleep 2700  # 45 min
done
echo "Watcher stopped (lockfile gone or PID mismatch)"
rm -f $LOCK
