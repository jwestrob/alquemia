#!/usr/bin/env bash
# cleanup_orca_scratch.sh
#
# Walks <workspace>/<stem>_qm/ dirs and deletes ORCA scratch files for
# candidates whose La AND Ca single-points have terminated normally.
# apo/water scratch is *also* deleted (cleanup is per-candidate, not per-job).
#
# Idempotent: re-running on already-cleaned dirs is a no-op.
#
# Usage:
#   bash scripts/cleanup_orca_scratch.sh [--dry-run]
#
# Safety:
#   - Only touches <workspace>/*_qm/ at the top level.
#   - SKIPS qmmm*, colin_handoff, b97_3c_panel, calexcitin_size_panel,
#     and anything outside the *_qm pattern.
#   - Requires sp_<stem>_La.out AND sp_<stem>_Ca.out to both contain
#     "ORCA TERMINATED NORMALLY" — otherwise leaves the dir untouched.
set -euo pipefail

WORKSPACE="/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs"

DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    -h|--help)
      echo "Usage: $0 [--dry-run]"
      echo "  --dry-run   Report what would be deleted without deleting."
      exit 0 ;;
    *)
      echo "Unknown arg: $arg" >&2
      exit 2 ;;
  esac
done

# Glob patterns to delete. Each is matched relative to the candidate dir.
# These are file *names*, not paths — safer than wildcard recursion.
DELETE_GLOBS=(
  "*.gbw"
  "*.densities"
  "*.densitiesinfo"
  "*.bas[0-9]"
  "*.bas[0-9][0-9]"
  "*.tmp"
  "*.tmp.[0-9]"
  "*.tmp.[0-9][0-9]"
  "*.bibtex"
  "*.cpcm"
  "*.cpcm_corr"
  "*.cpcmgrad*"
  "*.cpscfdata"
  "*.cpscfdata.tmp"
  "*.diisao.tmp.0"
  "*.diise.tmp.0"
  "*.diist.tmp.0"
  "*.EIJ.tmp"
  "*.en.tmp"
  "*.E.tmp"
  "*.FAO0.tmp"
  "*.FMO0.tmp"
  "*.G0.tmp"
  "*.gpot0.tmp"
  "*.grho0.tmp"
  "*.grid.tmp"
  "*.hostnames"
  "*.H.tmp"
  "*.int.tmp"
  "*.K.tmp"
  "*.opot0.tmp"
  "*.orho0.tmp"
  "*.P0.tmp"
  "*.PAUX.tmp"
  "*.PDAT.tmp"
  "*.PINP0.tmp"
  "*.POLD0.tmp"
  "*.propint.tmp.0"
  "*.shark_grid.tmp"
  "*.SHARKINP.tmp"
  "*.SHARK.KPRSCR.tmp"
  "*.SM12.tmp"
  "*.soscfv.tmp.0"
  "*.SP12.tmp"
  "*.SQRTVEXT.tmp"
  "*.SQRTVJ.tmp"
  "*.S.tmp"
  "*.T.tmp"
  "*.VAUXJ.tmp"
  "*.VCDJ.tmp"
  "*.vcpcm.tmp"
  "*.VEXT.tmp"
  "*.VM1EXT.tmp"
  "*.V.tmp"
  "*.VXC0.tmp"
  # ORCA also writes per-thread .B.<n>.tmp files
  "*.B.[0-9].tmp"
  "*.B.[0-9][0-9].tmp"
)

terminated_normally() {
  # Returns 0 if file exists and contains "ORCA TERMINATED NORMALLY"
  local f="$1"
  [[ -f "$f" ]] || return 1
  grep -q "ORCA TERMINATED NORMALLY" "$f"
}

human_mb() {
  # bytes -> "X.YZ MB"
  local b="$1"
  awk -v b="$b" 'BEGIN { printf "%.2f", b/1048576 }'
}

total_freed_bytes=0
total_dirs_cleaned=0
total_dirs_skipped=0
total_dirs_seen=0

shopt -s nullglob

for qm_dir in "$WORKSPACE"/workspaces/*_qm; do
  [[ -d "$qm_dir" ]] || continue
  base="$(basename "$qm_dir")"

  # Skip legacy qmmm* dirs (paths differ; out of scope)
  case "$base" in
    qmmm|qmmm_*|qmmm) continue ;;
  esac

  total_dirs_seen=$((total_dirs_seen + 1))
  stem="${base%_qm}"

  la_out="$qm_dir/sp_${stem}_La.out"
  ca_out="$qm_dir/sp_${stem}_Ca.out"

  if ! terminated_normally "$la_out" || ! terminated_normally "$ca_out"; then
    total_dirs_skipped=$((total_dirs_skipped + 1))
    continue
  fi

  # Build list of files to delete in this dir
  freed_bytes=0
  freed_count=0
  files_to_delete=()
  for pat in "${DELETE_GLOBS[@]}"; do
    for f in "$qm_dir"/$pat; do
      [[ -f "$f" ]] || continue
      sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
      freed_bytes=$((freed_bytes + sz))
      freed_count=$((freed_count + 1))
      files_to_delete+=("$f")
    done
  done

  if (( freed_count == 0 )); then
    # Already cleaned — silent (idempotent)
    continue
  fi

  if (( DRY_RUN )); then
    printf "[dry-run] %-55s freed=%s MB (%d files)\n" \
      "$stem" "$(human_mb "$freed_bytes")" "$freed_count"
  else
    for f in "${files_to_delete[@]}"; do
      rm -f -- "$f"
    done
    printf "%-55s freed=%s MB (%d files)\n" \
      "$stem" "$(human_mb "$freed_bytes")" "$freed_count"
  fi

  total_freed_bytes=$((total_freed_bytes + freed_bytes))
  total_dirs_cleaned=$((total_dirs_cleaned + 1))
done

echo "----------------------------------------------------------------"
if (( DRY_RUN )); then
  printf "DRY RUN: would clean %d / %d dirs, freeing %s MB total\n" \
    "$total_dirs_cleaned" "$total_dirs_seen" "$(human_mb "$total_freed_bytes")"
else
  printf "Cleaned %d / %d dirs, freed %s MB total\n" \
    "$total_dirs_cleaned" "$total_dirs_seen" "$(human_mb "$total_freed_bytes")"
fi
printf "Skipped %d dirs (La or Ca SP not terminated normally)\n" \
  "$total_dirs_skipped"
