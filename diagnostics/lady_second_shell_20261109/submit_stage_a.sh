#!/bin/bash
# Stage A submission (PLAN.md). Root runs this; it prints and runs each sbatch command.
# Per seed: Opt+Freq job, then an afterany restart job. Job IDs are appended to JOBS.txt.
set -euo pipefail
R=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
D=$R/diagnostics/lady_second_shell_20261109
W=$R/workspaces/lady_second_shell_20261109
cd "$R"
for n in Dy_ssaq_w21_s1 Dy_ssaq_w21_s2; do
  rd=$W/runs/$n
  if [ -e "$rd/job.out" ]; then echo "refusing: $rd/job.out already exists"; exit 1; fi
  grep -q "nprocs @NPROCS@ end" "$rd/job.inp" || { echo "refusing: $rd/job.inp lacks @NPROCS@"; exit 1; }
  c1=(sbatch --parsable -J ss_$n -o "$rd/slurm_%j.out" "$D/orca_job.sh" "$rd")
  echo "+ ${c1[*]}"; j1=$("${c1[@]}")
  c2=(sbatch --parsable --dependency=afterany:$j1 -J ssr_$n -o "$rd/slurm_restart_%j.out" "$D/restart_job.sh" "$rd")
  echo "+ ${c2[*]}"; j2=$("${c2[@]}")
  echo "$n $j1 $j2" | tee -a "$W/JOBS.txt"
done
echo "Next: scontrol show job <id> for each ID above; arm a completion watcher on all four IDs."
