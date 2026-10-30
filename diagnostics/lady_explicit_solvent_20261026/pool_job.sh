#!/bin/bash
#SBATCH -A banfield
#SBATCH -J lady_fep_pool
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: [NW=<workers>] [OMM_PLATFORM=CPU] sbatch -p ... [--gres=gpu:1] -c N -o log pool_job.sh <tasks.txt>
# Shared dynamic queue: every worker in every pool job walks the same task list and claims a
# window by atomic mkdir <outdir>/.claim; finished windows (meta.json) are skipped.
# GPU jobs: NW workers share the allocated GPUs (OpenCL). CPU jobs (OMM_PLATFORM=CPU): NW workers,
# each with SLURM_CPUS_ON_NODE/NW threads.
# tasks.txt lines: <la.prmtop> <dy.prmtop> <rst7> <lambda> <T> <prod_ns> <outdir> <seed>
P=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm
D=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lady_explicit_solvent_20261026
TASKS=$1
NW=${NW:-$SLURM_CPUS_ON_NODE}
export OMM_THREADS=$(( SLURM_CPUS_ON_NODE / NW ))
NG=1; [ "${OMM_PLATFORM:-}" != CPU ] && { nvidia-smi -L; NG=$(nvidia-smi -L | wc -l); }
echo "host $(hostname) workers $NW threads/worker $OMM_THREADS platform ${OMM_PLATFORM:-OpenCL}"
mapfile -t L < "$TASKS"
for ((w=0; w<NW; w++)); do
  (
    for line in "${L[@]}"; do
      set -- $line
      o=$7
      [ -f $o/meta.json ] && continue
      mkdir -p $(dirname $o) && mkdir $o 2>/dev/null || [ -d $o ] || continue
      mkdir $o/.claim 2>/dev/null || continue
      echo "$SLURM_JOB_ID $(hostname) w$w $(date +%FT%T)" > $o/.claim/by
      OMM_DEVICE=$((w % NG)) $P/bin/python $D/run_window.py "$@" > $o/run.log 2>&1
      echo "w$w rc=$? $o $(date +%T)"
    done
  ) &
done
wait
echo done
