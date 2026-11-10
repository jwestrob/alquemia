#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch --dependency=afterany:<main> -J name -o <rundir>/slurm_restart_%j.out restart_job.sh <rundir>
# Same as diagnostics/lady_relaxed_exchange_20261024/restart_job.sh plus a STOP-file check (PLAN.md cost
# rule: the watcher touches <rundir>/STOP before cancelling). Runs only if the main job left no Gibbs
# energy; restarts from the last geometry with the same Opt Freq input.
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
cd "$1"
if [ -f STOP ]; then echo "STOP present: no restart"; exit 0; fi
if grep -q "ORCA TERMINATED NORMALLY" job.out && grep -q "Final Gibbs free energy" job.out; then echo "already done"; exit 0; fi
mv job.out job_tight.out; cp job.xyz restart.xyz; cp job.inp job_tight.inp
sed -i -e 's/ TightOpt / Opt /' -e 's/start\.xyz/restart.xyz/' -e "s/nprocs [0-9@A-Z]* end/nprocs $SLURM_NTASKS end/" job.inp
grep -n "Opt\|xyzfile\|nprocs" job.inp
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
