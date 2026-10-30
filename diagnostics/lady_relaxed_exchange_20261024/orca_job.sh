#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch -J name -o dir/slurm_%j.out orca_job.sh <rundir>
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
cd "$1"
NP=$SLURM_NTASKS
sed -i "s/@NPROCS@/$NP/" job.inp
echo "host $(hostname) nprocs $NP"
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
