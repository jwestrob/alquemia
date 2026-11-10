#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch -J name -o <rundir>/slurm_%j.out orca_job.sh <rundir>
# Same header and MPI layout as diagnostics/lady_relaxed_exchange_20261024/orca_job.sh: 24 MPI ranks on
# 24 CPUs, the layout that completed the La partner jobs (e.g. 24975). ORCA starts its own mpirun,
# which needs one Slurm task per rank, so the allocation is --ntasks=24 -c 1.
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
cd "$1"
NP=$SLURM_NTASKS
sed -i "s/@NPROCS@/$NP/" job.inp
echo "host $(hostname) nprocs $NP cpus_on_node $SLURM_CPUS_ON_NODE job $SLURM_JOB_ID"
grep -n "nprocs" job.inp
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
