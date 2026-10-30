#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH -J ln_cpcm
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch -o log cpcm_pool.sh <ws>/cpcm   (claims by mkdir .claim; one 24-rank ORCA at a time)
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
C=$1
while read t; do o=$C/$t
  grep -q "ORCA TERMINATED NORMALLY" $o/job.out 2>/dev/null && continue
  mkdir $o/.claim 2>/dev/null || continue
  cd $o; sed -i "s/@NPROCS@/$SLURM_NTASKS/" job.inp
  $ORCA/orca job.inp > job.out 2> job.err < /dev/null; echo "$t rc=$? $(date +%T)"
  rm -f job.gbw job.densities job*.tmp
done < $C/tasks.txt
