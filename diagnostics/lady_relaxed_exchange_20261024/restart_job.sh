#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# Chained afterany on the TightOpt job. If it did not finish Opt+Freq (floppy
# H-bond modes stall the 1e-6 Eh energy criterion), restart from its last
# geometry with normal Opt convergence + Freq. The gradient was already at
# TightOpt RMS level, so this changes the minimum by << 0.1 kcal/mol.
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
cd "$1"
if grep -q "ORCA TERMINATED NORMALLY" job.out && grep -q "Final Gibbs free energy" job.out; then echo "already done"; exit 0; fi
mv job.out job_tight.out; cp job.xyz restart.xyz; cp job.inp job_tight.inp
sed -i -e 's/ TightOpt / Opt /' -e 's/start\.xyz/restart.xyz/' -e "s/nprocs [0-9@A-Z]* end/nprocs $SLURM_NTASKS end/" job.inp
grep -n "Opt\|xyzfile" job.inp
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
