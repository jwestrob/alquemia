#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -c 24
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: crest_job.sh <confdir> <charge>; confdir holds seed.xyz (snapshot of the
# partially DFT-relaxed geometry, so M-O distances are already metal-specific).
# CREST/GFN2/ALPB(water) iMTD-GC; default topology check keeps the coordination
# number fixed (waters leaving the first shell are discarded). CREST 3 native mode dies in GFN0 WBO setup
# (DSYGVD illegal value) on these Ln clusters; --legacy uses the xtb binary.
# Flexibility measure is 0 for pure aqua shells, so MTD length is set to 10 ps.
B=/home/jwestrob/miniconda3/envs/lanm_bench/bin
export PATH=$B:$PATH OMP_STACKSIZE=4G OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
ulimit -s unlimited
cd "$1"
# v2: legacy CREST aborts ("safety termination") when its initial xtb opt changes
# the detected topology (M-O right at the bond cutoff); a pre-opt alone does not
# fix it, so --noreftopo is used and the DFT array filters conformers by CN.
# v3: metal-donor distances softly constrained (constrain.inp, auto = current length)
# so CN-preserving H-bond rearrangements are sampled when GFN2 otherwise ejects a water.
xtb seed.xyz --opt --gfn 2 --alpb water --chrg $2 --uhf 0 -P $SLURM_CPUS_ON_NODE --input constrain.inp > preopt.out 2>&1 && cp xtbopt.xyz seed_pre.xyz
crest seed_pre.xyz --gfn2 --alpb water --chrg $2 --uhf 0 -T $SLURM_CPUS_ON_NODE --legacy -mdlen 10 --noreftopo -cinp constrain.inp > crest.out 2>&1
rc=$?; echo "crest gfn2 rc=$rc"
if [ $rc -ne 0 ] || [ ! -s crest_conformers.xyz ]; then
  mkdir -p gfn2_failed; mv crest.out gfn2_failed/ 2>/dev/null
  crest seed_pre.xyz --gfnff --alpb water --chrg $2 --uhf 0 -T $SLURM_CPUS_ON_NODE --legacy -mdlen 10 --noreftopo -cinp constrain.inp > crest.out 2>&1
  echo "crest gfnff fallback rc=$?"; echo gfnff > METHOD
else echo gfn2 > METHOD; fi
