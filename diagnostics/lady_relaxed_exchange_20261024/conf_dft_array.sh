#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --array=1-5
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: conf_dft_array.sh <confdir> <base job.inp>; task i re-optimizes CREST
# conformer i at the production level (normal Opt + Freq).
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
C="$1"; i=$SLURM_ARRAY_TASK_ID
[ -s "$C/crest_conformers.xyz" ] || { echo "no conformers"; exit 0; }
nat=$(head -1 "$C/crest_conformers.xyz"); nconf=$(( $(wc -l < "$C/crest_conformers.xyz") / (nat + 2) ))
[ $i -le $nconf ] || { echo "only $nconf conformers"; exit 0; }
R="$C/dft_c$i"; mkdir -p "$R"; cd "$R"
sed -n "$(( (i-1)*(nat+2) + 1 )),$(( i*(nat+2) ))p" ../crest_conformers.xyz > start.xyz
sed -e 's/ TightOpt / Opt /' -e "s/nprocs [0-9@A-Z]* end/nprocs $SLURM_NTASKS end/" "$2" > job.inp
grep -q "start.xyz" job.inp || { echo "bad base input"; exit 1; }
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
