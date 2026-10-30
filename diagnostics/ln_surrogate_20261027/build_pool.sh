#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -c 24
#SBATCH -J ln_build
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch -o log build_pool.sh <workspace>   (tasks.txt lines: lig_id smiles; claims by mkdir)
W=$1; D=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ln_surrogate_20261027
B=/groups/banfield/users/jwestrob/envs/lnsel/bin/python
M=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
TH=4; NW=$(( SLURM_CPUS_ON_NODE / TH )); export OMP_NUM_THREADS=$TH OMP_STACKSIZE=1G MKL_NUM_THREADS=$TH
mapfile -t L < $W/tasks.txt
for ((w=0; w<NW; w++)); do ( for line in "${L[@]}"; do set -- $line; o=$W/nist/$1
  [ -f $o/DONE ] && continue; mkdir -p $o; mkdir $o/.claim 2>/dev/null || continue
  echo "$SLURM_JOB_ID w$w $(date +%FT%T)" > $o/.claim/by
  ( cd $o && $B $D/build_complex.py $1 "$2" $o && $M $D/relax_mace.py $o $TH ) > $o/log 2>&1 && touch $o/DONE
  echo "w$w $1 rc=$? $(date +%T)"; done ) & done; wait
