#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -c 24
#SBATCH -J ln_series
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# usage: sbatch -o log series_pool.sh <workspace>   (every nist/<id> with La_cn*_mace.xyz; claims by mkdir .sclaim)
W=$1; D=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ln_surrogate_20261027
M=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
TH=4; NW=$(( SLURM_CPUS_ON_NODE / TH )); export OMP_NUM_THREADS=$TH OMP_STACKSIZE=1G MKL_NUM_THREADS=$TH
for ((w=0; w<NW; w++)); do ( for o in $W/nist/aqua $(ls -d $W/nist/*/ | grep -v aqua); do o=${o%/}
  [ -f $o/SDONE ] && continue; [ -f $o/DONE ] || continue
  mkdir $o/.sclaim 2>/dev/null || continue
  $M $D/relax_series.py $o $TH > $o/slog 2>&1 && touch $o/SDONE
  echo "w$w $(basename $o) rc=$? $(date +%T)"; done ) & done; wait
