#!/bin/bash
# Launch the Nd/Lu series Opt+Freq runs as dependency chains: N_SE standard_expand slots (24 ranks) and
# N_STD standard slots (12 ranks). Each state = TightOpt+Freq job + afterany restart; the next state in a slot
# waits for the previous restart. Size the slots to the queue load at launch time (squeue --me).
# Usage (after the La/Dy panel has finished):  python build_series.py && ./submit_series.sh N_SE N_STD
set -e
R=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
D=$R/diagnostics/lady_relaxed_exchange_20261024
W=$R/workspaces/lady_series_20261028
NSE=${1:?N_SE}; NSTD=${2:?N_STD}; NSLOT=$((NSE + NSTD))
ORDER=(ttha_w0 dtpa_w1 cdta_w3 cdta_w2 tmdta_w3 tmdta_w2 edtax_w3 edtax_w2 nta_w5 nta_w4 mal_w7 mal_w6 acbi_w8 acbi_w7 acbi_w6 acmono_w8 acmono_w7 aq9 aq8)
states=(); for s in "${ORDER[@]}"; do for m in Nd Lu; do [ -f $W/runs/${m}_$s/job.inp ] && states+=(${m}_$s); done; done
[ -s $W/JOBS.txt ] && mv $W/JOBS.txt $W/JOBS_$(date +%Y%m%d_%H%M%S).txt
declare -a last; i=0
for st in "${states[@]}"; do
  slot=$((i % NSLOT)); i=$((i + 1))
  if [ $slot -lt $NSE ]; then q="--qos=standard_expand --ntasks=24"; qn=standard_expand; else q="--qos=standard --ntasks=12"; qn=standard; fi
  dep=""; [ -n "${last[$slot]}" ] && dep="--dependency=afterany:${last[$slot]}"
  m=$(sbatch --parsable $q $dep -J sr_$st -o $W/runs/$st/slurm_%j.out $D/orca_job.sh $W/runs/$st)
  r=$(sbatch --parsable $q --dependency=afterany:$m -J srr_$st -o $W/runs/$st/slurm_restart_%j.out $D/restart_job.sh $W/runs/$st)
  last[$slot]=$r
  echo "$st $qn $m $r" | tee -a $W/JOBS.txt
done
echo "${#states[@]} states in $NSLOT slots; arm: $W/watch_series.sh (run_in_background)"
