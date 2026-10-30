#!/bin/bash
#SBATCH -A banfield
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH -J lady_fep
# usage: sbatch -o log window_job.sh <la.prmtop> <dy.prmtop> <rst7> <T> <prod_ns> <outroot> <seed> <lam1> [lam2 ...]
# Runs every listed lambda window concurrently on the one allocated GPU.
P=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm
D=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lady_explicit_solvent_20261026
la=$1 dy=$2 rst=$3 T=$4 ns=$5 root=$6 seed=$7; shift 7
nvidia-smi -L
for l in "$@"; do
  o=$root/T${T}/lam$(printf %.3f $l)
  mkdir -p $o
  $P/bin/python $D/run_window.py $la $dy $rst $l $T $ns $o $seed > $o/run.log 2>&1 &
done
wait
echo done
