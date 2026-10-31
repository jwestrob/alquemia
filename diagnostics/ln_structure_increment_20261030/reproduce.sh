#!/bin/bash
# Fresh output directory; reuse the pinned local assets and existing environments.
set -euo pipefail
ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$ROOT"
D="$ROOT/diagnostics/ln_structure_increment_20261030"
ORIGINAL="$ROOT/workspaces/ln_structure_increment_20261030"
export NIKASHA_WORKSPACE="$ROOT/workspaces/ln_structure_increment_reproduce_$(date -u +%Y%m%dT%H%M%SZ)"
export NIKASHA_RECORD_DIR="$NIKASHA_WORKSPACE/frozen"
mkdir "$NIKASHA_WORKSPACE"
ln -s "$ORIGINAL/sources" "$NIKASHA_WORKSPACE/sources"
ln -s "$ORIGINAL/upstream" "$NIKASHA_WORKSPACE/upstream"
ln -s "$ORIGINAL/upstream_logkpredict" "$NIKASHA_WORKSPACE/upstream_logkpredict"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/envs/lnsel/bin/python "$D/curate.py" > "$NIKASHA_WORKSPACE/curation.log" 2>&1
hist=$(sbatch --parsable -N 1 -n 1 -c 1 --qos=standard_expand --job-name=nik_hist_r --output="$NIKASHA_WORKSPACE/historical_%j.out" --error="$NIKASHA_WORKSPACE/historical_%j.err" "$D/historical.sh")
structure=$(sbatch --parsable -N 1 -n 1 -c 8 --qos=standard_expand --dependency="afterok:$hist" --job-name=nik_size_r --output="$NIKASHA_WORKSPACE/structure_%j.out" --error="$NIKASHA_WORKSPACE/structure_%j.err" "$D/structure.sh")
fit=$(sbatch --parsable -N 1 -n 1 -c 1 --qos=standard_expand --dependency="afterok:$structure" --job-name=nik_fit_r --output="$NIKASHA_WORKSPACE/evaluate_%j.out" --error="$NIKASHA_WORKSPACE/evaluate_%j.err" "$D/evaluate.sh")
checks=$(sbatch --parsable -N 1 -n 1 -c 1 --qos=standard_expand --dependency="afterok:$fit" --job-name=nik_check_r --output="$NIKASHA_WORKSPACE/checks_%j.out" --error="$NIKASHA_WORKSPACE/checks_%j.err" "$D/checks.sh")
pub=$(sbatch --parsable -N 1 -n 1 -c 4 --qos=standard_expand --job-name=nik_pub_r --output="$NIKASHA_WORKSPACE/published_%j.out" --error="$NIKASHA_WORKSPACE/published_%j.err" "$D/published.sh")
pep=$(sbatch --parsable -N 1 -n 1 -c 4 --qos=standard_expand --dependency="afterok:$pub" --job-name=nik_pep_r --output="$NIKASHA_WORKSPACE/peptide_%j.out" --error="$NIKASHA_WORKSPACE/peptide_%j.err" "$D/peptide.sh")
context=$(sbatch --parsable -N 1 -n 1 -c 1 --qos=standard_expand --dependency="afterok:$pub:$pep" --job-name=nik_ctx_r --output="$NIKASHA_WORKSPACE/enrich_%j.out" --error="$NIKASHA_WORKSPACE/enrich_%j.err" "$D/enrich_overlap.sh")
printf '%s\n' "Fresh workspace: $NIKASHA_WORKSPACE" "Jobs: $hist $structure $fit $checks $pub $pep"
printf '%s\n' "$hist $structure $fit $checks $pub $pep $context" > "$NIKASHA_WORKSPACE/jobs.txt"
sbatch --wait --parsable -N 1 -n 1 -c 1 --qos=standard_expand --dependency="afterok:$checks:$context" --job-name=nik_sum_r --output="$NIKASHA_WORKSPACE/summarize_%j.out" --error="$NIKASHA_WORKSPACE/summarize_%j.err" "$D/summarize.sh"
