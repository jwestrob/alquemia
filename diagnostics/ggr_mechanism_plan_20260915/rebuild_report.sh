#!/bin/bash
# Recollect approved immutable results; does not launch scientific calculations.
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: bash diagnostics/ggr_mechanism_plan_20260915/rebuild_report.sh NEW_ABSOLUTE_WORKSPACE_DIRECTORY' >&2
    exit 2
fi
GGR_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
GGR_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
GGR_REPORT=$1
case "$GGR_REPORT" in "$GGR_ROOT"/workspaces/*) ;; *) echo 'Output must be a new absolute workspace path' >&2; exit 2;; esac
if [ -e "$GGR_REPORT" ]; then
    echo 'Refusing existing output directory' >&2
    exit 2
fi
cd "$GGR_ROOT"
mkdir -p "$GGR_REPORT"
GGR_WORK=$GGR_ROOT/workspaces/ggr_mechanism_20260915
GGR_OLD_A=$GGR_ROOT/workspaces/baseline_benchmark_20260915/run_v1/collection_1199299.json
GGR_OLD_B=$GGR_ROOT/workspaces/benchmark_set_20260915/ready_tasks_v4/collection_1199508.json
for GGR_STAGE in a b; do
    "$GGR_PY" scripts/ggr_workflow.py collect --manifest "$GGR_WORK/stage_${GGR_STAGE}_tasks_v1/manifest.json" --output "$GGR_REPORT/collection_${GGR_STAGE}.json"
done
"$GGR_PY" scripts/ggr_sensitivity.py collect --manifest "$GGR_WORK/stage_c_tasks_v1/manifest.json" --output "$GGR_REPORT/collection_c.json"
GGR_COLLECTIONS=(--collection "$GGR_REPORT/collection_a.json" --collection "$GGR_REPORT/collection_b.json")
GGR_HISTORY=(--historical-collection "$GGR_OLD_A" --historical-collection "$GGR_OLD_B")
GGR_SENSITIVITY=(--sensitivity "$GGR_REPORT/collection_c.json")
GGR_MANIFESTS=(--manifest "$GGR_WORK/stage_a_tasks_v1/manifest.json" --manifest "$GGR_WORK/stage_b_tasks_v1/manifest.json" --manifest "$GGR_WORK/stage_c_tasks_v1/manifest.json")
if [ -f "$GGR_WORK/stage_c_half_tasks_v1/manifest.json" ]; then
    "$GGR_PY" scripts/ggr_sensitivity.py collect --manifest "$GGR_WORK/stage_c_half_tasks_v1/manifest.json" --output "$GGR_REPORT/collection_c_half.json"
    GGR_SENSITIVITY+=(--sensitivity "$GGR_REPORT/collection_c_half.json")
    GGR_MANIFESTS+=(--manifest "$GGR_WORK/stage_c_half_tasks_v1/manifest.json")
fi
"$GGR_PY" scripts/ggr_workflow.py compare "${GGR_COLLECTIONS[@]}" "${GGR_HISTORY[@]}" --output "$GGR_REPORT/comparison.json"
"$GGR_PY" scripts/ggr_component_audit.py "${GGR_COLLECTIONS[@]}" "${GGR_HISTORY[@]}" --output "$GGR_REPORT/components.json"
"$GGR_PY" scripts/ggr_execution_audit.py "${GGR_MANIFESTS[@]}" --output "$GGR_REPORT/execution_cost.json"
"$GGR_PY" scripts/ggr_benchmark_export.py --parent-release "$GGR_ROOT/workspaces/benchmark_set_20260915/scored_release_1199508/benchmark_manifest.json" "${GGR_COLLECTIONS[@]}" "${GGR_SENSITIVITY[@]}" --output "$GGR_REPORT/benchmark"
"$GGR_PY" scripts/ggr_mechanism_plot.py "${GGR_COLLECTIONS[@]}" "${GGR_HISTORY[@]}" "${GGR_SENSITIVITY[@]}" --comparison "$GGR_REPORT/comparison.json" --output-dir "$GGR_REPORT/plots"
echo "$GGR_REPORT"
