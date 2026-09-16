# Reproduce the GGR result figures

`scripts/ggr_mechanism_plot.py` reads existing collections and exports four
standalone PNG/PDF figures, a pinned plot manifest, and two TSV tables of the
plotted values. It runs no scientific executable and refuses existing output
directories. Missing results remain explicitly unavailable.

The completed A/B collections can be rendered from the repository root:

```bash
python scripts/ggr_mechanism_plot.py \
  --collection workspaces/ggr_mechanism_20260915/stage_a_tasks_v1/collection_1199770.json \
  --collection workspaces/ggr_mechanism_20260915/stage_b_tasks_v1/collection_1199802.json \
  --historical-collection workspaces/baseline_benchmark_20260915/run_v1/collection_1199299.json \
  --historical-collection workspaces/benchmark_set_20260915/ready_tasks_v4/collection_1199508.json \
  --comparison diagnostics/ggr_mechanism_plan_20260915/COMPARISON_AB.json \
  --output-dir workspaces/ggr_mechanism_20260915/stage_ab_plots_replay
```

Add `--sensitivity` with the actual nominal Stage C collection when available.
The option can be repeated for an actually executed conditional half-step
collection; nominal and half-step values remain separate. The figure shows
analytic `h × gradient`, actual `[E(+h) − E(−h)] / 2`, and the declared numerical
acceptance intervals. Those intervals are not statistical uncertainty.
In the exported gradient table, `units` and `gradient_units` describe the
projected derivative; `amplitude_units` separately describes the displacement
amplitude in Å or radians.

The comparison input is optional. When provided, its pinned current collections
must match the supplied collections and its alpha–GGR contrasts must agree with
the plotted direct algebra. Input/output energies and receipts are checked
against their archived records. Source/preparation protocols are never mixed
within an alpha–GGR contrast.

Actual A/B render verified visually:
`workspaces/ggr_mechanism_20260915/stage_ab_plots_v1/`.
The first three figures contain completed A/B results; its sensitivity figure
explicitly states that no sensitivity collection was supplied. Earlier Stage A
only renders are retained as technical render checks, not additional scientific
experiments.

All figures concern already-consumed development cases. They introduce no
classification threshold or relaxation correction. Structural replicates are
not independent biological observations.
