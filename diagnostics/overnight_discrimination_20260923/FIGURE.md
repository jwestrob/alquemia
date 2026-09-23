# Existing source-structure comparison, editable figure

`scripts/plot_discrimination_transfer.py` renders an existing comparison ledger.
It performs no new molecular work, fits no threshold and uses the existing
complete balanced aggregates. Every original source is exported to CSV, including
unavailable values. Counts and decisions are checked against the stored frozen
bands before plotting.

The first actual output uses the completed recovered minimal-adaptive comparison:
`workspaces/discrimination_transfer_figures_20260923/recovered_v1/` contains SVG,
PDF, PNG, source CSV and a hash-pinned receipt. The PNG was visually inspected:
labels, source conditioning, missing counts, wrong/inconclusive markers and
strict aggregate markers are legible. The real ledger has225 sources in25
protein groups, released203/2/2/18 versus adaptive204/1/1/19. These are consumed
structural repeats, not225 independent proteins.

Each panel subtracts its own frozen decision-band midpoint from R. This removes
only an arbitrary display offset; it does not rescale the energy, put different
methods on a common affinity scale or change their calibration. The grey band
shows that method's inconclusive region. Blue/purple protein labels denote the
known Ca/La classes; circle/triangle markers denote Ca/La-conditioned sources.
Black open diamonds are the existing equal means of the two complete arm
medians. Missing members prevent that aggregate from appearing.

Run from the repository root; choose a new output directory when repeating:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/plot_discrimination_transfer.py \
  --comparison workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json \
  --methods 'released=Released static' 'minimal_recovered=Original-context adaptive' \
  --output workspaces/discrimination_transfer_figures_20260923/recovered_v2
```

The running union/adaptive comparison will be plotted only after all four real
collections are terminal. No projected result is shown in the current figure.
