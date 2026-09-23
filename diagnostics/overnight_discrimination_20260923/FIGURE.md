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

## Completed union/adaptive comparison

All four real collections are now terminal. The actual three-method output is
`workspaces/discrimination_transfer_figures_20260923/union225_v1/`, again with
editable SVG/PDF, PNG, all-source CSV and receipt. Its PNG was visually inspected.
The union/adaptive result is205 correct/0 wrong/1 inconclusive/19 unavailable;
the plot retains C5AXV8's inconclusive call and both numerical coverage failures.
It shows that spread is not uniformly reduced. Use the source values and the
complete comparison ledger for counts; no unseen/missing source is inferred.

Exact rerender operation, using a new directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/plot_discrimination_transfer.py \
  --comparison workspaces/union_adaptive_20260923/transfer225_v1/COMPARISON_v1.json \
  --methods 'context_composite=Released static' \
    'native_minimal_recovered=Original-context adaptive' \
    'union_adaptive=Consistent pocket + adaptive' \
  --output workspaces/discrimination_transfer_figures_20260923/union225_v2
```

## Completed uniform-precision candidate

`workspaces/discrimination_transfer_figures_20260923/precision225_v1/` adds the
separately calibrated numerical candidate alongside released static and the
original union/adaptive method. The PNG was visually inspected; editable SVG/PDF,
all675 method/source rows and pinned receipt are present. The candidate has206
correct/0wrong/2inconclusive/17unavailable. Q4W6G0's new inconclusive and C5AXV8's
retained inconclusive are visible; recovered coverage is not concealed.

The numerical candidate uses a different frozen reference. NativeGFN2 component
sensitivity remains; neither identical numerical scores nor universally reduced
structural spread is implied. The figure is a completed development comparison,
not a production promotion or independent biological validation.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/plot_discrimination_transfer.py \
  --comparison workspaces/slsqp_precision_transfer_20260923/COMPARISON_v1.json \
  --methods 'context_composite=Released static' \
    'union_adaptive=Consistent pocket + adaptive' \
    'union_precision=Cheaper numerical candidate' \
  --output workspaces/discrimination_transfer_figures_20260923/precision225_v2
```
