# Lanthanide series: cubic polynomial per ligand plus one Gd|Tb gap (10 November 2026)

Read `REPORT.md` first. `PLAN.md` was frozen before any fit and `ADDENDUM_confirmation.md` before the 24
confirmation series were read; hashes and times are in `PLAN.sha256`.

The scripts expect a `data/` folder beside `scripts/` holding copies of three untracked tables:

- `workspaces/ln_surrogate_20261027/surrogate_nist.tsv`
- `workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv`
- `workspaces/ln_selectivity_calibration_20261027/features.tsv`

Order: `run_main.py` (frozen readings), `run_proteins.py`, `run_followup.py`, `fit_cubic_coordinate.py`,
`fit_step.py`, `explore_groups.py`, `yttrium_branch.py`, `run_confirm.py`, `robustness.py`, `make_figure.py`.
Outputs go to `results/` and `figures/`, which hold NIST-derived values and stay untracked.
