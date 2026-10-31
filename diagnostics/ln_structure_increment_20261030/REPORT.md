# Nikasha structural increment — 30 October 2026

1. **What predicts selectivity best?** Among the independently fitted ablations,
   donor counts plus ligand charge (M1): family-held-out RMSE **0.854 kcal/mol**
   and MAE **0.641**, on 115 ligands. Its family-bootstrap RMSE interval is
   0.694–0.981. Published-checkpoint errors are overlapping-data references,
   not a competing independent generalization result.
2. **Did topology add value?** No supported increment: full-data RMSE changes
   0.854 → 0.882. The paired RMSE improvement is −0.028 kcal/mol, with a 95%
   family-bootstrap interval of −0.057 to +0.026.
3. **Did 3D add value beyond topology?** Not supported. On the same 87 ligands,
   M2 → M3 changes RMSE 0.998 → 0.955: improvement +0.042, interval
   −0.042 to +0.087. With the predetermined softer restraint it instead worsens
   to 1.060. Neither numerical target is met, and the increment is unstable.
4. **Was the published-model comparison independent?** **No.** All 115 evaluation
   ligands have designated-family exposure in the released M1/M4 data; M4 has
   normalized matches for 103 and protonation-record matches for 96. Its SMILES-only
   ensemble has raw RMSE **0.563**, but that does not meet an independent validation
   gate. The clean architecture comparison remains unresolved.
5. **What single next action is justified?** **Target missed: preserve M1 as the
   honest baseline and end automatic expansion of this implementation.**
   This failure does not establish that structure contains no predictive information.

## Frozen test and data

[PLAN.md](PLAN.md), [groups.tsv](groups.tsv), [splits.json](splits.json) and
[FREEZE.json](FREEZE.json) specify target, groups, features, grids and source pins.
The current checkout remains at reviewed commit c3b250c; unrelated dirty changes
were preserved. No newer commits were present.

The target is directly tabulated, matched La(III)/Dy(III) K1=[ML]/([M][L]), 1:1,
25 °C, identical ionic strength and native species/condition annotations:
`y = -1.3642470130342386*(logK_Dy-logK_La)` kcal/mol. Negative means Dy preference.
No partial-series targets or measured affinity/series-anchor features enter.

There are **115 primary chemical examples**, **135 linked matched conditions**,
**11 outer chemical families**. One primary condition per ligand is chosen by
proximity to 0.1 M without using its affinity. All other conditions are retained
and assigned the same ligand/family fold. Ninety-four of the historical 101
ligands survive the stricter curation. Unresolved structures, named stereo missing
from historical SMILES, censored/flagged constants, duplicate native rows and
incompatible annotations remain explicit exclusions.

The 15 tetraacetate structures were recovered, including the cyclic cis/trans
assignments from systematic OPSIN diamines followed by two acetate arms per N,
with formula checks. Fourteen targets are admitted. TMDTA is excluded because La
and Dy use different native adjustments (`yv` versus `pv`); equivalence is not
assumed. TTHA is also excluded (different `yv` versus `v` adjustments), so its
historical extreme residual is not retested here. Every admitted ligand has at
most eight potential donors: the excess-donor feature is constant in this curated
iteration, and over-dentate behavior is not empirically qualified. Relative cis/trans stereochemistry is verified on every saved cyclic
candidate, including the 1,4 isomers. These ligands were a previously inspected
**development challenge**, not prospective validation.

A post-initial-results chemistry audit repaired the grouping code for nitrobenzoate,
placing it with benzoate relatives. This is disclosed in [DEVIATIONS.md](DEVIATIONS.md).
All initial freezes and scores are preserved under the workspace's `grouping_v1/`.
Only that ligand's primary group changes; no identities, targets, conditions,
features or numerical criteria change. The current grouping is therefore a
transparent metadata repair, not untouched prospective registration.

NIST-derived row-level records remain in the ignored workspace, following the
repository's existing redistribution policy. [Records](../../workspaces/ln_structure_increment_20261030/evaluation_records.tsv),
[primary rows](../../workspaces/ln_structure_increment_20261030/primary.tsv),
[exclusions](../../workspaces/ln_structure_increment_20261030/exclusions.tsv) preserve
exact native byte offsets, raw condition codes, ligand species, annotations and
structure provenance. Individual original-study bibliography links have not been
decoded; provenance is the critically selected NIST compilation and exact native
records. This is a limitation, not an assertion that compiled values are independent
replicates from the same original experiment.

## Ablations and uncertainty

Outer leave-one-family-out; inner leave-one-family-out selects ridge alpha from
0.01, 0.1, 1, 10, 100. Scaling is fitted inside each training fold. Explicit broad
aminocarboxylate grouping keeps all tetraacetates, IDA/NTA/DTPA/TTHA relatives out
together. No Murcko-only grouping, stacked model, feature search or checkpoint
fine-tuning is used. Confidence intervals resample whole families 2,000 times.

| Evaluation population | Model | n | RMSE | MAE |
|---|---|---:|---:|---:|
| Full | M0 mean | 115 | 1.339 | 1.009 |
| Full | M1 donors + charge | 115 | **0.854** | **0.641** |
| Full | M2 + topology | 115 | 0.882 | 0.646 |
| Same size-fit subset | M0 | 87 | 1.452 | 1.116 |
| Same size-fit subset | M1 | 87 | 1.031 | 0.738 |
| Same size-fit subset | M2 | 87 | 0.998 | 0.729 |
| Same size-fit subset | M3 + size fit | 87 | **0.955** | **0.703** |

All errors are kcal/mol. Full-data M1 and subset M3 have different denominators
and must not be ranked against one another. On the shared subset, topology's
+0.034 RMSE improvement has interval −0.021 to +0.145. The small data and broad
families leave substantial uncertainty; individual rows are not bootstrap units.

The **0.6 family-held-out target fails**. On the 14 tetraacetate development targets,
M1/M2 RMSE are 1.385/1.478, failing **0.8**. On the 11 with size-fit coverage,
M1/M2/M3 RMSE are 1.798/1.823/1.712, also failing. Stronger restraint or fewer
conformers are sensitivity results, not replacements for the primary protocol.

[Row-level predictions](../../workspaces/ln_structure_increment_20261030/predictions.tsv),
[paired errors](../../workspaces/ln_structure_increment_20261030/paired_errors.tsv),
[per-family errors and counts](../../workspaces/ln_structure_increment_20261030/per_family.tsv),
[inner fits/scalers](../../workspaces/ln_structure_increment_20261030/fit_metadata.json),
and [alternate conditions](../../workspaces/ln_structure_increment_20261030/alternate_condition_predictions.tsv)
are retained. The historical donor-only LOO result, **0.85 kcal/mol**, reproduces
with an exactly equal 101-row historical table; it is a different validation design
and dataset from this iteration's 0.854 result.

## Descriptor coverage and bounded sensitivity

M1 uses the existing 12 SMARTS counts, macrocycle indicator and L charge. M2 adds
plausible 5/6/7-member chelate paths, donor excess over eight, donor-path rotors,
ring constraints and macrocycle size. Counts refer to potential donors, not measured
coordination. Geometry uses only organic MMFF94 plus a dummy center; it supplies
no MMFF metal interaction parameters. Four deterministic conformers and up to eight
predetermined assignments are shared between metal sizes. Over-dentate candidates
bind at most eight selected donors. Only organic energies enter the difference;
artificial restraint energies are saved separately.

La/Dy CN9 Shannon radii 1.216/1.083 Å, with empirical O/N offsets 1.40/1.55 Å,
give O targets 2.616/2.483 and N targets 2.766/2.633 Å. These radius-sum distances
are documented descriptor assumptions, not validated bond parameters.
[Shannon (1976)](https://doi.org/10.1107/S0567739476001551).

Coverage is **87/115 (75.7%)**, including 11/14 tetraacetate targets. Failures:
12 disconnected assignments, 8 ambiguous nitroso/phenol tautomers, 4 relative-stereo
ambiguities, 3 ambiguous deprotonations/tautomers, 1 MMFF parameterization failure.
They are missing values, never zeros. Coverage by family and every candidate's
energies, coordinates, assignments, convergence and distance residuals are saved.

The fixed sensitivity changes k=20 to 10 or 40 kcal/mol/Å², or uses the first two
of the same four conformers. M3 RMSE is respectively **0.955, 1.060, 0.952, 0.960**.
The descriptor itself has RMS changes of 2.15/4.56/1.69 kcal/mol, with maximum
changes 14.67/33.85/11.39; many monodentate cases remain essentially zero.
This implementation is not stable enough to earn expansion. The stiffness giving
a positive interval is not selected after seeing its scores.

[Real-record checks](../../workspaces/ln_structure_increment_20261030/checks.json)
pass target sign/units, exact historical replay, all 2,088 candidate organic-energy
replays, restraint closure, cyclic relative stereo, acetate/EDTA finite-difference
gradients, ligand grouping and missing-value handling. No synthetic fixtures.

## Observable correction

The [series README](../lady_series_20261028/README.md) now has a dated correction.
The methods literature distinguishes intrinsic metal binding from CD conformational
response and discusses heavy-lanthanide divergence. Both series are preserved;
CD response, intrinsic affinity, adsorption and separation factors remain separate
observables. This validates neither the previous calculations nor Hans selectivity.
[Mattocks, Tirsch & Cotruvo (2021), section 4.1/Table 3](https://www.osti.gov/servlets/purl/1781338).

## Separate peptide illustration

The 2025 SI specifies **FIDTNNDGWYEGDELLA**, free N terminus and C-terminal amide;
N5D is **FIDTDNDGWYEGDELLA**. Tb-displacement measurements use pH6, 50mM MES,
100mM NaCl. The competition-assay temperature is not specified, so the comparison
is in log preference rather than an imputed experimental kcal/mol target.

Table S3 gives Lu/La preference 102.18-fold for WT and 2.45-fold for N5D: the
mutation changes `logK_Lu-logK_La` by **−1.621**. Released M4-S predicts WT/N5D
contrasts 0.673/0.669, hence a change of **−0.004** (ensemble spread 0.041 log).
It misses this total mutation effect. Added charge, coordination and hydrogen
bonding are not disentangled, and the result establishes no general protein
transfer capability. Whole-LanM repeated experimental values were never treated
as independently measured EF-hand labels.
[Primary paper and SI](https://doi.org/10.1002/chem.202501318),
[exact constructs/predictions](../../workspaces/ln_structure_increment_20261030/peptide_illustration.json).

## Published reference, cost and reproduction

Upstream release pinned to `d6099497416e5d3f69c8622fbef16bf17c08590f`,
Chemprop **1.5.2**, Python 3.10.20, CPU PyTorch 2.2.2, RDKit 2023.09.6.
[Full lock](chemprop-lock.txt), [analysis environment](analysis-environment.json),
[asset manifest](assets.json), [additional M4-S assets](additional_assets.json).
M1 and M4 use their published condition/metal features; the SMILES-only M4-S
ignores explicit conditions. Each La/Dy pair nevertheless uses the same experimental
conditions, and differences are formed within ensemble members before aggregation.
M4's 16-feature reference is restricted to the 16 primary ligands with explicitly
identified NaClO4 electrolyte; unreported electrolytes are not imputed.

| Released reference | n | RMSE | MAE | Interpretation |
|---|---:|---:|---:|---|
| M1, graph + temperature/I | 115 | 0.611 | 0.399 | overlapping-data reference |
| M4-S, graph only | 115 | 0.563 | 0.355 | overlapping-data reference |
| M4, full conditions/metal features | 16 | 0.612 | 0.309 | different subset; overlapping-data reference |

These are not out-of-fold predictions. Descriptive family intervals and per-family
errors are retained separately; they do not repair exposure. Do not interpret the
raw 0.563 versus ridge 0.854 difference as demonstrated superiority.

The audit inventories **121,078 released rows** across M1–M4 and the older
LOGKPREDICT training releases, including protonation constants, other metals,
conditions, stereo-agnostic representations and close Morgan analogues.
Normalized training/validation matches occur for **102/115 M1**, **99/115 M2**,
**99/115 M3**, and **103/115 M4** evaluation ligands. All designated families occur
in M1/M4 training/validation. Full tautomer normalization is required for every
heavy-element-composition-compatible exact match; other compositions cannot be
exact matches and receive charge normalization plus the family/analogue audit.
All rows remain inventoried. **215 parse failures** are preserved and prevent
using a negative lookup as proof of universal absence. The 115 evaluation identity
keys also agree between analysis and Chemprop RDKit versions.

[Overlap audit](../../workspaces/ln_structure_increment_20261030/training_overlap.tsv),
[other metals, conditions and protonation](../../workspaces/ln_structure_increment_20261030/training_overlap_conditions.tsv),
[complete inventory](../../workspaces/ln_structure_increment_20261030/released_training_conditions.tsv),
[published predictions](../../workspaces/ln_structure_increment_20261030/published_predictions.tsv),
and [published per-family results](published_per_family.tsv) provide source rows.
No verified external family is established. Family-excluded from-scratch published
architecture refits were **not completed**; this is the exact missing evidence
for a clean architecture comparison. No resource-impossibility claim is made.
Fine-tuning was not used as a substitute.

The original LOGKPREDICT Ni–ethylenediamine example yields **5.714140**, versus
bundled **5.408758**, failing the frozen 0.02-log-unit replay tolerance. The newer
Python API also has a 21-versus-48-feature checkpoint mismatch. Both failures are
preserved; legacy exact numerical reproduction remains unresolved. The separate
2025 Figure 8(b) M4-S application uses the exact released uranyl inputs: acetate
**6.0481**, oxalate **8.7006**, approximately matching the plotted red points.
The source does not tabulate exact ML outputs, so this is not claimed as exact
numeric reproduction. Its preprocessing and checkpoints are documented before
primary La/Dy inference; no settings were tuned against the application errors.
[Application replay](../../workspaces/ln_structure_increment_20261030/application_replay.json),
[2025 example](../../workspaces/ln_structure_increment_20261030/application_M4S/result.csv),
[upstream data/checkpoints](https://github.com/Critical-Materials-Institute/stabilityconstant-ml-models),
[2025 paper](https://doi.org/10.1021/acs.jcim.5c01546),
[LOGKPREDICT paper](https://doi.org/10.1063/5.0176000).

The exact amidated WT/N5D peptide graphs have no exact or stereo-agnostic identical
ligand graphs in the released inventory. Absence of related peptide-family exposure
has not been established; this illustration is not an independent protein-family
validation. [Peptide overlap record](../../workspaces/ln_structure_increment_20261030/peptide_illustration_with_overlap.json).

**Actual compute:** 19 Slurm jobs, **6,297 allocated CPU-seconds (1.749 CPU-hours)**,
**1,964.978 measured CPU-seconds (0.546 CPU-hours)**, zero GPU time. This includes
all cancelled audit attempts, failed examples/checks/collectors and recoveries.
All used `standard_expand`; no owned jobs remain. Source curation, downloads and
environment setup were not metered. [Accounting](COMPUTE_COST.json),
[full scheduler records](../../workspaces/ln_structure_increment_20261030/scheduler_accounting.tsv).

Reproduce on biotite using the pinned local assets/environments and NIST source:

```bash
bash diagnostics/ln_structure_increment_20261030/reproduce.sh
```

It creates a fresh workspace, queues the fixed stages through Slurm and leaves
these results intact. It reuses archived assets, without new scientific protocols.
No DFT, MD, alchemical free energies, force-field fits, ML-potential training,
reserved SpyCI-LAMBS/PLM outcomes, production changes, commits or pushes occurred.

**Decision: Target missed. Preserve M1 and end automatic expansion.**
