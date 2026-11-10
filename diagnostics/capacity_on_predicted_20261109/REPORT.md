# The four-carboxylate rule on predicted structures — 9 November 2026

**Result.** Frozen reading: **does not transfer.** On ESMFold2 models of the benchmark's own proteins
the detector still ranks the true lanthanide positions (top-five recall 0.552 [0.427, 0.682] against
0.60 on the stripped crystals), but the count of reachable Asp/Glu at those positions drops: 11.9% of
true positions reach four on the models against 25.7% on the crystals (ratio 0.46 [0.20, 0.72];
threshold 0.8). The loss sits at low-confidence sites. As a hard threshold on predicted structures,
"four or more" loses about half of the sites it finds on crystals; as a ranking it holds.

Three further observations:

1. **Low false-positive burden.** Centres with four or more reachable Asp/Glu away from any observed
   lanthanide: 0.09 per 100 residues on crystals and 0.13 on models; the median protein has none.
2. **Reference proteins.** Lanmodulin, lanpepsy and LanD sites reach four to six on crystals and keep
   that on models at 16 of 21 sites. The XoxF and PedH active sites reach three, so the rule misses
   cofactor-assisted sites on crystals as well.
3. **Repeat-protein candidates.** Twelve of 37 unannotated acidic-loop repeat families and 7 of 11
   FG-GAP / TSP_3 families carry period-matched centres of capacity four or more; none of 10
   TPR families does.

Authorized by Jacob on 9 November (exploratory mining tracks). Written by the lead agent from the
analysis subagent's reply and `results/stats.log`. `PLAN.md` sha256 `e16399ac…f4b2`, frozen before
any fold was scored. Candidate identifiers and sequences come from unpublished datasets and stay out
of this public repository; they are in `~/Documents/Berkeley/Lanthanides/capacity_on_predicted_2026-11-09/`
on Jacob's Mac.

## Test: the benchmark's proteins, folded

Fifty proteins, one per 30%-identity cluster, each with a single-chain lanthanide site of at least
three protein donors and a chain of at most 1,000 residues (seeded draw); 70 sites. ESMFold2 through
the Biohub API, single sequence. Models superposed locally on the crystal chain; sites transferred;
frozen v0 and B2 scored. Re-scoring the stripped crystals reproduced the benchmark's ranks at 70 of
70 sites.

| Tier | Sites / proteins | v0 recall at 1 / 3 / 5 | B2 recall at 1 / 3 / 5 | Mean reachable Asp/Glu at the true position | Fraction with ≥ 4 |
|---|---|---|---|---|---|
| Stripped crystal | 70 / 50 | 0.43 / 0.55 / 0.60 | 0.24 / 0.50 / 0.61 | 2.67 | 0.26 |
| Partner structure | 256 / 28 | 0.28 / 0.47 / 0.56 | 0.26 / 0.48 / 0.56 | 2.17 | 0.13 |
| Predicted model | 70 / 50 | 0.27 / 0.42 / 0.55 | 0.17 / 0.36 / 0.48 | 2.14 | 0.12 |
| Predicted, trimmed to crystal-observed residues | 70 / 50 | 0.26 / 0.44 / 0.55 | 0.19 / 0.38 / 0.50 | 2.14 | 0.12 |

- Frozen numbers: recall at 5 on models 0.552 [0.427, 0.682] (threshold 0.45, passed); fraction with
  ≥ 4 on models 0.119 [0.040, 0.208] against 0.257 [0.143, 0.380] on crystals, ratio 0.462
  [0.195, 0.720] (threshold 0.8, failed).
- Paired recall at 5, model minus crystal: −0.051 [−0.188, +0.083].
- Per site, the model's count is lower than the crystal's at 19 sites, equal at 51, higher at none.
- Predicted models match the lanthanide-free partner structures on every column.

**By model confidence at the site residues**

| Mean pLDDT at the site | Sites / proteins | Fraction ≥ 4: crystal → model | Recall at 5: crystal → model |
|---|---|---|---|
| Below 70 | 11 / 8 | 0.62 → 0.00 | 0.81 → 0.31 |
| 70 and above | 59 / 46 | 0.20 → 0.14 | 0.57 → 0.58 |

**Capacity curve on these proteins.** At four reachable Asp/Glu, 29 of 57 centres hold a lanthanide on
crystals (0.51) and 18 of 57 on models (0.32); at five, 9 of 12 (0.75) and 5 of 8 (0.62).

## False-positive burden

Merged centres (within 4 Å of each other) lying more than 8 Å from every observed lanthanide.

| Tier | Off-site centres with ≥ 4 | Residues | Per 100 residues | Off-site centres with ≥ 5 |
|---|---|---|---|---|
| Stripped crystals | 11 | 14,048 | 0.08 | 1 |
| Predicted models | 20 | 15,001 | 0.13 | 0 |
| Predicted, trimmed | 16 | 14,048 | 0.11 | |

The median protein has zero in every tier.

## Reference proteins with lanthanide-bound crystal structures

Reachable Asp/Glu at each known site, crystal → model (21 sites; precursor sequences folded).

| Protein (PDB) | Sites |
|---|---|
| XoxF (4MAE, 6OC5, 6DAM, 7O6Z) | 3→3, 3→3, 3→3, 4→4 |
| PedH (6ZCV) | 3→2 |
| Mex-LanM (8FNS) | 5→5, 4→3, 4→4, 3→3 |
| Hans-LanM (8DQ2) | 4→4, 4→4, 4→2 |
| Al-LanM (9PI2) | 4→4, 4→4, 4→3 |
| LanD (9C90, symmetry-bridging site) | 4→4 |
| Lanpepsy (9VY4) | 6→6, 5→5, 6→6, 4→3, 6→6 |

Sixteen equal, five lower, none higher. Fourteen of 21 sites reach four on the crystal and 10 on the
model.

## Repeat-protein candidates (descriptive)

One representative region per family from the acidic-loop table of the repeat-ladder scan; proteins
longer than 1,200 residues folded as the repeat region plus 50 residues each side.

| Group | Models | Median pLDDT | Centres with ≥ 4 per 100 residues (median) | Models with a period-matched centre of ≥ 4 | Period-matched centres, total |
|---|---|---|---|---|---|
| Unannotated families | 37 | 83.5 | 0.19 | 12 | 69 |
| FG-GAP / TSP_3 / TSP3_bac families | 11 | 83.2 | 0.95 | 7 | 45 |
| TPR families | 10 | 88.9 | 0.00 | 0 | 0 |

At most three repeat units supply donors to one centre (median 1.5–2). Whether period matching
exceeds chance spacing is untested.

## Compact EF-hand proteins (lead addition, descriptive)

Reachable Asp/Glu per EF-hand loop on models (PROSITE PS00018 loops; full-length sequences).

- Reference lanmodulins: 5, 5, 5 (Mex-LanM), 4 (Hans-LanM), 3 (Al-LanM); the pattern finds three, one
  and one loop in them.
- Five SR-VP lanmodulin homologs: all 20 loops reach four or more, 12 reach five.
- Four Omnitrophota proteins of the compact EF-hand family: 6 of 16 loops reach four, 10 reach three.
- A 6 Å window can take in residues from the paired loop.

## Limits

- Single-sequence folds; an MSA might raise confidence at the low-pLDDT sites.
- Fifty proteins and 70 sites; the low-confidence stratum has 11 sites.
- The crystals come from lanthanide-containing experiments, so "holds a lanthanide" is conditional on
  exposure.
- Lanthanide binding of every candidate is UNVERIFIED.

## Compute

127 of 150 sequences folded (50 benchmark, 10 reference, 58 candidates, 9 EF-hand); all succeeded. No
Slurm jobs; about 0.25 CPU-hours single-threaded on the login node.

## Files

- `PLAN.md`, `PLAN.sha256`; `scripts/` (`select_bench.py`, `score_pred.py`, `stats_pred.py`,
  `ef_score.py`).
- `results/`, tracked: `bench_eligible.tsv`, `bench_sample.tsv`, `bench_sample.faa`, `bench_sites.tsv`,
  `bench_fp.tsv`, `capacity.tsv`, `ref.faa`, `ref_sample.tsv`, `ref_sites.tsv`, `samples.sha256`.
- `results/`, untracked (unpublished identifiers and sequences): `stats.log`, `cand.faa`,
  `cand_sample.json`, `cand_models.tsv`, `cand_centres.tsv`, `ef.faa`, `ef_sample.tsv`, `ef_loops.tsv`,
  `failures.json`.
- `workspaces/capacity_on_predicted_20261109/` (untracked): frozen module copies, models, logs.
