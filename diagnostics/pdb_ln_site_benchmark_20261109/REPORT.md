# PDB lanthanide-site benchmark and the reachable-donor detector — 9 November 2026

**Result.** On lanthanide-free partner structures the reachable-donor detector v0 places 42% of
lanthanide sites with at least three protein donors within its top five centres (57% for single-chain
sites). Counting Asp and Glu Cβ atoms within 8 Å does as well: v0 − B2 = +0.025 [−0.033, +0.088].
**Frozen reading: not better than counting acidic residues.** The existing O/N density scanner reaches
about half the recall of the acidic count (0.19 against 0.40).

Two products stand apart from the detector's verdict:

1. **A labelled benchmark**: 879 lanthanide sites with at least two protein donors from 410 crystal
   entries, each classified as single-chain (621), symmetry-bridging (215) or inter-chain (43), with
   1,631 sites transferred onto 413 lanthanide-free deposits of the same proteins.
2. **A capacity curve**: the fraction of candidate centres that hold a lanthanide rises from 10% at
   three reachable Asp/Glu to 46% at four and 65% at five.

Authorized by Jacob on 9 November (exploratory mining tracks). Written by the lead agent from the
analysis subagent's reply and `results/stats.log`, which the lead read line by line for the numbers
below. `PLAN.md` sha256 `c9dc815c…7853`, frozen before any scoring. No function or native metal is
assigned to any protein.

## Benchmark (observed)

- 410 lanthanide entries (X-ray, ≤ 2.5 Å) hold 1,616 lanthanide ions; 879 have ≥ 2 protein donor
  atoms, in 327 entries, 126 clusters at 95% identity and 119 at 30%.
- **Stripped tier** (the same entry with metals, waters and ligands removed): 846 sites in 518 chains,
  125 proteins. Protein donors per site: 2: 327; 3: 123; 4: 135; 5: 62; 6: 108; 7: 38; ≥ 8: 53.
- **Site classes**, from the full crystal with symmetry images. Among stripped sites with ≥ 3 donors:
  299 single-chain, 193 symmetry-bridging, 27 inter-chain. 191 stripped sites have another lanthanide
  within 5 Å.
- **Partner tier** (primary): 413 lanthanide-free deposits at ≥ 95% identity, 768 chains, 1,631
  transferred sites (1,125 with ≥ 3 donors), 91 of 126 proteins. Ca-bound partners 145 entries (39
  proteins); metal-free 168 (58); other-metal 100 (45). Median local superposition RMSD 0.25 Å; 79
  transfers exceed 2 Å.

## Detector and baselines

- **v0**: candidate centres from the geometry of donor-capable residues; for each centre, the residues
  whose side chain could place a donor in the lanthanide first-shell band given each residue type's
  reach from Cβ and a backbone clash filter; weights anionic O 1.0, neutral O 0.5, N 0.25. No
  rotamer library, no tuning.
- **B1**: the existing O/N density scanner. **B2**: Asp and Glu Cβ atoms within 8 Å. **B3**: random
  solvent-accessible points.

## Recall

Sites with ≥ 3 protein donors; hit = a ranked centre within 4 Å of the lanthanide position; k = 1 / 3
/ 5; bootstrap over 30%-identity clusters.

| Tier, sites | v0 | B2 | B1 | B3 |
|---|---|---|---|---|
| Partner, all (1,125 sites, 72 proteins) | 0.22 / 0.36 / 0.42 | 0.14 / 0.31 / 0.40 | 0.09 / 0.17 / 0.19 | 0.00 / 0.01 / 0.01 |
| Partner, single-chain (580 sites, 45 proteins) | 0.29 / 0.49 / 0.57 | 0.22 / 0.43 / 0.54 | 0.12 / 0.24 / 0.28 | 0.01 |
| Stripped, all (519 sites, 103 proteins) | 0.30 / 0.43 / 0.48 | 0.21 / 0.39 / 0.46 | 0.20 / 0.27 / 0.32 | 0.01 |
| Stripped, single-chain (299 sites, 70 proteins) | 0.38 / 0.56 / 0.63 | 0.28 / 0.52 / 0.60 | 0.29 / 0.38 / 0.43 | 0.02 |

- **Frozen cell** (partner, ≥ 3 donors, k = 5): v0 0.424 [0.324, 0.524]; B2 0.399 [0.301, 0.492];
  v0 − B2 +0.025 [−0.033, +0.088]. Single-chain: v0 0.568 [0.457, 0.685]; v0 − B2 +0.026
  [−0.064, +0.116].
- **Where v0 − B2 excludes 0**: partner, all sites, k = 1: +0.086 [+0.025, +0.156]; metal-free
  partners, k = 1: +0.077 [+0.007, +0.161]; other-metal partners, k = 3: +0.112 [+0.025, +0.205].
- **By partner kind**, k = 5, all sites: Ca-bound 0.50 [0.36, 0.64]; metal-free 0.36 [0.23, 0.49];
  other-metal 0.44 [0.31, 0.58]. Single-chain: 0.63, 0.48, 0.59.
- **By site class**, raw hit rate at k = 5 on partners: single-chain 0.56, symmetry-bridging 0.14,
  inter-chain 0.03.
- The held-out split (40% of clusters set aside as development, unused because v0 has no tuning) and
  the ≥ 2-donor stratum give the same picture.

## Capacity curve

Stripped tier, all v0 centres. x = Asp/Glu residues of one chain that can reach the centre;
y = fraction of centres with an observed lanthanide within 4 Å.

| Reachable Asp/Glu | Centres | Fraction with a lanthanide [95%] | Single-chain sites only |
|---|---|---|---|
| 0 | 188,000 | 0.002 [0.002, 0.004] | 0.001 |
| 1 | 158,974 | 0.014 [0.011, 0.023] | 0.005 |
| 2 | 27,022 | 0.027 [0.015, 0.054] | 0.011 |
| 3 | 2,655 | 0.10 [0.07, 0.16] | 0.068 |
| 4 | 321 | 0.46 [0.36, 0.62] | 0.41 |
| 5 | 106 | 0.65 [0.49, 0.89] | 0.64 |
| 6 | 34 | 0.62 | 0.59 |
| 7 | 18 | 0.61 | 0.61 |
| ≥ 8 | 22 | 1.00 | 1.00 |

Intervals for six and above span most of the range because few clusters contribute. The observed
step sits between three and four reachable carboxylates. Every structure here comes from a crystal
that contains lanthanide, so these fractions are conditional on lanthanide exposure.

## Why v0 misses

- Median score at true positions is 2.0 (partner) and 2.5 (stripped) against 4.5 for the top centre;
  a chain has a median of 593 centres, many of them high-capacity points with no lanthanide.
- Ten sampled partner misses: six are symmetry-bridging or inter-chain sites where the scored chain
  supplies two to four donors; two are 1BIX sites transferred onto 5WN1/5WN3 that land beside
  Leu318 with a score of 0 (probable transfer error, UNVERIFIED); one is a metal-free partner whose
  side chains had rearranged; one is a single-chain site ranked eighth.
- The reach model is consistent: 1,491 of 1,505 same-chain donor residues are reachable at the true
  positions, and 12 of 519 true positions fail the clash filter (stripped tier).

## Composition-matched null

Fifty proteins, 200 shuffles each of donor-capable residue identities among exposure-matched
positions.

- Protein-wide maximum score: median percentile 0.66; 12% of proteins exceed the null's 95th
  percentile.
- Best score within 4 Å of an observed site: median percentile 0.98; 64% exceed the 95th; observed
  mean 3.87 against a null mean of 2.25.

Residue arrangement concentrates capacity at the real sites, and the protein-wide maximum carries
composition and little else (inference from the two lines above).

## What this means for mining (inference)

- Local acidic-residue density carries the positional signal; the reach geometry adds at k = 1 only.
- A protein-level score built on the maximum capacity would mostly rank composition. A count
  threshold at a centre (four or more reachable Asp/Glu) is the quantity the curve supports.
- Thirty-seven percent of the ≥ 3-donor sites in lanthanide crystals bridge symmetry mates. Any
  benchmark built from soaked crystals needs that label.

## Limits

- Unlabelled high-capacity centres count as negatives; some may bind a lanthanide.
- Copies within the biological assembly are unresolved: every donor from another crystal copy counts
  as symmetry-bridging.
- 223 transfers failed the eight-Cα-pair rule and were dropped; some retained transfers are suspect.
- Each chain is scored in isolation; partners extend to 3.0 Å resolution.
- BioMetAll, the published backbone-preorganization finder, was outside this run.

## Compute

Slurm jobs 69983 (scoring, 10 min) and 69990 (null, 28 min), `--qos=standard -c 12`, 7.7 CPU-hours
allocated. 651 new mmCIF downloads plus 218 reused from the imprint audit.

## Files

- `PLAN.md`, `PLAN.sha256`; `scripts/` (`discover_bench.py`, `fetch.py`, `common.py`, `sites.py`,
  `units.py`, `detect.py`, `score.py`, `stats.py`, `null.py`, `smoke.py`, two sbatch wrappers).
- `results/`: `stats.log` (every cell), `ln_sites.tsv` (labelled sites with donor origin),
  `site_ranks.tsv`, `null.jsonl`, `splits.tsv`, `partner_kinds.tsv`, `partner_selection.json`.
- `workspaces/pdb_ln_site_benchmark_20261109/` (untracked): structures, `units.json`, `scores.jsonl`.
