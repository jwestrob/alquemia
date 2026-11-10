# PDB lanthanide-site benchmark for blind "latent metalability" detection — frozen plan

Frozen 2026-11-09 before any detector, baseline or metric was computed on any benchmark structure.
Already done before freezing (benchmark construction only; no scoring): RCSB discovery
(`discover_bench.py`), download of all 410 Ln entries (72 reused from the La/Ca imprint audit
workspace, 338 new), Ln-site labelling (`sites.py` → `ln_sites.tsv`), partner selection
(`partner_selection.json`) and partner download (459 entries, ≤ 1,200 cap).
Execution: discovery, downloads, labelling and unit construction on the biotite login node
(single thread); scoring in at most 3 Slurm jobs (`standard`, `-c 12`); no GPU.

## 1. Benchmark

Ln entries: RCSB X-ray, resolution ≤ 2.5 Å, non-polymer component in {LA CE PR ND SM EU GD TB DY
HO ER TM YB LU} (410 entries). Protein = RCSB 95 % sequence-identity cluster (c95) of the entity;
supercluster = RCSB 30 % cluster (c30).

Site label (per Ln ion, first conformer: altloc blank or A). First shell = O/N atoms ≤ 3.0 Å
(`protein_shells.py` PRIMARY), searched on the full crystal with gemmi NeighborSearch and the
mmCIF cell/space group (symmetry images and lattice translations included). Protein donor = O/N of a
standard amino acid (incl. MSE) of a polymer entity. A **labelled site** has ≥ 2 protein donor
atoms. Recorded per site: protein donor count, donor types (carboxylate O, backbone O, neutral
side-chain O, His N, other O/N), CN (all O/N ≤ 3.0 Å incl. water/ligand), chelating carboxylates
(Asp/Glu with both carboxylate O ≤ 3.0 Å), occupancy, B and B / median protein Cα B.

Donor origin (added at the lead's request before freezing): donors are grouped by (chain,
symmetry image). The site is attributed to the group with most donors (ties: identity image, then
chain name) and its position is expressed in that chain's own coordinate frame. Counts: same chain,
other chain of the same asymmetric-unit copy, symmetry image / lattice translation. Classes:
`single-chain` (one group), `inter-chain` (≥ 2 chains, same copy), `symmetry-bridging` (any donor
from another copy). Biological-assembly membership of symmetry copies is not resolved: every
donor from another crystal copy counts as symmetry-bridging. Flag `cluster5`: another Ln within 5 Å.

Evaluation unit = one protein chain, scored **in isolation** (single chain, standard residues,
everything else removed). This matches the target use (a predicted monomer) and makes the per-chain
donor supply explicit; inter-chain and symmetry-bridging sites remain labels and are reported
separately.

- **Stripped tier** (optimistic): Ln-entry chains carrying ≥ 1 labelled site (≤ 4 chains per entry,
  most labelled sites first, then chain name). Labels = that chain's attributed sites.
- **Partner tier** (primary): for each c95 hosting labelled sites, Ln-free X-ray deposits (≤ 3.0 Å)
  of the same c95, up to 6 Ca-containing, 6 metal-free and 3 other-metal entries by resolution then
  ID (selection by component IDs; final kind from file atoms: Ca atom → `ca`; other metal atom →
  `other_metal`; none → `metal_free`; any Ln/Y/Sc/actinide atom → dropped). Partner chains with
  global sequence identity ≥ 0.95 to a source chain (≤ 4 per entry). Every labelled site of every
  Ln entry of that c95 is transferred: Kabsch fit of source Cα atoms within 12 Å of the site,
  mapped through the alignment (≥ 8 pairs required, RMSD recorded, unfiltered). Transferred sites
  within 2.0 Å are merged; representative = most donors, then single-chain first. Recorded per
  transferred site: metal of the partner within 2.5 Å (Ca / other / none).

Splits: c30 clusters shuffled with `random.Random(42)`; first 40 % = dev, rest = test
(`splits.tsv`). v0 has no tuning; all results are reported on all clusters and on test.

## 2. Detector v0 — reachable donor capacity (frozen, no tuning)

Input: one chain, first conformer. Cβ = observed Cβ or, if missing, ideal virtual Cβ from N/CA/C;
Gly anchor = CA.
Candidate points: 1.0 Å lattice at integer coordinates within 9.0 Å of a Cβ of a donor-capable
residue (D E N Q S T Y H) or within 2.9 Å of a backbone O; removed if any backbone N/CA/C < 3.0 Å,
any backbone O < 2.2 Å or any Cβ < 2.8 Å (backbone clash filter; side-chain atoms beyond Cβ are
ignored because they can move).
Donor band (Ln first shell): [2.2, 2.9] Å.
Side-chain reach from Cβ to the donor atom (ideal geometry, any rotamer): Ser OG, Thr OG1
[1.35, 1.55]; Asp OD, Asn OD1 [2.25, 2.55]; Glu OE, Gln OE1 [2.70, 3.90]; His ND1 [2.40, 2.60],
NE2 [3.50, 3.75]; Tyr OH [5.55, 5.85]. Residue can reach point c when D = |c − Cβ| satisfies
max(0, rmin − 2.9, 2.2 − rmax) ≤ D ≤ rmax + 2.9 and angle(Cβ − Cα, c − Cβ) ≤ 120°.
Backbone carbonyl O counts when |O − c| ∈ [2.2, 2.9] and angle C–O···c ≥ 120°.
Weights: anionic O (Asp, Glu) 1.0; neutral O (Asn, Gln, Ser, Thr, Tyr, backbone O) 0.5; His N 0.25.
Each residue contributes once, at its highest reachable weight. Score S(c) = Σ residues.
Also stored: n_anion(c) = number of Asp/Glu that can reach c.
Ranking: S descending, ties by mean distance from c to the anchors of contributing residues
(ascending). Greedy non-maximum suppression, radius 4.0 Å, over points with S > 0 → ranked centres.

## 3. Baselines

- B1: existing O/N density scanner `on_density_scanner/on_scanner/scanner.scan` unchanged, run on a
  file holding the same chain (standard residues, all altlocs, as the scanner reads them); its own
  ranking and NMS.
- B2: number of Asp/Glu Cβ (anchors as above) within 8.0 Å, evaluated on the same candidate points
  as v0; ties by mean distance to those Cβ; same NMS (4.0 Å) over points with count > 0.
- B3: random solvent-accessible points: 1.0 Å lattice points with no chain heavy atom within 3.0 Å
  and ≥ 1 within 5.0 Å, in random order (`numpy.random.default_rng([42, crc32(unit_id)])`), same NMS.

## 4. Metrics

Hit: a ranked centre within 4.0 Å of the site position. Site hit@k for k = 1, 3, 5 using the top k
centres of its chain. Aggregation: site → unit mean → entry mean → protein (c95) mean; reported
recall = mean over proteins. Partner tier: per partner kind (ca / metal_free / other_metal) and
pooled. Strata: sites with ≥ 3 and with ≥ 2 protein donors; all classes and single-chain only.
Intervals: bootstrap over c30 clusters, 2,000 resamples, `numpy.random.default_rng(42)`, 95 %
percentile; the detector − B2 difference is bootstrapped paired (same resamples).

**Frozen reading** (partner tier, pooled kinds, all site classes, ≥ 3 donors, k = 5): "useful" iff
detector recall ≥ 0.5 AND detector − B2 ≥ 0.10 AND the 95 % interval of the difference excludes 0;
otherwise "not better than counting acidic residues". The same reading is computed and reported
beside it for single-chain sites only (co-reported, lead addition).

## 5. Capacity curve

Stripped-tier units; all v0 centres (after NMS). x = n_anion(c): the number of Asp/Glu of that one
chain that can reach c (single-chain supply by construction), bins 0 … 7, ≥ 8. y = fraction of
centres with an observed Ln (any copy, any donor count) within 4.0 Å. Intervals: c30 bootstrap as
above (pooled fraction per resample). Single-chain version: positives = centres within 4.0 Å of a
labelled single-chain site of that chain; centres near any other Ln excluded.

## 6. Diagnostics (descriptive, no change to v0)

For each stripped-tier labelled site: S and n_anion at the true position (nearest candidate point,
or flagged when the true position fails the clash filter), and the fraction of same-chain donor
residues the reach model marks reachable there. Ten partner-tier ≥ 3-donor misses at k = 5
(sampled with seed 42) inspected for cause.

## 7. Optional — composition-matched null (if time)

50 proteins (seed 42) from partner-tier proteins with a ≥ 3-donor site; one stripped-tier chain
each. Exposure = number of Cα within 10 Å of the anchor, tertile bins within the chain. 200
permutations of residue identities among non-Gly/non-Pro positions within each bin (all residue
types, so donor-capable composition per bin is preserved); recompute the maximum S over the
candidate points. Report the observed max S percentile in its null and the fraction of proteins
above the null 95th percentile; same for S at the best centre within 4 Å of a labelled site.

## Deviations

Any departure is recorded with time in the final reply (no separate file).
