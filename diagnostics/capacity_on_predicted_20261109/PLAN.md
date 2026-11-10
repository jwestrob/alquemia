# Reachable acidic capacity on predicted structures — frozen plan (9 November 2026)

Question: does the crystal capacity rule (≥ 4 reachable Asp/Glu at a centre; detector v0 and baseline B2
of `pdb_ln_site_benchmark_20261109`) hold on single-sequence ESMFold2 models that carry no metal?

Frozen before any fold was requested or scored. Done before freezing: sample draws (below), the frozen
`sites.py` run on 6OC5 and 9C90 (labelling only), sequence extraction. Nothing has been scored.

## 0. Frozen code and environment

- `detect.py`, `common.py`, `units.py`, `sites.py` copied unchanged from the benchmark's `scripts/` into
  `workspaces/capacity_on_predicted_20261109/frozen/` (sha256 identical; recorded in the reply).
- `first_hit`, `residue_reach`, `clash` are copied verbatim from the benchmark's `score.py` (that file
  parses argv at import and cannot be imported).
- Interpreter: `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -I`, biotite login
  node, `OMP_NUM_THREADS=1`, one process (at most 3 threads in total).
- Folding: Mac, `Sharur/scripts/fold_biohub.py OUT FASTA` with defaults (endpoint `fold`, model
  `esmfold2-fast-2026-05`, single sequence, no MSA, no ligand). Budget ≤ 150 sequences, ≤ 1,500 aa each.
  pLDDT read from model CA B-factors (0–100 scale; checked against `summary.tsv`).

## 1. Samples (drawn before freezing; files in `results/`)

- **Part 2 benchmark sample** (`select_bench.py` → `bench_sample.tsv/.faa`): stripped-tier units with
  ≥ 1 single-chain site with ≥ 3 protein donors and an entity sequence ≤ 1,000 residues (178 units,
  69 c95, 64 c30). c30 clusters sorted, shuffled with `random.Random(42)`, first 50 taken; one unit per
  cluster by `rng.choice` over its eligible units sorted by uid (same RNG). Folded sequence = the full
  entity sequence (`_entity_poly_seq`, MSE→M). 50 proteins, 70 sites, 15,001 aa.
- **Part 4 reference proteins** (`ref_sample.tsv`, `ref.faa`): table rows whose PDB ids hold a
  lanthanide in an X-ray structure: 4MAE (I0JWN7), 6OC5 (C5B120), 6DAM (A0A3F2YLY8), 7O6Z (A0A5E6MBX8),
  6ZCV (Q88JH0), 8FNS (C5B164), 8DQ2 (A0A4V2JDD3), 9C90 (C5B159, Ho), 9VY4 (Q1H2G7), 9PI2 (no accession;
  entity sequence). 6MI5 is NMR without Ln; 9C8W has no Ln; both left out. The table gives no mature
  sequence, so the `ln_proteins.faa` sequence (precursor) is folded. Sites: every Ln site attributed
  to chain A with ≥ 2 protein donors (benchmark `ln_sites.tsv`; for 6OC5 and 9C90 the frozen `sites.py`
  run in `workspaces/.../ref_label/`). 10 sequences.
- **Part 5 candidates** (`cand_sample.json`, `cand.faa`): from `acidic_loop_families.tsv`, groups
  unannotated (family `pfam` = none; 37), comparators (family `pfam` contains FG-GAP, TSP_3, TSP3_bac
  or HemolysinCabind; 11, no HemolysinCabind family is present), TPR (families with TPR in `pfam`;
  sorted, `random.Random(42)` shuffle, first 10; exactly 10 exist). Representative region = highest
  `S`, then `region_id`. Sequence: Talea `record_id` → `sequence_sha256` → `protein_sequence_map` →
  smallest `protein_id` → Sharur `proteins.sequence` (sha256 verified, 58/58). Proteins > 1,200 aa:
  fold residues [reg_start − 50, reg_end + 50] (1-based, clipped), truncated to 1,200 from its start
  (16 proteins). Period = Talea `period` when Talea's registered start equals the region start, else
  the region's `unit_len`. 58 sequences, 27,667 aa.
- Total 118 sequences.

## 2. Predicted tier of the benchmark (the test)

Mapping: crystal chain residues (`common.chain_residues`) aligned to model residues with
`common.align_map`. **Site transfer** to the model as in `units.py`: crystal Cα within 12 Å of the site
with a mapped model residue, ≥ 8 pairs, `common.kabsch`, RMSD recorded. A site that fails transfer is
removed from every tier (same site set throughout); failures reported.

Tiers for the same sites: (a) stripped crystal chain, re-scored here with the frozen code (ranks must
equal the benchmark's `site_ranks.tsv`; mismatches reported); (b) partner: the benchmark's
`site_ranks.tsv` partner rows for these site ids (subset of proteins); (c) predicted model, full
folded sequence; (c′, secondary) model restricted to residues aligned to crystal-observed residues.

Per site: hit@k (k = 1, 3, 5; centre within 4 Å) for v0 and B2; `an_true` = reachable Asp/Glu at
the true position (`residue_reach`, no clash filter); ≥ 4 indicator. Aggregation site → unit → entry →
protein (c95) mean; bootstrap over c30 clusters, 2,000 resamples, `numpy.random.default_rng(42)`,
95% percentile (the benchmark's `boot`). Ratio F4(c)/F4(a) bootstrapped on the same resamples.

**Frozen reading**: "the count rule transfers to predicted structures" iff v0 recall@5 on (c) ≥ 0.45
AND F4(c) ≥ 0.8 × F4(a), where F4 = protein-mean fraction of true positions with an_true ≥ 4, on the
same sites; point estimates decide, intervals reported. Otherwise "does not transfer". B2 co-reported.

pLDDT stratification of (c): site pLDDT = mean model CA pLDDT over model residues mapped from crystal
same-chain residues with any O/N ≤ 3.0 Å of the site; strata < 70 and ≥ 70; recall@k and F4 per
stratum with bootstrap.

Descriptive: capacity curve on models (fraction of v0 centres within 4 Å of a transferred
any-Ln position, by n_anion bin) beside the same on the crystal chains (any Ln, any image).

## 3. False-positive burden

For t = 4, 5: candidate points with n_anion ≥ t, ordered by (−n_anion, −S, tiebreak), greedy NMS
radius 4 Å (`detect.nms`) → merged centres. Off-site = no Ln atom (any element in LN_EL, any copy or
symmetry image, any donor count) within 8 Å. Crystal chains: gemmi NeighborSearch on the full crystal
with its cell. Models: each centre mapped to the crystal frame by local superposition (model Cα within
12 Å of the centre with a mapped crystal residue, ≥ 8 pairs) and then queried as for crystals;
unmappable centres (residues absent from the crystal) count as off-site and are reported separately.
Per protein: counts and counts per 100 residues (crystal: observed residues; model: folded residues;
c′: retained residues). Medians, totals, and the paired model − crystal difference.

## 4. Reference proteins

Per site: class, protein donors, an_true and v0/B2 rank on the crystal chain A and on the model
(transfer as in section 2), site pLDDT, transfer RMSD, faa-to-crystal sequence identity.

## 5. Candidates (descriptive)

Per model: mean pLDDT (all folded residues); merged centres (section 3) with n_anion ≥ 4 and ≥ 5,
per 100 folded residues. For each merged centre: contributing Asp/Glu = Asp/Glu with reach weight 1.0
at the centre point (`residue_reach`); position = median protein position (1-based) of contributors;
the centre is in-region if that lies in [reg_start, reg_end]. Unit index of a residue =
floor((pos − reg_start)/period); units per centre = distinct unit indices among in-region contributors.
Period match: in-region centres (same threshold) sorted by position; a centre is period-matched if
its spacing to the previous or next in-region centre is within max(2, 0.15·period) of period or
2·period. Report per group (unannotated / comparator / TPR) counts and medians only; list the ten
models with most period-matched ≥ 4 centres (ties: ≥ 4 centres per 100 residues, then record id),
all UNVERIFIED. No function or bound metal is assigned to any candidate.

## Deviations

Recorded with time in the final reply.
