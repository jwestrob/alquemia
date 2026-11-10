# La/Ca imprint audit — 9 November 2026

**Result.** The general "imprint" hypothesis fails its frozen test. In native-Ca sites, a bound
lanthanide leaves the rule's call and the site geometry largely unchanged. La itself appears in one
native-Ca pair, and that pair does show the imprint. The La side of the crystal benchmark has too few
families for any difference between discriminators to be detectable.

Written by the lead agent from the analysis subagent's reply and the saved tables in `results/`
(see DEVIATIONS.md). `PLAN.md` sha256 `ed314201…0153` matches `PLAN.sha256`.

## Hypothesis under test

"Imprint": the discriminating features (coordination number, mean metal–oxygen distance, chelating
carboxylates) record which ion is bound, so the rule and the energetic discriminator answer "which
ion is in this structure". Alternative, "preference": the features are a property of the site that
persists whichever ion is bound.

## Part 1: what each benchmark source was scored on

- **Geometric rule, crystals.** 22 La sites in 5 deposits from 3 families (LanM, LanP, XoxF) and 125
  Ca sites in 34 deposits. Every site was scored on a crystal holding its native ion: 0
  counterfactuals.
- **Nikasha "25 canonical".** 25 PQQ-MDH proteins, 11 La-associated and 14 Ca-associated, one AF3
  fold each, with the metal set to La in all 25. The 14 Ca-labelled proteins were therefore scored on
  geometry made with La; all 14 were called correctly.
- **The "125" in the Nikasha ledger** is those 25 proteins × 5 Ca-conditioned AF3 folds, a different
  set from the 125 Ca crystal sites. It puts 55 La-labelled sources on Ca-made geometry:
  accommodated 51 correct / 1 wrong / 3 inconclusive; static 49 / 1 / 5.
- **Total.** 65 of 69 counterfactual-ion sources were called correctly.
- **Crystal controls.** 1H4I (Ca), 1KB0 (Ca), 4MAE (Ce), each with its native ion.

Inference: on the PQQ-MDH benchmark the energetic calls follow the site encoded in the sequence, the
D+2-Asp / core-charge confound the project already documents.

Open points: no document read measures how much AF3 ion conditioning changes site geometry, so these
sources may be counterfactual in label only. A pre-registered α-lactalbumin 2×2 crossing (6IP9 La /
1F6S Ca) exists; its result was not found. No energetic run on the 147 crystal sites was found, so
"performs about the same" there is UNVERIFIED.

## Part 2: paired same-site test (PDB; cluster login node)

Survey: 410 lanthanide X-ray entries at ≤ 2.5 Å in 166 clusters at 95% identity; 46 clusters have a Ca
partner (255 files, no sampling; 8 clusters reached the frozen cap of 8 entries per side). 333 site
pairs: median metal offset 0.24 Å, median superposition RMSD 0.14 Å, minimum chain identity 0.969.
194 of the 333 lanthanide ions have occupancy below 1.

**Primary stratum** (curated native-Ca site; Ca deposit with ≥ 4 protein donors; 12 proteins, all
Swiss-Prot reviewed; bootstrap over proteins):

| Measure | Value [95% CI] |
|---|---|
| p, lanthanide deposits called La | 0.175 |
| q, Ca deposits called La | 0.083 |
| p − q | 0.092 [−0.083, 0.300] |
| Frozen reading | **Imprint refuted** (criterion p − q ≤ 0.10 with ≥ 10 proteins) |
| ΔCN, lanthanide − Ca | −0.28 [−1.05, 0.43] |
| Δ mean M–O | +0.009 Å [−0.059, 0.091] |
| Δ chelating carboxylates | +0.12 [0.00, 0.33] |

The interval on p − q reaches 0.30. Unweighted by protein, p = 0.296 and q = 0.179.

**Split by lanthanide (added after the primary result; descriptive):**

| Group | Proteins | p − q [95% CI] | ΔCN |
|---|---|---|---|
| Heavy (Gd, Tb, Ho, Yb, Lu) | 8 | −0.050 [−0.188, 0.037] | −0.66 |
| Light (La, Pr, Sm, Eu) | 4 | 0.375 [0.000, 0.750], q = 0 | +0.47 [−0.48, 1.70] |

- La occurs in one native-Ca pair, bovine α-lactalbumin: 6IP9 La203 has CN 9, mean La–O 2.517 Å and
  is called La; the three Ca deposits (1F6S, 1HFZ, 9OZW) have CN 6–7, mean Ca–O 2.389–2.450 Å and are
  called Ca.
- Proteinase K: Pr in 5AVK and 6K2X has CN 9 (called La); Pr in 5WRC and Eu in 9KW1 have CN 8 (called
  Ca); all 7 Ca deposits have CN 8.

Inference: ions near or below the Ca²⁺ radius keep the Ca shell; La and Pr can add a donor and flip
the rule. The general imprint is refuted at n = 12. For La the evidence is one protein, which shows
the imprint.

**Other strata.** Sites with no curated native-Ca annotation (UNVERIFIED, 6 proteins): p − q 0.368
[0.125, 0.667], ΔCN +0.83. Reverse pairs (curated native-lanthanide site with a Ca-bound deposit): 0.
The frozen keyword rule excluded 6 pairs from 3 proteins whose source organism is "synthetic
construct".

## Part 3: power of the crystal benchmark

α = 0.05 two-sided, exact tests, 80% power. The paired column is a best case (every disagreement
favours one method).

| Unit | Smallest detectable accuracy difference, paired | Unpaired |
|---|---|---|
| La family (n = 3) | none at any effect size (best p = 0.25) | none (best p = 0.10) |
| La site, clustering ignored (n = 22) | 0.34 | 0.36 |
| Ca site (n = 125) | 0.063 | 0.10 |
| Ca deposit (n = 34) | 0.22 | 0.255 |

Any detectable difference between two discriminators has to come from the Ca calls.

## Limits

- Native-Ca status rests on UniProt annotations; primary citations were not read.
- Light-lanthanide occupancies are 0.70–1.0, so some shells may mix ions.
- Soaked and co-crystallised structures are pooled.
- The light-lanthanide split rests on 4 proteins and the La result on 1.

## Files

- `PLAN.md`, `PLAN.sha256`; `scripts/` (`rcsb.py`, `discover.py`, `download.py`, `pair_sites.py`,
  `stats.py`); `results/` (`pairs.tsv`, `per_protein.tsv`, `stats_summary.json`, `candidates.json`,
  `pair_log.json`, `stats.log`).
- `workspaces/laca_imprint_audit_20261109/`: 255 mmCIFs, API cache, logs.
- Reproduce: `discover.py`, `download.py`, then `pair_sites.py` under `conda_envs/lanm_qmmm` with
  `COD=../../cod_reference_geometry_2026-10-13/scripts`, then `stats.py`.
