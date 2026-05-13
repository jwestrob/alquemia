# alchemical_bvs — vertical-swap DFT discriminator for Ca²⁺ vs Ln³⁺ pocket preference

![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![ORCA](https://img.shields.io/badge/ORCA-6.1.1-orange)
![License](https://img.shields.io/badge/license-TBD-lightgrey)

A DFT cluster-level discriminator for asking, of a given candidate metal-binding
pocket: **does this site prefer Ca²⁺ or Ln³⁺ at the coordination-chemistry level?**

For each candidate, the pipeline carves a QM cluster (metal + first-shell residues
+ cofactor if present), runs four r²SCAN-3c / CPCM(Water) single points
(La form, Ca form, apo, bulk water), and computes a single number ΔΔE in
kcal/mol. The aquo reference (DIFF_AQUO = −646.064556 Ha) cancels absolute-energy
noise — protonation non-determinism, ECP-vs-AE offset, basis incompleteness —
leaving a clean coordination-chemistry signal. See [METHODS.md](METHODS.md).

---

## Headline result — the XoxF/MxaF gold standard

The cleanest in-vivo Ca²⁺/Ln³⁺ comparison in the literature is the bacterial
paralog pair XoxF (Ln³⁺-MDH) and MxaF (Ca²⁺-MDH): same fold, same PQQ cofactor,
different active-site donors, different metal in vivo. The discriminator
splits them cleanly:

- **XoxF** (4MAE): ΔΔE = **+24.13** kcal/mol → **Ln-evolved**
- **MxaF** (1H4I): ΔΔE = **−6.48** kcal/mol → **Ca-evolved**

A **31-kcal/mol gap** between two paralogs, with the sign correct in both
cases. This is the methodology gold standard; reproducibility on rerun is
within ~1 kcal/mol (the PDBFixer noise floor). See
[VALIDATION.md §1](VALIDATION.md).

---

## Current scale (2026-05-12)

| tier                                          | count |
|-----------------------------------------------|------:|
| **Ln-evolved** (ΔΔE ≥ +20)                    | **17** |
| **Ln-preferring** (+5..+20)                   | **278** |
| marginal (0..+5)                              | 216 |
| ambiguous (−5..0)                             | 201 |
| Ca-evolved (< −5)                             | 418 |
| OUTLIER (|ΔΔE| > 30, n_atoms < 50)            | 96 |
| EXCLUDE (under-carved, n_atoms < 15)          | 90 |
| pending                                       | 137 |
| **TOTAL**                                     | **1,453** |

Full top-hit tables and the eukaryotic PQQ-ADH discovery story are in
[RESULTS.md](RESULTS.md); class definitions in [TIERS.md](TIERS.md).

---

## Quick start

You have a CIF with La placed in a candidate metal-binding pocket.

```bash
# 1. Drop the CIF in the inbox
cp my_protein.cif inbox/

# 2. Trigger the pipeline (or let the 45-min watcher do it)
bash scripts/process_inbox.sh

# 3. After SLURM finishes (~20–60 min per candidate), refresh the manifest
python scripts/update_results_jsonl.py .
```

Full how-to: [HOWTO.md](HOWTO.md). Pipeline internals: [PIPELINE.md](PIPELINE.md).
Conda env spec: [environment.yml](environment.yml).

---

## Confirmed Ln-binder families

Beyond XoxF, the discriminator has been validated against — or independently
discovered — the following families. Full catalogue in
`results/confirmed_ln_binders.tsv`; details in [VALIDATION.md §6](VALIDATION.md):

- **XoxF / PQQ-MDH** — gold standard (+24.13)
- **Tannase (Marco β-roll)** — `tannasepara_NZ_JYMT` (+27.46), canonical `tannase` (+15.68)
- **SilE family** — `SilE_La_sample_0` (+21.33), `SilE_like_LanM_associated_sample_0` (+22.40); convergent
- **Thermolysin Ca₂ site** — `A0A226QCS9` (+24.42); independently rediscovers Holmquist & Vallee 1974
- **Colin's 11 La-verified PQQ-MDH controls** — `colinpqq_la_*`, 11/11 score +11 to +23
- **Eukaryotic PQQ-ADHs (2026-05-12 discovery)** — 21+/50 of Colin's `03_Euk` batch score Ln-class; top hits `A0A2V0PRL1` (+23.45), `A0A2V0NW76` (+22.35), `A0A423VTM4` (+18.70)
- **AFDB Acidobacterium/Actinobacteriota clade** — `A0A7V5CTA7` (+17.87), a porin-family Ln-preferring hit
- **`small_unannotated_candidate_marco`** — 31-protein novel family in +7 to +13 band

---

## Documentation map

| document                                       | what's there                                                |
|------------------------------------------------|-------------------------------------------------------------|
| [README.md](README.md)                         | this file — top-level entry point                           |
| [METHODS.md](METHODS.md)                       | theory and chemistry — ΔΔE equation, DIFF_AQUO, carve rules, PQQ, electronic structure, class thresholds, what the discriminator does not measure |
| [PIPELINE.md](PIPELINE.md)                     | runtime architecture — data flow, per-script reference, idempotency, failure modes, SLURM integration, state files |
| [TIERS.md](TIERS.md)                           | classification reference — Ln-evolved, Ln-preferring, marginal, ambiguous, Ca-evolved, OUTLIER, EXCLUDE, SOLVENT_EXCLUDE, CARVE_AMBIGUOUS, SEVERE asymmetric |
| [VALIDATION.md](VALIDATION.md)                 | validation panel + reproducibility — XoxF/MxaF, Tier-0/Tier-1 controls, B97-3c cross-check, cluster-size sensitivity, known false-positive families, limitations |
| [RESULTS.md](RESULTS.md)                       | current findings — class breakdown, top hits, Euk PQQ-ADH discovery, cross-domain pattern, wet-lab handoff list |
| [HOWTO.md](HOWTO.md)                           | hands-on pipeline operation — quick start, walkthroughs, troubleshooting, ORCA settings reference |
| [CONTEXT.md](CONTEXT.md)                       | current project state — counts, families, open work, methodology revisions log |
| [CHANGELOG.md](CHANGELOG.md)                   | methodology revisions, chronologically                      |
| [CONTRIBUTING.md](CONTRIBUTING.md)             | how to extend, submit from other accounts, add a sidechain  |
| [scripts/README.md](scripts/README.md)         | per-script compact reference                                |
| [examples/walkthrough.md](examples/walkthrough.md) | fully worked single-protein example                       |
| [CITATION.cff](CITATION.cff)                   | citation file                                               |
| [environment.yml](environment.yml)             | conda environment spec for the primary `lanm_qmmm` env      |
| [LICENSE](LICENSE)                             | TBD — choose before public release                          |

Legacy docs at the repo root: [HANDOFF.md](HANDOFF.md) (pre-discriminator,
2026-05-03) and [REPORT.md](REPORT.md) (8DQ2 BVS pre-check). Both predate
the methodology pivot to a vertical Ca/Ln swap.

---

## License

See [LICENSE](LICENSE). Currently a placeholder; choose an OSS license (MIT,
BSD-3-Clause, or Apache-2.0 are reasonable defaults for a methodology
pipeline) before public release.

---

## Citation

See [CITATION.cff](CITATION.cff). Until a manuscript lands the suggested
form is:

> Westrob J. et al. (in prep). *A vertical-swap DFT discriminator for Ca²⁺
> vs Ln³⁺ pocket preference in metalloproteins.* Target: JCTC.

---

## Affiliation

Banfield Lab, UC Berkeley.
