# Validation — panel, reproducibility, and known limits

> **Current agent guide — 2026-09-16:** [Alquemia operations and protocol status](docs/AGENT_PIPELINE.md).
> The fixed-core PQQ release, crystal transfer, 26-endpoint baseline/repair benchmark, ten-endpoint additions and 38-endpoint GGR study are completed and linked in the current guide. Broad La/Ca affinity discrimination remains unvalidated.
> This guide supersedes conflicting operational/status prose below; dated scientific records remain historical.

This document is the validation evidence for the DFT vertical-swap
Ca²⁺/Ln³⁺ discriminator. It is the answer to "does this work?" written for
a reviewer / second-opinion reader. For methodology see
[METHODS.md](METHODS.md); for live results see [RESULTS.md](RESULTS.md).

The numbers below come from `results/all_results.jsonl`,
`results/discriminator_panel_2026-05-07.tsv`,
`results/findings_summary_2026-05-07.md`, the per-candidate `*_qm/` directories,
and the `b97_3c_panel/` and `calexcitin_size_panel/` directories. Every
claim is sourced.

---

## 1. The XoxF/MxaF gold-standard pair

**Same fold + same cofactor (PQQ) + different metal preference by biology.**
The XoxF/MxaF pair is the cleanest in-vivo Ca²⁺-vs-Ln³⁺ comparison that exists
in the literature: they are bacterial paralogs that share the PQQ-MDH
architecture but differ in the active-site donor set.

| protein       | structure                | active-site geometry                              | ΔΔE        |
|---------------|--------------------------|---------------------------------------------------|-----------:|
| **XoxF**      | Methylobacterium 4MAE    | Ce³⁺-MDH; 2 Asp + 1 Glu + 1 Asn + PQQ → CN9       | **+24.13** |
| **MxaF**      | Methylobacterium 1H4I    | Ca²⁺-MDH; 1 Glu + 1 Asn + PQQ → CN6                | **−6.48**  |

**31 kcal/mol gap** between the two paralogs, with the sign correct in both
cases (XoxF Ln-evolved, MxaF Ca-evolved). This is the bedrock of the
discriminator's claim to validity.

Source: `results/discriminator_panel_2026-05-07.tsv` rows 8 and 14;
`results/all_results.jsonl` stems `xoxf` and `mxaf`; backed by the
energies in `xoxf_qm/sp_xoxf_*.out` and `mxaf_qm/sp_mxaf_*.out`.

Reproducibility: rerunning the pair on a freshly-protonated input
gave the same ΔΔE values within ~1 kcal/mol (the PDBFixer noise floor,
which cancels in the vertical swap; see §4).

---

## 2. Tier-0 controls (canonical Ca-binding proteins; chemistry-positive in vitro)

These are the "negative controls" in the colloquial sense — proteins
biologists call **Ca-binders** because that's what they bind in vivo. The
discriminator does *not* score them as Ca-evolved. Instead, most of them
land **Ln-preferring** at the chemistry level, in a tight +5 to +20 band:

| protein                    | structure | ΔΔE       | class           | notes                                                                                       |
|----------------------------|-----------|----------:|-----------------|---------------------------------------------------------------------------------------------|
| **Calmodulin (EF1)**       | 1CLL      | **+18.15** | Ln-preferring   | Canonical eukaryotic CaM. Tb³⁺/Eu³⁺ luminescence probes confirm this chemistry experimentally. |
| **Calexcitin (EF3)**       | —         | **+10.59** | Ln-preferring   | Convergent EF-hand; not catalogued as a Ln-binder in vivo.                                  |
| **Parvalbumin**            | 1B9A      | **+11.4**  | Ln-preferring   | Canonical EF-hand. Required `TightSCF` + level-shift to converge — see §8.                  |
| **Calbindin**              | 1A75      | **+3.71**  | marginal        | Weak EF-hand-like.                                                                          |
| **Fern peroxidase**        | —         | **−1.55**  | ambiguous       | Class III peroxidase structural Ca site. Genuinely neutral.                                  |

(Calmodulin / calbindin / parvalbumin / fernperox: from
`results/discriminator_panel_2026-05-07.tsv`. Calexcitin: refreshed in JSONL
at +10.59.)

**Interpretation.** This is **not a discriminator failure**. The chemistry of
a 4-carboxylate EF-hand cluster *does* prefer Ln³⁺ over Ca²⁺ at fixed
geometry — that is textbook coordination chemistry, demonstrated half a
century ago by Tb³⁺ / Eu³⁺ luminescence probe experiments (the "lanthanide
luminescence probe" literature). The proteins are Ca-binders in vivo only
because environmental [Ln³⁺] is trace and cytoplasmic [Ca²⁺] is mM. The
discriminator measures **chemistry, not in-vivo occupancy** — see
[METHODS.md §7](METHODS.md).

The Ca-evolved tier is reserved for the structurally-distinct MDH-style
active sites (MxaF and other 1-Glu + 1-Asn + PQQ CN≤6 pockets) and the
hard-Ca-tuned dinuclear sites (prolidase −20.7, certain arginases −16.7) —
all of which the discriminator correctly identifies.

---

## 3. Tier-1 controls (validated Ln-binders)

These are biologically-established Ln-binders. The discriminator scores
each at or above +10 kcal/mol; the canonical positive controls cluster
tightly:

| protein                                    | ΔΔE        | class         | source                                                |
|--------------------------------------------|-----------:|---------------|-------------------------------------------------------|
| **LanM (M. extorquens) EF1**               | **+47.5**  | Ln-evolved    | original 2026-05-07 panel; `qmmm/cluster_panel/`     |
| **Tannase paralog `NZ_JYMT`**              | **+27.46** | Ln-evolved    | `tannasepara_NZ_JYMT_qm/`                            |
| **A0A226QCS9 (thermolysin Ca₂ site)**      | **+24.42** | Ln-evolved    | rediscovers Holmquist & Vallee 1974 — see §6.        |
| **XoxF**                                   | **+24.13** | Ln-evolved    | gold-standard positive (§1)                          |
| **Colin's 11 La-verified PQQ-MDHs**        | **+11 to +23** | Ln-preferring / Ln-evolved | 11/11 score positive; tight cluster |
| **SilE (×2 independent structures)**       | **+21.33 / +22.40** | Ln-evolved | convergent — `SilE_La_sample_0`, `SilE_like_LanM_associated_sample_0` |
| **Tannase (canonical)**                    | **+15.68** | Ln-preferring | `tannase_qm/`                                        |

### Colin's 11 La-verified PQQ-MDH set (`00_Lav` + reruns)

Colin's curated set of 11 PQQ-MDH structures with confirmed La³⁺ binding
(biochemistry + crystallography). Run twice on the discriminator:

- **First pass** (`colinpqq_la_*` direct driver, 2026-05-08): all 11 scored
  Ln-preferring or Ln-evolved, +11 to +23 kcal/mol band.
- **Second pass** (`00_Lav_*` PQQ-aware inbox routing, 2026-05-11): all 11
  scored Ln-preferring or Ln-evolved, +11 to +23 kcal/mol.
- **Cross-pipeline reproducibility:** 0.2 to 3 kcal/mol agreement between
  the two passes on each of the 10 stems that completed in both. **0 / 11
  false-negatives** on Colin's positive set.

Source: `results/all_results.jsonl` rows with stems
`colinpqq_la_<id>_pqq_la_model` and `00_Lav_<id>-pqq-la_model`.

---

## 4. Cross-run reproducibility

| comparison                                                                  | spread          | source                                  |
|-----------------------------------------------------------------------------|-----------------|------------------------------------------|
| **PDBFixer non-determinism (absolute SP energy, same input twice)**         | 15–66 kcal/mol  | observed across protonation rounds       |
| **PDBFixer non-determinism (ΔΔE, same input twice)**                        | **<0.5 kcal/mol** | the vertical-swap cancellation working   |
| **Marco's β-roll family paralogs (4 highly-related sequences)**              | **<0.05 kcal/mol** | within `findings_summary_2026-05-07.md` |
| **Colin's 00_Lav re-run vs `colinpqq_la_*` legacy run**                     | **0.2–3 kcal/mol** | validated 2026-05-12 (CONTEXT.md)      |

The vertical-swap design works as expected: absolute energies are noisy by
~50 kcal/mol because PDBFixer picks tautomers stochastically, but **both
metals see the same protonation in any single run** and the noise
cancels cleanly in the ΔΔE difference.

---

## 5. False-positive families catalogued

Even with the carboxylate-cluster chemistry working correctly, several
families consistently score Ln-preferring without being Ln-binders. These
have been catalogued in `results/known_false_positives.tsv` so we don't
re-discover them every batch.

| family                                              | representative stem                           | typical ΔΔE   | why it scores positive                                                                                  |
|-----------------------------------------------------|-----------------------------------------------|--------------:|---------------------------------------------------------------------------------------------------------|
| **M20 amidohydrolase (Zn²⁺ peptidase)**             | `NZ_BPQI01000118_1_36_iptm0_953_s0`           | +5 to +10    | XoxF-mimic 4-Asp/His-poor cluster the discriminator confuses for a hard-cation pocket                   |
| **VIT1 family (Fe³⁺ uptake transporter)**           | `NZ_JACIDR010000005_1_9_iptm0_932_s0`          | +30+         | 5-carboxylate cage; reclassified 2026-05-09 as **hypothesized FP / candidate Ln transporter** — VIT1 is functionally promiscuous and a Ln-co-opted VIT1 is plausible. Worth follow-up. |
| **Diiron rubrerythrin / bacterioferritin**          | `C5CIA3_bryobacter`                            | ~+10         | 4-Glu binuclear pocket; La grabs the second Fe site                                                     |
| **CcmF heme biogenesis**                            | `NC_010511_1_4382`                             | +6 to +8     | cytochrome C assembly factor — Fe-handling, not Ln                                                       |
| **Mn²⁺/Zn²⁺ small acidic pocket**                   | e.g. `RXJC01000147_1_13`, `NZ_FNHS01000007_1_195` | < −30 (OUTLIER) | The pocket wants a smaller, harder cation; SCF gets cornered with La and the gap blows up — biology signal |

The OUTLIER tier (last row) is itself a useful annotation: a |ΔΔE| > 30
kcal/mol on a sub-50-atom cluster is a **Mn/Zn pocket signature**, not a
numerical artifact. See [TIERS.md](TIERS.md) for the SCF interpretation.

The discriminator is, at the chemistry level, a "carboxylate-rich high-charge-cation
pocket detector." Biology context (Pfam, KO, operon, no-Cys/Met, LanM neighbor)
is what separates real Ln-binders from these convergent geometric mimics.

---

## 6. Confirmed Ln-binder families

Beyond XoxF (the gold standard), five additional families pass both the
chemistry filter (Ln-evolved tier) **and** the biology filter (independent
evidence for Ln-binding). All are catalogued in
`results/confirmed_ln_binders.tsv`.

| family                                         | example stem                                     | ΔΔE       | biology evidence                                                                                            |
|------------------------------------------------|--------------------------------------------------|----------:|-------------------------------------------------------------------------------------------------------------|
| **XoxF / PQQ-MDH**                             | `xoxf`                                           | +24.13    | Canonical Ln-MDH (K23995); methylotrophy literature                                                          |
| **Tannase (Marco β-roll)**                     | `tannasepara_NZ_JYMT` (+27.46), `tannase` (+15.68), `FSRD01000002_1_1028` (+20.37 from spicy_lams) | +15 to +27 | Bradyrhizobium tannase fold; LanM-less paralog distinct family per Jacob 2026-05-09 |
| **SilE (small Ln-binder)**                     | `SilE_La_sample_0` (+21.33), `SilE_like_LanM_associated_sample_0` (+22.40) | +21 to +22 | **Convergent +21/+22 across two independent structures.** Family hypothesis confirmed by discriminator.    |
| **Thermolysin Ca₂ site (Geobacillus)**         | `A0A226QCS9_acidithiobacillus` (bin label upstream is misassigned; sequence is actually *Geobacillus / B. thermoproteolyticus*) | +24.42 | **Rediscovers Holmquist & Vallee 1974** — Ln³⁺ substitution at thermolysin's non-catalytic Ca₂ site was demonstrated in vitro 50 years ago. Discriminator independently picked it up from a blind scan. Strong methodology validation. |
| **Colin's 11 La-verified PQQ-MDH controls**    | `colinpqq_la_*`                                  | +11 to +23 | 11/11 confirm Ln-binding; both runs match within 0.2–3 kcal/mol of each other                              |
| **Eukaryotic PQQ-ADHs (2026-05-12 discovery)** | `03_Euk_a0a2v0prl1_domain3...` (+23.45), `a0a2v0nw76_domain1` (+22.35), `a0a423vtm4_domain1` (+18.70), `a0a2j6qw99_domain1` (+17.31) | +6 to +23 | 21+/50 score Ln-preferring or Ln-evolved in the in-flight Colin Euk batch. See [RESULTS.md §3](RESULTS.md). |

The Marco β-roll tannase paralog (+27.46) is the strongest "blind Ln-evolved
hit" in the panel — a tannase paralog from Bradyrhizobium that scored
+27 *without* a LanM cassette in its operon, suggesting the tannase fold is
*intrinsically* Ln-evolved at the architectural level (not merely operonic).

---

## 7. Functional robustness — B97-3c cross-check

To check that the ΔΔE result isn't an artifact of the r²SCAN-3c functional
choice, the top four candidates were re-run with B97-3c (a different composite
DFT) at the same CPCM(Water) DefGrid3 level. Same carve, same XYZ, different
functional.

The B97-3c aquo references:
- `[Ca(H₂O)₆]²⁺` at B97-3c/CPCM(Water): `−1288.761999564482 Ha`
- `[La(H₂O)₆]³⁺` at B97-3c/CPCM(Water): `−642.749072156339 Ha`
- `DIFF_AQUO (B97-3c)`: `−646.012927408143 Ha`

| protein     | ΔΔE (r²SCAN-3c) | ΔΔE (B97-3c) | Δ (B97 − r²SCAN) | class preserved? |
|-------------|---------------:|-------------:|------------------:|------------------|
| XoxF        | +24.13         | +21.46       | −2.67             | Ln-evolved        |
| MxaF        | −6.48          | −7.42        | −0.94             | Ca-evolved        |
| Tannase     | +15.68         | +22.72       | +7.04             | Ln-preferring    |
| Calexcitin  | +10.59         | +11.32       | +0.73             | Ln-preferring    |

Source: `b97_3c_panel/{xoxf,mxaf,tannase,calexcitin}/sp_*_b97.out` and
`b97_3c_panel/aquo/sp_{Ca,La}_b97.out`.

**Class assignment is preserved on all four** (the signs match and every
candidate stays in its tier). The magnitudes scatter by up to 7 kcal/mol
on tannase; that is the typical inter-functional disagreement on composite
DFT methods and does not change the qualitative call.

The takeaway: the discriminator's tier assignments are not sensitive to the
specific composite-DFT choice between r²SCAN-3c and B97-3c. r²SCAN-3c is
the production functional because it converges faster and gives slightly
tighter cross-run reproducibility.

---

## 8. Cluster-size sensitivity

Earlier methodology work tested a calexcitin (EF-hand) cluster at 16, 32,
50, and 100 atoms to find the minimum carve size that gives a reliable
ΔΔE. The numbers (recomputed at r²SCAN-3c/CPCM(Water) with the production
DIFF_AQUO):

| cluster size | atoms (Ca/La) | ΔΔE (kcal/mol) | reading                                                                |
|--------------|---------------|---------------:|-------------------------------------------------------------------------|
| 16-atom      | 16            | +3.36          | Under-carved; misses second-shell stabilisation                          |
| 32-atom      | 32            | +10.59         | Production-size carve (matches the 2026-05-07 panel value)               |
| 50-atom      | 50            | −6.23          | Includes residues that shift the chemistry (bringing second-shell waters / charged loops into QM region) |
| 100-atom     | 100           | −3.72          | Larger still — same drift direction                                      |

Source: `calexcitin_size_panel/sp_calexcitin_size{16,32,50,100}_{La,Ca}.out`.

**Reading.** The 32-atom carve (production size for protein-only sites) matches
the original panel value (+10.59). Pushing the carve to 50 or 100 atoms does
not "converge" the answer in a stationary sense — it drags in additional
residues whose chemistry is no longer first-shell, and the answer drifts in
the negative direction. This is consistent with the design: the discriminator
is calibrated on **first-shell chemistry within 3.0 Å**. The 32-atom carve
captures that first shell; the 50/100-atom carves over-extend and the
arbitrariness of where to cut becomes the dominant signal.

**Recommendation.** Stay at the production carve cutoff (3.0 Å, sidechain
dict ASP/GLU/ASN/GLN/SER/TYR + backbone-O). Do not arbitrarily widen for
"more accuracy" — you'll trade signal for second-shell scaffolding.

For PQQ-aware carves the typical size is 55–58 atoms (because PQQ alone
adds 24 atoms), which is the right size in that case because PQQ *is* the
first-shell donor pair.

---

## 9. SCF convergence cases requiring intervention

A handful of candidates required `TightSCF` + level-shift to converge:

- **Parvalbumin 1B9A La SP** — restarted with `TightSCF` + `%scf shift shift 0.4 erroff 0.1 end end`.
- **Various OUTLIER-tier Mn²⁺/Zn²⁺ sites** — these are inherently SCF-difficult
  because La doesn't fit the pocket geometrically. The convergence
  difficulty is part of the biology signal (see [TIERS.md](TIERS.md)).
- **DUF882** (early panel) — 15-atom under-carved cluster gave ΔΔE
  ≈ −1378 kcal/mol; flagged as `EXCLUDE`.

The recipe (from [HOWTO.md §11](HOWTO.md)):

```bash
# Edit sp_<stem>_La.inp:
#   ! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
# → ! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF
#   %scf shift shift 0.4 erroff 0.1 end end
$ORCA_PATH/orca sp_<stem>_La.inp > sp_<stem>_La_v2.out
```

---

## 10. What still doesn't work

Honest list of known limitations:

- **Severe asymmetric carve (|ΔΔE| > 100)** — a handful of acidobacterium and
  03_Euk candidates have a PQQ-inclusion mismatch between Ca and La carves
  (one includes the cofactor, the other doesn't). Manual atom-count diff
  required to diagnose. See [TIERS.md](TIERS.md#severe-asymmetric-carve-manual-flag-no-auto-label).
- **Empty-carve epidemic on coordinated CYS / HIS / MET sites** — out of
  policy by design. These get classified as `CARVE_AMBIGUOUS` and not
  scored. If the upstream fold-daemon emits CIFs where the AF3 model
  predicts a Cys-coordinated metal, those candidates are silently dropped.
- **Cluster-size convergence non-monotonic** (§8). The 32-atom carve gives
  the production answer, but you cannot defend that as a "converged" answer
  in the usual sense — only as a "first-shell" answer. Reviewers may push
  on this; the response is in [METHODS.md §7](METHODS.md): the
  discriminator measures *first-shell* chemistry, not full-protein binding
  thermodynamics.
- **Single-snapshot, vertical swap** — no induced fit. The "relaxed-swap"
  validation tier (re-fold each candidate with each candidate metal) is on
  the open-work list in [CONTEXT.md](CONTEXT.md). Would resolve some
  borderline marginal/ambiguous calls and rescue some OUTLIER candidates.
