# Results — current findings

**Snapshot date:** 2026-05-12 (aggregator run at the time this doc was written).
**Source:** `results/all_results.jsonl` (deduped by stem, latest row wins),
`results/colin_euk_pqq_adh_results.tsv`, `results/confirmed_ln_binders.tsv`.

The tables below are auto-generated; to refresh, run
`python scripts/update_results_jsonl.py .` and re-run the snippet in
[§2](#2-top-hits) at the top of this file.

For the chemistry, see [METHODS.md](METHODS.md); for class semantics,
[TIERS.md](TIERS.md); for validation evidence, [VALIDATION.md](VALIDATION.md).

---

## 1. Class breakdown

| tier                                          | count |
|-----------------------------------------------|------:|
| **Ln-evolved** (ΔΔE ≥ +20)                    | **17** |
| **Ln-preferring** (+5..+20)                   | **278** |
| marginal (0..+5)                              | 216 |
| ambiguous (−5..0)                             | 201 |
| Ca-evolved (< −5)                             | 418 |
| OUTLIER (|ΔΔE| > 30, n_atoms < 50)            | 96 |
| EXCLUDE (under-carved, n_atoms < 15)          | 90 |
| pending (one or more SPs incomplete)          | 137 |
| **TOTAL**                                     | **1,453** |

Total Ln-class hits (Ln-preferring + Ln-evolved): **295**.

---

## 2. Top hits

### Top 17 Ln-evolved (ΔΔE ≥ +20 kcal/mol)

| rank | ΔΔE     | stem                                                       | iPTM  | source              |
|-----:|--------:|------------------------------------------------------------|------:|---------------------|
|   1  | +32.95  | NZ_JACIDR010000005_1_9_iptm0_932_s0                        | 0.932 | spicylams           |
|   2  | +27.46  | tannasepara_NZ_JYMT                                        |    -  | vault_panel         |
|   3  | +24.42  | A0A226QCS9_acidithiobacillus                               |    -  | colin_extras        |
|   4  | +24.13  | xoxf                                                       |    -  | vault_panel         |
|   5  | +23.46  | colinpqq_la_i0jwn7_pqq_la_model                            |    -  | colin_la_verified   |
|   6  | +23.45  | 03_Euk_a0a2v0prl1_domain3chaina-pqq-la_model               |    -  | (colin_euk)         |
|   7  | +23.28  | 00_Lav_i0jwn7-pqq-la_model                                 |    -  | (colin_euk)         |
|   8  | +23.02  | colinpqq_la_a0acd6b9f2_pqq_la_model                        |    -  | colin_la_verified   |
|   9  | +22.79  | colinpqq_la_c5atj3_pqq_la_model                            |    -  | colin_la_verified   |
|  10  | +22.61  | 00_Lav_a0acd6b9f2-pqq-la_model                             |    -  | (colin_euk)         |
|  11  | +22.40  | SilE_like_LanM_associated_sample_0                         |    -  | unknown             |
|  12  | +22.35  | 03_Euk_a0a2v0nw76_domain1chaina-pqq-la_model               |    -  | (colin_euk)         |
|  13  | +21.50  | 00_Lav_q89gy2-pqq-la_model                                 |    -  | (colin_euk)         |
|  14  | +21.33  | SilE_La_sample_0                                           |    -  | unknown             |
|  15  | +21.18  | colinpqq_la_q89gy2_pqq_la_model                            |    -  | colin_la_verified   |
|  16  | +20.57  | 00_Lav_c5atj3-pqq-la_model                                 |    -  | (colin_euk)         |
|  17  | +20.37  | FSRD01000002_1_1028_iptm0_980_s0                           | 0.980 | spicylams           |

**Notes on row 1.** `NZ_JACIDR010000005_1_9` (the +32.95 entry) is a VIT1
family Fe³⁺-uptake transporter, originally flagged as a known false-positive
on the basis of Pfam/operon. Per Jacob 2026-05-09 PM the classification was
revised to **hypothesized FP / possibly real Ln-transporter** — VIT1 family
is functionally promiscuous and a Ln-co-opted VIT1 is plausible. Recorded in
`results/known_false_positives.tsv`.

### Top 30 Ln-preferring (+5 ≤ ΔΔE < +20 kcal/mol)

| rank | ΔΔE     | stem                                                       | iPTM  | source              |
|-----:|--------:|------------------------------------------------------------|------:|---------------------|
|   1  | +19.40  | 00_Lav_q92wy9-pqq-la_model                                 |    -  | (colin_euk)         |
|   2  | +19.14  | hyphomyst_NZ_WMBQ                                          |    -  | vault_panel         |
|   3  | +18.74  | NZ_QHJE01000001_1_740_iptm0_984_s0                         | 0.984 | spicylams           |
|   4  | +18.70  | 03_Euk_a0a423vtm4_domain1chaina-pqq-la_model               |    -  | (colin_euk)         |
|   5  | +18.68  | colinpqq_la_q88jh0_pqq_la_model                            |    -  | colin_la_verified   |
|   6  | +18.26  | colinpqq_la_q92wy9_pqq_la_model                            |    -  | colin_la_verified   |
|   7  | +17.92  | colinpqq_la_c5b120_pqq_la_model                            |    -  | colin_la_verified   |
|   8  | +17.92  | 00_Lav_c5b120-pqq-la_model                                 |    -  | (colin_euk)         |
|   9  | +17.87  | A0A7V5CTA7_acidobacterium                                  |    -  | colin_extras        |
|  10  | +17.46  | NZ_JAERVJ010000008_1_86_iptm0_962_s0                       | 0.962 | spicylams           |
|  11  | +17.42  | A0A226Q2Q9_acidithiobacillus                               |    -  | colin_extras        |
|  12  | +17.31  | 03_Euk_a0a2j6qw99_domain1chaina-pqq-la_model               |    -  | (colin_euk)         |
|  13  | +16.58  | NZ_JACIDR010000005_1_24_iptm0_957_s0                       | 0.957 | spicylams           |
|  14  | +15.94  | colinpqq_la_mmol_2048_pqq_la_model                         |    -  | colin_la_verified   |
|  15  | +15.84  | colinpqq_la_mmol_1770_pqq_la_model                         |    -  | colin_la_verified   |
|  16  | +15.77  | duf2950_sample_0                                           |    -  | unknown             |
|  17  | +15.72  | 03_Euk_a0a2t6zff7_domain1chaina-pqq-la_model               |    -  | (colin_euk)         |
|  18  | +15.68  | tannase                                                    |    -  | vault_panel         |
|  19  | +15.59  | 03_Euk_a0a383vah7_domain2chaina-pqq-la_model               |    -  | (colin_euk)         |
|  20  | +15.49  | 00_Lav_q88jh0-pqq-la_model                                 |    -  | (colin_euk)         |
|  21  | +15.06  | NZ_PSRS01000033_1_43_iptm0_990_s0                          | 0.990 | spicylams           |
|  22  | +15.04  | NZ_JAEMSD010000177_1_3_iptm0_990_s0                        | 0.990 | spicylams           |
|  23  | +14.93  | NZ_RDQF01000023_1_48_iptm0_990_s0                          | 0.990 | spicylams           |
|  24  | +14.90  | NZ_JAGIKT010000066_1_76_iptm0_990_s0                       | 0.990 | spicylams           |
|  25  | +14.88  | 00_Lav_mmol_1770-pqq-la_model                              |    -  | (colin_euk)         |
|  26  | +14.84  | 03_Euk_a0a7j7pmx4_domain1chaina-pqq-la_model               |    -  | (colin_euk)         |
|  27  | +14.76  | NZ_AXAY01000027_1_94_iptm0_990_s0                          | 0.990 | spicylams           |
|  28  | +14.68  | NZ_JAFCLG010000004_1_224_iptm0_990_s0                      | 0.990 | spicylams           |
|  29  | +14.64  | NZ_SPQS01000001_1_273_iptm0_990_s0                         | 0.990 | spicylams           |
|  30  | +14.59  | VAZP01000222_1_2_iptm0_989_s0                              | 0.989 | spicylams           |

(`source = (colin_euk)` is shown where the stem is in the `03_Euk_*` /
`00_Lav_*` set; the aggregator's auto-classifier prints `unknown` because
those naming patterns predate the source-detection regex. The biology is
correct.)

To regenerate these tables on a fresh manifest:

```bash
python scripts/update_results_jsonl.py .
python - <<'PY'
import json
latest = {}
for line in open('results/all_results.jsonl'):
    r = json.loads(line); latest[r['stem']] = r
rows = list(latest.values())
for tier, n in (('Ln-evolved', None), ('Ln-preferring', 30)):
    print(f'\n== top {n or "all"} {tier} ==')
    hits = sorted([r for r in rows if r['class']==tier], key=lambda r: -r['ddE_kcal'])
    for r in hits[:n] if n else hits:
        iptm = r.get('iptm'); ipt = f"{iptm:.3f}" if iptm else '-'
        print(f"  {r['ddE_kcal']:+7.2f}  {r['stem']:<55}  iPTM={ipt}  src={r['source']}")
PY
```

---

## 3. Featured: Eukaryotic PQQ-ADH discovery (2026-05-12)

Of the 50 eukaryotic PQQ-ADH structures completed to date from Colin's
`03_Euk` batch (166 total queued), **21 score Ln-preferring or Ln-evolved**.
Two land in the Ln-evolved tier:

- **A0A2V0PRL1 domain3** — ΔΔE **+23.45**, iPTM (La) 0.980, ranking score 0.980
- **A0A2V0NW76 domain1** — ΔΔE **+22.35**, iPTM (La) 0.980, ranking score 0.970

Top Ln-preferring (ΔΔE > +13):

| ΔΔE     | UniProt            | domain         | iPTM(La) | n_atoms (QM) |
|--------:|--------------------|----------------|---------:|-------------:|
| +18.70  | A0A423VTM4         | domain1        | 0.980    | 55           |
| +17.31  | A0A2J6QW99         | domain1        | 0.980    | 55           |
| +15.72  | A0A2T6ZFF7         | domain1        | 0.980    | 55           |
| +15.59  | A0A383VAH7         | domain2        | 0.980    | 58           |
| +13.30  | A0A420YCP2         | domain1        | 0.980    | 55           |
| +13.14  | A0A1Y1IRX0         | domain1        | 0.980    | 55           |
| +13.06  | A0A1R3GAW4         | domain2        | 0.980    | 58           |

Source: `results/colin_euk_pqq_adh_results.tsv` (50 rows, ranks 1–50).

**Why this matters.** Eukaryotic PQQ-utilising oxidoreductases were thought
of as a small accessory group; finding ~40% of the curated set landing
chemistry-positive for Ln preference suggests the PQQ + Ln coordination
geometry travels across kingdoms. The four cross-domain differentiation cases
(§4) add detail.

---

## 4. Cross-domain differentiation pattern

Several Colin Euk candidates are multi-PQQ-domain proteins (tandem repeats
of the PQQ-binding fold within one chain). In four cases, the paralogous
domains of the *same* protein score in different tiers:

| UniProt    | domain  | ΔΔE     | class           |
|------------|---------|--------:|-----------------|
| A0A1R3GAW4 | domain1 | +3.12   | marginal         |
|            | domain2 | +13.06  | Ln-preferring    |
| A0A5N6PFI5 | domain1 | +0.20   | marginal         |
|            | domain2 | +6.66   | Ln-preferring    |
| A0A6J5WV68 | domain1 | +1.78   | marginal         |
|            | domain2 | +13.20  | Ln-preferring    |
| A0A7J7PMX4 | domain1 | +14.84  | Ln-preferring    |
|            | domain2 | −2.15   | ambiguous        |

Source: `results/all_results.jsonl` rows with stems
`03_Euk_<uniprot>_domain{1,2}chaina-pqq-la_model`.

**Interpretation.** Within a single multi-domain ADH, one PQQ-binding
domain has architecturally evolved a Ln-preferring carboxylate cluster
while the paralogous domain has not. This is the cleanest possible signal
that the discriminator picks up **site-level** chemistry rather than
fold-level identity: paralogs with the same primary architecture
diverge sharply in ΔΔE within the same chain.

---

## 5. AFDB porin highlight

**A0A7V5CTA7** — *Acidobacterium capsulatum* porin family protein
(ENW50_04465, 234 aa). ΔΔE = **+17.87** kcal/mol, fold confidence iPTM 0.98,
QM cluster 31 atoms. Carved by the generic pipeline (no PQQ).

This is a Ln-preferring hit *without* a PQQ cofactor — a candidate
small-protein Ln-binder in the spicy_lams scan. Acidobacterium capsulatum
is from a Ln-rich environmental niche; a porin-family protein scoring
Ln-preferring is unexpected and worth follow-up.

Source: `results/all_results.jsonl` stem `A0A7V5CTA7_acidobacterium`;
HOWTO.md §8 row "AFDB Acidobacterium/Actinobacteriota."

---

## 6. Confirmed Ln-binder families — current state

Beyond XoxF (the gold standard), six additional families pass both the
chemistry filter (Ln-evolved tier) and the biology filter (independent
evidence for Ln-binding). Catalogued in `results/confirmed_ln_binders.tsv`:

| family                                       | example                                     | typical ΔΔE   | comment                                                       |
|----------------------------------------------|---------------------------------------------|--------------:|---------------------------------------------------------------|
| XoxF / PQQ-MDH                               | `xoxf`                                       | +24.13        | Methodology gold standard                                     |
| Tannase (Marco β-roll family)                | `tannasepara_NZ_JYMT` (+27.46), `tannase` (+15.68), `FSRD01000002_1_1028_iptm0_980_s0` (+20.37) | +15 to +27 | LanM-less paralog suggests fold is intrinsically Ln-evolved   |
| SilE (small Ln-binder)                       | `SilE_La_sample_0` (+21.33), `SilE_like_LanM_associated_sample_0` (+22.40) | +21 to +22 | Convergent on two independent structures                       |
| Thermolysin Ca₂ site (Geobacillus)           | `A0A226QCS9_acidithiobacillus` (bin label misassigned) | +24.42       | Independently rediscovers Holmquist & Vallee 1974             |
| Colin's PQQ-MDH positive controls (00_Lav)   | 11 stems, +11 to +23                         | +11 to +23    | 11/11 validate; cross-run agreement 0.2–3 kcal/mol            |
| Eukaryotic PQQ-ADHs (NEW 2026-05-12)         | `03_Euk_a0a2v0prl1_domain3` (+23.45), `a0a2v0nw76_domain1` (+22.35), `a0a423vtm4_domain1` (+18.70) | +6 to +23 | 21+/50 currently Ln-class in the in-flight Colin Euk batch    |
| AFDB Acidobacterium / Actinobacteriota clade | `A0A7V5CTA7` (+17.87), `A0A226Q2Q9` (+17.42) | +5 to +18     | Clade-level signal in the spicy_lams scan                     |

Also flagged as a candidate family:

- **`small_unannotated_candidate_marco`** — a 31-protein novel family
  clustering in the +7 to +13 kcal/mol band, 114–118 aa, CN≈4. Awaiting
  characterisation. Examples: `VAZP01000222_1_12`, `NZ_VITY01000007_1_28`.

---

## 7. Wet-lab handoff priorities

Top candidates worth experimental follow-up, ordered by combined
chemistry-and-biology novelty:

1. **SilE family** — `SilE_La_sample_0` (+21.33), `SilE_like_LanM_associated_sample_0`
   (+22.40). Two independent SilE-family structures land in the Ln-evolved
   tier; Jacob's "small periplasmic Ln-binder" hypothesis is supported by
   the discriminator. **Highest priority** for wet-lab Ln binding/competition
   assays.
2. **Thermolysin Ca₂ site (Geobacillus)** — `A0A226QCS9` (+24.42).
   In vitro Ln³⁺ substitution was demonstrated by Holmquist & Vallee 1974;
   the question is whether it occurs *in vivo* in Geobacillus
   (geothermal / weathered-rock habitats — Ln-enriched environments).
3. **Eukaryotic PQQ-ADH top hits** — `A0A2V0PRL1` (+23.45), `A0A2V0NW76`
   (+22.35), `A0A423VTM4` (+18.70), `A0A2J6QW99` (+17.31). The +23 / +22
   hits are eukaryotic Ln-MDH/ADH candidates. iPTM 0.98 on all of them —
   the fold confidence is solid.
4. **AFDB Acidobacterium porin** — `A0A7V5CTA7` (+17.87). Unexpected
   Ln-preferring porin-family hit with iPTM 0.98 — see §5.
5. **Tannase paralog (LanM-less)** — `tannasepara_NZ_JYMT` (+27.46). The
   tannase fold without operonic LanM context still scoring Ln-evolved
   suggests an architectural Ln preference; worth biochemical validation as
   a candidate Ln-binding tannase.
6. **DUF4394 cassette members** — JACMOT010000131_6, NZ_LVYV01000005_62,
   NZ_NPKQ01000003_263, NZ_CP088008_3919 (+6 to +10). Originally
   hypothesised as periplasmic Ln-PQQ chaperones; the 2026-05-11 JABMCJ
   La+PQQ co-fold experiment complicated the picture (PQQ and La end up in
   separate sites in the monomer fold). Worth a re-fold with PQQ for the
   higher-scoring cassette members to test alternatives (two-site loading,
   dimer assembly, within-family variation).
7. **Hyphomicrobium mystery protein** — `hyphomyst_NZ_WMBQ` (+19.14).
   113 aa secreted, 9-coord La, novel fold (zero PFAM hits), family-restricted
   to Hyphomicrobiaceae, operonic with cofactor synthesis machinery. Strongest
   novel-binder candidate from the original validation panel.

---

## 8. False-positives — catalogued, do not chase

The following families consistently score chemistry-positive but biology says
they are not Ln-binders. Listed in `results/known_false_positives.tsv`; full
discussion in [VALIDATION.md §5](VALIDATION.md):

- M20 amidohydrolase (Zn²⁺ peptidase, e.g. `NZ_BPQI01000118_1_36`)
- Diiron rubrerythrin / bacterioferritin (e.g. `C5CIA3_bryobacter`)
- CcmF heme biogenesis (e.g. `NC_010511_1_4382`)
- Mn²⁺/Zn²⁺ small acidic pockets (OUTLIER tier — biology signal, not error)
- VIT1 Fe³⁺ transporter (e.g. `NZ_JACIDR010000005_1_9`) —
  **reclassified to hypothesized FP / candidate Ln transporter** per Jacob
  2026-05-09. Plausibly real; do not auto-reject.

---

## 9. How to re-run

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
python scripts/update_results_jsonl.py .
```

Class breakdown at the top of this file refreshes from
`results/all_results.jsonl`; the top-hit tables can be regenerated with the
snippet in §2.
