# Tiers — classification reference

> **Current agent guide — 2026-09-16:** [Alquemia operations and protocol status](docs/AGENT_PIPELINE.md).
> Only exact compatible fixed-core PQQ preparations inherit the released Ca/La bands. Generic v2, repaired v3 and the global challenger do not inherit those bands or a universal zero. Older tier names below remain historical heuristics.
> This guide supersedes conflicting operational/status prose below; dated scientific records remain historical.

Every candidate scored by `update_results_jsonl.py` gets a `class` label
based on ΔΔE(Ca − La) in kcal/mol and the QM-cluster atom count. This file
explains what each label means. The thresholds are locked by the
2026-05-07 validation panel; see [VALIDATION.md](VALIDATION.md) for the panel
and [METHODS.md](METHODS.md) for what ΔΔE actually measures.

ΔΔE is computed by [`scripts/update_results_jsonl.py`](scripts/update_results_jsonl.py):

```
ΔΔE_kcal = ((E_Ca − E_La) − DIFF_AQUO) × 627.5095
DIFF_AQUO = −646.064555551745 Ha
```

---

## Primary tiers (by ΔΔE)

| ΔΔE (kcal/mol) | class label              | what biology it represents                                  | representative stem                                          |
|---------------:|--------------------------|-------------------------------------------------------------|--------------------------------------------------------------|
| ≥ +20          | **Ln-evolved**           | architecturally amplified Ln³⁺ preference; XoxF-tier        | `xoxf` (+24.13), `tannasepara_NZ_JYMT` (+27.46)              |
| +5 to +20      | **Ln-preferring**        | intrinsic carboxylate-cluster Ln preference                 | `tannase` (+15.68), `03_Euk_a0a423vtm4_domain1...` (+18.70)  |
| 0 to +5        | **marginal**             | weakly Ln-preferring; often a paralog of a Ln-evolved hit   | `03_Euk_a0a1r3gaw4_domain1...` (+3.12)                       |
| −5 to 0        | **ambiguous**            | within the noise band; treat as neutral                     | `rifoxy` (−3.29)                                             |
| < −5           | **Ca-evolved**           | active-site geometry tuned for Ca²⁺; MxaF-tier              | `mxaf` (−6.48)                                               |
| \|ΔΔE\| > 30   | **OUTLIER** (verify SCF) | typically Mn²⁺/Zn²⁺/Fe³⁺ pocket signature — see below       | (Mn/Zn) small-acidic-pocket entries in `ln_class_hits.tsv`   |

The ±5 kcal/mol band around 0 spans the **noise floor** of the method
(PDBFixer non-determinism ≈ 0.5 kcal/mol; cross-pipeline reproducibility on
Colin's 00_Lav set was 0.2–3 kcal/mol). The thresholds at ±5 and ±20 were
chosen at panel lock so that the XoxF (+24)/MxaF (−6.5) gold-standard pair
falls cleanly on opposite sides.

The OUTLIER tier is a **useful annotation, not an error**: when ΔΔE goes
sharply negative on a small acidic pocket, the SCF is telling you the pocket
wants a smaller, harder cation (Mn²⁺, Zn²⁺, Fe³⁺) and the discriminator can't
accommodate La geometrically. Treat as biology signal — cross-reference with
Pfam/KO.

---

## Status flags (not chemistry tiers)

These are returned by the aggregator instead of a primary tier:

### `EXCLUDE (under-carved)`

Set when the La QM cluster has **fewer than 15 atoms**. Almost always
means the carver couldn't find enough donors in the first shell (cutoff
3.0 Å) and produced a near-bare-La cluster. ΔΔE is unreliable below this
atom count.

Recovery: the entry is logged to `results/recarve_queue.tsv` and reprocessed
by `scripts/process_recarve_queue.py`, which retries at 3.5 Å and routes
PQQ-bearing structures through `carve_with_pqq.py`. Still-failing entries are
re-classified as one of the two carve-status labels below.

### `SOLVENT_EXCLUDE`

Re-carve outcome. The closest non-metal O/N/S donor is farther than 4 Å from
the metal — the pocket is genuinely empty (e.g. the AF3/Protenix model
placed La in solvent, or the candidate is a true negative). Not submitted
for ORCA SPs.

### `CARVE_AMBIGUOUS`

Re-carve outcome. The closest donor is a residue the carve dict does not
include — typically CYS (SG), HIS (NE2/ND1), MET (SD), LYS (NZ), or ARG
(NH1/NH2). Per the explicit policy (Jacob, 2026-05-11; see
[CHANGELOG.md](CHANGELOG.md#2026-05-11-0110--carve-dict-tightening-tyr-and-backbone-o-carving)),
the carve dict is restricted to **ASP, GLU, ASN, GLN, SER, TYR** and is not
expanded to soft S-donors or N-aromatic donors. Such entries are flagged
for manual review rather than auto-included.

### `pending`

One or more of the four single-points (La, Ca, apo, water) has not yet
completed. Run `python scripts/update_results_jsonl.py .` after the SLURM
jobs finish to refresh.

### `SEVERE asymmetric carve` (manual flag, no auto-label)

`|ΔΔE| > 100 kcal/mol` indicates the Ca and La carves used different ligand
sets (typically a PQQ-inclusion mismatch between the two metal placements,
or one of the SPs localized SCF differently). Six entries from the
acidobacterium + 03_Euk batches have this signature. Diagnose by diffing
atom counts:

```bash
head -1 <stem>_qm/<stem>_{Ca,La}_qm.xyz
diff <(awk 'NR>2{print $1}' <stem>_qm/<stem>_Ca_qm.xyz | sort | uniq -c) \
     <(awk 'NR>2{print $1}' <stem>_qm/<stem>_La_qm.xyz | sort | uniq -c)
```

These need manual investigation; the ΔΔE value should not be trusted.

---

## How to interpret a hit

A positive ΔΔE means **the chemistry** of the carved cluster — the
combination of carboxylates, amides, hydroxyls, backbone Os, and (if
present) PQQ — favours Ln³⁺ over Ca²⁺. That is **necessary but not
sufficient** evidence for in-vivo Ln binding:

- **Chemistry favours Ln.** A positive ΔΔE says you have a hard-O cluster
  with enough donors that the larger Ln³⁺ ion is comfortable. Required for
  a real Ln-binder; the discriminator filters out Ca-evolved pockets.
- **Biology determines function.** Cross-reference with Pfam/KO, genome
  context (LanM neighbor? PQQ machinery? ECF-σ? VIT1 operon?), protein
  family, and organism. Several known false-positive families
  (M20 amidohydrolase, VIT1 Fe-transporter, diiron rubrerythrin, CcmF heme
  biogenesis, small Mn/Zn acidic pockets) score chemistry-positive but are
  not Ln-binders. See `results/known_false_positives.tsv` and
  [VALIDATION.md](VALIDATION.md).
- **Selectivity is concentration-dominated in vivo.** Most canonical
  Ca-binding proteins (calmodulin, calexcitin, calbindin, parvalbumin) also
  score Ln-preferring at the chemistry level; in vivo they bind Ca²⁺ only
  because cytoplasmic [Ca²⁺] ≈ 1 mM vs trace [Ln³⁺]. This is the same
  effect that powers Tb³⁺/Eu³⁺ luminescence-probe biochemistry. See
  [METHODS.md §7](METHODS.md) for what the discriminator does and does
  not measure.

The decision tree for a fresh Ln-evolved hit: **(1) classify**
(this file), **(2) check known FPs** (`results/known_false_positives.tsv`),
**(3) look at Pfam/KO/operon** (Sharur DB annotations are auto-joined in
`results/ln_class_hits.tsv`), then **(4) escalate** (wet-lab handoff list
in [RESULTS.md](RESULTS.md)).
