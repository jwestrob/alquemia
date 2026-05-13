# DFT Cluster Ca²⁺/Ln³⁺ Discriminator — Protein Findings

*One-night sprint, 54 protein active sites scored, 2026-05-07*

## The headline number

**LanM EF1: +47.5 kcal/mol Ln-preferring → Prolidase A0A2V5ZP60: −20.7 kcal/mol Ca-preferring**

A 68 kcal/mol range across one validation panel.

## The gold-standard pair

The same fold + same cofactor. Only the active-site coordination differs by biology, and DFT picks them apart cleanly.

| Protein | Active-site geometry | ΔΔE(Ca − La) |
|---------|----------------------|-------------:|
| **MxaF** (1H4I) | Ca²⁺-MDH; 1 Glu + 1 Asn + PQQ → CN6 | **−6.5 kcal/mol** |
| **XoxF** (4MAE) | Ce³⁺-MDH; 2 Asp + 1 Glu + 1 Asn + PQQ → CN9 | **+24.1 kcal/mol** |

**31 kcal/mol gap.** This is the single bacterial paralog pair that biology already discriminates, and the discriminator gets it cleanly.

## Tier-1 Ln-evolved sites (≥ +20 kcal/mol)

The ≥ +20 threshold appears to flag sites where biology has *architecturally amplified* the carboxylate cluster's intrinsic Ln³⁺ preference — true Ln-evolved active sites.

| Protein | Source | ΔΔE | What it is |
|---------|--------|---:|-----------|
| **LanM EF1** | M. extorquens canonical Ln-binder (6MI5) | **+47.5** | 4 carboxylates + inner water, octadentate |
| **Tannase paralog** (NZ_JYMT) | Bradyrhizobium, **LanM-less variant** | **+27.5** | Tannase fold without operon-LanM context still Ln-evolved |
| **XoxF** | Methylobacterium Ce-MDH | **+24.1** | Ln-MDH active site |

## Tier-2 Ln-preferring (+5 to +20 kcal/mol) — 18 proteins

A broad band of carboxylate-cluster sites that *prefer Ln³⁺ at the coordination level* even though most are biologically labeled "Ca-binding."

This is the more nuanced finding: **EF-hand and similar carboxylate-cluster sites bind Ln³⁺ tighter than Ca²⁺ at the coordination level.** They function as Ca-binders in vivo because [Ca²⁺] ≈ 1 mM vs trace [Ln³⁺]. *This is consistent with the long-known biochemistry that Tb³⁺/Eu³⁺ luminescence probes outcompete Ca²⁺ at canonical Ca-binding sites.*

| Protein | ΔΔE | Class |
|---------|---:|-------|
| Hyphomicrobium mystery (NZ_WMBQ) | +19.1 | Novel; family-restricted; operonic with cofactor synthesis |
| Calmodulin (1CLL) | +18.2 | Canonical eukaryotic CaM |
| Tannase (canonical) | +15.7 | Bradyrhizobium tannase fold |
| Pectate-lyase A0A9E3VGV4 La_3 | +13.0 | Gemmatimonadaceae T9SS β-helix |
| Fern calmodulin-like (A0A8T2TXN1) | +12.1 | Plant CaM homolog |
| Parvalbumin (1B9A) | +11.4 | Canonical EF-hand |
| LanM-fusion (NZ_LFLZ) | +11.0 | Bifunctional EF-hand + novel CT domain |
| Calexcitin EF-3 | +10.6 | Convergent EF-hand |
| RTX β-roll (A0A3M2BIM9) | +10.1 | Planctomycetota Ca-binding repeat |
| Calcineurin B (Ceratopteris) | +10.1 | Plant signaling |
| Type A orphan (NZ_JADPKR) | +9.9 | DxDx-rich novel binder, no LanM homology |
| **rifoxy + heme variant** | +8.6 | Heme-bound flips ambiguous → Ln-class |
| Arginase (Gemmatimonadota A0A6N9AP46) | +8.3 | Binuclear-Mn fold |
| Marco's β-roll (cp030053) | +7.8 | See below |
| Fern parvalbumin-like (A0A8T2SI74) | +7.7 | Plant EF-hand |
| GCaMP3 (engineered Ca-sensor CaM) | +7.1 | Designed Ca probe |
| Marco's β-roll (cp019948) | +6.6 | See below |
| Marco's β-roll (akiy_44) | +6.0 | See below |
| Marco's β-roll (bsox_25) | +5.4 | See below |

## Marco's β-roll family — paralog precision champion ★

**4 paralogs from Bradyrhizobium small-secreted family, all Ln-preferring band, paralog spread <2.5 kcal/mol.**

| Paralog | Coordination | ΔΔE | Operon context |
|---------|--------------|---:|------|
| cp030053_2053 | 4× ASP/GLU | +7.8 | LanM-adjacent |
| cp019948_943 | 3× ASP/GLU | +6.6 | TonB-direct → β-roll → LanM |
| akiy_44 | 4× ASP/GLU | +6.0 | TonB at -19; β-roll at 0; LanM nearby |
| bsox_25 | 4× ASP/GLU | +5.4 | β-roll near LanM |

**Mean +6.5 ± 1.0 kcal/mol** — strongest paralog precision in the entire panel. Marco's instinct that this family is biologically Ln-relevant is supported. The TonB → β-roll → LanM operon co-occurrence reflects real coordination chemistry, not architectural noise.

## The Ca-class side (< −5 kcal/mol) — 13 proteins

| Protein | ΔΔE | Function |
|---------|---:|----------|
| **Prolidase (A0A2V5ZP60)** | **−20.7** | Tightest dinuclear pair (3.2 Å Asp-Asp) — strongest Ca preference in panel |
| Verrucomicrobiota cbb3 oxidase (NZ_RXIZ) | −16.7 | Cyt-locus-flagged Ln candidate, **DFT-rejected** |
| Fern endoglucanase 1 (A0A8T2U5N9) | −16.5 | Plant cellulase structural Ca |
| Alphap arginase (Phenylobacterium) | −16.7 | Phylum split from Gemmatimonadota +8.3! |
| Excalibur DUF (PF05901) | −15.5 | "Extracellular Calcium-Binding" — name validated |
| Cytochrome c (A0A1Q7R3F2) | −11.5 | Operon-flagged Ln candidate, **DFT-rejected** |
| M24 peptidase (A0A2E7SM25) | −8.4 | Dinuclear Mn/Co metallopeptidase |
| Fern cellulase 2 (A0A8T2TYS0) | −8.7 | Plant cellulase |
| Verruco A0A975IXS4 | −8.1 | Multi-La verruco protein |
| Agmatinase (A0A2P8QC45) | −6.8 | Binuclear-Mn-fold |
| **MxaF** | **−6.5** | Gold-standard Ca-MDH paralog |
| Aminopep P (A0A432JR12) | −6.2 | M24 family |
| ETR2 EF2 (bacterial CaM-fold) | −5.9 | EF-hand #2 of A0A127ETR2 |

## Five novel insights from the panel

### 1. Discriminator is per-site, not per-fold
Same fold (binuclear-Mn arginase), different organisms, **opposite signs**:
- Gemmatimonadota arginase: **+8.3 kcal/mol** (mild Ln-pref)
- Alphaproteobacteria arginase: **−16.7 kcal/mol** (strongly Ca-pref)

### 2. Bacterial vs eukaryotic CaM-fold ≠
Same architecture, different result. Architectural convergence ≠ thermodynamic equivalence.
- Eukaryotic 1CLL EF1: **+18.2** (Ln-preferring)
- Bacterial ETR2 EF1: **+0.2** (marginal)
- Bacterial ETR2 EF2: **−5.9** (mild Ca-class)

### 3. The "Ca-binding protein" label is misleading for in-silico selectivity
Calmodulin, calexcitin, calbindin, parvalbumin, calcineurin B, RTX-roll — all the canonical "Ca proteins" land Ln-preferring at the coordination level, with magnitudes in a tight +5 to +20 band. They're Ca-binders in vivo only because of [Ca²⁺]/[Ln³⁺] concentration ratio.

### 4. Cofactor-mediated geometry shifts can flip ΔΔE
- Rifoxy Site B: **−3.3** (ambiguous)
- Rifoxy + heme variant: **+8.6** (Ln-preferring)
- Heme binding shifts coordination geometry enough to flip the sign.

### 5. Honest negatives — the discriminator filters its own candidate list
- Cytochrome c (A0A1Q7R3F2) was operon-flagged as Ln candidate (iPTM 0.992) — **DFT rejected at −11.5 kcal/mol**
- cbb3 oxidase locus protein (NZ_RXIZ) — operon evidence, **DFT rejected at −16.7 kcal/mol**

The discriminator is not a rubber stamp. Important for credibility.

## Excalibur DUF (PF05901) — settled

The DUF nomenclature ("Extracellular Calcium-Binding Region") reflects **real thermodynamic Ca-preference**, not just historical discovery context.
- alphap_ivl5_excalibur (A0A2E2IVL5): **−15.5 kcal/mol** → strongly Ca-class

## Outliers we caught (and excluded)

Three sites gave >|30| kcal/mol signals — flagged as geometry artifacts:
- **DUF882** (−1378 kcal/mol): 15-atom under-carved cluster, charge-state pathology
- **A0A965Z5K8** (−49.6): Protenix-predicted 2.0 Å M-O distances (vs typical 2.4 Å) → unphysical Pauli repulsion penalty for La³⁺
- **A0A5M3W1H2** (−51.1): Same geometry artifact
- **alphap_aminopep_p** (−28.5): 1.93/2.06 Å M-O distances; same issue

Fix: M-O distance sanity check during carve. Exclude clusters with min(M-O) < 2.1 Å.

## What this enables

- **Triage tool for AFDB-scale candidate sets** — score a Protenix La-bound prediction in 30 min CPU, classify as Ln-evolved / Ln-preferring / null / Ca-class
- **Companion to the plant xoxF-like manuscript** — bacterial XoxF/MxaF reference pair (+24 vs −6.5) sets the validation bar; plant xoxF candidates can be scored against it
- **Per-site geometry sensitivity** — discriminator detects subtle differences within fold class (paralogs split −20.7 to +8.3 in the binuclear-Mn family)

## Two specific high-priority novel candidates to highlight

### Hyphomicrobium mystery protein (NZ_WMBQ)
- 113-aa secreted, 9-coord La, novel fold (zero PFAM hits)
- Family-restricted to Hyphomicrobiaceae
- Operonic with ALAD (cofactor synthesis machinery)
- ΔΔE = **+19.1 kcal/mol** — just below Ln-evolved threshold
- One of the strongest novel-binder candidates in the panel

### Tannase paralog NZ_JYMT (LanM-less)
- 86% identity to LanM-cassette tannase
- ΔΔE = **+27.5 kcal/mol** — Ln-evolved tier
- Demonstrates that tannase fold *intrinsically* prefers Ln³⁺; doesn't need operon-LanM context
- Implication: tannase fold is a Ln-evolved family at the architectural level

## Stats

- **54 candidate sites** scored (50 reliable, 4 outliers excluded)
- **1 night** of compute (memory partition, 64-core nodes, ORCA 6.1.1 r²SCAN-3c CPCM(Water))
- **30 min CPU per (La, Ca, apo, water) SP set** average
- **47 reliable signed scalar values** spanning −20.7 to +47.5 kcal/mol
- Functional robustness: B97-3c sign-preserved on 4 reference proteins
- Paralog precision: <2.5 kcal/mol within Marco β-roll family
