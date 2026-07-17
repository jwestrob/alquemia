# Paper 2 — Methodology Paper (In Preparation)

**Authors:** West-Roberts, Robinson, Banfield, et al. (order TBD)
**Working title:** Cluster-DFT Ca²⁺/Ln³⁺ Metal-Selectivity Scoring for AlphaFold Co-Folded Structures
**Target venue:** *J. Chem. Theory Comput.* (primary) or *J. Phys. Chem. B* (backup)
**Status:** Outline sketch, unwritten. This document holds the pitch, hooks, and open questions to package into a full manuscript.

---

## Core contribution

A physically grounded, computationally cheap Ca²⁺/Ln³⁺ selectivity scoring method for co-folded metal-binding sites that:

1. Uses single-point cluster DFT with an aquo reference to isolate site-specific selectivity
2. Achieves clean separation between biochemically-validated Ln-standards and Ca-standards (p ~ 3×10⁻⁵)
3. Scales to thousands of AF3 co-folded structures (~1 CPU-hr per protein for a Ca/La pair)
4. Provides interpretable structural covariates (CN, bidentate count) that flag near-boundary calibration cases
5. Correctly identifies why FEP is the wrong tool for this quantity, and points at DFT-quality sampling (ML potentials, QM/MM) as the honest path forward for reorganization free energy

## Story arcs the paper should address

### Arc 1: Why not FEP?

The reflex response to "compute Ca vs Ln binding" is alchemical FEP. The paper should articulate the 5-point case for why FEP is mis-specified for this quantity (see CONTEXT.md "why-not-FEP" section for the full argument). Key claims:

1. Charge change +2→+3 makes PME/finite-size corrections dominate
2. Fixed-charge FFs cannot render differential Ln-O polarization — the physics FF throws away IS the discriminating physics
3. Signal (few kcal/mol) sits on solvation magnitudes (hundreds of kcal/mol)
4. Coordination number changes (6→8-9): endpoints don't share a first shell
5. Reorganization slow relative to λ-window sampling

The relative cycle (Ca→La in site minus Ca→La in water) cancels #1 and #3 but NOT #2 or #5. Precision was never the problem; specification is.

**Actionable claim for the paper:** running FEP for this quantity produces precise-looking numbers that are systematically biased in ways that inflate confidence rather than accuracy. Cluster DFT with an aquo reference is the correct tool because Ln/Ca discrimination is fundamentally coordination chemistry, not conformational entropy.

### Arc 2: Benchmarking

- **Standards panel:** 11 Ln-standards (colinpqq_la_* set) + 14 Ca-standards (pqq_ca_* set), all biochemically validated
- **Result:** Mann-Whitney p ~ 3×10⁻⁵ on ΔΔE alone; no overlap between groups
- **Median ΔΔE:** +18.3 kcal/mol (Ln-standards), −8.2 kcal/mol (Ca-standards)
- **Cluster-size sensitivity panel:** completed on calexcitin — need to write up
- **Functional-robustness panel:** B97-3c vs r²SCAN-3c — completed, need to write up
- **Ca-class negative controls:** calmodulin, parvalbumin, calbindin — completed, need to write up (parvalbumin 1B9A La SP needed TightSCF + level shift; document this as a known convergence gotcha)

### Arc 3: CN and bidentate structural covariates

Cross-table of CN @3Å × n_bidentate_asp_glu reveals structural patterns that inform calibration:

- **CN=5 with 1 bidentate** (77% of CN=5): canonical Ln-preferring non-cofactor geometry
- **CN=8 with 1 bidentate** (84% of CN=8): classic XoxF geometry (2D-1E-1N + 4 PQQ-O with one Asp bidentate)
- **CN=9 with 2 bidentate** (98% of CN=9): PQQ-bound signature; PQQ contributes 2 bidentate carboxylates
- **CN=4 with 1 bidentate** (61% of CN=4): "effectively 3-residue" sites — sub-saturated

Paper should propose using CN + bidentate as pre-filter and calibration diagnostic. Add cross-table figure.

**Filter policy:** `carve_generic.py` skips CN ≤ 3 by default (18% of raw candidates). Motivated by cross-table.

### Arc 4: The ExaF calibration-boundary case

ExaF (C5AXV8), a validated strict-Ln-dependent enzyme (Good et al. 2016), scores ΔΔE = +10 kcal/mol — well below expected. Two independent AF3 co-folds converged on the same CN=7 (plus one donor at 3.10 Å). Diagnosis: **AF3 systematically displaces one donor ~0.3–0.5 Å from crystallographic Ln-O distance.**

**Verification (needed for paper):** enumerate 8 first-shell donor distances in PDB 6OC6 (XoxF crystal, La-bound), compare to AF3 XoxF co-fold. If crystal is at ≤2.9 Å and AF3 at ~3.1 Å for the same donor position, that's a methodological finding: **AF3's PQQ-Ln geometry has sub-Angstrom systematic bias** that our carve+DFT protocol propagates into ΔΔE.

Options to mitigate:
- Multi-seed averaging (cheap)
- Brief QM geometry optimization of first shell (medium cost, most principled)
- Larger cluster + geometry constraint from crystal template (expensive)

### Arc 5: PQQ vs alternative quinone cofactors

Structural + coordination + residue-signature arguments for why our eukaryotic Ln-PQQ hits are PQQ-binders and NOT TPQ/TTQ/LTQ/CTQ enzymes:

- 8-bladed β-propeller fold is unique to PQQ enzymes
- 8-oxygen coordination (2D + 1E + 1N + 4 PQQ-O) is diagnostic of PQQ
- Absence of TPQ/TTQ/LTQ/CTQ-forming residue signatures in active site (Tyr for TPQ, 2×Trp for TTQ, Lys+Tyr for LTQ, Cys+Trp for CTQ)
- iPTM ≥ 0.8 with exogenous PQQ specifically

This is more of a discussion-section item than a core methods contribution, but worth including as a sanity check.

### Arc 6: Forward pointer to sampled free energy

The static cluster DFT ΔΔE does not include reorganization free energy. The honest way to recover it:

- **ML potential trained on cluster DFT reference data** — gives sampling at DFT-quality electronics; can then run MBAR-like free-energy estimators
- **QM/MM with the QM region matching the cluster carve** — pay full DFT cost but with explicit sampling
- **NOT FEP** — see Arc 1

This is a "future work" section, not something to complete in Paper 2. But naming it correctly matters — future readers should know the discovery paper used static ΔΔE for coordination chemistry (the dominant term) and reorganization sampling is a separate methodological effort.

## Data and code to include

- Full cluster-DFT recipe (carve_generic.py, carve_with_pqq.py, ORCA input template, submit script template)
- Aquo reference calculation protocol
- Standards panel (25 stems: 11 Ln + 14 Ca) as supplementary data
- Cluster XYZ files for the standards panel
- CN/bidentate cross-table (all 3985 entries or a representative subset)
- Sensitivity panels (cluster size, functional)
- Alquemia repo tag or DOI at publication time (Zenodo)

## Open methodological questions to close before writing

1. AF3 geometry-bias magnitude (crystal vs AF3 donor distance comparison for XoxF, MxaF, and 2-3 other structurally-characterized PQQ enzymes)
2. Does the bias persist at multi-seed AF3 averaging, or does it converge to crystal?
3. Cluster-size sensitivity on the standards panel — do medians shift when we go from 3.0 Å to 3.5 Å to 4.0 Å carve?
4. Do outliers (|ΔΔE| > 30) share structural signatures that our CN/bidentate covariates don't capture?
5. What's the equivalent number for the Ca side reference? DIFF_AQUO is the Ca-La aquo swap; what would the corresponding number be for Ca-Nd, Ca-Sm, Ca-Eu, etc.? (Enables extending the discriminator beyond La to the full Ln series.)

## Anti-goals (things NOT to attempt in Paper 2)

- Wet-lab validation of any specific candidate (that's Paper 1 territory and future biochemistry work)
- Comparison to other software (Gaussian, Q-Chem, Psi4) unless a reviewer specifically demands
- Detailed comparison to bond-valence sum (BVS) methods — the DFT results have superseded BVS in this project
- Extension to actinides (project scope has been strictly Ln)
