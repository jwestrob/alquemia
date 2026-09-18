# Khoury 2025 author-model La/Ca challenge — declared before scores

2026-09-18. Jacob requested this paper benchmark in parallel with discriminator
improvements. Contained pilot authorized by the active goal and explicit pilot
permission. Baseline/default and historical references remain unchanged.

## Evidence and fixed selection

Use all modeled sites in author-supplied AF3 Ca-conditioned domain structures:
A0A7 (s001, 87 residues, six Ca), HEW5 (s003, 117 residues, eight Ca), and
RTX (s008, 152 residues, eight Ca). Original Ca-chain order defines site order.
All seven supplied PDBs remain in the inventory; the other four are not scored.
These are predictions, not crystal structures. Table S3 confirms their domain
sequences. Experimental ITC constructs retain an N-terminal MPVP cleavage scar
absent from these models. This pilot tests author domain models, not exact
experimental constructs; no scar modeling or new fold is permitted here.

Table S6 (SI p20) reports La ITC Kd A0A7 17±2, HEW5 5.2±1.7, RTX40±4 µM
(95% confidence intervals). Calcium values ≈2000, ≈750, ≈250 µM are CD folding
concentrations, not fitted Kd. ITC: 50mM MES,50mM NaCl,pH6,25°C,30µM protein.
CD: 1mM MES,50mM KCl,pH6,25°C,60µM protein. This is qualified protein-level
multisite/cross-readout La-favoring evidence; no site-level affinity labels.
No weak-XO negatives, FRET La labels, or HEW5>A0A7 ordering are inferred.

## Preparation and method

Preserve every source heavy coordinate and all Ca sites. Normalize metal-chain
identifiers only through an exact atom map, without changing assembly. Native
OpenMM hydrogen preparation uses ff19SB, pH6, Reference platform, fixed seed;
retain the generated output and exact source mapping. Existing terminal-OXT
completion is allowed if needed. Missing internal atoms, broken peptides or
unsupported chemistry fail explicitly. No water is present or added. Paired
Ca/La endpoints share one protonation/H geometry; all other ions remain Ca2+.
The old multisite policy/default is unchanged. New policy ID:
`omol_author_AF3_domain_multisite_pH6_fixed_Ca_ff19sb_v1`.

Use the existing charge-feature-masked MACE OMOL checkpoint/software and
qualified exact two-call disconnected-atom factorization. For each of22 sites,
calculate bound Ca and La:44 endpoint forwards total. No new DFT, folding,
AMOEBA, geometry relaxation or learned-model fitting. Do not reuse all-Ca
energies unless exact prepared/cache equivalence is demonstrated; initial
manifest retains all44 tasks. Existing one-A5000/16CPU/64474MiB host allocation.
Record actual wall/CPU/GPU allocation cost and failures; no project time cap.

## Declared interpretation

Report all ordered site vectors, including null/invalid entries. Primary
protein summary is the unweighted arithmetic mean over ALL modeled sites;
any unavailable site makes that domain mean unavailable. No selecting best
sites or post-hoc aggregation. The primary comparison is each of three domain
means minus each consumed GGR descriptor (1GLG,2FW0,2FVY):nine comparisons,
grouped by three domain families versus one GGR family. Direction passes only
if the difference exceeds the previously used0.02 model-kcal screen. Larger
descriptor is La-like. These are exploratory cross-family comparisons, not
nine independent observations or direct predictions of measured Kd. No PQQ
decision bands or universal zero. Exact La Kd rank and within-Ln preferences
are separate questions and not acceptance criteria here.

The models/labels were inspected during selection, while their scores were
unseen before this plan. This is a newly scored exploratory challenge to a
frozen existing descriptor, not a blind clinical-style validation. Search
existing project records for prior accessions/scores before reporting exposure.

Sources: DOI10.1039/D5SC02315G; main XML, SI PDF and original author PDBs pinned
under workspaces/mace_omol_20260917/d5sc02315g_reading_v1. Raw products and
finite manifests belong under workspaces/mace_omol_20260917/khoury_*.
