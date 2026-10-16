# A small consumed monomer benchmark for La/Dy

16 October 2026. Evidence and source inventory only; no scientific calculations,
molecular preparation, or jobs. This note follows Jacob's request to start with
monomers and retains protein-level LanM discrimination as the larger goal.

## Recommendation

**Start with Mex-LanM and the already-consumed LBT3 peptide.** They provide
opposing within-series directions without requiring a dimer model: Mex's modest
light-lanthanide preference versus LBT3's directly measured Dy-over-La binding.
This is one natural protein plus one engineered peptide, not two natural LanMs.
It is a useful small mechanism/control panel, not sufficient validation of
protein-level LanM prediction. No second natural monomer with both a verified
La/Dy label and an actual deposited structure was established in this bounded
consumed inventory.

| Candidate | Defensible measured target | Actual structural starting points | Main limitation |
|---|---|---|---|
| Mex-LanM WT | Modest light-over-heavy apparent affinity; qualitative La/Dy direction from the consumed pH-5 comparison | Nd-bound monomer 8FNS; Y-bound NMR 6MI5 | Cooperative multisite folding/occupancy; no La- or Dy-conditioned deposit in this set |
| LBT3 | Dy binds more favorably than La in the same ITC experiment | La-bound 7CCO and Lu-bound 7CCN, both monomeric NMR ensembles | A 17-residue peptide; no Dy deposit; different assay from Mex |

Mex's evidence is the primary [2023 Hans/Mex comparison, Fig. 1e](https://www.nature.com/articles/s41586-023-05945-5),
using earlier [Mex measurements](https://doi.org/10.1021/acs.inorgchem.0c01303).
Exact matched Mex La/Dy values and uncertainties were not recovered here; retain
the qualitative target without digitizing a quantitative benchmark. The article
describes Mex as monomeric. Its EF1–3 high-affinity response and weak EF4 are not
four identical sites. See `EXPERIMENTAL_TARGET.md` for the fuller observable and
construct distinctions.

## Direct LBT3 evidence, already consumed

LBT3 sequence is `FIDTNNDGWIEGDELLA`; LBTv is a different peptide.
The original [Hatanaka et al. paper](https://doi.org/10.1038/s41598-020-76527-y)
and [Supplementary Table S1 / Fig. S5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-020-76527-y/MediaObjects/41598_2020_76527_MOESM1_ESM.pdf)
report:

| ITC quantity | La | Dy |
|---|---:|---:|
| Binding free energy, kcal/mol | −6.82 ± 0.03 | −9.18 ± 0.11 |
| Fitted stoichiometry N | 1.01 ± 0.05 | 0.97 ± 0.08 |

Conditions: 10 °C, pH 6.0, 50 mM MES, 100 mM NaCl; 15 µM peptide,
300 µM lanthanide nitrate injectant. Means use at least three experiments;
the corresponding main figure specifies standard-error bars, while Table S1
does not separately restate the uncertainty convention. The methods allow
stoichiometry to be fixed at one for weak nonsigmoidal titrations.
This is a binding-thermodynamic measurement, not fluorescence brightness,
capture loading, or a separation factor. No new ratios or free energies were
calculated for this note.

Gel filtration supports monomers for apo, La-bound and Lu-bound LBT3; no separate
Dy gel-filtration trace was verified. [7CCO](https://www.rcsb.org/structure/7CCO)
and [7CCN](https://www.rcsb.org/structure/7CCN) are author-assigned monomers,
each with 20 submitted solution-NMR conformers and model 1 designated the best
representative. Their coordination restraints limit independent interpretation
of the metal geometry. Retain both endpoints; a Lu-shaped source is not an
observed Dy-bound structure.

## Exact local sources and preparation state

All paths below are repository-relative unless absolute.

**Mex:**

- Actual Nd crystal: `workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif`.
- Actual Y NMR: `/groups/banfield/projects/environmental/sr/srvp2020/protenix/mmcif/6mi5.cif`.
- Existing accepted source contract:
  `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`.
- Accepted La/Dy collection:
  `workspaces/metal_environment_response_20260926/lady_repaired_origins_Mex8FNS_v1/FINAL_COLLECTION.json`.
- Matched classical terms:
  `workspaces/metal_environment_response_20260926/scaffold_origin_transfer_v1/Mex8FNS/{La,Dy}_A.json`.
- Current prepared state: isolated chain A, residues 29–133; protein charge −10,
  EF3 QM charge 0, seven QM waters / 167 total waters, three Nd spectators.
  Electronic spectators are +3 charges; the classical ledger also has
  ion-specific pure 12–6 interactions. This is one conditional occupied state,
  not the measured full-protein binding equilibrium.

**LBT3:**

- Actual deposits: `lanthanide_series_capability/lbt3/raw/{7CCO,7CCN}.pdb`.
- Preserved model-1 sources:
  `lanthanide_series_capability/lbt3/first_models/{7CCO,7CCN}_model1.pdb`.
- Exact experimental transcription:
  `lanthanide_series_capability/lbt3/experimental/itc_table_s1.tsv`;
  prior consumed declaration: `lanthanide_series_capability/LBT3_CAPABILITY_PREREGISTRATION.md`.
- Historical preparation:
  `lanthanide_series_capability/lbt3/benchmark_manifest.json`,
  `build_benchmark.py`, and `clusters/{7CCO_model1_La,7CCN_model1_Lu}_cluster.xyz`.
- Those legacy clusters contain 47 ligand atoms plus one trivalent metal,
  ligand charge −4 / total charge −1. They cap D3/N5/D7/E11/E14 side chains
  and W9 backbone carbonyl; they are not intact-peptide preparations.
  The model-1 cluster hashes are respectively
  `66407a94d81d4d0da873e4fab1d8bd9a1148fe5540ad34f343d4befd8647112e` and
  `48242c2411761651484d37137e02844ea34d0b6db812973c8d88ad7a5c606695`.

The old LBT3 series calculation is not a qualified reusable affinity method:
`diagnostics/lanm_series_followup_20260923/READINESS.md` records geometry/aquo
dependence and open-shell failures. Do not rerun that campaign. A profitable
new pilot first needs a declared intact-peptide charge/terminus policy,
backbone treatment, water inventory and metal-state contract. Absence of
deposited NMR waters does not establish a dry coordination site. La/Dy exchange
must have a consistent reference or balanced relative cycle; raw cross-element
total-energy subtraction is not affinity. Use the existing Mex collection for
compatible diagnostics rather than restarting its solved SCF calculation.

**Smallest useful question:** does one declared response/hydration treatment
retain LBT3's Dy-over-La direction from both consumed endpoint geometries while
allowing a different response in Mex? Preserve assay-specific signs; do not fit
one numerical threshold to ITC at pH 6 and cooperative CD at pH 5, and do not
assign the Mex protein label to isolated EF3. Full-protein site/occupancy and
folding treatment remains required for the larger LanM goal.

## Bounded alternatives

- **Mex 4P2A:** already consumed; four substitutions P36A/P60A/P85A/P109A.
  The exact [2018 JACS supplement](https://ndownloader.figshare.com/files/13330367)
  supports altered response and a monomer-like La-loaded SEC size, but no
  direct mutant La/Dy pair was verified. No deposited mutant structure was
  identified. Do not assign a Dy label from generic lanthanide language.
  No direct La/Dy pair for a Mex D9 mutant was confirmed in the bounded check.
- **Hans R100K:** strongest existing protein-level extension: measured La/Dy
  apparent affinities and increased monomerization, but no deposited mutant
  structure in the inspected set. Its state is condition-dependent, so a
  mutated isolated-chain model would remain a hypothesis. Exact measurements
  and assembly limits are in `EXPERIMENTAL_TARGET.md`.
- **Lanpepsy / LanP, 9VY7:** biological C2 dimer with shared coordination;
  excluded from this monomer panel. Local source:
  `diagnostics/laca_benchmark_expansion_20260915/sources/pdb/9VY7.cif`.
- **Landiscernin / LanD, 9C8X:** measured dilute monomer binding exists, but
  deposited holo structures are metal-centered dimers; no direct quantitative
  La/Dy pair was verified here. [Primary paper](https://doi.org/10.1073/pnas.2410926121).
- **Xan-LanM / Al-LanM:** the new [primary article](https://doi.org/10.1021/jacs.6c08525)
  explicitly studies dimerizing LanMs; not simpler monomer alternatives.

## Reserved structures: metadata only

`workspaces/spicy_lams_inventory_20260923/inventory_v2/STRUCTURES.json` indexes
32 one-chain, protein-only AF2 predictions for these 16 reserved sentinels:

| Metadata pair | Accessions |
|---|---|
| o382 / o383 | WP_198355877.1 / WP_216634451.1 |
| o398 / o399 | WP_109971873.1 / WP_236990760.1 |
| o199 / o208 | WP_050405229.1 / WP_076859038.1 |
| o126 / o128 | WP_196407746.1 / WP_036293838.1 |
| o372 / o373 | WP_092869117.1 / MBY0561954.1 |
| o171 / o416 | WP_012589863.1 / WP_102844174.1 |
| o516 / o520 | MBE7203409.1 / MBE7198988.1 |
| o15 / o27 | MBE7249430.1 / MBE7247691.1 |

A one-chain prediction does not establish a biological monomer, occupied sites,
or a deposited holo structure. These are not selected by preference or opened
for outcomes. Nonreserved FSRD-LanM/o526/SIO63612.1 has existing four-La predicted
models, including
`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/candidate_bundle/lanm_calibration/FSRD_LanM_gene1031_4La_iPTM0p951.cif`;
no verified deposited structure or matched monomer-state evidence was found.

The reserved numerical files and PLM holdouts were not opened. Broad web searches
unavoidably returned snippets from the reserved-family paper despite exclusions;
those results were not opened or used, and broad searching stopped. Earlier
documented exposure also remains recorded. This panel is consumed developmental
evidence and makes no blanket-blindness claim.
