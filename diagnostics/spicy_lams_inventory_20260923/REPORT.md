# Spicy-Lams input inventory and reserved-outcome boundaries

The local workbook contains **616 unique accessions and616 unique mature
sequences** (81–137aa). All are exported with stable joins and source paths.
**18 proteins have exact mature-chain structures in the inspected locations**:
16 reserved sentinel proteins, Hans-LanM and FSRD-LanM. Mex-LanM has related
experimental constructs. No structure was located for the remaining597 in these
bounded locations; this is not a claim that none exists elsewhere on biotite.

This inventory performed no molecular calculation, new prediction, job, label
assignment, source move or reserved-outcome lookup. Existing selection,
preregistration and structure checksums all pass. Only workbook `Data!A:S` and
`Data!BM:BS` were accessed; ion, normalized-logD, cluster and UMAP columns remain
unread. **Whole-project blindness is incomplete:** prior parent-agent article
browsing exposed published cluster/qualitative information for reserved member
o-27 (`MBE7247691.1`, pair p08): C7, poor expression and weak binding. Numeric
reserved workbook profiles remain unopened. This limited exposure is recorded
without changing panel membership or reading further outcomes; see [evidence](EVIDENCE.md).

## Usable files and joins

Repository-relative workspace root:
`workspaces/spicy_lams_inventory_20260923/`.

| Product | Exact relative path | Purpose |
|---|---|---|
| All616 metadata rows | `inventory_v2/ORTHOLOGS.tsv` | workbook row,ID4,ID6/o-ID,accession,sequence hash,taxonomy,EF annotation,reserved membership,source assembly/files |
| Exact mature sequences | `inventory_v2/MATURE_SEQUENCES.faa` | all616, no tag/proton/sequence normalization |
| Reserved16 folding inputs/results | `inventory_v2/RESERVED_FOLDS.tsv` | exact FASTA/MSA pins, both existing ranks, confidence-file paths |
| Structure composition/sequence pins | `inventory_v2/STRUCTURES.json` | full-chain versus construct identity, metal/water inventory, exact file hashes |
| Per-protein structure coverage | `delivery_v1/STRUCTURE_COVERAGE.tsv` | exact sequence joins; missing/related-construct statuses explicit |
| Global/local resource index | `delivery_v1/GLOBAL_AND_LOCAL_RESOURCES.tsv` | existing models/MD/local cores, status and evidence paths |
| Existing files | `inventory_v2/{PROJECT_FILES,REPRESENTATION_FILES}.tsv`, `delivery_v1/ADDITIONAL_REPRESENTATION_FILES.tsv` | paths/types/sizes; trajectories and energies not reanalysed |
| Delivery pins | `delivery_v1/DELIVERY.json` | joinable artifact hashes/counts;[compact pointer](ARTIFACTS.json) |

The dataset has92three-EF and524four-EF annotations; motifs are not occupied-site
measurements. The paper describes a621-ortholog library, while this exact local
sheet has616rows acrossID6range0–621, with371,425,432,544,595,620 absent. Retain
these actual denominators; do not invent missing records or infer exclusion causes.

## Structures and constructs

- **Reserved16:** all16fold directories completed,32unrelaxed AF2 models total,
  two ranks per protein. Each model exactly preserves its full mature sequence.
  All32 are protein-only: no metals or waters are present. They are conditioned
  by an8FNS template; they do not establish occupancy, apo/holo equilibrium or
  a four-ion folded state. Rank order and EF1–3 primary policy remain frozen.
  No post-fold site transfer/core manifest was located in the sentinel directory;
  do not describe the models as ready15-ion energetic inputs.
- **Hans:**8DQ2's four110-aa chains exactly match workbooko-180,
  `WP_131004249.1`. It contains12La and4Na in total; native crystal waters remain.
- **FSRD:** workbooko-526,`SIO63612.1`, matches the110-aa five-sample Protenix
  four-La trial and the separate calibration-bundle model. These are the same
  biological sequence, not six independent proteins. Four requested/modelled
  ions do not validate four-site occupancy.
- **Mex:** workbooko-621,`WP_003601797.1`, is113aa. Local8FNS is105aa and exactly
  omits its leading`MAPTTTTK`;6MI5 equals the workbook sequence minus initial`MA`
  plus a C-terminal`HHHHHH`. Equivalently6MI5=`PTTTTK`+8FNS+`HHHHHH` (117aa).
  Preserve these source constructs; no tags were stripped or added.8FNS contains
  fourNd;6MI5 has12NMRmodels, each with threeY. The old statement about native
  sequence minus one Lys concerns a different normalization than raw deposited
  or workbook constructs and must not silently govern a new preparation.
- The `lanm_benchmark/structures` folder also contains PQQ enzymes6DAM/6OC5/6OC6;
  directory membership does not make them LanM structures.

No two inspected structure files have identical whole-file hashes. Multiple
ranks, seeds, chains and repeated sequence copies remain linked, not deduplicated
into fictitious biological replicates.

## Sharur discovery data stays separate

`Jacob/spicy_lams/metadata/lanm_download_summary.tsv` links595assay accessions to
375unique assemblies. All595 associated genome/proteome paths were located.
Their stored paths omit`/Jacob` and are stale; ORTHOLOGS retains both stored and
verified located paths without editing old metadata.21rows have no mapped local
assembly/proteome; their mature sequences still exist in the workbook/FASTA.

The GTDB discovery extension currently has2349genome entries but1856proteome and
1856annotation entries, including177working symlinks in each category. These
actual counts differ from the older briefing's complete-download claim. Whole
proteomes are discovery backgrounds, not experimentally characterized LanMs.
Their precursor proteins have not been newly extracted or aligned here.

`candidate_bundle/spici_lams_other/` holds four separate predicted proteins
(DUF882/TonB,NlpC,TerB,cbb3-locus;150–671aa), each with one modelledLa and no exact
mature-sequence match to the616. Preserve these as Sharur discovery candidates,
not SpyCI-LAMBS ortholog labels. The FSRD exact join is separately explicit.

## Existing global representations and their limits

- Mex full117-aa FEP preparations exist at `lanthanide_binding/fep/prep/`:
  apo protein and **three separate one-La EF1/EF2/EF3 models**, plus aqueousLa.
  They do not represent simultaneous three- or four-metal loading. Explicit-
  solvent trajectories and topology files remain indexed. Historical notes report
  poor series magnitudes/direction; no new analysis or reuse qualification here.
- Hans Amber12-6-4 has20completed trajectories and actual parameterized solvated
  models. One inspected source retains the110-residue scaffold (109standard
  residues plus protonation-labelledHIE), oneLa and7302waters. Its frozen test
  failed the required La/Dy coordination discrimination and weak-EF4 challenge.
- The four Hans Dy outer-water15ns walkers and later QM/MM correction artifacts
  remain available. Their reported outcomes are respectivelyNO_CALL and no final
  qualified production transfer; neither establishes whole-protein folding
  thermodynamics. Trajectory-file existence alone is not validation.
- Existing Hans monomer/dimer Protenix predictions failed their frozen geometry
  test. Recent Hans/Mex MACE+GFN2 work used44–50atom EF1–3 cores; its initial
  La/Dy signal failed transfer to the Dy-conditioned Hans crystal. That result
  is local-core evidence, not a failed global model.

[Assay and folding clarification](EVIDENCE.md) explains why a future global
representation must declare sequence/construct, assembly, metal/site occupancy,
waters, protonation and assay condition explicitly. No universal four-ion recipe
or new global simulation is selected by this inventory.

## Reproduction and integrity

[Commands](COMMANDS.md) reconstruct metadata only. Source seals and616unique IDs/
sequences,16reserved memberships,32exact full-chain matches, file joins and
18+1+597coverage accounting were verified. `VALIDATION.json` records these checks.
An initial scan was stopped after recognizing that the four-La run's`mmcif/`
directory is a general template cache; completedv2 inventories only its actual
`predict_out/`. Partialv1metadata is preserved. Neither version read reserved
outcomes. No originals or shared workflow files changed.
