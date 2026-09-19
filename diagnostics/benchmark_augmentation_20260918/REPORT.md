# Delegated benchmark extension: one ready PQQ structural control

## Delivered

**8GY2/O05542 is prepared for the existing DFT protocol:** one additional
Ca-associated PQQ enzyme, two La/Ca endpoint inputs. No energy calculation,
MACE inference, fold, classifier fit or cluster job ran. Baseline/default and
earlier records remain unchanged. No new protocol was introduced.

The curation acquired 19 coordinate sets from a 51-entry deposited-PQQ search,
crosschecked their sequences against the current panel, and retained 13 earlier
gated candidates. These counts are inventories, not independent labelled tests.
No additional direct-affinity group or matched DFT/MACE row became available.
Acquisition selected unfamiliar enzyme/accession deposits and explicit duplicate
or unsupported-metal checks from entry titles/citations, before any scoring.
This is a bounded deposited-structure inventory, not a complete literature census.

## Ready control and its actual meaning

Gluconobacter oxydans membrane alcohol dehydrogenase, **8GY2 author chain A**,
contains deposited PQQ802/Ca803. It supplies a **Ca-associated enzyme/structural
transfer control**, not a measured La/Ca affinity direction or proof that La
cannot support activity. The deposition links the biochemical ethanol-oxidation
study [Adachi et al., 2023](https://doi.org/10.1021/acscatal.3c01962).
Coordinate and citation evidence came directly from the [deposition](https://www.rcsb.org/structure/8GY2);
the publisher full text returned HTTP 403 and was not independently reread here.
This evidence qualification is carried into the machine manifest.

Maximum identity to current source sequences is **43.7068%** under the frozen
global-alignment/longer-sequence rule. It therefore adds a group at the existing
50% cutoff; this is not an independent-fold claim. Two newer membrane-region
variant deposits, 9X0Q/9X0R, have identical catalytic-chain coordinate sequences
and remain the same group. No score was generated or inspected. Absence of this
identifier from the inspected diagnostic records is not an exhaustive historical
blindness audit.

The unchanged `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3` preparation
passes: **Glu215, Asn297, Asp342, nonacidic Thr344 and Lys369** map to its frozen
roles. The actual deposited assembly is a heme-bearing heterotrimer; the fixed
QM core retains its established selected-chain/core policy. Heme and other
chains are explicitly excluded by that existing carve policy. **Do not reuse
this normalized chain as a heme-free whole-protein MACE preparation.** A matched
MACE input remains unsupported, not silently substituted.

Preparation: 71 atoms/endpoint; Ca charge −2, La −1; singlets; 7 typed contacts;
no crystal/synthetic water. All 42 retained fragment heavy atoms match the
source exactly (maximum displacement 0 Å); no missing heavy atoms were rebuilt.
Nonmetal coordinates are byte-identical between endpoints. No geometry search,
metal relocation, cofactor substitution or altered decision bands occurred.

## Useful gated additions

| Candidate | What was recovered | Remaining decision/gate |
|---|---|---|
| **4CVB/Q93RE9** | Published Ca–PQQ enzyme structure; only 28.95% maximum identity to current cases | Its GLN220/ASP333/GLU335/TYR476 donor arrangement differs from the canonical core. Existing fixed-core preparation cannot be inherited. |
| **7WMK DepA** | PQQ/Ca complex and linked functional enzyme study | Noncanonical preparation and exact metal-use evidence crosswalk; not a direct-affinity label. |
| **4MH1 L-sorbose dehydrogenase** | Ca/PQQ coordinates and primary citation | Noncanonical site and full assay/metal-state crosswalk. |
| **3DAS, 3A9H** | Additional sugar-dehydrogenase structures | Different core/water chemistry; 3A9H has no metal within 4 Å of its first PQQ. No invented metal placement. |
| **9M2J/9M2K** | Additional deposited PQQ structures | Citation remains “to be published”; no verified functional label. |

4CVB is especially useful as a future *composition challenge*: an acidic residue
occurs at catalytic Asp+2 despite its observed Ca-associated state. The primary
study explicitly describes the different ligand arrangement. It does **not**
measure La/Ca affinity. [Rozeboom et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4815231/).
4CVB, 7WMK, 4MH1 and 9M2J/K connect under the frozen 50% rule and must not become
five independent family observations. 7WMK is [Yang et al., 2022](https://doi.org/10.1021/acs.jafc.2c01083).

1KV9/1YIQ/6DAM/6ZCV/6ZCW/7O6Z overlap existing biological/homology groups.
2D0V is 90.2% identical to an existing case. 5XM3 contains **Mg**, not deposited
Ca. These findings prevent an inflated new-control denominator.

Earlier seven functional PQQ sequence candidates still lack verified holo
structures. B2/C5 remain apo-model/cofactor-placement gates. PqqT/aqualysin
remain assay-to-site gates; the new PqqT Y161W deposition 9OOZ contains no metal
and does not resolve the issue. Existing records were reused, not relabelled.

## Checks, cost and handoff

Six read-only real-artifact checks passed: acquisition hashes, existing runner
task loading/hashes, preserved source coordinates, paired atom/charge invariants,
evidence-stratum separation and explicit absence of new scoring. These are
preparation/integrity checks, not energy integration or accuracy results.

The actual single preparation took **10.863996157422662 s wall /
10.759353711 CPU s**. Downloads/curation were not instrumented as an aggregate
compute benchmark; no GPU or Slurm allocation. No timing estimate is a result.

Full expanded manifest and prepared products are pinned in [RELEASE.json](RELEASE.json).
[READINESS.tsv](READINESS.tsv) lists every acquired structure and exclusion.
[COMMANDS.md](COMMANDS.md) gives exact preparation and verification operations.
The requested vault note is
`agent-captures/2026-09-18_laca-benchmark-augmentation.md`.
The existing manifest runner can consume the two prepared tasks when energy
evaluation is included in an agreed analysis. This curation did not execute it.

**Practical next handoff:** 8GY2 is ready for a DFT association-transfer test;
4CVB is a concrete, already deposited noncanonical chemistry challenge for a
separately agreed representation policy. Neither completes the missing balanced,
independent direct-affinity panel, and neither demonstrates an accuracy gain.
