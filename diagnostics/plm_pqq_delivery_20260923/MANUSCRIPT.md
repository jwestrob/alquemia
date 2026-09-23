# Nikasha in the integrated PLM paper — initial manuscript material

**Status, 23 September 2026:** completed method-development results; two PLM
proteins have been run with the practical accommodated candidate. The full PLM
candidate scan remains to be performed. The original DFT results and biological
joins are retained. The following text separates those completed findings from
the application section that will be written after the scan.

## Suggested methods text

We developed Nikasha to assess the compatibility of PQQ-containing protein
structures with Ca and La coordination chemistry. The accommodated protocol
evaluates three predeclared La-conditioned structural predictions per protein.
Each structure is prepared using the same PQQ state and homologous donor-role
mapping. A common set of complete polar fragments is retained across the three
structures, using a 4.3 Å contact envelope. Within each structure, Ca and La
share initial nonmetal coordinates and chemical composition, with the total
charge adjusted for metal substitution.

A pretrained MACE-OMOL potential supplies vacuum energies and forces. We select
up to four independent, source-connected donor torsions using normalized
individual-metal and differential force projections. Each metal independently
generates a bounded accommodation proposal. Both metals are then evaluated at
the same candidate geometries: the original structure and the Ca- and La-derived
proposals. Selection uses a composite energy consisting of the MACE vacuum
energy plus the difference between native GFN2-xTB energies in ALPB water and
vacuum. This solvent difference is an approximate correction evaluated with
GFN2-xTB, not a self-consistent reaction field of the MACE model. The MACE model
was not fine-tuned on PLM metal labels.

The Ca-minus-La contrast is calculated from the operationally selected energies,
retaining the origin when an improvement is below 0.1 model kcal mol−1. We
report its median across all three structures, individual contrasts, structural
spread, and protocol-specific Ca-supported, La-supported or inconclusive calls.
Decision bands are fixed using the designated 25 canonical reference structures.
Incomplete preparations or required energy calculations yield unavailable
results. These contrasts are electronic-energy descriptors for classification;
they are not binding free energies, dissociation constants or probabilities of
physiological metal occupancy.

## Suggested results text: completed development evidence

Bounded accommodation improved discrimination across previously examined
alternative structures of known-class PQQ proteins. In the most extensively
tested context policy, which fixes pocket membership across ten reference folds,
207 of 208 scorable alternative structures received the expected class, one was
inconclusive and none received the opposing class; 17 of 225 attempted structures
were unavailable. On the same supported structures, a static calculation using
the same pocket membership and numerical solver gave 199 correct calls, one
opposing call and eight inconclusive calls. This comparison supports a useful
contribution from local structural response beyond changing pocket membership
or numerical precision. Structural replicas represent 25 protein groups and
were used during development; they are not independent prospective biological
tests.

The practical three-structure protocol retained the expected median class for
all 91 complete predefined reference triples in its primary evaluation. Six
triples retained prior preparation exclusions, and three were unavailable
because they shared one nonconverged energy calculation. A separately reported
two-start recovery of that calculation restored all three expected calls without
changing geometry or decision bands. Individual structures could remain
inconclusive, and accommodation did not consistently narrow structural spread
across the reference population.

In two PLM development examples, accommodation relieved compressed metal–Asp
contacts and reduced the ranges across the three structures from 40.14 to 2.31
and from 72.42 to 17.23 model kcal mol−1. Both median predictions remained
Ca-supported; their experimental metal preferences are unknown. These examples
demonstrate a response to input geometry and greater consistency for those
inputs, without treating an expected XoxF annotation as a biochemical label.

## Application paragraph to complete after the cohort scan

Do not insert an invented cohort count or imply that the two examples constitute
the completed PLM scan. The final paragraph should report:

1. The actual declared PQQ cohort and available three-structure groups, including
   which XoxF/Exa/Ped/ADH-like groups were represented and why others were absent.
2. Numbers Ca-supported, La-supported, inconclusive and unavailable; distinguish
   preparation exclusions from numerical failures and report discordant folds.
3. The actual joined motif and phylogenetic distribution, preserving unresolved
   functional annotations and the original DFT calls separately.
4. Associations with the existing sample- and genome-context transcript measures,
   using their original normalization definitions and missingness. Do not label
   the existing read-depth-normalized profiles TPM or infer a replicated depth
   effect from confounded sampling.

No substrate assignment, new environmental contrast or transcriptional
normalization is proposed by this packaging task.

## Suggested figure / table placement

- **Main text:** a reference-decision panel with proteins as groups; a compact
  before/after geometry example; the final PQQ phylogeny linked to actual Nikasha
  predictions and measured transcript profiles.
- **Supplement:** source-level result table, exclusion ledger, canonical bands,
  static-versus-accommodated comparison and all correlated structural replicas.
- Existing editable geometry-response figure:
  `workspaces/plm_envelope_response_20260923/figures_v2/` (SVG/PDF and preview).
  Its caption must identify the two proteins as unlabeled development examples.

## Sources and reproducibility

Model provenance: [official MACE-OMOL release](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_omol_0).
Native solvent implementation: [ORCA 6.1 semiempirical-method documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb).
These sources describe the components; the Nikasha discrimination and accommodation
claims derive from the actual project results linked below.

[Exact settings and evidence](SUPPLEMENT.md) · [PLM SOP](../../docs/PLM_PQQ_SOP.md) ·
[strongest transfer](../strict_native_transfer_20260923/REPORT.md) ·
[matched static ablation](../strict_static_ablation_20260923/REPORT.md) ·
[practical transfer](../motion_envelope_transfer_20260923/REPORT.md) ·
[separate recovery](../native_failed_cell_recovery_20260923/POOLED_SENSITIVITY.md) ·
[actual PLM response](../plm_envelope_response_20260923/REPORT.md).
