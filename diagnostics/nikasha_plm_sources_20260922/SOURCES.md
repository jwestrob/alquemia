# Existing PLM sources for the Nikasha manuscript

**The existing joins cover all 176 scored-or-excluded proteins, 200 source genes
and 170 selected PQQ domains.** All 170 domains occur in the completed updated
4,028-tip tree. Existing RNA profiles cover 77 source genes (74 distinct full
proteins); the current within-genome screen covers 28 genes (27 full proteins)
in 25 named bins. These are source-coverage counts, not new biological results.

The integrated PLM manuscript introduces Nikasha. This inventory supplies existing
inputs for that delivery; it makes no new matrix, normalization, phylogeny,
correlation, substrate assignment, chemical score or cohort selection. Paths,
SHA256 hashes, schemas, key uniqueness and exact join checks are in
[INVENTORY.json](INVENTORY.json). All hash comparisons inspected here pass.

## Locations and join keys

Root **B** is `/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM`.
**P** is `B/revision_analysis/2026-09-11_PQQ_ADH`.
**I** is `P/energetics_queue/xoxf_all/inventory`.
**T** is `P/functional_reference_update/tree`.
The machine-readable inventory gives absolute paths for every table below.

| Existing source | Size | Reliable join / purpose |
|---|---:|---|
| `B/PLM_XoxF_energetics/results.tsv` | 176 proteins | `target_id`; 137 scored, 37 unsupported, 2 source-residue failures. Byte-identical to `P/energetics_queue/xoxf_all/final/results.tsv`. Older attempt tables are not the consolidated delivery. |
| `I/proteins.tsv` | 176 | `target_id` plus full `sequence_sha256`; source genes, domains, FASTA, retained scope. |
| `I/original_gene_provenance.tsv` | 200 | Original `gene_id` → `target_id`, original FASTA and verified protein hash. |
| `I/source_aliases.tsv` | 210 rows / 200 genes | Preserve `(target_id, gene_id, catalog, genome_id)`, exact protein/domain hashes, translation status and source provenance. Multiple catalog rows are not extra genes. |
| `I/selected_leaves.tsv` | 170 | Frozen `domain_id` ↔ `leaf_id`; use this frozen cohort rather than all current reference/curation rows. |
| `I/role_mappings.tsv` | 704 | `(sequence_id, role)` gives all four source-sequence/alignment-verified role positions for each protein. |
| `P/energetics_queue/coordination_review/curated_gene_sites.tsv` | 200 | `gene_id`, `unique_sequence_id`, `domain_id`, four-site motif and observed positions. Motif is not a metal/substrate label. |
| `T/domain_entity_crosswalk.tsv` | 7,795 | `(domain_id, entity_role, entity_id)` links exact full proteins to extracted domains; preserve domain boundaries and entity role. |
| `T/panel_tip_ledger.tsv`, `T/leaf_name_crosswalk.tsv` | 4,028 each | Exact `domain_id` → current tree identity/readable name. The original inference is `T/main.treefile`; current named delivery is `B/PLM_itol/PQQ_broad.treefile`. |
| `P/quantification/candidate_MAG_membership.tsv` | 15,349 broad candidate rows | `(candidate_id, catalog, genome_id)`, exact source/scaffold coordinates, translation evidence and genome-file hash. Join from aliases; do not parse a genome out of a gene or tree-tip name. |
| `B/PLM_XoxF_expression/genes.tsv` | 77 | `(reference, gene_id)`, domain/tree tip, full protein hash, available genome/taxonomy and explicit membership evidence. |
| `B/PLM_XoxF_expression/plotted_values.tsv`, `samples.tsv` | 1,925 cells / 25 samples | `(reference, gene_id, RNA_sample)`; retain assembly reference, original counts/lengths/denominators and sampling metadata. |
| `B/PLM_XoxF_correlations/within_genome_membership.tsv` | 28 | Exact source-CDS/genome checks, `reference`, `gene_id`, `genome_id`, full protein hash. |
| `.../within_genome_XoxF_profiles.tsv` | 700 | `(reference, genome_id, xoxF_gene, RNA_sample)`; original read-depth-normalized profiles for 28 × 25. These values are not genome-adjusted expression. |
| `.../within_genome_correlations.tsv.gz` | 40,426 pairs | `(reference, genome_id, xoxF_gene, partner_gene)`; existing raw, campaign-only and genome/campaign-adjusted correlations. Methods and target eligibility remain alongside. |

The broad full-protein alias source at
`P/energetics_queue/candidates/exact_full_protein_source_aliases.tsv` is useful
for tracing discovery, but its 16,912 rows do not define this frozen cohort.
The current curation source table has 210 rows but the same **200 unique genes**;
the difference is repeated catalog ownership, not an enlarged protein panel.
Likewise, expression contains 77 genes but only 73 tree domains and 74 full
proteins. Exact-domain aliases do not make identical full proteins or independent
biological observations.

## Three distinct existing transcript measures

1. **Sample-wide TPM:** `P/quantification/candidate_source_MAG_RNA.tsv` stores
   `samplewide_TPM`, `rna_readname_gene_count`, source CDS length, library and
   availability. Its historical arithmetic is gene RPK divided by the
   all-reference RPK sum × 10⁶, as recorded in `join_provenance.py` and
   `B/revision_analysis/2026-09-11_ureolysis_ammonia/library_audit.tsv`.
   This source-library observation is not the 25-sample expression matrix.
   The 29 older `B/Counts_Files/normalized/*.tpm.tsv` files are listed separately;
   they are not genome-normalized merely because their directory says normalized.

2. **Expression figure / correlation input:** the 77 × 25 export uses
   `unique_QNAME_gene_assignments × 10⁹ / (CDS_nt × input_RNA_read_ends)`.
   These units are **not TPM**. All values were copied from existing exports;
   none were recalculated. The 123 remaining cohort genes have no exported
   profiles and remain missing, not zero. Independent assembly mappings must
   not be summed into a community total.

3. **Genome/campaign adjustment:** the completed September 15 within-genome
   screen ranks background gene profiles across samples, forms a pair-specific
   shared-genome transcription score excluding the tested two genes, and
   residualizes both tested rank vectors on campaign and that score. Its result
   is a partial Spearman correlation, not DNA-normalized or per-cell expression.
   The earlier `genome_adjusted_profiles.tsv` has four selected pairs × 25
   samples and explicit residuals; it is preserved history, superseded in scope
   by the 28-target screen. A single general genome-normalized expression matrix
   for all 176 proteins was **not located in the reviewed exports**.

`P/integration/PQQ_domain_site_RNA_evidence.tsv` also retains a DNA recruitment
descriptor per scaffold kilobase per million input pairs. Keep its declared
`DNA_metric_scope`; it is not an existing gene-level RNA/DNA ratio or permission
to invent one. Historical ambiguous read assignments and missing BAMs limit
gene-specific recruitment claims, as the existing expression report explains.

## Practical assembly order and limits

Start from the frozen 176 `target_id` rows, retain each result status, then expand
through the **200 original source genes** only when attaching source-specific
genome/RNA evidence. Verify full protein hashes; join motif/domain positions by
the recorded IDs. Attach current tree names through exact `domain_id`. Attach
RNA by original gene **and mapping reference**, followed by `RNA_sample`. Keep
catalog membership, sequence-verified genome ownership and archived taxonomy
as separate evidence; 24 cohort genes have explicit sequence-audited ownership
in the frozen aliases, while the later within-genome study verifies its own 28
targets against additional archived bins.

All 77 expression genes and all 28 within-genome targets join the frozen source
genes without protein-hash disagreement. No missing profile or unsupported
chemical score should remove a cohort member from the denominator. PLM bands
are predictions; tree proximity and motifs do not turn them into affinity,
physiological occupancy or biochemical substrate labels.

The old `P/phylogeny/nearest_experimental_anchors.tsv` belongs to the original
3,998-tip tree. Its distances are not recomputed for the 4,028-tip update here.
Use the newer tree/crosswalk and its source reference-evidence table for current
placement, with any new distance/biological analysis planned separately. A
manuscript methods outline follows protocol selection, rather than freezing a
research challenger prematurely.

## Repeat this source inventory

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/nikasha_plm_sources_20260922/inventory_sources.py \
  --plm-root /groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM \
  --output /tmp/nikasha_plm_inventory_recheck.json
```

This reads existing files and reports schemas, hashes and identifier coverage;
it does not launch biological, structural or molecular calculations. No source
data were edited. Authorization is the September 22 delegated handoff recorded
in [the parent plan](../nikasha_shared_pool_20260922/PLAN.md).
