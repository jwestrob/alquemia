# PLM evidence export ready — 22 September 2026

The existing PLM scan is now joined to its gene, genome, motif, tree and
transcript evidence in a reproducible export. **All 176 proteins remain present**,
including 39 unscored proteins. No scientific score or transcript value changed.

Start with
`workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv`.
`EXPORT.json` records the table schemas, source/output hashes, original method,
reference gauge, decision bands and transcript definitions. The complete export
is approximately 7.3 MB. [Commands](COMMANDS.md) reproduce it in a new directory.

## Actual contents

| Table | Rows | Meaning |
|---|---:|---|
| `proteins.tsv` | 176 | Original native DFT outcome plus exact sequence, gene/genome evidence IDs, current domains/tree names, motifs and RNA availability. |
| `genes.tsv` | 200 | Original source provenance, motif and explicit profile availability. |
| `gene_aliases.tsv` | 210 | Original exact-protein/domain aliases and catalog records, copied without deduplication across evidence types. |
| `protein_domains.tsv` | 176 | Protein-to-domain links covering 170 unique domains in the current 4,028-tip tree. |
| `protein_roles.tsv` | 704 | Four original aligned role observations per protein. |
| `gene_genomes.tsv` | 97 | Separate frozen catalog/translation, expression metadata and later exact-CDS genome evidence; full original evidence rows retained as JSON. |
| `source_library_rna.tsv` | 200 | Original source-library RNA status, counts and samplewide TPM. |
| `expression_genes.tsv`, `rna_samples.tsv` | 77 / 25 | Original gene/reference and RNA-sample metadata. |
| `rna_profiles.tsv` | 1,925 | All original 77 × 25 read-depth-normalized expression cells. |
| `within_genome_membership.tsv`, `within_genome_profiles.tsv` | 28 / 700 | Existing verified target membership and original expression profiles. |
| `within_genome_correlations.tsv.gz` | 40,426 | Existing pair-specific genome/campaign-adjusted coefficients and their original columns. |
| `historical_genome_adjusted_pair_profiles.tsv` | 100 | The earlier four-pair residual-profile records, with their distinct original method. |

Protein-level list columns use JSON arrays. They summarize child-table links;
they do not collapse source genes or upgrade weak genome evidence. Join RNA on
`reference`, `gene_id` and `RNA_sample`, retaining full-protein `target_id` as the
parent key. Source gene and tree domain are different units.

## Method and missingness preserved

The score method is the existing
`pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`: native ORCA
r2SCAN-3c/CPCM(Water), fixed core. It is **not** the released fast scorer or a
research accommodation/pool result. Each original result field is checked
against the authoritative consolidated JSON and copied as its original string.
The existing 137 scored / 37 unsupported / 2 source-residue-failed outcomes remain
unchanged. Calibration and reporting-gauge caveats are copied in `EXPORT.json`.

The 123 genes without existing sample profiles have explicit missing status and
blank profile fields. No 25-sample zero vectors were manufactured. The separate
source-library RNA table has 124 unavailable records: that measurement describes
the gene's own source library, while cross-sample mapping profiles have a
different coverage rule. These counts need not agree.

Source-library samplewide TPM, read-end/CDS-length-normalized expression, and
genome/campaign-adjusted correlations remain distinct. `within_genome_profiles`
contains the original expression values, not a general genome-normalized
abundance. Historical residuals remain a separately named pair-specific table.
There is no new general genome-normalized matrix or candidate-score overlay.

Unknown PLM metal labels remain unknown. This export adds usable joins, not
affinity validation, substrate assignments, occupancy estimates or independent
biological replicates. No tree distances, expression normalization, correlations,
geometry or molecular energies were recomputed.

## Verification and execution

**Nine real-fixture tests pass, zero skips** ([log](TESTS_v2.txt)). They check all
176 outcomes, exact transcript strings, source TPM/availability, all source
genome evidence, current tree IDs and output preservation. Explicitly corrupted
copies of real records test hash failure, wrong gene identity, missing samples,
changed scores and duplicate cohort rows. Original sources remain untouched.

The actual CLI export completed in **6.110343515872955 seconds wall time**, with
zero molecular calls and no cluster allocation. All 14 output table hashes and
row counts were independently re-read and checked; [RESULT.json](RESULT.json)
records the result. The initial test invocation in the scientific CPU environment
could not run because it lacks pytest ([original log](TESTS.txt)); the existing
normal Python environment ran the passing suite. No dependency was installed.

The input inventory and source hashes fail closed if an input changes. Existing
output directories cannot be overwritten. [PLAN.md](PLAN.md) records the approved
read-only scope; the vault note records this handoff. Shared scoring code,
production data and existing scientific outputs were not modified.
