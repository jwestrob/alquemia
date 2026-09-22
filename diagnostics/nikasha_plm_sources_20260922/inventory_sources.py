#!/usr/bin/env python3
"""Inventory existing PLM source keys; do not recompute scientific values."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


P = "revision_analysis/2026-09-11_PQQ_ADH/"
I = P + "energetics_queue/xoxf_all/inventory/"
T = P + "functional_reference_update/tree/"
SOURCES = {
    "energetics": ("PLM_XoxF_energetics/results.tsv", ["target_id"]),
    "energetics_origin": (P + "energetics_queue/xoxf_all/final/results.tsv", ["target_id"]),
    "proteins": (I + "proteins.tsv", ["target_id"]),
    "aliases": (I + "source_aliases.tsv", ["target_id", "gene_id", "catalog", "genome_id"]),
    "original_genes": (I + "original_gene_provenance.tsv", ["gene_id"]),
    "selected_leaves": (I + "selected_leaves.tsv", ["domain_id"]),
    "roles": (I + "role_mappings.tsv", ["sequence_id", "role"]),
    "motifs": (P + "energetics_queue/coordination_review/curated_gene_sites.tsv", ["gene_id"]),
    "all_full_protein_aliases": (P + "energetics_queue/candidates/exact_full_protein_source_aliases.tsv", ["target_id", "gene_id"]),
    "current_curation_genes": (P + "curation/xoxF/source_genes.tsv", ["gene_id"]),
    "current_curation_leaves": (P + "curation/xoxF/selected_leaves.tsv", ["domain_id"]),
    "tree_tip_ledger": (T + "panel_tip_ledger.tsv", ["domain_id"]),
    "tree_entities": (T + "domain_entity_crosswalk.tsv", ["domain_id", "entity_role", "entity_id"]),
    "tree_names": (T + "leaf_name_crosswalk.tsv", ["original_leaf_id"]),
    "download_tree_names": ("PLM_itol/PQQ_broad.Leaf_names.tsv", ["original_leaf_id"]),
    "old_nearest_anchors": (P + "phylogeny/nearest_experimental_anchors.tsv", ["domain_id"]),
    "MAG_membership": (P + "quantification/candidate_MAG_membership.tsv", ["candidate_id", "catalog", "genome_id"]),
    "RNA_MAG_join": (P + "quantification/candidate_source_MAG_RNA.tsv", ["candidate_id"]),
    "domain_RNA_DNA": (P + "integration/PQQ_domain_site_RNA_evidence.tsv", ["candidate_id", "domain_id"]),
    "expression_genes": ("PLM_XoxF_expression/genes.tsv", ["reference", "gene_id"]),
    "expression_samples": ("PLM_XoxF_expression/samples.tsv", ["RNA_sample"]),
    "expression_values": ("PLM_XoxF_expression/plotted_values.tsv", ["reference", "gene_id", "RNA_sample"]),
    "genome_targets": ("PLM_XoxF_correlations/within_genome_targets.tsv", ["reference", "genome_id", "xoxF_gene"]),
    "genome_membership": ("PLM_XoxF_correlations/within_genome_membership.tsv", ["reference", "gene_id"]),
    "genome_profiles": ("PLM_XoxF_correlations/within_genome_XoxF_profiles.tsv", ["reference", "genome_id", "xoxF_gene", "RNA_sample"]),
    "genome_correlations": ("PLM_XoxF_correlations/within_genome_correlations.tsv.gz", ["reference", "genome_id", "xoxF_gene", "partner_gene"]),
    "old_genome_profiles": ("PLM_XoxF_correlations/genome_adjusted_profiles.tsv", ["reference", "genome_id", "xoxF_gene", "partner_gene", "RNA_sample"]),
}
DOCUMENTS = [
    I + "README.md", I + "manifest.json",
    "PLM_XoxF_energetics/README.md", "PLM_XoxF_expression/README.md",
    "PLM_XoxF_expression/provenance.json", "PLM_XoxF_expression/validation.json",
    "PLM_XoxF_correlations/within_genome_report.md",
    "PLM_XoxF_correlations/within_genome_methods.json",
    "PLM_XoxF_correlations/within_genome_validation.json",
    "PLM_XoxF_correlations/genome_adjusted_methods.json",
    T + "README.md", T + "export_validation.json", T + "main.treefile",
    "PLM_itol/PQQ_broad.treefile", "PLM_itol/PQQ_broad.Reference_evidence.tsv",
    P + "quantification/join_provenance.py",
    "revision_analysis/2026-09-11_ureolysis_ammonia/library_audit.tsv",
]


def pin(path):
    if not path.is_file():
        return {"path": str(path), "status": "missing"}
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return {"path": str(path), "status": "present", "sha256": h.hexdigest(), "bytes": path.stat().st_size}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plm-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    tables, records = {}, {}
    for name, (relative, keys) in SOURCES.items():
        path = args.plm_root / relative
        record = pin(path)
        if record["status"] == "present":
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                fields = reader.fieldnames
            tables[name] = rows
            key_counts = Counter(tuple(row.get(key) for key in keys) for row in rows)
            record.update(row_count=len(rows), columns=fields, inspected_key=keys,
                          missing_key_columns=sorted(set(keys)-set(fields)),
                          distinct_inspected_keys=len(key_counts),
                          duplicated_inspected_key_groups=sum(n > 1 for n in key_counts.values()))
        records[name] = record

    def values(name, key):
        return {row[key] for row in tables[name] if row[key]}

    targets = values("proteins", "target_id")
    genes = values("original_genes", "gene_id")
    domains = values("selected_leaves", "domain_id")
    proteins = {r["target_id"]: r for r in tables["proteins"]}
    original = {r["gene_id"]: r for r in tables["original_genes"]}
    aliases = tables["aliases"]
    expression = tables["expression_genes"]
    genome = tables["genome_membership"]
    counts = {
        "frozen_full_protein_count": len(targets),
        "frozen_original_gene_count": len(genes),
        "frozen_selected_domain_count": len(domains),
        "energetics_status_counts": dict(Counter(r["status"] for r in tables["energetics"])),
        "energetics_exact_target_set_match": values("energetics", "target_id") == targets,
        "download_origin_table_bytes_equal": records["energetics"]["sha256"] == records["energetics_origin"]["sha256"],
        "frozen_alias_rows": len(aliases),
        "frozen_alias_distinct_genes": len({r["gene_id"] for r in aliases}),
        "frozen_alias_targets_not_in_cohort": sorted({r["target_id"] for r in aliases}-targets),
        "frozen_alias_protein_hash_mismatches": [r["gene_id"] for r in aliases if r["protein_sha256"] != proteins[r["target_id"]]["sequence_sha256"]],
        "current_curation_unique_genes": len(values("current_curation_genes", "gene_id")),
        "current_curation_genes_not_in_frozen_cohort": sorted(values("current_curation_genes", "gene_id")-genes),
        "genes_with_motif_rows": len(genes & values("motifs", "gene_id")),
        "proteins_with_four_declared_roles": sum(n == 4 for n in Counter(r["sequence_id"] for r in tables["roles"]).values()),
        "frozen_domains_in_updated_tree": len(domains & values("tree_names", "original_leaf_id")),
        "frozen_domains_missing_updated_tree": sorted(domains-values("tree_names", "original_leaf_id")),
        "updated_tree_name_rows": len(tables["tree_names"]),
        "frozen_proteins_in_updated_tree_entities": len(targets & {r["entity_id"] for r in tables["tree_entities"] if r["entity_role"] == "PLM"}),
        "expression_source_genes": len(expression),
        "expression_genes_in_frozen_cohort": len({r["gene_id"] for r in expression} & genes),
        "cohort_genes_without_exported_expression_profiles": len(genes-values("expression_genes", "gene_id")),
        "expression_distinct_proteins": len({original[r["gene_id"]]["target_id"] for r in expression if r["gene_id"] in original}),
        "expression_protein_hash_mismatches": [r["gene_id"] for r in expression if r["gene_id"] in original and r["source_protein_sha256"] != original[r["gene_id"]]["normalized_sequence_sha256"]],
        "expression_cells": len(tables["expression_values"]),
        "expression_sample_count": len(tables["expression_samples"]),
        "expression_domains_in_updated_tree": len(values("expression_genes", "domain_id") & values("tree_names", "original_leaf_id")),
        "within_genome_target_genes": len(genome),
        "within_genome_targets_in_frozen_cohort": len(values("genome_membership", "gene_id") & genes),
        "within_genome_named_bins": len(values("genome_membership", "genome_id")),
        "within_genome_distinct_proteins": len({original[r["gene_id"]]["target_id"] for r in genome if r["gene_id"] in original}),
        "within_genome_protein_hash_mismatches": [r["gene_id"] for r in genome if r["gene_id"] in original and r["protein_sha256"] != original[r["gene_id"]]["normalized_sequence_sha256"]],
        "within_genome_profile_rows": len(tables["genome_profiles"]),
        "within_genome_pair_rows": len(tables["genome_correlations"]),
    }
    by_gene = defaultdict(list)
    for row in aliases:
        by_gene[row["gene_id"]].append(row)
    counts["frozen_genes_with_declared_sequence_audited_genome_alias"] = sum(any(r["genome_id"] and r["translation_status"] in {"exact_CDS_translation", "identical_body_initiator_methionine_only"} for r in by_gene[g]) for g in genes)
    counts["frozen_alias_translation_status_counts"] = dict(Counter(r["translation_status"] for r in aliases))
    result = {
        "schema": "nikasha_plm_existing_source_inventory_v1", "date": "2026-09-22",
        "scope": "Read-only source schemas, exact identifier/hash joins and coverage; no new matrix, scientific score, normalization, correlation or biological label.",
        "plm_root": str(args.plm_root.resolve()), "sources": records,
        "documents": [pin(args.plm_root / rel) for rel in DOCUMENTS],
        "archived_TPM_files": [str(p) for p in sorted((args.plm_root / "Counts_Files/normalized").glob("*.tpm.tsv"))],
        "join_checks": counts,
        "unavailable_products": [
            "No cohort-wide 176-protein genome-normalized expression matrix located in the reviewed exports.",
            "Within-genome adjustment is a pair-specific correlation covariate; raw normalized_expression is not a genome-adjusted abundance.",
            "No new phylogenetic distances were computed for the 4,028-tip update; old nearest-anchor distances retain their 3,998-tip origin."
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
