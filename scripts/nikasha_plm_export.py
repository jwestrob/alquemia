#!/usr/bin/env python3
"""Export existing PLM evidence through exact, hash-pinned relational joins."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


SCHEMA = "nikasha_plm_relational_export_v1"
INPUT_SCHEMA = "nikasha_plm_export_inputs_v1"
DFT_PROTOCOL = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
REQUIRED = (
    "energetics", "proteins", "original_genes", "aliases", "selected_leaves",
    "roles", "motifs", "tree_names", "expression_genes", "expression_samples",
    "expression_values", "genome_membership", "genome_profiles",
    "genome_correlations", "old_genome_profiles",
)
DOC_SUFFIXES = {
    "expression_methods": "PLM_XoxF_expression/README.md",
    "genome_methods": "PLM_XoxF_correlations/within_genome_methods.json",
    "historical_genome_methods": "PLM_XoxF_correlations/genome_adjusted_methods.json",
}


class ExportError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ExportError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_pin(pin):
    path = Path(pin["path"])
    require(path.is_absolute(), f"Input path must be absolute: {path}")
    raw = path.read_bytes()
    require(digest(raw) == pin["sha256"], f"Input SHA256 mismatch: {path}")
    return raw


def parse_tsv(raw, compressed=False):
    if compressed:
        raw = gzip.decompress(raw)
    reader = csv.DictReader(io.StringIO(raw.decode()), delimiter="\t")
    rows = list(reader)
    require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames), "Missing or duplicate TSV header")
    require(all(None not in row and all(v is not None for v in row.values()) for row in rows), "Malformed TSV row width")
    return {"fields": reader.fieldnames, "rows": rows}


def unique(rows, key, name):
    out = {}
    for row in rows:
        value = row[key]
        require(value and value not in out, f"Duplicate or empty {name} {key}: {value}")
        out[value] = row
    return out


def keyed(rows, keys, name):
    seen = set()
    for row in rows:
        key = tuple(row[k] for k in keys)
        require(all(key) and key not in seen, f"Duplicate or empty {name} key: {key}")
        seen.add(key)
    return seen


def only(rows, key, name):
    values = {row[key] for row in rows}
    require(len(values) == 1, f"Conflicting {name} field {key}: {sorted(values)}")
    return next(iter(values))


def jarray(values):
    return json.dumps(sorted(set(values)), separators=(",", ":"))


def load_inputs(inventory):
    inventory = Path(inventory).resolve()
    raw = inventory.read_bytes()
    config = json.loads(raw)
    require(config.get("schema") == INPUT_SCHEMA, "Unsupported export input inventory")
    source = json.loads(read_pin(config["source_inventory"]))
    require(source.get("schema") == "nikasha_plm_existing_source_inventory_v1", "Unsupported source inventory")
    tables, pins = {}, {"source_inventory": config["source_inventory"], "score_record": config["score_record"]}
    for name in REQUIRED:
        pin = source["sources"][name]
        require(pin.get("status") == "present", f"Unavailable input: {name}")
        tables[name] = parse_tsv(read_pin(pin), pin["path"].endswith(".gz"))
        require(len(tables[name]["rows"]) == pin["row_count"], f"Row count changed: {name}")
        pins[name] = pin
    documents = {}
    for name, suffix in DOC_SUFFIXES.items():
        matches = [p for p in source["documents"] if p["path"].endswith("/" + suffix)]
        require(len(matches) == 1, f"Document pin missing or ambiguous: {suffix}")
        pins[name] = matches[0]
        documents[name] = read_pin(matches[0]).decode()
    score = json.loads(read_pin(config["score_record"]))
    return {"config": config, "tables": tables, "score": score, "documents": documents,
            "pins": pins, "inventory": {"path": str(inventory), "sha256": digest(raw)}}


def build_export(data):
    """Copy values; validate exact identities. No new score or normalization."""
    tables, score = data["tables"], data["score"]
    rows = {name: table["rows"] for name, table in tables.items()}
    proteins = unique(rows["proteins"], "target_id", "protein")
    original = unique(rows["original_genes"], "gene_id", "original gene")
    energies = unique(rows["energetics"], "target_id", "energy result")
    score_rows = unique(score["per_target"], "target_id", "JSON energy result")
    require(score.get("protocol_id") == DFT_PROTOCOL, "Expected the archived native r2SCAN-3c/CPCM fixed-core scan")
    require(set(energies) == set(proteins) == set(score_rows), "Protein/score cohort mismatch")
    require(len(proteins) == data["config"]["expected_proteins"] == score["total_targets"], "Protein denominator mismatch")
    require(len(original) == data["config"]["expected_genes"], "Gene denominator mismatch")
    for target, row in energies.items():
        for key, value in row.items():
            expected = score_rows[target][key]
            require(value == ("" if expected is None else str(expected)), f"TSV/JSON result mismatch: {target}/{key}")
        numeric = [row[k] for k in ("E_La_hartree", "E_Ca_hartree", "R_hartree", "S_kcal_mol")]
        require(all(numeric) if row["status"] == "scored" else not any(numeric), f"Score availability mismatch: {target}")
    aliases = defaultdict(list)
    for row in rows["aliases"]:
        gene, target = row["gene_id"], row["target_id"]
        require(gene in original and target in proteins, f"Unknown alias identity: {gene}")
        require(original[gene]["target_id"] == target, f"Alias protein mismatch: {gene}")
        require(row["protein_sha256"] == proteins[target]["sequence_sha256"] == original[gene]["normalized_sequence_sha256"], f"Alias protein hash mismatch: {gene}")
        aliases[gene].append(row)
    require(set(aliases) == set(original), "Missing source gene aliases")
    motifs = unique(rows["motifs"], "gene_id", "motif gene")
    require(set(motifs) == set(original), "Motif gene cohort mismatch")
    tree = unique(rows["tree_names"], "original_leaf_id", "tree domain")
    selected_domains = {r["domain_id"] for r in rows["selected_leaves"]}
    require(len(selected_domains) == data["config"]["expected_domains"] and selected_domains <= set(tree), "Selected tree coverage mismatch")
    role_counts = Counter(row["sequence_id"] for row in rows["roles"])
    keyed(rows["roles"], ["sequence_id", "role"], "protein role")
    require(set(role_counts) == set(proteins) and set(role_counts.values()) == {4}, "Four-role protein coverage mismatch")
    expression = unique(rows["expression_genes"], "gene_id", "expression gene")
    genome_members = unique(rows["genome_membership"], "gene_id", "within-genome target")
    for name, records, hash_field in [("expression", expression, "source_protein_sha256"), ("within-genome", genome_members, "protein_sha256")]:
        require(set(records) <= set(original), f"Unknown {name} gene")
        for gene, row in records.items():
            require(row[hash_field] == original[gene]["normalized_sequence_sha256"], f"{name} protein hash mismatch: {gene}")
            require(row["domain_id"] == only(aliases[gene], "domain_id", gene), f"{name} domain mismatch: {gene}")
            require(row["tree_leaf_id"] == tree[row["domain_id"]]["new_leaf_id"], f"{name} current tree name mismatch: {gene}")
    sample_ids = set(unique(rows["expression_samples"], "RNA_sample", "RNA sample"))
    expr_keys = keyed(rows["expression_values"], ["reference", "gene_id", "RNA_sample"], "expression cell")
    expected_keys = {(r["reference"], g, s) for g, r in expression.items() for s in sample_ids}
    require(expr_keys == expected_keys, "Expression sample coverage mismatch")
    for row in rows["genome_profiles"]:
        gene = row["xoxF_gene"]
        require(gene in genome_members and row["reference"] == genome_members[gene]["reference"] and row["genome_id"] == genome_members[gene]["genome_id"], f"Within-genome profile membership mismatch: {gene}")
    genome_keys = keyed(rows["genome_profiles"], ["reference", "genome_id", "xoxF_gene", "RNA_sample"], "within-genome profile")
    require(genome_keys == {(r["reference"], r["genome_id"], g, s) for g, r in genome_members.items() for s in sample_ids}, "Within-genome sample coverage mismatch")
    keyed(rows["genome_correlations"], ["reference", "genome_id", "xoxF_gene", "partner_gene"], "within-genome pair")
    for row in rows["genome_correlations"] + rows["old_genome_profiles"]:
        gene = row["xoxF_gene"]
        require(gene in genome_members and row["reference"] == genome_members[gene]["reference"] and row["genome_id"] == genome_members[gene]["genome_id"], f"Correlation membership mismatch: {gene}")

    output = {}
    def add(name, fields, records, description, source_names):
        require(len(fields) == len(set(fields)), f"Duplicate output columns: {name}")
        require(all(set(r) == set(fields) for r in records), f"Output row schema mismatch: {name}")
        output[name] = {"fields": fields, "rows": records, "description": description, "sources": source_names}

    def copy(name, source, gene_field=None):
        fields = tables[source]["fields"][:]
        records = [dict(r) for r in rows[source]]
        if gene_field:
            require("target_id" not in fields, f"Ambiguous target column: {source}")
            fields = ["target_id"] + fields
            for row in records:
                row["target_id"] = original[row[gene_field]]["target_id"]
        add(name, fields, records, "Original strings copied; target_id is an exact source-gene join where added.", [source])

    for name, source, gene_field in [
        ("gene_aliases.tsv", "aliases", None), ("protein_roles.tsv", "roles", None),
        ("expression_genes.tsv", "expression_genes", "gene_id"),
        ("rna_samples.tsv", "expression_samples", None), ("rna_profiles.tsv", "expression_values", "gene_id"),
        ("within_genome_membership.tsv", "genome_membership", "gene_id"),
        ("within_genome_profiles.tsv", "genome_profiles", "xoxF_gene"),
        ("within_genome_correlations.tsv.gz", "genome_correlations", "xoxF_gene"),
        ("historical_genome_adjusted_pair_profiles.tsv", "old_genome_profiles", "xoxF_gene"),
    ]:
        copy(name, source, gene_field)

    gene_rows, source_rna, domain_rows, genome_rows = [], [], {}, []
    by_target = defaultdict(list)
    gene_fields = tables["original_genes"]["fields"] + ["library_id", "domain_id", "current_tree_leaf_id", "four_site_pattern", "Glu_site_position", "Asn_site_position", "catalytic_Asp_site_position", "additional_site_position", "all_four_observed", "motif_interpretation", "rna_profile_status", "rna_profile_reference", "rna_profile_sample_count", "within_genome_profile_status", "genome_evidence_rows"]
    source_rna_fields = ["target_id", "gene_id", "library_id", "source_RNA_status", "source_RNA_counts", "source_RNA_TPM", "measure_id"]
    genome_fields = ["target_id", "gene_id", "genome_id", "catalog", "evidence_source", "evidence_kind", "source_row_json"]
    for gene in sorted(original):
        base = original[gene]
        target = base["target_id"]
        by_target[target].append(gene)
        alias_rows, motif = aliases[gene], motifs[gene]
        domain = only(alias_rows, "domain_id", gene)
        require(domain in selected_domains and motif["domain_id"] == domain, f"Motif domain mismatch: {gene}")
        require(motif["protein_sha256"] == base["normalized_sequence_sha256"] and motif["unique_sequence_id"] == target, f"Motif protein mismatch: {gene}")
        library = only(alias_rows, "library_id", gene)
        for record in alias_rows:
            if record["genome_id"]:
                genome_rows.append(dict(target_id=target, gene_id=gene, genome_id=record["genome_id"], catalog=record["catalog"], evidence_source="aliases", evidence_kind="frozen_alias_catalog_and_translation_record", source_row_json=json.dumps(record, sort_keys=True, separators=(",", ":"))))
        if gene in expression and expression[gene]["taxonomy_genome_id"]:
            record = expression[gene]
            genome_rows.append(dict(target_id=target, gene_id=gene, genome_id=record["taxonomy_genome_id"], catalog="", evidence_source="expression_genes", evidence_kind=record["genome_membership_evidence"], source_row_json=json.dumps(record, sort_keys=True, separators=(",", ":"))))
        if gene in genome_members:
            record = genome_members[gene]
            genome_rows.append(dict(target_id=target, gene_id=gene, genome_id=record["genome_id"], catalog=record["source_genome_catalog"], evidence_source="genome_membership", evidence_kind=record["membership_status"], source_row_json=json.dumps(record, sort_keys=True, separators=(",", ":"))))
        source_rna.append(dict(target_id=target, gene_id=gene, library_id=library, **{k: only(alias_rows, k, gene) for k in ["source_RNA_status", "source_RNA_counts", "source_RNA_TPM"]}, measure_id="source_library_samplewide_TPM"))
        gene_rows.append(dict(base, library_id=library, domain_id=domain, current_tree_leaf_id=tree[domain]["new_leaf_id"], **{k: motif[k] for k in ["four_site_pattern", "Glu_site_position", "Asn_site_position", "catalytic_Asp_site_position", "additional_site_position", "all_four_observed"]}, motif_interpretation=motif["residue_interpretation"], rna_profile_status="present" if gene in expression else "missing_from_existing_export", rna_profile_reference=expression[gene]["reference"] if gene in expression else "", rna_profile_sample_count=str(len(sample_ids)) if gene in expression else "", within_genome_profile_status="present" if gene in genome_members else "outside_existing_screen", genome_evidence_rows=""))
        domain_rows[(target, domain)] = dict(target_id=target, domain_id=domain, current_tree_leaf_id=tree[domain]["new_leaf_id"], tree_naming_status=tree[domain]["naming_status"], source_gene_ids_json="")
    evidence_counts = Counter(r["gene_id"] for r in genome_rows)
    for row in gene_rows:
        row["genome_evidence_rows"] = str(evidence_counts[row["gene_id"]])
    for (target, domain), row in domain_rows.items():
        row["source_gene_ids_json"] = jarray(g for g in by_target[target] if only(aliases[g], "domain_id", g) == domain)
    add("genes.tsv", gene_fields, gene_rows, "All original genes; absent sample profiles remain explicitly missing.", ["original_genes", "aliases", "motifs", "tree_names", "expression_genes", "genome_membership"])
    add("source_library_rna.tsv", source_rna_fields, source_rna, "Original source-library TPM/count strings; unavailable source RNA remains blank.", ["aliases"])
    add("gene_genomes.tsv", genome_fields, genome_rows, "Evidence records, not unique biological genomes; catalog, historical assignment and exact-CDS evidence stay separate.", ["aliases", "expression_genes", "genome_membership"])
    add("protein_domains.tsv", ["target_id", "domain_id", "current_tree_leaf_id", "tree_naming_status", "source_gene_ids_json"], [domain_rows[k] for k in sorted(domain_rows)], "Exact domain links to the current 4,028-tip tree; no new placement or distances.", ["aliases", "tree_names"])

    extra_fields = ["sequence_sha256", "full_length_aa", "protocol_id", "score_method", "source_gene_ids_json", "genome_evidence_ids_json", "domain_ids_json", "current_tree_leaf_ids_json", "four_site_patterns_json", "source_gene_count", "rna_profile_source_gene_count", "rna_profile_missing_source_gene_count", "rna_profile_status", "within_genome_source_gene_count", "prediction_evidence_status"]
    protein_rows = []
    require(set(by_target) == set(proteins), "Protein without original source gene")
    for target in sorted(proteins):
        genes = by_target[target]
        require(set(energies[target]["source_gene_ids"].split(";")) == set(genes), f"Existing result source genes changed: {target}")
        domains = {only(aliases[g], "domain_id", g) for g in genes}
        present = sum(g in expression for g in genes)
        protein_rows.append(dict(energies[target], sequence_sha256=proteins[target]["sequence_sha256"], full_length_aa=proteins[target]["full_length_aa"], protocol_id=score["protocol_id"], score_method="native_ORCA_r2SCAN-3c_CPCM_water_fixed_core_DFT", source_gene_ids_json=jarray(genes), genome_evidence_ids_json=jarray(r["genome_id"] for r in genome_rows if r["target_id"] == target), domain_ids_json=jarray(domains), current_tree_leaf_ids_json=jarray(tree[d]["new_leaf_id"] for d in domains), four_site_patterns_json=jarray(motifs[g]["four_site_pattern"] for g in genes), source_gene_count=str(len(genes)), rna_profile_source_gene_count=str(present), rna_profile_missing_source_gene_count=str(len(genes)-present), rna_profile_status="all_source_genes_present" if present == len(genes) else "some_source_genes_missing" if present else "missing_from_existing_export", within_genome_source_gene_count=str(sum(g in genome_members for g in genes)), prediction_evidence_status="unknown_PLM_metal_labels_protocol_prediction_only"))
    add("proteins.tsv", tables["energetics"]["fields"] + extra_fields, protein_rows, "All archived native DFT outcomes, including39 unscored; no score, band or model substitution.", ["energetics", "proteins", "original_genes", "aliases", "motifs", "tree_names", "expression_genes", "genome_membership"])
    summary = {"proteins": len(protein_rows), "genes": len(gene_rows), "domains": len(selected_domains), "protein_domain_links": len(domain_rows), "genome_evidence_rows": len(genome_rows), "rna_profile_genes": len(expression), "rna_missing_genes": len(original)-len(expression), "rna_samples": len(sample_ids), "rna_cells": len(rows["expression_values"]), "within_genome_genes": len(genome_members), "within_genome_profiles": len(rows["genome_profiles"]), "within_genome_pairs": len(rows["genome_correlations"]), "score_status_counts": dict(Counter(r["status"] for r in protein_rows))}
    return {"tables": output, "summary": summary}


def write_export(data, built, output, elapsed_start):
    output = Path(output).resolve()
    require(not output.exists(), f"Output already exists; choose a new directory: {output}")
    output.mkdir(parents=True)
    output_files = {}
    for name, table in built["tables"].items():
        handle = io.StringIO(newline="")
        writer = csv.DictWriter(handle, fieldnames=table["fields"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table["rows"])
        raw = handle.getvalue().encode()
        if name.endswith(".gz"):
            raw = gzip.compress(raw, mtime=0)
        (output / name).write_bytes(raw)
        output_files[name] = {"sha256": digest(raw), "bytes": len(raw), "rows": len(table["rows"]), "columns": table["fields"], "description": table["description"], "sources": table["sources"]}
    metadata = {"schema": SCHEMA, "status": "PASS", "scope": "Existing evidence export only; no scientific values recomputed.", "inventory": data["inventory"], "input_pins": data["pins"], "outputs": output_files, "summary": built["summary"], "score_record": {k: v for k, v in data["score"].items() if k != "per_target"}, "transcript_methods": {"source_library_rna": "source-library samplewide TPM; original strings/availability retained", "rna_profiles": "assigned_QNAME_counts * 1e9 / (CDS_nt * mapper_input_RNA_read_ends); not TPM or genome-adjusted expression", "within_genome_profiles": "same read-depth-normalized expression; genome background enters pair-specific correlation, not this column", "within_genome_correlations": json.loads(data["documents"]["genome_methods"]), "historical_genome_adjusted_pair_profiles": json.loads(data["documents"]["historical_genome_methods"])}, "missing_policy": "Blank scientific value remains missing; measured zero is copied as zero. No fabricated rows for unavailable RNA profiles. All200 genes and176 proteins retained.", "candidate_results": "No candidate/research score overlay implemented; this delivery remains the archived native DFT scan."}
    (output / "EXPORT.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output / "README.md").write_text("# Nikasha: existing PLM evidence export\n\nStart with `proteins.tsv` (all176 outcomes) and `EXPORT.json` (input hashes, original DFT protocol/bands and table schemas). `genes.tsv` retains all200 source genes, including123 without exported25-sample RNA profiles. Exact gene-to-genome evidence is in `gene_genomes.tsv`; current domain/tree links are in `protein_domains.tsv`. Child tables keep their original values and measurement meanings.\n\n`rna_profiles.tsv` and `within_genome_profiles.tsv` are read-depth-normalized transcript recruitment, not TPM or genome-adjusted abundance. `source_library_rna.tsv` contains source-library samplewide TPM. `within_genome_correlations.tsv.gz` contains the existing pair-specific genome/campaign-adjusted correlations. The separately named historical pair-profile table preserves its narrower earlier method.\n\nScores are the existing native r2SCAN-3c/CPCM fixed-core DFT results, not current fast-mode or research-pool scores. All39 unscored proteins remain explicit. These unknown PLM labels are predictions, not affinity or physiological occupancy measurements. Independent mappings/genome aliases are not biological replicates. No new calculations, normalization, tree distances, correlations or biological interpretation were generated.\n")
    receipt = {"schema": "nikasha_plm_export_receipt_v1", "completed_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv, "python": sys.executable, "wall_seconds": time.monotonic()-elapsed_start, "scientific_calls": 0, "implementation": {"path": str(Path(__file__).resolve()), "sha256": digest(Path(__file__).read_bytes())}, "export_sha256": digest((output / "EXPORT.json").read_bytes())}
    (output / "EXECUTION.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return metadata


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", type=Path, required=True, help="Hash-pinned export input inventory; no implicit data paths")
    ap.add_argument("--output", type=Path, required=True, help="New output directory; existing paths are never overwritten")
    args = ap.parse_args(argv)
    started = time.monotonic()
    try:
        data = load_inputs(args.inventory)
        built = build_export(data)
        metadata = write_export(data, built, args.output, started)
    except (ExportError, OSError, KeyError, json.JSONDecodeError) as exc:
        ap.exit(2, f"Export unavailable: {exc}\n")
    print(json.dumps({"status": metadata["status"], "output": str(args.output.resolve()), "summary": metadata["summary"]}, indent=2))


if __name__ == "__main__":
    main()
