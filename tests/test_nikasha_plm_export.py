"""Real pinned PLM exports; malformed cases are corrupted copies of those records."""
import importlib.util
import json
import time
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("nikasha_plm_export", ROOT / "scripts/nikasha_plm_export.py")
export = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(export)
INPUTS = ROOT / "diagnostics/nikasha_plm_export_20260922/INPUTS.json"


@pytest.fixture(scope="module")
def actual():
    data = export.load_inputs(INPUTS)
    return data, export.build_export(data)


def changed_table(data, name, change):
    """Corrupt only a copied real table; retain all scientific source files."""
    result = dict(data, tables=dict(data["tables"]))
    table = data["tables"][name]
    rows = [dict(row) for row in table["rows"]]
    change(rows)
    result["tables"][name] = dict(table, rows=rows)
    return result


def test_complete_cohort_and_actual_method(actual):
    data, built = actual
    summary = built["summary"]
    assert (summary["proteins"], summary["genes"], summary["domains"]) == (176, 200, 170)
    assert summary["score_status_counts"] == {"scored": 137, "unsupported": 37, "unscored_source_residue": 2}
    source = {r["target_id"]: r for r in data["tables"]["energetics"]["rows"]}
    for row in built["tables"]["proteins.tsv"]["rows"]:
        assert row["protocol_id"] == export.DFT_PROTOCOL
        assert row["score_method"] == "native_ORCA_r2SCAN-3c_CPCM_water_fixed_core_DFT"
        assert {key: row[key] for key in source[row["target_id"]]} == source[row["target_id"]]
        if row["status"] != "scored":
            assert row["S_kcal_mol"] == row["R_hartree"] == ""


def test_all_existing_transcript_strings_preserved(actual):
    data, built = actual
    pairs = {
        "rna_profiles.tsv": "expression_values",
        "within_genome_profiles.tsv": "genome_profiles",
        "within_genome_correlations.tsv.gz": "genome_correlations",
        "historical_genome_adjusted_pair_profiles.tsv": "old_genome_profiles",
    }
    for output_name, source_name in pairs.items():
        rows = built["tables"][output_name]["rows"]
        assert [{k: v for k, v in r.items() if k != "target_id"} for r in rows] == data["tables"][source_name]["rows"]
    assert len(built["tables"]["rna_profiles.tsv"]["rows"]) == 1925
    assert len(built["tables"]["within_genome_profiles.tsv"]["rows"]) == 700
    assert len(built["tables"]["within_genome_correlations.tsv.gz"]["rows"]) == 40426


def test_missing_profiles_remain_missing_and_source_TPM_is_separate(actual):
    data, built = actual
    genes = built["tables"]["genes.tsv"]["rows"]
    missing = [r for r in genes if r["rna_profile_status"] == "missing_from_existing_export"]
    assert len(missing) == 123
    assert all(r["rna_profile_sample_count"] == r["rna_profile_reference"] == "" for r in missing)
    represented = {r["gene_id"] for r in built["tables"]["rna_profiles.tsv"]["rows"]}
    assert not represented & {r["gene_id"] for r in missing}
    aliases = {}
    for row in data["tables"]["aliases"]["rows"]:
        aliases.setdefault(row["gene_id"], row)
    rna = built["tables"]["source_library_rna.tsv"]["rows"]
    for row in rna:
        for field in ["source_RNA_status", "source_RNA_counts", "source_RNA_TPM"]:
            assert row[field] == aliases[row["gene_id"]][field]
    assert sum(r["source_RNA_status"] == "RNA_unavailable_for_source_library" for r in rna) == 124
    assert all(r["source_RNA_TPM"] == r["source_RNA_counts"] == "" for r in rna if r["source_RNA_status"] == "RNA_unavailable_for_source_library")


def test_genome_evidence_and_current_tree_ids_are_original(actual):
    data, built = actual
    source_rows = {name: [json.dumps(r, sort_keys=True, separators=(",", ":")) for r in data["tables"][name]["rows"]] for name in ["aliases", "expression_genes", "genome_membership"]}
    for row in built["tables"]["gene_genomes.tsv"]["rows"]:
        assert row["source_row_json"] in source_rows[row["evidence_source"]]
    tree = {r["original_leaf_id"]: r["new_leaf_id"] for r in data["tables"]["tree_names"]["rows"]}
    domains = built["tables"]["protein_domains.tsv"]["rows"]
    assert len({r["domain_id"] for r in domains}) == 170
    assert all(r["current_tree_leaf_id"] == tree[r["domain_id"]] for r in domains)


def test_corrupted_actual_input_bytes_fail_hash(actual, tmp_path):
    data, _ = actual
    pin = data["pins"]["energetics"]
    corrupted = tmp_path / "explicitly_corrupted_real_results.tsv"
    corrupted.write_bytes(Path(pin["path"]).read_bytes() + b"corrupted fixture\n")
    with pytest.raises(export.ExportError, match="SHA256 mismatch"):
        export.read_pin(dict(pin, path=str(corrupted)))


def test_corrupted_actual_gene_hash_fails(actual):
    data, _ = actual
    def corrupt(rows):
        old = rows[0]["source_protein_sha256"]
        rows[0]["source_protein_sha256"] = ("0" if old[0] != "0" else "1") + old[1:]
    with pytest.raises(export.ExportError, match="expression protein hash mismatch"):
        export.build_export(changed_table(data, "expression_genes", corrupt))


def test_corrupted_actual_missing_cell_is_not_silently_zero_filled(actual):
    data, _ = actual
    with pytest.raises(export.ExportError, match="sample coverage mismatch"):
        export.build_export(changed_table(data, "expression_values", lambda rows: rows.pop()))


def test_corrupted_actual_score_and_duplicate_cohort_fail(actual):
    data, _ = actual
    def corrupt(rows):
        rows[0]["S_kcal_mol"] = ""
    with pytest.raises(export.ExportError, match="TSV/JSON result mismatch"):
        export.build_export(changed_table(data, "energetics", corrupt))
    with pytest.raises(export.ExportError, match="Duplicate or empty protein"):
        export.build_export(changed_table(data, "proteins", lambda rows: rows.append(dict(rows[0]))))


def test_existing_output_is_never_overwritten(actual, tmp_path):
    data, built = actual
    marker = tmp_path / "preserved.txt"
    marker.write_text("existing output")
    with pytest.raises(export.ExportError, match="Output already exists"):
        export.write_export(data, built, tmp_path, time.monotonic())
    assert marker.read_text() == "existing output"
    assert list(tmp_path.iterdir()) == [marker]
