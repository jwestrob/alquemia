#!/usr/bin/env python3
"""Materialize the frozen PQQ external-fold sequence panel with provenance."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import urllib.parse
import urllib.request


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[2]
EXPANSION_FASTA = PROJECT / "testset_expansion" / "expansion_sequences_v1.faa"
OUT = ROOT / "inputs"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

TARGETS = (
    ("H19_FAM1_XoxF5", "WP_024300827.1", "strict_Ln", "local_expansion"),
    ("H19_FAM1_ExaF", "WP_019918664.1", "strict_Ln", "local_expansion"),
    ("H19_Rkho_XoxF5", "WP_108028715.1", "strict_Ln", "local_expansion"),
    ("H19_Tcon_XoxF5", "WP_085126452.1", "strict_Ln", "local_expansion"),
    ("H19_Gmar_XoxF5", "WP_002539484.1", "strict_Ln", "local_expansion"),
    ("M107_PedH", "WP_245258612.1", "promiscuous_boundary", "local_expansion"),
    ("5GB1C_MxaF", "QCW83947.1", "Ca_functional", "ncbi"),
    ("LW13_XoxF", "QBC28924.1", "Ln_functional_mapping_M2", "ncbi"),
    ("LW13_MxaF", "QBC27574.1", "Ca_functional_mapping_M2", "ncbi"),
    ("DAMO_XoxF1", "CBE67239.1", "unresolved_boundary", "ncbi"),
    ("Msil_XoxF1", "WP_012592127.1", "unresolved_boundary", "ncbi"),
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def atomic_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_bytes(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def parse_fasta(data: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in data.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError("sequence precedes FASTA header")
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    for header, sequence in records:
        if not re.fullmatch(r"[A-Za-z*]+", sequence):
            raise ValueError(f"invalid sequence for {header}")
    return [(header, sequence.upper().rstrip("*")) for header, sequence in records]


def local_by_accession() -> dict[str, tuple[str, str]]:
    records = parse_fasta(EXPANSION_FASTA.read_text())
    output: dict[str, tuple[str, str]] = {}
    for header, sequence in records:
        accession = header.split("|", 1)[0]
        if accession in output:
            raise ValueError(f"duplicate local accession {accession}")
        output[accession] = (header, sequence)
    return output


def fetch_ncbi(accession: str) -> tuple[str, str, str, bytes]:
    query = urllib.parse.urlencode(
        {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
    )
    url = f"{NCBI_EFETCH}?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "alchemical-bvs/1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        raw = response.read()
    records = parse_fasta(raw.decode())
    if len(records) != 1:
        raise ValueError(f"NCBI returned {len(records)} records for {accession}")
    header, sequence = records[0]
    returned = header.split()[0].split("|", 1)[0]
    if returned.upper() != accession.upper():
        raise ValueError(f"NCBI returned {returned}, expected {accession}")
    return header, sequence, url, raw


def wrap_record(header: str, sequence: str) -> str:
    lines = [f">{header}"]
    lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    return "\n".join(lines) + "\n"


def main() -> None:
    if (ROOT / "sequence_provenance.json").exists():
        raise RuntimeError("sequence panel is already frozen; refusing overwrite")
    local = local_by_accession()
    source_hash = sha256_file(EXPANSION_FASTA)
    records: list[dict[str, object]] = []
    combined: list[str] = []
    for target_id, accession, label, source in TARGETS:
        if source == "local_expansion":
            if accession not in local:
                raise KeyError(f"{accession} absent from {EXPANSION_FASTA}")
            original_header, sequence = local[accession]
            source_record = {
                "kind": source,
                "path": str(EXPANSION_FASTA.resolve()),
                "sha256": source_hash,
                "original_header": original_header,
            }
        else:
            original_header, sequence, url, raw = fetch_ncbi(accession)
            raw_path = OUT / "raw_ncbi" / f"{accession}.fasta"
            atomic_bytes(raw_path, raw)
            source_record = {
                "kind": source,
                "url": url,
                "retrieved_record": str(raw_path.resolve()),
                "retrieved_record_sha256": sha256_bytes(raw),
                "original_header": original_header,
            }
        if len(sequence) < 300:
            raise ValueError(f"unexpectedly short PQQ enzyme {accession}: {len(sequence)} aa")
        canonical_header = f"{target_id}|{accession}|{label}"
        text = wrap_record(canonical_header, sequence)
        path = OUT / "per_target" / f"{target_id}.faa"
        atomic_bytes(path, text.encode())
        combined.append(text)
        records.append(
            {
                "target_id": target_id,
                "accession": accession,
                "biological_label": label,
                "length": len(sequence),
                "sequence_sha256": sha256_bytes(sequence.encode()),
                "canonical_fasta": str(path.resolve()),
                "canonical_fasta_sha256": sha256_bytes(text.encode()),
                "source": source_record,
            }
        )
    combined_path = OUT / "pqq_external_panel.faa"
    atomic_bytes(combined_path, "".join(combined).encode())
    atomic_json(
        ROOT / "sequence_provenance.json",
        {
            "schema_version": "alchemical_bvs.pqq_external_sequences.v1",
            "target_count": len(records),
            "combined_fasta": {
                "path": str(combined_path.resolve()),
                "sha256": sha256_file(combined_path),
            },
            "records": records,
        },
    )
    print(f"frozen {len(records)} sequences in {combined_path}")


if __name__ == "__main__":
    main()
