#!/usr/bin/env python3
"""Snapshot all experimental-coordinate sources named by the master ledger."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import urllib.request


ROOT = Path(__file__).resolve().parent
PDB_CACHE = Path("/groups/banfield/projects/environmental/sr/srvp2020/protenix/mmcif")
PDB_IDS = (
    "1KB0", "1GLG", "4CPV", "1SL8", "1KSM", "2BCB", "4DZT", "1HCH",
    "1W4W", "3CNA", "6MI5", "6IP9", "1F6S", "9B1U", "9B1V", "1C9U",
    "1CQ1", "6JWF", "9C8X", "9VY7",
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def main() -> None:
    manifest_path = ROOT / "structure_sources.json"
    if manifest_path.exists():
        raise RuntimeError("structure snapshot already exists; refusing overwrite")
    output = ROOT / "sources" / "pdb"
    records = []
    for pdb_id in PDB_IDS:
        cached = PDB_CACHE / f"{pdb_id.lower()}.cif"
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        if cached.is_file():
            data = cached.read_bytes()
            origin = {"kind": "local_shared_cache", "path": str(cached.resolve())}
        else:
            request = urllib.request.Request(url, headers={"User-Agent": "alchemical-bvs/1"})
            with urllib.request.urlopen(request, timeout=180) as response:
                data = response.read()
            origin = {"kind": "rcsb_download", "url": url}
        first_line = data.splitlines()[0].decode(errors="replace").lower()
        if first_line != f"data_{pdb_id.lower()}":
            raise ValueError(f"unexpected mmCIF payload for {pdb_id}")
        destination = output / f"{pdb_id}.cif"
        atomic_write(destination, data)
        records.append(
            {
                "pdb_id": pdb_id,
                "path": str(destination.resolve()),
                "sha256": digest(data),
                "bytes": len(data),
                "origin": origin,
            }
        )
    payload = {
        "schema_version": "alchemical_bvs.external_structure_snapshot.v1",
        "record_count": len(records),
        "records": records,
    }
    atomic_write(manifest_path, (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode())
    print(f"snapshotted {len(records)} structures under {output}")


if __name__ == "__main__":
    main()
