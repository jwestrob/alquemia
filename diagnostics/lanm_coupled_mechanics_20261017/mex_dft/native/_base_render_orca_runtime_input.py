#!/usr/bin/env python3
"""Render an immutable ORCA template with allocation-aware MPI provenance.

The scientific input written by a carver remains hashed and unchanged.  On an
allocated node this helper inserts only a ``%pal nprocs`` block, writes the
concrete runtime input, and records both hashes plus the SLURM allocation.  The
result parser verifies this sidecar before accepting an output generated from a
runtime input.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any


RUNTIME_SCHEMA = "alchemical_bvs.orca_runtime_input.v1"
PAL_BLOCK_RE = re.compile(r"^\s*%pal\b", re.IGNORECASE | re.MULTILINE)
PAL_SIMPLE_RE = re.compile(r"(?:^|\s)pal(?:\d+|\()", re.IGNORECASE)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def render_runtime_input(
    template_path: Path,
    runtime_path: Path,
    *,
    nprocs: int,
    sidecar_path: Path | None = None,
) -> dict[str, Any]:
    if nprocs < 1:
        raise ValueError(f"nprocs must be positive, got {nprocs}")
    allocated = os.environ.get("SLURM_CPUS_ON_NODE")
    if allocated and allocated.isdigit() and nprocs > int(allocated):
        raise ValueError(
            f"nprocs={nprocs} exceeds SLURM_CPUS_ON_NODE={allocated}"
        )

    template_path = template_path.resolve()
    runtime_path = runtime_path.resolve()
    text = template_path.read_text()
    if PAL_BLOCK_RE.search(text) or any(
        PAL_SIMPLE_RE.search(line)
        for line in text.splitlines()
        if line.lstrip().startswith("!")
    ):
        raise ValueError(
            f"template already declares PAL parallelism; refusing a second policy: {template_path}"
        )
    lines = text.splitlines(keepends=True)
    coordinate_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.lstrip().lower().startswith("* xyzfile")
        ),
        None,
    )
    if coordinate_index is None:
        raise ValueError(f"template lacks a '* xyzfile' line: {template_path}")
    pal_block = f"%pal\n  nprocs {nprocs}\nend\n\n"
    rendered = "".join(lines[:coordinate_index]) + pal_block + "".join(
        lines[coordinate_index:]
    )
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = runtime_path.with_name(runtime_path.name + ".tmp")
    temporary.write_text(rendered)
    temporary.replace(runtime_path)

    record: dict[str, Any] = {
        "schema_version": RUNTIME_SCHEMA,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "template": {
            "path": str(template_path),
            "sha256": sha256_file(template_path),
        },
        "runtime_input": {
            "path": str(runtime_path),
            "sha256": sha256_file(runtime_path),
        },
        "parallelism": {
            "nprocs": nprocs,
            "policy": "allocation_aware_mpi_rank_partition",
            "slurm_cpus_on_node": int(allocated) if allocated and allocated.isdigit() else None,
        },
        "execution_context": {
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_job_partition": os.environ.get("SLURM_JOB_PARTITION"),
            "host": os.uname().nodename,
        },
    }
    if sidecar_path is None:
        sidecar_path = runtime_path.with_suffix(".json")
    sidecar_path = sidecar_path.resolve()
    sidecar_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path)
    parser.add_argument("runtime", type=Path)
    parser.add_argument("--nprocs", type=int, required=True)
    parser.add_argument("--sidecar", type=Path)
    args = parser.parse_args()
    record = render_runtime_input(
        args.template,
        args.runtime,
        nprocs=args.nprocs,
        sidecar_path=args.sidecar,
    )
    print(
        f"rendered {record['runtime_input']['path']} with "
        f"nprocs={record['parallelism']['nprocs']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
