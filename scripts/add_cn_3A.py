#!/usr/bin/env python3
"""
add_cn_3A.py — compute first-shell donor coordination number @3.0 Å (N+O only)
for every entry in results/all_results.jsonl and write a `cn_3A` field back to
the manifest in-place.

CN @3.0 Å counts heavy-atom donors (N or O) within 3.0 Å of the metal in the
La cluster XYZ. This is the standard Ln-O first-shell cutoff for r²SCAN-3c
clusters; it's generous on the Ln side but tight on the Ca side.

Usage:
    python3 scripts/add_cn_3A.py
"""
from __future__ import annotations
import json
import math
from pathlib import Path

ALCH = Path(__file__).resolve().parent.parent
JSONL = ALCH / "results" / "all_results.jsonl"
CUTOFF = 3.0  # Å


def count_donors_3A(xyz_path: Path) -> int | None:
    """Count N + O atoms within CUTOFF of the metal (first atom)."""
    if not xyz_path.exists():
        return None
    try:
        with open(xyz_path) as f:
            lines = f.readlines()
        if len(lines) < 3:
            return None
        # Line 0: n_atoms; line 1: comment; line 2+: atoms
        atoms = []
        for ln in lines[2:]:
            parts = ln.split()
            if len(parts) < 4:
                continue
            try:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            except ValueError:
                continue
        if not atoms:
            return None
        # Metal is first atom (La or Ca)
        metal = atoms[0]
        if metal[0] not in ("La", "Ca"):
            # not a metal cluster
            return None
        mx, my, mz = metal[1], metal[2], metal[3]
        cn = 0
        for sym, x, y, z in atoms[1:]:
            if sym not in ("N", "O"):
                continue
            d = math.sqrt((x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2)
            if d <= CUTOFF:
                cn += 1
        return cn
    except Exception:
        return None


def main():
    entries = []
    with open(JSONL) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            try:
                entries.append(json.loads(ln))
            except json.JSONDecodeError:
                continue

    n_total = len(entries)
    n_updated = 0
    n_no_xyz = 0
    for e in entries:
        qmdir = e.get("cluster_panel_path")
        if not qmdir:
            continue
        stem = e.get("stem")
        if not stem:
            continue
        # PQQ-refold stems may have different name; the La XYZ is at qmdir/stem_La_qm.xyz
        xyz = Path(qmdir) / f"{stem}_La_qm.xyz"
        cn = count_donors_3A(xyz)
        if cn is None:
            n_no_xyz += 1
            e["cn_3A"] = None
        else:
            e["cn_3A"] = cn
            n_updated += 1

    # Write back
    with open(JSONL, "w") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")

    print(f"Manifest updated: {JSONL}")
    print(f"  Total entries:      {n_total}")
    print(f"  CN @3Å assigned:    {n_updated}")
    print(f"  No XYZ / unparseable: {n_no_xyz}")


if __name__ == "__main__":
    main()
