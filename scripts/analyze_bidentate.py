#!/usr/bin/env python3
"""
analyze_bidentate.py — for each cluster, count bidentate Asp/Glu κ²-O,O
coordinations to the metal.

Detection: two O atoms are first-shell (≤3.0 Å of metal) AND share a common
carboxylate C neighbor (each O within 1.5 Å of the same C). That's the
unambiguous Asp/Glu carboxylate signature; Asn/Gln amides have only one O
per C, so they don't appear here.

Output: extends results/all_results.jsonl with `n_bidentate_asp_glu` field.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

ALCH = Path(__file__).resolve().parent.parent
JSONL = ALCH / "results" / "all_results.jsonl"

FIRST_SHELL = 3.0     # Å; metal-O cutoff
CO_BOND = 1.5         # Å; C-O covalent cutoff for carboxylate detection


def parse_xyz(path: Path):
    """Return list of (symbol, x, y, z)."""
    if not path.exists():
        return None
    try:
        with open(path) as f:
            lines = f.readlines()
        if len(lines) < 3:
            return None
        atoms = []
        for ln in lines[2:]:
            parts = ln.split()
            if len(parts) < 4:
                continue
            try:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            except ValueError:
                continue
        return atoms
    except Exception:
        return None


def dist(a, b):
    return math.sqrt((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)


def count_bidentate(atoms) -> int | None:
    """Count Asp/Glu carboxylate bidentate κ²-O,O coordinations."""
    if not atoms:
        return None
    # Metal is first atom
    metal = atoms[0]
    if metal[0] not in ("La", "Ca"):
        return None

    # Index O atoms in first shell and all C atoms
    first_shell_O_idx = []
    C_idx = []
    for i, a in enumerate(atoms[1:], start=1):
        d = dist(a, metal)
        if a[0] == "O" and d <= FIRST_SHELL:
            first_shell_O_idx.append(i)
        if a[0] == "C":
            C_idx.append(i)

    if len(first_shell_O_idx) < 2:
        return 0

    # For each O in first shell, find its bonded C (≤1.5 Å)
    O_to_C = {}
    for oi in first_shell_O_idx:
        oa = atoms[oi]
        for ci in C_idx:
            ca = atoms[ci]
            if dist(oa, ca) <= CO_BOND:
                O_to_C[oi] = ci
                break

    # Count pairs of first-shell O sharing the same C
    n_bidentate = 0
    seen_C = set()
    o_by_c = {}
    for oi, ci in O_to_C.items():
        o_by_c.setdefault(ci, []).append(oi)
    for ci, ois in o_by_c.items():
        if len(ois) >= 2:
            n_bidentate += 1
    return n_bidentate


def main():
    entries = []
    with open(JSONL) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            try:
                entries.append(json.loads(ln))
            except:
                continue

    n_total = len(entries)
    n_updated = 0
    n_missing = 0
    for e in entries:
        qmdir = e.get("cluster_panel_path")
        stem = e.get("stem")
        if not qmdir or not stem:
            continue
        xyz = Path(qmdir) / f"{stem}_La_qm.xyz"
        atoms = parse_xyz(xyz)
        nb = count_bidentate(atoms) if atoms else None
        if nb is None:
            n_missing += 1
            e["n_bidentate_asp_glu"] = None
        else:
            e["n_bidentate_asp_glu"] = nb
            n_updated += 1

    with open(JSONL, "w") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")

    print(f"Manifest updated: {JSONL}")
    print(f"  Total:              {n_total}")
    print(f"  n_bidentate set:    {n_updated}")
    print(f"  Missing XYZ:        {n_missing}")


if __name__ == "__main__":
    main()
