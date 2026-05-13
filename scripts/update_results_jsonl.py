#!/usr/bin/env python3
"""
Incremental aggregator for the Ca/Ln³⁺ DFT discriminator panel.

Maintains `results/all_results.jsonl` (one JSON object per line, append-only).
Per row we capture rich metadata: stem, source label, iPTM, atom count,
charge, energies, ΔΔE(Ca − La), classification, ORCA termination flags, etc.

Idempotent:
  - Skips dirs whose row in the JSONL is newer than the dir's mtime.
  - Use --rebuild to wipe the file and re-process every dir from scratch.

At the end, regenerates `results/discriminator_panel_LATEST.tsv` (sorted by
ΔΔE desc) for backward compatibility with `analyze_panel.py` consumers.

Usage:
  python scripts/update_results_jsonl.py <workspace>
  python scripts/update_results_jsonl.py <workspace> --rebuild
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants — must match analyze_panel.py
# ---------------------------------------------------------------------------
HA2KCAL = 627.5095
E_AQUO_CA = -1288.921041427522
E_AQUO_LA = -642.856485875777
DIFF_AQUO = E_AQUO_CA - E_AQUO_LA  # -646.064555551745 Ha

# Dirs we never touch (legacy P3 work, handoffs, panels with different layout)
SKIP_DIRS = {
    "qmmm",
    "qmmm_lc",
    "qmmm_md",
    "colin_handoff",
    "b97_3c_panel",
    "calexcitin_size_panel",
    "inbox",
    "inputs",
    "logs",
    "params",
    "protenix_jobs",
    "results",
    "scripts",
}


# ---------------------------------------------------------------------------
# Source detection
# ---------------------------------------------------------------------------
_RE_SPICYLAMS = re.compile(
    r"^(?:NC_|NZ_|CP\d|JA[A-Z]+\d|FSRD|MUNW|DYYG|NCCZ)"
    r".*_iptm[\d_.]+_s\d+$"
)
_RE_IPTM = re.compile(r"_iptm0[._](\d+)_s\d+")

_VAULT_PREFIXES = (
    "marco_", "alphap_", "hyphomyst_", "typeAorphan_", "lanmfusion_",
    "terbdiv_", "lanthhyd_", "cbb3oxid_", "marcobetar_", "rtxbroll_",
    "arginase_", "agmatinase_", "cytoc_", "midas_", "m24pep_",
    "adendeam_", "aminopepP_", "fernperox", "ferncellul", "fernparv_",
    "ferncalmod_", "gcamp3_", "calcineurinB_", "verruco_", "rifoxy_LaH",
    "novelfold_", "opaque_", "phase1_", "tannasepara_", "duf882_",
    "fern_endogluc", "tannase", "rifoxy_", "calexcitin", "lanm",
    "mxaf", "prolidase_",
)

# Curated vault-panel stems whose names don't fit a prefix rule.
_VAULT_EXACT = {
    "xoxf", "rifoxy", "calexcitin", "tannase", "lanm", "mxaf",
}

_RE_PDB_STYLE = re.compile(r"^1[A-Z0-9]{3}$")  # 1CLL, 1B9A, 1A75


def detect_source(stem: str) -> str:
    """Assign a source label for a per-candidate stem."""
    s = stem
    sl = s.lower()

    if sl.startswith("colinpqq_la_"):
        return "colin_la_verified"
    if sl.startswith("colinpqq_ca_"):
        return "colin_ca_verified"
    if sl.startswith("droideka"):
        return "colin_droideka"

    # Methylococcales inbox: anything mentioning the order or its genera
    if "_methylococcales" in sl:
        return "methylococcales_inbox"
    if any(g in sl for g in (
        "_methylotenera", "_methylocaldum_", "_methylobacillus_",
        "_methylomicrobium_",
    )):
        return "methylococcales_inbox"

    # Spicylams: NCBI-accession-like with _iptm_s# tail
    if _RE_SPICYLAMS.match(s):
        return "spicylams"
    # Generic fallback for anything ending in _iptm…_s# (Spicylams-like batch)
    if _RE_IPTM.search(s):
        return "spicylams"

    # Vault panel (curated controls + designed scaffolds)
    if s in _VAULT_EXACT or s.startswith(_VAULT_PREFIXES) \
            or _RE_PDB_STYLE.match(s):
        return "vault_panel"

    # Colin's PQQ controls (AFDB-style ID at start or end of stem)
    afdb = r"(?:A0A[A-Z0-9]+|Q\d[A-Z0-9]+|P\d[A-Z0-9]+|" \
           r"H8[A-Z0-9]+|C5[A-Z0-9]+|C6[A-Z0-9]+)"
    if re.search(rf"{afdb}$", s) or re.match(rf"^{afdb}_", s):
        return "colin_extras"

    return "unknown"


def extract_iptm(stem: str) -> float | None:
    m = _RE_IPTM.search(stem)
    if not m:
        return None
    digits = m.group(1)
    # "980" -> 0.980, "9520" -> 0.9520, etc.
    return float("0." + digits)


# ---------------------------------------------------------------------------
# ORCA output parsing
# ---------------------------------------------------------------------------
def grab_energy(out: Path) -> float | None:
    """Return last FINAL SINGLE POINT ENERGY from an ORCA .out file."""
    if not out.exists():
        return None
    try:
        for line in reversed(out.read_text(errors="ignore").splitlines()):
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])
    except Exception:
        return None
    return None


def is_terminated(out: Path) -> bool:
    if not out.exists():
        return False
    try:
        return "ORCA TERMINATED NORMALLY" in out.read_text(errors="ignore")
    except Exception:
        return False


def parse_charge_from_inp(inp: Path) -> int | None:
    """Parse '* xyzfile <CHARGE> <MULT> ...' line."""
    if not inp.exists():
        return None
    try:
        for line in inp.read_text(errors="ignore").splitlines():
            line = line.strip()
            if line.startswith("*") and "xyzfile" in line:
                parts = line.split()
                # Format: "* xyzfile -2 1 stem_La_qm.xyz"
                idx = parts.index("xyzfile")
                return int(parts[idx + 1])
    except Exception:
        return None
    return None


def n_atoms_from_xyz(xyz: Path) -> int | None:
    if not xyz.exists():
        return None
    try:
        first = xyz.read_text(errors="ignore").splitlines()[0].strip()
        return int(first)
    except Exception:
        return None


_RE_QM_REGION = re.compile(
    r"QM region:\s*(\d+)\s*atoms?\s*\((\d+)\s*carboxylates?",
    re.IGNORECASE,
)
_RE_FIRST_SHELL_LIST = re.compile(
    r"First-shell residues(?:\s*\(protein\))?:\s*\[([^\]]*)\]",
    re.IGNORECASE,
)


def parse_carve_log(qm_dir: Path) -> tuple[int | None, list[str]]:
    """Best-effort scan of any carve-log-like file in `qm_dir` for n_carbox
    and first_shell_residues. Returns (n_carbox, residues_list)."""
    n_carbox: int | None = None
    residues: list[str] = []

    # Candidates: any .log, carve_*.txt, or slurm logs (some pipelines pipe
    # carve stdout into the slurm.out before the QM banner).
    candidates: list[Path] = []
    candidates.extend(sorted(qm_dir.glob("*.log")))
    candidates.extend(sorted(qm_dir.glob("carve*.txt")))
    candidates.extend(sorted(qm_dir.glob("carve*.out")))
    candidates.extend(sorted(qm_dir.glob("slurm_*.out")))

    for f in candidates:
        try:
            text = f.read_text(errors="ignore")
        except Exception:
            continue

        if n_carbox is None:
            m = _RE_QM_REGION.search(text)
            if m:
                try:
                    n_carbox = int(m.group(2))
                except ValueError:
                    pass

        if not residues:
            m = _RE_FIRST_SHELL_LIST.search(text)
            if m:
                items = [x.strip().strip("'\"") for x in m.group(1).split(",")]
                residues = [x for x in items if x]

        if n_carbox is not None and residues:
            break

    return n_carbox, residues


# ---------------------------------------------------------------------------
# Classification (matches analyze_panel.py)
# ---------------------------------------------------------------------------
def classify(ddE_kcal: float | None, n_atoms: int | None) -> str:
    if ddE_kcal is None:
        return "pending"
    n = n_atoms if (n_atoms is not None) else 0
    if n < 15:
        return "EXCLUDE (under-carved)"
    if abs(ddE_kcal) > 30 and n < 50:
        return f"OUTLIER ({'Ca' if ddE_kcal < 0 else 'Ln'}, verify SCF)"
    if ddE_kcal >= 20:
        return "Ln-evolved"
    if ddE_kcal >= 5:
        return "Ln-preferring"
    if ddE_kcal >= 0:
        return "marginal"
    if ddE_kcal >= -5:
        return "ambiguous"
    return "Ca-evolved"


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------
def extract_row(qm_dir: Path) -> dict[str, Any]:
    base = qm_dir.name
    stem = base[:-3] if base.endswith("_qm") else base

    la_out = qm_dir / f"sp_{stem}_La.out"
    ca_out = qm_dir / f"sp_{stem}_Ca.out"
    apo_out = qm_dir / f"sp_{stem}_apo.out"
    water_out = qm_dir / f"sp_{stem}_water.out"
    if not water_out.exists():
        alt = qm_dir / "sp_bulk_water.out"
        if alt.exists():
            water_out = alt

    la_inp = qm_dir / f"sp_{stem}_La.inp"
    la_xyz = qm_dir / f"{stem}_La_qm.xyz"

    E_La = grab_energy(la_out)
    E_Ca = grab_energy(ca_out)
    E_apo = grab_energy(apo_out)
    E_water = grab_energy(water_out)

    term_la = is_terminated(la_out)
    term_ca = is_terminated(ca_out)
    term_apo = is_terminated(apo_out)
    term_water = is_terminated(water_out)

    n_atoms_la = n_atoms_from_xyz(la_xyz)
    charge_la = parse_charge_from_inp(la_inp)

    n_carbox, first_shell_residues = parse_carve_log(qm_dir)

    if E_La is not None and E_Ca is not None:
        ddE_kcal = ((E_Ca - E_La) - DIFF_AQUO) * HA2KCAL
        ddE_kcal = round(ddE_kcal, 4)
    else:
        ddE_kcal = None

    cls = classify(ddE_kcal, n_atoms_la)

    completed_at = _dt.datetime.fromtimestamp(
        qm_dir.stat().st_mtime,
    ).astimezone().isoformat(timespec="seconds")

    return {
        "stem": stem,
        "source": detect_source(stem),
        "iptm": extract_iptm(stem),
        "n_atoms_la": n_atoms_la,
        "n_carbox": n_carbox,
        "charge_la": charge_la,
        "first_shell_residues": first_shell_residues,
        "E_La_Ha": E_La,
        "E_Ca_Ha": E_Ca,
        "E_apo_Ha": E_apo,
        "E_water_Ha": E_water,
        "ddE_kcal": ddE_kcal,
        "class": cls,
        "term_la": term_la,
        "term_ca": term_ca,
        "term_apo": term_apo,
        "term_water": term_water,
        "completed_at": completed_at,
        "cluster_panel_path": str(qm_dir),
    }


# ---------------------------------------------------------------------------
# JSONL persistence
# ---------------------------------------------------------------------------
def load_existing(jsonl: Path) -> dict[str, dict[str, Any]]:
    """Return latest row keyed by stem (later rows override earlier)."""
    if not jsonl.exists():
        return {}
    out: dict[str, dict[str, Any]] = {}
    with jsonl.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            stem = row.get("stem")
            if stem:
                out[stem] = row
    return out


def row_is_fresh(row: dict[str, Any], qm_dir: Path) -> bool:
    """True if `row.completed_at` is at or after the dir's mtime."""
    ts = row.get("completed_at")
    if not ts:
        return False
    try:
        row_dt = _dt.datetime.fromisoformat(ts)
    except ValueError:
        return False
    if row_dt.tzinfo is None:
        row_dt = row_dt.astimezone()
    # Truncate dir mtime to whole seconds to match `timespec="seconds"` storage
    dir_mtime = int(qm_dir.stat().st_mtime)
    dir_dt = _dt.datetime.fromtimestamp(dir_mtime).astimezone()
    return row_dt >= dir_dt


def append_row(jsonl: Path, row: dict[str, Any]) -> None:
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    with jsonl.open("a") as fh:
        fh.write(json.dumps(row, separators=(",", ":")) + "\n")


# ---------------------------------------------------------------------------
# TSV regeneration (matches analyze_panel.py format)
# ---------------------------------------------------------------------------
def write_tsv(rows: list[dict[str, Any]], out_tsv: Path) -> None:
    rows_sorted = sorted(
        rows,
        key=lambda r: (r["ddE_kcal"] is None, -(r["ddE_kcal"] or 0)),
    )
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as f:
        f.write("# DFT Ca²⁺/Ln³⁺ discriminator — auto-generated\n")
        f.write(
            f"# E_aquo_Ca − E_aquo_La = {DIFF_AQUO:.6f} Ha "
            f"(= {DIFF_AQUO * HA2KCAL:.2f} kcal/mol)\n"
        )
        f.write(
            "# ΔΔE(Ca − La) = (E_holo_Ca − E_holo_La) − "
            "(E_aquo_Ca − E_aquo_La)\n"
        )
        f.write(
            "# Class: Ln-evolved ≥ +20; Ln-preferring +5 to +20; "
            "marginal 0 to +5; ambiguous -5 to 0; Ca-evolved < -5\n"
        )
        f.write("stem\tn_atoms\tE_La\tE_Ca\tE_apo\tE_water\tddE_kcal\tclass\n")
        for r in rows_sorted:
            n_atoms = (
                str(r["n_atoms_la"]) if r.get("n_atoms_la") is not None
                else "?"
            )
            ddE_s = (
                f"{r['ddE_kcal']:+.2f}" if r["ddE_kcal"] is not None
                else "pending"
            )
            E_La = (
                f"{r['E_La_Ha']:.6f}" if r.get("E_La_Ha") is not None
                else "-"
            )
            E_Ca = (
                f"{r['E_Ca_Ha']:.6f}" if r.get("E_Ca_Ha") is not None
                else "-"
            )
            E_apo = (
                f"{r['E_apo_Ha']:.6f}" if r.get("E_apo_Ha") is not None
                else "-"
            )
            E_w = (
                f"{r['E_water_Ha']:.6f}" if r.get("E_water_Ha") is not None
                else "-"
            )
            f.write(
                f"{r['stem']}\t{n_atoms}\t{E_La}\t{E_Ca}\t{E_apo}\t{E_w}"
                f"\t{ddE_s}\t{r['class']}\n"
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "workspace",
        help="Workspace root (contains the *_qm/ candidate dirs).",
    )
    p.add_argument(
        "--rebuild",
        action="store_true",
        help="Delete the existing JSONL and re-process every dir.",
    )
    args = p.parse_args(argv)

    workspace = Path(args.workspace).resolve()
    if not workspace.is_dir():
        print(f"ERROR: workspace not found: {workspace}", file=sys.stderr)
        return 2

    results_dir = workspace / "results"
    jsonl = results_dir / "all_results.jsonl"
    tsv = results_dir / "discriminator_panel_LATEST.tsv"

    if args.rebuild and jsonl.exists():
        jsonl.unlink()

    existing = load_existing(jsonl)

    # Walk candidate dirs
    n_new = 0
    n_updated = 0
    n_skipped = 0
    for qm_dir in sorted(workspace.glob("*_qm")):
        if not qm_dir.is_dir():
            continue
        if qm_dir.name in SKIP_DIRS:
            continue
        # Belt-and-suspenders: skip qmmm-prefixed
        if qm_dir.name.startswith("qmmm"):
            continue

        stem = qm_dir.name[:-3]
        prev = existing.get(stem)
        if prev and row_is_fresh(prev, qm_dir):
            n_skipped += 1
            continue

        try:
            row = extract_row(qm_dir)
        except Exception as exc:
            print(f"  [warn] {qm_dir.name}: {exc}", file=sys.stderr)
            continue

        append_row(jsonl, row)
        existing[stem] = row
        if prev is None:
            n_new += 1
        else:
            n_updated += 1

    # Write TSV from current latest-by-stem rows
    rows = list(existing.values())
    write_tsv(rows, tsv)

    # Summary
    by_class: dict[str, int] = {}
    for r in rows:
        by_class[r["class"]] = by_class.get(r["class"], 0) + 1

    print(f"Workspace:  {workspace}")
    print(f"JSONL:      {jsonl}")
    print(f"TSV:        {tsv}")
    print(f"New rows added:     {n_new}")
    print(f"Existing rows refreshed: {n_updated}")
    print(f"Up-to-date (skipped):    {n_skipped}")
    print(f"Total rows in JSONL:     {len(rows)}")
    print()
    print("Class breakdown:")
    for cls in sorted(by_class):
        print(f"  {cls:35s} {by_class[cls]}")

    # Auto-rebuild enriched Ln-class tracker (results/ln_class_hits.tsv)
    try:
        import subprocess, sys
        rebuild = Path(__file__).parent / "rebuild_ln_tracker.py"
        if rebuild.exists():
            subprocess.run(["python3", str(rebuild)], check=False)
    except Exception as e:
        import sys as _sys
        print(f"warn: ln_class_hits.tsv auto-rebuild failed: {e}", file=_sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
