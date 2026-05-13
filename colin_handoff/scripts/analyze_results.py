#!/usr/bin/env python3
"""
Compute ΔΔE(Ca − La) for completed candidates in a workspace dir.
Walks through `<workspace>/*_qm/sp_*_La.out` and `sp_*_Ca.out` pairs,
extracts FINAL SINGLE POINT ENERGY values, applies the discriminator cycle.

Usage: python analyze_results.py <workspace_dir>
Outputs: <workspace_dir>/discriminator_results.tsv
"""
from __future__ import annotations
import sys
import re
from pathlib import Path

HA2KCAL = 627.5095

# Universal aquo references (idealized [M(H2O)8]^q+, r²SCAN-3c CPCM(Water) DefGrid3)
# DO NOT change these — must stay constant across all panel members for ΔΔE comparability.
E_AQUO_CA = -1288.921041427522
E_AQUO_LA =  -642.856485875777
DIFF_AQUO = E_AQUO_CA - E_AQUO_LA  # -646.064556 Ha


def grab_final_sp(out_path: Path) -> float | None:
    if not out_path.exists():
        return None
    try:
        for line in reversed(out_path.read_text(errors="ignore").splitlines()):
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])
    except Exception:
        return None
    return None


def classify(ddE_kcal: float, n_atoms: int) -> str:
    if n_atoms < 15:
        return "EXCLUDE (under-carved)"
    if abs(ddE_kcal) > 30 and n_atoms < 50:
        return f"OUTLIER ({'Ca' if ddE_kcal < 0 else 'Ln'}, verify SCF/geometry)"
    if ddE_kcal >= 20:
        return "Ln-evolved"
    elif ddE_kcal >= 5:
        return "Ln-preferring"
    elif ddE_kcal >= 0:
        return "marginal"
    elif ddE_kcal >= -5:
        return "ambiguous"
    else:
        return "Ca-evolved"


def main(workspace: Path):
    rows = []
    for qm_dir in sorted(workspace.glob("*_qm")):
        if not qm_dir.is_dir():
            continue
        la_outs = sorted(qm_dir.glob("sp_*_La.out"))
        ca_outs = sorted(qm_dir.glob("sp_*_Ca.out"))
        if not la_outs or not ca_outs:
            continue
        # Recover stem from filename
        m = re.match(r"sp_(.+)_La\.out", la_outs[0].name)
        stem = m.group(1) if m else qm_dir.name.removesuffix("_qm")

        E_La = grab_final_sp(la_outs[0])
        E_Ca = grab_final_sp(ca_outs[0])

        # Get atom count from xyz
        xyz = qm_dir / f"{stem}_La_qm.xyz"
        n_atoms = 0
        if xyz.exists():
            try:
                n_atoms = int(xyz.read_text().splitlines()[0].strip())
            except Exception:
                pass

        if E_La is not None and E_Ca is not None:
            ddE = (E_Ca - E_La) - DIFF_AQUO
            ddE_kcal = ddE * HA2KCAL
            cls = classify(ddE_kcal, n_atoms)
        else:
            ddE_kcal = None
            cls = "pending"

        rows.append({
            "stem": stem,
            "atoms": n_atoms,
            "E_La": E_La,
            "E_Ca": E_Ca,
            "ddE_kcal": ddE_kcal,
            "class": cls,
        })

    rows.sort(key=lambda r: (r["ddE_kcal"] is None, -(r["ddE_kcal"] or 0)))

    out_tsv = workspace / "discriminator_results.tsv"
    with open(out_tsv, "w") as f:
        f.write("# DFT Ca²⁺/Ln³⁺ discriminator results — auto-generated\n")
        f.write(f"# E_aquo_Ca − E_aquo_La = {DIFF_AQUO:.6f} Ha (universal reference)\n")
        f.write("# ΔΔE(Ca − La) = (E_holo_Ca − E_holo_La) − (E_aquo_Ca − E_aquo_La)\n")
        f.write("# Class thresholds: ≥+20 Ln-evolved; +5 to +20 Ln-preferring; 0 to +5 marginal;\n")
        f.write("#                   −5 to 0 ambiguous; <−5 Ca-evolved; |ΔΔE|>30 with <50 atoms = outlier\n")
        f.write("stem\tn_atoms\tE_La_Ha\tE_Ca_Ha\tddE_kcal\tclass\n")
        for r in rows:
            ddE_s = f"{r['ddE_kcal']:+.2f}" if r['ddE_kcal'] is not None else "pending"
            E_La = f"{r['E_La']:.6f}" if r['E_La'] is not None else "-"
            E_Ca = f"{r['E_Ca']:.6f}" if r['E_Ca'] is not None else "-"
            f.write(f"{r['stem']}\t{r['atoms']}\t{E_La}\t{E_Ca}\t{ddE_s}\t{r['class']}\n")

    # Print summary
    print(f"\n{'stem':50s} {'atoms':>6s} {'ddE_kcal':>9s}  class")
    print("-" * 100)
    for r in rows:
        ddE_s = f"{r['ddE_kcal']:+.2f}" if r["ddE_kcal"] is not None else "pending"
        print(f"{r['stem']:50s} {r['atoms']:>6d} {ddE_s:>9s}  {r['class']}")

    n_done = sum(1 for r in rows if r["ddE_kcal"] is not None)
    print(f"\nTotal: {len(rows)} candidates, {n_done} done, {len(rows)-n_done} pending")
    print(f"Wrote: {out_tsv}")


if __name__ == "__main__":
    workspace = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    main(workspace)
