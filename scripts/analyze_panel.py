#!/usr/bin/env python3
"""
Scan all *_qm/ directories under alchemical_bvs/ for completed La/Ca SP outputs
and compute ΔΔE(Ca − La) using idealized [M(H2O)8]^q+ aquo references.

Outputs:
  results/discriminator_panel_LATEST.tsv
"""
from __future__ import annotations
from pathlib import Path
import re
import sys

ROOT = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs")
HA2KCAL = 627.5095

# r²SCAN-3c CPCM(Water) DefGrid3 idealized [M(H2O)8] references
E_AQUO_CA = -1288.921041427522
E_AQUO_LA =  -642.856485875777
DIFF_AQUO = E_AQUO_CA - E_AQUO_LA  # -646.064556 Ha


def grab(out: Path) -> float | None:
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
        text = out.read_text(errors="ignore")
        return "ORCA TERMINATED NORMALLY" in text
    except Exception:
        return False


def classify(ddE_kcal: float, n_atoms_str: str) -> str:
    # Sanity-flag outliers — if magnitude > 30 kcal/mol, double-check carve
    n_atoms = int(n_atoms_str) if n_atoms_str.isdigit() else 0
    if n_atoms < 15:
        return "EXCLUDE (under-carved)"
    if abs(ddE_kcal) > 30 and n_atoms < 50:
        return f"OUTLIER ({'Ca' if ddE_kcal < 0 else 'Ln'}, verify SCF)"
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


def main():
    rows = []
    for qm_dir in sorted(ROOT.glob("*_qm")):
        if qm_dir.is_dir() and qm_dir.name not in ("qmmm",):
            # Find La/Ca/apo/water by glob (case-insensitive stem)
            la_outs = sorted(qm_dir.glob("sp_*_La.out"))
            ca_outs = sorted(qm_dir.glob("sp_*_Ca.out"))
            apo_outs = sorted(qm_dir.glob("sp_*_apo.out"))
            water_outs = sorted(qm_dir.glob("sp_*_water.out")) + sorted(qm_dir.glob("sp_bulk_water.out"))
            la_out = la_outs[0] if la_outs else qm_dir / "MISSING_La.out"
            ca_out = ca_outs[0] if ca_outs else qm_dir / "MISSING_Ca.out"
            apo_out = apo_outs[0] if apo_outs else qm_dir / "MISSING_apo.out"
            water_out = water_outs[0] if water_outs else qm_dir / "MISSING_water.out"
            # Recover stem from filename for xyz lookup
            if la_outs:
                stem_match = re.match(r"sp_(.+)_La\.out", la_outs[0].name)
                stem = stem_match.group(1) if stem_match else qm_dir.name.removesuffix("_qm")
            else:
                stem = qm_dir.name.removesuffix("_qm")
            E_La = grab(la_out)
            E_Ca = grab(ca_out)
            E_apo = grab(apo_out)
            E_water = grab(water_out)
            term_La = is_terminated(la_out)
            term_Ca = is_terminated(ca_out)

            # Heuristic protein/site description from xyz title
            xyz = qm_dir / f"{stem}_La_qm.xyz"
            n_atoms = "?"
            if xyz.exists():
                try:
                    n_atoms = xyz.read_text().splitlines()[0].strip()
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
                "stem": qm_dir.name.removesuffix("_qm"),
                "atoms": n_atoms,
                "E_La": E_La,
                "E_Ca": E_Ca,
                "E_apo": E_apo,
                "E_water": E_water,
                "term_La": term_La,
                "term_Ca": term_Ca,
                "ddE_kcal": ddE_kcal,
                "class": cls,
            })

    # Inject LanM EF1 from earlier P3 cluster panel (lives outside *_qm/)
    lanm_la = ROOT / "qmmm" / "cluster_panel" / "La" / "opt_cluster_La.out"
    lanm_ca = ROOT / "qmmm" / "cluster_panel" / "Ca" / "opt_cluster_Ca.out"
    if lanm_la.exists() and lanm_ca.exists():
        E_La = grab(lanm_la); E_Ca = grab(lanm_ca)
        if E_La is not None and E_Ca is not None:
            ddE_kcal = ((E_Ca - E_La) - DIFF_AQUO) * HA2KCAL
            rows = [r for r in rows if r["stem"] != "lanm"]
            rows.append({
                "stem": "lanm_EF1",
                "atoms": "163",
                "E_La": E_La, "E_Ca": E_Ca,
                "E_apo": None, "E_water": None,
                "term_La": True, "term_Ca": True,
                "ddE_kcal": ddE_kcal,
                "class": classify(ddE_kcal, "163"),
            })

    # Sort by ddE_kcal, with pending at bottom
    rows.sort(key=lambda r: (r["ddE_kcal"] is None, -(r["ddE_kcal"] or 0)))

    out_tsv = ROOT / "results" / "discriminator_panel_LATEST.tsv"
    out_tsv.parent.mkdir(exist_ok=True)
    with open(out_tsv, "w") as f:
        f.write("# DFT Ca²⁺/Ln³⁺ discriminator — auto-generated\n")
        f.write(f"# E_aquo_Ca − E_aquo_La = {DIFF_AQUO:.6f} Ha (= {DIFF_AQUO*HA2KCAL:.2f} kcal/mol)\n")
        f.write("# ΔΔE(Ca − La) = (E_holo_Ca − E_holo_La) − (E_aquo_Ca − E_aquo_La)\n")
        f.write("# Class: Ln-evolved ≥ +20; Ln-preferring +5 to +20; marginal 0 to +5; ambiguous -5 to 0; Ca-evolved < -5\n")
        f.write("stem\tn_atoms\tE_La\tE_Ca\tE_apo\tE_water\tddE_kcal\tclass\n")
        for r in rows:
            ddE_s = f"{r['ddE_kcal']:+.2f}" if r["ddE_kcal"] is not None else "pending"
            E_La = f"{r['E_La']:.6f}" if r['E_La'] is not None else "-"
            E_Ca = f"{r['E_Ca']:.6f}" if r['E_Ca'] is not None else "-"
            E_apo = f"{r['E_apo']:.6f}" if r['E_apo'] is not None else "-"
            E_w = f"{r['E_water']:.6f}" if r['E_water'] is not None else "-"
            f.write(f"{r['stem']}\t{r['atoms']}\t{E_La}\t{E_Ca}\t{E_apo}\t{E_w}\t{ddE_s}\t{r['class']}\n")

    # Print summary
    print(f"\n{'stem':35s} {'atoms':>6s} {'ddE (kcal/mol)':>15s}  class")
    print("-" * 80)
    for r in rows:
        ddE_s = f"{r['ddE_kcal']:+.2f}" if r["ddE_kcal"] is not None else "pending"
        print(f"{r['stem']:35s} {r['atoms']:>6s} {ddE_s:>15s}  {r['class']}")

    n_done = sum(1 for r in rows if r["ddE_kcal"] is not None)
    n_pending = len(rows) - n_done
    print(f"\nTotal: {len(rows)} candidates, {n_done} done, {n_pending} pending")
    print(f"\nWrote: {out_tsv}")


if __name__ == "__main__":
    main()
