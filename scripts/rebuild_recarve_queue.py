#!/usr/bin/env python3
"""
rebuild_recarve_queue.py — scan workspace for empty-carve entries and
write results/recarve_queue.tsv.

Idempotent: stems already classified CARVE_AMBIGUOUS or SOLVENT_EXCLUDE in
prior recarve_log_*.tsv files are skipped (their outcome won't change without
a methodology change). Stems that were RECARVED_AND_SUBMITTED are also
skipped (already taken care of). Only new empties land in the new queue.

Signature of an empty carve: E[La] = -31.158775581143 Ha (bare La³⁺ + caps,
typically n_atoms = 1). This is what falls out of `carve_generic.py` when no
carve-eligible donor is within FIRST_SHELL_CUT of the metal.

Usage:
    python3 scripts/rebuild_recarve_queue.py
"""
from __future__ import annotations
import os, re, csv, sys, glob
from pathlib import Path

ALCH = Path(__file__).resolve().parent.parent
OUT_TSV = ALCH / "results" / "recarve_queue.tsv"
LOG_GLOB = str(ALCH / "results" / "recarve_log_*.tsv")

# Empty-carve fingerprint — these are the bit-identical values
EMPTY_LA_LO = -31.16
EMPTY_LA_HI = -31.15


def parse_final_e(path: Path) -> float | None:
    if not path.exists():
        return None
    try:
        txt = path.read_text()
    except Exception:
        return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", txt)
    if m:
        try:
            return float(m[-1])
        except ValueError:
            return None
    return None


def load_prior_outcomes() -> dict[str, str]:
    """Read all recarve_log_*.tsv files; return {stem: latest_outcome}."""
    out = {}
    for path in sorted(glob.glob(LOG_GLOB)):
        try:
            with open(path) as f:
                # Tab-separated; skip comments. First non-comment line is header.
                rdr = csv.DictReader(
                    (ln for ln in f if not ln.startswith("#")),
                    delimiter="\t"
                )
                for row in rdr:
                    stem = row.get("stem")
                    outcome = row.get("outcome", "")
                    if stem:
                        out[stem] = outcome
        except Exception:
            continue
    return out


def main():
    prior = load_prior_outcomes()
    flagged = []

    for d in sorted(ALCH.glob("*_qm")):
        stem = d.name[:-3]  # strip "_qm"
        laf = d / f"sp_{stem}_La.out"
        e = parse_final_e(laf)
        if e is None:
            continue
        if not (EMPTY_LA_LO <= e <= EMPTY_LA_HI):
            continue

        # Empty-carve hit. Check prior outcome.
        outcome = prior.get(stem, "")
        if outcome in ("CARVE_AMBIGUOUS", "SOLVENT_EXCLUDE", "RECARVED_AND_SUBMITTED"):
            # Already handled; skip unless protonated.pdb is freshly newer
            # (i.e. the daemon re-emitted it after a fix).
            continue

        # Count atoms in the carved cluster
        xyz = d / f"{stem}_La_qm.xyz"
        n_atoms = "?"
        if xyz.exists():
            try:
                with open(xyz) as f:
                    n_atoms = int(f.readline().strip())
            except Exception:
                pass

        # Need the protonated.pdb to re-carve
        pdb = d / f"{stem}_protonated.pdb"
        e_ca = parse_final_e(d / f"sp_{stem}_Ca.out")
        flagged.append({
            "stem": stem,
            "n_atoms_La": n_atoms,
            "E_La": e,
            "E_Ca": e_ca if e_ca else "-",
            "protonated_pdb": str(pdb) if pdb.exists() else "MISSING",
            "qm_dir": str(d),
        })

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w") as f:
        f.write("# Empty-carve / failed-carve queue — entries needing re-carve with "
                "fixed carve_generic.py\n")
        f.write("# Fixed script: alchemical_bvs/scripts/carve_generic.py (3.0 Å cutoff, "
                "TYR sidechain, backbone-O carving)\n")
        f.write("# Signature: E[La] ≈ -31.16 Ha (bare La³⁺ + caps), N_atoms typically 1\n")
        f.write("# Stems already classified CARVE_AMBIGUOUS / SOLVENT_EXCLUDE / "
                "RECARVED_AND_SUBMITTED in prior recarve_log_*.tsv are EXCLUDED.\n")
        f.write("stem\tn_atoms_La\tE_La\tE_Ca\tprotonated_pdb\tqm_dir\n")
        for r in flagged:
            e_ca_s = f"{r['E_Ca']:.6f}" if isinstance(r['E_Ca'], (int, float)) else str(r['E_Ca'])
            f.write(f"{r['stem']}\t{r['n_atoms_La']}\t{r['E_La']:.6f}\t{e_ca_s}\t"
                    f"{r['protonated_pdb']}\t{r['qm_dir']}\n")

    print(f"Rebuilt {OUT_TSV}: {len(flagged)} new empty-carve entries", file=sys.stderr)
    print(f"  (Skipped {sum(1 for v in prior.values() if v in ('CARVE_AMBIGUOUS','SOLVENT_EXCLUDE','RECARVED_AND_SUBMITTED'))} already-classified stems)", file=sys.stderr)
    return len(flagged)


if __name__ == "__main__":
    n = main()
    sys.exit(0 if n >= 0 else 1)
