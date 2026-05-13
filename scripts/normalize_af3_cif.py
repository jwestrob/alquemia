#!/usr/bin/env python3
"""
Normalize AlphaFold3 mmCIF naming for our discriminator pipeline.
Detects La-containing residues and PQQ-like ligands by content (not residue name)
so it works for single-La / multi-La / La+PQQ structures alike.

Renames:
  LIG_* residue containing a La atom → "LA"   (atom name LA1 → LA)
  LIG_* residue with PQQ-like atom composition (≥14 C + ≥6 O + ≥1 N) → "PQQ"

Output: a PDB file (since downstream pipeline uses PDBFixer→PDB anyway).
"""
import sys
from pathlib import Path
import gemmi


def is_pqq_like(residue) -> bool:
    """Heuristic for PQQ: 24 atoms, 14 C + 8 O + 2 N."""
    counts = {"C": 0, "O": 0, "N": 0, "H": 0}
    for atom in residue:
        e = atom.element.name
        if e in counts:
            counts[e] += 1
    # Allow some flexibility (Hs may or may not be present, ring tautomers etc.)
    return counts["C"] >= 12 and counts["O"] >= 6 and counts["N"] >= 1


def normalize(src_cif: Path, out_pdb: Path):
    st = gemmi.read_structure(str(src_cif))
    st.setup_entities()
    n_la = 0
    n_pqq = 0
    for model in st:
        for chain in model:
            for residue in chain:
                rname = residue.name
                if not rname.startswith("LIG_") and rname not in ("LA1",):
                    continue
                # Element-based detection
                has_la = any(a.element.name == "La" for a in residue)
                if has_la:
                    residue.name = "LA"
                    for atom in residue:
                        if atom.name == "LA1":
                            atom.name = "LA"
                    n_la += 1
                elif is_pqq_like(residue):
                    residue.name = "PQQ"
                    n_pqq += 1
                else:
                    # Other unidentified ligand — leave alone (could be a cofactor
                    # we don't need to carve)
                    pass

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(out_pdb))
    print(f"normalized: {out_pdb}  La={n_la}  PQQ={n_pqq}")


if __name__ == "__main__":
    src = Path(sys.argv[1])
    out = Path(sys.argv[2])
    normalize(src, out)
