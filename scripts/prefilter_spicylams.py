#!/usr/bin/env python3
"""
Pre-filter spicy_lams CIFs for the discriminator pipeline.

Sharur's "CN ≥ 7" filter counts backbone-O donors. Our carve_generic.py
only includes sidechains for ASP/GLU/ASN/GLN/SER. So candidates whose
first-shell donors are all backbone-O (e.g., GLY/ALA/VAL coordination)
produce 1-atom QM clusters and are useless to score.

This script reads each CIF, finds La, counts ASP/GLU/ASN/GLN/SER
within 3.2 Å, and outputs candidates with ≥2 such "carve-friendly"
residues.

Usage:
  python prefilter_spicylams.py <inbox_dir> [min_donors]
  -> writes inbox_dir/_filter_passed.txt and inbox_dir/_filter_rejected.txt
"""
from __future__ import annotations
import sys
from pathlib import Path
import gemmi

CARVE_DICT = {"ASP", "GLU", "ASN", "GLN", "SER"}
CUTOFF = 3.2


def first_shell_dict_donors(cif_path: Path) -> tuple[int, int, list[str]]:
    """Return (n_dict_donors, n_total_donors, donor_summary)."""
    try:
        st = gemmi.read_structure(str(cif_path))
        st.setup_entities()
    except Exception as e:
        return -1, -1, [f"READ_ERROR: {e}"]

    model = st[0]
    la_atoms = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.name == "La":
                    la_atoms.append(atom)

    if not la_atoms:
        return 0, 0, ["NO_LA"]

    # Use first La (multi-La rare in this set)
    la_pos = la_atoms[0].pos

    dict_residues = set()
    total_residues = set()
    for chain in model:
        for residue in chain:
            if residue.name in ("LA", "PQQ"):
                continue
            for atom in residue:
                if atom.element.name not in ("O", "N"):
                    continue
                if la_pos.dist(atom.pos) <= CUTOFF:
                    rkey = (chain.name, residue.seqid.num, residue.name)
                    total_residues.add(rkey)
                    if residue.name in CARVE_DICT:
                        dict_residues.add(rkey)

    summary = sorted(total_residues, key=lambda x: x[1])
    return len(dict_residues), len(total_residues), [f"{r[2]}{r[1]}" for r in summary]


def main(inbox: Path, min_donors: int = 2):
    cifs = sorted(inbox.glob("*.cif"))
    print(f"Found {len(cifs)} CIFs in inbox", file=sys.stderr)

    pass_p = inbox / "_filter_passed.txt"
    reject_p = inbox / "_filter_rejected.txt"
    with open(pass_p, "w") as fp, open(reject_p, "w") as fr:
        n_pass = n_rej = 0
        for cif in cifs:
            n_dict, n_total, donors = first_shell_dict_donors(cif)
            short = cif.name
            line = f"{short}\tn_dict={n_dict}\tn_total={n_total}\tdonors={','.join(donors)}\n"
            if n_dict >= min_donors:
                fp.write(line)
                n_pass += 1
            else:
                fr.write(line)
                n_rej += 1
            if (n_pass + n_rej) % 100 == 0:
                print(f"  scanned {n_pass + n_rej} ...", file=sys.stderr)

    print(f"\nPASS: {n_pass}  (≥{min_donors} dict-eligible donors)", file=sys.stderr)
    print(f"REJ:  {n_rej}  (<{min_donors} dict-eligible donors)", file=sys.stderr)
    print(f"Wrote: {pass_p}  and  {reject_p}", file=sys.stderr)


if __name__ == "__main__":
    inbox = Path(sys.argv[1])
    min_donors = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    main(inbox, min_donors)
