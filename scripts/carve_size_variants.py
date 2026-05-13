#!/usr/bin/env python3
"""
Build cluster-size variants for the sensitivity figure.
Takes a protonated PDB + a metal site, produces xyz files at multiple QM region sizes.

Usage: python carve_size_variants.py <pdb> <out_dir> --stem calexcitin

Sizes produced:
  size16 (~16 atoms): metal + 4 carboxylate "COO⁻" groups only (4 sidechain cuts at Cβ)
  size32 (current default): full sidechains + link Hs (carve_generic equivalent)
  size50 (~50 atoms): include Cα + backbone N for first-shell residues
  size100 (~100 atoms): include 2nd-shell residues (any with sidechain heavy atom within 5 Å of metal)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import gemmi

ORCA_PATH = "/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
FIRST_SHELL_CUT = 3.2
SECOND_SHELL_CUT = 5.5

SIDECHAIN_QM_ATOMS = {
    "ASP": ["CB", "HB2", "HB3", "CG", "OD1", "OD2"],
    "GLU": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"],
    "ASN": ["CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"],
    "SER": ["CB", "HB2", "HB3", "OG", "HG"],
}
# Just the carboxylate (smallest)
COO_ONLY = {
    "ASP": ["CG", "OD1", "OD2"],
    "GLU": ["CD", "OE1", "OE2"],
    "ASN": ["CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CD", "OE1", "NE2", "HE21", "HE22"],
    "SER": ["OG", "HG"],
}
BACKBONE = ["CA", "HA", "N", "H", "C", "O"]


def find_first_shell(model, la_pos, cut=FIRST_SHELL_CUT):
    donors = []
    for chain in model:
        for residue in chain:
            if residue.name in ("LA", "LA3", "CA", "CE"):
                continue
            for atom in residue:
                if atom.element.name not in ("O", "N"):
                    continue
                if la_pos.dist(atom.pos) <= cut:
                    donors.append((chain.name, residue.seqid.num, residue.name, atom.name))
    return sorted(set((c, r, n) for c, r, n, _ in donors))


def find_second_shell(model, la_pos, first_shell):
    """Residues with any sidechain heavy atom within SECOND_SHELL_CUT, not in first shell."""
    fs_keys = {(c, r) for c, r, n in first_shell}
    second = []
    for chain in model:
        for residue in chain:
            if (chain.name, residue.seqid.num) in fs_keys:
                continue
            if residue.name in ("LA", "LA3", "CA", "CE", "WAT", "HOH"):
                continue
            min_d = 999
            for atom in residue:
                if atom.name in ("CA", "C", "N", "O", "H"):
                    continue
                if atom.element.name == "H":
                    continue
                min_d = min(min_d, la_pos.dist(atom.pos))
            if FIRST_SHELL_CUT < min_d <= SECOND_SHELL_CUT:
                second.append((chain.name, residue.seqid.num, residue.name))
    return second


def link_h(cb, ca):
    """Place an H between Cβ-Cα at 1.09 Å from Cβ."""
    v_cb = np.array([cb.pos.x, cb.pos.y, cb.pos.z])
    v_ca = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
    v_dir = v_ca - v_cb
    v = v_cb + 1.09 * v_dir / np.linalg.norm(v_dir)
    return ("H", *v.tolist())


def link_h_at(parent_pos, neighbor_pos):
    """Generic: place H along parent-neighbor at 1.09 Å from parent."""
    v_p = np.array([parent_pos.x, parent_pos.y, parent_pos.z])
    v_n = np.array([neighbor_pos.x, neighbor_pos.y, neighbor_pos.z])
    d = v_n - v_p
    v = v_p + 1.09 * d / np.linalg.norm(d)
    return ("H", *v.tolist())


def carve(pdb_path, out_dir, stem):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()
    model = st[0]

    # Find all La sites, pick one with most carboxylate donors (matches carve_generic)
    la_candidates = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.name == "La":
                    la_candidates.append((chain.name, residue.seqid.num, atom))
    if not la_candidates:
        raise RuntimeError(f"No La in {pdb_path}")
    best = None; best_n = -1
    for c in la_candidates:
        la_pos = c[2].pos
        n = 0
        for chain in model:
            for residue in chain:
                if residue.name in ("ASP", "GLU"):
                    for at in residue:
                        if at.name in ("OD1", "OD2", "OE1", "OE2") and la_pos.dist(at.pos) <= FIRST_SHELL_CUT:
                            n += 1
                            break
        if n > best_n:
            best_n = n
            best = c
    la_atom = best[2]
    la_pos = la_atom.pos
    print(f"  picked best site: chain {best[0]} res {best[1]} ({best_n} carbox donors)")
    print(f"  La position: ({la_pos.x:.3f}, {la_pos.y:.3f}, {la_pos.z:.3f})")

    fsr = find_first_shell(model, la_pos)
    print(f"First-shell residues: {fsr}")
    second = find_second_shell(model, la_pos, fsr)
    print(f"Second-shell residues: {second}")

    # Build size variants
    def build(sidechain_dict, include_backbone=False, include_second=False, label="size"):
        atoms = [("La", la_pos.x, la_pos.y, la_pos.z)]
        n_carbox = 0
        seen_atoms = set()  # by (chain, resnum, atomname)
        # First-shell residues
        for ch_id, ri, rn in fsr:
            if rn in ("ASP", "GLU"):
                n_carbox += 1
            chain = model[ch_id]
            res = next((r for r in chain if r.seqid.num == ri), None)
            if not res: continue
            sc_list = sidechain_dict.get(res.name, [])
            for atom in res:
                aname = atom.name
                if aname in sc_list or (include_backbone and aname in BACKBONE):
                    key = (ch_id, ri, aname)
                    if key in seen_atoms: continue
                    seen_atoms.add(key)
                    atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
            # Link H placement
            if include_backbone:
                # backbone N-Cα link H, C-O link H
                ca = next((a for a in res if a.name == "CA"), None)
                n = next((a for a in res if a.name == "N"), None)
                c = next((a for a in res if a.name == "C"), None)
                if ca and n:
                    atoms.append(link_h_at(n.pos, ca.pos))  # actually N→prev-C would be more correct; approximate
                if ca and c:
                    atoms.append(link_h_at(c.pos, ca.pos))
            else:
                cb = next((a for a in res if a.name == "CB"), None)
                ca = next((a for a in res if a.name == "CA"), None)
                if cb and ca and any(aname in sc_list for aname in [a.name for a in res]):
                    if "CB" in sc_list:
                        atoms.append(link_h(cb, ca))  # H at Cβ-Cα for sidechain-only carve
                    else:
                        # CG-only carve: link H at CG-CB
                        cg = next((a for a in res if a.name == "CG"), None)
                        if cg and cb:
                            atoms.append(link_h_at(cg.pos, cb.pos))
        # Second-shell residues
        if include_second:
            for ch_id, ri, rn in second:
                chain = model[ch_id]
                res = next((r for r in chain if r.seqid.num == ri), None)
                if not res: continue
                sc_list = SIDECHAIN_QM_ATOMS.get(res.name, [])
                if not sc_list:
                    continue
                if res.name in ("ASP", "GLU"):
                    n_carbox += 1
                for atom in res:
                    if atom.name in sc_list:
                        atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
                cb = next((a for a in res if a.name == "CB"), None)
                ca = next((a for a in res if a.name == "CA"), None)
                if cb and ca:
                    atoms.append(link_h(cb, ca))

        return atoms, n_carbox

    sizes = [
        ("size16", COO_ONLY, False, False),
        ("size32", SIDECHAIN_QM_ATOMS, False, False),
        ("size50", SIDECHAIN_QM_ATOMS, True, False),
        ("size100", SIDECHAIN_QM_ATOMS, True, True),
    ]

    for label, scd, bb, second_shell in sizes:
        atoms, n_carbox = build(scd, bb, second_shell, label)
        net = 3 - n_carbox
        print(f"  {label}: {len(atoms)} atoms, {n_carbox} carbox, charge={net:+d}")
        # Write La and Ca xyz
        for elem, charge in [("La", net), ("Ca", net - 1)]:
            xyz = out_dir / f"{stem}_{label}_{elem}.xyz"
            with open(xyz, "w") as f:
                f.write(f"{len(atoms)}\n")
                f.write(f"{stem} {label} {elem} cluster (charge={charge})\n")
                for at in atoms:
                    el = elem if at[0] == "La" else at[0]  # swap metal
                    f.write(f"{el:<3s} {at[1]:>14.6f} {at[2]:>14.6f} {at[3]:>14.6f}\n")
            inp = out_dir / f"sp_{stem}_{label}_{elem}.inp"
            la_block = '%basis\n  NewECP La "def2-ECP" end\n  NewGTO La "def2-TZVP" end\nend' if elem == "La" else ""
            with open(inp, "w") as f:
                f.write(f"""! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
{la_block}
* xyzfile {charge} 1 {xyz.name}
""")

    # Single submit script for ALL size variants
    sh = out_dir / f"submit_size_variants.sh"
    with open(sh, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p memory
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J size_{stem}
#SBATCH -o {out_dir}/slurm_%j.out
set -uo pipefail
ORCA_PATH={ORCA_PATH}
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:${{LD_LIBRARY_PATH:-}}
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
cd {out_dir}
export ORCA_TMPDIR=$PWD; export TMPDIR=$PWD
echo "host=$(hostname) jobid=$SLURM_JOB_ID started=$(date -Iseconds)"
""")
        for label, _, _, _ in sizes:
            for elem in ("La", "Ca"):
                f.write(f"""echo "=== {label} {elem} ==="
$ORCA_PATH/orca sp_{stem}_{label}_{elem}.inp > sp_{stem}_{label}_{elem}.out 2>&1 || true
grep "FINAL SINGLE POINT" sp_{stem}_{label}_{elem}.out | tail -1 || echo "(no FINAL SP)"
""")
        f.write('echo "finished=$(date -Iseconds)"\n')
    sh.chmod(0o755)
    print(f"\nWrote outputs in {out_dir}/")
    print(f"Submit with: sbatch {sh}")


if __name__ == "__main__":
    pdb = Path(sys.argv[1])
    out = Path(sys.argv[2])
    stem = sys.argv[3] if len(sys.argv) > 3 else "test"
    carve(pdb, out, stem)
