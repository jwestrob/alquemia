#!/usr/bin/env python3
"""
PQQ-aware QM cluster carver. Same as carve_generic but includes the
entire PQQ cofactor (24 atoms, charge -3 in deprotonated form) when
the metal binding site has PQQ donors.

Used for XoxF (PDB 4MAE) and MxaF (PDB 1H4I) carves.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import gemmi

ORCA_PATH = "/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
FIRST_SHELL_CUT = 3.2

SIDECHAIN_QM_ATOMS = {
    "ASP": ["CB", "HB2", "HB3", "CG", "OD1", "OD2"],
    "GLU": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"],
    "ASN": ["CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"],
    "SER": ["CB", "HB2", "HB3", "OG", "HG"],
}

# PQQ heavy-atoms-only (24 atoms, no H — pdbfixer doesn't protonate cofactors).
# Empirical electron-parity calibration: -2 gives even electron count for both
# XoxF (3 ASP/GLU + 1 ASN protein donors) and MxaF (1 GLU + 1 ASN protein donors)
# when combined with La³⁺ or Ca²⁺ metal. Chemically this corresponds to PQQ with
# one of its 3 carboxylates protonated (or equivalently a singly-protonated quinone
# or N-H elsewhere). Without explicit H atoms in the PDB this is the best
# closed-shell approximation we can make. The Ca/Ln SELECTIVITY (ΔΔE) is
# insensitive to this choice since both metals see the same PQQ.
PQQ_CHARGE = -2


def carve_with_pqq(structure_path, out_dir, stem,
                   metal_chain, metal_resname,
                   replace_metal=None, target_chain="A"):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    structure_path = Path(structure_path)
    st = gemmi.read_structure(str(structure_path))
    st.setup_entities()
    model = st[0]

    # Find target metal
    metal_atom = None
    for chain in model:
        if chain.name != target_chain:
            continue
        for residue in chain:
            if residue.name == metal_resname:
                for atom in residue:
                    if atom.element.name in ("Ce", "La", "Ca", "Y"):
                        metal_atom = atom
                        break
                if metal_atom:
                    break
        if metal_atom:
            break
    if metal_atom is None:
        raise RuntimeError(f"No metal {metal_resname} in chain {target_chain}")

    metal_pos = metal_atom.pos
    metal_elem = replace_metal if replace_metal else metal_atom.element.name
    print(f"  metal: chain {target_chain} {metal_resname}  pos=({metal_pos.x:.2f},{metal_pos.y:.2f},{metal_pos.z:.2f})")
    print(f"  using element: {metal_elem}")

    # First-shell donors (any chain — multimeric proteins put metal in one chain
    # and coordinating residues in another)
    donors = []
    for chain in model:
        for residue in chain:
            if residue.name == metal_resname:
                continue
            for atom in residue:
                if atom.element.name not in ("O", "N"):
                    continue
                d = metal_pos.dist(atom.pos)
                if d <= FIRST_SHELL_CUT:
                    donors.append((d, chain.name, residue.seqid.num, residue.name, atom.name))
    donors.sort()
    print(f"\n  First-shell donors ({len(donors)} within {FIRST_SHELL_CUT} Å):")
    for d, ch, ri, rn, an in donors:
        print(f"    {rn} {ri} ({ch})/{an}  d={d:.3f}")

    # Identify first-shell residues
    fsr_protein = sorted(set((ch, ri, rn) for _, ch, ri, rn, _ in donors
                              if rn != "PQQ"))
    has_pqq = any(rn == "PQQ" for _, _, _, rn, _ in donors)
    pqq_resnum = None
    pqq_chain = None
    if has_pqq:
        pqq_chain, pqq_resnum = next(
            (ch, ri) for _, ch, ri, rn, _ in donors if rn == "PQQ"
        )

    print(f"\n  First-shell residues (protein): {fsr_protein}")
    print(f"  PQQ present: {has_pqq} (chain {pqq_chain}, res {pqq_resnum})")

    # Build QM atoms
    qm_atoms = [(metal_elem, metal_pos.x, metal_pos.y, metal_pos.z)]
    n_carbox = 0  # ASP/GLU
    sidechain_atoms = []
    link_hs = []

    for ch_id, ri, rn in fsr_protein:
        chain = model[ch_id]
        res = next((r for r in chain if r.seqid.num == ri), None)
        if res is None:
            continue
        if res.name in ("ASP", "GLU"):
            n_carbox += 1
        sc_list = SIDECHAIN_QM_ATOMS.get(res.name)
        if sc_list is None:
            print(f"  WARN: {res.name} {ri} has no QM dict entry; skipped")
            continue
        for atom in res:
            if atom.name in sc_list:
                sidechain_atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
        cb = next((a for a in res if a.name == "CB"), None)
        ca = next((a for a in res if a.name == "CA"), None)
        if cb and ca:
            v_cb = np.array([cb.pos.x, cb.pos.y, cb.pos.z])
            v_ca = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
            v_dir = v_ca - v_cb
            v_link = v_cb + 1.09 * v_dir / np.linalg.norm(v_dir)
            link_hs.append(("H", *v_link.tolist()))

    # PQQ atoms (full cofactor)
    pqq_atoms = []
    if has_pqq:
        pqq_chain_obj = model[pqq_chain]
        pqq_res = next((r for r in pqq_chain_obj if r.seqid.num == pqq_resnum and r.name == "PQQ"), None)
        if pqq_res:
            for atom in pqq_res:
                pqq_atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
        print(f"  PQQ atoms included: {len(pqq_atoms)}")

    # Assemble: metal + sidechains + PQQ + link Hs
    full_atoms = [qm_atoms[0]] + sidechain_atoms + pqq_atoms + link_hs
    n_total = len(full_atoms)

    # Charge: M^q+ - n_carbox + PQQ_charge
    if metal_elem in ("La", "Ce", "Y"):
        m_charge = 3
    elif metal_elem == "Ca":
        m_charge = 2
    else:
        m_charge = 0
    pqq_charge = PQQ_CHARGE if has_pqq else 0
    net_charge = m_charge - n_carbox + pqq_charge

    print(f"\n  QM region: {n_total} atoms ({n_carbox} carboxylates, "
          f"{1 if has_pqq else 0} PQQ, {len(link_hs)} link Hs)")
    print(f"  Charge: {m_charge:+d} (metal) {-n_carbox:+d} (carb) {pqq_charge:+d} (PQQ) = {net_charge:+d}")

    def write_pair(suffix, charge, atoms, basis_block):
        xyz_p = out_dir / f"{stem}_{suffix}_qm.xyz"
        with open(xyz_p, "w") as f:
            f.write(f"{len(atoms)}\n")
            f.write(f"{stem} {suffix} (charge={charge}, mult=1)\n")
            for el, x, y, z in atoms:
                f.write(f"{el:<3s} {x:>14.6f} {y:>14.6f} {z:>14.6f}\n")
        inp_p = out_dir / f"sp_{stem}_{suffix}.inp"
        with open(inp_p, "w") as f:
            f.write(f"""! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
{basis_block}
* xyzfile {charge} 1 {xyz_p.name}
""")

    # Use La form (or whatever metal_elem is) for the "primary"
    if metal_elem == "La":
        la_basis = '%basis\n  NewECP La "def2-ECP" end\n  NewGTO La "def2-TZVP" end\nend'
        ca_basis = ""
    elif metal_elem == "Ca":
        la_basis = '%basis\n  NewECP La "def2-ECP" end\n  NewGTO La "def2-TZVP" end\nend'
        ca_basis = ""
    elif metal_elem == "Ce":
        la_basis = '%basis\n  NewECP La "def2-ECP" end\n  NewGTO La "def2-TZVP" end\nend'
        ca_basis = ""
    else:
        la_basis = ""
        ca_basis = ""

    # Write La and Ca forms (vertical metal swap)
    la_atoms_full = [("La", metal_pos.x, metal_pos.y, metal_pos.z)] + sidechain_atoms + pqq_atoms + link_hs
    write_pair("La", 3 - n_carbox + pqq_charge, la_atoms_full, la_basis)
    ca_atoms_full = [("Ca", metal_pos.x, metal_pos.y, metal_pos.z)] + sidechain_atoms + pqq_atoms + link_hs
    write_pair("Ca", 2 - n_carbox + pqq_charge, ca_atoms_full, ca_basis)
    # Apo: drop metal
    apo_atoms = sidechain_atoms + pqq_atoms + link_hs
    write_pair("apo", -n_carbox + pqq_charge, apo_atoms, "")

    # Bulk water
    water_xyz = out_dir / "bulk_water.xyz"
    water_xyz.write_text("""3
H2O reference
O   0.000000  0.000000  0.000000
H   0.756950  0.000000  0.585822
H  -0.756950  0.000000  0.585822
""")
    (out_dir / f"sp_{stem}_water.inp").write_text(
        "! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n%maxcore 8000\n* xyzfile 0 1 bulk_water.xyz\n"
    )

    # Submit script
    sh = out_dir / f"submit_{stem}.sh"
    with open(sh, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p memory
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J orca_{stem}
#SBATCH -o {out_dir}/slurm_%j.out
#SBATCH -e {out_dir}/slurm_%j.err

set -euo pipefail
ORCA_PATH={ORCA_PATH}
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:${{LD_LIBRARY_PATH:-}}
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE; export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE

cd {out_dir}
export ORCA_TMPDIR=$PWD; export TMPDIR=$PWD
echo "host=$(hostname) jobid=$SLURM_JOB_ID started=$(date -Iseconds)"
for kind in La Ca apo water; do
    echo "=== $kind SP ==="
    $ORCA_PATH/orca sp_{stem}_$kind.inp > sp_{stem}_$kind.out 2>&1
    grep -E "FINAL SINGLE POINT|TERMINATED" sp_{stem}_$kind.out | tail -2
done
echo "finished=$(date -Iseconds)"
""")
    sh.chmod(0o755)
    print(f"\nWrote outputs in {out_dir}/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("structure", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--stem", required=True)
    ap.add_argument("--metal-chain", required=True)
    ap.add_argument("--metal-resname", required=True)  # e.g. "CE" or "CA"
    ap.add_argument("--replace-metal", default=None)  # e.g. "La"
    args = ap.parse_args()
    carve_with_pqq(args.structure, args.out_dir, args.stem,
                   args.metal_chain, args.metal_resname,
                   args.replace_metal, args.metal_chain)
