#!/usr/bin/env python3
"""
Generic QM cluster carver for the Ca/Ln DFT discriminator panel.
Outputs ORCA inputs for {La, Ca, apo} + a SLURM submit script.
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


def carve(structure_path: Path, out_dir: Path, stem: str,
          site_chain=None, site_resnum=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    st = gemmi.read_structure(str(structure_path))
    st.setup_entities()
    model = st[0]

    la_candidates = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.name in ("La", "Y", "Ca", "Ce"):
                    la_candidates.append((chain.name, residue.seqid.num, residue.name, atom))

    if not la_candidates:
        raise RuntimeError(f"No La/Y/Ca/Ce in {structure_path}")

    if site_chain and site_resnum:
        target = next((c for c in la_candidates
                       if c[0] == site_chain and c[1] == site_resnum), None)
    else:
        # Pick site with most carboxylate donors
        best = None
        best_n = -1
        for c in la_candidates:
            la_pos = c[3].pos
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
        target = best
        print(f"  picked best site: chain {target[0]} res {target[1]} ({best_n} carbox donors)")

    target_chain, target_resnum, target_resname, target_atom = target
    la_pos = target_atom.pos
    print(f"  metal at chain {target_chain} {target_resname} {target_resnum}: ({la_pos.x:.3f},{la_pos.y:.3f},{la_pos.z:.3f})")

    # First-shell donors
    donors = []
    for chain in model:
        for residue in chain:
            if residue.name in ("LA", "LA3", "YT3", "Y", "CA", "CE"):
                continue
            for atom in residue:
                if atom.element.name not in ("O", "N"):
                    continue
                d = la_pos.dist(atom.pos)
                if d <= FIRST_SHELL_CUT:
                    donors.append((d, chain.name, residue.seqid.num, residue.name, atom.name))
    donors.sort()
    print(f"\n  First-shell donors ({len(donors)} within {FIRST_SHELL_CUT} Å):")
    for d, ch, ri, rn, an in donors:
        print(f"    {rn} {ri} ({ch})/{an}  d={d:.3f} Å")

    # Inner-shell water by O position
    inner_w_O = []
    for chain in model:
        for residue in chain:
            if residue.name not in ("WAT", "HOH"):
                continue
            for atom in residue:
                if atom.name == "O":
                    if la_pos.dist(atom.pos) <= FIRST_SHELL_CUT:
                        inner_w_O.append(atom.pos)
    print(f"  Inner-shell waters: {len(inner_w_O)}")

    fsr = sorted(set((ch, ri) for _, ch, ri, _, _ in donors))
    print(f"  First-shell residues: {fsr}")

    # Build QM list, tracking sidechain count for apo carve
    n_sc = 0
    n_carbox = 0
    sidechain_atoms = []
    link_hs = []
    for ch_id, ri in fsr:
        chain = model[ch_id]
        res = next((r for r in chain if r.seqid.num == ri), None)
        if res is None:
            continue
        if res.name in ("ASP", "GLU"):
            n_carbox += 1
        sc_list = SIDECHAIN_QM_ATOMS.get(res.name)
        if sc_list is None:
            print(f"  WARN: residue {res.name} {ri} not in dict; skipped")
            continue
        for atom in res:
            if atom.name in sc_list:
                sidechain_atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
                n_sc += 1
        # Link H at Cβ-Cα
        cb = next((a for a in res if a.name == "CB"), None)
        ca = next((a for a in res if a.name == "CA"), None)
        if cb and ca:
            v_cb = np.array([cb.pos.x, cb.pos.y, cb.pos.z])
            v_ca = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
            v_dir = v_ca - v_cb
            v_link = v_cb + 1.09 * v_dir / np.linalg.norm(v_dir)
            link_hs.append(("H", *v_link.tolist()))

    # Inner water atoms
    water_atoms = []
    for o_pos in inner_w_O:
        water_atoms.append(("O", o_pos.x, o_pos.y, o_pos.z))
        h_count = 0
        for chain in model:
            for residue in chain:
                if residue.name not in ("WAT", "HOH"):
                    continue
                for atom in residue:
                    if atom.element.name == "H" and o_pos.dist(atom.pos) <= 1.2:
                        water_atoms.append(("H", atom.pos.x, atom.pos.y, atom.pos.z))
                        h_count += 1
                        if h_count == 2: break
                if h_count >= 2: break
            if h_count >= 2: break

    # Construct full holo: [La, sidechains..., waters..., link_Hs...]
    n_total = 1 + len(sidechain_atoms) + len(water_atoms) + len(link_hs)
    net_charge = 3 - n_carbox
    print(f"\n  QM region: {n_total} atoms ({n_carbox} carboxylates, {len(inner_w_O)} water, {len(link_hs)} link H)")
    print(f"  Net charge (La): {net_charge:+d}")

    def write_pair(suffix, charge, atoms, metal_basis_block):
        xyz_p = out_dir / f"{stem}_{suffix}_qm.xyz"
        with open(xyz_p, "w") as f:
            f.write(f"{len(atoms)}\n")
            f.write(f"{stem} {suffix} cluster (charge={charge}, mult=1)\n")
            for el, x, y, z in atoms:
                f.write(f"{el:<3s} {x:>14.6f} {y:>14.6f} {z:>14.6f}\n")
        inp_p = out_dir / f"sp_{stem}_{suffix}.inp"
        with open(inp_p, "w") as f:
            f.write(f"""# {stem} {suffix}
! r2SCAN-3c CPCM(Water) DefGrid3

%maxcore 8000

{metal_basis_block}
* xyzfile {charge} 1 {xyz_p.name}
""")

    # La holo
    la_atoms = [("La", la_pos.x, la_pos.y, la_pos.z)] + sidechain_atoms + water_atoms + link_hs
    write_pair("La", net_charge, la_atoms,
               '%basis\n  NewECP La "def2-ECP" end\n  NewGTO La "def2-TZVP" end\nend')
    # Ca holo (metal swap)
    ca_atoms = [("Ca", la_pos.x, la_pos.y, la_pos.z)] + sidechain_atoms + water_atoms + link_hs
    write_pair("Ca", net_charge - 1, ca_atoms, "")
    # Apo (no metal, no inner water)
    apo_atoms = sidechain_atoms + link_hs
    write_pair("apo", net_charge - 3, apo_atoms, "")

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
for kind in La Ca apo; do
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
    ap.add_argument("--site-chain", default=None)
    ap.add_argument("--site-resnum", type=int, default=None)
    args = ap.parse_args()
    carve(args.structure, args.out_dir, args.stem, args.site_chain, args.site_resnum)
