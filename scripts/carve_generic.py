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
FIRST_SHELL_CUT = 3.0  # 2026-05-11: tightened from 3.2 → 3.0 Å (per Jacob)

# Sidechain atom lists for the carve. Only residues here have sidechain contacts
# pulled into the QM region (along with a link H at Cβ to cap Cβ-Cα).
# Note: TYR added 2026-05-11 to fix carve failures on Tyr-OH-coordinated sites.
# Note: 2026-05-11 (later) — an agent attempted to add THR/CYS/HIS/MET/LYS/ARG,
#       but per Jacob's explicit direction those are reverted out. Only TYR
#       was authorized in addition to the original carboxylate/amide/Ser set.
#       If a re-carve still produces a 1-atom cluster after the TYR+backbone-O
#       fix, classify it as SOLVENT_EXCLUDE / CARVE_AMBIGUOUS rather than
#       expanding the dict.
SIDECHAIN_QM_ATOMS = {
    "ASP": ["CB", "HB2", "HB3", "CG", "OD1", "OD2"],
    "GLU": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"],
    "ASN": ["CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"],
    "SER": ["CB", "HB2", "HB3", "OG", "HG"],
    "TYR": ["CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1",
            "CZ", "OH", "HH", "CE2", "HE2", "CD2", "HD2"],
}


def carve(structure_path: Path, out_dir: Path, stem: str,
          site_chain=None, site_resnum=None, first_shell_cut=None):
    # Allow caller to override the module-level default per-call without
    # mutating module state. Falls back to FIRST_SHELL_CUT if not given.
    global FIRST_SHELL_CUT
    if first_shell_cut is not None:
        FIRST_SHELL_CUT = float(first_shell_cut)
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

    # First-shell donors.
    # 2026-05-11: include S (Cys-SG, Met-SD) in donor detection — was O/N only,
    # which silently dropped all Cys-coordinated sites to bare-La carves.
    donors = []
    for chain in model:
        for residue in chain:
            if residue.name in ("LA", "LA3", "YT3", "Y", "CA", "CE"):
                continue
            for atom in residue:
                if atom.element.name not in ("O", "N", "S"):
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

    # Identify which residues have a backbone-O contact (for carbonyl carving).
    # Added 2026-05-11 to capture GLY/PRO and other non-sidechain-coordinated sites.
    # A residue gets a backbone-O carve if the *backbone* O atom is within FIRST_SHELL_CUT.
    backbone_o_residues = set()
    for d, ch_id, ri, rn, an in donors:
        if an == "O":   # backbone carbonyl oxygen
            backbone_o_residues.add((ch_id, ri))

    # Build QM list, tracking sidechain count for apo carve
    n_sc = 0
    n_carbox = 0
    n_backbone = 0
    sidechain_atoms = []
    link_hs = []
    for ch_id, ri in fsr:
        chain = model[ch_id]
        res = next((r for r in chain if r.seqid.num == ri), None)
        if res is None:
            continue
        if res.name in ("ASP", "GLU"):
            n_carbox += 1

        # --- (A) sidechain carve, if residue has a dict entry ---
        sc_list = SIDECHAIN_QM_ATOMS.get(res.name)
        if sc_list is not None:
            for atom in res:
                if atom.name in sc_list:
                    sidechain_atoms.append((atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z))
                    n_sc += 1
            # Link H at Cβ-Cα (caps the broken Cβ-Cα bond on the dropped backbone)
            cb = next((a for a in res if a.name == "CB"), None)
            ca = next((a for a in res if a.name == "CA"), None)
            if cb and ca:
                v_cb = np.array([cb.pos.x, cb.pos.y, cb.pos.z])
                v_ca = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
                v_dir = v_ca - v_cb
                v_link = v_cb + 1.09 * v_dir / np.linalg.norm(v_dir)
                link_hs.append(("H", *v_link.tolist()))
        else:
            print(f"  WARN: residue {res.name} {ri} not in sidechain dict; sidechain skipped")

        # --- (B) backbone-O carve, if this residue's backbone O is in first shell ---
        # Independent of sidechain dict — captures GLY/PRO/etc. carbonyl coordination.
        if (ch_id, ri) in backbone_o_residues:
            c_atom = next((a for a in res if a.name == "C"), None)
            o_atom = next((a for a in res if a.name == "O"), None)
            ca_atom = next((a for a in res if a.name == "CA"), None)
            if c_atom and o_atom and ca_atom:
                # Skip if we already added these via sidechain dict (shouldn't happen — backbone C/O not in any SIDECHAIN_QM_ATOMS list)
                sidechain_atoms.append(("C", c_atom.pos.x, c_atom.pos.y, c_atom.pos.z))
                sidechain_atoms.append(("O", o_atom.pos.x, o_atom.pos.y, o_atom.pos.z))
                n_backbone += 1
                # Link H from C → toward Cα (caps the broken C-Cα peptide bond)
                v_c = np.array([c_atom.pos.x, c_atom.pos.y, c_atom.pos.z])
                v_ca = np.array([ca_atom.pos.x, ca_atom.pos.y, ca_atom.pos.z])
                v_dir = v_ca - v_c
                v_link = v_c + 1.09 * v_dir / np.linalg.norm(v_dir)
                link_hs.append(("H", *v_link.tolist()))
                # Link H from C → away from O (caps the broken C-N(next) peptide bond)
                # Use the direction perpendicular to C=O and C-Cα, roughly in-plane
                v_o = np.array([o_atom.pos.x, o_atom.pos.y, o_atom.pos.z])
                v_co = v_o - v_c
                v_cca = v_ca - v_c
                # Bisector of the (-O, -Cα) directions gives the "where N should be" direction
                v_n_dir = -(v_co / np.linalg.norm(v_co)) - (v_cca / np.linalg.norm(v_cca))
                v_n_dir = v_n_dir / np.linalg.norm(v_n_dir)
                v_link_n = v_c + 1.09 * v_n_dir
                link_hs.append(("H", *v_link_n.tolist()))

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

    # ─── CN @3.0 Å filter (N+O donors) ─────────────────────────────────────
    # Sites with too few first-shell donors are not credible Ln binders;
    # running DFT on them wastes cluster time. Threshold: skip if CN ≤ 3.
    # See VALIDATION.md / RESULTS.md for the cross-table that motivated this.
    CN_MIN = 4
    cn_3A = 0
    for el, x, y, z in sidechain_atoms + water_atoms:
        if el not in ("O", "N"):
            continue
        dx = x - la_pos.x; dy = y - la_pos.y; dz = z - la_pos.z
        if (dx*dx + dy*dy + dz*dz) ** 0.5 <= 3.0:
            cn_3A += 1
    print(f"  CN @3Å (N+O): {cn_3A}")
    if cn_3A < CN_MIN:
        skip_marker = out_dir / f"{stem}_SKIPPED_CN{cn_3A}.txt"
        with open(skip_marker, "w") as f:
            f.write(f"# Skipped DFT pre-submit (CN @3Å = {cn_3A} < {CN_MIN})\n")
            f.write(f"# stem={stem}\n# n_carbox={n_carbox}\n# n_total_qm={n_total}\n")
        print(f"  ⏭️  Skipping ORCA writes — CN={cn_3A} < {CN_MIN}. Marker: {skip_marker.name}")
        return  # exit before write_pair / submit script generation
    # ───────────────────────────────────────────────────────────────────────

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
for kind in La Ca apo water; do
    echo "=== $kind SP ==="
    $ORCA_PATH/orca sp_{stem}_$kind.inp > sp_{stem}_$kind.out 2>&1
    grep -E "FINAL SINGLE POINT|TERMINATED" sp_{stem}_$kind.out | tail -2
done
echo "finished=$(date -Iseconds)"

# Auto-cleanup on success: if all 4 SPs terminated normally, drop ORCA scratch + slurm logs.
# Keeps .out (energies + diagnostics), .inp (recipe), .xyz (carves), .pdb (protonated).
ALL_OK=1
for kind in La Ca apo water; do
    grep -q "ORCA TERMINATED NORMALLY" sp_{stem}_$kind.out 2>/dev/null || ALL_OK=0
done
if [ "$ALL_OK" = "1" ]; then
    rm -f *.gbw *.densities *.densitiesinfo *.bas[0-9] *.bas[0-9][0-9] *.tmp *.tmp.[0-9]* \\
          *.bibtex *.cpcm *.cpcm_corr *.cpcm_achol *.SHARK*.tmp *.shark_grid.tmp *.SHARKINP.tmp \\
          *.ges *.property.txt *.K.tmp *.J.tmp *.E.tmp *.S.tmp *.T.tmp *.V.tmp *.H.tmp \\
          *.P0.tmp *.G0.tmp *.PINP[0-9].tmp *.POLD[0-9].tmp *.PAUX.tmp *.PDAT.tmp \\
          *.VAUXJ.tmp *.VCDJ.tmp *.VEXT.tmp *.VM1EXT.tmp *.VXC[0-9].tmp \\
          *.SP12.tmp *.SM12.tmp *.SQRTVEXT.tmp *.SQRTVJ.tmp \\
          *.FAO[0-9].tmp *.FMO[0-9].tmp *.EIJ.tmp *.cpscfdata.tmp.* \\
          *.diisao.tmp.* *.diise.tmp.* *.diist.tmp.* \\
          *.gpot[0-9].tmp *.grho[0-9].tmp *.opot[0-9].tmp *.orho[0-9].tmp \\
          *.grid.tmp *.hostnames *.int.tmp *.propint.tmp.* *.soscfv.tmp.* \\
          *.vcpcm.tmp *.en.tmp *.C[0-9].tmp *.citations.tmp 2>/dev/null
    rm -f slurm_*.err slurm_*.out 2>/dev/null
fi
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
    ap.add_argument("--first-shell-cut", type=float, default=None,
                    help=f"Override first-shell distance cutoff in Å (module default {FIRST_SHELL_CUT}).")
    args = ap.parse_args()
    carve(args.structure, args.out_dir, args.stem,
          args.site_chain, args.site_resnum,
          first_shell_cut=args.first_shell_cut)
