#!/usr/bin/env python3
"""Recon: enumerate La sites + first-shell donors in 8DQ2 chain A.

Purpose: hand off donor lists to the pre-check pipeline so it knows
which residues to free during sidechain relaxation per metal.
"""
from pathlib import Path
import gemmi
import json

PDB = "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/calibration/pdb_controls/8DQ2.pdb"
CHAIN = "A"
FIRST_SHELL_CUT = 3.2  # Å
FREE_RADIUS = 6.0      # Å — sidechains within this of metal are released

st = gemmi.read_structure(PDB)
st.setup_entities()
model = st[0]

# Build chain A atom list: protein heavy atoms + La ions
chain_a_atoms = []
la_atoms = []  # (resname, seqid, atom)
chain = model[CHAIN]
for residue in chain:
    rname = residue.name
    if rname == "HOH":
        continue
    for atom in residue:
        if atom.element.name == "H":
            continue
        if rname == "LA":
            la_atoms.append((rname, residue.seqid.num, atom))
        else:
            chain_a_atoms.append((rname, residue.seqid.num, atom))

print(f"Chain A: {len(chain_a_atoms)} protein heavy atoms, {len(la_atoms)} La ions")
for rn, sn, atom in la_atoms:
    print(f"  La {sn}: pos=({atom.pos.x:.2f}, {atom.pos.y:.2f}, {atom.pos.z:.2f})")

sites = {}
for rn, la_seq, la_atom in la_atoms:
    metal_pos = la_atom.pos
    # First-shell donors: O/N within 3.2 Å
    donors = []
    free_residues = set()
    for prn, pseq, patom in chain_a_atoms:
        elem = patom.element.name
        d = metal_pos.dist(patom.pos)
        if d <= FREE_RADIUS:
            # Any heavy atom within FREE_RADIUS marks the residue for freeing
            free_residues.add((prn, pseq))
        if elem in ("O", "N") and d <= FIRST_SHELL_CUT:
            donors.append({
                "resname": prn,
                "resid": pseq,
                "atom": patom.name,
                "element": elem,
                "dist_A": round(d, 3),
            })
    donors.sort(key=lambda d: d["dist_A"])
    sites[f"LA{la_seq}"] = {
        "site_id": la_seq,
        "metal_pos_A": [metal_pos.x, metal_pos.y, metal_pos.z],
        "n_first_shell_donors": len(donors),
        "first_shell_donors": donors,
        "free_residues": sorted(list(free_residues), key=lambda x: x[1]),
        "n_free_residues": len(free_residues),
    }
    print(f"\nLA {la_seq}: {len(donors)} first-shell donors (≤{FIRST_SHELL_CUT} Å), {len(free_residues)} residues within {FREE_RADIUS} Å")
    for d in donors:
        print(f"  {d['resname']:3s} {d['resid']:3d} {d['atom']:4s} {d['element']} {d['dist_A']:.2f} Å")

# Save to JSON for downstream pipeline
out_path = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/inputs/8dq2_chainA_sites.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    json.dump(sites, f, indent=2)
print(f"\nSaved to {out_path}")
