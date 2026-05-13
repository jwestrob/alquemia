"""Convert Protenix CIF → PDB, add hydrogens at pH 7, write back as PDB."""
import sys
from pathlib import Path

def protonate(src_cif: Path, out_pdb: Path):
    import gemmi
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, Modeller

    # CIF → temp PDB
    st = gemmi.read_structure(str(src_cif))
    tmp_pdb = out_pdb.with_suffix(".tmp.pdb")
    st.write_pdb(str(tmp_pdb))

    # Add Hs
    fixer = PDBFixer(str(tmp_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    with open(out_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    tmp_pdb.unlink(missing_ok=True)
    print(f"protonated: {out_pdb}")

if __name__ == "__main__":
    src = Path(sys.argv[1])
    out = Path(sys.argv[2])
    protonate(src, out)
