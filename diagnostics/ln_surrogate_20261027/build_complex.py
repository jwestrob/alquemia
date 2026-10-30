"""Build La(III) complexes of one ligand with Architector (GFN2-xTB via tblite), at total CN 8 and 9
(waters fill CN - denticity). Fully deprotonated acids; donors = one O per carboxylate, aliphatic amine
N, ether O, alcohol O (capped at CN, dropping carboxylates last-first, then alcohols/ethers).
Also builds the aqua ions when lig_id == "aqua".

usage: build_complex.py <lig_id> <smiles_neutral> <outdir>
writes <outdir>/La_cn<CN>.xyz (comment line: charge)
"""
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from architector import build_complex

ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def deprotonate(smi):
    m = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol(m)
    for match in m.GetSubstructMatches(ACID):
        o = rw.GetAtomWithIdx(match[2]); o.SetFormalCharge(-1); o.SetNumExplicitHs(0); o.SetNoImplicit(True)
    mol = rw.GetMol(); Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def donors(smi):
    m = Chem.MolFromSmiles(smi)
    carbox = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() == -1]
    amine = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N" and not a.GetIsAromatic()
             and not any(n.GetSymbol() in "CS" and any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() in "OS"
                                                       for b in n.GetBonds()) for n in a.GetNeighbors())]
    ether = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() == 0 and a.GetDegree() == 2
             and not a.GetIsAromatic() and all(n.GetSymbol() == "C" and not any(b.GetBondTypeAsDouble() == 2 for b in n.GetBonds()) for n in a.GetNeighbors())]
    alcohol = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 1
               and a.GetNeighbors()[0].GetHybridization() == Chem.HybridizationType.SP3]
    return carbox, amine, ether, alcohol


def chelatable(smi, cand):
    """Greedy donor set in priority order: a donor joins if it is <= 5 bonds from an accepted donor
    (chelate ring <= 7 members including the metal); the first candidate always joins."""
    D = Chem.GetDistanceMatrix(Chem.MolFromSmiles(smi))
    acc = []
    for k in cand:
        if not acc or min(D[k][j] for j in acc) <= 5: acc.append(k)
    return acc


def main():
    lig, smi, out = sys.argv[1], sys.argv[2], Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
    for cn in (8, 9):
        f = out / f"La_cn{cn}.xyz"
        if f.exists(): continue
        if lig == "aqua":
            ligs, q = [{"smiles": "O", "coordList": [0]}] * cn, 3
        else:
            s = deprotonate(smi)
            c, a, e, h = donors(s)
            dl = chelatable(s, a + c + e + h[:2])[:cn]  # priority: amines, carboxylates, ethers, <=2 alcohols
            q = 3 + Chem.GetFormalCharge(Chem.MolFromSmiles(s))
            ligs = [{"smiles": s, "coordList": dl}] + [{"smiles": "O", "coordList": [0]}] * (cn - len(dl))
            if cn - len(dl) < 0: continue
        inp = {"core": {"metal": "La", "coreCN": cn}, "ligands": ligs,
               "parameters": {"metal_ox": 3, "full_method": "GFN2-xTB", "assemble_method": "GFN2-xTB",
                              "n_conformers": 1, "return_only_1": True, "relax": True}}
        try:
            res = build_complex(inp)
        except Exception as ex:
            print(f"{lig} cn{cn} FAILED {ex!r}"); continue
        if not res:
            print(f"{lig} cn{cn} no structure"); continue
        at = res[list(res)[0]]["ase_atoms"]
        at.write(f, format="xyz", comment=f"charge={q}")
        print(f"{lig} cn{cn} ok q={q} n={len(at)}")


if __name__ == "__main__":
    main()
