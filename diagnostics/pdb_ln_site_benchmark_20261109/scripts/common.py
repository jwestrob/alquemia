"""Shared helpers: file lookup, chain extraction, sequences, alignment, Kabsch."""
import os
import re

import gemmi
import numpy as np

LN_EL = ('La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu')
AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
      'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MSE'}
CIF_DIRS = []


def cif_path(e):
    for d in CIF_DIRS:
        p = os.path.join(d, e.upper() + '.cif')
        if os.path.exists(p):
            return p
    raise FileNotFoundError(e)


def alt_ok(a):
    return a.altloc in ('\0', 'A')


def read(e):
    st = gemmi.read_structure(cif_path(e))
    st.setup_entities()
    st.remove_hydrogens()
    return st


def chain_residues(model, chname):
    """Ordered list of (resname, seqnum, icode, {atomname: xyz}) for standard AA residues of one chain."""
    out = []
    for ch in model:
        if ch.name != chname:
            continue
        for res in ch:
            if res.name not in AA or res.het_flag == 'H' and res.name != 'MSE':
                continue
            atoms = {}
            for a in res:
                if alt_ok(a) and a.name not in atoms:
                    atoms[a.name] = np.array(a.pos.tolist())
            if 'CA' in atoms:
                out.append((res.name, res.seqid.num, res.seqid.icode, atoms))
    return out


def seq_of(residues):
    s = []
    for r in residues:
        t = gemmi.find_tabulated_residue(r[0])
        s.append(t.one_letter_code.upper() if t else 'X')
    return ''.join(s)


def align_map(a, b):
    al = gemmi.align_string_sequences(list(a), list(b), [])
    ident = al.match_count / max(1, min(len(a), len(b)))
    mp, i, j = {}, 0, 0
    for op in re.findall(r'(\d+)([MID])', al.cigar_str()):
        n, t = int(op[0]), op[1]
        for _ in range(n):
            if t == 'M':
                mp[i] = j; i += 1; j += 1
            elif t == 'I':
                i += 1
            else:
                j += 1
    return ident, mp


def kabsch(P, Q):  # R,t minimising |R Q + t - P|
    pc, qc = P.mean(0), Q.mean(0)
    H = (Q - qc).T @ (P - pc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = pc - R @ qc
    rmsd = float(np.sqrt(np.mean(np.sum((Q @ R.T + t - P) ** 2, axis=1))))
    return R, t, rmsd
