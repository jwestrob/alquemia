"""Addition after PLAN.md (lead request): EF-hand loop capacity on compact EF-hand models and reference
lanmodulins. Usage: python -I ef_score.py WS DIAG BENCH_WS AUDIT_CIF
Loops: PROSITE PS00018 as regex (overlapping). Per loop: max n_anion over v0 ranked centres (after NMS)
and over all candidate points within 6 A of any loop atom; loop mean pLDDT (models); at the best ranked
centre, how many loop positions 1,3,5,9,12 that are Asp/Glu reach it."""
import csv
import os
import re
import sys

import gemmi
import numpy as np

WS, DIAG, BWS, ACIF = sys.argv[1:5]
sys.path.insert(0, os.path.join(WS, 'frozen'))
sys.path.insert(0, os.path.join(DIAG, 'scripts'))
import common as C  # noqa: E402
import detect as Dt  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(BWS, 'cif'), ACIF, os.path.join(WS, 'cif')]
PAT = re.compile(r'(?=(D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE]))')
POS = (0, 2, 4, 8, 11)


def reach_w(P, i, x):  # Asp/Glu weight-1.0 reach test, same rule as residue_reach in score.py
    rn = P['names'][i]
    if rn not in ('ASP', 'GLU'):
        return False
    cb, ca = P['anchor'][i], P['ca'][i]
    v = x - cb
    D = np.linalg.norm(v)
    ok = np.dot(cb - ca, v) / (np.linalg.norm(cb - ca) * max(D, 1e-9)) >= np.cos(np.radians(Dt.MAX_ANGLE))
    rmin, rmax, _ = Dt.REACH[rn][0]
    return bool(ok and max(0, rmin - Dt.BAND[1], Dt.BAND[0] - rmax) <= D <= rmax + Dt.BAND[1])


def run(name, res, pl, dataset, source):
    seq = C.seq_of(res)
    P = Dt.prepare(res)
    pts = Dt.candidate_points(P)
    S, an, tb = Dt.v0_scores(P, pts)
    k0 = Dt.rank(S, tb, pts)
    rows = []
    for m in PAT.finditer(seq):
        a = m.start()
        idx = list(range(a, a + 12))
        atoms = np.array([p for i in idx for p in res[i][3].values()])
        d_all = np.min(np.linalg.norm(pts[:, None, :] - atoms[None, :, :], axis=2), axis=1)
        near_all = np.where(d_all <= 6.0)[0]
        near_c = [i for i in k0 if d_all[i] <= 6.0]
        best = max(near_c, key=lambda i: (an[i], S[i])) if near_c else None
        own = [idx[q] for q in POS if res[idx[q]][0] in ('ASP', 'GLU')]
        rows.append({'name': name, 'dataset': dataset, 'source': source, 'len': len(res), 'loop_start': res[a][1],
                     'loop': seq[a:a + 12], 'max_an_centre': int(an[best]) if best is not None else 0,
                     'S_centre': float(S[best]) if best is not None else 0.0,
                     'max_an_point': int(an[near_all].max()) if len(near_all) else 0,
                     'own_acidic_15912': len(own),
                     'own_reaching': sum(reach_w(P, i, pts[best]) for i in own) if best is not None else 0,
                     'loop_plddt': round(float(np.mean([pl[i] for i in idx])), 1) if pl is not None else ''})
    if not rows:
        rows.append({'name': name, 'dataset': dataset, 'source': source, 'len': len(res), 'loop': 'none found'})
    return rows


out = []
for f in sorted(os.listdir(os.path.join(WS, 'models_ef'))) + ['R_MexLanM.pdb', 'R_HansLanM.pdb', 'R_AlLanM.pdb']:
    if not f.endswith('.pdb'):
        continue
    d = 'models_ef' if f.startswith('E_') else 'models'
    st = gemmi.read_structure(os.path.join(WS, d, f))
    st.remove_hydrogens()
    ch = st[0][0].name
    res = C.chain_residues(st[0], ch)
    pl = [r.find_atom('CA', '*').b_iso for r in st[0][ch] if r.name in C.AA]
    nm = f[:-4]
    ds = nm.split('_')[1] if nm.startswith('E_') else 'reference'
    out += run(nm, res, pl, ds, 'model')
for nm, e in (('R_MexLanM', '8FNS'), ('R_HansLanM', '8DQ2'), ('R_AlLanM', '9PI2')):
    st = C.read(e)
    out += run(nm, C.chain_residues(st[0], 'A'), None, 'reference', f'crystal {e}:A')
cols = ['name', 'dataset', 'source', 'len', 'loop_start', 'loop', 'max_an_centre', 'S_centre', 'max_an_point',
        'own_acidic_15912', 'own_reaching', 'loop_plddt']
with open(os.path.join(DIAG, 'results', 'ef_loops.tsv'), 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in out:
        fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
        print('\t'.join(str(r.get(c, '')) for c in cols))
