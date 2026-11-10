"""Build evaluation units (stripped + partner tiers) per PLAN.md section 1.
Usage: python -I units.py WS SCRIPTS AUDIT_CIF  -> WS/out/units.json, WS/out/partner_kinds.tsv, splits.tsv"""
import csv
import json
import os
import random
import sys
from collections import defaultdict

import numpy as np
import gemmi

WS, SCR, ACIF = sys.argv[1:4]
sys.path.insert(0, SCR)
import common as C  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(WS, 'cif'), ACIF]
LNLIKE = set(C.LN_EL) | {'Y', 'Sc', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'}
rows = [r for r in csv.DictReader(open(os.path.join(WS, 'out', 'ln_sites.tsv')), delimiter='\t')
        if int(r['n_prot_don']) >= 2]
sel = json.load(open(os.path.join(WS, 'out', 'partner_selection.json')))
by_entry = defaultdict(list)
for r in rows:
    by_entry[r['pdb']].append(r)

_cache = {}


def chains_of(e, keep_st=False):
    if e in _cache and not keep_st:
        return None, _cache[e]
    st = C.read(e)
    ch = {}
    for c in st[0]:
        res = C.chain_residues(st[0], c.name)
        if len(res) >= 20:
            ch[c.name] = {'seq': C.seq_of(res), 'ca': np.array([x[3]['CA'] for x in res])}
    _cache[e] = ch
    return (st if keep_st else None), ch


units = []
# ---- stripped tier
for e, rs in sorted(by_entry.items()):
    st, ch = chains_of(e)
    per = defaultdict(list)
    for r in rs:
        if r['chain'] in ch:
            per[r['chain']].append(r)
    for c in sorted(per, key=lambda c: (-len(per[c]), c))[:4]:
        units.append({'uid': f'S|{e}|{c}', 'tier': 'stripped', 'entry': e, 'chain': c, 'kind': 'self',
                      'c95': per[c][0]['c95'], 'c30': per[c][0]['c30'],
                      'sites': [{'sid': f"{e}|{r['site']}", 'xyz': [float(r['qx']), float(r['qy']), float(r['qz'])],
                                 'n_don': int(r['n_prot_don']), 'n_chain': int(r['n_same_chain']), 'cls': r['cls'],
                                 'cluster5': int(r['cluster5']), 'rmsd': 0.0, 'metal_at_site': r['metal']}
                                for r in per[c]]})
# ---- partner tier
src_by_c95 = defaultdict(list)
for r in rows:
    src_by_c95[r['c95']].append(r)
kinds = {}
log = defaultdict(int)
aln = {}
for c95, partners in sorted(sel.items()):
    srcs = src_by_c95[c95]
    for p in sorted(partners):
        try:
            pst, pch = chains_of(p, keep_st=True)
        except Exception:  # noqa: BLE001
            log['partner_read_fail'] += 1; continue
        els = {a.element.name for c in pst[0] for r in c for a in r}
        if els & LNLIKE:
            kinds[p] = 'ln_like_dropped'; log['ln_like_dropped'] += 1; continue
        k = 'ca' if 'Ca' in els else ('other_metal' if any(gemmi.Element(x).is_metal for x in els) else 'metal_free')
        kinds[p] = k
        ns = gemmi.NeighborSearch(pst[0], pst.cell, 5).populate()
        pu = defaultdict(list)
        for r in srcs:
            sst, sch = chains_of(r['pdb'])
            if r['chain'] not in sch:
                continue
            S = sch[r['chain']]
            q = np.array([float(r['qx']), float(r['qy']), float(r['qz'])])
            for pc, PC in pch.items():
                key = (r['pdb'], r['chain'], p, pc)
                if key not in aln:
                    aln[key] = C.align_map(S['seq'], PC['seq'])
                ident, mp = aln[key]
                if ident < 0.95:
                    continue
                near = np.where(np.linalg.norm(S['ca'] - q, axis=1) <= 12.0)[0]
                pairs = [(i, mp[i]) for i in near if i in mp]
                if len(pairs) < 8:
                    log['transfer_lt8'] += 1; continue
                Pp = PC['ca'][[j for _, j in pairs]]
                Qs = S['ca'][[i for i, _ in pairs]]
                R, t, rmsd = C.kabsch(Pp, Qs)
                x = R @ q + t
                met = 'none'
                for mk in ns.find_atoms(gemmi.Position(*x), '\0', radius=2.5):
                    a = mk.to_cra(pst[0]).atom
                    if a.element.is_metal:
                        met = 'Ca' if a.element.name == 'Ca' else 'other'
                        break
                pu[pc].append({'sid': f"{r['pdb']}|{r['site']}", 'xyz': x.tolist(), 'n_don': int(r['n_prot_don']),
                               'n_chain': int(r['n_same_chain']), 'cls': r['cls'], 'cluster5': int(r['cluster5']),
                               'rmsd': round(rmsd, 3), 'metal_at_site': met})
        # merge within 2 A, cap 4 chains
        for pc in sorted(pu, key=lambda c: (-len(pu[c]), c))[:4]:
            cand = sorted(pu[pc], key=lambda s: (-s['n_don'], s['cls'] != 'single-chain', s['sid']))
            kept = []
            for s in cand:
                if all(np.linalg.norm(np.array(s['xyz']) - np.array(k['xyz'])) > 2.0 for k in kept):
                    kept.append(s)
            units.append({'uid': f'P|{p}|{pc}|{c95}', 'tier': 'partner', 'entry': p, 'chain': pc, 'kind': k,
                          'c95': c95, 'c30': srcs[0]['c30'], 'sites': kept})
    print(c95, len(units), flush=True)
json.dump(units, open(os.path.join(WS, 'out', 'units.json'), 'w'))
with open(os.path.join(WS, 'out', 'partner_kinds.tsv'), 'w') as fh:
    for p, k in sorted(kinds.items()):
        fh.write(f'{p}\t{k}\n')
c30 = sorted({u['c30'] for u in units})
random.Random(42).shuffle(c30)
nd = int(round(0.4 * len(c30)))
with open(os.path.join(WS, 'out', 'splits.tsv'), 'w') as fh:
    for i, c in enumerate(c30):
        fh.write(f"{c}\t{'dev' if i < nd else 'test'}\n")
print(dict(log), 'units', len(units), 'stripped', sum(u['tier'] == 'stripped' for u in units),
      'partner', sum(u['tier'] == 'partner' for u in units))
