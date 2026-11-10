"""Score every unit with v0, B1, B2, B3 (PLAN.md sections 2-3, 5-6).
Usage: python -I score.py WS SCRIPTS AUDIT_CIF SCANNER_ROOT NPROC [UNIT_LIMIT]
Writes WS/out/scores.jsonl (one line per unit)."""
import json
import os
import sys
import tempfile
import time
import traceback
from multiprocessing import Pool

import numpy as np
import gemmi

WS, SCR, ACIF, SROOT, NPROC = sys.argv[1:6]
LIMIT = int(sys.argv[6]) if len(sys.argv) > 6 else None
sys.path.insert(0, SCR)
sys.path.insert(0, SROOT)
import common as C  # noqa: E402
import detect as Dt  # noqa: E402
from on_scanner.scanner import scan  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(WS, 'cif'), ACIF]
TMP = os.path.join(WS, 'tmp')
os.makedirs(TMP, exist_ok=True)
HIT = 4.0


def first_hit(centres, x):
    if len(centres) == 0:
        return None
    d = np.linalg.norm(np.asarray(centres) - np.asarray(x), axis=1)
    w = np.where(d <= HIT)[0]
    return int(w[0]) if len(w) else None


def residue_reach(P, x):
    """Per-residue v0 weight at point x (no clash filter)."""
    S, an, _ = Dt.v0_scores(P, x[None, :])
    out = {}
    for i, rn in enumerate(P['names']):
        w = 0.0
        if rn in Dt.REACH and not P['is_gly'][i]:
            cb, ca = P['anchor'][i], P['ca'][i]
            v = x - cb
            D = np.linalg.norm(v)
            ax = cb - ca
            ok_ang = np.dot(ax, v) / (np.linalg.norm(ax) * max(D, 1e-9)) >= np.cos(np.radians(Dt.MAX_ANGLE))
            for (rmin, rmax, wt) in Dt.REACH[rn]:
                if ok_ang and max(0, rmin - Dt.BAND[1], Dt.BAND[0] - rmax) <= D <= rmax + Dt.BAND[1]:
                    w = max(w, wt)
        O, Cc = P['bbO'][i], P['bbC'][i]
        if np.isfinite(O).all() and np.isfinite(Cc).all():
            d = np.linalg.norm(x - O)
            u1, u2 = Cc - O, x - O
            if Dt.BAND[0] <= d <= Dt.BAND[1] and np.dot(u1, u2) / (np.linalg.norm(u1) * np.linalg.norm(u2)) <= np.cos(np.radians(120)):
                w = max(w, Dt.BB_W)
        out[i] = w
    return float(S[0]), float(an[0]), out


def clash(P, x):
    okO = np.isfinite(P['bbO'][:, 0])
    return bool((np.linalg.norm(P['bbNCAC'] - x, axis=1) < Dt.CLASH['bbNCAC']).any()
                or (np.linalg.norm(P['bbO'][okO] - x, axis=1) < Dt.CLASH['bbO']).any()
                or (np.linalg.norm(P['anchor'][~P['is_gly']] - x, axis=1) < Dt.CLASH['CB']).any())


def run_b1(st, chain, uid):
    s2 = gemmi.Structure()
    s2.cell = st.cell
    m = gemmi.Model('1')
    for ch in st[0]:
        if ch.name != chain:
            continue
        nc = gemmi.Chain(ch.name)
        for r in ch:
            if r.name in C.AA:
                nc.add_residue(r)
        m.add_chain(nc)
    s2.add_model(m)
    s2.setup_entities()
    fd, path = tempfile.mkstemp(suffix='.pdb', dir=TMP)
    os.close(fd)
    s2.write_pdb(path)
    try:
        sites = scan(path)
    finally:
        os.remove(path)
    return [s.probe.tolist() for s in sites]


def work(u):
    t0 = time.time()
    try:
        st = C.read(u['entry'])
        res = C.chain_residues(st[0], u['chain'])
        P = Dt.prepare(res)
        pts = Dt.candidate_points(P)
        S, an, tb = Dt.v0_scores(P, pts)
        k0 = Dt.rank(S, tb, pts)
        c0 = pts[k0]
        n2, tb2 = Dt.b2_scores(P, pts)
        k2 = Dt.rank(n2, tb2, pts)
        c2 = pts[k2]
        g3, k3 = Dt.b3_points(P, u['uid'])
        c3 = g3[k3]
        try:
            c1 = run_b1(st, u['chain'], u['uid'])
            b1_err = None
        except Exception as ex:  # noqa: BLE001
            c1, b1_err = [], repr(ex)
        out = {'uid': u['uid'], 'n_res': len(res), 'n_pts': int(len(pts)),
               'n_cent': {'v0': int(len(c0)), 'b1': len(c1), 'b2': int(len(c2)), 'b3': int(len(c3))},
               'b1_err': b1_err, 'top': {'v0': c0[:10].tolist(), 'b2': c2[:10].tolist()}, 'sites': []}
        for s in u['sites']:
            x = np.array(s['xyz'])
            Sx, anx, contrib = residue_reach(P, x)
            # same-chain residues with an O/N <= 3.0 A of the site position
            near = [i for i, r in enumerate(res)
                    if any(np.linalg.norm(p - x) <= 3.0 for n, p in r[3].items() if n[0] in 'ON')]
            out['sites'].append({'sid': s['sid'], 'rank': {'v0': first_hit(c0, x), 'b1': first_hit(c1, x),
                                                           'b2': first_hit(c2, x), 'b3': first_hit(c3, x)},
                                 'S_true': Sx, 'an_true': anx, 'clash_true': clash(P, x),
                                 'n_near_res': len(near), 'n_near_reach': sum(contrib[i] > 0 for i in near),
                                 'near_res': [f'{res[i][0]}{res[i][1]}:{contrib[i]}' for i in near],
                                 'S_top1': float(S[k0[0]]) if len(k0) else None})
        if u['tier'] == 'stripped':
            ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate()
            sc = [np.array(s['xyz']) for s in u['sites'] if s['cls'] == 'single-chain']
            cap = []
            for j, i in enumerate(k0):
                x = pts[i]
                anyln = False
                for mk in ns.find_atoms(gemmi.Position(*x), '\0', radius=HIT):
                    a = mk.to_cra(st[0]).atom
                    if a.element.name in C.LN_EL and C.alt_ok(a):
                        anyln = True
                        break
                scp = any(np.linalg.norm(x - y) <= HIT for y in sc)
                cap.append([round(float(S[i]), 2), int(an[i]), int(anyln), int(scp)])
            out['cap'] = cap
        out['sec'] = round(time.time() - t0, 2)
        return out
    except Exception:  # noqa: BLE001
        return {'uid': u['uid'], 'error': traceback.format_exc()}


if __name__ == '__main__':
    units = json.load(open(os.path.join(WS, 'out', 'units.json')))
    if LIMIT:
        units = units[:LIMIT]
    outp = os.path.join(WS, 'out', 'scores.jsonl' if not LIMIT else 'scores_test.jsonl')
    # biggest first for load balance
    units.sort(key=lambda u: -os.path.getsize(C.cif_path(u['entry'])))
    with Pool(int(NPROC)) as pool, open(outp, 'w') as fh:
        for i, r in enumerate(pool.imap_unordered(work, units, chunksize=1)):
            fh.write(json.dumps(r) + '\n')
            fh.flush()
            if i % 50 == 0:
                print(i, r.get('uid'), r.get('sec'), 'ERR' if 'error' in r else '', flush=True)
    print('done', len(units))
